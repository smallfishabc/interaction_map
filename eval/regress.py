"""Regression harness: the combined branch must reproduce the August 2023 p53 output.

Reference data (never written to, only read):
    All_atom_p53_interaction/p53-summary/BB/S_0/p53_0.8_contact_df_1201.csv
    All_atom_p53_interaction/p53-summary/BB/S_0/interaction_0317.csv

Two levels:

    python regress.py              fast. Replays normalization and the arc selection on the
                                   stored contact probabilities. No mdtraj, no trajectory.
                                   Run this after every refactoring step.

    python regress.py --full       also recomputes the contact map from the 5 repeats, on a
                                   COPY of the trajectory in scratch, and compares it too.
                                   Slow. Run before declaring a step finished.

Exit code 0 means every checked column matched exactly.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
DATA = os.path.join(PROJECT, 'All_atom_p53_interaction', 'p53-summary')
REF_DIR = os.path.join(DATA, 'BB', 'S_0')
# Working copies live under scratch/_work, which is not tracked. The reference folder is
# only ever read: --full copies the trajectory out before recomputing, because building a
# contact map writes a csv next to the trajectory it came from.
WORK = os.path.join(HERE, '_work', 'regress')
COPY_DIR = os.path.join(WORK, 'p53-summary', 'BB', 'S_0')

sys.path.insert(0, os.path.join(PROJECT, 'src'))

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd

from ensemble_map import profiles
from ensemble_map import normalization as nl

FAIL = []
FULL = '--full' in sys.argv

# The stored reference was produced in August 2023 with a different numpy/libm build. Its
# pow() results differ from this machine's in the last couple of digits: gs_standard is
# bit-identical on only 9.8% of rows, with a maximum RELATIVE difference of 1.1e-14, and
# the file cannot even regenerate its own relative_strength column from its own
# gs_standard column (max 5.3e-13). Every formulation of the power -- row-wise apply,
# Series **, array **, a Python loop -- reproduces the file equally well, so this is
# environmental and not something the code controls.
#
# So: continuous columns are compared with a tolerance far below anything that could
# change a result, while the things that decide the figure -- interaction level, arc
# count, arc endpoints -- must match exactly.
RTOL = 1e-11


def check(label, ok, detail=''):
    print('  {0} {1}{2}'.format('ok  ' if ok else 'FAIL', label,
                                '   ' + str(detail) if detail else ''))
    if not ok:
        FAIL.append(label)


def section(title):
    print('\n== {0} =='.format(title))


def read_sequence():
    with open(os.path.join(DATA, 'seq.fasta')) as handle:
        return handle.read().strip()


SEQ = read_sequence()
REF = pd.read_csv(os.path.join(REF_DIR, 'interaction_0317.csv'), index_col=0)
REF_CONTACT = pd.read_csv(os.path.join(REF_DIR, 'p53_0.8_contact_df_1201.csv'), index_col=0)

print('reference: {0} pairs, sequence {1} residues'.format(len(REF), len(SEQ)))


# ------------------------------------------------------------------ numbering
section('residue numbering')
from ensemble_map import contact_map_generation as cg

check('r_1 min is 1', int(REF['r_1'].min()) == 1, int(REF['r_1'].min()))
check('r_2 max is len(seq)', int(REF['r_2'].max()) == len(SEQ),
      '{0} vs {1}'.format(int(REF['r_2'].max()), len(SEQ)))
check('check_residue_numbering accepts the reference',
      cg.check_residue_numbering(REF, SEQ) is True)
# A deliberately shifted table must be rejected, so the guard is not vacuous.
shifted = REF.copy()
shifted['r_1'] = shifted['r_1'] - 1
shifted['r_2'] = shifted['r_2'] - 1
check('check_residue_numbering rejects a shifted table',
      cg.check_residue_numbering(shifted, SEQ) is False)


# ------------------------------------------------------------------ normalization
section('normalization replayed on the stored contact probabilities')
replay = REF[['r_1', 'r_2', 'cont_prob', 'distance']].copy()
out = nl.normalization(replay, profile='all_atom')

for col in ('gs_standard', 'relative_strength'):
    got, want = out[col].values, REF[col].values
    rel = np.abs(got - want) / np.maximum(np.abs(want), 1e-300)
    check('{0} within {1:g} relative'.format(col, RTOL), rel.max() <= RTOL,
          'max rel diff {0:.3e}'.format(rel.max()))

same = out['plot_value'].values == REF['plot_value'].values
check('plot_value identical', same.all(), '{0} of {1} rows differ'.format((~same).sum(), len(same)))
if not same.all():
    bad = out.loc[~same, ['r_1', 'r_2', 'distance', 'cont_prob', 'relative_strength', 'plot_value']].copy()
    bad['reference'] = REF.loc[~same, 'plot_value'].values
    print(bad.head(20).to_string())

print('  interaction level counts (replay vs reference):')
for lvl in (2, 1, 0, -1, -2):
    print('    {0:>2}: {1:5d}  {2:5d}'.format(lvl, int((out['plot_value'] == lvl).sum()),
                                              int((REF['plot_value'] == lvl).sum())))

# The reference curve constants must still be the calibrated ones.
aa = profiles.get_profile('all_atom')
check('reference curve is a=1.64 b=-1.33', (aa.reference_a, aa.reference_b) == (1.64, -1.33),
      (aa.reference_a, aa.reference_b))
check('cutoffs are (1.5, 0.5, -1, -2)', tuple(aa.inter_cutoff) == (1.5, 0.5, -1, -2), aa.inter_cutoff)
check('contact cutoff is 0.8', aa.cutoff == 0.8, aa.cutoff)


# ------------------------------------------------------------------ arc selection
section('arc selection and line widths')
from ensemble_map import interaction_plot as ip

graphg = ip.create_network(SEQ, len(SEQ))
pos, layout = ip.create_position(graphg)

# Reproduce, independently of the plotting code, which arcs the old code would draw and how
# thick. The old all-atom code filtered on `distance <= 2: continue`, which never fired
# because mdtraj already drops |i-j| < 3, and used a 2x line width scale.
expected = {}
for lvl in (2, 1, -1, -2):
    sel = REF[(REF['plot_value'] == lvl) & (REF['distance'] > 2)]
    widths = []
    for _, row in sel.iterrows():
        if lvl > 0:
            widths.append(2 * row['cont_prob'] * (1 + row['relative_strength']))
        else:
            widths.append(-2 * row['cont_prob'] * (row['relative_strength'] - 1))
    expected[lvl] = (list(zip(sel['r_1'].astype(int), sel['r_2'].astype(int))), widths)


class RecordingAxis(object):
    """Captures annotate() calls instead of rendering, so this runs without a GUI stack."""

    def __init__(self):
        self.calls = []

    def annotate(self, _text, xy=None, xytext=None, arrowprops=None, **_kw):
        self.calls.append((xy, xytext, arrowprops['lw'], arrowprops['color']))


for lvl in (2, 1, -1, -2):
    rec = RecordingAxis()
    ip.interaction_plotting(out, layout, rec, lvl, profile='all_atom')
    pairs, widths = expected[lvl]
    check('level {0:>2}: arc count'.format(lvl), len(rec.calls) == len(pairs),
          '{0} vs {1}'.format(len(rec.calls), len(pairs)))
    if len(rec.calls) == len(pairs):
        # Endpoint x coordinates must land on the residues the pair names.
        xs_ok = True
        for (xy, xytext, _lw, _c), (r1, r2) in zip(rec.calls, pairs):
            if abs(xy[0] - (layout[r1][0] - 0.2)) > 1e-9 or abs(xytext[0] - (layout[r2][0] + 0.2)) > 1e-9:
                xs_ok = False
                break
        check('level {0:>2}: endpoints on the named residues'.format(lvl), xs_ok)
        got_w = np.array([c[2] for c in rec.calls])
        want_w = np.array(widths)
        rel = np.abs(got_w - want_w) / np.maximum(np.abs(want_w), 1e-300)
        check('level {0:>2}: line widths within {1:g} relative'.format(lvl, RTOL),
              rel.max() <= RTOL, 'max rel diff {0:.3e}'.format(rel.max()))


# ------------------------------------------------------------------ rendered figure
section('rendered figure vs p53.png')
# The strongest check available: render from the stored interaction values and compare
# against the committed PNG pixel by pixel. Requires the letter residue style, which is
# what produced the published figures.
import tempfile
import matplotlib.image as mpimg

fig_dir = tempfile.mkdtemp(prefix='imfig_')
here = os.getcwd()
os.chdir(fig_dir)
letter = profiles.Profile(**dict(profiles.ALL_ATOM.__dict__, residue_style='letter'))
ip.interaction_map(SEQ, len(SEQ), REF.copy(), 'p53', profile=letter, show=False, save=True)
os.chdir(here)

ref_png = mpimg.imread(os.path.join(REF_DIR, 'p53.png'))
new_png = mpimg.imread(os.path.join(fig_dir, 'p53.png'))
check('same image size', ref_png.shape == new_png.shape,
      '{0} vs {1}'.format(ref_png.shape, new_png.shape))
if ref_png.shape == new_png.shape:
    diff = np.abs(ref_png.astype(float) - new_png.astype(float))
    check('pixel identical to p53.png', diff.max() == 0.0,
          '{0} pixels differ'.format(int((diff.max(axis=2) > 0).sum())))


# ------------------------------------------------------------------ full recompute
def ensure_trajectory_copy():
    """Copy the topology and the five repeats into the working folder.

    Recomputing writes a contact csv beside the trajectory, so this runs on a copy and the
    reference folder is left untouched. About 100 MB, copied once.
    """
    import shutil
    if os.path.isdir(COPY_DIR) and os.path.exists(os.path.join(COPY_DIR, '__traj_4.xtc')):
        return True
    if not os.path.exists(os.path.join(REF_DIR, '__traj_0.xtc')):
        return False
    os.makedirs(COPY_DIR, exist_ok=True)
    print('  copying the trajectory to {0} (once, about 100 MB)'.format(COPY_DIR))
    for name in os.listdir(REF_DIR):
        if name.startswith('__START_0') or name.startswith('__traj_'):
            shutil.copy(os.path.join(REF_DIR, name), COPY_DIR)
    return True


if FULL:
    section('contact map recomputed from the 5 repeats')
    if not ensure_trajectory_copy():
        print('  SKIP: trajectory not present in {0}'.format(REF_DIR))
    else:
        os.chdir(COPY_DIR)
        contact = cg.generate_contact('p53', '__START_0.pdb', 5, profile='all_atom')
        new = contact.contact
        check('same row count', len(new) == len(REF_CONTACT),
              '{0} vs {1}'.format(len(new), len(REF_CONTACT)))
        if len(new) == len(REF_CONTACT):
            for col in ('r_1', 'r_2', 'distance'):
                check('{0} identical'.format(col),
                      (new[col].values == REF_CONTACT[col].values).all())
            # cont_prob is a mean over ~28000 booleans. The count is integral, but the
            # summation order inside numpy's mean differs between builds, so this lands
            # one ULP away rather than exactly.
            got, want = new['cont_prob'].values, REF_CONTACT['cont_prob'].values
            rel = np.abs(got - want) / np.maximum(np.abs(want), 1e-300)
            check('cont_prob within {0:g} relative'.format(RTOL), rel.max() <= RTOL,
                  'max rel diff {0:.3e}'.format(rel.max()))
else:
    print('\n(skipping the trajectory recompute; pass --full to include it)')


# ------------------------------------------------------------------ summary
section('summary')
if FAIL:
    print('{0} FAILURE(S): {1}'.format(len(FAIL), ', '.join(FAIL)))
    sys.exit(1)
print('matches the August 2023 reference')
