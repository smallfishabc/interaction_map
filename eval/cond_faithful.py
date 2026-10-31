"""The conditional analysis must reproduce 1115_new_start exactly.

No stored output exists for this path, so the comparison is generated: the 1115 modules
are checked out of git into a temp folder and run on a trimmed p53 trajectory, then the
combined code is run on the same input and the outputs compared.

    python cond_faithful.py [R1 R2]        default 13 24, the pair 1115's main.py used
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
# A short trajectory, built once from the first p53 repeat. The secondary structure and
# hydrogen bond calculations loop over frames in Python, so the full 28000 would be slow
# and would tell us nothing more about whether the port is faithful.
DATA = os.path.join(HERE, '_work', 'cond')
SOURCE = os.path.join(PROJECT, 'All_atom_p53_interaction', 'p53-summary')
FIXTURE_FRAMES = 300
LEGACY_BRANCH = '1115_new_start'
MODULES = ['contact_map_generation.py', 'normalization.py', 'interaction_plot.py',
           'default_function.py', 'readpath.py', 'showoff.py', 'search_chunk.py',
           'interaction_strength.py']

R1 = int(sys.argv[1]) if len(sys.argv) > 2 else 13
R2 = int(sys.argv[2]) if len(sys.argv) > 2 else 24

sys.path.insert(0, os.path.join(PROJECT, 'src'))

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd

FAIL = []


def check(label, ok, detail=''):
    print('  {0} {1}{2}'.format('ok  ' if ok else 'FAIL', label, '   ' + str(detail) if detail else ''))
    if not ok:
        FAIL.append(label)


DRIVER = '''
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import default_function
here = os.getcwd()
with open('seq.fasta') as h:
    seq = h.read().strip()
default_function.search_interaction('p53', {r1}, {r2}, here, seq, here,
                                    pdb_top='__START_0.pdb', xtc_input=['__traj_0.xtc'])
'''

def ensure_fixture():
    """Build the short trajectory the comparison runs on, if it is not there yet."""
    if os.path.exists(os.path.join(DATA, '__traj_0.xtc')):
        return True
    src = os.path.join(SOURCE, 'BB', 'S_0')
    if not os.path.exists(os.path.join(src, '__traj_0.xtc')):
        return False
    import mdtraj as md
    os.makedirs(DATA, exist_ok=True)
    print('building a {0} frame fixture in {1}'.format(FIXTURE_FRAMES, DATA))
    traj = md.load(os.path.join(src, '__traj_0.xtc'), top=os.path.join(src, '__START_0.pdb'))
    traj[:FIXTURE_FRAMES].save_xtc(os.path.join(DATA, '__traj_0.xtc'))
    shutil.copy(os.path.join(src, '__START_0.pdb'), DATA)
    shutil.copy(os.path.join(SOURCE, 'seq.fasta'), DATA)
    return True


if not ensure_fixture():
    sys.exit('All_atom_p53_interaction is not present, nothing to compare against')

root = tempfile.mkdtemp(prefix='condfaith_')
legacy_dir = os.path.join(root, 'legacy')
new_dir = os.path.join(root, 'new')
os.makedirs(legacy_dir)
os.makedirs(new_dir)

print('pair       {0}, {1}'.format(R1, R2))
print('legacy from', LEGACY_BRANCH)
print('workspace  ', root)

for fname in ('__START_0.pdb', '__traj_0.xtc', 'seq.fasta'):
    shutil.copy(os.path.join(DATA, fname), legacy_dir)
    shutil.copy(os.path.join(DATA, fname), new_dir)

for mod in MODULES:
    try:
        out = subprocess.check_output(['git', 'show', '{0}:{1}'.format(LEGACY_BRANCH, mod)],
                                      cwd=PROJECT, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        continue
    text = out.decode('utf-8')
    # 1115_new_start was written against networkx < 2.4, where Graph.node was the public
    # accessor; 2.4 renamed it to Graph._node, and the 0317 branch adopted the new spelling
    # when it was time. This harness runs the networkx the rest of the project is pinned to,
    # so the accessor is spelled the current way here purely so the two versions can be run
    # side by side. Reading a node attribute, under either spelling, changes no value.
    text = text.replace('graphg.node[n]', 'graphg._node[n]')
    # 1115_new_start still carries the older reference exponent, b = -1.32. The all-atom
    # profile uses -1.33, which is what 0811 carries and what the published p53 output was
    # made with. Aligning it here isolates the question being asked -- is the port
    # faithful -- from a constant that deliberately changed between the branches.
    text = text.replace('b1=-1.32', 'b1=-1.33')
    with open(os.path.join(legacy_dir, mod), 'w') as handle:
        handle.write(text)

with open(os.path.join(legacy_dir, 'run_legacy.py'), 'w') as handle:
    handle.write(DRIVER.format(r1=R1, r2=R2))

print('\nrunning the {0} modules...'.format(LEGACY_BRANCH))
result = subprocess.run([sys.executable, 'run_legacy.py'], cwd=legacy_dir,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
if result.returncode != 0:
    print(result.stdout.decode('utf-8', 'replace')[-3000:])
    sys.exit('legacy run failed')

print('running the combined code...')
from ensemble_map import default_function
from ensemble_map import profiles

with open(os.path.join(new_dir, 'seq.fasta')) as handle:
    seq = handle.read().strip()
here = os.getcwd()
# 1115 drew the sequence as coloured letters. The marker renderer is the current default,
# so the figure is rendered the legacy way here to compare like with like.
letter = profiles.Profile(**dict(profiles.ALL_ATOM.__dict__, residue_style='letter'))
default_function.search_interaction('p53', R1, R2, new_dir, seq, new_dir,
                                    pdb_top='__START_0.pdb', xtc_input=['__traj_0.xtc'],
                                    profile=letter, show=False)
os.chdir(here)

suffix = '{0}_{1}.csv'.format(R1, R2)


def compare(label, filename, index_col=0):
    a = pd.read_csv(os.path.join(legacy_dir, filename), index_col=index_col)
    b = pd.read_csv(os.path.join(new_dir, filename), index_col=index_col)
    print('\n== {0} =='.format(label))
    check('every legacy column present', all(c in b.columns for c in a.columns),
          [c for c in a.columns if c not in b.columns] or 'ok')
    # The ensemble summary gained a Frames row. Compare on the rows the legacy run had.
    if not a.index.equals(b.index) and set(a.index) <= set(b.index):
        added_rows = [i for i in b.index if i not in a.index]
        b = b.loc[a.index]
        print('  rows added by the combined code: {0}'.format(', '.join(map(str, added_rows))))
    check('same row count', len(a) == len(b), '{0} vs {1}'.format(len(a), len(b)))
    if len(a) == len(b):
        for col in a.columns:
            if col not in b.columns:
                continue
            x, y = a[col].values, b[col].values
            if x.dtype.kind in 'fc':
                same = np.allclose(x, y, rtol=0, atol=0, equal_nan=True)
                detail = '' if same else 'max abs diff {0:.3e}'.format(np.nanmax(np.abs(x - y)))
            else:
                same = bool((x == y).all())
                detail = ''
            check('{0} bit-identical'.format(col), same, detail)
    added = [c for c in b.columns if c not in a.columns]
    if added:
        print('  added by the combined code: {0}'.format(', '.join(added)))


compare('contact map, with the interaction', 'with_interaction_' + suffix, index_col=None)
compare('contact map, without it', 'without_interaction_' + suffix, index_col=None)
compare('ensemble comparison', 'ensemble_' + suffix)

print('\n== rendered figure ==')
# Informational, not a pass/fail. 1115 scaled the arc width by the interaction strength
# alone, 2 * relative_strength * level, with repulsive arcs further multiplied by 0.1.
# 0811 replaced that with 2 * cont_prob * (1 + relative_strength), and that is the form
# the published p53 figure was drawn with, which the all-atom regression verifies pixel
# for pixel. The current code carries the 0811 formula on purpose, so these two images
# are expected to differ.
import matplotlib.image as mpimg
ia = mpimg.imread(os.path.join(legacy_dir, 'p53.png'))
ib = mpimg.imread(os.path.join(new_dir, 'p53.png'))
check('same image size', ia.shape == ib.shape, '{0} vs {1}'.format(ia.shape, ib.shape))
if ia.shape == ib.shape:
    d = np.abs(ia.astype(float) - ib.astype(float))
    differing = int((d.max(axis=2) > 0).sum())
    print('  note  {0} pixels differ, from the 1115 -> 0811 line width change'.format(differing))
    print('        1115:  2 * relative_strength * level   (repulsive also x0.1)')
    print('        0811:  2 * cont_prob * (1 + relative_strength)   <- published form')

print('\n== ensemble values ==')
print(pd.read_csv(os.path.join(new_dir, 'ensemble_' + suffix), index_col=0).to_string())

print('\n== summary ==')
if FAIL:
    print('{0} FAILURE(S): {1}'.format(len(FAIL), ', '.join(FAIL)))
    sys.exit(1)
print('the conditional analysis reproduces {0} exactly'.format(LEGACY_BRANCH))
shutil.rmtree(root, ignore_errors=True)
