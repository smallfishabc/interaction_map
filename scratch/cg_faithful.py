"""The CG profile must reproduce the 0519_CG branch exactly.

There is no stored CG output to regress against, so this regenerates one: it checks the
original branch's modules out of git into a temp folder, runs them on a CALVADOS mutant,
runs the combined code on the same mutant, and compares the contact map, the interaction
table and the rendered figure.

    python cg_faithful.py [MUTANT]        default A39G_K24Y

Run it after any change that touches the CG path.
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
CG_DATA = os.path.join(PROJECT, 'CG_p53_random')
LEGACY_BRANCH = '0519_CG'
MODULES = ['contact_map_generation.py', 'normalization.py', 'interaction_plot.py',
           'default_function.py', 'readpath.py', 'showoff.py']

NAME = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('-') else 'A39G_K24Y'

sys.path.insert(0, PROJECT)

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.image as mpimg

FAIL = []


def check(label, ok, detail=''):
    print('  {0} {1}{2}'.format('ok  ' if ok else 'FAIL', label, '   ' + str(detail) if detail else ''))
    if not ok:
        FAIL.append(label)


def read_fasta(path):
    with open(path) as handle:
        return ''.join(l.strip() for l in handle if not l.startswith('>'))


DRIVER = '''
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import default_function
NAME = {name!r}
with open(NAME + '.fasta') as h:
    seq = ''.join(l.strip() for l in h if not l.startswith('>'))
here = os.getcwd()
default_function.interaction_map_pairwise(NAME, here, seq, here,
                                          pdb_top=NAME + '.pdb',
                                          xtc_input=NAME + '.dcd',
                                          read_from_file=False)
'''

root = tempfile.mkdtemp(prefix='cgfaith_')
legacy_dir = os.path.join(root, 'legacy')
new_dir = os.path.join(root, 'new')
os.makedirs(legacy_dir)
os.makedirs(new_dir)

print('mutant     ', NAME)
print('legacy from', LEGACY_BRANCH)
print('workspace  ', root)

for ext in ('.pdb', '.dcd', '.fasta'):
    src = os.path.join(CG_DATA, NAME + ext)
    if not os.path.exists(src):
        print('missing input:', src)
        sys.exit(2)
    shutil.copy(src, legacy_dir)
    shutil.copy(src, new_dir)

for mod in MODULES:
    out = subprocess.check_output(['git', 'show', '{0}:{1}'.format(LEGACY_BRANCH, mod)],
                                  cwd=PROJECT)
    with open(os.path.join(legacy_dir, mod), 'wb') as handle:
        handle.write(out)

with open(os.path.join(legacy_dir, 'run_legacy.py'), 'w') as handle:
    handle.write(DRIVER.format(name=NAME))

print('\nrunning the {0} modules...'.format(LEGACY_BRANCH))
subprocess.check_call([sys.executable, 'run_legacy.py'], cwd=legacy_dir,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print('running the combined code...')
import profiles
import contact_map_generation as cg
import normalization as nl
import interaction_plot as ip

seq = read_fasta(os.path.join(new_dir, NAME + '.fasta'))
prof = profiles.get_profile('cg')
here = os.getcwd()
os.chdir(new_dir)
contact = cg.generate_contact(NAME, NAME + '.pdb', NAME + '.dcd', profile=prof)
inter = nl.normalization(contact.contact, profile=prof)
inter.to_csv('interaction_0317.csv')
ip.interaction_map(seq, len(seq), inter, NAME, profile=prof, show=False, save=True)
os.chdir(here)

contact_csv = '{0}_1.2_contact_df_1201.csv'.format(NAME)

def compare_tables(label, legacy_path, new_path):
    """Every column the legacy run produced must still be there and unchanged.

    The combined code adds columns (residue identity, repeat spread). Those are new by
    design, so they are listed rather than compared.
    """
    a = pd.read_csv(legacy_path, index_col=0)
    b = pd.read_csv(new_path, index_col=0)
    print('\n== {0} =='.format(label))
    check('every legacy column still present',
          all(col in b.columns for col in a.columns),
          [col for col in a.columns if col not in b.columns] or 'ok')
    check('same row count', len(a) == len(b), '{0} vs {1}'.format(len(a), len(b)))
    if len(a) == len(b):
        for col in a.columns:
            if col in b.columns:
                check('{0} bit-identical'.format(col), bool((a[col].values == b[col].values).all()))
    added = [col for col in b.columns if col not in a.columns]
    if added:
        print('  added by the combined code: {0}'.format(', '.join(added)))


compare_tables('contact map',
               os.path.join(legacy_dir, contact_csv),
               os.path.join(new_dir, contact_csv))
compare_tables('interaction table',
               os.path.join(legacy_dir, 'interaction_0317.csv'),
               os.path.join(new_dir, 'interaction_0317.csv'))

print('\n== rendered figure ==')
ia = mpimg.imread(os.path.join(legacy_dir, NAME + '.png'))
ib = mpimg.imread(os.path.join(new_dir, NAME + '.png'))
check('same image size', ia.shape == ib.shape, '{0} vs {1}'.format(ia.shape, ib.shape))
if ia.shape == ib.shape:
    d = np.abs(ia.astype(float) - ib.astype(float))
    check('pixel identical', d.max() == 0.0, '{0} pixels differ'.format(int((d.max(axis=2) > 0).sum())))

print('\n== summary ==')
if FAIL:
    print('{0} FAILURE(S): {1}'.format(len(FAIL), ', '.join(FAIL)))
    sys.exit(1)
print('the CG profile reproduces {0} exactly'.format(LEGACY_BRANCH))
shutil.rmtree(root, ignore_errors=True)
