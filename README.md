# Interaction Map

Author: Feng Yu

Interaction maps for intrinsically disordered regions. The method takes a simulation
ensemble, measures how often each pair of residues is in contact, compares that against
what a generic flexible chain of the same sequence separation would do, and draws the
excess and the deficit as arcs along the sequence.

Four operations, and everything here is scaffolding around them:

1. trajectory → Cα–Cα contact probability per residue pair, with sequence separation
   `s = j − i`
2. divide by a reference polymer curve `p_ref(s) = a·s^b` and take the log, giving the
   interaction strength `I = ln(p_obs / p_ref)`
3. bin that into five levels and draw a sequence track with arcs, attractive above,
   repulsive below
4. optionally condition the ensemble on one contact being formed and compare the two
   sub-ensembles

## Profiles

All-atom and coarse-grained analyses use different calibrated constants. Both live in
`profiles.py` and are selected with `--profile`; neither can be changed by editing the
other.

| | `all_atom` | `cg` |
|---|---|---|
| simulation | CAMPARI / ABSINTH | CALVADOS |
| contact cutoff | 0.8 nm | 1.2 nm |
| pair selection | `compute_contacts(contacts='all', scheme='CA')`, drops \|i−j\| < 3 | `select_pairs('all','all')`, keeps every pair |
| reference curve | 1.64·s^−1.33 | 13.12·s^−2.32 |
| interaction cutoffs | (1.5, 0.5, −1, −2) | (2, 1, −1, −2) |
| minimum separation drawn | 3 | 5 |

The coarse-grained reference falls off faster, which is why its exponent is steeper.

## Layout

```
src/interaction_map/   the package
scripts/main.py        run it without installing
tests/                 the test suite, with a small coarse-grained fixture
eval/                  checks that the code still reproduces published output
```

## Installing

```
conda env create -f environment.yml
conda activate interaction-map
pip install -e . --no-deps
```

`--no-deps` because mdtraj and the rest come from conda; mdtraj publishes no Windows
wheels for 1.9.x, so pip cannot supply them. Installing gives you the `interaction-map`
command. Skipping the install is fine too — `python scripts/main.py` takes the same
arguments and puts `src/` on the path itself.

mdtraj publishes no Windows wheels for 1.9.x, so this has to be a conda environment
rather than a pip virtualenv. **On Windows the environment must be activated**, not just
invoked by path: running `<env>/python.exe` directly leaves the environment's
`Library/bin` off the DLL search path, and the compiled extensions then fail to load,
killing the interpreter with no traceback.

## Running

Multi-repeat simulation database layout:

```
<name>-summary/
  seq.fasta              protein sequence
  BB/                    residue group whose transfer free energy is altered (--restype)
    S_0/                 psi value of the transfer free energy (--psi)
      __START_0.pdb      topology
      __traj_0.xtc ...   one file per repeat
```

```
python scripts/main.py -dir "path/to/<name>-summary" --profile all_atom
```

A single topology and trajectory instead:

```
python scripts/main.py --single_traj -dir "path/to/folder" -p protein.pdb -x protein.dcd -name protein --profile cg
```

Several trajectories of the same construct are treated as repeats, and the spread between
them is reported as `p_sem`:

```
python scripts/main.py -p top.pdb -x run1.xtc run2.xtc run3.xtc -name protein
```

Useful flags: `--psi`, `--restype`, `--repeat`, `--name`, `--out`, `--seq`, `--no_show`.

## Conditional analysis

Split the ensemble on whether one interaction is formed, and compare the two halves:

```
python scripts/main.py -p top.pdb -x traj.xtc -name protein --search 13 24
```

A frame counts as containing the interaction when any residue pair in the 3x3 window
around (13, 24) is within the contact cutoff. `--window N` changes the half width; `0`
narrows it to the single pair. Writes `ensemble_13_24.csv` with Rg, end-to-end distance,
helicity, beta content, hydrogen bonds and the frame count for each side, plus a contact
table for each.

Outputs land next to the trajectory: `<name>_<cutoff>_contact_df_1201.csv` with the
contact probabilities, `interaction_0317.csv` with the interaction strengths and levels,
and `<name>.png` / `.svg`.

Folders used for previous figures are listed in `RUN_PATHS.md`.

## Residue numbering

`r_1` and `r_2` in every output are **1-based sequence positions**. The two profiles get
there differently, and the difference matters:

- All-atom CAMPARI topologies carry an ACE cap at residue index 0 and an NME cap at the
  end. `ignore_nonprotein=True` drops the caps from the pair list but not from the
  numbering, so the first real residue is already index 1 and no shift is applied.
- CALVADOS topologies have no caps, so one is added.

An **uncapped all-atom topology would therefore come out one residue low**, with every arc
drawn against the wrong residue and no error raised.
`contact_map_generation.check_residue_numbering` warns when the numbering does not span
`1..len(sequence)`, and `check_sequence_against_topology` warns when the sequence file
does not match the topology. Both are diagnostic and alter nothing.

## Modules

| `src/interaction_map/` | Contents |
|---|---|
| `cli.py` | command line entry point |
| `default_function.py` | the pipeline: contacts → normalization → figure |
| `profiles.py` | the two sets of calibrated constants |
| `contact_map_generation.py` | trajectory loading, contact probabilities, numbering guards |
| `normalization.py` | reference curve and the interaction levels |
| `interaction_plot.py` | the sequence track and the arcs |
| `schema.py` | shape checks for the tables passed between stages |
| `search_chunk.py` | conditional analysis: split the ensemble on one interaction |
| `readpath.py` | sequence files |

## Tests

```
pytest
```

153 tests. Most run on synthetic data; those comparing against the real p53 output or a
CALVADOS mutant skip when the data folders are absent.

The suite pins the calibrated constants, the bin edges, the arc geometry and line widths,
and — when the reference data is present — that the rendered figure is **pixel identical**
to the published `p53.png`. If one of these fails, a published number has moved.

## Known behaviour

The coarse-grained reference curve returns values above 1 below separation 4 (13.12 at
s=1, 2.63 at s=2, 1.03 at s=3). Bonded neighbours are inside the 1.2 nm cutoff in every
frame, so they come out as strongly repulsive. `min_sep = 5` keeps them out of the figure,
but they remain in the interaction CSV. The count is exactly *N−1* for a chain of length
*N*, so it cancels when comparing equal-length sequences and does not cancel across
different lengths.

## Further reading

- `SUMMARY.md` — every change on this branch against the branches it was combined from
- `CHANGES_AND_TESTS.md` — verification status and what still needs testing
- `RUN_PATHS.md` — data folders used by previous runs
