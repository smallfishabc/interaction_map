# Ensemble Map

Code for **"Ensemble Map Enabled Intrinsically Disordered Protein Design by
Understanding Intramolecular Interactions."**

Feng Yu, Shahar Sukenik

---

An intrinsically disordered region has no fixed structure, but it is not a random coil
either. Pairs of residues touch more often, or less often, than a featureless chain of the
same length would, and those preferences shape what the region does.

A raw contact map cannot show this. Long-range contacts have low absolute probability
whatever their meaning, so a map of raw probabilities buries the interactions that matter.
The Ensemble Map divides that probability by what an interaction-free chain would give
at the same sequence separation, which puts short-range and long-range contacts on one
scale.

## The method

For a pair of residues separated by `s = j - i` along the sequence:

**1. Contact probability.** The fraction of simulation frames in which the two residues
are within a cutoff, measured Cα–Cα for all-atom.

**2. Reference.** What a glycine–serine linker gives at the same separation, a power law
fitted on GS linkers of 8 to 32 repeats:

```
p_ref(s) = a * s ** b                                  Eq. 1
```

**3. Interaction strength.** How far the pair departs from that reference:

```
I = ln( p_obs / p_ref )                                Eq. 2
```

Positive means the pair touches more than an interaction-free chain would — attraction,
from opposite charges, hydrophobic contact or π–π stacking. Negative means less —
repulsion, from like charges or from interactions drawn elsewhere in the sequence.

**4. The map.** `I` is binned into five levels and drawn as arcs along the sequence:
attractive above, repulsive below, residues coloured by type, arc thickness following the
contact probability and the strength.

Cutoffs are applied to both `I` and `p_obs` to suppress noise, and held constant across
every map.

## What is here

| | |
|---|---|
| `src/ensemble_map/` | the analysis |
| `src/ensemble_map/util/` | tools beside it, currently the reference-curve fit |
| `scripts/main.py` | run without installing |
| `tutorials/` | worked notebooks |
| `tests/` | 213 tests, with a small coarse-grained example |
| `eval/` | checks that the code still reproduces the published output |

## Install

```
conda env create -f environment.yml
conda activate ensemble-map
pip install -e . --no-deps
```

mdtraj publishes no Windows wheels for 1.9.x, so this has to be a conda environment rather
than a pip virtualenv, and `--no-deps` keeps pip from trying to supply them. On Windows the
environment must be **activated**, not called by path, or the compiled extensions cannot
find their DLLs.

Installing gives the `ensemble-map` command. Without installing,
`python scripts/main.py` takes the same arguments.

## Try it without any data

One CALVADOS simulation of a p53 N-terminal double mutant ships with the repository:
61 residues, 1000 frames, 805 KB.

```
ensemble-map -p A39G_K24Y.pdb -x A39G_K24Y.dcd -dir tests/fixtures \
                -name A39G_K24Y --profile cg --out out
```

`TUTORIAL.md` walks through this end to end, including conditioning the ensemble on chosen
interactions. Its numbers are pinned by the test suite, so the two cannot drift apart.

## Running your own

All-atom, five repeats in the simulation database layout
`<name>-summary/<restype>/S_<psi>/`:

```
ensemble-map --database -dir path/to/p53-summary --psi 0 --restype BB --repeat 5
```

A topology and one or more trajectories directly. Several trajectories of the same
construct are treated as repeats, and the spread between them is reported:

```
ensemble-map -p top.pdb -x run1.xtc run2.xtc run3.xtc -name protein
```

Outputs land beside the trajectory unless `--out` says otherwise:

| File | Contents |
|---|---|
| `<name>_<cutoff>_contact_df_1201.csv` | contact probability per pair, with residue identity and the spread across repeats |
| `interaction_0317.csv` | the above plus the reference, the interaction strength and the level |
| `<name>.png`, `.svg` | the map |

`r_1` and `r_2` are 1-based sequence positions.

## Conditioning on an interaction

Resolving an ensemble into subpopulations by which contacts are formed. Split the frames
on one interaction and compare the halves:

```
ensemble-map -p top.pdb -x traj.xtc -name protein --search 13 24
```

A frame counts when any pair in the 3×3 window around (13, 24) is within the cutoff.
`--window` changes the half width. Several interactions can be combined:

```
ensemble-map -p top.pdb -x traj.xtc -name protein --search 13 24 37 44 --combine all
```

`any` selects frames where at least one is formed, `all` where every one is formed
together. The output records Rg, end-to-end distance, helicity, β content, hydrogen bonds,
and **the frame count on each side** — a split where one side holds a handful of frames
means something different from an even one.

## Two force fields

The paper uses all-atom simulation for the maps and coarse-grained simulation to screen
mutant ensembles, because CALVADOS contact probability falls off about an order of
magnitude faster than CAMPARI and so loses sensitivity to mid- and long-range pairs
(Fig. 8A). Both are supported through `--profile`, and their constants are declared
together in `profiles.py` so that neither can be altered by editing the other:

| | `all_atom` | `cg` |
|---|---|---|
| force field | ABSINTH / CAMPARI v2_09052017 | CALVADOS2 (2023) |
| contact cutoff | 0.8 nm | 1.2 nm |
| distance measured | Cα–Cα | between beads |
| reference `a`, `b` | 1.64, −1.33 | 13.12, −2.32 |
| interaction cutoffs | (1.5, 0.5, −1, −2) | (2, 1, −1, −2) |
| minimum separation drawn | 3 | 5 |

### Fitting the reference curve

`util/calibrate.py` performs the GS linker fit of Eq. 1:

```
python -m ensemble_map.util.calibrate GS16/__START_0.pdb GS16/__traj_*.xtc \
    --cutoff 0.8 --sep 3 60 --out references/gs16.json
```

It records the cutoff, the scheme, the separation range, the trajectories, the standard
errors and the residuals alongside `a` and `b`. A curve fitted at one cutoff does not
describe another, and beyond its fitted range it is an extrapolation.

The constants in the table above were inherited from earlier work and do not trace to a
fit recorded in this repository. `tutorials/calibration.ipynb` shows what re-deriving them
involves, and what the GS16 notebooks give on their own.

## Simulation settings

From the paper's Methods, for reference when preparing input.

**All-atom** — ABSINTH implicit solvent, CAMPARI v2_09052017, 310 K, 10⁷ equilibration
steps, conformations written every 12,500 steps. Five independent simulations of ~5,600
conformations each per protein per condition, ~28,000 in total. PUMA scrambles used three
repeats, ~16,800 conformations.

**Coarse-grained** — CALVADOS2 (2023), 310 K, 1,000 frames per trajectory.

## Verification

This code was reorganised from several older branches. Three checks confirm the published
output is still reproduced, each run against the original code rather than against a
description of it:

| Check | Against |
|---|---|
| `eval/regress.py` | the August 2023 p53 output; the rendered figure matches `p53.png` pixel for pixel |
| `eval/cg_faithful.py` | the `0519_CG` branch, taken from git and run side by side |
| `eval/cond_faithful.py` | the `1115_new_start` branch, likewise |

They need the simulation data, so they are not part of `pytest`. `eval/README.md` explains
what each establishes.

## Not yet in this package

Analyses described in the paper that currently live in notebooks outside the repository:

- PCA clustering of conformations by interaction pattern (Fig. 3D, E)
- Mutation scanning with the 3-residue sliding window of Eq. 3 (Fig. 5, 6)
- Hydrophobic scaling per interaction chunk (Fig. 5C/D, 6C/D)
- Solution-sensitivity regression against attractive interaction count (Fig. 4D, E)

## ColabMap

A Colab notebook combining CALVADOS simulation with the Ensemble Map, for use without a
local installation: give it an IDR sequence and the residue types to mutate, and it returns
Interaction Maps and a list of mutants as CSV.

*Link to be added.*

## Further reading

- `TUTORIAL.md` — a complete run on the bundled example
- `tutorials/calibration.ipynb` — fitting the reference curve
- `RUN_PATHS.md` — the simulation folders behind the published figures
- `SUMMARY.md`, `CHANGES_AND_TESTS.md` — what changed in the reorganisation, and how it was checked

## Citing

Please cite the paper. A `CITATION.cff` will be added on publication.

## License

MIT. See `LICENSE`.

Work still outstanding before release is listed in `TODO.md`.
