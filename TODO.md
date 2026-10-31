# TODO

Work deferred, roughly in the order it would need doing. Nothing here blocks using the
code; all of it blocks publishing it properly.

## Before release

**Citation.** Add `CITATION.cff` at the root, with the paper's authors, title and DOI once
it has one. GitHub reads this file and shows a "Cite this repository" button, and reference
managers pick it up. The README currently says a citation file will be added.

**Archive and DOI.** Connect the repository to Zenodo *before* tagging a release. Zenodo
then mints a DOI for the tagged snapshot, and that DOI is what goes in the paper's data
availability statement. Reviewers look for this more than they look for a PyPI entry. The
order matters: connect first, then tag, or the release is not captured.

**Continuous integration.** A GitHub Actions workflow running `pytest` on ubuntu across
Python 3.10 to 3.12. The 253 tests run on the committed coarse-grained fixture and need
none of the large data folders, so this works without any further setup. Two things it
buys: a visible green badge, and the first evidence that the code runs on anything other
than this one Windows machine.

**Merge to master.** All of this work is on `refactor_20261024`; `master` is still the 2021
code. The Colab notebook installs the package from the branch, so that line needs updating
to a tag once merged.

**Repository name.** Rename `inteactionmap` to `interaction_map` in GitHub settings, then
locally:

```
git remote set-url origin https://github.com/smallfishabc/interaction_map.git
```

GitHub redirects the old name, so nothing breaks before this is done.

## Documentation

**Tutorial.** `TUTORIAL.md` exists and every command in it has been run, but it was written
before the mutation scanning and the calibration util landed, so it covers the map and the
conditional analysis only. It needs a pass to cover the full path from simulation to
mutant list, and to agree with the Colab notebook rather than duplicate it.

**ColabMap.** `Google Colab/InteractionMap_2026.ipynb` has been run only for its analysis
half, against the committed fixture. The simulation half has never been executed in this
project. Run it once end to end on a Colab GPU before linking it anywhere, then add the
badge to the README.

## Analyses still outside the package

Described in the paper, currently only in notebooks under `temp_examples/`:

- PCA clustering of conformations by interaction pattern (Fig. 3D, E)
- Hydrophobic scaling per interaction chunk (Fig. 5C/D, 6C/D) — partly covered, since
  `util/mutation_scan.py` computes per-chunk Kyte-Doolittle, charge and aromatic counts
- Solution-sensitivity regression against attractive interaction count (Fig. 4D, E), which
  depends on an external `simulationdatabase` package not in this repository

## Open scientific questions

**The reference constants have no recorded provenance.** `profiles.py` carries
a = 1.64, b = −1.33 for all-atom and a = 13.12, b = −2.32 for coarse-grained. Neither
traces to a fit kept here. `util/calibrate.py` performs the fit and
`tutorials/calibration.ipynb` shows what re-deriving them involves; it needs the GS linker
trajectories, 8 to 32 repeats, as the Methods describe.

**Which cutoffs are the published ones.** The contact-probability cutoff differs between
the sources: 0.1 with a minimum separation of 10 in the 2024 scan notebook, 0.01 with a
minimum separation of 5 in the 2025 Colab notebook, and none at all in the shipped
analysis. They are user-facing parameters now, but the paper reports one pair.

**Which interaction thresholds are the published ones.** Four sets are in circulation
across the old branches. The current pairwise set reproduces the published p53 figure, so
that one is settled; whether the chunk set is a separate calibrated choice is not.
