# TODO

Blocks releasing alongside the paper. Nothing here blocks using the code.

- [ ] **Merge to `master`.** This work is all on `refactor_20261024`; `master` is still the
      2021 code.

- [ ] **`CITATION.cff`** at the root, with the paper's authors, title and DOI. GitHub turns
      it into a "Cite this repository" button.

- [ ] **Zenodo DOI.** Connect the repository to Zenodo *before* tagging a release, or the
      release is not captured. The DOI goes in the paper's data availability statement.

- [ ] **CI.** GitHub Actions running `pytest` on ubuntu. The 253 tests use the committed
      fixture and need none of the large data folders, so it works as-is. It is also the
      first evidence the code runs anywhere other than this machine.

## Known gap

The reference constants in `profiles.py` do not trace to a fit recorded here.
`util/calibrate.py` does the fit and `tutorials/calibration.ipynb` shows what re-deriving
them involves; it needs the GS linker trajectories.
