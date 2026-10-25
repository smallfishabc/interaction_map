# Change summary — `refactor_20261024`

Branch cut from `0811_final_all_atom_version` (04f2679).

One purpose: let the all-atom and the coarse-grained analysis run from a single code path
without either one's calibrated numbers being disturbed. Nothing was restructured beyond
what that required, and nothing was merged in from `1115_new_start`, `chunkreal` or
`clean_up_structure`.

Both paths are verified to reproduce their origin branch exactly. Details at the end.

```
 __init__.py               |   1 +
 contact_map_generation.py |  84 +++++++++++++++------
 default_function.py       |  25 +++++---
 interaction_plot.py       | 154 ++++++++++++++++++++++++---------------
 main.py                   |  90 ++++++++++++-----------
 normalization.py          |  43 +++++++++--
 profiles.py               | new
```

---

## 1. The two profiles

All-atom and CG were the same code with nine forked constants on two branches. Those
constants now sit together in `profiles.py` instead of being whichever branch you had
checked out.

| Constant | `all_atom` (from `0811`) | `cg` (from `0519_CG`) |
|---|---|---|
| `cutoff` | 0.8 nm | 1.2 nm |
| `pair_selection` | `mdtraj_ca` — `compute_contacts(contacts='all', scheme='CA', ignore_nonprotein=True)`, which drops \|i−j\| < 3 | `all_pairs` — `top.select_pairs('all','all')`, which keeps every pair |
| `slice_protein` | True | False |
| `residue_index_offset` | 0 | 1 |
| `reference_a`, `reference_b` | 1.64, −1.33 | 13.12, −2.32 |
| `inter_cutoff` | (1.5, 0.5, −1, −2) | (2, 1, −1, −2) |
| `value_list` | (2, 1, 0, −1, −2) | same |
| `min_sep` | 3 | 5 |
| `line_width_scale` | 2.0 | 4.0 |
| `prob_floor` | 0.0 | 0.0 |
| `residue_style` | `marker` | `marker` |

### On `residue_index_offset`

This is the one field whose name carries an argument, so it is worth stating plainly.

`r_1` / `r_2` are **1-based sequence positions in both profiles**, which is what the plot
layer and every existing CSV expect. They get there differently:

- All-atom: the CAMPARI topology carries an **ACE cap at residue index 0** and an NME cap
  at the end. `ignore_nonprotein=True` removes the caps from the pair list but not from
  the numbering, so the first real residue is already index 1. Offset 0.
- CG: CALVADOS topologies have no caps, so the first residue is index 0. Offset 1. The
  `0519_CG` branch did this with a hardcoded `pairs[:, 0] + 1`.

Confirmed against the data: the p53 protein slice has **63 residues** (ACE + 61 + NME) and
`compute_contacts` returns pair indices spanning exactly **1…61**.

The consequence is that all-atom correctness depends on the cap being present. An
uncapped all-atom topology would come out one residue low with no error, which is why
`check_residue_numbering()` exists.

### `min_sep` and `line_width_scale` are not new behaviour

Both replace hardcoded values and are numerically identical to what each branch did:

- all-atom `distance >= 3` replaces `if distance <= 2: continue`, which could never fire
  because mdtraj had already dropped those pairs.
- CG `distance >= 5` replaces `(interaction['distance'] > 4)`.
- all-atom 2.0 and CG 4.0 are each branch's own multiplier.

---

## 2. `profiles.py` — new

A frozen dataclass `Profile` with the twelve fields above, the two instances `ALL_ATOM`
and `CG`, a `PROFILES` lookup, and `get_profile()` which accepts a name, a `Profile`, or
`None` (meaning `ALL_ATOM`) and raises `ValueError` on anything else.

---

## 3. `contact_map_generation.py`

| Change | Before | After |
|---|---|---|
| `ContactProbData.__init__` | `(name, cutoff, traj, read_from_file=0)` | `+ profile=None`, stored as `self.profile` |
| `compute_contact` pair selection | one hardcoded `compute_contacts` call | branches on `profile.pair_selection`; unknown value raises |
| residue numbering | `pairs[:, 0]` | `pairs[:, 0] + profile.residue_index_offset` |
| debug output | `print(cont_df.dtypes)` on every run | removed |
| `load_traj_protein` | always `atom_slice(select('protein'))` | slices only when `profile.slice_protein`; gains `profile=` |
| `generate_contact` cutoff | `cutoff=0.8` | `cutoff=None`, taken from the profile |
| `generate_contact` read path | `ContactProbData(protein_name, cutoff, read_from_file)` | `ContactProbData(protein_name, cutoff, traj=None, read_from_file=1, profile=profile)` |

**New function `check_residue_numbering(contact_df, sequence)`.** Prints a warning when
the numbering does not span `1..len(sequence)` and returns a bool. Diagnostic only; it
alters nothing. It guards the uncapped-topology case described above.

**Bug fixed.** The `read_from_file=True` branch passed the flag into the `traj` argument
slot, leaving `read_from_file` at its default of 0, so the object tried to compute
contacts from the boolean `True`. Reading a saved contact map had never worked on any
branch. It now does. This is new behaviour rather than a restored behaviour, so it has no
regression baseline.

`load_traj_protein` returning a three-tuple matches `0811`. The `0519_CG` branch returned
only the trajectory; any CG caller expecting one value needs updating, and
`generate_contact` already was.

---

## 4. `normalization.py`

| Change | Before | After |
|---|---|---|
| signature | `(target_map, a1=1.64, b1=-1.33, inter_cutoff=(1.5,0.5,-1,-2), value_list=(2,1,0,-1,-2))` | `(target_map, profile=None, a1=None, b1=None, inter_cutoff=None, value_list=None, prob_floor=None)` |
| constants | baked into the defaults | taken from the profile; every one still overridable by keyword |
| debug output | `print(target_map.columns)` | removed |
| `prob_floor` | did not exist | optional; at the default 0.0 it is a no-op |

`gs_standard` is still computed with the original row-wise
`target_map.apply(fitting_function, axis=1, args=(a1, b1))`. A vectorised form was tried
and reverted: it gives bit-identical results, so there was nothing to gain from changing
the form of a calibrated step.

**New helper `reference_probability(separation, a, b)`** evaluates the same curve on a
plain number or array rather than a dataframe row. It is not used by `normalization`; it
exists for the reference-curve fitting work.

**`prob_floor`** forces `plot_value` to 0 where `cont_prob <= prob_floor`. The published
Methods state a cutoff was applied to p₁ as well as to I, and the old `chunkreal` code did
this as `i > 0.001`, but neither current branch has it. Left at 0.0, which reproduces both
branches, until the value used is confirmed.

---

## 5. `interaction_plot.py`

**Three new module-level tables.** `RESIDUE_TYPE` and `RESIDUE_COLOR` are lifted verbatim
from `0519_CG`. `INTERACTION_STYLE` maps level → (colour, curvature) using exactly the
values both branches used; it replaces a chain of four separate `if` statements that left
`colorset` and `connect` undefined for any level outside {2, 1, −1, −2}.

| Function | Change |
|---|---|
| `color_text` | gains `residue_style='marker'`. `'letter'` is the `0811` renderer unchanged; `'marker'` is the `0519_CG` renderer unchanged. An unrecognised residue now raises `ValueError` naming the residue and position, instead of `KeyError` |
| `create_color_coding` | gains `residue_style` and threads it to `color_text` |
| `interaction_plotting` | gains `profile`; colour/curvature from the table; `min_sep` filter; `line_width_scale` from the profile; `r_1`/`r_2` cast with `int()`; per-row `print(data)` removed; the unused `full_strength = []` assignments removed |
| `interaction_map` | gains `profile`, `show=True`, `save=True`; returns `(fig, ax)`; `graphg.add_edge(1, 5)` removed |

`graphg.add_edge(1, 5)` ran after `nx.draw` had already rendered, so it never appeared in
any figure. Its removal is confirmed by the pixel-identical comparison.

networkx was still used at this point. It was removed in the cleanup phase below, and the
figure stayed pixel identical.

---

## 6. `default_function.py`

`interaction_map_pairwise` gains `profile=None`, `show=True` and
`interaction_csv='interaction_0317.csv'`, and threads the profile into
`generate_contact`, `normalization` and `interaction_map`. The status line now prints the
profile name.

The output filename is unchanged. An earlier version of this branch renamed it per profile
to avoid an all-atom and a CG run colliding in one folder; that was reverted to stay
faithful. The collision is real but pre-existing — run the two profiles in one directory
and the second overwrites the first, unless `interaction_csv=` is passed.

---

## 7. `main.py`

Rewritten, but the default behaviour is unchanged: with no flags it runs the
multi-trajectory path, exactly as `0811` did.

| Change | Before | After |
|---|---|---|
| hardcoded paths | ~20 commented `path = ...` lines with one active | one `DEFAULT_PATH` constant; the full list moved to `RUN_PATHS.md` |
| `--profile` | did not exist | `all_atom` (default) or `cg` |
| `--no_show` | did not exist | suppresses the figure window; png and svg still written |
| `--single_traj` | declared with `default=1`, then overwritten by a hardcoded `single_traj = 0`, so the flag did nothing | a real `store_true` flag |
| single-trajectory code | inline in `__main__` | `single_trajectory(args)` |
| protein name | `path.split("\\")[-1].split("-")[0]` inside a bare `try/except` printing `'err'` | `os.path.basename(os.path.normpath(path)).split('-')[0]`, no bare except |
| `test = 1` dead branch | present | removed |
| `--repeat` | passed through as a string | cast with `int()` |

Note for CG users: `0519_CG` left the `single_traj = 0` line commented out, so that branch
defaulted to the **single-trajectory** path. The combined script defaults to
multi-trajectory for both profiles; pass `--single_traj` for the CG style of input.

---

## 8. `__init__.py`

`import profiles` added. The pre-existing `import main` is untouched.

---

## 9. Changed, then reverted

Recorded because the reasoning is the kind that looks right and is not.

**An "off-by-one fix" in the all-atom plot path.** The argument was that
`compute_contacts` returns 0-based indices while node *k* draws `seq[k-1]`, so
`layout[r_1]` must be one residue off. The data says otherwise: `interaction_0317.csv`
has `r_1` starting at 1 and `r_2` reaching 61 for a 61-residue sequence, because of the
ACE cap. `layout[r_1]` was correct, and the change would have displaced every all-atom arc
by one residue. Reverted.

**An `n_frames` column on the contact dataframe.** Useful for recovering the detection
floor, but it changes the saved CSV format and nothing asked for it. Reverted.

**A per-profile interaction CSV filename.** See section 6. Reverted.

**A vectorised `gs_standard`.** Bit-identical, therefore pointless. Reverted.

---

## 10. Verification

### All-atom, against `All_atom_p53_interaction` (August 2023)

| Check | Result |
|---|---|
| rendered figure vs `p53.png` | **pixel identical, 0 of 1,000,000 pixels differ** |
| `plot_value`, all 1711 pairs | identical |
| level counts (2/1/0/−1/−2) | 60 / 129 / 625 / 166 / 731, exact |
| arc counts and endpoints, all four levels | exact |
| contact map recomputed from the 5 repeats — `r_1`, `r_2`, `distance` | exact |
| `cont_prob`, `gs_standard`, `relative_strength`, line widths | within 1e-11 relative |

The pixel-identical figure used `residue_style='letter'`.

The residual ~1e-14 is not attributable to this branch. The stored reference cannot
regenerate its own `relative_strength` from its own `gs_standard` to better than 5.3e-13,
and row-wise `apply`, `Series **`, `array **` and a Python loop all reproduce the file
equally well — 9.76% of rows bit-identical in every case. It is the 2023 numpy/libm
against the current one.

### CG, against the `0519_CG` modules taken from git

Verified on two mutants from `CG_p53_random`, `A39G_K24Y` and `D21E_W53D`:

| Check | Result |
|---|---|
| contact map, 4 columns, 1830 pairs | **bit-identical** |
| interaction table, 7 columns incl. `plot_value` | **bit-identical** |
| rendered figure | **pixel identical, 0 pixels differ** |

### Harnesses

- `regress.py` — all-atom. Fast mode replays normalization, arc selection and the figure
  from the stored contact probabilities. `--full` also recomputes from the 5 repeats, on a
  copy. Reads the data folder, never writes to it.
- `cg_faithful.py [MUTANT]` — CG. Checks the `0519_CG` modules out of git into a temp
  folder, runs both versions on the same mutant, compares everything including the PNG.

Both live in the session scratch folder. Run them after every refactoring step.

---

## 11. Known and accepted

**The CG reference curve returns probabilities above 1 below separation 4** — 13.12 at
s=1, 2.63 at s=2, 1.03 at s=3. Bonded neighbours sit inside the 1.2 nm cutoff in every
frame, so `ln(1.0/13.12) = −2.57` puts all 60 of them at level −2. The all-atom path never
sees this because `contacts='all'` drops \|i−j\| < 3; the CG path keeps every pair.
`min_sep=5` removes them from the figure, so figures are unaffected, but the interaction
CSV carries them. The count is exactly *N−1* for a chain of length *N*, so it cancels when
comparing equal-length mutants and does not cancel across different lengths. Left as is,
to stay faithful to `0519_CG`.

---

## 12. Supporting files added

| File | Contents |
|---|---|
| `profiles.py` | the two profiles |
| `SUMMARY.md` | this document |
| `CHANGES_AND_TESTS.md` | verification status, how to run the environment, what still needs testing |
| `RUN_PATHS.md` | every hardcoded data path from every branch, grouped by purpose, with the folder layout |
| `.gitignore` | environments, `__pycache__`, and the three data folders, which total 1.2 GB |

An environment lives in `.conda-env/` (ignored): Python 3.10, mdtraj 1.9.9, numpy 1.23.5,
pandas 1.5.3, scipy 1.9.3, matplotlib 3.6.3, networkx 2.8.8, setuptools 68.2.2. On Windows
its library directories must be on `PATH` or the compiled extensions cannot resolve their
DLLs and the interpreter exits with no traceback. See `CHANGES_AND_TESTS.md`.

---

## 13. Phases after the profile merge

Sections 1 to 12 describe the first step, combining the two profiles. Three further phases
followed on the same branch.

**Cleanup.** networkx removed, with the figure verified pixel identical before and after,
so the package no longer depends on it. `readpath.subdir` removed, file handles closed.
`schema.py` added, checking the shape of the tables passed between stages and wired into
the contact csv read as a warning. `check_sequence_against_topology` added, which warns
when the sequence file does not match the simulated topology. `environment.yml` and
`pyproject.toml` added.

**Data and interface.** The contact table gained `r_1_aa`, `r_2_aa`, `r_1_chain`,
`r_2_chain` from the topology, and `n_frames`, `n_rep`, `p_sem` describing the spread
between repeats, with `cont_prob` itself untouched. `classify()` split out of
`normalization` so the five levels can be recut without recomputing contacts. Residue
indices pinned to int64, which a test caught drifting between int32 and int64 across a csv
round trip. `--out` added for the results directory; the date stamps in the file names stay,
since they identify the version of the analysis. The entry point now treats a single
topology and trajectory as the primary path, with `--database` as the adapter for the
`<name>-summary` layout, and `-x` accepting several trajectories as repeats.

**Conditional analysis.** `search_chunk.py` ported from `1115_new_start`, reproducing it
bit for bit on both contact maps and all five ensemble observables. The 3x3 frame selector
is now a `ChunkContact` object, so the window, the cutoff and the residue numbering are
explicit; the numbering matters because the two profiles centre the window differently. The
ensemble comparison now records how many frames fell on each side of the split. Reachable
as `--search R1 R2`.

Still on other branches, untouched: chunk aggregation and chunk plotting (`chunkreal`), the
mutation-target sweep (`clean_up_structure:pick_ensemble.py`), `interaction_strength.py`,
which cannot run on any branch, and the reference-curve fit, which exists on no branch.

With no code anywhere in the repository: PCA clustering, hydrophobic scaling, the
solution-sensitivity regression, and the weighted 3-residue chunk of Methods Eq. 3. Their
outputs are in `All_atom_p53_interaction`.

Open decisions: the p₁ noise cutoff value; which of the four threshold sets in circulation
are the published ones; whether `residue_style` should default to `letter` so published
figures stay reproducible by default; and how a CG sequence should be located, since
`readsequence_single` reads `seq.txt` while `CG_p53_random` carries `<name>.fasta`.
