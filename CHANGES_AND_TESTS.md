# Step 1 — two profiles in one branch

Branch: `refactor_20261024`, cut from `0811_final_all_atom_version` (04f2679).

Scope of this step only: the all-atom and coarse-grained analyses now run from one code
path, selected with `--profile`. Nothing was merged in from `1115_new_start`,
`chunkreal` or `clean_up_structure` yet. No restructuring beyond what was needed to let
the two parameter sets coexist.

## Why a profile object

All-atom and CG were the same code with nine forked constants. Picking either set as
"the" default would silently change the other analysis. They are now declared side by
side in `profiles.py`:

| | `all_atom` | `cg` |
|---|---|---|
| source branch | `0811_final_all_atom_version` | `0519_CG` |
| contact cutoff | 0.8 nm | 1.2 nm |
| pair selection | `compute_contacts(contacts='all', scheme='CA')`, drops \|i−j\| < 3 | `top.select_pairs('all','all')`, keeps every pair |
| slice to `protein` | yes | no |
| residue index written to `r_1`/`r_2` | 0-based | 1-based |
| reference curve `a`, `b` | 1.64, −1.33 | 13.12, −2.32 |
| interaction cutoffs | (1.5, 0.5, −1, −2) | (2, 1, −1, −2) |
| minimum separation drawn | 3 | 5 |
| line width scale | 2× | 4× |

The steeper CG exponent matches Fig. 8A of the manuscript, so both sets are real.

## What changed, file by file

**`profiles.py`** (new) — the table above, plus `get_profile()`.

**`contact_map_generation.py`**
- `ContactProbData` and `generate_contact` take a `profile`; cutoff defaults from it.
- Pair selection switches on `profile.pair_selection`.
- `load_traj_protein` skips the `protein` atom slice for CG and still returns the
  `(full_traj, sliced_traj, n_frames)` triple. The CG branch returned only the trajectory.
- The contact dataframe gains an **`n_frames`** column, so the detection floor
  (1/n_frames) is recoverable when the CSV is read back.
- Fixed: `generate_contact(..., read_from_file=True)` passed the flag into the `traj`
  argument, so reading a saved contact map had never worked.
- Removed `print(cont_df.dtypes)`.

**`normalization.py`**
- Takes a `profile`; `a1`, `b1`, `inter_cutoff`, `value_list` still override it.
- `reference_probability(separation, a, b)` added so the curve can be evaluated on a
  plain number or array, not only on a dataframe row. Needed by the fitting work later.
- `gs_standard` is now computed vectorised instead of row-wise `apply`. Same arithmetic.
- Optional `prob_floor`: contact probabilities at or below it are forced to level 0.
  Default 0.0, which reproduces the current branches. The published Methods apply a
  cutoff to p₁ as well as to I; the old `chunkreal` code did this as `i > 0.001`. Left
  off until you confirm the value.
- Removed `print(target_map.columns)`.

**`interaction_plot.py`**
- Marker style (coloured dots, 5 residue classes) is now the default for both profiles.
  The letter style is still there, selected by `profile.residue_style = 'letter'`.
- Arc colour and curvature come from one table; previously `colorset` and `connect` were
  left undefined for any level outside {2, 1, −1, −2}.
- **Off-by-one fix on the all-atom path — see below.**
- `min_sep` and the line width scale come from the profile.
- Removed the per-row `print(data)` and the `graphg.add_edge(1, 5)` that ran after
  drawing and had no effect.
- `interaction_map` returns `(fig, ax)` and takes `show=` and `save=`.

**`default_function.py`** — threads the profile through. The output CSV is now
`interaction_all_atom.csv` / `interaction_cg.csv` instead of `interaction_0317.csv`, so
the two profiles do not overwrite each other in the same folder.

**`main.py`**
- `--profile {all_atom,cg}`, default `all_atom`.
- `--single_traj` is a real flag now. It used to be overridden by a hardcoded
  `single_traj = 0` on the all-atom branch, and left at `1` on the CG branch, so the two
  branches silently defaulted to opposite modes.
- `--no_show` for batch runs.
- One `DEFAULT_PATH` constant at the top; the stack of commented paths moved to
  `RUN_PATHS.md`.
- Protein name derived with `os.path.basename`, which also works on Linux paths, and the
  bare `except` around it is gone.

## Residue numbering — a change I made and then reverted

I briefly "fixed" an off-by-one in the all-atom plotting path. **That was wrong and has
been reverted.** Recording it because the reasoning matters for anyone who looks at this
code later and reaches for the same fix.

The argument was: `compute_contacts` returns 0-based residue indices, network node *k*
draws `seq[k-1]`, so looking up `layout[r_1]` directly must be one residue off.

The real data says otherwise. `All_atom_p53_interaction/.../interaction_0317.csv` has
`r_1` starting at **1** and `r_2` reaching **61** for a 61-residue sequence, i.e. already
1-based. The reason is in the topology: the CAMPARI PDB carries an **ACE cap as residue
index 0** and an NME cap at the end. `ignore_nonprotein=True` removes the caps from the
pair list but not from the numbering, so the first real residue is already index 1.
`layout[r_1]` was correct all along, and my change would have displaced every all-atom arc
by one residue.

`profiles.residue_index_offset` now records this explicitly: 0 for all-atom because the
cap already provides the shift, 1 for CG because CALVADOS topologies have no caps. Both
profiles end up writing 1-based sequence positions, which is what the plot layer and every
existing CSV expect.

The fragile part is that all-atom correctness depends on the cap being present. An
uncapped all-atom topology would come out one residue low with no error.
`contact_map_generation.check_residue_numbering()` prints a warning when the numbering
does not span `1..len(sequence)`. It is diagnostic only and changes nothing.

## Changes that can move numbers or figures

After the revert, the only deliberate visual change is the residue style:

1. **Residue markers instead of letters.** Coloured dots with 5 residue classes, the CG
   branch's renderer, now the default for both profiles. Set
   `profiles.ALL_ATOM.residue_style = 'letter'` to get the old letters back. The paper
   figures use letters.

2. Reading a saved contact map now works rather than raising. This path could never have
   run before, so it is new behaviour, not a regression.

Everything else was checked to be numerically identical:

- the vectorised reference curve against the old row-wise `apply`
- `distance >= 3` (all-atom) in place of a `distance <= 2: continue` that could never fire
- `distance >= 5` (CG) in place of `distance > 4`
- arc colours and curvatures, line width scales (2× all-atom, 4× CG)
- the contact CSV and interaction CSV keep their original names and columns

The two separation filters were rewritten but are numerically identical: all-atom
`distance >= 3` replaces a `distance <= 2: continue` that could never fire, and CG
`distance >= 5` replaces `distance > 4`.

## Verification status

Everything below has been run against `All_atom_p53_interaction`, in the pinned
environment described further down.

| Check | Result |
|---|---|
| `plot_value`, all 1711 pairs | **identical** |
| interaction level counts (2/1/0/−1/−2) | 60 / 129 / 625 / 166 / 731, exact |
| arc counts and endpoints, all four levels | exact |
| rendered figure vs `p53.png` | **pixel identical, 0 of 1,000,000 pixels differ** |
| contact map recomputed from the 5 repeats | `r_1`, `r_2`, `distance` exact |
| `cont_prob`, `gs_standard`, `relative_strength`, line widths | within 1e-11 relative |

The residual ~1e-14 differences are not ours. The stored reference cannot regenerate its
own `relative_strength` from its own `gs_standard` column to better than 5.3e-13, and
row-wise `apply`, `Series **`, `array **` and a plain Python loop all reproduce the file
equally well (9.76% of rows bit-identical in every case). It is the 2023 numpy/libm
against this one.

The pixel-identical figure was produced with `residue_style='letter'`. That is the
renderer behind the published figures.

### CG profile

There is no stored CG output to regress against, so the check regenerates one: the
`0519_CG` modules are taken straight out of git, run on a CALVADOS mutant from
`CG_p53_random`, and compared against the combined code on the same input. Verified on
two independent mutants, `A39G_K24Y` and `D21E_W53D`:

| Check | Result |
|---|---|
| contact map, all 4 columns, 1830 pairs | **bit-identical** |
| interaction table, all 7 columns incl. `plot_value` | **bit-identical** |
| rendered figure | **pixel identical, 0 pixels differ** |

Known and accepted: the CG reference curve gives `p_ref > 1` below separation 4
(13.12 at s=1, 2.63 at s=2, 1.03 at s=3). Bonded neighbours sit inside the 1.2 nm cutoff
in every frame, so `ln(1.0/13.12) = -2.57` puts all 60 of them at level -2. The all-atom
path never sees this because mdtraj's `contacts='all'` drops `|i-j| < 3`, while the CG
path keeps every pair. `min_sep=5` removes them from the figure, so the figures are
unaffected; the interaction CSV still carries them. The count is exactly N-1 for a chain
of length N, so it cancels when comparing equal-length mutants and does not cancel across
different lengths. Left as it is, to stay faithful to `0519_CG`.

## Regression harness

`regress.py` lives in the session scratch folder. Run it after every refactoring step.

```
regress.py            fast: replays normalization, arc selection and the figure from the
                      stored contact probabilities. No trajectory needed.
regress.py --full     also recomputes the contact map from the 5 repeats.
```

It reads `All_atom_p53_interaction` and never writes to it; `--full` works on a copy of
the trajectory in scratch.

`cg_faithful.py [MUTANT]` is the CG counterpart. It checks the `0519_CG` modules out of
git into a temp folder, runs both versions on the same mutant and compares everything
including the rendered PNG. Run it after any change that touches the CG path.

## Running the environment

`.conda-env/` in the project folder: Python 3.10, mdtraj 1.9.9, numpy 1.23.5,
pandas 1.5.3, scipy 1.9.3, matplotlib 3.6.3, networkx 2.8.8, setuptools 68.2.2.

**Windows needs the env's library directories on `PATH`**, otherwise the compiled
extensions cannot resolve their DLLs and the interpreter dies with no traceback. Calling
`.conda-env/python.exe` directly is not enough. Either activate the env, or:

```
E=.conda-env
PATH="$E:$E/Library/bin:$E/Library/mingw-w64/bin:$E/Library/usr/bin:$E/Scripts:$PATH"
```

Two notes on versions. mdtraj publishes no Windows wheels for any 1.9.x, so pip cannot
install it and this has to be a conda env. mdtraj 1.9.6 has no Python 3.10 build at all;
1.9.9 is from April 2023 and predates the August 2023 reference, and is also what the
anaconda base env carries. setuptools must stay below 81, since mdtraj 1.9.x imports
`pkg_resources`.

## What I would like tested

Still untested by me, needs your data or your judgement:

**1. The CG profile end to end.** `CG_p53_random/` holds ~400 CALVADOS mutants as
`.dcd` / `.pdb` / `.fasta` triples. The CG path has not been exercised since the merge,
and there is no stored CG output to compare against, so this needs your eye on the result
rather than a regression.

```
python main.py --single_traj --profile cg -dir CG_p53_random -p A39G_K24Y.pdb -x A39G_K24Y.dcd -name A39G_K24Y --no_show
```

One thing will fail here: the single-trajectory path reads the sequence from `seq.txt`
via `readpath.readsequence_single`, and these folders carry `.fasta` files named after
each mutant instead. I have not changed that function, because deciding how a CG sequence
is located is a real choice, not a cleanup. Options: read `<name>.fasta`, accept a
`--seq` argument, or take the sequence from the topology.

**2. The 2024 SARS-CoV-2 linker set.** `2024_simulation/SARS_COV_2_LINKER-summary` has
the full `BB/S_-3 … S_3` condition spread and no prior analysis output, so it is a clean
check that the multi-condition path still works.

```
python main.py -dir "2024_simulation/SARS_COV_2_LINKER-summary" --psi 0 --no_show
```

**3. Read back a saved map.** `read_from_file=True` now loads the contact CSV instead of
raising. That path had never run before, so it is new behaviour rather than a regression
check, and deserves a look at whether the loaded frame is what you expect.

**4. Residue style for the paper.** The marker renderer is the default, as you asked. The
letter renderer is what reproduces `p53.png` pixel for pixel. If the published figures
should stay reproducible by default, this is a one-word change in `profiles.py`.
