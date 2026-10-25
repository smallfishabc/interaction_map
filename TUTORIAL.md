# Tutorial

A complete run on a small example that ships with the repository, so nothing has to be
downloaded or simulated first.

The example is `tests/fixtures/A39G_K24Y` — one CALVADOS coarse-grained simulation of a
p53 N-terminal double mutant. 61 residues, 1000 frames, 805 KB.

```
tests/fixtures/
  A39G_K24Y.pdb      topology, one bead per residue
  A39G_K24Y.dcd      1000 frames
  A39G_K24Y.fasta    the sequence
```

Everything below takes a few seconds. The same numbers are pinned in
`tests/test_pipeline_cg.py`, so if a step here gives you something different, the test
suite will say so too.

## Setup

```
conda env create -f environment.yml
conda activate interaction-map
```

On Windows the environment has to be **activated**, not just called by path, or the
compiled extensions cannot find their DLLs.

## 1. Draw the map

```
python main.py -p A39G_K24Y.pdb -x A39G_K24Y.dcd -dir tests/fixtures \
               -name A39G_K24Y --profile cg --out out --no_show
```

`--profile cg` selects the coarse-grained constants: a 1.2 nm contact cutoff, the
reference curve 13.12·s^−2.32, and interaction cutoffs (2, 1, −1, −2). Use
`--profile all_atom` for CAMPARI/ABSINTH trajectories; that is the default.

The sequence is found automatically. There is no `seq.txt` here, so `<name>.fasta` is
picked up and the run says which file it used:

```
sequence read from tests/fixtures\A39G_K24Y.fasta (61 residues)
```

Four files land in `out/`:

| File | Contents |
|---|---|
| `A39G_K24Y_1.2_contact_df_1201.csv` | contact probability per residue pair |
| `interaction_0317.csv` | the same table plus the reference curve, the interaction strength and the level |
| `A39G_K24Y.png` / `.svg` | the map |

The date stamps in those names are deliberate: they identify the version of the analysis
that produced the file.

## 2. Read the table

```
r_1  r_2  cont_prob  distance  r_1_aa  r_2_aa  gs_standard  relative_strength  plot_value
```

- `r_1`, `r_2` are **1-based sequence positions**, so `r_1 = 1` is the first residue.
- `r_1_aa`, `r_2_aa` are the residues themselves, read from the topology.
- `distance` is the sequence separation, j − i.
- `gs_standard` is what a generic flexible chain would do at that separation.
- `relative_strength` is `ln(cont_prob / gs_standard)`: positive means the pair touches
  more than a generic chain would, negative means less.
- `plot_value` bins that into 2, 1, 0, −1, −2 — the five arc styles.

For this mutant: 1830 pairs, 22 weakly attractive, 33 weakly repulsive, 60 strongly
repulsive, the rest neutral.

Those 60 are all bonded neighbours at separation 1. The coarse-grained reference curve
returns a value above 1 below separation 4, and adjacent beads are always within 1.2 nm,
so they come out strongly repulsive. They are kept out of the figure by `min_sep = 5`, but
they are still in the table — worth knowing if you ever count repulsive pairs.

## 3. Condition on an interaction

Split the ensemble on whether one interaction is formed and compare the two halves:

```
python main.py -p A39G_K24Y.pdb -x A39G_K24Y.dcd -dir tests/fixtures \
               -name A39G_K24Y --profile cg --out out --no_show --search 20 40
```

```
selecting frames with ChunkContact(20, 40, window=1, cutoff=1.2)
  46 of 1000 frames contain the interaction (4.6%)
```

A frame counts when **any** residue pair in the 3×3 window around (20, 40) is within the
cutoff. `--window 2` widens it, `--window 0` narrows it to the single pair.

This writes `ensemble_20_40.csv` with Rg, end-to-end distance, helicity, β content,
hydrogen bonds and the frame count for each side, plus a contact table for each.

Read the frame count first. 46 against 954 is a lopsided split, and the properties of the
46 carry a lot less weight than the properties of the 954.

**Helicity, β content and hydrogen bonds come out blank here.** DSSP needs a protein
backbone and the hydrogen bond search needs hydrogens, and a one-bead-per-residue model
has neither. They are left empty rather than reported as zero, because they are not
measurable rather than absent. Run the same command on an all-atom trajectory and all five
are filled in.

## 4. Condition on several interactions

Give more than one pair to compare the ensemble with and without a *group* of
interactions:

```
python main.py -p A39G_K24Y.pdb -x A39G_K24Y.dcd -dir tests/fixtures \
               -name A39G_K24Y --profile cg --out out --no_show \
               --search 20 40 30 50 --combine any
```

`--combine` decides how they are joined:

| | selects | this example |
|---|---|---|
| `any` | frames where **at least one** is formed | 78 of 1000 |
| `all` | frames where **every one** is formed together | 2 of 1000 |

The counts are a union and an intersection, and behave like it: (20,40) alone gives 46,
(30,50) alone gives 34, and 46 + 34 − 2 = 78.

`all` is the stricter question — what the chain looks like when these contacts are formed
at the same time — and will usually leave very few frames. Two frames out of a thousand is
not something to draw conclusions from; the run warns you if either side is empty, but it
cannot tell you when a side is merely too small.

Add `--with_maps` to also rebuild a contact map and figure for each side.

Output is `ensemble_any_20_40__30_50.csv`, named after the condition, with the mode and
the pairs recorded inside the file as well:

```
,interaction,non_interaction
Rg,2.0252417786838970,2.4283933295774083
Re,4.9985683774336790,6.0712451327778000
Frames,78.0,922.0
Heli,,
Beta,,
HB,,
mode,any,any
interactions,20-40; 30-50,20-40; 30-50
```

The chain is more compact when either contact is formed, 2.03 against 2.43, and more
compact still when both are formed at once, 1.67 — though that last figure rests on two
frames.

## 5. Several repeats

Give several trajectories of the same construct and they are treated as repeats:

```
python main.py -p top.pdb -x run1.xtc run2.xtc run3.xtc -name protein
```

The contact probability is unchanged — still the average over every frame — but the table
gains `n_rep` and `p_sem`, the standard error across repeats. With a single trajectory, as
in this tutorial, `n_rep` is 1 and `p_sem` is empty, because there is no spread to report.

This matters more than it sounds. On the all-atom p53 set, the strongest attractive pair
is 0.737 with a standard error of 0.186 across its five repeats.

## 6. The simulation database layout

For the `<name>-summary/<restype>/S_<psi>/` convention:

```
python main.py --database -dir "path/to/p53-summary" --psi 0 --restype BB --repeat 5
```

This resolves the folder convention and then runs exactly the same analysis.

## Where to go next

- `README.md` — the profiles, the residue numbering, and what each module does
- `scratch/README.md` — the checks that verify this code still reproduces the published
  output
- `SUMMARY.md` — what changed against the branches this was combined from
