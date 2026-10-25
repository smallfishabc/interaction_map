# scratch

Throwaway checks, kept because they are the only thing standing between a future edit and
a silently changed published number.

They are not part of `pytest`: each one needs one of the data folders, which are too large
to commit, and two of them check old branches out of git and run them. Run them by hand
after any change that touches the analysis.

Working copies and fixtures are written to `scratch/_work/`, which is not tracked. The data
folders themselves are only ever read.

| Script | Compares against | Needs |
|---|---|---|
| `regress.py` | the August 2023 p53 output | `All_atom_p53_interaction/` |
| `cg_faithful.py` | the `0519_CG` branch, taken from git | `CG_p53_random/` |
| `cond_faithful.py` | the `1115_new_start` branch, taken from git | `All_atom_p53_interaction/` |

```
python scratch/regress.py                 replay normalization, arc selection, the figure
python scratch/regress.py --full          also recompute the contact map from the 5 repeats
python scratch/cg_faithful.py [MUTANT]    default A39G_K24Y
python scratch/cond_faithful.py [R1 R2]   default 13 24
```

Exit code 0 means everything matched.

## What each one establishes

`regress.py` renders the map from the stored interaction values and compares it against
`p53.png` pixel by pixel, then checks `plot_value` on all 1711 pairs and the four
continuous columns. `--full` recomputes the contact map from all five repeats.

`cg_faithful.py` and `cond_faithful.py` have no stored output to compare against, so they
generate one: the old branch's modules are checked out into a temporary folder, run on the
same input, and every column of every table compared.

Both of them align two things first, so that the question being asked is whether the port
is faithful rather than whether the branches differ:

- `1115_new_start` carries the older reference exponent `b = -1.32`; the all-atom profile
  uses `-1.33`, which is what `0811` carries and what the published output was made with.
- `1115_new_start` was written against networkx below 2.4, where `Graph.node` was the
  public accessor. The project is pinned to a later networkx, so the accessor is spelled
  the current way in the copy. Reading a node attribute computes nothing.

`cond_faithful.py` reports the figure difference rather than failing on it: the arc width
formula changed between `1115` and `0811`, and the `0811` form is the published one, which
`regress.py` verifies pixel for pixel.

## Tolerances

`regress.py` compares the continuous columns with a relative tolerance of 1e-11 rather than
demanding equality. The stored reference was produced on a different numpy build: it cannot
regenerate its own `relative_strength` from its own `gs_standard` to better than 5.3e-13,
and row-wise `apply`, `Series **`, `array **` and a plain Python loop all reproduce it
equally well. What decides the figure — the interaction level, the arc count, the arc
endpoints — is compared exactly.
