# Data folders used by previous runs

`main.py` used to carry these as a stack of commented-out `path = ...` lines, with one
uncommented. They are collected here so the script keeps only a single `DEFAULT_PATH`
constant, and nothing is lost.

All of them follow the simulation database layout:

```
<name>-summary/
  seq.fasta              protein sequence, read by readpath.readsequence
  BB/                    residue group whose transfer free energy is altered (--restype)
    S_0/                 psi value of the transfer free energy (--psi)
      __START_0.pdb      topology
      __traj_0.xtc ...   one file per repeat (5 repeats, 3 for the PUMA scrambles)
```

Run one of them with:

```
python main.py -dir "F:\DATA_F\Interaction_map_target\DP00086_p53Nterm_P04637\p53-summary" --profile all_atom
```

## GS linkers — reference curve calibration

These are the ideal-polymer runs the `a`, `b` reference curve is fitted to.

| Path | Note |
|---|---|
| `F:\globus\simulation_contactmap_validation\GS18-summary` | |
| `F:\globus\simulation_contactmap_validation\GS22-summary` | |
| `F:\globus\simulation_contactmap_validation\GS44-summary` | |
| `F:\globus\simulation_contactmap_validation\GS56-summary` | |
| `F:\globus\simulation_contactmap_validation\ChargeP_30-summary` | charge patterning variant |
| `F:\globus\simulation_contactmap_validation` | holds `database_entry.csv`, read by `standard_fitting.directory_location` |
| `F:\globus\simulation_sticker_spacer\F1_GS_20-summary` | |
| `F:\globus\simulation_sticker_spacer\F1_GS_40-summary` | |
| `F:\DATA_F\GSlinker\GS56-summary` | |
| `F:\DATA_F\GSlinker\GS64-summary` | |

The manuscript Methods state the fit uses GS linkers of 8 to 32 repeats, i.e. 16 to 64
residues, which matches the GS18 - GS64 spread above.

## Target IDRs

| Path | Protein |
|---|---|
| `F:\DATA_F\Interaction_map_target\DP00086_p53Nterm_P04637\p53-summary` | p53 NTAD (Fig. 3) |
| `F:\DATA_F\Interaction_map_target\Individual Proteins\p53-summary` | p53, current `DEFAULT_PATH` |
| `F:\DATA_F\PDBsumreal_0814\p53_new-summary` | p53, earlier run |
| `F:\DATA_F\Interaction_map_target\DP00949_UPF2_Q9HAU5\DP00949-summary` | UPF2 MIF4G-2 NTD (Fig. 5) |
| `F:\DATA_F\Interaction_map_target\DP00719_TDG_Q13569\DP00719-summary` | TDG N-terminal IDR (Fig. 6) |
| `F:\DATA_F\Interaction_map_target\DP00007_APEX1_P27695\DP00007-summary` | APEX1 (Fig. 2C) |
| `F:\DATA_F\Interaction_map_target\DP00356_CITED2_Q99967\DP00356-summary` | CITED2 |
| `F:\DATA_F\Interaction_map_target\DP00549_zincfinger593_O00488\DP00549-summary` | zinc finger 593 |
| `F:\DATA_F\Interaction_map_target\DP00011_sulfotransferase_P0DMM9\DP00011-summary` | sulfotransferase |
| `F:\DATA_F\Interaction_map_target\DP01128_DP00018_CDKN1B_p27_P46527\DP00018-summary` | p27 / CDKN1B |
| `F:\DATA_F\Interaction_map_target\SAXS_Cov2_Linker_P0DTC9\SARS_COV_2_LINKER-summary` | SARS-CoV-2 linker |
| `F:\DATA_F\UGDH0.5_simplified_interaction_map` | UGDH |
| `F:\DATA_F\PDBsum\Ash1-summary` (Linux: `/media/lemoncatboy/WD_BLACK/DATA_F/PDBsum/Ash1-summary`) | Ash1 |

## PUMA scrambles

| Path | Note |
|---|---|
| `F:\DATA_F\puma_scramble_new\puma123\puma_wildfull-summary` | wild type |
| `F:\DATA_F\puma_scramble_new\puma_scrammble_sum\puma_wildfull-summary` | wild type, summary set |
| `F:\DATA_F\puma_scrammble_sum\puma_scramble_20-summary` | scramble 20 |
| `F:\DATA_F\puma_scramble_new\puma123\puma_scramble_2-summary` | scramble 2 |
| `F:\DATA_F\puma_scramble_new\puma123\puma_scramble_3-summary` | scramble 3 |
| `F:\DATA_F\Eli_lilly_figure\puma_wildfull-summary` | figure run |

PUMA scrambles were simulated with **3** repeats, not 5, so these need `--repeat 3`.

## Mutation and validation runs

| Path | Note |
|---|---|
| `F:\DATA_F\contact_mutation\p53_W53G-summary` | p53 point mutant |
| `F:\DATA_F\puma_scramble_new\Mutation\PUMAWT_13toS_23toS-summary` | PUMA double mutant |
| `F:\DATA_F\E1A_pat-summary` | E1A, the active path on the CG branch |
| `F:\DATA_F\All_analysis_test_github` | the folder `1115_new_start` pointed at for the conditional search |

## Coarse-grained / Colab

`CG_google_colab` ran from `D:\Test_colab` on a single trajectory pair rather than the
database layout:

```
D:\Test_colab\p53NTAD.pdb
D:\Test_colab\p53NTAD.dcd
sequence: MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPD
```

Equivalent command with the combined script:

```
python main.py --single_traj --profile cg -dir "D:\Test_colab" -p p53NTAD.pdb -x p53NTAD.dcd -name p53NTAD
```

Note that the single-trajectory path reads the sequence from `seq.txt` in that folder
(`readpath.readsequence_single`), not from `seq.fasta`.

## Linux paths

The same data was also mounted on Linux during the `pick_ensemble` work:

```
/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma123/puma_wildfull-summary/BB/S_0
/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum/puma_wildfull-summary/BB/S_0
```
