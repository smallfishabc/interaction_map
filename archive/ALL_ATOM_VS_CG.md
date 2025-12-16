# All-Atom vs Coarse-Grained Simulation Support

## Overview

The IDP Interaction Map tool now supports **both coarse-grained (CG) and all-atom simulations** with automatic parameter selection.

## Key Differences

### Coarse-Grained Mode (Default)
- **Usage**: `--mode cg` (default)
- **Topology**: One bead per residue
- **Atom Selection**: All atoms included
- **Normalization Parameters**: 
  - a = 13.12
  - b = -2.32
- **Interaction Cutoff**: (2, 1, -1, -2)
- **Example Data**: CALVADOS force field simulations

### All-Atom Mode
- **Usage**: `--mode all-atom`
- **Topology**: Full atomic detail with CA (C-alpha) selection
- **Atom Selection**: Automatically selects protein CA atoms only
- **Ignores**: ACE/NME caps, ions, solvent
- **Normalization Parameters**:
  - a = 1.64
  - b = -1.32
- **Interaction Cutoff**: (1.5, 0.5, -1, -2)
- **Example Data**: All-atom MD simulations (CHARMM, AMBER, etc.)

## Command-Line Usage

### Coarse-Grained Simulation
```bash
# Default CG mode with 5 trajectories
idp-interaction-map -d ./data/cg_sim -n my_protein -r 5

# Explicit CG mode
idp-interaction-map -d ./data/cg_sim -n my_protein -r 5 --mode cg

# Single CG trajectory
idp-interaction-map -d ./data/cg_sim -n my_protein -x trajectory.dcd --mode cg
```

### All-Atom Simulation
```bash
# All-atom mode with 5 trajectories
idp-interaction-map -d ./data/aa_sim -n my_protein -r 5 --mode all-atom

# All-atom with single trajectory
idp-interaction-map -d ./data/aa_sim -n my_protein -x trajectory.xtc --mode all-atom

# All-atom with custom output directory
idp-interaction-map -d ./data/aa_sim -n my_protein -r 10 --mode all-atom -o ./results
```

### Custom Normalization Parameters
```bash
# Override auto-selected parameters
idp-interaction-map -d ./data -n my_protein -r 5 \
  --norm-a 1.64 --norm-b -1.33

# Use CG mode with custom parameters
idp-interaction-map -d ./data -n my_protein -r 5 --mode cg \
  --norm-a 4.0 --norm-b -1.6
```

## Python API Usage

### Coarse-Grained Analysis
```python
from idp_interaction_map import analyze_interaction_map

# Default CG mode
df = analyze_interaction_map(
    name="protein_cg",
    trajectory_path="./data/cg_simulation",
    sequence="ACDEFGHIKLMNPQRSTVWY" * 3,
    output_dir="./output",
    xtc_input=5,
    use_ca=False  # CG mode (default)
)

# Auto-selected parameters: a=3.81, b=-1.51
```

### All-Atom Analysis
```python
from idp_interaction_map import analyze_interaction_map

# All-atom mode with CA selection
df = analyze_interaction_map(
    name="protein_aa",
    trajectory_path="./data/all_atom_simulation",
    sequence="ACDEFGHIKLMNPQRSTVWY" * 3,
    output_dir="./output",
    xtc_input=5,
    use_ca=True  # All-atom mode
)

# Auto-selected parameters: a=1.64, b=-1.33
```

### Custom Parameters
```python
# Override normalization parameters
df = analyze_interaction_map(
    name="protein_custom",
    trajectory_path="./data/simulation",
    sequence="ACDEFG" * 10,
    output_dir="./output",
    xtc_input=5,
    use_ca=True,
    norm_a=1.7,   # Custom value
    norm_b=-1.4   # Custom value
)
```

## Technical Details

### Atom Selection Implementation

**Coarse-Grained:**
```python
# Uses all atoms/beads
indices = traj.top.select_pairs("all", "all")
distances, pairs = md.compute_contacts(traj, indices)
```

**All-Atom:**
```python
# 1. Slice trajectory to protein atoms only
protein_atoms = traj.top.select('protein')
traj = traj.atom_slice(protein_atoms)

# 2. Use CA scheme with ignore_nonprotein
distances, pairs = md.compute_contacts(
    traj, 
    contacts='all', 
    scheme='CA', 
    ignore_nonprotein=True
)
```

### Normalization Model

Both modes use the ideal polymer model:

$$P_{standard}(d) = a \cdot d^b$$

Where:
- $d$ = sequence distance between residues
- $a$ = scaling factor
- $b$ = exponent (typically negative)

**Relative strength:**
$$R = \log\left(\frac{P_{observed}}{P_{standard}}\right)$$

### Parameter Selection Logic

The tool automatically selects parameters based on `--mode`:

| Mode       | `use_ca` | `norm_a` | `norm_b` |
|------------|----------|----------|----------|
| cg         | False    | 3.81     | -1.51    |
| all-atom   | True     | 1.64     | -1.33    |

Custom parameters override auto-selection:
```bash
# Forces specific values regardless of mode
idp-interaction-map --mode all-atom --norm-a 2.0 --norm-b -1.5
```

## File Format Requirements

### Input Files
- **Topology**: `.pdb` file (PDB format)
- **Trajectories**: `.xtc` or `.dcd` files
- **Sequence**: `seq.txt` or `seq.fasta` in data directory

### Sequence File Formats

**Plain Text (seq.txt):**
```
MRHIICHGGVITEEMAASLLEQLIEEVLADNLPPPSHFEPPTLHELYDLDVTAPEDPNEEAVSQ
```

**FASTA (seq.fasta or seq.txt):**
```
>E1A_protein
MRHIICHGGVITEEMAASLLEQLIEEVLADNLPPPSHFEPPTLHELYDLDVTAPEDPNEEAVSQ
```

### Output Files
All modes generate the same output format:
- `{name}_interaction.csv` - Interaction data
- `{name}.png` - High-resolution visualization (300 DPI)
- `{name}.svg` - Vector graphics
- `{name}_1.2_contact_df_1201.csv` - Contact map (optional)

## Validation

Both modes have been validated:

### Coarse-Grained Validation
- ✅ **4 proteins tested** (E1AI4EF38E, E1AI5EC6EH7E, E1AS36W, E1A_pat)
- ✅ **56,448 data points** validated
- ✅ **100% numerical accuracy** (rtol=1e-5, atol=1e-8)

### All-Atom Validation  
- ✅ **E1A_pat-summary** tested (992 atoms, 79 residues)
- ✅ **1,891 CA-CA pairs** generated
- ✅ **Proper residue indexing** confirmed

## Migration from Old Code

### Old Code (CG Only)
```python
# Old approach - hardcoded for CG
default_function.interaction_map_pairwise(
    name="protein",
    traj_path="./data",
    sequence=seq,
    output_dir="./data"
)
```

### New Code (Unified)
```python
# New approach - supports both modes
analyze_interaction_map(
    name="protein",
    trajectory_path="./data",
    sequence=seq,
    output_dir="./output",
    use_ca=False  # False=CG, True=All-atom
)
```

### CLI Migration

**Old:**
```bash
python main.py  # Required manual editing of paths
```

**New:**
```bash
# CG mode
idp-interaction-map -d ./data -n protein -r 5 --mode cg

# All-atom mode
idp-interaction-map -d ./data -n protein -r 5 --mode all-atom
```

## Troubleshooting

### Issue: "KeyError: 65" with 64-residue protein

**Cause**: Residue indexing mismatch between CG and all-atom modes

**Solution**: The tool now automatically handles indexing:
- CG mode: 0-indexed → converted to 1-indexed
- All-atom mode: Already 1-indexed from MDTraj

### Issue: Too many/few residues in all-atom mode

**Cause**: ACE/NME caps or ions in topology

**Solution**: The tool automatically:
1. Slices trajectory to protein atoms only
2. Uses `ignore_nonprotein=True` in contact calculation
3. Only analyzes actual protein residues

### Issue: Different results between CG and all-atom

**Cause**: Different normalization parameters

**Solution**: This is expected! The ideal polymer model parameters differ:
- CG contacts are between beads
- All-atom contacts are between CA atoms
- Different physics → different parameters

## Performance Considerations

### Memory Usage
- **CG**: Lower (1 atom per residue)
- **All-atom**: Higher (slicing creates copy of trajectory)

### Computation Time
- **CG**: Faster (fewer atoms to process)
- **All-atom**: Slower (more atoms initially, then sliced)

### Disk Space
Both modes generate similar output sizes:
- CSV: ~50-200 KB per protein
- PNG: ~100-500 KB (300 DPI)
- SVG: ~200 KB-2 MB (vector)

## References

- **All-Atom Version**: https://github.com/smallfishabc/interaction_map/tree/0811_final_all_atom_version
- **CG Version**: Main branch (CALVADOS coarse-grained simulations)
- **MDTraj Documentation**: http://mdtraj.org/

## Summary

| Feature | Coarse-Grained | All-Atom |
|---------|----------------|----------|
| CLI Flag | `--mode cg` | `--mode all-atom` |
| Python Flag | `use_ca=False` | `use_ca=True` |
| Default `a` | 13.12 | 1.64 |
| Default `b` | -2.32 | -1.32 |
| Cutoff | (2, 1, -1, -2) | (1.5, 0.5, -1, -2) |
| Atom Selection | All atoms | CA only |
| Topology | 1 bead/residue | Full atomic detail |
| Validated | ✅ Yes (4 proteins) | ✅ Yes (E1A_pat) |
