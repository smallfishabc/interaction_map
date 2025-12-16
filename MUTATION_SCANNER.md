# Mutation Scanner - User Guide

The Mutation Scanner is a powerful tool for generating targeted mutations based on interaction map analysis. It helps identify key residues driving protein interactions and suggests mutations to modulate these interactions.

## Overview

The mutation scanner analyzes interaction maps to:
1. **Identify strong interactions** between residue pairs or chunks
2. **Calculate chunk strength** - weighted interaction strength for groups of adjacent residues  
3. **Generate mutations** to either enhance attractive interactions or disrupt repulsive ones
4. **Export mutation libraries** ready for experimental validation

## Installation

The mutation scanner is included in the IDP Interaction Map package:

```bash
pip install -e .
```

After installation, two CLI commands are available:
- `idp-interaction-map` - Generate interaction maps
- `idp-mutation-scan` - Generate mutations from interaction maps

## Quick Start

### Command Line Usage

```bash
# Generate attractive mutations
idp-mutation-scan \
  -i protein_interaction.csv \
  -s ACDEFGHIKLMNPQRSTVWY \
  -n MyProtein \
  -o ./mutations \
  --type attractive

# Generate repulsive mutations with forbidden regions
idp-mutation-scan \
  -i protein_interaction.csv \
  -f seq.txt \
  -n MyProtein \
  -o ./mutations \
  --type repulsive \
  --forbidden 1-5,40-45

# Custom chunk strength threshold
idp-mutation-scan \
  -i protein_interaction.csv \
  -s ACDEFGHIKLMNPQRSTVWY \
  -n MyProtein \
  -o ./mutations \
  --min-strength 2.0
```

### Python API Usage

```python
from idp_interaction_map import scan_mutations_from_csv

# Quick scan
mutations = scan_mutations_from_csv(
    interaction_csv="protein_interaction.csv",
    sequence="ACDEFGHIKLMNPQRSTVWY",
    protein_name="MyProtein",
    output_dir="./mutations",
    interaction_type='attractive',
    forbidden_regions=[1, 2, 3],
    min_chunk_strength=1.0
)

# Returns dict of DataFrames with generated mutations
print(f"Generated {sum(len(df) for df in mutations.values())} mutations")
```

## Detailed Usage

### 1. Interaction Analysis

First, generate an interaction map using the main tool:

```bash
idp-interaction-map \
  -d ./trajectory_data \
  -n MyProtein \
  -r 5 \
  --mode cg
```

This creates `MyProtein_interaction.csv` with columns:
- `r_1`, `r_2`: Residue pair
- `cont_prob`: Contact probability
- `distance`: Sequence separation
- `relative_strength`: Normalized interaction strength
- `plot_value`: Interaction classification (2, 1, -1, -2)

### 2. Mutation Generation

#### Attractive Mutations
Enhance favorable interactions by:
- Adding complementary charges (E/K to attract opposite partners)
- Introducing polar residues (Q for hydrogen bonding)
- Increasing hydrophobicity (L for hydrophobic core)

```python
from idp_interaction_map import MutationScanner
import pandas as pd

# Load interaction data
df = pd.read_csv("protein_interaction.csv", index_col=0)

# Create scanner
scanner = MutationScanner(
    interaction_df=df,
    sequence="ACDEFG...",
    protein_name="MyProtein",
    forbidden_regions=[1, 2, 3]  # Protect specific regions
)

# Filter significant interactions
selected = scanner.filter_interactions(
    min_contact_prob=0.01,
    min_distance=4,
    min_relative_strength=0.1
)

# Calculate chunk strength
chunk_data = scanner.calculate_chunk_strength(selected)

# Identify attractive candidates
candidates = scanner.identify_candidates(
    chunk_data,
    interaction_type='attractive',
    min_chunk_strength=1.0
)

# Generate mutations
mutations = scanner.generate_single_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    interaction_type='attractive'
)
```

#### Repulsive Mutations
Disrupt unfavorable interactions by:
- Neutralizing charges (to Ala)
- Reducing hydrophobicity (to Ser)
- Removing polar interactions

```python
# Identify repulsive candidates
candidates = scanner.identify_candidates(
    chunk_data,
    interaction_type='repulsive',
    min_chunk_strength=1.0
)

# Generate neutralizing mutations
mutations = scanner.generate_single_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    interaction_type='repulsive'
)
```

### 3. Mutation Types

#### Single Mutations
Mutate one residue in each interaction:

```python
# Left side (r_1) charge mutations
left_charge = scanner.generate_single_mutations(
    candidates, position='left', mutation_type='charge'
)

# Right side (r_2) polar mutations
right_polar = scanner.generate_single_mutations(
    candidates, position='right', mutation_type='polar'
)
```

#### Pair Mutations
Mutate both residues in an interaction:

```python
# Both residues to same charge
pair_charge = scanner.generate_pair_mutations(
    candidates, mutation_type='charge'
)

# Both to polar
pair_polar = scanner.generate_pair_mutations(
    candidates, mutation_type='polar'
)
```

#### Chunk Mutations
Mutate 3-residue chunks (residue ± 1):

```python
# Mutate left chunk
chunk_muts = scanner.generate_chunk_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    chunk_size=3
)
```

### 4. Forbidden Regions

Protect specific regions from mutation:

```python
# Protect N-terminus (1-5) and C-terminus (60-64)
forbidden = list(range(1, 6)) + list(range(60, 65))

scanner = MutationScanner(
    interaction_df=df,
    sequence=sequence,
    protein_name=name,
    forbidden_regions=forbidden
)
```

CLI format:
```bash
# Comma-separated
--forbidden 1,2,3,60,61,62

# Ranges
--forbidden 1-5,60-64

# Mixed
--forbidden 1-5,10,20-25,60-64
```

### 5. Full Mutation Scan

Generate all mutation types at once:

```python
results = scanner.generate_full_scan(
    output_dir="./mutations",
    interaction_type='attractive',
    min_chunk_strength=1.0
)

# Results contains:
# - left_single_charge
# - left_single_polar
# - right_single_charge
# - right_single_polar
# - pair_charge
# - pair_polar
# - left_chunk_charge
# - right_chunk_charge
```

## Output Format

Mutations are saved as CSV files with two columns:
1. **Mutation name**: `ProteinName_R10E` (R at position 10 → E)
2. **Sequence**: Full mutated sequence

Example:
```
MyProtein_R10E,ACDEFGHIKE...
MyProtein_K15A,ACDEFGHIKA...
```

Multiple mutations in one sequence:
```
MyProtein_R10E_K15E,ACDEFGHIEE...
```

## Chunk Strength Calculation

Chunk strength quantifies how strongly two chunks interact:

1. **Define chunks**: 3 adjacent residues centered on each residue
2. **Calculate pairwise**: All 9 interactions between chunks (3x3)
3. **Weight by distance**: Closer to center = higher weight
   - Center-center: weight = 1.0
   - One away: weight = 0.5
   - Two away: weight = 0.25
4. **Sum weighted strengths**

Example for residues 10-15 interaction:
```
Chunk 1: [9, 10, 11]
Chunk 2: [14, 15, 16]

Strength = Σ (interaction_strength / 2^distance_from_center)
```

## Mutation Strategy Guide

### For Attractive Interactions (Enhancing)

| Current Type | Target | Mutation | Rationale |
|--------------|--------|----------|-----------|
| Positive charge | Negative partner | → E | Attract opposite |
| Negative charge | Positive partner | → K | Attract opposite |
| Hydrophobic | Hydrophobic partner | → L | Enhance core |
| Polar | Polar partner | → Q | H-bonding |

### For Repulsive Interactions (Disrupting)

| Current Type | Mutation | Rationale |
|--------------|----------|-----------|
| Charged | → A | Neutralize |
| Hydrophobic | → S | Reduce hydrophobicity |
| Polar | → A | Remove H-bonding |

## Advanced Features

### Custom Filtering

```python
# Stricter filtering
selected = scanner.filter_interactions(
    min_contact_prob=0.05,  # Higher threshold
    min_distance=5,         # More sequence separation
    min_relative_strength=0.2,  # Stronger interactions
    exclude_termini=True    # Skip terminal residues
)
```

### Chunk Size Variation

```python
# Larger chunks (5 residues)
chunk_data = scanner.calculate_chunk_strength(
    selected,
    chunk_size=5
)

# Generate 5-residue chunk mutations
mutations = scanner.generate_chunk_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    chunk_size=5
)
```

### Custom Properties

Access calculated properties:
```python
# Chunk data includes:
chunk_data['r_1_hydro']     # Hydrophobicity
chunk_data['r_1_charge']    # Net charge
chunk_data['r_1_aromatic']  # Aromatic count
chunk_data['chunk_strength'] # Overall strength

# Filter by custom criteria
high_charge = chunk_data[abs(chunk_data['r_1_charge']) >= 2]
```

## Integration with Workflow

Typical workflow:

```bash
# 1. Generate interaction map
idp-interaction-map -d ./data -n Protein -r 5

# 2. Generate attractive mutations
idp-mutation-scan \
  -i ./output/Protein_interaction.csv \
  -f ./data/seq.txt \
  -n Protein \
  -o ./mutations_attractive \
  --type attractive

# 3. Generate repulsive mutations
idp-mutation-scan \
  -i ./output/Protein_interaction.csv \
  -f ./data/seq.txt \
  -n Protein \
  -o ./mutations_repulsive \
  --type repulsive
```

## Examples

See `examples/generate_mutations.py` for comprehensive examples including:
- Quick mutation scanning
- Custom mutation workflows
- Analyzing generated mutations
- Filtering and selection strategies

## Troubleshooting

**No mutations generated?**
- Lower `min_chunk_strength` threshold
- Check `plot_value` distribution in interaction CSV
- Verify sequence matches interaction data
- Check forbidden regions aren't excluding everything

**Too many mutations?**
- Increase `min_chunk_strength` threshold
- Add forbidden regions
- Filter by interaction distance or contact probability
- Focus on specific mutation types

**Import errors?**
- Ensure package is installed: `pip install -e .`
- Check Python version >= 3.8
- Verify all dependencies installed

## Citation

If you use the mutation scanner in your research, please cite:

```
Yu, F. et al. (2024). IDP Interaction Map: Analyzing and Modulating
Interactions in Intrinsically Disordered Proteins.
```

## Support

For questions or issues:
- GitHub Issues: https://github.com/yourusername/interaction_map
- Email: your.email@example.com
