"""
Example: Generate mutations from interaction map analysis.

This example demonstrates how to use the mutation scanner to identify
and generate targeted mutations based on interaction analysis.
"""

from idp_interaction_map import scan_mutations_from_csv, MutationScanner
import pandas as pd

# Example 1: Quick mutation scan from CSV
print("Example 1: Quick mutation scan from CSV file")
print("=" * 60)

# Assuming you've already run interaction map analysis
interaction_csv = "protein_interaction.csv"
sequence = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"
protein_name = "MyProtein"
output_dir = "./mutations_attractive"

# Generate attractive mutations
mutations = scan_mutations_from_csv(
    interaction_csv=interaction_csv,
    sequence=sequence,
    protein_name=protein_name,
    output_dir=output_dir,
    interaction_type='attractive',
    forbidden_regions=[1, 2, 3],  # Don't mutate first 3 residues
    min_chunk_strength=1.0
)

print(f"\nGenerated {sum(len(df) for df in mutations.values())} mutations")
for name, df in mutations.items():
    if len(df) > 0:
        print(f"  - {name}: {len(df)} mutations")


# Example 2: Custom mutation generation
print("\n" + "=" * 60)
print("Example 2: Custom mutation generation workflow")
print("=" * 60)

# Load interaction data
df = pd.read_csv(interaction_csv, index_col=0)

# Create scanner
scanner = MutationScanner(
    interaction_df=df,
    sequence=sequence,
    protein_name=protein_name,
    forbidden_regions=[1, 2, 3, 40, 41, 42]  # Protect termini and specific region
)

# Filter significant interactions
selected = scanner.filter_interactions(
    min_contact_prob=0.01,
    min_distance=4,
    min_relative_strength=0.1
)

print(f"\nFound {len(selected)} significant interactions")

# Calculate chunk strength
chunk_data = scanner.calculate_chunk_strength(selected)
print(f"Calculated chunk strength for {len(chunk_data)} interactions")

# Identify candidates for attractive mutations
candidates = scanner.identify_candidates(
    chunk_data,
    interaction_type='attractive',
    min_chunk_strength=1.5  # Higher threshold
)

print(f"Identified {len(candidates)} strong attractive candidates")

# Generate specific types of mutations
print("\nGenerating targeted mutations...")

# Single charge mutations on left side
left_charge = scanner.generate_single_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    interaction_type='attractive'
)
print(f"  - Left side charge mutations: {len(left_charge)}")

# Pair mutations for both residues
pair_charge = scanner.generate_pair_mutations(
    candidates,
    mutation_type='charge',
    interaction_type='attractive'
)
print(f"  - Pair charge mutations: {len(pair_charge)}")

# Chunk mutations (3 residues)
chunk_mutations = scanner.generate_chunk_mutations(
    candidates,
    position='left',
    mutation_type='charge',
    interaction_type='attractive',
    chunk_size=3
)
print(f"  - Chunk mutations (3 residues): {len(chunk_mutations)}")

# Save mutations
output_dir2 = "./custom_mutations"
scanner.save_mutations(left_charge, f"{output_dir2}/left_charge.csv")
scanner.save_mutations(pair_charge, f"{output_dir2}/pair_charge.csv")
scanner.save_mutations(chunk_mutations, f"{output_dir2}/chunk_mutations.csv")

print(f"\nMutations saved to {output_dir2}/")


# Example 3: Repulsive mutations to disrupt interactions
print("\n" + "=" * 60)
print("Example 3: Generate repulsive mutations")
print("=" * 60)

# Generate mutations to disrupt unfavorable interactions
repulsive_mutations = scan_mutations_from_csv(
    interaction_csv=interaction_csv,
    sequence=sequence,
    protein_name=protein_name,
    output_dir="./mutations_repulsive",
    interaction_type='repulsive',
    min_chunk_strength=1.0
)

print(f"\nGenerated {sum(len(df) for df in repulsive_mutations.values())} repulsive mutations")


# Example 4: Analyze generated mutations
print("\n" + "=" * 60)
print("Example 4: Analyze generated mutations")
print("=" * 60)

if len(left_charge) > 0:
    print("\nFirst 5 left-side charge mutations:")
    print(left_charge.head())
    
    print("\nMutation statistics:")
    print(f"  - Total mutations: {len(left_charge)}")
    print(f"  - Unique positions: {left_charge['mutation_name'].str.extract(r'(\d+)')[0].nunique()}")
    
    # Example: Get mutations for a specific position
    # Extract position from mutation name (format: ProteinName_R10E)
    left_charge['position'] = left_charge['mutation_name'].str.extract(r'[A-Z](\d+)[A-Z]')
    
    print("\nMutations targeting position 10:")
    pos_10_muts = left_charge[left_charge['position'] == '10']
    if len(pos_10_muts) > 0:
        print(pos_10_muts[['mutation_name', 'sequence']])


print("\n" + "=" * 60)
print("Mutation generation complete!")
print("=" * 60)
print("\nNext steps:")
print("1. Review generated mutation CSV files")
print("2. Select mutations for experimental validation")
print("3. Use sequences for protein expression and testing")
