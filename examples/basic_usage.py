"""Example script demonstrating basic usage of IDP Interaction Map."""

from pathlib import Path

from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt


def main():
    """Run example analysis."""
    # Set up paths
    data_dir = Path("./data")
    output_dir = Path("./output")
    output_dir.mkdir(exist_ok=True)

    # Read protein sequence
    print("Reading protein sequence...")
    sequence = read_sequence_from_txt(data_dir)
    print(f"Sequence length: {len(sequence)}")

    # Run analysis
    print("\nAnalyzing interactions...")
    interaction_df = analyze_interaction_map(
        name="example_protein",
        trajectory_path=data_dir,
        sequence=sequence,
        output_dir=output_dir,
        pdb_top="__START_0.pdb",
        xtc_input=5,  # Use 5 trajectory files
        read_from_file=False,
    )

    # Print summary statistics
    print("\n=== Analysis Summary ===")
    print(f"Total residue pairs: {len(interaction_df)}")
    print(f"Strong favorable interactions: {(interaction_df['plot_value'] == 2).sum()}")
    print(f"Weak favorable interactions: {(interaction_df['plot_value'] == 1).sum()}")
    print(f"Weak unfavorable interactions: {(interaction_df['plot_value'] == -1).sum()}")
    print(f"Strong unfavorable interactions: {(interaction_df['plot_value'] == -2).sum()}")

    # Show top interactions
    print("\n=== Top 5 Strongest Contacts ===")
    top_contacts = interaction_df.nlargest(5, "cont_prob")
    for _, row in top_contacts.iterrows():
        print(
            f"Residues {int(row['r_1'])}-{int(row['r_2'])}: "
            f"P={row['cont_prob']:.3f}, Strength={row['relative_strength']:.2f}"
        )

    print("\nAnalysis complete! Check the output directory for visualizations.")


if __name__ == "__main__":
    main()
