"""Example: Comparing multiple protein variants."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt


def compare_proteins(protein_configs):
    """
    Compare interaction patterns across multiple proteins.

    Args:
        protein_configs: List of dicts with protein configuration
    """
    results = {}

    for config in protein_configs:
        print(f"\nAnalyzing {config['name']}...")

        sequence = read_sequence_from_txt(config["data_dir"])

        interaction_df = analyze_interaction_map(
            name=config["name"],
            trajectory_path=config["data_dir"],
            sequence=sequence,
            output_dir=config["output_dir"],
            xtc_input=config.get("xtc_input", 5),
        )

        results[config["name"]] = interaction_df

    return results


def plot_comparison(results):
    """Plot comparison of interaction strengths."""
    fig, ax = plt.subplots(figsize=(10, 6))

    for name, df in results.items():
        favorable = (df["plot_value"] > 0).sum()
        unfavorable = (df["plot_value"] < 0).sum()
        neutral = (df["plot_value"] == 0).sum()

        print(f"\n{name}:")
        print(f"  Favorable: {favorable}")
        print(f"  Unfavorable: {unfavorable}")
        print(f"  Neutral: {neutral}")

    plt.tight_layout()
    plt.savefig("protein_comparison.png", dpi=300)
    print("\nComparison plot saved to protein_comparison.png")


def main():
    """Run comparison analysis."""
    # Define proteins to compare
    configs = [
        {
            "name": "wildtype",
            "data_dir": Path("./data/wildtype"),
            "output_dir": Path("./output/wildtype"),
            "xtc_input": 5,
        },
        {
            "name": "mutant_1",
            "data_dir": Path("./data/mutant1"),
            "output_dir": Path("./output/mutant1"),
            "xtc_input": 5,
        },
    ]

    # Run comparison
    results = compare_proteins(configs)
    plot_comparison(results)


if __name__ == "__main__":
    main()
