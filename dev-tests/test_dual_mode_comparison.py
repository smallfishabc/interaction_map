#!/usr/bin/env python3
"""
Integration test: Compare CG and All-Atom modes

This test runs the same protein through both CG and all-atom modes
to demonstrate the differences in analysis results.

Author: Feng Yu
Date: December 15, 2025
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt


def test_dual_mode_comparison():
    """Compare CG vs all-atom analysis on the same data."""
    
    print("=" * 70)
    print("DUAL-MODE COMPARISON TEST")
    print("=" * 70)
    
    # Test data
    data_path = Path("BB/S_0")
    if not data_path.exists():
        print(f"⚠️  Test data not found: {data_path}")
        print("This test requires E1A_pat-summary data in BB/S_0/")
        return
    
    sequence = read_sequence_from_txt(str(data_path))
    print(f"\nProtein: E1A_pat")
    print(f"Sequence length: {len(sequence)} residues")
    print(f"Trajectories: 5 files (__traj_0.xtc through __traj_4.xtc)")
    
    # Run CG analysis
    print("\n" + "-" * 70)
    print("RUNNING COARSE-GRAINED ANALYSIS")
    print("-" * 70)
    
    df_cg = analyze_interaction_map(
        name="E1A_pat_cg",
        trajectory_path=str(data_path),
        sequence=sequence,
        output_dir="test_output_comparison/cg",
        xtc_input=5,
        use_ca=False  # CG mode
    )
    
    print(f"\n✅ CG Analysis Complete")
    print(f"   Pairs analyzed: {len(df_cg)}")
    print(f"   Strong favorable: {len(df_cg[df_cg['plot_value'] == 2])}")
    print(f"   Moderate favorable: {len(df_cg[df_cg['plot_value'] == 1])}")
    print(f"   Moderate unfavorable: {len(df_cg[df_cg['plot_value'] == -1])}")
    print(f"   Strong unfavorable: {len(df_cg[df_cg['plot_value'] == -2])}")
    
    # Run All-Atom analysis
    print("\n" + "-" * 70)
    print("RUNNING ALL-ATOM ANALYSIS")
    print("-" * 70)
    
    df_aa = analyze_interaction_map(
        name="E1A_pat_aa",
        trajectory_path=str(data_path),
        sequence=sequence,
        output_dir="test_output_comparison/all_atom",
        xtc_input=5,
        use_ca=True  # All-atom mode
    )
    
    print(f"\n✅ All-Atom Analysis Complete")
    print(f"   Pairs analyzed: {len(df_aa)}")
    print(f"   Strong favorable: {len(df_aa[df_aa['plot_value'] == 2])}")
    print(f"   Moderate favorable: {len(df_aa[df_aa['plot_value'] == 1])}")
    print(f"   Moderate unfavorable: {len(df_aa[df_aa['plot_value'] == -1])}")
    print(f"   Strong unfavorable: {len(df_aa[df_aa['plot_value'] == -2])}")
    
    # Comparison
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    
    print(f"\n1. Contact Pair Counts:")
    print(f"   CG:        {len(df_cg):5d} pairs")
    print(f"   All-Atom:  {len(df_aa):5d} pairs")
    print(f"   Difference: {abs(len(df_cg) - len(df_aa)):5d} pairs")
    
    print(f"\n2. Favorable Interaction Counts:")
    cg_favorable = len(df_cg[df_cg['relative_strength'] > 0.8])
    aa_favorable = len(df_aa[df_aa['relative_strength'] > 0.8])
    print(f"   CG:        {cg_favorable:5d} favorable")
    print(f"   All-Atom:  {aa_favorable:5d} favorable")
    print(f"   Difference: {abs(cg_favorable - aa_favorable):5d}")
    
    print(f"\n3. Unfavorable Interaction Counts:")
    cg_unfavorable = len(df_cg[df_cg['relative_strength'] < -0.8])
    aa_unfavorable = len(df_aa[df_aa['relative_strength'] < -0.8])
    print(f"   CG:        {cg_unfavorable:5d} unfavorable")
    print(f"   All-Atom:  {aa_unfavorable:5d} unfavorable")
    print(f"   Difference: {abs(cg_unfavorable - aa_unfavorable):5d}")
    
    print(f"\n4. Mean Relative Strength:")
    print(f"   CG:        {df_cg['relative_strength'].mean():7.3f}")
    print(f"   All-Atom:  {df_aa['relative_strength'].mean():7.3f}")
    
    print(f"\n5. Top 5 Strongest Interactions (CG mode):")
    top_cg = df_cg.nlargest(5, 'relative_strength')
    for idx, row in top_cg.iterrows():
        print(f"   {int(row['r_1']):3d}-{int(row['r_2']):3d}: "
              f"R = {row['relative_strength']:6.3f}")
    
    print(f"\n6. Top 5 Strongest Interactions (All-Atom mode):")
    top_aa = df_aa.nlargest(5, 'relative_strength')
    for idx, row in top_aa.iterrows():
        print(f"   {int(row['r_1']):3d}-{int(row['r_2']):3d}: "
              f"R = {row['relative_strength']:6.3f}")
    
    # Analysis notes
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    print("""
Expected Differences:
- CG mode uses all atoms/beads (one per residue)
- All-atom mode uses CA atoms only (from full atomic detail)
- Different normalization parameters reflect different physics
- Number of pairs should be similar (~2016 for 64 residues)
- Interaction categorization may differ due to different models

Output Files Generated:
- test_output_comparison/cg/E1A_pat_cg_interaction.csv
- test_output_comparison/cg/E1A_pat_cg.png
- test_output_comparison/cg/E1A_pat_cg.svg
- test_output_comparison/all_atom/E1A_pat_aa_interaction.csv
- test_output_comparison/all_atom/E1A_pat_aa.png
- test_output_comparison/all_atom/E1A_pat_aa.svg

You can visually compare the .png files to see differences in
interaction networks between CG and all-atom analyses.
    """)
    
    print("=" * 70)
    print("✅ DUAL-MODE COMPARISON TEST COMPLETED")
    print("=" * 70)
    
    return df_cg, df_aa


if __name__ == "__main__":
    try:
        df_cg, df_aa = test_dual_mode_comparison()
        print("\n✅ Test completed successfully!")
        print(f"   CG results: {len(df_cg)} pairs")
        print(f"   All-atom results: {len(df_aa)} pairs")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
