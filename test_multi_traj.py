#!/usr/bin/env python3
"""
Test script for E1A_pat-summary with multiple trajectories.

NOTE: This test currently fails due to topology mismatch:
- E1A_pat-summary has all-atom structure (992 atoms, 79 residues including ACE/NME caps and ions)
- Our code was validated with coarse-grained models (64 atoms, 64 residues)
- MDTraj's compute_contacts has issues with the all-atom topology
  Error: "contacts requests a residue that is not in the permitted range"

The core functionality works correctly with CG models (validated on 4 proteins, 56k+ data points).
This test is kept for future reference when handling all-atom topologies.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Import new functions
from idp_interaction_map import analyze_interaction_map


def test_multi_trajectory():
    """Test E1A_pat-summary with 5 trajectories."""
    
    print("\n" + "="*80)
    print("Testing Multi-Trajectory Analysis: E1A_pat-summary")
    print("="*80)
    
    # Test data path
    data_path = "/Users/fengyu/interaction_map_test/E1A_pat-summary/BB/S_0"
    
    if not Path(data_path).exists():
        print(f"❌ Test data not found at {data_path}")
        print("Skipping multi-trajectory test")
        return False
    
    # Check trajectory files
    traj_files = list(Path(data_path).glob("__traj_*.xtc"))
    print(f"\n📁 Found {len(traj_files)} trajectory files:")
    for traj in sorted(traj_files):
        size_mb = traj.stat().st_size / (1024*1024)
        print(f"   - {traj.name}: {size_mb:.1f} MB")
    
    # Output directory
    output_new = Path("test_output_multi_new")
    output_new.mkdir(exist_ok=True)
    
    print("\n🚀 Running NEW implementation with 5 trajectories...")
    try:
        # Read sequence
        seq_file_path = Path(data_path) / "../../seq.fasta"
        with open(seq_file_path) as f:
            lines = f.read().strip().split('\n')
            seq = lines[-1] if len(lines) > 1 else lines[0]  # Handle FASTA format
        
        df_new = analyze_interaction_map(
            name="E1A_pat_summary",
            trajectory_path=data_path,
            sequence=seq,
            output_dir=str(output_new),
            pdb_top="__START_0.pdb",
            xtc_input=5,
            read_from_file=False
        )
        print(f"   ✓ Completed interaction analysis: {len(df_new)} interactions")
        
    except Exception as e:
        print(f"   ❌ Error in new implementation: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Validate results
    print("\n📊 Validating Results...")
    print("-" * 80)
    
    # Check expected columns
    expected_columns = ['r_1', 'r_2', 'cont_prob', 'distance', 
                       'gs_standard', 'relative_strength', 'plot_value']
    
    missing_cols = set(expected_columns) - set(df_new.columns)
    if missing_cols:
        print(f"❌ Missing columns: {missing_cols}")
        return False
    
    print(f"✓ All expected columns present: {len(expected_columns)} columns")
    
    # Data quality checks
    print(f"\n📈 Data Quality Checks:")
    print(f"   Total interaction pairs: {len(df_new)}")
    print(f"   Residue range: {int(df_new['r_1'].min())}-{int(df_new['r_2'].max())}")
    print(f"   Contact probability range: {df_new['cont_prob'].min():.3f} - {df_new['cont_prob'].max():.3f}")
    print(f"   Distance range: {df_new['distance'].min():.3f} - {df_new['distance'].max():.3f} nm")
    
    # Check for NaN values
    nan_count = df_new.isna().sum().sum()
    if nan_count > 0:
        print(f"   ⚠️  Warning: {nan_count} NaN values found")
    else:
        print(f"   ✓ No NaN values")
    
    # Check plot_value distribution
    plot_value_counts = df_new['plot_value'].value_counts().sort_index()
    print(f"\n   Interaction classification:")
    for val, count in plot_value_counts.items():
        labels = {-2: "Strong unfavorable", -1: "Weak unfavorable", 
                  1: "Weak favorable", 2: "Strong favorable"}
        label = labels.get(int(val), f"Unknown ({val})")
        print(f"     {label:20s}: {count:4d} pairs ({count/len(df_new)*100:5.1f}%)")
    
    # Check output files
    print(f"\n📁 Output Files:")
    csv_file = output_new / "E1A_pat_summary_interaction.csv"
    png_file = output_new / "E1A_pat_summary.png"
    svg_file = output_new / "E1A_pat_summary.svg"
    
    files_exist = []
    if csv_file.exists():
        print(f"   ✓ CSV: {csv_file.name} ({csv_file.stat().st_size} bytes)")
        files_exist.append(True)
    else:
        print(f"   ❌ CSV file not found")
        files_exist.append(False)
        
    if png_file.exists():
        print(f"   ✓ PNG: {png_file.name} ({png_file.stat().st_size / 1024:.1f} KB)")
        files_exist.append(True)
    else:
        print(f"   ❌ PNG file not found")
        files_exist.append(False)
        
    if svg_file.exists():
        print(f"   ✓ SVG: {svg_file.name} ({svg_file.stat().st_size / 1024:.1f} KB)")
        files_exist.append(True)
    else:
        print(f"   ❌ SVG file not found")
        files_exist.append(False)
    
    # Summary
    print("\n" + "="*80)
    if all(files_exist) and len(df_new) > 0 and nan_count == 0:
        print("✅ MULTI-TRAJECTORY TEST PASSED")
        print(f"   - Processed {len(traj_files)} trajectory files successfully")
        print(f"   - Generated {len(df_new)} interaction pairs")
        print(f"   - All output files created")
        print("="*80)
        return True
    else:
        print("❌ MULTI-TRAJECTORY TEST FAILED")
        print("="*80)
        return False


if __name__ == "__main__":
    try:
        success = test_multi_trajectory()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
