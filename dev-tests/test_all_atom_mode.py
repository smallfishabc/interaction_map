#!/usr/bin/env python3
"""
Test script to verify all-atom mode works with E1A_pat-summary data.
"""

from pathlib import Path
from idp_interaction_map import analyze_interaction_map

# Test data path
data_path = "/Users/fengyu/interaction_map_test/E1A_pat-summary/BB/S_0"

if not Path(data_path).exists():
    print(f"❌ Test data not found at {data_path}")
    exit(1)

# Read sequence
seq_file = Path(data_path) / "../../seq.fasta"
with open(seq_file) as f:
    lines = f.read().strip().split('\n')
    sequence = lines[-1] if len(lines) > 1 else lines[0]

print(f"Sequence: {sequence}")
print(f"Length: {len(sequence)}")

# Output directory
output_dir = Path("test_all_atom_output")
output_dir.mkdir(exist_ok=True)

print("\n" + "="*80)
print("Testing All-Atom Mode with CA Selection")
print("="*80)

try:
    df = analyze_interaction_map(
        name="E1A_pat_all_atom",
        trajectory_path=data_path,
        sequence=sequence,
        output_dir=str(output_dir),
        pdb_top="__START_0.pdb",
        xtc_input=5,
        read_from_file=False,
        use_ca=True,  # All-atom mode with CA selection
        norm_a=1.64,  # All-atom normalization parameters
        norm_b=-1.33
    )
    
    print(f"\n✅ Analysis completed successfully!")
    print(f"   Generated {len(df)} interaction pairs")
    print(f"   Columns: {list(df.columns)}")
    print(f"\n📊 Sample data:")
    print(df.head())
    
    print(f"\n📁 Output files:")
    for file in output_dir.glob("E1A_pat_all_atom*"):
        print(f"   - {file.name}")
    
    print("\n" + "="*80)
    print("✅ ALL-ATOM MODE TEST PASSED")
    print("="*80)
    
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
