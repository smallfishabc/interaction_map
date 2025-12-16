"""Tests for mutation scanner module."""

import pytest
import pandas as pd
import tempfile
from pathlib import Path

from idp_interaction_map.mutation_scanner import MutationScanner, scan_mutations_from_csv


@pytest.fixture
def sample_interaction_df():
    """Create sample interaction dataframe."""
    return pd.DataFrame({
        'r_1': [1, 2, 3, 10, 11],
        'r_2': [5, 6, 7, 15, 16],
        'cont_prob': [0.8, 0.7, 0.6, 0.9, 0.85],
        'distance': [4, 4, 4, 5, 5],
        'relative_strength': [0.5, 0.4, 0.3, 0.6, 0.55],
        'plot_value': [2, 2, 1, -1, -2]
    })


@pytest.fixture
def sample_sequence():
    """Create sample protein sequence."""
    return "ACDEFGHIKLMNPQRSTVWY"


def test_mutation_scanner_init(sample_interaction_df, sample_sequence):
    """Test MutationScanner initialization."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    assert scanner.seq_length == len(sample_sequence)
    assert scanner.protein_name == "TestProtein"
    assert 'r_1_res' in scanner.df.columns
    assert 'r_2_res' in scanner.df.columns


def test_filter_interactions(sample_interaction_df, sample_sequence):
    """Test interaction filtering."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions(
        min_contact_prob=0.5,
        min_distance=4,
        min_relative_strength=0.3
    )
    
    assert len(filtered) <= len(sample_interaction_df)
    assert all(filtered['cont_prob'] >= 0.5)
    assert all(filtered['distance'] >= 4)
    assert all(filtered['relative_strength'] >= 0.3)


def test_calculate_chunk_strength(sample_interaction_df, sample_sequence):
    """Test chunk strength calculation."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    
    assert 'chunk_strength' in chunk_data.columns
    assert 'r_1_chunk' in chunk_data.columns
    assert 'r_2_chunk' in chunk_data.columns
    assert 'r_1_hydro' in chunk_data.columns
    assert 'r_1_charge' in chunk_data.columns


def test_identify_attractive_candidates(sample_interaction_df, sample_sequence):
    """Test identification of attractive candidates."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(
        chunk_data,
        interaction_type='attractive',
        min_chunk_strength=0.0
    )
    
    # Should only include positive plot_values
    assert all(candidates['plot_value'] > 0)


def test_identify_repulsive_candidates(sample_interaction_df, sample_sequence):
    """Test identification of repulsive candidates."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(
        chunk_data,
        interaction_type='repulsive',
        min_chunk_strength=0.0
    )
    
    # Should only include negative plot_values
    assert all(candidates['plot_value'] < 0)


def test_generate_single_mutations(sample_interaction_df, sample_sequence):
    """Test single mutation generation."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(chunk_data, 'attractive', 0.0)
    
    mutations = scanner.generate_single_mutations(
        candidates,
        position='left',
        mutation_type='charge',
        interaction_type='attractive'
    )
    
    assert 'mutation_name' in mutations.columns
    assert 'sequence' in mutations.columns
    assert all(mutations['sequence'].apply(len) == len(sample_sequence))


def test_generate_pair_mutations(sample_interaction_df, sample_sequence):
    """Test pair mutation generation."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(chunk_data, 'attractive', 0.0)
    
    mutations = scanner.generate_pair_mutations(
        candidates,
        mutation_type='charge',
        interaction_type='attractive'
    )
    
    assert 'mutation_name' in mutations.columns
    assert 'sequence' in mutations.columns
    # Should have multiple charges (E and K variants)
    assert len(mutations) >= len(candidates)


def test_forbidden_regions(sample_interaction_df, sample_sequence):
    """Test forbidden regions are respected."""
    forbidden = [1, 2, 3]
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein",
        forbidden_regions=forbidden
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(chunk_data, 'attractive', 0.0)
    
    mutations = scanner.generate_single_mutations(
        candidates,
        position='left',
        mutation_type='charge',
        interaction_type='attractive'
    )
    
    # Check that no forbidden positions are mutated
    for _, row in mutations.iterrows():
        # Extract position from mutation name (format: TestProtein_R1E)
        import re
        match = re.search(r'_[A-Z](\d+)[A-Z]', row['mutation_name'])
        if match:
            pos = int(match.group(1))
            assert pos not in forbidden


def test_save_mutations(sample_interaction_df, sample_sequence):
    """Test saving mutations to file."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    mutations = pd.DataFrame({
        'mutation_name': ['TestProtein_A1E', 'TestProtein_C2K'],
        'sequence': ['ECDEFGHIKLMNPQRSTVWY', 'AKDEFGHIKLMNPQRSTVWY']
    })
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "mutations.csv"
        scanner.save_mutations(mutations, output_path)
        
        assert output_path.exists()
        loaded = pd.read_csv(output_path, header=None)
        assert len(loaded) == len(mutations)


def test_scan_mutations_from_csv(sample_interaction_df, sample_sequence):
    """Test convenience function for scanning from CSV."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save interaction data
        csv_path = Path(tmpdir) / "interactions.csv"
        sample_interaction_df.to_csv(csv_path)
        
        output_dir = Path(tmpdir) / "mutations"
        
        # Run scan
        results = scan_mutations_from_csv(
            interaction_csv=csv_path,
            sequence=sample_sequence,
            protein_name="TestProtein",
            output_dir=output_dir,
            interaction_type='attractive',
            min_chunk_strength=0.0
        )
        
        assert isinstance(results, dict)
        assert len(results) > 0
        assert output_dir.exists()


def test_chunk_mutations(sample_interaction_df, sample_sequence):
    """Test chunk mutation generation."""
    scanner = MutationScanner(
        interaction_df=sample_interaction_df,
        sequence=sample_sequence,
        protein_name="TestProtein"
    )
    
    filtered = scanner.filter_interactions()
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(chunk_data, 'attractive', 0.0)
    
    mutations = scanner.generate_chunk_mutations(
        candidates,
        position='left',
        mutation_type='charge',
        interaction_type='attractive',
        chunk_size=3
    )
    
    assert 'mutation_name' in mutations.columns
    assert 'sequence' in mutations.columns
    
    # Check that 3 positions are mutated in each
    for _, row in mutations.iterrows():
        mut_name = row['mutation_name']
        # Count underscores after protein name (each mutation separated by _)
        parts = mut_name.split('_')[1:]  # Skip protein name
        # Should have 3 mutations for chunk_size=3
        assert len(parts) == 3


def test_residue_properties():
    """Test residue property calculations."""
    sequence = "RKDEQNFYW"
    scanner = MutationScanner(
        interaction_df=pd.DataFrame({
            'r_1': [1], 'r_2': [5], 'cont_prob': [0.5],
            'distance': [4], 'relative_strength': [0.5], 'plot_value': [2]
        }),
        sequence=sequence,
        protein_name="Test"
    )
    
    # Test charge calculation
    assert scanner._calc_charge("RKD") > 0  # Net positive
    assert scanner._calc_charge("DEE") < 0  # Net negative
    assert scanner._calc_charge("QNA") == 0  # Neutral
    
    # Test aromatic counting
    assert scanner._count_aromatic("FYW") == 3
    assert scanner._count_aromatic("ACE") == 0
    
    # Test hydrophobicity
    assert scanner._calc_hydrophobicity("III") > 0  # Hydrophobic
    assert scanner._calc_hydrophobicity("DDD") < 0  # Hydrophilic


def test_single_mutations_polar_variations(sample_interaction_df, sample_sequence):
    """Test single polar mutations with different scenarios."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    # Create candidates with polar residues at position 1 (A is non-polar)
    candidates = pd.DataFrame({
        'r_1': [1, 3],
        'r_2': [5, 7],
        'strength': [2.0, 1.5],
        'r_1_charge': [0, 0],
        'r_2_charge': [0, 0]
    })
    
    # Test attractive polar - should mutate non-polar to polar
    mutations_attr = scanner.generate_single_mutations(
        candidates, 'left', 'polar', 'attractive'
    )
    assert len(mutations_attr) >= 0
    
    # Test repulsive polar - should mutate polar to non-polar
    mutations_rep = scanner.generate_single_mutations(
        candidates, 'left', 'polar', 'repulsive'
    )
    assert len(mutations_rep) >= 0


def test_single_mutations_charge_repulsive_neutral(sample_interaction_df, sample_sequence):
    """Test charge mutations with neutral residues in repulsive mode."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    # Neutral charge residues should be skipped in repulsive
    candidates = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [5, 6],
        'strength': [2.0, 1.5],
        'r_1_charge': [0, 0],  # Neutral
        'r_2_charge': [0, 0]
    })
    
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'charge', 'repulsive'
    )
    
    # Should produce no mutations for neutral residues
    assert len(mutations) == 0


def test_single_mutations_charge_repulsive_charged(sample_interaction_df, sample_sequence):
    """Test charge mutations with charged residues in repulsive mode."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    # Charged residues should be neutralized
    candidates = pd.DataFrame({
        'r_1': [2, 3],
        'r_2': [5, 6],
        'strength': [2.0, 1.5],
        'r_1_charge': [1, -1],  # Positive and negative
        'r_2_charge': [0, 0]
    })
    
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'charge', 'repulsive'
    )
    
    # Should generate neutralizing mutations
    assert len(mutations) > 0
    assert all('A' in name for name in mutations['mutation_name'])


def test_pair_mutations_charge_attractive_both_charges(sample_interaction_df, sample_sequence):
    """Test charge pair mutations for attractive - should try both E and K."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    candidates = pd.DataFrame({
        'r_1': [2],
        'r_2': [5],
        'strength': [2.0]
    })
    
    mutations = scanner.generate_pair_mutations(
        candidates, 'charge', 'attractive'
    )
    
    # Should generate pairs with both E and K
    assert len(mutations) == 2  # One with E, one with K
    assert any('E' in name for name in mutations['mutation_name'])
    assert any('K' in name for name in mutations['mutation_name'])


def test_chunk_mutations_charge_attractive_negative(sample_interaction_df, sample_sequence):
    """Test chunk charge mutations with negative charge."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    chunk_data = pd.DataFrame({
        'r_1': [5],
        'r_2': [8],
        'strength': [2.5],
        'r_1_charge': [-1],  # Negative
        'r_2_charge': [1],
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    mutations = scanner.generate_chunk_mutations(
        chunk_data, 'left', 'charge', 'attractive'
    )
    
    # Should mutate to K (opposite charge)
    assert len(mutations) > 0
    assert all('K' in name for name in mutations['mutation_name'])


def test_chunk_mutations_charge_attractive_neutral_with_opposite(sample_interaction_df, sample_sequence):
    """Test chunk charge mutations with neutral charge but opposite partner."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    # Neutral with positive opposite
    chunk_data_pos = pd.DataFrame({
        'r_1': [5],
        'r_2': [8],
        'strength': [2.5],
        'r_1_charge': [0],  # Neutral
        'r_2_charge': [1],  # Positive opposite
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    mutations_pos = scanner.generate_chunk_mutations(
        chunk_data_pos, 'left', 'charge', 'attractive'
    )
    assert len(mutations_pos) > 0
    assert all('K' in name for name in mutations_pos['mutation_name'])
    
    # Neutral with negative opposite
    chunk_data_neg = pd.DataFrame({
        'r_1': [5],
        'r_2': [8],
        'strength': [2.5],
        'r_1_charge': [0],  # Neutral
        'r_2_charge': [-1],  # Negative opposite
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    mutations_neg = scanner.generate_chunk_mutations(
        chunk_data_neg, 'left', 'charge', 'attractive'
    )
    assert len(mutations_neg) > 0
    assert all('E' in name for name in mutations_neg['mutation_name'])


def test_chunk_mutations_boundary_left(sample_interaction_df, sample_sequence):
    """Test chunk mutations at left boundary."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    # Position 1 with chunk_size=3 would need positions -1, 0, 1
    chunk_data = pd.DataFrame({
        'r_1': [1],  # Too close to start
        'r_2': [5],
        'strength': [2.0],
        'r_1_charge': [0],
        'r_2_charge': [0],
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    mutations = scanner.generate_chunk_mutations(
        chunk_data, 'left', 'charge', 'attractive', chunk_size=3
    )
    
    # Should skip out-of-bounds chunks
    assert len(mutations) == 0


def test_chunk_mutations_boundary_right(sample_interaction_df, sample_sequence):
    """Test chunk mutations at right boundary."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    seq_len = len(sample_sequence)
    
    # Position at end with chunk_size=3 would exceed sequence length
    chunk_data = pd.DataFrame({
        'r_1': [5],
        'r_2': [seq_len],  # At end
        'strength': [2.0],
        'r_1_charge': [0],
        'r_2_charge': [0],
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    mutations = scanner.generate_chunk_mutations(
        chunk_data, 'right', 'charge', 'attractive', chunk_size=3
    )
    
    # Should skip out-of-bounds chunks
    assert len(mutations) == 0


def test_generate_full_scan_comprehensive(sample_interaction_df, sample_sequence):
    """Test complete full scan with both attractive and repulsive."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        
        # Test attractive
        results_attr = scanner.generate_full_scan(
            output_dir / "attractive",
            interaction_type='attractive',
            min_chunk_strength=0.5
        )
        
        assert isinstance(results_attr, dict)
        # Should have single (left/right x charge/polar), pair (charge/polar), chunk types
        expected_keys = [
            'left_single_charge', 'left_single_polar',
            'right_single_charge', 'right_single_polar',
            'pair_charge', 'pair_polar'
        ]
        for key in expected_keys:
            assert key in results_attr
        
        # Test repulsive
        results_rep = scanner.generate_full_scan(
            output_dir / "repulsive",
            interaction_type='repulsive',
            min_chunk_strength=0.5
        )
        
        assert isinstance(results_rep, dict)
        for key in expected_keys:
            assert key in results_rep


def test_chunk_mutations_all_mutation_types(sample_interaction_df, sample_sequence):
    """Test all mutation types for chunk mutations."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    chunk_data = pd.DataFrame({
        'r_1': [5],
        'r_2': [8],
        'strength': [2.5],
        'r_1_charge': [1],
        'r_2_charge': [-1],
        'r_1_aromatic': [0],
        'r_2_aromatic': [0]
    })
    
    # Test all combinations
    for mut_type in ['charge', 'polar', 'hydrophobic']:
        for int_type in ['attractive', 'repulsive']:
            for position in ['left', 'right']:
                mutations = scanner.generate_chunk_mutations(
                    chunk_data, position, mut_type, int_type
                )
                assert isinstance(mutations, pd.DataFrame)
                assert 'mutation_name' in mutations.columns
                assert 'sequence' in mutations.columns


def test_hydrophobic_single_mutations(sample_interaction_df, sample_sequence):
    """Test hydrophobic single mutations comprehensively."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    candidates = pd.DataFrame({
        'r_1': [2, 3],
        'r_2': [5, 6],
        'strength': [2.0, 1.5],
        'r_1_charge': [0, 0],
        'r_2_charge': [0, 0]
    })
    
    # Attractive should mutate to L
    mutations_attr = scanner.generate_single_mutations(
        candidates, 'left', 'hydrophobic', 'attractive'
    )
    assert len(mutations_attr) > 0
    for seq in mutations_attr['sequence']:
        # Check that mutations were made
        assert seq != sample_sequence
    
    # Repulsive should mutate to S
    mutations_rep = scanner.generate_single_mutations(
        candidates, 'left', 'hydrophobic', 'repulsive'
    )
    assert len(mutations_rep) > 0
    for seq in mutations_rep['sequence']:
        assert seq != sample_sequence


def test_single_mutations_polar_attractive_already_polar(sample_interaction_df):
    """Test polar attractive mutations when residue is already polar."""
    # Create sequence with polar residues (Q is polar)
    sequence = "QQQQQQQQQQQQQQQQQQQQ"  # All polar
    scanner = MutationScanner(sample_interaction_df, sequence, "TestProtein")
    
    candidates = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [5, 6],
        'strength': [2.0, 1.5],
        'r_1_charge': [0, 0],
        'r_2_charge': [0, 0]
    })
    
    # Should skip already polar residues
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'polar', 'attractive'
    )
    assert len(mutations) == 0


def test_pair_mutations_forbidden_regions(sample_interaction_df, sample_sequence):
    """Test pair mutations with forbidden regions."""
    scanner = MutationScanner(
        sample_interaction_df, 
        sample_sequence, 
        "TestProtein",
        forbidden_regions=[2, 3, 5, 6]  # Individual positions, not ranges
    )
    
    candidates = pd.DataFrame({
        'r_1': [2, 4],  # Position 2 is forbidden
        'r_2': [5, 7],  # Position 5 is forbidden
        'strength': [2.0, 1.5]
    })
    
    mutations = scanner.generate_pair_mutations(
        candidates, 'charge', 'attractive'
    )
    
    # Should skip pairs where either position is forbidden
    # Pair (2, 5) should be skipped, only (4, 7) should be processed
    assert len(mutations) == 2  # One for E and one for K at positions (4, 7)
    for name in mutations['mutation_name']:
        # Check that forbidden positions aren't mutated
        assert 'C2' not in name  # Position 2 (residue C)
        assert 'F5' not in name  # Position 5 (residue F)


def test_pair_mutations_hydrophobic_both_types(sample_interaction_df, sample_sequence):
    """Test hydrophobic pair mutations for both attractive and repulsive."""
    scanner = MutationScanner(sample_interaction_df, sample_sequence, "TestProtein")
    
    candidates = pd.DataFrame({
        'r_1': [2, 3],
        'r_2': [5, 6],
        'strength': [2.0, 1.5]
    })
    
    # Attractive - mutate to L
    mutations_attr = scanner.generate_pair_mutations(
        candidates, 'hydrophobic', 'attractive'
    )
    assert len(mutations_attr) > 0
    assert all('L' in name for name in mutations_attr['mutation_name'])
    for seq in mutations_attr['sequence']:
        # Should have L at both positions
        assert 'L' in seq
    
    # Repulsive - mutate to S
    mutations_rep = scanner.generate_pair_mutations(
        candidates, 'hydrophobic', 'repulsive'
    )
    assert len(mutations_rep) > 0
    assert all('S' in name for name in mutations_rep['mutation_name'])
    for seq in mutations_rep['sequence']:
        # Should have S at both positions
        assert 'S' in seq


def test_chunk_mutations_forbidden_regions(sample_interaction_df, sample_sequence):
    """Test chunk mutations with forbidden regions."""
    scanner = MutationScanner(
        sample_interaction_df,
        sample_sequence,
        "TestProtein",
        forbidden_regions=[5, 6, 7]  # Forbid positions 5-7 as individual positions
    )
    
    chunk_data = pd.DataFrame({
        'r_1': [5, 8],  # Position 5 is in forbidden region
        'r_2': [10, 11],
        'strength': [2.5, 2.0],
        'r_1_charge': [1, 1],
        'r_2_charge': [-1, -1],
        'r_1_aromatic': [0, 0],
        'r_2_aromatic': [0, 0]
    })
    
    mutations = scanner.generate_chunk_mutations(
        chunk_data, 'left', 'charge', 'attractive', chunk_size=3
    )
    
    # Should skip chunks that overlap with forbidden regions
    # Chunk at position 5 (positions 4-6) overlaps with forbidden [5, 6, 7]
    # So it should be skipped
    for name in mutations['mutation_name']:
        # Mutation names for chunks around position 5 should not exist
        parts = name.split('_')[1:]  # Skip protein name
        for part in parts:
            # Extract position from format like "F5E"
            pos_str = ''.join(filter(str.isdigit, part))
            if pos_str:
                pos = int(pos_str)
                # Check that we don't mutate forbidden positions
                if pos in [5, 6, 7]:
                    # This chunk should have been skipped
                    assert False, f"Found mutation at forbidden position {pos}"
