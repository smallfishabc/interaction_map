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
