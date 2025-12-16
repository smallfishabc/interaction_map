"""
Additional tests for contact_map module to improve coverage.
"""

import pytest
import numpy as np
import mdtraj as md
import pandas as pd
from pathlib import Path
from unittest.mock import Mock, patch
from idp_interaction_map.contact_map import (
    ContactProbData,
    load_xtc,
    load_traj_protein,
    generate_contact,
)


def test_contact_prob_data_read_from_file(tmp_path):
    """Test loading ContactProbData from file."""
    # Create a mock contact CSV file
    contact_file = tmp_path / "test_1.2_contact_df_1201.csv"
    df = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [2, 3],
        'cont_prob': [0.5, 0.3],
        'distance': [1, 1]
    })
    df.to_csv(contact_file)
    
    # Change to temp directory
    import os
    old_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        
        # Load from file
        contact_data = ContactProbData("test", cutoff=1.2, read_from_file=True)
        
        assert contact_data.name == "test"
        assert contact_data.cutoff == 1.2
        assert len(contact_data.contact) == 2
        assert 'r_1' in contact_data.contact.columns
        
    finally:
        os.chdir(old_cwd)


def test_contact_prob_data_with_trajectory():
    """Test ContactProbData with actual trajectory."""
    # Create a minimal mock trajectory
    mock_traj = Mock(spec=md.Trajectory)
    mock_top = Mock()
    
    # Mock select_pairs to return some atom pairs
    mock_top.select_pairs.return_value = np.array([[0, 1], [0, 2]])
    mock_traj.top = mock_top
    
    # Mock compute_contacts to return distances
    with patch('idp_interaction_map.contact_map.md.compute_contacts') as mock_compute:
        # Return distances for 2 pairs across 10 frames
        # Shape should be (n_frames, n_pairs) = (10, 2)
        distances = np.array([[0.5, 0.7], [0.6, 0.8], [0.5, 0.7], [0.6, 0.8], [0.5, 0.7],
                             [0.6, 0.8], [0.5, 0.7], [0.6, 0.8], [0.5, 0.7], [0.6, 0.8]])
        pairs = np.array([[0, 1], [0, 2]])
        mock_compute.return_value = (distances, pairs)
        
        contact_data = ContactProbData("test", cutoff=0.6, traj=mock_traj)
        
        assert contact_data.name == "test"
        assert len(contact_data.contact) > 0
        assert 'cont_prob' in contact_data.contact.columns


def test_contact_prob_data_invalid():
    """Test ContactProbData with invalid arguments."""
    with pytest.raises(ValueError, match="Must provide either trajectory"):
        ContactProbData("test", cutoff=1.2)


def test_load_xtc_single_file(tmp_path):
    """Test loading single XTC file."""
    pdb_file = tmp_path / "test.pdb"
    xtc_file = tmp_path / "test.xtc"
    
    pdb_content = """ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00
END
"""
    pdb_file.write_text(pdb_content)
    
    with patch('idp_interaction_map.contact_map.md.load') as mock_load:
        mock_traj = Mock(spec=md.Trajectory)
        mock_load.return_value = mock_traj
        
        result = load_xtc(str(xtc_file), str(pdb_file))
        
        assert result == mock_traj


def test_load_xtc_multiple_files(tmp_path):
    """Test loading multiple XTC files."""
    pdb_file = tmp_path / "test.pdb"
    xtc_files = [tmp_path / f"test_{i}.xtc" for i in range(2)]
    
    pdb_content = """ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00
END
"""
    pdb_file.write_text(pdb_content)
    
    with patch('idp_interaction_map.contact_map.md.load') as mock_load:
        combined_traj = Mock(spec=md.Trajectory)
        mock_load.return_value = combined_traj
        
        result = load_xtc([str(f) for f in xtc_files], str(pdb_file))
        
        # Verify md.load was called with the list of files
        mock_load.assert_called_once_with([str(f) for f in xtc_files], top=str(pdb_file))
        assert result == combined_traj


def test_load_traj_protein_xtc_int(tmp_path):
    """Test load_traj_protein with integer (number of XTC files)."""
    pdb_file = "test.pdb"
    
    with patch('idp_interaction_map.contact_map.Path') as mock_path_class:
        mock_path = Mock()
        mock_path_class.return_value = mock_path
        
        def glob_side_effect(pattern):
            if "*.xtc" in pattern:
                return [Path(f"traj_{i}.xtc") for i in range(2)]
            return []
        
        mock_path.glob = Mock(side_effect=glob_side_effect)
        
        with patch('idp_interaction_map.contact_map.load_xtc') as mock_load_xtc:
            mock_traj = Mock(spec=md.Trajectory)
            mock_load_xtc.return_value = mock_traj
            
            result = load_traj_protein(pdb_file, 2)
            
            assert result == mock_traj


def test_load_traj_protein_list_xtc():
    """Test load_traj_protein with list of XTC files."""
    pdb_file = "test.pdb"
    xtc_files = ["traj_0.xtc", "traj_1.xtc"]
    
    with patch('idp_interaction_map.contact_map.load_xtc') as mock_load:
        mock_traj = Mock(spec=md.Trajectory)
        mock_load.return_value = mock_traj
        
        result = load_traj_protein(pdb_file, xtc_files)
        
        assert result == mock_traj
        mock_load.assert_called_once_with(xtc_files, pdb_file)


def test_generate_contact_read_from_file(tmp_path):
    """Test generate_contact with read_from_file=True."""
    contact_file = tmp_path / "protein_1.2_contact_df_1201.csv"
    df = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [2, 3],
        'cont_prob': [0.5, 0.3],
        'distance': [1, 1]
    })
    df.to_csv(contact_file)
    
    import os
    old_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        
        result = generate_contact("protein", read_from_file=True)
        
        assert result.name == "protein"
        assert len(result.contact) == 2
        
    finally:
        os.chdir(old_cwd)


def test_generate_contact_with_trajectory():
    """Test generate_contact with trajectory data."""
    with patch('idp_interaction_map.contact_map.load_traj_protein') as mock_load:
        mock_traj = Mock(spec=md.Trajectory)
        mock_top = Mock()
        mock_top.select_pairs.return_value = np.array([[0, 1]])
        mock_traj.top = mock_top
        mock_load.return_value = mock_traj
        
        with patch('idp_interaction_map.contact_map.md.compute_contacts') as mock_compute:
            # Shape: (n_frames, n_pairs) = (10, 1)
            distances = np.array([[0.5], [0.7], [0.5], [0.7], [0.5], [0.7], [0.5], [0.7], [0.5], [0.7]])
            pairs = np.array([[0, 1]])
            mock_compute.return_value = (distances, pairs)
            
            result = generate_contact("test", "test.pdb", 1, cutoff=0.6)
            
            assert result.name == "test"
            assert result.cutoff == 0.6
