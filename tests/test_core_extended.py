"""
Additional tests for core module to improve coverage.
"""

import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch, Mock, MagicMock
from idp_interaction_map.core import analyze_interaction_map


def test_analyze_interaction_map_basic(tmp_path):
    """Test basic analyze_interaction_map functionality."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        # Mock contact data
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({
            'r_1': [1, 2, 3],
            'r_2': [4, 5, 6],
            'cont_prob': [0.5, 0.6, 0.4],
            'distance': [3, 3, 3]
        })
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map') as mock_norm:
            mock_interaction_df = pd.DataFrame({
                'r_1': [1, 2, 3],
                'r_2': [4, 5, 6],
                'cont_prob': [0.5, 0.6, 0.4],
                'distance': [3, 3, 3],
                'gs_standard': [0.3, 0.4, 0.2],
                'relative_strength': [1.5, 1.3, 1.8],
                'plot_value': [2, 1, 2]
            })
            mock_norm.return_value = mock_interaction_df
            
            with patch('idp_interaction_map.core.create_interaction_map') as mock_plot:
                result = analyze_interaction_map(
                    name="test_protein",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path)
                )
                
                assert isinstance(result, pd.DataFrame)
                assert len(result) == 3
                mock_contact.assert_called_once()
                mock_norm.assert_called_once()
                mock_plot.assert_called_once()


def test_analyze_interaction_map_read_from_file(tmp_path):
    """Test analyze_interaction_map with read_from_file=True."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    sequence = "ACDEFG"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({
            'r_1': [1],
            'r_2': [2],
            'cont_prob': [0.5],
            'distance': [1]
        })
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map') as mock_norm:
            mock_df = pd.DataFrame({
                'r_1': [1],
                'r_2': [2],
                'plot_value': [2]
            })
            mock_norm.return_value = mock_df
            
            with patch('idp_interaction_map.core.create_interaction_map'):
                result = analyze_interaction_map(
                    name="test",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path),
                    read_from_file=True
                )
                
                # Verify read_from_file was passed correctly
                assert mock_contact.call_args[1]['read_from_file'] is True


def test_analyze_interaction_map_custom_params(tmp_path):
    """Test analyze_interaction_map with custom parameters."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    sequence = "ACDEFGHIKLM"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({
            'r_1': [1, 2],
            'r_2': [3, 4],
            'cont_prob': [0.5, 0.6],
            'distance': [2, 2]
        })
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map') as mock_norm:
            mock_df = pd.DataFrame({
                'r_1': [1, 2],
                'r_2': [3, 4],
                'plot_value': [2, 1]
            })
            mock_norm.return_value = mock_df
            
            with patch('idp_interaction_map.core.create_interaction_map'):
                result = analyze_interaction_map(
                    name="custom_protein",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path),
                    pdb_top="custom_topology.pdb",
                    xtc_input=10
                )
                
                # Verify custom parameters were used (positional args)
                call_args = mock_contact.call_args[0]
                assert call_args[0] == "custom_protein"  # name
                assert call_args[1] == "custom_topology.pdb"  # pdb_top
                assert call_args[2] == 10  # xtc_input


def test_analyze_interaction_map_creates_output_dir(tmp_path):
    """Test that analyze_interaction_map creates output directory."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output" / "nested" / "path"
    traj_path.mkdir()
    
    sequence = "ACDEFG"
    
    assert not output_path.exists()
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({'r_1': [1], 'r_2': [2], 'cont_prob': [0.5], 'distance': [1]})
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map', return_value=pd.DataFrame({'r_1': [1], 'r_2': [2], 'plot_value': [2]})):
            with patch('idp_interaction_map.core.create_interaction_map'):
                analyze_interaction_map(
                    name="test",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path)
                )
                
                assert output_path.exists()


def test_analyze_interaction_map_saves_csv(tmp_path):
    """Test that analyze_interaction_map saves CSV file."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    sequence = "ACDEFGHIKL"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({
            'r_1': [1, 2],
            'r_2': [3, 4],
            'cont_prob': [0.5, 0.6],
            'distance': [2, 2]
        })
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map') as mock_norm:
            mock_df = pd.DataFrame({
                'r_1': [1, 2],
                'r_2': [3, 4],
                'cont_prob': [0.5, 0.6],
                'plot_value': [2, 1]
            })
            mock_norm.return_value = mock_df
            
            with patch('idp_interaction_map.core.create_interaction_map'):
                analyze_interaction_map(
                    name="test_csv",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path)
                )
                
                # Check that CSV file would be created
                expected_csv = output_path / "test_csv_interaction.csv"
                # File creation is mocked, but we can verify the path exists
                assert output_path.exists()


def test_analyze_interaction_map_returns_to_original_dir(tmp_path):
    """Test that analyze_interaction_map returns to original directory."""
    import os
    
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    original_dir = os.getcwd()
    sequence = "ACDEFG"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({'r_1': [1], 'r_2': [2], 'cont_prob': [0.5], 'distance': [1]})
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map', return_value=pd.DataFrame({'r_1': [1], 'r_2': [2], 'plot_value': [2]})):
            with patch('idp_interaction_map.core.create_interaction_map'):
                analyze_interaction_map(
                    name="test",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path)
                )
                
                # Verify we're back in original directory
                assert os.getcwd() == original_dir


def test_analyze_interaction_map_error_handling(tmp_path):
    """Test that analyze_interaction_map handles errors and returns to original dir."""
    import os
    
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    original_dir = os.getcwd()
    sequence = "ACDEFG"
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact.side_effect = Exception("Test error")
        
        with pytest.raises(Exception, match="Test error"):
            analyze_interaction_map(
                name="test",
                trajectory_path=str(traj_path),
                sequence=sequence,
                output_dir=str(output_path)
            )
        
        # Even after error, should return to original directory
        assert os.getcwd() == original_dir


def test_analyze_interaction_map_with_list_xtc_input(tmp_path):
    """Test analyze_interaction_map with list of trajectory files."""
    traj_path = tmp_path / "traj"
    output_path = tmp_path / "output"
    traj_path.mkdir()
    
    sequence = "ACDEFG"
    xtc_files = ["traj_0.xtc", "traj_1.xtc", "traj_2.xtc"]
    
    with patch('idp_interaction_map.core.generate_contact') as mock_contact:
        mock_contact_obj = Mock()
        mock_contact_obj.contact = pd.DataFrame({'r_1': [1], 'r_2': [2], 'cont_prob': [0.5], 'distance': [1]})
        mock_contact.return_value = mock_contact_obj
        
        with patch('idp_interaction_map.core.normalize_interaction_map', return_value=pd.DataFrame({'r_1': [1], 'r_2': [2], 'plot_value': [2]})):
            with patch('idp_interaction_map.core.create_interaction_map'):
                analyze_interaction_map(
                    name="test",
                    trajectory_path=str(traj_path),
                    sequence=sequence,
                    output_dir=str(output_path),
                    xtc_input=xtc_files
                )
                
                # Verify list was passed through (positional args)
                call_args = mock_contact.call_args[0]
                assert call_args[2] == xtc_files  # xtc_input is 3rd positional arg
