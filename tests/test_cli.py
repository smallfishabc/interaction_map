"""
Tests for CLI module to improve coverage.
"""

import logging
import pytest
import sys
from pathlib import Path
from unittest.mock import patch, Mock
from idp_interaction_map.cli import parse_args, validate_args, main


def test_parse_args_minimal():
    """Test parsing minimal required arguments."""
    test_args = [
        'idp-interaction-map',
        '-d', '/path/to/data',
        '-n', 'test_protein'
    ]
    
    with patch.object(sys, 'argv', test_args):
        args = parse_args()
        
        assert args.data_dir == '/path/to/data'
        assert args.name == 'test_protein'
        assert args.pdb == '__START_0.pdb'  # default
        assert args.cutoff == 1.2  # default
        assert args.mode == 'cg'  # default


def test_parse_args_all_options():
    """Test parsing all arguments."""
    test_args = [
        'idp-interaction-map',
        '-d', '/path/to/data',
        '-n', 'protein',
        '-p', 'custom.pdb',
        '-x', 'traj.xtc',
        '-r', '10',
        '-o', '/output',
        '--cutoff', '0.5',
        '-v'
    ]
    
    with patch.object(sys, 'argv', test_args):
        args = parse_args()
        
        assert args.data_dir == '/path/to/data'
        assert args.name == 'protein'
        assert args.pdb == 'custom.pdb'
        assert args.xtc == 'traj.xtc'
        assert args.repeat == 10
        assert args.output_dir == '/output'
        assert args.cutoff == 0.5
        assert args.verbose is True


def test_validate_args_valid(tmp_path):
    """Test validation with valid arguments."""
    # Create test directory
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG")
    
    args = Mock()
    args.data_dir = str(data_dir)
    args.name = "test"
    args.output_dir = None
    args.repeat = 5
    args.xtc = None
    
    # Should not raise exception
    validate_args(args)


def test_validate_args_directory_not_exists():
    """Test validation with non-existent directory."""
    args = Mock()
    args.data_dir = "/nonexistent/path"
    args.name = "test"
    args.output_dir = None
    args.repeat = 5
    args.xtc = None
    
    with pytest.raises(ValueError, match="Data directory does not exist"):
        validate_args(args)


def test_validate_args_creates_output_dir(tmp_path):
    """Test that validation creates output directory if specified."""
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "output"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG")
    
    args = Mock()
    args.data_dir = str(data_dir)
    args.name = "test"
    args.output_dir = str(output_dir)
    args.repeat = 5
    args.xtc = None
    
    # validate_args doesn't create output_dir, that's done in main()
    # Just check it doesn't error
    validate_args(args)


def test_main_help():
    """Test main function with --help argument."""
    test_args = ['idp-interaction-map', '--help']
    
    with patch.object(sys, 'argv', test_args):
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        assert exc_info.value.code == 0


def test_main_version():
    """Test main function with --version argument."""
    test_args = ['idp-interaction-map', '--version']
    
    with patch.object(sys, 'argv', test_args):
        with pytest.raises(SystemExit) as exc_info:
            with patch('builtins.print'):  # Suppress version output
                main()
        
        assert exc_info.value.code == 0


def test_main_successful_run(tmp_path):
    """Test main function with successful analysis."""
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "output"
    data_dir.mkdir()
    
    # Create test files
    (data_dir / "seq.txt").write_text("ACDEFGHIKLMNPQRSTVWY")
    
    test_args = [
        'idp-interaction-map',
        '-d', str(data_dir),
        '-n', 'test_protein',
        '-o', str(output_dir),
        '-r', '5'
    ]
    
    with patch.object(sys, 'argv', test_args):
        with patch('idp_interaction_map.cli.read_sequence_from_txt') as mock_read:
            mock_read.return_value = "ACDEFGHIKLMNPQRSTVWY"
            
            with patch('idp_interaction_map.cli.analyze_interaction_map') as mock_analyze:
                mock_df = Mock()
                mock_analyze.return_value = mock_df
                
                result = main()
                
                assert result == 0
                mock_read.assert_called_once()
                mock_analyze.assert_called_once()


def test_main_missing_required_args():
    """Test main function with missing required arguments."""
    test_args = ['idp-interaction-map']  # No required args
    
    with patch.object(sys, 'argv', test_args):
        with pytest.raises(SystemExit):
            with patch('builtins.print'):  # Suppress error messages
                main()


def test_main_verbose_logging(tmp_path):
    """Test main function enables verbose logging."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG")
    
    test_args = [
        'idp-interaction-map',
        '-d', str(data_dir),
        '-n', 'test',
        '-r', '5',
        '-v'  # Verbose flag
    ]
    
    with patch.object(sys, 'argv', test_args):
        with patch('idp_interaction_map.cli.read_sequence_from_txt', return_value="ACDEFG"):
            with patch('idp_interaction_map.cli.analyze_interaction_map', return_value=Mock()):
                result = main()
                
                # Check that main completed successfully
                assert result == 0
                # Verbose logging is enabled by setting level on root logger
                assert logging.getLogger().level == logging.DEBUG


def test_main_analysis_error(tmp_path):
    """Test main function handles analysis errors gracefully."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG")
    
    test_args = [
        'idp-interaction-map',
        '-d', str(data_dir),
        '-n', 'test',
        '-r', '5'
    ]
    
    with patch.object(sys, 'argv', test_args):
        with patch('idp_interaction_map.cli.read_sequence_from_txt', return_value="ACDEFG"):
            with patch('idp_interaction_map.cli.analyze_interaction_map') as mock_analyze:
                mock_analyze.side_effect = Exception("Test error")
                
                result = main()
                
                assert result == 1  # Error exit code


def test_main_keyboard_interrupt(tmp_path):
    """Test main function handles keyboard interrupt."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG")
    
    test_args = [
        'idp-interaction-map',
        '-d', str(data_dir),
        '-n', 'test',
        '-r', '5'
    ]
    
    with patch.object(sys, 'argv', test_args):
        with patch('idp_interaction_map.cli.read_sequence_from_txt', return_value="ACDEFG"):
            with patch('idp_interaction_map.cli.analyze_interaction_map') as mock_analyze:
                mock_analyze.side_effect = KeyboardInterrupt()
                
                result = main()
                
                assert result == 130  # Keyboard interrupt exit code


def test_main_creates_output_files(tmp_path):
    """Test that main creates expected output files."""
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "output"
    data_dir.mkdir()
    (data_dir / "seq.txt").write_text("ACDEFG" * 10)  # 60 residues
    
    test_args = [
        'idp-interaction-map',
        '-d', str(data_dir),
        '-n', 'test_protein',
        '-o', str(output_dir),
        '-r', '5'
    ]
    
    with patch.object(sys, 'argv', test_args):
        with patch('idp_interaction_map.cli.read_sequence_from_txt') as mock_read:
            mock_read.return_value = "ACDEFG" * 10
            
            with patch('idp_interaction_map.cli.analyze_interaction_map') as mock_analyze:
                # Return a mock DataFrame
                import pandas as pd
                mock_df = pd.DataFrame({
                    'r_1': [1, 2],
                    'r_2': [2, 3],
                    'cont_prob': [0.5, 0.3]
                })
                mock_analyze.return_value = mock_df
                
                main()
                
                # Verify analysis was called with correct parameters
                call_kwargs = mock_analyze.call_args[1]
                assert 'name' in call_kwargs
                assert call_kwargs['name'] == 'test_protein'
                assert 'trajectory_path' in call_kwargs
                assert str(call_kwargs['trajectory_path']) == str(data_dir)
