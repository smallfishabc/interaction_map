"""
Additional tests for plotting module to improve coverage.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, Mock
from idp_interaction_map.plotting import create_interaction_map


def test_create_interaction_map_minimal():
    """Test create_interaction_map with minimal data."""
    sequence = "ACDEFG"
    interaction_df = pd.DataFrame({
        'r_1': [1, 2, 3],
        'r_2': [4, 5, 6],
        'plot_value': [2, -2, 1]
    })
    
    with patch('idp_interaction_map.plotting.create_network') as mock_network:
        import networkx as nx
        # Create a real graph to avoid mock issues
        real_graph = nx.MultiDiGraph()
        for i in range(1, len(sequence)+1):
            real_graph.add_node(i, residue=sequence[i-1], pos=(i, 10))
        mock_network.return_value = real_graph
        
        with patch('idp_interaction_map.plotting.create_sequence_visualization') as mock_viz:
            mock_fig = Mock()
            mock_ax = Mock()
            mock_viz.return_value = (mock_fig, mock_ax)
            
            with patch('idp_interaction_map.plotting.plot_interactions') as mock_plot:
                with patch('idp_interaction_map.plotting.plt.savefig'):
                    with patch('idp_interaction_map.plotting.plt.show'):
                        create_interaction_map(sequence, len(sequence), interaction_df, "test_output")
                    
                        mock_network.assert_called_once()
                        # plot_interactions is called 4 times (for each interaction type)
                        assert mock_plot.call_count == 4


def test_create_interaction_map_filters_interactions():
    """Test that create_interaction_map filters by plot_value."""
    sequence = "ACDEFGHIKLM"
    interaction_df = pd.DataFrame({
        'r_1': [1, 2, 3, 4],
        'r_2': [5, 6, 7, 8],
        'plot_value': [2, 1, 0, -1]  # 0 should be filtered out
    })
    
    with patch('idp_interaction_map.plotting.create_network') as mock_network:
        import networkx as nx
        real_graph = nx.MultiDiGraph()
        for i in range(1, len(sequence)+1):
            real_graph.add_node(i, residue=sequence[i-1], pos=(i, 10))
        mock_network.return_value = real_graph
        
        with patch('idp_interaction_map.plotting.create_sequence_visualization') as mock_viz:
            mock_viz.return_value = (Mock(), Mock())
            
            with patch('idp_interaction_map.plotting.plot_interactions') as mock_plot:
                with patch('idp_interaction_map.plotting.plt.savefig'):
                    with patch('idp_interaction_map.plotting.plt.show'):
                        create_interaction_map(sequence, len(sequence), interaction_df, "test")
                    
                        # plot_interactions is called 4 times for the 4 interaction types
                        assert mock_plot.call_count == 4


def test_create_interaction_map_custom_dpi():
    """Test create_interaction_map saves with correct DPI."""
    sequence = "ACDEFG"
    interaction_df = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [4, 5],
        'plot_value': [2, -2]
    })
    
    with patch('idp_interaction_map.plotting.create_network') as mock_network:
        import networkx as nx
        real_graph = nx.MultiDiGraph()
        for i in range(1, len(sequence)+1):
            real_graph.add_node(i, residue=sequence[i-1], pos=(i, 10))
        mock_network.return_value = real_graph
        
        with patch('idp_interaction_map.plotting.create_sequence_visualization') as mock_viz:
            mock_viz.return_value = (Mock(), Mock())
            
            with patch('idp_interaction_map.plotting.plot_interactions'):
                with patch('idp_interaction_map.plotting.plt.savefig') as mock_save:
                    with patch('idp_interaction_map.plotting.plt.show'):
                        create_interaction_map(
                            sequence, len(sequence), interaction_df, "test"
                        )
                    
                        # Verify savefig was called with dpi=300 (default)
                        assert mock_save.call_count == 2  # PNG and SVG
                        png_call = mock_save.call_args_list[0]
                        assert png_call[1]['dpi'] == 300


def test_create_interaction_map_large_sequence():
    """Test create_interaction_map with large sequence."""
    sequence = "A" * 200  # 200 residues
    interaction_df = pd.DataFrame({
        'r_1': list(range(1, 51)),
        'r_2': list(range(51, 101)),
        'plot_value': [2] * 25 + [-2] * 25
    })
    
    with patch('idp_interaction_map.plotting.create_network') as mock_network:
        import networkx as nx
        real_graph = nx.MultiDiGraph()
        for i in range(1, len(sequence)+1):
            real_graph.add_node(i, residue='A', pos=(i, 10))
        mock_network.return_value = real_graph
        
        with patch('idp_interaction_map.plotting.create_sequence_visualization') as mock_viz:
            mock_viz.return_value = (Mock(), Mock())
            
            with patch('idp_interaction_map.plotting.plot_interactions') as mock_plot:
                with patch('idp_interaction_map.plotting.plt.savefig'):
                    with patch('idp_interaction_map.plotting.plt.show'):
                        create_interaction_map(sequence, len(sequence), interaction_df, "test")
                    
                        mock_network.assert_called_once()
                        # Should handle large sequence without error


def test_create_network_with_interaction_data():
    """Test create_network creates graph with correct nodes."""
    from idp_interaction_map.plotting import create_network
    
    sequence = "ACDEFG"
    length = len(sequence)
    
    G = create_network(sequence, length)
    
    assert G.number_of_nodes() == length
    # Verify nodes have correct attributes
    for i in range(1, length + 1):
        assert i in G.nodes
        assert 'pos' in G.nodes[i]
        assert 'residue' in G.nodes[i]


def test_plot_interactions_creates_figure():
    """Test that plot_interactions plots interactions correctly."""
    from idp_interaction_map.plotting import plot_interactions
    
    interaction_df = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [6, 7],
        'plot_value': [2, -2],
        'cont_prob': [0.5, 0.3],
        'relative_strength': [0.1, 0.2],
        'distance': [5, 6]
    })
    
    layout = {1: (1, 10), 2: (2, 10), 6: (6, 10), 7: (7, 10)}
    mock_ax = Mock()
    
    with patch('idp_interaction_map.plotting.plt.annotate'):
        # Test with interaction_type 2 (strong favorable)
        result = plot_interactions(interaction_df, layout, mock_ax, 2)
        assert result == 1  # Favorable interaction returns 1
        
        # Test with interaction_type -2 (strong unfavorable)
        result = plot_interactions(interaction_df, layout, mock_ax, -2)
        assert result == -1  # Unfavorable interaction returns -1
