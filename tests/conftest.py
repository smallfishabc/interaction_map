"""Test configuration and shared fixtures."""

import tempfile
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_sequence():
    """Provide sample protein sequence for testing."""
    # Small test sequence with various residue types
    return "MDEYKLPPFGESYRGDERKRFQNVPVDYFLPSDGRPRPIVTPG"


@pytest.fixture
def sample_sequence_file(temp_dir, sample_sequence):
    """Create sample sequence file."""
    seq_file = temp_dir / "seq.txt"
    seq_file.write_text(sample_sequence)
    return seq_file


@pytest.fixture
def sample_contact_data():
    """Provide sample contact probability data."""
    import pandas as pd

    data = {
        "r_1": [1, 1, 2, 2, 3],
        "r_2": [5, 10, 6, 11, 8],
        "cont_prob": [0.8, 0.3, 0.7, 0.2, 0.6],
        "distance": [4, 9, 4, 9, 5],
    }
    return pd.DataFrame(data)


@pytest.fixture
def mock_trajectory(monkeypatch):
    """Mock MDTraj trajectory object."""
    import mdtraj as md

    class MockTopology:
        def select_pairs(self, sel1, sel2):
            # Return indices for 5 atoms
            return np.array([[i, j] for i in range(5) for j in range(i + 1, 5)])

    class MockTrajectory:
        def __init__(self):
            self.top = MockTopology()
            self.n_frames = 100
            self.n_atoms = 5

    return MockTrajectory()
