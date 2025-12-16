"""Tests for utility functions."""

import pytest
from pathlib import Path

from idp_interaction_map.utils import (
    read_sequence,
    read_sequence_from_txt,
    read_sequence_from_fasta,
)


def test_read_sequence_from_txt(temp_dir, sample_sequence):
    """Test reading sequence from text file."""
    seq_file = temp_dir / "seq.txt"
    seq_file.write_text(sample_sequence)

    result = read_sequence(seq_file)
    assert result == sample_sequence


def test_read_sequence_with_whitespace(temp_dir):
    """Test reading sequence with trailing whitespace."""
    seq_file = temp_dir / "seq.txt"
    seq_file.write_text("MDEYK  \n\n")

    result = read_sequence(seq_file)
    assert result == "MDEYK"


def test_read_sequence_fasta_format(temp_dir):
    """Test reading sequence from FASTA format."""
    seq_file = temp_dir / "seq.fasta"
    seq_file.write_text(">test_protein\nMDEYKLPPFG\nESYRGDERKR\n")

    result = read_sequence(seq_file)
    assert result == "MDEYKLPPFGESYRGDERKR"


def test_read_sequence_file_not_found(temp_dir):
    """Test error handling for missing file."""
    with pytest.raises(FileNotFoundError):
        read_sequence(temp_dir / "nonexistent.txt")


def test_read_sequence_from_directory(temp_dir, sample_sequence):
    """Test reading seq.txt from directory."""
    seq_file = temp_dir / "seq.txt"
    seq_file.write_text(sample_sequence)

    result = read_sequence_from_txt(temp_dir)
    assert result == sample_sequence


def test_read_sequence_from_fasta_directory(temp_dir, sample_sequence):
    """Test reading seq.fasta from directory."""
    seq_file = temp_dir / "seq.fasta"
    seq_file.write_text(f">protein\n{sample_sequence}")

    result = read_sequence_from_fasta(temp_dir)
    # Remove FASTA header
    assert result == sample_sequence
