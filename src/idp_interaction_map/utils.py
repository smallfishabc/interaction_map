"""Utility functions for file I/O and sequence processing."""

import logging
import re
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


def read_sequence(sequence_path: Union[str, Path]) -> str:
    """
    Read and clean protein sequence from file.

    Args:
        sequence_path: Path to sequence file (FASTA or plain text)

    Returns:
        Cleaned protein sequence string

    Raises:
        FileNotFoundError: If sequence file doesn't exist
    """
    path = Path(sequence_path)

    if not path.exists():
        raise FileNotFoundError(f"Sequence file not found: {path}")

    logger.info(f"Reading sequence from {path}")

    with open(path, "r") as f:
        seq = f.read()

    # Remove whitespace from end
    while seq and re.search(r"\s", seq[-1]) is not None:
        seq = seq[:-1]

    # Remove FASTA header if present
    if seq.startswith(">"):
        seq = "".join(seq.split("\n")[1:])

    logger.info(f"Read sequence of length {len(seq)}")
    return seq


def read_sequence_from_fasta(trajectory_path: Union[str, Path]) -> str:
    """
    Read sequence from seq.fasta file in trajectory directory.

    Args:
        trajectory_path: Path to trajectory directory

    Returns:
        Protein sequence string
    """
    path = Path(trajectory_path) / "seq.fasta"
    return read_sequence(path)


def read_sequence_from_txt(data_path: Union[str, Path]) -> str:
    """
    Read sequence from seq.txt file.

    Args:
        data_path: Path to data directory

    Returns:
        Protein sequence string
    """
    path = Path(data_path) / "seq.txt"
    return read_sequence(path)
