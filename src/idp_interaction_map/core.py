"""Core analysis functions."""

import logging
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from idp_interaction_map.contact_map import generate_contact
from idp_interaction_map.normalization import normalize_interaction_map
from idp_interaction_map.plotting import create_interaction_map

logger = logging.getLogger(__name__)


def analyze_interaction_map(
    name: str,
    trajectory_path: Union[str, Path],
    sequence: str,
    output_dir: Union[str, Path],
    pdb_top: str = "__START_0.pdb",
    xtc_input: Union[int, list] = 5,
    read_from_file: bool = False,
    use_ca: bool = False,
    norm_a: float = None,
    norm_b: float = None,
) -> pd.DataFrame:
    """
    Main analysis function to generate interaction map.

    This function:
    1. Loads or computes contact map from trajectory
    2. Normalizes interactions against ideal polymer model
    3. Generates and saves visualization
    4. Returns interaction strength data

    Args:
        name: Protein identifier
        trajectory_path: Path to simulation trajectory files
        sequence: Protein sequence string
        output_dir: Directory for output files
        pdb_top: PDB topology filename
        xtc_input: Number of trajectory files or list of filenames
        read_from_file: Whether to load existing contact map
        use_ca: If True, use CA-only mode for all-atom simulations. If False, use CG mode.
        norm_a: Normalization parameter a (auto-selected based on use_ca if None)
        norm_b: Normalization parameter b (auto-selected based on use_ca if None)

    Returns:
        DataFrame with normalized interaction data

    Example:
        >>> # Coarse-grained simulation (default)
        >>> df = analyze_interaction_map(
        ...     name="protein_cg",
        ...     trajectory_path="./data/cg_traj",
        ...     sequence="ACDEFG...",
        ...     output_dir="./output",
        ...     use_ca=False
        ... )
        >>> 
        >>> # All-atom simulation with CA selection
        >>> df = analyze_interaction_map(
        ...     name="protein_aa",
        ...     trajectory_path="./data/aa_traj",
        ...     sequence="ACDEFG...",
        ...     output_dir="./output",
        ...     use_ca=True
        ... )
    """
    # Auto-select normalization parameters based on simulation type
    if norm_a is None:
        norm_a = 1.64 if use_ca else 13.12
    if norm_b is None:
        norm_b = -1.32 if use_ca else -2.32
    
    # Auto-select inter_cutoff based on simulation type
    inter_cutoff = (1.5, 0.5, -1, -2) if use_ca else (2, 1, -1, -2)
    
    logger.info(f"Analysis mode: {'All-atom (CA)' if use_ca else 'Coarse-grained'}")
    logger.info(f"Normalization parameters: a={norm_a}, b={norm_b}")
    logger.info(f"Interaction cutoff: {inter_cutoff}")
    traj_path = Path(trajectory_path).absolute()
    out_path = Path(output_dir).absolute()
    out_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Analyzing {name} from {traj_path}")

    # Change to trajectory directory (needed for relative file paths)
    import os
    original_dir = Path.cwd()
    try:
        os.chdir(traj_path)
        logger.debug(f"Changed to directory: {traj_path}")

        # Generate or load contact map
        if not read_from_file:
            contact_data = generate_contact(name, pdb_top, xtc_input, use_ca=use_ca)
            interaction_df = normalize_interaction_map(
                contact_data.contact, a1=norm_a, b1=norm_b, inter_cutoff=inter_cutoff
            )

            # Save raw interaction data
            output_csv = out_path / f"{name}_interaction.csv"
            interaction_df.to_csv(output_csv)
            logger.info(f"Saved interaction data to {output_csv}")
        else:
            # Load from existing file
            contact_data = generate_contact(name, read_from_file=True, use_ca=use_ca)
            interaction_df = normalize_interaction_map(
                contact_data.contact, a1=norm_a, b1=norm_b, inter_cutoff=inter_cutoff
            )

        # Generate visualization
        output_name = str(out_path / name)
        create_interaction_map(sequence, len(sequence), interaction_df, output_name)

        logger.info(f"Analysis complete for {name}")
        return interaction_df
    
    finally:
        # Always change back to original directory
        os.chdir(original_dir)
        logger.debug(f"Returned to directory: {original_dir}")
