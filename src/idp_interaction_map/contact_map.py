"""Contact map generation module for molecular dynamics trajectories."""

import logging
from pathlib import Path
from typing import List, Optional, Union

import mdtraj as md
import pandas as pd

logger = logging.getLogger(__name__)


class ContactProbData:
    """
    Contact probability data container.

    Computes and stores contact probabilities between residue pairs from MD trajectories.

    Attributes:
        name: Protein and contact map identifier
        cutoff: Distance cutoff for contacts (in nm)
        contact: DataFrame with columns ['r_1', 'r_2', 'cont_prob', 'distance']
        use_ca: Whether to use CA-only mode (all-atom) or all atoms (CG)
    """

    def __init__(
        self,
        name: str,
        cutoff: float,
        traj: Optional[md.Trajectory] = None,
        read_from_file: bool = False,
        use_ca: bool = False,
    ):
        """
        Initialize ContactProbData.

        Args:
            name: Protein identifier
            cutoff: Contact distance cutoff in nm
            traj: MDTraj trajectory object (required if not reading from file)
            read_from_file: Whether to load existing contact map from file
            use_ca: If True, use CA atoms only (all-atom). If False, use all atoms (CG)
        """
        self.name = name
        self.cutoff = cutoff
        self.use_ca = use_ca

        if read_from_file:
            self.read_contact()
        elif traj is not None:
            self.contact = self.compute_contact(traj, use_ca=use_ca)
            self.update_csv()
        else:
            raise ValueError("Must provide either trajectory or set read_from_file=True")

    def compute_contact(self, traj: md.Trajectory, use_ca: bool = False) -> pd.DataFrame:
        """
        Calculate contact probability from trajectory.

        Args:
            traj: MDTraj trajectory object
            use_ca: If True, use CA atoms only (all-atom mode). If False, use all atoms (CG mode)

        Returns:
            DataFrame with contact probabilities and residue pair information
        """
        logger.info(f"Computing contacts for {self.name} with cutoff {self.cutoff} nm (CA-only: {use_ca})")

        # Calculate pairwise distances
        if use_ca:
            # All-atom mode: use CA scheme and ignore non-protein
            # Note: scheme='CA' returns 1-indexed residue numbers already
            distances, pairs = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
            # pairs are already 1-indexed when using scheme='CA'
            r1_values = pairs[:, 0]
            r2_values = pairs[:, 1]
        else:
            # CG mode: use all atoms (returns 0-indexed atom/residue numbers)
            indices = traj.top.select_pairs("all", "all")
            distances, pairs = md.compute_contacts(traj, indices)
            # Convert to 1-indexed
            r1_values = pairs[:, 0] + 1
            r2_values = pairs[:, 1] + 1

        # Determine contacts based on cutoff
        is_contact = distances < self.cutoff
        contact_probability = is_contact.mean(axis=0)

        # Create DataFrame with contact information
        contact_df = pd.DataFrame(
            {
                "r_1": r1_values,
                "r_2": r2_values,
                "cont_prob": contact_probability,
            }
        )

        contact_df["r_1"] = contact_df["r_1"].astype("int")
        contact_df["r_2"] = contact_df["r_2"].astype("int")
        contact_df["distance"] = contact_df["r_2"] - contact_df["r_1"]

        logger.info(f"Computed {len(contact_df)} residue pairs")
        return contact_df

    def read_contact(self) -> None:
        """Load contact probability data from CSV file."""
        filename = f"{self.name}_{self.cutoff}_contact_df_1201.csv"
        logger.info(f"Reading contact map from {filename}")
        self.contact = pd.read_csv(filename, index_col=0)

    def save_contact(self) -> None:
        """Save contact probability data to CSV file."""
        filename = f"{self.name}_{self.cutoff}_contact_df_1201.csv"
        logger.info(f"Saving contact map to {filename}")
        self.contact.to_csv(filename)

    def update_csv(self) -> None:
        """Save contact data to file."""
        self.save_contact()


def load_xtc(xtc: Union[str, List[str]], pdb: Union[str, Path]) -> md.Trajectory:
    """
    Load XTC trajectory file(s) with PDB topology.

    Args:
        xtc: XTC filename or list of filenames
        pdb: PDB topology filename

    Returns:
        MDTraj trajectory object
    """
    return md.load(xtc, top=str(pdb))


def load_traj_protein(
    pdb_top: str = "__START_0.pdb", 
    xtc_input: Union[int, List[str]] = "__traj_0.xtc",
    use_ca: bool = False
) -> md.Trajectory:
    """
    Load protein trajectory from PDB topology and XTC file(s).

    Args:
        pdb_top: PDB topology filename
        xtc_input: Either number of trajectory files (generates names as __traj_{i}.xtc)
                   or list of XTC filenames
        use_ca: If True, slice trajectory to protein atoms only (all-atom mode)

    Returns:
        MDTraj trajectory object (sliced to protein if use_ca=True)
    """
    if isinstance(xtc_input, int):
        xtc_list = [f"__traj_{i}.xtc" for i in range(xtc_input)]
    else:
        xtc_list = xtc_input if isinstance(xtc_input, list) else [xtc_input]

    traj = load_xtc(xtc_list, pdb_top)
    
    # Slice to protein atoms if in all-atom mode
    if use_ca:
        logger.info("Slicing trajectory to protein atoms only (all-atom mode)")
        protein_atoms = traj.top.select('protein')
        traj = traj.atom_slice(protein_atoms)
    
    return traj
    return load_xtc(xtc_list, pdb_top)


def generate_contact(
    protein_name: str,
    pdb_top: str = "__START_0.pdb",
    xtc_input: Union[int, List[str]] = 5,
    cutoff: float = 1.2,
    read_from_file: bool = False,
    use_ca: bool = False,
) -> ContactProbData:
    """
    Generate or load contact map for a protein.

    Args:
        protein_name: Protein identifier
        pdb_top: PDB topology filename
        xtc_input: Number of trajectory files or list of XTC filenames
        cutoff: Contact distance cutoff in nm
        read_from_file: Whether to load existing contact map
        use_ca: If True, use CA-only mode for all-atom simulations. If False, use CG mode

    Returns:
        ContactProbData object with contact information
    """
    if read_from_file:
        return ContactProbData(protein_name, cutoff, read_from_file=True, use_ca=use_ca)

    traj = load_traj_protein(pdb_top, xtc_input, use_ca=use_ca)
    return ContactProbData(protein_name, cutoff, traj, use_ca=use_ca)
