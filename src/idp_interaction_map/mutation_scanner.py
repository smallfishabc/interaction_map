"""
Mutation scanning module for identifying and generating mutations based on interaction maps.

This module provides functionality to:
- Analyze interaction maps to identify residue chunks driving interactions
- Generate mutations to enhance attractive or disrupt repulsive interactions
- Create mutation libraries for experimental validation
"""

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import logging

logger = logging.getLogger(__name__)


class MutationScanner:
    """
    Scan interaction maps and generate mutation candidates.
    
    This class analyzes interaction strength between residue pairs and generates
    targeted mutations to modulate protein interactions.
    """
    
    # Amino acid properties
    RESIDUE_TYPES = {
        'F': 'aromatic', 'Y': 'aromatic', 'W': 'aromatic',
        'R': 'positive', 'K': 'positive', 'H': 'positive',
        'D': 'negative', 'E': 'negative',
        'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar',
        'A': 'hydrophobic', 'V': 'hydrophobic', 'I': 'hydrophobic',
        'L': 'hydrophobic', 'M': 'hydrophobic', 'C': 'hydrophobic',
        'G': 'hydrophobic', 'P': 'hydrophobic'
    }
    
    # Kyte-Doolittle hydrophobicity scale
    HYDROPHOBICITY = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    
    def __init__(
        self,
        interaction_df: pd.DataFrame,
        sequence: str,
        protein_name: str,
        forbidden_regions: Optional[List[int]] = None
    ):
        """
        Initialize mutation scanner.
        
        Args:
            interaction_df: DataFrame with interaction map data
            sequence: Protein sequence
            protein_name: Name of the protein
            forbidden_regions: List of residue positions that should not be mutated
        """
        self.df = interaction_df.copy()
        self.sequence = sequence
        self.protein_name = protein_name
        self.forbidden_regions = forbidden_regions or []
        self.seq_length = len(sequence)
        
        # Create sequence reference dictionary
        self.seq_ref_dic = dict(zip(range(1, len(sequence) + 1), list(sequence)))
        
        # Add residue information to dataframe
        self._add_residue_info()
        
        logger.info(f"Initialized MutationScanner for {protein_name} ({self.seq_length} residues)")
    
    def _add_residue_info(self):
        """Add residue type information to interaction dataframe."""
        # Map residue identities
        self.df['r_1_res'] = self.df['r_1'].map(self.seq_ref_dic)
        self.df['r_2_res'] = self.df['r_2'].map(self.seq_ref_dic)
        
        # Map residue types
        self.df['r_1_type'] = self.df['r_1_res'].map(self.RESIDUE_TYPES)
        self.df['r_2_type'] = self.df['r_2_res'].map(self.RESIDUE_TYPES)
    
    def filter_interactions(
        self,
        min_contact_prob: float = 0.01,
        min_distance: int = 4,
        min_relative_strength: float = 0.1,
        exclude_termini: bool = True
    ) -> pd.DataFrame:
        """
        Filter interaction map for significant interactions.
        
        Args:
            min_contact_prob: Minimum contact probability
            min_distance: Minimum sequence separation
            min_relative_strength: Minimum relative interaction strength
            exclude_termini: Whether to exclude terminal residues
        
        Returns:
            Filtered DataFrame
        """
        selected = self.df[
            (self.df['cont_prob'] >= min_contact_prob) &
            (self.df['distance'] >= min_distance) &
            (self.df['relative_strength'] >= min_relative_strength)
        ]
        
        if exclude_termini:
            selected = selected[
                ~(selected['r_1'] == 1) &
                ~(selected['r_2'] == self.seq_length)
            ]
        
        logger.info(f"Filtered to {len(selected)} significant interactions")
        return selected
    
    def calculate_chunk_strength(
        self,
        selected_df: pd.DataFrame,
        chunk_size: int = 3
    ) -> pd.DataFrame:
        """
        Calculate interaction strength for residue chunks.
        
        A chunk is defined as adjacent residues (default: 3 residues).
        Chunk strength accounts for all pairwise interactions between
        residues in two chunks, weighted by distance from chunk center.
        
        Args:
            selected_df: DataFrame with filtered interactions
            chunk_size: Size of residue chunks (default: 3)
        
        Returns:
            DataFrame with chunk information added
        """
        result = selected_df.copy()
        
        # Define chunk indices (centered on each residue)
        offset = chunk_size // 2
        result['r_1_chunk_index'] = result['r_1'].apply(
            lambda x: list(range(max(1, x - offset), min(self.seq_length, x + offset) + 1))
        )
        result['r_2_chunk_index'] = result['r_2'].apply(
            lambda x: list(range(max(1, x - offset), min(self.seq_length, x + offset) + 1))
        )
        
        # Define chunk sequences
        result['r_1_chunk'] = result['r_1'].apply(
            lambda x: self.sequence[max(0, x - 1 - offset):min(self.seq_length, x + offset)]
        )
        result['r_2_chunk'] = result['r_2'].apply(
            lambda x: self.sequence[max(0, x - 1 - offset):min(self.seq_length, x + offset)]
        )
        
        # Calculate chunk strength
        def calc_strength(row):
            left = row['r_1']
            right = row['r_2']
            chunk_strength = 0
            
            for i in row['r_1_chunk_index']:
                for j in row['r_2_chunk_index']:
                    # Weight by distance from chunk center
                    dist = abs(i - left) + abs(j - right)
                    
                    # Find interaction strength
                    match = selected_df[
                        (selected_df['r_1'] == i) &
                        (selected_df['r_2'] == j)
                    ]
                    
                    if not match.empty:
                        strength = match['relative_strength'].values[0]
                        chunk_strength += strength / (2 ** dist)
            
            return chunk_strength
        
        result['chunk_strength'] = result.apply(calc_strength, axis=1)
        
        # Calculate chunk properties
        result['r_1_hydro'] = result['r_1_chunk'].apply(self._calc_hydrophobicity)
        result['r_2_hydro'] = result['r_2_chunk'].apply(self._calc_hydrophobicity)
        result['r_1_charge'] = result['r_1_chunk'].apply(self._calc_charge)
        result['r_2_charge'] = result['r_2_chunk'].apply(self._calc_charge)
        result['r_1_aromatic'] = result['r_1_chunk'].apply(self._count_aromatic)
        result['r_2_aromatic'] = result['r_2_chunk'].apply(self._count_aromatic)
        
        logger.info(f"Calculated chunk strength for {len(result)} interactions")
        return result
    
    def _calc_hydrophobicity(self, chunk: str) -> float:
        """Calculate total hydrophobicity of a chunk."""
        return sum(self.HYDROPHOBICITY.get(aa, 0) for aa in chunk)
    
    def _calc_charge(self, chunk: str) -> int:
        """Calculate net charge of a chunk."""
        charge = 0
        for aa in chunk:
            if aa in ['K', 'R', 'H']:
                charge += 1
            elif aa in ['D', 'E']:
                charge -= 1
        return charge
    
    def _count_aromatic(self, chunk: str) -> int:
        """Count aromatic residues in a chunk."""
        return sum(1 for aa in chunk if aa in ['F', 'Y', 'W'])
    
    def identify_candidates(
        self,
        chunk_df: pd.DataFrame,
        interaction_type: str = 'attractive',
        min_chunk_strength: float = 1.0
    ) -> pd.DataFrame:
        """
        Identify mutation candidates based on interaction type.
        
        Args:
            chunk_df: DataFrame with chunk information
            interaction_type: 'attractive' or 'repulsive'
            min_chunk_strength: Minimum chunk strength threshold
        
        Returns:
            DataFrame with mutation candidates
        """
        if interaction_type == 'attractive':
            # Target favorable interactions (plot_value > 0)
            candidates = chunk_df[chunk_df['plot_value'] > 0]
        else:  # repulsive
            # Target unfavorable interactions (plot_value < 0)
            candidates = chunk_df[chunk_df['plot_value'] < 0]
        
        # Filter by chunk strength
        candidates = candidates[candidates['chunk_strength'] > min_chunk_strength]
        candidates = candidates.sort_values(by='chunk_strength', ascending=False)
        
        logger.info(f"Identified {len(candidates)} {interaction_type} candidates")
        return candidates
    
    def generate_single_mutations(
        self,
        candidates: pd.DataFrame,
        position: str = 'left',
        mutation_type: str = 'charge',
        interaction_type: str = 'attractive'
    ) -> pd.DataFrame:
        """
        Generate single point mutations.
        
        Args:
            candidates: DataFrame with mutation candidates
            position: 'left' (r_1) or 'right' (r_2)
            mutation_type: 'charge', 'polar', or 'hydrophobic'
            interaction_type: 'attractive' or 'repulsive'
        
        Returns:
            DataFrame with mutation sequences
        """
        seq_list = []
        seq_name_list = []
        
        pos_col = f'r_{1 if position == "left" else 2}'
        charge_col = f'{pos_col}_charge'
        
        for residue_idx in candidates[pos_col].unique():
            if residue_idx in self.forbidden_regions:
                logger.debug(f"Skipping forbidden position {residue_idx}")
                continue
            
            wt_res = self.seq_ref_dic[residue_idx]
            mut_seq = list(self.sequence)
            
            if mutation_type == 'charge':
                if interaction_type == 'attractive':
                    # Enhance attraction: opposite charges
                    charge_status = candidates[candidates[pos_col] == residue_idx][charge_col].iloc[0]
                    mut_res = 'E' if charge_status >= 0 else 'K'
                else:  # repulsive
                    # Disrupt repulsion: neutralize or reverse charges
                    charge_status = candidates[candidates[pos_col] == residue_idx][charge_col].iloc[0]
                    if charge_status > 0:
                        mut_res = 'A'  # Neutralize positive
                    elif charge_status < 0:
                        mut_res = 'A'  # Neutralize negative
                    else:
                        continue
            
            elif mutation_type == 'polar':
                if interaction_type == 'attractive':
                    # Enhance polarity
                    if self.RESIDUE_TYPES.get(wt_res) != 'polar':
                        mut_res = 'Q'
                    else:
                        continue
                else:  # repulsive
                    # Reduce polarity
                    if self.RESIDUE_TYPES.get(wt_res) == 'polar':
                        mut_res = 'A'
                    else:
                        continue
            
            elif mutation_type == 'hydrophobic':
                if interaction_type == 'attractive':
                    # Enhance hydrophobicity
                    mut_res = 'L'
                else:  # repulsive
                    # Reduce hydrophobicity
                    mut_res = 'S'
            
            mut_seq[residue_idx - 1] = mut_res
            mut_name = f"{self.protein_name}_{wt_res}{residue_idx}{mut_res}"
            
            seq_name_list.append(mut_name)
            seq_list.append("".join(mut_seq))
        
        result_df = pd.DataFrame({
            'mutation_name': seq_name_list,
            'sequence': seq_list
        })
        
        logger.info(f"Generated {len(result_df)} single {mutation_type} mutations ({position} side)")
        return result_df
    
    def generate_pair_mutations(
        self,
        candidates: pd.DataFrame,
        mutation_type: str = 'charge',
        interaction_type: str = 'attractive'
    ) -> pd.DataFrame:
        """
        Generate double mutations (both residues in a pair).
        
        Args:
            candidates: DataFrame with mutation candidates
            mutation_type: 'charge', 'polar', or 'hydrophobic'
            interaction_type: 'attractive' or 'repulsive'
        
        Returns:
            DataFrame with mutation sequences
        """
        seq_list = []
        seq_name_list = []
        
        for _, row in candidates.iterrows():
            i, j = int(row['r_1']), int(row['r_2'])
            
            if i in self.forbidden_regions or j in self.forbidden_regions:
                logger.debug(f"Skipping forbidden pair ({i}, {j})")
                continue
            
            mut_seq = list(self.sequence)
            wt_res_i = self.seq_ref_dic[i]
            wt_res_j = self.seq_ref_dic[j]
            
            if mutation_type == 'charge':
                if interaction_type == 'attractive':
                    # Both to same charge to enhance attraction with opposite partner
                    for charge in ['E', 'K']:
                        mut_seq = list(self.sequence)
                        mut_seq[i - 1] = charge
                        mut_seq[j - 1] = charge
                        mut_name = f"{self.protein_name}_{wt_res_i}{i}{charge}_{wt_res_j}{j}{charge}"
                        seq_name_list.append(mut_name)
                        seq_list.append("".join(mut_seq))
                else:  # repulsive
                    # Neutralize both
                    mut_seq[i - 1] = 'A'
                    mut_seq[j - 1] = 'A'
                    mut_name = f"{self.protein_name}_{wt_res_i}{i}A_{wt_res_j}{j}A"
                    seq_name_list.append(mut_name)
                    seq_list.append("".join(mut_seq))
            
            elif mutation_type == 'polar':
                if interaction_type == 'attractive':
                    mut_seq[i - 1] = 'Q'
                    mut_seq[j - 1] = 'Q'
                    mut_name = f"{self.protein_name}_{wt_res_i}{i}Q_{wt_res_j}{j}Q"
                else:
                    mut_seq[i - 1] = 'A'
                    mut_seq[j - 1] = 'A'
                    mut_name = f"{self.protein_name}_{wt_res_i}{i}A_{wt_res_j}{j}A"
                
                seq_name_list.append(mut_name)
                seq_list.append("".join(mut_seq))
            
            elif mutation_type == 'hydrophobic':
                if interaction_type == 'attractive':
                    mut_seq[i - 1] = 'L'
                    mut_seq[j - 1] = 'L'
                    mut_name = f"{self.protein_name}_{wt_res_i}{i}L_{wt_res_j}{j}L"
                else:
                    mut_seq[i - 1] = 'S'
                    mut_seq[j - 1] = 'S'
                    mut_name = f"{self.protein_name}_{wt_res_i}{i}S_{wt_res_j}{j}S"
                
                seq_name_list.append(mut_name)
                seq_list.append("".join(mut_seq))
        
        result_df = pd.DataFrame({
            'mutation_name': seq_name_list,
            'sequence': seq_list
        })
        
        logger.info(f"Generated {len(result_df)} pair {mutation_type} mutations")
        return result_df
    
    def generate_chunk_mutations(
        self,
        candidates: pd.DataFrame,
        position: str = 'left',
        mutation_type: str = 'charge',
        interaction_type: str = 'attractive',
        chunk_size: int = 3
    ) -> pd.DataFrame:
        """
        Generate chunk mutations (mutate entire chunk).
        
        Args:
            candidates: DataFrame with mutation candidates
            position: 'left' (r_1) or 'right' (r_2)
            mutation_type: 'charge', 'polar', or 'hydrophobic'
            interaction_type: 'attractive' or 'repulsive'
            chunk_size: Size of chunk to mutate
        
        Returns:
            DataFrame with mutation sequences
        """
        seq_list = []
        seq_name_list = []
        
        pos_col = f'r_{1 if position == "left" else 2}'
        charge_col = f'{pos_col}_charge'
        opposite_charge_col = f'r_{2 if position == "left" else 1}_charge'
        
        offset = chunk_size // 2
        
        for _, row in candidates.iterrows():
            center = int(row[pos_col])
            chunk_start = center - offset
            chunk_end = center + offset
            
            # Check forbidden regions
            if any(i in self.forbidden_regions for i in range(chunk_start, chunk_end + 1)):
                logger.debug(f"Skipping forbidden chunk at {center}")
                continue
            
            # Ensure within sequence bounds
            if chunk_start < 1 or chunk_end > self.seq_length:
                continue
            
            mut_seq = list(self.sequence)
            
            if mutation_type == 'charge':
                if interaction_type == 'attractive':
                    charge_status = row[charge_col]
                    opposite_charge = row[opposite_charge_col]
                    
                    # Choose mutation based on charge
                    if charge_status > 0:
                        mut_aa = 'E'
                    elif charge_status < 0:
                        mut_aa = 'K'
                    elif opposite_charge >= 0:
                        mut_aa = 'K'
                    else:
                        mut_aa = 'E'
                else:  # repulsive
                    mut_aa = 'A'  # Neutralize
                
                for i in range(chunk_start, chunk_end + 1):
                    mut_seq[i - 1] = mut_aa
                
                # Generate mutation name
                mut_parts = [f"{self.seq_ref_dic[i]}{i}{mut_aa}" 
                            for i in range(chunk_start, chunk_end + 1)]
                mut_name = f"{self.protein_name}_{'_'.join(mut_parts)}"
            
            elif mutation_type == 'polar':
                mut_aa = 'Q' if interaction_type == 'attractive' else 'A'
                for i in range(chunk_start, chunk_end + 1):
                    mut_seq[i - 1] = mut_aa
                
                mut_parts = [f"{self.seq_ref_dic[i]}{i}{mut_aa}" 
                            for i in range(chunk_start, chunk_end + 1)]
                mut_name = f"{self.protein_name}_{'_'.join(mut_parts)}"
            
            elif mutation_type == 'hydrophobic':
                mut_aa = 'L' if interaction_type == 'attractive' else 'S'
                for i in range(chunk_start, chunk_end + 1):
                    mut_seq[i - 1] = mut_aa
                
                mut_parts = [f"{self.seq_ref_dic[i]}{i}{mut_aa}" 
                            for i in range(chunk_start, chunk_end + 1)]
                mut_name = f"{self.protein_name}_{'_'.join(mut_parts)}"
            
            seq_name_list.append(mut_name)
            seq_list.append("".join(mut_seq))
        
        result_df = pd.DataFrame({
            'mutation_name': seq_name_list,
            'sequence': seq_list
        })
        
        logger.info(f"Generated {len(result_df)} chunk {mutation_type} mutations ({position} side)")
        return result_df
    
    def save_mutations(self, mutations_df: pd.DataFrame, output_path: Union[str, Path]):
        """
        Save mutations to CSV file.
        
        Args:
            mutations_df: DataFrame with mutation_name and sequence columns
            output_path: Path to output CSV file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        mutations_df.to_csv(output_path, header=False, index=False)
        logger.info(f"Saved {len(mutations_df)} mutations to {output_path}")
    
    def generate_full_scan(
        self,
        output_dir: Union[str, Path],
        interaction_type: str = 'attractive',
        min_chunk_strength: float = 1.0
    ) -> Dict[str, pd.DataFrame]:
        """
        Generate complete mutation scan (all types).
        
        Args:
            output_dir: Directory to save mutation files
            interaction_type: 'attractive' or 'repulsive'
            min_chunk_strength: Minimum chunk strength threshold
        
        Returns:
            Dictionary of mutation DataFrames
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Generating full mutation scan ({interaction_type})")
        
        # Filter and process interactions
        selected = self.filter_interactions()
        chunk_data = self.calculate_chunk_strength(selected)
        candidates = self.identify_candidates(chunk_data, interaction_type, min_chunk_strength)
        
        results = {}
        
        # Single mutations
        for position in ['left', 'right']:
            for mut_type in ['charge', 'polar']:
                mutations = self.generate_single_mutations(
                    candidates, position, mut_type, interaction_type
                )
                key = f"{position}_single_{mut_type}"
                results[key] = mutations
                
                if len(mutations) > 0:
                    filename = f"{key}_{interaction_type}.csv"
                    self.save_mutations(mutations, output_dir / filename)
        
        # Pair mutations
        for mut_type in ['charge', 'polar']:
            mutations = self.generate_pair_mutations(
                candidates, mut_type, interaction_type
            )
            key = f"pair_{mut_type}"
            results[key] = mutations
            
            if len(mutations) > 0:
                filename = f"{key}_{interaction_type}.csv"
                self.save_mutations(mutations, output_dir / filename)
        
        # Chunk mutations
        for position in ['left', 'right']:
            for mut_type in ['charge']:
                mutations = self.generate_chunk_mutations(
                    candidates, position, mut_type, interaction_type
                )
                key = f"{position}_chunk_{mut_type}"
                results[key] = mutations
                
                if len(mutations) > 0:
                    filename = f"{key}_{interaction_type}.csv"
                    self.save_mutations(mutations, output_dir / filename)
        
        logger.info(f"Full scan complete. Generated {sum(len(df) for df in results.values())} total mutations")
        return results


def scan_mutations_from_csv(
    interaction_csv: Union[str, Path],
    sequence: str,
    protein_name: str,
    output_dir: Union[str, Path],
    interaction_type: str = 'attractive',
    forbidden_regions: Optional[List[int]] = None,
    min_chunk_strength: float = 1.0
) -> Dict[str, pd.DataFrame]:
    """
    Convenience function to scan mutations from interaction CSV file.
    
    Args:
        interaction_csv: Path to interaction map CSV file
        sequence: Protein sequence
        protein_name: Name of the protein
        output_dir: Directory to save mutation files
        interaction_type: 'attractive' or 'repulsive'
        forbidden_regions: List of positions that should not be mutated
        min_chunk_strength: Minimum chunk strength threshold
    
    Returns:
        Dictionary of mutation DataFrames
    
    Example:
        >>> mutations = scan_mutations_from_csv(
        ...     "protein_interaction.csv",
        ...     "ACDEFGHIKLMNPQRSTVWY",
        ...     "MyProtein",
        ...     "./mutations",
        ...     interaction_type='attractive'
        ... )
    """
    # Load interaction data
    df = pd.read_csv(interaction_csv, index_col=0)
    
    # Create scanner and run full scan
    scanner = MutationScanner(df, sequence, protein_name, forbidden_regions)
    return scanner.generate_full_scan(output_dir, interaction_type, min_chunk_strength)
