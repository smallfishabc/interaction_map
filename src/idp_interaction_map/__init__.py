"""IDP Interaction Map package for analyzing intrinsically disordered proteins."""

__version__ = "1.0.0"
__author__ = "Feng Yu"

from idp_interaction_map.contact_map import ContactProbData, generate_contact, load_traj_protein
from idp_interaction_map.core import analyze_interaction_map
from idp_interaction_map.normalization import normalize_interaction_map

__all__ = [
    "ContactProbData",
    "generate_contact",
    "load_traj_protein",
    "analyze_interaction_map",
    "normalize_interaction_map",
]
