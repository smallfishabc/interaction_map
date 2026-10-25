"""Interaction maps for intrinsically disordered regions.

The analysis in four steps: contact probability per residue pair, divided by a reference
polymer curve and logged, binned into five levels, drawn as arcs along the sequence.
Optionally, condition the ensemble on one or more interactions and compare the two halves.
"""

from . import profiles
from . import readpath
from . import schema
from . import normalization
from . import contact_map_generation
from . import interaction_plot
from . import search_chunk
from . import default_function
from . import showoff

__version__ = '0.1.0'

__all__ = ['profiles', 'readpath', 'schema', 'normalization', 'contact_map_generation',
           'interaction_plot', 'search_chunk', 'default_function', 'showoff']
