"""Interaction maps for intrinsically disordered regions.

The analysis in four steps: contact probability per residue pair, divided by a reference
polymer curve and logged, binned into five levels, drawn as arcs along the sequence.
Optionally, condition the ensemble on one or more interactions and compare the two halves.

Only the first step reads trajectories, and only that step needs mdtraj. Everything from
the contact table onwards -- the reference curve, the interaction strength, the levels and
the figure -- is numpy, pandas and matplotlib.

So the modules that need mdtraj are imported on demand rather than when the package is
imported. Starting from a stored contact table, this works with no mdtraj installed:

    from interaction_map import normalization, interaction_plot

while these pull mdtraj in at the moment they are first used:

    from interaction_map import contact_map_generation, search_chunk, default_function
"""

import importlib

from . import profiles
from . import readpath
from . import schema
from . import normalization
from . import interaction_plot
from . import showoff

# Imported on first use, because each of these reads trajectories.
_ON_DEMAND = ('contact_map_generation', 'search_chunk', 'default_function')


def __getattr__(name):
    if name in _ON_DEMAND:
        module = importlib.import_module('.' + name, __name__)
        globals()[name] = module
        return module
    raise AttributeError('module {0!r} has no attribute {1!r}'.format(__name__, name))


def __dir__():
    return sorted(list(globals()) + list(_ON_DEMAND))


__version__ = '0.1.0'

__all__ = ['profiles', 'readpath', 'schema', 'normalization', 'interaction_plot',
           'showoff', 'contact_map_generation', 'search_chunk', 'default_function']
