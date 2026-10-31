#!/usr/bin/env python
"""Run the interaction map without installing the package.

    python scripts/main.py -p protein.pdb -x traj.xtc -name protein

Once the package is installed, `interaction-map` with the same arguments does the same
thing and this file is not needed.
"""
import os
import sys

# src/ onto the path, so the package imports without being installed first.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from ensemble_map.cli import main

if __name__ == '__main__':
    main()
