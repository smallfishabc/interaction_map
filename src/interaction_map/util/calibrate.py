# -*- coding: utf-8 -*-
"""
Fit the reference contact probability curve, p_ref(s) = a * s ** b.

This is the curve the interaction strength is measured against: the behaviour of a chain
with no interactions of its own, so that a real protein's departure from it can be read as
attraction or repulsion. It is fitted on glycine-serine linkers, which are taken to be
that reference chain.

The method is the one from the GS notebooks, kept as it was:

    per-pair contact probability  = fraction of frames in which the pair is inside cutoff
    per-separation probability    = mean of those, over all pairs at that separation
    fit                           = curve_fit of a * s ** b against separation

Deliberately independent of the rest of the package. It takes a cutoff and a scheme as
arguments rather than reading a profile, because its whole purpose is to produce the
numbers a profile is later given, and a calibration that quietly inherited the constants
it is meant to measure would be worth nothing.

Run it when you have reference trajectories:

    python -m interaction_map.util.calibrate GS16.pdb GS16.dcd --cutoff 1.2 --all-pairs
    python -m interaction_map.util.calibrate __START_0.pdb __traj_*.xtc --out gs.json

@author: ShaharGroup-fyu
"""
from __future__ import print_function

import argparse
import datetime
import glob
import json
import os
import sys

import numpy as np
import pandas as pd


# The reference curve itself. Kept here, rather than imported, so this module can be run
# against a checkout that has no package installed.
def power_law(separation, a, b):
    return a * np.asarray(separation, dtype=float) ** b


def pair_separations_and_distances(traj, scheme='ca', all_pairs=False, min_separation=3):
    """
    Residue pair distances for every frame.

    :param traj: mdtraj trajectory
    :param scheme: how a residue-residue distance is measured. 'ca' for all-atom,
                   'closest-heavy' for one bead per residue.
    :param all_pairs: select every residue pair explicitly, which is what a bead model
                      needs. Otherwise mdtraj's contacts='all' is used, which already
                      drops pairs closer than 3 along the sequence.
    :param min_separation: pairs closer than this along the sequence are discarded.
    :return: (separations, distances) where distances is pairs x frames
    """
    import mdtraj as md

    if all_pairs:
        indices = traj.top.select_pairs('all', 'all')
        indices = indices[indices[:, 1] - indices[:, 0] >= min_separation]
        cont, pairs = md.compute_contacts(traj, indices, scheme=scheme)
    else:
        cont, pairs = md.compute_contacts(traj, contacts='all', scheme=scheme,
                                          ignore_nonprotein=True)
    separations = pairs[:, 1] - pairs[:, 0]
    keep = separations >= min_separation
    return separations[keep], cont.T[keep]


def contact_probability_by_separation(separations, distances, cutoff):
    """
    Mean contact probability at each sequence separation.

    Each pair contributes the fraction of frames in which it is inside the cutoff; those
    are then averaged within each separation.

    :return: (probability, n_pairs) both indexed by separation
    """
    frame_count = distances.shape[1]
    per_pair = (distances < cutoff).sum(axis=1) / float(frame_count)
    table = pd.DataFrame({'distance': separations, 'p': per_pair})
    grouped = table.groupby('distance')['p']
    return grouped.mean(), grouped.size()


def fit_power_law(probability, p0=(1.0, -1.0)):
    """
    Least squares fit of a * s ** b.

    :param probability: contact probability indexed by separation
    :return: dict with a, b, their standard errors, r squared, and the residuals
    """
    from scipy.optimize import curve_fit

    separations = np.asarray(probability.index, dtype=float)
    observed = np.asarray(probability.values, dtype=float)
    usable = observed > 0
    if usable.sum() < 3:
        raise ValueError('need at least three separations with a non-zero contact '
                         'probability to fit, got {0}'.format(int(usable.sum())))

    popt, pcov = curve_fit(power_law, separations[usable], observed[usable], p0=p0,
                           maxfev=10000)
    predicted = power_law(separations, *popt)
    residuals = observed - predicted
    # R squared on the points that were actually fitted.
    ss_res = float((residuals[usable] ** 2).sum())
    ss_tot = float(((observed[usable] - observed[usable].mean()) ** 2).sum())
    errors = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else [float('nan')] * 2
    return {'a': float(popt[0]),
            'b': float(popt[1]),
            'a_stderr': float(errors[0]),
            'b_stderr': float(errors[1]),
            'r_squared': float(1 - ss_res / ss_tot) if ss_tot > 0 else float('nan'),
            'separations': [int(s) for s in separations],
            'observed': [float(v) for v in observed],
            'predicted': [float(v) for v in predicted],
            'residuals': [float(v) for v in residuals],
            'n_separations_fitted': int(usable.sum()),
            'n_separations_dropped': int((~usable).sum())}


def calibrate(topology, trajectories, cutoff=0.8, scheme='ca', all_pairs=False,
              min_separation=3, max_separation=None, source=None):
    """
    Fit the reference curve on one or more reference trajectories.

    :param topology: topology file
    :param trajectories: one trajectory file, or a list of them
    :param cutoff: contact cutoff in nm. The curve is only valid for this cutoff.
    :param scheme: 'ca' for all-atom, 'closest-heavy' for a bead model
    :param all_pairs: True for a bead model, where contacts='all' finds nothing useful
    :param min_separation: shortest separation included in the fit
    :param max_separation: longest separation included, or None for all of them
    :param source: a note recording what was fitted, free text
    :return: dict holding the fit and enough provenance to know what it describes
    """
    import mdtraj as md

    if isinstance(trajectories, str):
        trajectories = [trajectories]
    trajectories = list(trajectories)

    traj = md.load(trajectories, top=topology)
    separations, distances = pair_separations_and_distances(
        traj, scheme=scheme, all_pairs=all_pairs, min_separation=min_separation)
    probability, pair_counts = contact_probability_by_separation(separations, distances, cutoff)

    if max_separation is not None:
        probability = probability[probability.index <= max_separation]
        pair_counts = pair_counts[pair_counts.index <= max_separation]

    fit = fit_power_law(probability)
    fit.update({
        'cutoff': float(cutoff),
        'scheme': scheme,
        'all_pairs': bool(all_pairs),
        # The range the curve was fitted over. Outside it the curve is an extrapolation,
        # and a power law with a negative exponent climbs without limit as s falls, so it
        # can return values above 1, which no probability can be.
        'separation_min': int(probability.index.min()),
        'separation_max': int(probability.index.max()),
        'pairs_per_separation': {int(k): int(v) for k, v in pair_counts.items()},
        'n_frames': int(traj.n_frames),
        'n_residues': int(traj.n_residues),
        'n_trajectories': len(trajectories),
        'topology': os.path.basename(topology),
        'trajectories': [os.path.basename(t) for t in trajectories],
        'source': source or '',
        'fitted': datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
    })
    return fit


def probability_ceiling(fit):
    """
    The largest separation at which the fitted curve returns a value above 1.

    Reported so the range can be seen, not acted on. What it means for a map is a question
    about the method, not about this function.

    :return: the largest integer separation where p_ref > 1, or None if there is none
    """
    a, b = fit['a'], fit['b']
    if a <= 1.0 or b >= 0:
        return None
    # a * s**b > 1  <=>  s < a ** (-1/b)
    boundary = a ** (-1.0 / b)
    largest = int(np.floor(boundary))
    if power_law(largest, a, b) <= 1.0:
        largest -= 1
    return largest if largest >= 1 else None


def summarise(fit):
    """A short human readable report of a fit."""
    lines = [
        'p_ref(s) = {0:.4g} * s ** {1:.4g}'.format(fit['a'], fit['b']),
        '  standard error      a {0:.3g}, b {1:.3g}'.format(fit['a_stderr'], fit['b_stderr']),
        '  R squared           {0:.5f}'.format(fit['r_squared']),
        '  fitted over         separation {0} to {1}, {2} points'.format(
            fit['separation_min'], fit['separation_max'], fit['n_separations_fitted']),
        '  cutoff              {0} nm, scheme {1!r}'.format(fit['cutoff'], fit['scheme']),
        '  from                {0} frames, {1} residues, {2} trajectory file(s)'.format(
            fit['n_frames'], fit['n_residues'], fit['n_trajectories']),
    ]
    if fit['n_separations_dropped']:
        lines.append('  dropped             {0} separation(s) with no contacts at all'.format(
            fit['n_separations_dropped']))
    ceiling = probability_ceiling(fit)
    if ceiling:
        lines.append('  p_ref above 1       at separation {0} and below'.format(ceiling))
    worst = int(np.argmax(np.abs(fit['residuals'])))
    lines.append('  largest residual    {0:+.4g} at separation {1}'.format(
        fit['residuals'][worst], fit['separations'][worst]))
    return '\n'.join(lines)


def save(fit, path):
    with open(path, 'w') as handle:
        json.dump(fit, handle, indent=2, sort_keys=True)
    return path


def load(path):
    with open(path) as handle:
        fit = json.load(handle)
    # json has no integer keys, so the separations come back as strings. Put them back,
    # otherwise a fit that has been through a file indexes differently from one that has
    # not, which is the kind of difference that is noticed far from where it was caused.
    if 'pairs_per_separation' in fit:
        fit['pairs_per_separation'] = {int(k): int(v)
                                       for k, v in fit['pairs_per_separation'].items()}
    return fit


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='interaction_map.util.calibrate',
        description='Fit the reference contact probability curve on reference '
                    'trajectories, usually glycine-serine linkers.')
    parser.add_argument('topology', help='topology file')
    parser.add_argument('trajectory', nargs='+',
                        help='one or more trajectory files; wildcards are expanded')
    parser.add_argument('--cutoff', type=float, default=0.8,
                        help='contact cutoff in nm. The fitted curve belongs to this '
                             'cutoff and to no other. Default 0.8.')
    parser.add_argument('--scheme', default=None,
                        help="how a residue-residue distance is measured. Defaults to "
                             "'ca', or 'closest-heavy' with --all-pairs.")
    parser.add_argument('--all-pairs', action='store_true',
                        help='select every residue pair explicitly. Needed for a bead '
                             'model, where the beads are not named CA.')
    parser.add_argument('--sep', nargs=2, type=int, metavar=('MIN', 'MAX'), default=None,
                        help='separation range to fit over. Default: 3 upwards.')
    parser.add_argument('--source', default=None,
                        help='free text recording what was fitted, stored with the result')
    parser.add_argument('--out', default=None, help='write the fit to this json file')
    args = parser.parse_args(argv)

    paths = []
    for pattern in args.trajectory:
        matched = sorted(glob.glob(pattern))
        paths.extend(matched if matched else [pattern])

    scheme = args.scheme or ('closest-heavy' if args.all_pairs else 'ca')
    sep_min, sep_max = (args.sep if args.sep else (3, None))

    fit = calibrate(args.topology, paths, cutoff=args.cutoff, scheme=scheme,
                    all_pairs=args.all_pairs, min_separation=sep_min,
                    max_separation=sep_max, source=args.source)
    print(summarise(fit))
    if args.out:
        save(fit, args.out)
        print('\nwritten to {0}'.format(args.out))
    return fit


if __name__ == '__main__':
    main()
