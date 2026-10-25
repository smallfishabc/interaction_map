# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 14:37:47 2021

@author: ShaharGroup-fyu
"""
import argparse
import os
import profiles
import readpath
import showoff
import default_function

# Trajectory folder used by --database when --protein_directory is not given. Edit this
# line to run the simulation database from the IDE. The folders used for previous figures
# are listed in RUN_PATHS.md.
DEFAULT_PATH = r'F:\DATA_F\Interaction_map_target\Individual Proteins\p53-summary'



# Results go next to the trajectory unless --out says otherwise, in which case the folder
# is created if it is not there yet.
def output_directory(args, default):
    if not args.out:
        return default
    if not os.path.isdir(args.out):
        os.makedirs(args.out)
    return os.path.abspath(args.out)


# --search takes a flat list of residue numbers, read as pairs.
def search_pairs(args):
    values = args.search
    if len(values) % 2:
        raise SystemExit('--search needs an even number of residues, read as pairs')
    return [(values[i], values[i + 1]) for i in range(0, len(values), 2)]


# One topology plus one or more trajectory files. This is the published entry point.
def run_single(args):
    directory = args.protein_directory or '.'
    name = args.name
    if not name:
        # Name the map after the topology file when nothing better is given.
        name = os.path.splitext(os.path.basename(args.pdb))[0]
    sequence, source = readpath.find_sequence(directory, name, args.seq)
    print('sequence read from {0} ({1} residues)'.format(source, len(sequence)))

    # One trajectory stays a plain name; several are passed through as a list, and each
    # one counts as a repeat when the spread between repeats is reported.
    xtc_input = args.xtc[0] if len(args.xtc) == 1 else list(args.xtc)
    if args.search:
        pairs = search_pairs(args)
        if len(pairs) == 1:
            return default_function.search_interaction(
                name, pairs[0][0], pairs[0][1], directory, sequence, output_directory(args, directory),
                pdb_top=args.pdb, xtc_input=xtc_input,
                profile=args.profile, show=not args.no_show,
                window=args.window, chunk=args.window > 0)
        return default_function.compare_interactions(
            name, pairs, directory, sequence, output_directory(args, directory),
            pdb_top=args.pdb, xtc_input=xtc_input, profile=args.profile,
            mode=args.combine, window=args.window, chunk=args.window > 0,
            with_maps=args.with_maps, show=not args.no_show)
    return default_function.interaction_map_pairwise(
        name, directory, sequence, output_directory(args, directory),
        pdb_top=args.pdb, xtc_input=xtc_input,
        profile=args.profile, show=not args.no_show)


# Adapter for the simulation database layout:
#     <name>-summary/seq.fasta
#     <name>-summary/<restype>/S_<psi>/__START_0.pdb, __traj_0.xtc ...
# It resolves the folder convention and then calls the same single trajectory path.
def run_database(args):
    path = args.protein_directory or DEFAULT_PATH
    name = args.name or os.path.basename(os.path.normpath(path)).split('-')[0]
    psi = args.psi if args.psi else '0'
    residue = args.restype if args.restype else 'BB'
    repeat = int(args.repeat) if args.repeat else 5

    traj_path, sequence = default_function.multi_traj_pre(name, path, psi, residue)
    print(traj_path)
    if args.search:
        pairs = search_pairs(args)
        if len(pairs) == 1:
            return default_function.search_interaction(
                name, pairs[0][0], pairs[0][1], traj_path, sequence, output_directory(args, traj_path),
                xtc_input=repeat, profile=args.profile, show=not args.no_show,
                window=args.window, chunk=args.window > 0)
        return default_function.compare_interactions(
            name, pairs, traj_path, sequence, output_directory(args, traj_path),
            xtc_input=repeat, profile=args.profile, mode=args.combine,
            window=args.window, chunk=args.window > 0,
            with_maps=args.with_maps, show=not args.no_show)
    return default_function.interaction_map_pairwise(
        name, traj_path, sequence, output_directory(args, traj_path),
        xtc_input=repeat, read_from_file=False,
        profile=args.profile, show=not args.no_show)


if __name__ == "__main__":
    showoff.welcome()
    parser = argparse.ArgumentParser(
        description='Interaction map for an intrinsically disordered region.')
    parser.add_argument("--pdb", "-p", help="Topology file")
    parser.add_argument("--xtc", "-x", nargs='+',
                        help="One or more trajectory files. Several files are treated as "
                             "repeats of the same construct and the spread between them "
                             "is reported.")
    parser.add_argument("--protein_directory", "-dir",
                        help="Folder holding the trajectory, or the -summary folder with "
                             "--database. Defaults to the working directory.")
    parser.add_argument("--name", "-name", help="Name of the protein and of the output files")
    parser.add_argument("--seq", help="Sequence file. Plain text or FASTA. When not given, "
                                      "<name>.fasta, seq.txt and seq.fasta are tried in "
                                      "that order.")
    parser.add_argument("--out", "-o", help="Directory for the results. Defaults to the "
                                            "trajectory folder, which is where they have "
                                            "always been written.")
    parser.add_argument("--profile", default="all_atom", choices=sorted(profiles.PROFILES),
                        help="Analysis profile. 'all_atom' for CAMPARI/ABSINTH trajectories, "
                             "'cg' for CALVADOS trajectories. Decides the contact cutoff, the "
                             "reference curve and the interaction strength cutoffs.")
    parser.add_argument("--no_show", action="store_true",
                        help="Do not open the figure window. The png and svg are still written.")
    parser.add_argument("--database", "-db", action="store_true",
                        help="Read the simulation database layout, <name>-summary/<restype>/"
                             "S_<psi>, instead of a topology and trajectory given directly.")
    parser.add_argument("--repeat", "-repeat", help="Number of repeats (--database only)")
    parser.add_argument("--restype", "-restype",
                        help="Residue group whose transfer free energy is altered (--database only)")
    parser.add_argument("--psi", "-psi",
                        help="Psi value of the transfer free energy (--database only)")
    parser.add_argument("--search", nargs='+', metavar='R', type=int,
                        help="Conditional analysis. Give residue pairs: --search 13 24 "
                             "splits the ensemble on that one interaction and also rebuilds "
                             "the contact maps. More than one pair, --search 13 24 37 44, "
                             "compares the ensemble properties with and without the group, "
                             "combined according to --combine. Residues are numbered as in "
                             "the contact table.")
    parser.add_argument("--combine", choices=('any', 'all'), default='any',
                        help="How several interactions are combined. 'any' selects frames "
                             "where at least one is formed, 'all' only frames where every "
                             "one is formed together. Default any.")
    parser.add_argument("--with_maps", action="store_true",
                        help="With several interactions, also rebuild a contact map and "
                             "figure for each side of the split.")
    parser.add_argument("--window", type=int, default=1,
                        help="Half width of the residue window used to decide that the "
                             "interaction is present. 1, the default, is the published "
                             "3x3 window. 0 narrows it to the single pair.")
    parser.add_argument("--single_traj", "-straj", action="store_true",
                        help=argparse.SUPPRESS)  # accepted and ignored; now the default
    args = parser.parse_args()

    if args.database:
        run_database(args)
    elif args.pdb and args.xtc:
        run_single(args)
    else:
        parser.error('Give a topology and a trajectory with -p and -x, '
                     'or use --database for the simulation database layout.')
