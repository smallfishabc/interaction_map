"""Command-line interface for IDP Interaction Map."""

import argparse
import logging
import sys
from pathlib import Path

from idp_interaction_map import __version__
from idp_interaction_map.core import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)


def print_banner() -> None:
    """Print welcome banner."""
    banner = """
#################################################################

.................................................................

.................IDP Interaction Map.........................

.................................................................

#################################################################

IDP Interaction Map v{version}
Analyzing intramolecular interactions in intrinsically disordered proteins

""".format(
        version=__version__
    )
    print(banner)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description="Analyze intramolecular interactions in IDPs from MD simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Coarse-grained simulation with multiple trajectories (default)
  idp-interaction-map -d ./data -n my_protein -r 5 --mode cg

  # All-atom simulation with CA selection
  idp-interaction-map -d ./data -n my_protein -r 5 --mode all-atom

  # Single trajectory file (coarse-grained)
  idp-interaction-map -d ./data -n my_protein -x trajectory.xtc

  # Single trajectory file (all-atom)
  idp-interaction-map -d ./data -n my_protein -x trajectory.dcd --mode all-atom

  # Custom normalization parameters
  idp-interaction-map -d ./data -n my_protein -r 10 --norm-a 1.64 --norm-b -1.33
        """,
    )

    parser.add_argument(
        "-p",
        "--pdb",
        type=str,
        default="__START_0.pdb",
        help="PDB topology file (default: __START_0.pdb)",
    )

    parser.add_argument(
        "-x",
        "--xtc",
        type=str,
        help="XTC trajectory file (required for single trajectory mode)",
    )

    parser.add_argument(
        "-d",
        "--data-dir",
        type=str,
        required=True,
        help="Data directory containing trajectory files and seq.txt",
    )

    parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=True,
        help="Protein name (used for output files)",
    )

    parser.add_argument(
        "-r",
        "--repeat",
        type=int,
        default=None,
        help="Number of trajectory files (generates __traj_0.xtc, __traj_1.xtc, etc.)",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: same as data-dir)",
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.2,
        help="Contact distance cutoff in nm (default: 1.2)",
    )

    parser.add_argument(
        "--mode",
        type=str,
        choices=["cg", "all-atom"],
        default="cg",
        help="Simulation type: 'cg' for coarse-grained, 'all-atom' for all-atom with CA selection (default: cg)",
    )

    parser.add_argument(
        "--norm-a",
        type=float,
        default=None,
        help="Normalization parameter a (auto-selected if not specified: CG=3.81, all-atom=1.64)",
    )

    parser.add_argument(
        "--norm-b",
        type=float,
        default=None,
        help="Normalization parameter b (auto-selected if not specified: CG=-1.51, all-atom=-1.33)",
    )

    parser.add_argument(
        "--read-from-file",
        action="store_true",
        help="Read existing contact map instead of recalculating",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"IDP Interaction Map {__version__}",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate command-line arguments.

    Args:
        args: Parsed arguments

    Raises:
        ValueError: If arguments are invalid
    """
    data_dir = Path(args.data_dir)

    if not data_dir.exists():
        raise ValueError(f"Data directory does not exist: {data_dir}")

    # Check for sequence file
    seq_file = data_dir / "seq.txt"
    if not seq_file.exists():
        raise ValueError(
            f"Sequence file not found: {seq_file}\n"
            "Please create seq.txt in the data directory with the protein sequence"
        )

    # Validate trajectory input
    if args.repeat is None and args.xtc is None:
        raise ValueError(
            "Must specify either --xtc (single trajectory) or --repeat (multiple trajectories)"
        )

    if args.repeat is not None and args.xtc is not None:
        logger.warning(
            "Both --repeat and --xtc specified. Using --xtc for single trajectory mode"
        )


def main() -> int:
    """
    Main entry point for CLI.

    Returns:
        Exit code (0 for success, 1 for error)
    """
    try:
        print_banner()

        args = parse_args()

        # Set logging level
        if args.verbose:
            logging.getLogger().setLevel(logging.DEBUG)

        # Validate arguments
        validate_args(args)

        # Set paths
        data_dir = Path(args.data_dir)
        output_dir = Path(args.output_dir) if args.output_dir else data_dir

        # Read sequence
        logger.info("Reading protein sequence")
        sequence = read_sequence_from_txt(data_dir)
        logger.info(f"Sequence length: {len(sequence)}")

        # Determine trajectory input
        if args.xtc:
            xtc_input = args.xtc
            logger.info(f"Using single trajectory: {args.xtc}")
        else:
            xtc_input = args.repeat
            logger.info(f"Using {args.repeat} trajectory files")

        # Determine simulation mode
        use_ca = (args.mode == "all-atom")
        logger.info(f"Simulation mode: {args.mode} (CA-only: {use_ca})")

        # Run analysis
        logger.info("Starting interaction map analysis")
        analyze_interaction_map(
            name=args.name,
            trajectory_path=data_dir,
            sequence=sequence,
            output_dir=output_dir,
            pdb_top=args.pdb,
            xtc_input=xtc_input,
            read_from_file=args.read_from_file,
            use_ca=use_ca,
            norm_a=args.norm_a,
            norm_b=args.norm_b,
        )

        logger.info("Analysis completed successfully!")
        return 0

    except KeyboardInterrupt:
        logger.info("\nAnalysis interrupted by user")
        return 130

    except Exception as e:
        logger.error(f"Error: {e}", exc_info=args.verbose if "args" in locals() else False)
        return 1


if __name__ == "__main__":
    sys.exit(main())
