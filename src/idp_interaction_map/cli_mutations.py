"""Command-line interface for mutation scanning."""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from idp_interaction_map import __version__, scan_mutations_from_csv
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

...............Mutation Scanner for IDP Interaction Maps..........

#################################################################

IDP Interaction Map - Mutation Scanner v{version}
Generate targeted mutations based on interaction analysis

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
        description="Generate mutations to modulate protein interactions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate attractive mutations
  idp-mutation-scan -i protein_interaction.csv -s ACDEFG... -n MyProtein -o ./mutations --type attractive

  # Generate repulsive mutations with forbidden regions
  idp-mutation-scan -i protein_interaction.csv -s ACDEFG... -n MyProtein -o ./mutations \\
      --type repulsive --forbidden 1,2,3,10-20

  # Load sequence from file
  idp-mutation-scan -i protein_interaction.csv -f seq.txt -n MyProtein -o ./mutations

  # Custom chunk strength threshold
  idp-mutation-scan -i protein_interaction.csv -s ACDEFG... -n MyProtein -o ./mutations \\
      --min-strength 2.0
        """,
    )

    parser.add_argument(
        "-i",
        "--interaction-csv",
        type=str,
        required=True,
        help="Input interaction CSV file from IDP interaction map analysis",
    )

    parser.add_argument(
        "-s",
        "--sequence",
        type=str,
        help="Protein sequence (one-letter amino acid codes)",
    )

    parser.add_argument(
        "-f",
        "--sequence-file",
        type=str,
        help="File containing protein sequence (alternative to --sequence)",
    )

    parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=True,
        help="Protein name (used in mutation naming)",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for mutation files",
    )

    parser.add_argument(
        "--type",
        type=str,
        choices=["attractive", "repulsive"],
        default="attractive",
        help="Type of mutations to generate (default: attractive)",
    )

    parser.add_argument(
        "--forbidden",
        type=str,
        help="Forbidden regions (comma-separated or ranges, e.g., '1,2,3' or '1-10,15-20')",
    )

    parser.add_argument(
        "--min-strength",
        type=float,
        default=1.0,
        help="Minimum chunk strength threshold (default: 1.0)",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"IDP Mutation Scanner v{__version__}",
    )

    return parser.parse_args()


def parse_forbidden_regions(forbidden_str: str) -> List[int]:
    """
    Parse forbidden regions string into list of positions.
    
    Args:
        forbidden_str: String like '1,2,3' or '1-10,15-20'
    
    Returns:
        List of forbidden positions
    """
    if not forbidden_str:
        return []
    
    forbidden = []
    parts = forbidden_str.split(',')
    
    for part in parts:
        part = part.strip()
        if '-' in part:
            # Range like '1-10'
            start, end = map(int, part.split('-'))
            forbidden.extend(range(start, end + 1))
        else:
            # Single number
            forbidden.append(int(part))
    
    return sorted(set(forbidden))


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate command-line arguments.

    Args:
        args: Parsed arguments

    Raises:
        ValueError: If arguments are invalid
    """
    # Check that either sequence or sequence_file is provided
    if not args.sequence and not args.sequence_file:
        raise ValueError("Either --sequence or --sequence-file must be provided")
    
    if args.sequence and args.sequence_file:
        raise ValueError("Provide only one of --sequence or --sequence-file, not both")
    
    # Check that interaction CSV exists
    interaction_path = Path(args.interaction_csv)
    if not interaction_path.exists():
        raise ValueError(f"Interaction CSV file not found: {args.interaction_csv}")
    
    # Check sequence file if provided
    if args.sequence_file:
        seq_path = Path(args.sequence_file)
        if not seq_path.exists():
            raise ValueError(f"Sequence file not found: {args.sequence_file}")


def main() -> int:
    """
    Main entry point for mutation scanner CLI.

    Returns:
        Exit code (0 for success, non-zero for error)
    """
    # Parse arguments
    args = parse_args()

    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    # Print banner
    print_banner()

    try:
        # Validate arguments
        validate_args(args)

        # Get sequence
        if args.sequence_file:
            logger.info(f"Reading sequence from {args.sequence_file}")
            sequence = read_sequence_from_txt(args.sequence_file)
        else:
            sequence = args.sequence
        
        logger.info(f"Protein: {args.name}")
        logger.info(f"Sequence length: {len(sequence)}")
        logger.info(f"Interaction type: {args.type}")

        # Parse forbidden regions
        forbidden_regions = parse_forbidden_regions(args.forbidden) if args.forbidden else []
        if forbidden_regions:
            logger.info(f"Forbidden regions: {len(forbidden_regions)} positions")

        # Run mutation scan
        logger.info(f"Reading interaction data from {args.interaction_csv}")
        mutations = scan_mutations_from_csv(
            interaction_csv=args.interaction_csv,
            sequence=sequence,
            protein_name=args.name,
            output_dir=args.output_dir,
            interaction_type=args.type,
            forbidden_regions=forbidden_regions,
            min_chunk_strength=args.min_strength,
        )

        # Report results
        total_mutations = sum(len(df) for df in mutations.values())
        logger.info("")
        logger.info("=" * 60)
        logger.info("MUTATION SCAN COMPLETE")
        logger.info("=" * 60)
        logger.info(f"Generated {total_mutations} total mutations")
        logger.info(f"Output directory: {Path(args.output_dir).absolute()}")
        
        for name, df in mutations.items():
            if len(df) > 0:
                logger.info(f"  - {name}: {len(df)} mutations")
        
        logger.info("")
        logger.info("Mutation files saved with naming format:")
        logger.info(f"  {args.name}_<ResidueID><MutationCode>")
        logger.info("")

        return 0

    except KeyboardInterrupt:
        logger.info("\nMutation scan interrupted by user")
        return 130

    except Exception as e:
        logger.error(f"Error during mutation scan: {e}", exc_info=args.verbose)
        return 1


if __name__ == "__main__":
    sys.exit(main())
