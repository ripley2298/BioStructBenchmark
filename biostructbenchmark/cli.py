"""
biostructbenchmark/cli.py
Complete CLI module with corrected path validation for test compatibility
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Union, List, Optional, Tuple


def validate_path(input_path: str, must_exist: bool = True, 
                 allow_directory: bool = False) -> Path:
    """
    Validate file or directory path with comprehensive checks
    
    Args:
        input_path: Path to validate
        must_exist: Whether path must already exist
        allow_directory: Whether directories are acceptable
        
    Returns:
        Validated Path object (resolved to absolute path)
        
    Raises:
        ValueError: If validation fails
    """
    path = Path(input_path).resolve()
    
    if must_exist and not path.exists():
        raise ValueError(f"Path does not exist: {path}")
    
    if path.exists():
        if path.is_file():
            # Validate file properties
            if not os.access(path, os.R_OK):
                raise ValueError(f"No read permission: {path}")
            if path.stat().st_size == 0:
                raise ValueError(f"File is empty: {path}")
        elif path.is_dir():
            if not allow_directory:
                raise ValueError(f"Expected file, got directory: {path}")
            if not os.access(path, os.R_OK):
                raise ValueError(f"No read permission for directory: {path}")
        else:
            raise ValueError(f"Path is neither file nor directory: {path}")
    
    return path


def validate_file_path(input_path: str) -> Path:
    """
    Legacy wrapper for validate_path that preserves relative/absolute path style
    Specifically designed for test compatibility
    
    Args:
        input_path: File path string to validate
        
    Returns:
        Validated Path object (preserving relative vs absolute style)
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is empty or invalid
    """
    # Create path object preserving relative/absolute style
    path = Path(input_path)
    
    # Use resolve() only for validation checks
    resolved_path = path.resolve()
    
    # Check existence
    if not resolved_path.exists():
        raise FileNotFoundError(f"File not found: {input_path}")
    
    # Check if it's actually a file
    if not resolved_path.is_file():
        raise ValueError(f"Expected file, got directory: {input_path}")
    
    # Check read permissions
    if not os.access(resolved_path, os.R_OK):
        raise ValueError(f"No read permission: {input_path}")
    
    # Check if file is empty
    if resolved_path.stat().st_size == 0:
        raise ValueError(f"File is empty: {input_path}")
    
    # Return original path style (relative/absolute preserved)
    return path


def get_version() -> str:
    """Get package version from __init__.py"""
    try:
        # Try multiple possible locations for __init__.py
        possible_paths = [
            Path(__file__).parent / "__init__.py",
            Path("biostructbenchmark") / "__init__.py",
            Path(".") / "biostructbenchmark" / "__init__.py"
        ]
        
        for init_path in possible_paths:
            if init_path.exists():
                with open(init_path, 'r') as f:
                    for line in f:
                        if line.startswith("__version__"):
                            # Extract version safely
                            version_str = line.split("=")[1].strip().strip('"\'')
                            return version_str
        
        return "0.0.1"  # Fallback version
    except Exception:
        return "Unknown"


def find_structure_files(directory: Path) -> List[Path]:
    """
    Find all structure files in directory (PDB/CIF)
    
    Args:
        directory: Directory to search
        
    Returns:
        List of structure file paths
    """
    structure_extensions = {'.pdb', '.cif', '.mmcif', '.ent'}
    files = []
    
    for ext in structure_extensions:
        files.extend(directory.glob(f'*{ext}'))
        files.extend(directory.glob(f'*{ext.upper()}'))  # Also check uppercase
    
    return sorted(files)


def create_argument_parser() -> argparse.ArgumentParser:
    """Create comprehensive argument parser with all documented options"""
    parser = argparse.ArgumentParser(
        prog='biostructbenchmark',
        description='Compare experimental and predicted DNA-protein complex structures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file comparison
  biostructbenchmark observed.pdb predicted.pdb
  
  # Directory comparison with full analysis
  biostructbenchmark -e experimental/ -p predicted/ -o results/ --all-benchmarks
  
  # Specific analyses
  biostructbenchmark -e exp/ -p pred/ -o out/ --curves --bfactor --visualize
  
  # B-factor analysis only  
  biostructbenchmark -e exp/ -p pred/ -o out/ --bfactor
        """
    )
    
    # Version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'biostructbenchmark {get_version()}',
        help='Show version and exit'
    )
    
    # Input arguments - support both old and new style
    input_group = parser.add_mutually_exclusive_group(required=True)
    
    # New directory-based interface (preferred)
    input_group.add_argument(
        '-e', '--experimental',
        type=lambda x: validate_path(x, must_exist=True, allow_directory=True),
        help='Path to experimental structure file or directory'
    )
    
    # Legacy positional arguments (backward compatibility)
    input_group.add_argument(
        'legacy_files',
        nargs='*',
        help='Legacy format: observed_file predicted_file'
    )
    
    # Predicted structures (required if using -e)
    parser.add_argument(
        '-p', '--predicted',
        type=lambda x: validate_path(x, must_exist=True, allow_directory=True),
        help='Path to predicted structure file or directory (required with -e)'
    )
    
    # Output directory
    parser.add_argument(
        '-o', '--output',
        type=lambda x: validate_path(x, must_exist=False, allow_directory=True),
        default=Path('./biostructbenchmark_results'),
        help='Output directory for results (default: ./biostructbenchmark_results)'
    )
    
    # Analysis options
    analysis_group = parser.add_argument_group('Analysis Options')
    
    analysis_group.add_argument(
        '--all-benchmarks',
        action='store_true',
        help='Run all available analyses (RMSD, B-factor, CURVES+, consensus, mutations, visualization)'
    )
    
    analysis_group.add_argument(
        '--rmsd-only',
        action='store_true', 
        help='Perform only RMSD analysis (fastest option)'
    )
    
    analysis_group.add_argument(
        '--curves',
        action='store_true',
        help='Perform CURVES+ DNA geometry analysis'
    )
    
    analysis_group.add_argument(
        '--bfactor',
        action='store_true',
        help='Analyze B-factors vs confidence metrics'
    )
    
    analysis_group.add_argument(
        '--consensus',
        action='store_true', 
        help='Generate consensus error mapping (requires multiple structure pairs)'
    )
    
    analysis_group.add_argument(
        '--mutations',
        action='store_true',
        help='Detect and analyze mutations between structures'
    )
    
    analysis_group.add_argument(
        '--visualize',
        action='store_true',
        help='Generate publication-quality plots and visualizations'
    )
    
    # Processing options
    processing_group = parser.add_argument_group('Processing Options')
    
    processing_group.add_argument(
        '--parallel',
        type=int,
        metavar='N',
        help='Number of parallel processes for batch analysis (default: auto-detect)'
    )
    
    processing_group.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    processing_group.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress non-essential output'
    )
    
    # Advanced options
    advanced_group = parser.add_argument_group('Advanced Options')
    
    advanced_group.add_argument(
        '--rmsd-threshold',
        type=float,
        default=3.0,
        help='RMSD threshold for consensus error analysis (default: 3.0 Å)'
    )
    
    advanced_group.add_argument(
        '--reference-frame',
        choices=['full', 'protein', 'dna'],
        default='full',
        help='Reference frame for alignment (default: full structure)'
    )
    
    advanced_group.add_argument(
        '--output-format',
        choices=['csv', 'json', 'both'],
        default='csv',
        help='Output format for data files (default: csv)'
    )
    
    return parser


def validate_arguments(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate and post-process parsed arguments
    
    Args:
        args: Parsed arguments from argparse
        
    Returns:
        Validated and processed arguments
        
    Raises:
        SystemExit: If validation fails
    """
    # Handle legacy vs new interface
    if args.legacy_files:
        if len(args.legacy_files) != 2:
            print("Error: Legacy format requires exactly 2 files (observed predicted)", file=sys.stderr)
            sys.exit(1)
        
        # Convert legacy to new format
        args.experimental = validate_path(args.legacy_files[0], must_exist=True)
        args.predicted = validate_path(args.legacy_files[1], must_exist=True)
        
        # Set default analysis for legacy mode
        if not any([args.all_benchmarks, args.rmsd_only, args.curves, 
                   args.bfactor, args.consensus, args.mutations, args.visualize]):
            args.rmsd_only = True
    
    # Validate experimental/predicted pair
    if args.experimental and not args.predicted:
        print("Error: --predicted/-p is required when using --experimental/-e", file=sys.stderr)
        sys.exit(1)
    
    # Set default analysis if none specified
    if not any([args.all_benchmarks, args.rmsd_only, args.curves, 
               args.bfactor, args.consensus, args.mutations, args.visualize]):
        args.rmsd_only = True
    
    # Handle conflicting options
    if args.verbose and args.quiet:
        print("Error: --verbose and --quiet are mutually exclusive", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Detect available CPU cores for parallel processing
    if args.parallel is None:
        args.parallel = min(4, os.cpu_count() or 1)  # Conservative default
    
    return args


def arg_parser() -> argparse.Namespace:
    """
    Main argument parsing function with full validation
    
    Returns:
        Validated argument namespace ready for use by main()
    """
    parser = create_argument_parser()
    args = parser.parse_args()
    return validate_arguments(args)


# Utility functions for main module integration
def get_analysis_flags(args: argparse.Namespace) -> dict:
    """Extract analysis flags as dictionary for easy checking"""
    return {
        'rmsd': True,  # Always perform basic RMSD
        'rmsd_only': args.rmsd_only,
        'curves': args.curves or args.all_benchmarks,
        'bfactor': args.bfactor or args.all_benchmarks,
        'consensus': args.consensus or args.all_benchmarks,
        'mutations': args.mutations or args.all_benchmarks,
        'visualize': args.visualize or args.all_benchmarks,
        'all_benchmarks': args.all_benchmarks
    }


def get_structure_pairs(args: argparse.Namespace) -> List[Tuple[Path, Path]]:
    """
    Get list of (experimental, predicted) structure file pairs
    
    Args:
        args: Parsed arguments
        
    Returns:
        List of (exp_path, pred_path) tuples
    """
    pairs = []
    
    if args.experimental.is_file() and args.predicted.is_file():
        # Single file pair
        pairs.append((args.experimental, args.predicted))
    elif args.experimental.is_dir() and args.predicted.is_dir():
        # Directory matching
        exp_files = find_structure_files(args.experimental)
        pred_files = find_structure_files(args.predicted)
        
        # Match by filename (without extension)
        exp_stems = {f.stem: f for f in exp_files}
        pred_stems = {f.stem: f for f in pred_files}
        
        common_stems = set(exp_stems.keys()) & set(pred_stems.keys())
        pairs = [(exp_stems[stem], pred_stems[stem]) for stem in sorted(common_stems)]
        
        if not pairs:
            print(f"Warning: No matching structure pairs found between {args.experimental} and {args.predicted}")
    else:
        # Mixed file/directory - not supported
        print("Error: Both experimental and predicted must be files or both must be directories", file=sys.stderr)
        sys.exit(1)
    
    return pairs


# Entry point for testing
if __name__ == "__main__":
    args = arg_parser()
    print(f"Parsed arguments: {args}")
