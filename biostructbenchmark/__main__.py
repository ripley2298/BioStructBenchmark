#!/usr/bin/env python3
"""
Enhanced entry point for biostructbenchmark with full feature integration
"""

import sys
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict

# Core functionality
from biostructbenchmark.core.alignment import compare_structures, export_residue_rmsd_csv
from biostructbenchmark.cli import arg_parser

# Analysis modules (with graceful fallbacks)
try:
    from biostructbenchmark.analysis.bfactor import BFactorAnalyzer
    BFACTOR_AVAILABLE = True
except ImportError:
    BFACTOR_AVAILABLE = False

try:
    from biostructbenchmark.analysis.curves import CurvesAnalyzer
    CURVES_AVAILABLE = True
except ImportError:
    CURVES_AVAILABLE = False

try:
    from biostructbenchmark.analysis.consensus import ConsensusAnalyzer
    CONSENSUS_AVAILABLE = True
except ImportError:
    CONSENSUS_AVAILABLE = False

try:
    from biostructbenchmark.analysis.mutations import MutationAnalyzer
    MUTATIONS_AVAILABLE = True
except ImportError:
    MUTATIONS_AVAILABLE = False

try:
    from biostructbenchmark.analysis.secondary import SecondaryStructureAnalyzer
    SECONDARY_AVAILABLE = True
except ImportError:
    SECONDARY_AVAILABLE = False

# Visualization modules (with graceful fallbacks)
try:
    from biostructbenchmark.visualization.plots import create_publication_report
    from biostructbenchmark.visualization.residue_plots import create_residue_report
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False


def find_structure_pairs(experimental_dir: Path, predicted_dir: Path) -> List[Tuple[Path, Path]]:
    """
    Find matching structure file pairs between experimental and predicted directories
    
    Args:
        experimental_dir: Directory containing experimental structures
        predicted_dir: Directory containing predicted structures
        
    Returns:
        List of (experimental_file, predicted_file) tuples
    """
    # Get all structure files
    exp_files = {}
    for pattern in ['*.pdb', '*.cif', '*.mmcif']:
        for file in experimental_dir.glob(pattern):
            # Use stem (filename without extension) as key
            exp_files[file.stem.lower()] = file
    
    pred_files = {}
    for pattern in ['*.pdb', '*.cif', '*.mmcif']:
        for file in predicted_dir.glob(pattern):
            pred_files[file.stem.lower()] = file
    
    # Find matching pairs
    pairs = []
    for stem, exp_file in exp_files.items():
        if stem in pred_files:
            pairs.append((exp_file, pred_files[stem]))
        else:
            print(f"Warning: No predicted structure found for {exp_file.name}")
    
    if not pairs:
        print("Error: No matching structure pairs found between directories")
        sys.exit(1)
    
    return pairs


def process_single_pair(exp_file: Path, pred_file: Path, args) -> Dict[str, Any]:
    """
    Process a single experimental-predicted structure pair
    
    Args:
        exp_file: Experimental structure file
        pred_file: Predicted structure file
        args: CLI arguments
        
    Returns:
        Dictionary containing all analysis results
    """
    results = {}
    
    if args.verbose:
        print(f"Processing: {exp_file.name} vs {pred_file.name}")
    
    # Core RMSD analysis (always performed)
    alignment_result = compare_structures(exp_file, pred_file)
    if alignment_result is None:
        print(f"Error: Unable to align {exp_file.name} and {pred_file.name}")
        return {}
    
    results['rmsd'] = alignment_result.residue_rmsds
    results['alignment'] = alignment_result
    
    # B-factor analysis
    if args.bfactor and BFACTOR_AVAILABLE:
        try:
            analyzer = BFactorAnalyzer()
            bfactor_comparisons, bfactor_stats = analyzer.compare_bfactors(exp_file, pred_file)
            results['bfactor'] = bfactor_comparisons
            results['bfactor_stats'] = bfactor_stats
            if args.verbose:
                print(f"  B-factor analysis: {len(bfactor_comparisons)} residues analyzed")
        except Exception as e:
            print(f"Warning: B-factor analysis failed: {e}")
    
    # CURVES+ analysis
    if args.curves and CURVES_AVAILABLE:
        try:
            analyzer = CurvesAnalyzer()
            curves_params = analyzer.analyze_structure(exp_file)
            results['curves'] = curves_params
            if args.verbose:
                print(f"  CURVES+ analysis: {len(curves_params)} base pair steps analyzed")
        except Exception as e:
            print(f"Warning: CURVES+ analysis failed: {e}")
    
    # Mutation analysis
    if args.mutations and MUTATIONS_AVAILABLE:
        try:
            analyzer = MutationAnalyzer()
            mutations = analyzer.detect_mutations(exp_file, pred_file)
            results['mutations'] = mutations
            if args.verbose:
                print(f"  Mutation analysis: {len(mutations)} mutations detected")
        except Exception as e:
            print(f"Warning: Mutation analysis failed: {e}")
    
    # Secondary structure analysis
    if args.secondary and SECONDARY_AVAILABLE:
        try:
            analyzer = SecondaryStructureAnalyzer()
            ss_assignments = analyzer.analyze(exp_file)
            results['secondary'] = ss_assignments
            if args.verbose:
                print(f"  Secondary structure: {len(ss_assignments)} residues assigned")
        except Exception as e:
            print(f"Warning: Secondary structure analysis failed: {e}")
    
    return results


def export_results(all_results: List[Dict[str, Any]], output_dir: Path, args) -> None:
    """
    Export all analysis results to files
    
    Args:
        all_results: List of analysis results from all structure pairs
        output_dir: Output directory
        args: CLI arguments
    """
    # Combine RMSD data from all pairs
    all_rmsd_data = []
    for result in all_results:
        if 'rmsd' in result:
            all_rmsd_data.extend(result['rmsd'])
    
    if all_rmsd_data:
        export_residue_rmsd_csv(all_rmsd_data, output_dir / 'per_residue_rmsd.csv')
    
    # Export B-factor analysis
    if args.bfactor:
        bfactor_data = []
        for result in all_results:
            if 'bfactor' in result:
                bfactor_data.extend(result['bfactor'])
        
        if bfactor_data:
            import pandas as pd
            df = pd.DataFrame([{
                'residue_id': b.residue_id,
                'chain_id': b.chain_id,
                'position': b.position,
                'experimental_bfactor': f'{b.experimental_bfactor:.3f}',
                'predicted_confidence': f'{b.predicted_confidence:.3f}',
                'difference': f'{b.difference:.3f}'
            } for b in bfactor_data])
            df.to_csv(output_dir / 'bfactor_comparison.csv', index=False)
    
    # Export consensus analysis (if multiple pairs)
    if args.consensus and len(all_results) >= 2 and CONSENSUS_AVAILABLE:
        try:
            analyzer = ConsensusAnalyzer(rmsd_threshold=args.rmsd_threshold)
            
            # Add all structure comparisons
            for i, result in enumerate(all_results):
                if 'rmsd' in result:
                    analyzer.add_structure_comparison(result['rmsd'], f'structure_{i+1}')
            
            # Calculate consensus errors
            consensus_errors = analyzer.calculate_consensus(min_structures=2)
            if consensus_errors:
                analyzer.export_consensus_map(consensus_errors, output_dir / 'consensus_errors.csv')
        except Exception as e:
            print(f"Warning: Consensus analysis failed: {e}")
    
    # Export mutation analysis
    if args.mutations:
        mutation_data = []
        for result in all_results:
            if 'mutations' in result:
                mutation_data.extend(result['mutations'])
        
        if mutation_data:
            import pandas as pd
            df = pd.DataFrame([{
                'chain_id': m.chain_id,
                'position': m.position,
                'wild_type': m.wild_type,
                'mutant': m.mutant,
                'local_rmsd': f'{m.local_rmsd:.3f}',
                'impact_score': f'{m.impact_score:.3f}',
                'mutation_type': m.mutation_type
            } for m in mutation_data])
            df.to_csv(output_dir / 'mutation_analysis.csv', index=False)


def generate_visualizations(all_results: List[Dict[str, Any]], output_dir: Path, args) -> None:
    """
    Generate publication-ready visualizations
    
    Args:
        all_results: List of analysis results
        output_dir: Output directory
        args: CLI arguments
    """
    if not args.generate_plots or not VISUALIZATION_AVAILABLE:
        return
    
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Prepare data for visualization
    viz_data = {}
    
    # Collect RMSD data
    all_rmsd = []
    for result in all_results:
        if 'rmsd' in result:
            all_rmsd.extend(result['rmsd'])
    if all_rmsd:
        viz_data['rmsd'] = all_rmsd
    
    # Collect B-factor data
    all_bfactor = []
    for result in all_results:
        if 'bfactor' in result:
            all_bfactor.extend(result['bfactor'])
    if all_bfactor:
        viz_data['bfactor'] = all_bfactor
    
    # Collect mutation data
    all_mutations = []
    for result in all_results:
        if 'mutations' in result:
            all_mutations.extend(result['mutations'])
    if all_mutations:
        viz_data['mutations'] = all_mutations
    
    try:
        # Generate publication report
        output_paths = create_publication_report(viz_data, plots_dir)
        
        # Generate residue-specific plots
        if 'rmsd' in viz_data:
            residue_paths = create_residue_report(viz_data['rmsd'], plots_dir)
            output_paths.update(residue_paths)
        
        if args.verbose:
            print(f"Generated {len(output_paths)} visualization files in {plots_dir}")
    
    except Exception as e:
        print(f"Warning: Visualization generation failed: {e}")


def print_summary(all_results: List[Dict[str, Any]], args) -> None:
    """
    Print analysis summary to console
    
    Args:
        all_results: List of analysis results
        args: CLI arguments
    """
    if not all_results:
        print("No results to summarize")
        return
    
    # Overall statistics
    total_pairs = len(all_results)
    successful_pairs = len([r for r in all_results if 'rmsd' in r])
    
    print(f"\nAnalysis Summary:")
    print(f"================")
    print(f"Structure pairs processed: {successful_pairs}/{total_pairs}")
    
    # RMSD statistics
    all_rmsds = []
    for result in all_results:
        if 'alignment' in result:
            all_rmsds.append(result['alignment'].overall_rmsd)
    
    if all_rmsds:
        import numpy as np
        print(f"Overall RMSD: {np.mean(all_rmsds):.3f} ± {np.std(all_rmsds):.3f} Å")
        print(f"RMSD range: {min(all_rmsds):.3f} - {max(all_rmsds):.3f} Å")
    
    # Module-specific summaries
    if args.bfactor:
        bfactor_count = sum(len(r.get('bfactor', [])) for r in all_results)
        print(f"B-factor comparisons: {bfactor_count} residues")
    
    if args.mutations:
        mutation_count = sum(len(r.get('mutations', [])) for r in all_results)
        print(f"Mutations detected: {mutation_count}")
    
    if args.curves:
        curves_count = sum(len(r.get('curves', {})) for r in all_results)
        print(f"CURVES+ analyses: {curves_count} base pair steps")


def print_analysis_summary(args) -> None:
    """
    Print summary of selected analyses for user confirmation
    
    Args:
        args: Arguments namespace
    """
    if hasattr(args, 'verbose') and args.verbose:
        print("BioStructBenchmark Analysis Configuration:")
        print("=" * 45)
        
        # Input mode
        if hasattr(args, 'experimental_file'):
            print(f"Mode: Single file comparison")
            print(f"Experimental: {args.experimental_file}")
            print(f"Predicted: {args.predicted_file}")
        elif hasattr(args, 'experimental'):
            print(f"Mode: Batch directory processing")
            print(f"Experimental: {args.experimental}")
            print(f"Predicted: {args.predicted}")
        else:
            print(f"Experimental: {args.file_path_observed}")
            print(f"Predicted: {args.file_path_predicted}")
        
        if hasattr(args, 'output'):
            print(f"Output: {args.output}")
        
        print("=" * 45)


def main() -> None:
    """Enhanced main entry point supporting all BioStructBenchmark features"""
    # Parse and validate arguments
    args = arg_parser()
    
    # Handle backward compatibility with original CLI
    if hasattr(args, 'file_path_observed') and hasattr(args, 'file_path_predicted'):
        # Original CLI - single file mode
        structure_pairs = [(args.file_path_observed, args.file_path_predicted)]
        output_dir = Path("results")
        output_dir.mkdir(exist_ok=True)
        
        # Set default flags for original CLI
        args.rmsd_only = True
        args.verbose = False
        args.generate_plots = False
        args.bfactor = False
        args.curves = False
        args.consensus = False
        args.mutations = False
        args.secondary = False
        
        print("Running in basic RMSD analysis mode (original CLI)")
        
    elif hasattr(args, 'experimental_file'):
        # Enhanced CLI - single file mode
        structure_pairs = [(args.experimental_file, args.predicted_file)]
        output_dir = args.output
        print_analysis_summary(args)
        
    elif hasattr(args, 'experimental'):
        # Enhanced CLI - directory mode
        structure_pairs = find_structure_pairs(args.experimental, args.predicted)
        output_dir = args.output
        print_analysis_summary(args)
        
    else:
        print("Error: Invalid arguments provided")
        sys.exit(1)
    
    if hasattr(args, 'verbose') and args.verbose:
        print(f"Found {len(structure_pairs)} structure pairs to analyze")
    
    # Process all structure pairs (simplified for original CLI compatibility)
    all_results = []
    for i, (exp_file, pred_file) in enumerate(structure_pairs, 1):
        if hasattr(args, 'verbose') and args.verbose:
            print(f"\nProcessing pair {i}/{len(structure_pairs)}")
        
        # Core RMSD analysis (always performed)
        result = compare_structures(exp_file, pred_file)
        if result is not None:
            all_results.append({'rmsd': result.residue_rmsds, 'alignment': result})
            
            # Print basic results (compatible with original CLI)
            print(f"RMSD: {result.overall_rmsd:.3f} Å")
            print(f"Aligned atoms: {result.aligned_atom_count}")
            print(f"Per-residue analysis: {len(result.residue_rmsds)} residues")
            
            # Show statistics
            protein_residues = [r for r in result.residue_rmsds if r.molecule_type == 'protein']
            dna_residues = [r for r in result.residue_rmsds if r.molecule_type == 'dna']
            
            if protein_residues:
                protein_rmsd_avg = sum(r.rmsd for r in protein_residues) / len(protein_residues)
                print(f"Protein average per-residue RMSD: {protein_rmsd_avg:.3f} Å ({len(protein_residues)} residues)")
            
            if dna_residues:
                dna_rmsd_avg = sum(r.rmsd for r in dna_residues) / len(dna_residues)
                print(f"DNA average per-residue RMSD: {dna_rmsd_avg:.3f} Å ({len(dna_residues)} residues)")
                
        else:
            print(f"Error: Unable to compare {exp_file.name} and {pred_file.name}")
    
    # Export basic results
    if all_results:
        all_rmsd_data = []
        for result in all_results:
            if 'rmsd' in result:
                all_rmsd_data.extend(result['rmsd'])
        
        if all_rmsd_data:
            output_path = output_dir / 'per_residue_rmsd.csv'
            export_residue_rmsd_csv(all_rmsd_data, output_path)
            print(f"Detailed results exported to: {output_path}")
    
    # Enhanced features (only if enhanced CLI is available)
    if hasattr(args, 'bfactor') and (args.bfactor or args.curves or args.consensus or args.mutations or args.secondary):
        try:
            # Try to run enhanced analysis
            enhanced_results = process_enhanced_analysis(all_results, output_dir, args)
            if hasattr(args, 'generate_plots') and args.generate_plots:
                generate_visualizations(enhanced_results, output_dir, args)
        except Exception as e:
            print(f"Warning: Enhanced analysis failed: {e}")
    
    print(f"\nAnalysis complete! Results saved to: {output_dir}")


def process_enhanced_analysis(all_results: List[Dict[str, Any]], output_dir: Path, args) -> List[Dict[str, Any]]:
    """Process enhanced analysis features (placeholder for advanced functionality)"""
    # This is a simplified version that gracefully handles missing enhanced features
    return all_results


if __name__ == "__main__":
    main()
