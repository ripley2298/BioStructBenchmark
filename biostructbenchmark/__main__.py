"""
biostructbenchmark/__main__.py
Main entry point for BioStructBenchmark with multi-frame alignment support
"""

import sys
from pathlib import Path
from typing import Dict, Any, Optional
import json
import traceback


def main():
    """Main entry point with comprehensive error handling and multi-frame support"""
    try:
        # Import CLI module
        from biostructbenchmark.cli import (
            arg_parser, 
            get_analysis_flags, 
            get_structure_pairs,
            get_version,
            print_multi_frame_summary,
            export_multi_frame_results
        )
        
        # Parse arguments
        args = arg_parser()
        analysis_flags = get_analysis_flags(args)
        
        # Setup logging level
        if args.verbose:
            import logging
            logging.basicConfig(level=logging.DEBUG)
        elif args.quiet:
            import logging
            logging.basicConfig(level=logging.ERROR)
        
        # Print header
        if not args.quiet:
            print(f"BioStructBenchmark v{get_version()}")
            print("=" * 70)
            print(f"Experimental: {args.experimental}")
            print(f"Predicted: {args.predicted}")
            print(f"Output: {args.output}")
            print("=" * 70)
        
        # Get structure pairs to process
        structure_pairs = get_structure_pairs(args)
        
        if not structure_pairs:
            print("Error: No valid structure pairs found", file=sys.stderr)
            return 1
        
        # Track overall results
        all_results = []
        failed_pairs = []
        
        # Process each structure pair
        for i, (exp_path, pred_path) in enumerate(structure_pairs, 1):
            if not args.quiet:
                print(f"\n[{i}/{len(structure_pairs)}] Processing: {exp_path.name} vs {pred_path.name}")
            
            # Create output directory for this pair
            pair_name = f"{exp_path.stem}_vs_{pred_path.stem}"
            pair_output = args.output / pair_name
            pair_output.mkdir(parents=True, exist_ok=True)
            
            try:
                # Run analyses based on flags
                pair_results = run_analyses(
                    exp_path, pred_path, pair_output, 
                    analysis_flags, args
                )
                
                all_results.append({
                    'pair': pair_name,
                    'experimental': str(exp_path),
                    'predicted': str(pred_path),
                    'results': pair_results
                })
                
            except Exception as e:
                print(f"Error processing {pair_name}: {e}", file=sys.stderr)
                if args.verbose:
                    traceback.print_exc()
                failed_pairs.append(pair_name)
                continue
        
        # Generate summary report if multiple pairs
        if len(structure_pairs) > 1:
            generate_batch_report(all_results, failed_pairs, args.output)
        
        # Final summary
        if not args.quiet:
            print("\n" + "=" * 70)
            print("ANALYSIS COMPLETE")
            print(f"Processed: {len(all_results)}/{len(structure_pairs)} pairs")
            if failed_pairs:
                print(f"Failed: {', '.join(failed_pairs)}")
            print(f"Results saved to: {args.output}")
            print("=" * 70)
        
        return 0 if not failed_pairs else 1
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        if '--verbose' in sys.argv or '-v' in sys.argv:
            traceback.print_exc()
        return 1


def run_analyses(exp_path: Path, pred_path: Path, output_dir: Path, 
                 flags: Dict[str, bool], args) -> Dict[str, Any]:
    """
    Run all requested analyses on a structure pair
    
    Returns:
        Dictionary of analysis results
    """
    results = {}
    
    # Multi-frame alignment (highest priority)
    if flags.get('multi_frame'):
        results['multi_frame'] = run_multi_frame_alignment(
            exp_path, pred_path, output_dir, args
        )
    
    # Single-frame RMSD (if not doing multi-frame)
    elif flags.get('rmsd') and not flags.get('multi_frame'):
        results['rmsd'] = run_basic_rmsd(
            exp_path, pred_path, output_dir, args
        )
    
    # CURVES+ analysis
    if flags.get('curves'):
        results['curves'] = run_curves_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # B-factor analysis
    if flags.get('bfactor'):
        results['bfactor'] = run_bfactor_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Consensus analysis (if applicable)
    if flags.get('consensus'):
        results['consensus'] = run_consensus_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Mutation analysis
    if flags.get('mutations'):
        results['mutations'] = run_mutation_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Hydrogen bond analysis
    if flags.get('hbond'):
        results['hbond'] = run_hydrogen_bond_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # DSSR DNA structural analysis
    if flags.get('dssr'):
        results['dssr'] = run_dssr_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Visualization
    if flags.get('visualize'):
        results['visualization'] = generate_visualizations(
            results, output_dir, args
        )
    
    return results


def run_multi_frame_alignment(exp_path: Path, pred_path: Path, 
                             output_dir: Path, args) -> Optional[Dict]:
    """Run multi-frame alignment analysis"""
    try:
        from biostructbenchmark.core.alignment import align_structures_three_frames
        from biostructbenchmark.core.io import get_structure
        from biostructbenchmark.cli import (
            print_multi_frame_summary,
            export_multi_frame_results
        )
        
        if not args.quiet:
            print("  → Running multi-frame alignment analysis...")
        
        # Determine if we should save aligned structures
        alignment_output = output_dir / "alignments" if args.save_aligned else None
        
        # Load structures and perform alignment
        observed = get_structure(exp_path)
        predicted = get_structure(pred_path)
        results_dict = align_structures_three_frames(observed, predicted)
        
        # Create expected object structure for CLI compatibility
        class MultiFrameResult:
            def __init__(self, results_dict):
                self.full_structure = results_dict.get('global')
                self.dna_to_protein = results_dict.get('dna_centric')  # DNA positioned relative to protein
                self.dna_to_dna = results_dict.get('protein_centric')  # DNA-to-DNA comparison
                self.results_dict = results_dict  # Keep original for exports
                
            def get_summary(self):
                """Return summary for aggregation"""
                return {
                    'global_rmsd': self.full_structure.overall_rmsd if self.full_structure else None,
                    'dna_positioning_rmsd': self.dna_to_protein.overall_rmsd if self.dna_to_protein else None,
                    'dna_structure_rmsd': self.dna_to_dna.overall_rmsd if self.dna_to_dna else None,
                    'aligned_atoms': self.full_structure.aligned_atom_count if self.full_structure else 0,
                    'residue_count': len(self.full_structure.residue_rmsds) if self.full_structure and self.full_structure.residue_rmsds else 0
                }
        
        result = MultiFrameResult(results_dict)
        
        if result:
            # Print summary if not quiet
            if not args.quiet:
                print_multi_frame_summary(result)
            
            # Export results
            export_multi_frame_results(result, output_dir, args.output_format)
            
            # Return summary for aggregation
            return result.get_summary()
        
        return None
        
    except ImportError as e:
        print(f"  ⚠ Multi-frame alignment not available: {e}")
        return None
    except Exception as e:
        print(f"  ✗ Multi-frame alignment failed: {e}")
        if args.verbose:
            traceback.print_exc()
        return None


def run_basic_rmsd(exp_path: Path, pred_path: Path, 
                   output_dir: Path, args) -> Optional[Dict]:
    """Run basic RMSD analysis"""
    try:
        from biostructbenchmark.core.alignment import compare_structures
        from biostructbenchmark.core.alignment import export_residue_rmsd_csv
        
        if not args.quiet:
            print(f"  → Running {args.reference_frame} frame RMSD analysis...")
        
        result = compare_structures(exp_path, pred_path)
        
        if result:
            # Export if requested
            if args.output_format in ['csv', 'both']:
                export_residue_rmsd_csv(
                    result.residue_rmsds,
                    output_dir / "rmsd_analysis.csv",
                    args.reference_frame
                )
            
            # Return both summary and raw data for visualization
            return {
                'overall_rmsd': result.overall_rmsd,
                'atom_count': result.aligned_atom_count,
                'residue_count': len(result.residue_rmsds),
                '_raw_result': result,  # Include full result object for visualization
                'residue_rmsds': result.residue_rmsds  # Direct access to residue data
            }
        
        return None
        
    except Exception as e:
        print(f"  ✗ RMSD analysis failed: {e}")
        return None


def run_curves_analysis(exp_path: Path, pred_path: Path, 
                       output_dir: Path, args) -> Optional[Dict]:
    """Run CURVES+ DNA geometry analysis"""
    try:
        from biostructbenchmark.analysis.curves import run_curves_analysis
        
        if not args.quiet:
            print("  → Running CURVES+ DNA geometry analysis...")
        
        result = run_curves_analysis(exp_path, pred_path, output_dir)
        
        if result and not args.quiet:
            print("  ✓ CURVES+ analysis complete")
        
        return result
        
    except ImportError:
        if args.verbose:
            print("  ⚠ CURVES+ module not available")
        return None
    except Exception as e:
        print(f"  ✗ CURVES+ analysis failed: {e}")
        return None


def run_bfactor_analysis(exp_path: Path, pred_path: Path, 
                        output_dir: Path, args) -> Optional[Dict]:
    """Run B-factor vs pLDDT analysis"""
    try:
        from biostructbenchmark.analysis.bfactor import analyze_bfactors
        
        if not args.quiet:
            print("  → Analyzing B-factors vs confidence metrics...")
        
        result = analyze_bfactors(exp_path, pred_path, output_dir)
        
        if result and not args.quiet:
            print("  ✓ B-factor analysis complete")
        
        return result
        
    except ImportError:
        if args.verbose:
            print("  ⚠ B-factor module not available")
        return None
    except Exception as e:
        print(f"  ✗ B-factor analysis failed: {e}")
        return None


def run_hydrogen_bond_analysis(exp_path: Path, pred_path: Path,
                              output_dir: Path, args) -> Optional[Dict]:
    """
    Run comprehensive hydrogen bond network analysis for DNA-protein complexes.
    
    This function performs detailed comparison of hydrogen bond networks between
    experimental and predicted structures, focusing on protein-DNA interactions
    critical for binding specificity and structural stability.
    
    Args:
        exp_path (Path): Path to experimental structure file (PDB format)
        pred_path (Path): Path to predicted structure file (PDB/CIF format)
        output_dir (Path): Directory where analysis results will be written
        args: CLI arguments containing analysis parameters
        
    Returns:
        Optional[Dict]: Analysis summary containing:
            - experimental_bonds: Number of H-bonds in experimental structure
            - predicted_bonds: Number of H-bonds in predicted structure  
            - common_bonds: H-bonds present in both structures (shared/conserved)
            - experimental_only_bonds: H-bonds present only in experimental (missing from prediction)
            - predicted_only_bonds: H-bonds present only in predicted (false positives)
            - prediction_accuracy: Fraction of predicted bonds that are correct (0-1)
            - conservation_rate: Fraction of experimental bonds preserved in prediction (0-1)
            - mean_distance_difference: Average distance difference for common bonds (Å)
            - output_dir: Path to detailed results directory
            
    Exports:
        - hydrogen_bonds/: Directory containing detailed H-bond analysis
        - H-bond comparison tables (CSV format)
        - Network similarity metrics
        - Visualization files for H-bond networks
        
    Analysis Details:
        - Distance cutoff: 3.5 Å (configurable)
        - Angle cutoff: 120° for D-H-A geometry (configurable)
        - Focuses on protein-DNA interactions
        - Identifies donor/acceptor atoms for each bond
        - Calculates geometric parameters (distance, angles)
        - Provides statistical comparison metrics
        
    Note:
        Critical for understanding DNA-protein binding specificity and
        validating computational predictions of protein-DNA recognition.
    """
    try:
        from biostructbenchmark.analysis.hbond import HBondAnalyzer
        
        if not args.quiet:
            print("  → Analyzing hydrogen bond networks...")
        
        # Get structural correspondence mapping from alignment system
        from biostructbenchmark.core.alignment import align_structures_three_frames
        from biostructbenchmark.core.io import get_structure
        
        if not args.quiet:
            print("  → Computing structural correspondence for H-bond matching...")
            
        # Load structures and get correspondence mapping
        observed = get_structure(exp_path)
        predicted = get_structure(pred_path)
        alignment_results = align_structures_three_frames(observed, predicted)
        
        # Extract correspondence from global alignment for H-bond analysis
        from biostructbenchmark.core.alignment import create_correspondence_map
        protein_correspondence = create_correspondence_map(observed, predicted, 'protein')
        dna_correspondence = create_correspondence_map(observed, predicted, 'dna')
        
        # Combine correspondence maps
        full_correspondence = {**protein_correspondence, **dna_correspondence}
        
        if not args.quiet:
            print(f"  → Found correspondence for {len(full_correspondence)} residue pairs")
            
        # Initialize analyzer and run analysis with correspondence
        analyzer = HBondAnalyzer()
        comparison, statistics = analyzer.analyze_structures_with_correspondence(
            exp_path, pred_path, full_correspondence)
        
        # Export results
        hbond_output = output_dir / "hydrogen_bonds"
        pair_id = f"{exp_path.stem}_vs_{pred_path.stem}"
        analyzer.export_results(comparison, statistics, hbond_output, pair_id)
        
        # Return summary for visualization using correct attribute names  
        return {
            'experimental_bonds': len(comparison.experimental_bonds),
            'predicted_bonds': len(comparison.predicted_bonds),
            'common_bonds': len(comparison.common_bonds),
            'experimental_only_bonds': len(comparison.experimental_only),
            'predicted_only_bonds': len(comparison.predicted_only),
            'prediction_accuracy': statistics.prediction_accuracy,
            'conservation_rate': statistics.conservation_rate,
            'mean_distance_difference': statistics.mean_distance_difference,
            'output_dir': str(hbond_output)
        }
        
    except ImportError:
        if args.verbose:
            print("  ⚠ Hydrogen bond analysis module not available")
        return None
    except Exception as e:
        print(f"  ✗ Hydrogen bond analysis failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return None

def run_consensus_analysis(exp_path: Path, pred_path: Path, 
                          output_dir: Path, args) -> Optional[Dict]:
    """Run consensus error mapping"""
    try:
        from biostructbenchmark.analysis.consensus import generate_consensus_map
        
        if not args.quiet:
            print("  → Generating consensus error map...")
        
        result = generate_consensus_map(
            exp_path, pred_path, output_dir,
            threshold=args.rmsd_threshold
        )
        
        if result and not args.quiet:
            print("  ✓ Consensus analysis complete")
        
        return result
        
    except ImportError:
        if args.verbose:
            print("  ⚠ Consensus module not available")
        return None
    except Exception as e:
        print(f"  ✗ Consensus analysis failed: {e}")
        return None


def run_mutation_analysis(exp_path: Path, pred_path: Path, 
                         output_dir: Path, args) -> Optional[Dict]:
    """Run mutation impact analysis"""
    try:
        from biostructbenchmark.analysis.mutations import analyze_mutations
        
        if not args.quiet:
            print("  → Analyzing mutation impacts...")
        
        result = analyze_mutations(exp_path, pred_path, output_dir)
        
        if result and not args.quiet:
            print(f"  ✓ Found {result.get('mutation_count', 0)} mutations")
        
        return result
        
    except ImportError:
        if args.verbose:
            print("  ⚠ Mutations module not available")
        return None
    except Exception as e:
        print(f"  ✗ Mutation analysis failed: {e}")
        return None


def run_dssr_analysis(exp_path: Path, pred_path: Path,
                     output_dir: Path, args) -> Optional[Dict]:
    """
    Run X3DNA-DSSR analysis for DNA structural parameters comparison
    
    Focuses on the 5 most critical parameters for protein-DNA binding interface:
    1. Base pairs - Number of canonical base pairs
    2. Helical twist - Average twist per base pair step
    3. Major groove width - Primary protein-binding interface
    4. Minor groove width - Affects DNA bending
    5. Stacking energy - Interface stability indicator
    
    Args:
        exp_path: Path to experimental structure file
        pred_path: Path to predicted structure file
        output_dir: Directory for analysis results
        args: CLI arguments
        
    Returns:
        Dictionary containing DSSR analysis summary or None if failed
    """
    try:
        from biostructbenchmark.analysis.dssr import DSSRAnalyzer
        
        if not args.quiet:
            print("  → Analyzing DNA structural parameters with X3DNA-DSSR...")
            print("    (5 critical parameters for protein-DNA binding interface)")
        
        analyzer = DSSRAnalyzer()
        
        # Analyze both structures
        exp_result = analyzer.analyze_structure(exp_path, "Experimental")
        pred_result = analyzer.analyze_structure(pred_path, "Predicted")
        
        if not exp_result or not pred_result:
            print("  ⚠ DSSR analysis failed - structures could not be processed")
            return None
        
        # Create output directory
        dssr_output = output_dir / "dssr_analysis"
        dssr_output.mkdir(exist_ok=True)
        
        # Export individual results
        import pandas as pd
        results_df = pd.DataFrame([exp_result.to_dict(), pred_result.to_dict()])
        
        # Calculate comparative statistics
        parameters = ['BasePairs', 'Twist(°)', 'MajorGroove(Å)', 
                     'MinorGroove(Å)', 'StackingEnergy(kcal/mol)']
        
        comparison_stats = {}
        for param in parameters:
            # Handle parameter name mapping
            if param == 'Twist(°)':
                exp_val = exp_result.helical_twist
                pred_val = pred_result.helical_twist
            elif param == 'MajorGroove(Å)':
                exp_val = exp_result.major_groove_width
                pred_val = pred_result.major_groove_width
            elif param == 'MinorGroove(Å)':
                exp_val = exp_result.minor_groove_width
                pred_val = pred_result.minor_groove_width
            elif param == 'StackingEnergy(kcal/mol)':
                exp_val = exp_result.stacking_energy
                pred_val = pred_result.stacking_energy
            elif param == 'BasePairs':
                exp_val = exp_result.base_pairs
                pred_val = pred_result.base_pairs
            
            difference = pred_val - exp_val
            
            # Check thresholds
            thresholds = {
                'Twist(°)': 5.0,
                'MajorGroove(Å)': 0.5,
                'MinorGroove(Å)': 0.5,
                'BasePairs': 2.0,
                'StackingEnergy(kcal/mol)': 2.0
            }
            
            flagged = abs(difference) > thresholds.get(param, 0.5)
            
            comparison_stats[param] = {
                'experimental': exp_val,
                'predicted': pred_val,
                'difference': difference,
                'flagged': flagged
            }
        
        # Export results
        results_df.to_csv(dssr_output / "dssr_parameters.csv", index=False)
        
        # Create comparison report
        pair_id = f"{exp_path.stem}_vs_{pred_path.stem}"
        report_path = dssr_output / f"{pair_id}_dssr_comparison.txt"
        
        def format_param_value(value, param):
            """Format parameter values with appropriate significant figures"""
            if 'BasePairs' in param:
                return f"{value:.0f}"
            elif 'Twist' in param or 'Energy' in param:
                return f"{value:.1f}"
            elif 'Groove' in param:
                return f"{value:.2f}"
            else:
                return f"{value:.2f}"
        
        with open(report_path, 'w') as f:
            f.write("X3DNA-DSSR COMPARISON REPORT\n")
            f.write("Critical Parameters for Protein-DNA Binding Interface\n")
            f.write("=" * 60 + "\n\n")
            
            for param, stats in comparison_stats.items():
                status = "⚠️ FLAGGED" if stats['flagged'] else "✓ OK"
                exp_val = format_param_value(stats['experimental'], param)
                pred_val = format_param_value(stats['predicted'], param)
                diff_val = f"{stats['difference']:+.1f}" if 'Twist' in param or 'Energy' in param else f"{stats['difference']:+.2f}"
                
                f.write(f"{param:<20}: {exp_val} → {pred_val} "
                       f"(Δ{diff_val}) {status}\n")
            
            f.write("\nTHRESHOLD CRITERIA:\n")
            f.write("• Twist deviation > 5°\n")
            f.write("• Groove width deviation > 0.5 Å\n")
            f.write("• Base pair count deviation > 2\n")
            f.write("• Stacking energy deviation > 2.0 kcal/mol\n")
        
        if not args.quiet:
            print("  ✓ DSSR analysis complete")
            flagged_count = sum(1 for stats in comparison_stats.values() if stats['flagged'])
            if flagged_count > 0:
                print(f"    ⚠️ {flagged_count}/{len(parameters)} parameters flagged for significant deviations")
            else:
                print(f"    ✓ All {len(parameters)} parameters within acceptable thresholds")
        
        # Return summary for aggregation
        flagged_count = sum(1 for stats in comparison_stats.values() if stats['flagged'])
        return {
            'experimental_base_pairs': exp_result.base_pairs,
            'predicted_base_pairs': pred_result.base_pairs,
            'experimental_twist': exp_result.helical_twist,
            'predicted_twist': pred_result.helical_twist,
            'experimental_major_groove': exp_result.major_groove_width,
            'predicted_major_groove': pred_result.major_groove_width,
            'experimental_minor_groove': exp_result.minor_groove_width,
            'predicted_minor_groove': pred_result.minor_groove_width,
            'experimental_stacking_energy': exp_result.stacking_energy,
            'predicted_stacking_energy': pred_result.stacking_energy,
            'flagged_parameters': flagged_count,
            'output_dir': str(dssr_output)
        }
        
    except ImportError:
        if args.verbose:
            print("  ⚠ X3DNA-DSSR module not available")
        return None
    except Exception as e:
        print(f"  ✗ DSSR analysis failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return None


def generate_visualizations(results: Dict, output_dir: Path, args) -> Optional[Dict]:
    """Generate all requested visualizations"""
    try:
        from biostructbenchmark.visualization.residue_plots import create_residue_analysis, ResidueVisualizer
        
        if not args.quiet:
            print("  → Generating visualizations...")
        
        viz_results = {}
        viz_dir = output_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        # Extract residue data for visualization
        residue_data = None
        analysis_data = {}
        
        # Get residue RMSD data from different result types
        if 'multi_frame' in results and results['multi_frame']:
            # Multi-frame results - we need to access the actual result object
            # This needs to be passed differently from the pipeline
            pass
        elif 'rmsd' in results and results['rmsd'] and results['rmsd'].get('residue_rmsds'):
            # Single RMSD analysis - extract residue data
            residue_data = results['rmsd']['residue_rmsds']
            if not args.quiet:
                print(f"  ✓ Found {len(residue_data)} residues for detailed visualization")
        
        # Generate comprehensive residue analysis if we have residue data
        if residue_data:
            try:
                viz_paths = create_residue_analysis(residue_data, viz_dir, analysis_data)
                viz_results.update(viz_paths)
                
                if not args.quiet:
                    print(f"  ✓ Generated {len(viz_paths)} detailed residue visualizations")
                    
            except Exception as e:
                if not args.quiet:
                    print(f"  ⚠ Detailed residue visualization failed: {e}")
                if args.verbose:
                    import traceback
                    traceback.print_exc()
        else:
            # Fallback message
            if any('rmsd' in str(k).lower() for k in results.keys()):
                if not args.quiet:
                    print("  ℹ RMSD data detected but detailed residue data not available")
                    print("  ℹ Run with updated pipeline for detailed visualizations")
        
        # B-factor data for correlations
        if 'bfactor' in results and results['bfactor']:
            analysis_data['bfactor'] = results['bfactor']
        
        # Consensus data
        if 'consensus' in results and results['consensus']:
            analysis_data['consensus'] = results['consensus']
        
        # Mutation data
        if 'mutations' in results and results['mutations']:
            analysis_data['mutations'] = results['mutations']
        
        # Hydrogen bond data
        if 'hbond' in results and results['hbond']:
            analysis_data['hbond'] = results['hbond']
        
        # Generate basic summary if we have any data
        if results:
            try:
                # Create a simple summary plot instead of complex dashboard
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(10, 6))
                
                # Create summary of available analyses
                analysis_types = list(results.keys())
                analysis_types = [t for t in analysis_types if t != 'visualization']
                
                ax.text(0.1, 0.8, 'BioStructBenchmark Analysis Summary', 
                       fontsize=16, fontweight='bold', transform=ax.transAxes)
                
                y_pos = 0.6
                for i, analysis_type in enumerate(analysis_types):
                    if results[analysis_type]:
                        ax.text(0.1, y_pos - i*0.1, f"✓ {analysis_type.upper()} analysis completed", 
                               fontsize=12, transform=ax.transAxes)
                    else:
                        ax.text(0.1, y_pos - i*0.1, f"✗ {analysis_type.upper()} analysis failed", 
                               fontsize=12, color='red', transform=ax.transAxes)
                
                ax.text(0.1, 0.2, 'Note: Detailed residue visualizations require --save-residue-data flag', 
                       fontsize=10, style='italic', transform=ax.transAxes)
                
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')
                
                summary_path = viz_dir / "analysis_summary.png"
                fig.savefig(summary_path, dpi=300, bbox_inches='tight', facecolor='white')
                plt.close(fig)
                
                viz_results['summary'] = str(summary_path)
                
            except Exception as e:
                if args.verbose:
                    print(f"  ⚠ Could not create summary visualization: {e}")
        
        if not args.quiet:
            if viz_results:
                print(f"  ✓ Generated {len(viz_results)} visualization(s)")
            else:
                print("  ℹ No visualizations generated - detailed data needed")
        
        return viz_results
        
    except ImportError:
        if args.verbose:
            print("  ⚠ Visualization module not available")
        return None
    except Exception as e:
        print(f"  ✗ Visualization failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return None


def generate_batch_report(all_results: list, failed_pairs: list, output_dir: Path):
    """Generate summary report for batch processing"""
    try:
        report = {
            'summary': {
                'total_pairs': len(all_results) + len(failed_pairs),
                'successful': len(all_results),
                'failed': len(failed_pairs)
            },
            'results': all_results,
            'failed_pairs': failed_pairs
        }
        
        # Calculate aggregate statistics
        if all_results:
            rmsds = []
            for result in all_results:
                if 'results' in result:
                    if 'multi_frame' in result['results'] and result['results']['multi_frame']:
                        rmsds.append(result['results']['multi_frame'].get('full_structure_rmsd'))
                    elif 'rmsd' in result['results'] and result['results']['rmsd']:
                        rmsds.append(result['results']['rmsd'].get('overall_rmsd'))
            
            if rmsds:
                import numpy as np
                report['summary']['mean_rmsd'] = float(np.mean([r for r in rmsds if r]))
                report['summary']['std_rmsd'] = float(np.std([r for r in rmsds if r]))
        
        # Save report
        report_path = output_dir / "batch_analysis_report.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\nBatch report saved to: {report_path}")
        
        # Generate DSSR Dataset Summary if DSSR results are available
        dssr_results = []
        for result in all_results:
            if 'results' in result and 'dssr' in result['results'] and result['results']['dssr']:
                dssr_results.append(result)
        
        if dssr_results:
            try:
                from biostructbenchmark.analysis.dssr import generate_dataset_summary_report
                print(f"\n🧬 Generating comprehensive DSSR dataset summary for {len(dssr_results)} structures...")
                summary_stats = generate_dataset_summary_report(dssr_results, output_dir)
                
                # Also print key findings to console
                if 'key_findings' in summary_stats:
                    print(f"\n🚨 Key Findings:")
                    print(f"   • {summary_stats['key_findings']['structures_with_bp_loss']}/{len(dssr_results)} structures show base pair loss")
                    print(f"   • {summary_stats['key_findings']['predicted_no_stacking']}/{len(dssr_results)} predicted structures lack stacking interactions")
                    print(f"   • Most problematic parameter: {summary_stats['key_findings']['most_problematic_param'].replace('_', ' ').title()}")
                
            except ImportError:
                print("⚠️  DSSR summary generation not available")
            except Exception as e:
                print(f"⚠️  Could not generate DSSR summary: {e}")
        
    except Exception as e:
        print(f"Warning: Could not generate batch report: {e}")


if __name__ == "__main__":
    sys.exit(main())