"""
Analysis runner functions for larger/complex analyses
Contains runners that would make their target modules too large
"""

from pathlib import Path
from typing import Dict, Optional
import traceback


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
                self.dna_to_protein = results_dict.get('dna_centric')
                self.dna_to_dna = results_dict.get('protein_centric')
                self.results_dict = results_dict
                
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
                '_raw_result': result,
                'residue_rmsds': result.residue_rmsds
            }
        
        return None
        
    except Exception as e:
        print(f"  ✗ RMSD analysis failed: {e}")
        return None



def run_hydrogen_bond_analysis(exp_path: Path, pred_path: Path,
                              output_dir: Path, args) -> Optional[Dict]:
    """Run comprehensive hydrogen bond network analysis for DNA-protein complexes"""
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

        # Convert tuple-based correspondence to DSSR-format strings for H-bond analysis
        # Tuple format: (('A', (' ', 8, ' ')), ('A', (' ', 8, ' ')))
        # DSSR format: "A.PRO8" -> "A.PRO8"
        dssr_correspondence = {}
        for obs_key, pred_key in full_correspondence.items():
            obs_chain_id, obs_res_id = obs_key
            pred_chain_id, pred_res_id = pred_key

            # Get residue objects to extract residue names
            obs_res = observed[0][obs_chain_id][obs_res_id]
            pred_res = predicted[0][pred_chain_id][pred_res_id]

            # Create DSSR-format keys: "CHAIN.RESNAMERESNUM"
            obs_dssr_key = f"{obs_chain_id}.{obs_res.get_resname().strip()}{obs_res_id[1]}"
            pred_dssr_key = f"{pred_chain_id}.{pred_res.get_resname().strip()}{pred_res_id[1]}"

            dssr_correspondence[obs_dssr_key] = pred_dssr_key

        if not args.quiet:
            print(f"  → Found correspondence for {len(dssr_correspondence)} residue pairs")

        # Initialize analyzer and run analysis with correspondence
        analyzer = HBondAnalyzer()
        comparison, statistics = analyzer.analyze_structures_with_correspondence(
            exp_path, pred_path, dssr_correspondence)
        
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
            traceback.print_exc()
        return None


def run_dssr_analysis(exp_path: Path, pred_path: Path,
                     output_dir: Path, args) -> Optional[Dict]:
    """Run X3DNA-DSSR analysis for DNA structural parameters comparison"""
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
            traceback.print_exc()
        return None


def run_mutant_pca_analysis(exp_path: Path, pred_path: Path, 
                           output_dir: Path, args) -> Optional[Dict]:
    """Run mutant PCA error analysis for systematic modeling errors"""
    try:
        from biostructbenchmark.analysis.mutant_error_pipeline import MutantErrorPipeline
        
        if not args.quiet:
            print("  → Running mutant PCA error analysis...")
        
        # For mutant PCA, we need multiple mutant structures
        # Use the experimental as wildtype reference and predicted as mutant
        wildtype_pdb = str(exp_path)
        mutant_pdb_dict = {"predicted": str(pred_path)}
        
        # Check if we have multiple structure pairs to work with
        if hasattr(args, 'structure_pairs') and len(args.structure_pairs) > 1:
            # Use first pair as wildtype, rest as mutants
            pairs = args.structure_pairs
            wildtype_pdb = str(pairs[0][0])  # First experimental as wildtype
            mutant_pdb_dict = {
                f"mutant_{i}": str(pair[1]) 
                for i, pair in enumerate(pairs)
            }
        
        # Initialize and run pipeline
        pipeline = MutantErrorPipeline(wildtype_pdb, mutant_pdb_dict)
        
        # Create output subdirectory
        mutant_output_dir = output_dir / "mutant_pca_analysis"
        mutant_output_dir.mkdir(exist_ok=True)
        
        # Change to output directory for pipeline execution
        import os
        original_cwd = os.getcwd()
        os.chdir(mutant_output_dir)
        
        try:
            pipeline.run_pipeline()
            
            # Return summary for aggregation
            return {
                'analysis_type': 'mutant_pca',
                'wildtype_structure': wildtype_pdb,
                'mutant_count': len(mutant_pdb_dict),
                'output_files': [
                    'error_hotspots.csv',
                    'pca_scatter_plot.png', 
                    'top_residues_pc1_loading.png'
                ],
                'success': True
            }
        finally:
            os.chdir(original_cwd)
            
    except Exception as e:
        if not args.quiet:
            print(f"    ⚠ Mutant PCA analysis failed: {e}")
        if args.verbose:
            print(f"    Traceback: {traceback.format_exc()}")
        return None