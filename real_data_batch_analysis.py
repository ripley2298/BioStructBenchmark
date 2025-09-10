#!/usr/bin/env python3
"""
Real Data Batch Analysis for BioStructBenchmark
Processes experimental vs AlphaFold3 predicted structures
"""

import sys
import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Tuple
import traceback

# Add project root to path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

def find_structure_pairs(experimental_dir: Path, predicted_dir: Path) -> List[Tuple[Path, Path, str]]:
    """Find matching experimental and predicted structure pairs"""
    pairs = []
    
    # Get experimental files
    exp_files = list(experimental_dir.glob("*.pdb"))
    exp_files.extend(list(experimental_dir.glob("*.cif")))
    
    print(f"Found {len(exp_files)} experimental structures")
    
    for exp_file in exp_files:
        # Extract base identifier (e.g., p456_02 from p456_02_experimental.pdb)
        base_name = exp_file.stem.replace('_experimental', '')
        
        # Look for matching predicted file with different naming patterns
        predicted_candidates = [
            predicted_dir / f"{base_name}_alphafold3.cif",
            predicted_dir / f"{base_name}_alphafold.cif", 
            predicted_dir / f"{base_name}_predicted.cif",
            predicted_dir / f"{base_name}.cif"
        ]
        
        pred_file = None
        for candidate in predicted_candidates:
            if candidate.exists():
                pred_file = candidate
                break
        
        if pred_file:
            pairs.append((exp_file, pred_file, base_name))
            print(f"✓ Paired: {exp_file.name} <-> {pred_file.name}")
        else:
            print(f"⚠ No match for {exp_file.name}")
    
    print(f"\nTotal structure pairs found: {len(pairs)}")
    return pairs

def run_batch_analysis(structure_pairs: List[Tuple[Path, Path, str]], output_dir: Path):
    """Run comprehensive batch analysis on structure pairs"""
    
    from biostructbenchmark.core.alignment import align_structures_three_frames
    from biostructbenchmark.core.io import get_structure
    from biostructbenchmark.analysis.pca import PCAAnalyzer
    from biostructbenchmark.visualization.residue_plots import PublicationPlotter, create_publication_report
    
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Results storage
    batch_results = []
    detailed_results = {}
    failed_analyses = []
    
    # Process each structure pair
    print(f"\n=== PROCESSING {len(structure_pairs)} STRUCTURE PAIRS ===")
    
    for i, (exp_path, pred_path, pair_id) in enumerate(structure_pairs, 1):
        print(f"\n[{i}/{len(structure_pairs)}] Analyzing {pair_id}...")
        
        try:
            # Load structures
            exp_structure = get_structure(exp_path)
            pred_structure = get_structure(pred_path)
            
            if not exp_structure or not pred_structure:
                failed_analyses.append((pair_id, "Structure loading failed"))
                print(f"  ✗ Failed to load structures")
                continue
            
            # Perform three-frame alignment
            alignment_results = align_structures_three_frames(exp_structure, pred_structure)
            
            if not alignment_results:
                failed_analyses.append((pair_id, "Alignment failed"))
                print(f"  ✗ Alignment failed")
                continue
            
            # Extract results for batch summary
            result_summary = {
                'pair_id': pair_id,
                'experimental_file': exp_path.name,
                'predicted_file': pred_path.name,
                'analysis_timestamp': datetime.now().isoformat()
            }
            
            # Add alignment results
            for frame_name, frame_result in alignment_results.items():
                result_summary[f'{frame_name}_rmsd'] = round(float(frame_result.overall_rmsd), 4)
                result_summary[f'{frame_name}_atoms'] = int(getattr(frame_result, 'aligned_atom_count', 0))
                
                # Include residue data for visualization (especially for global frame)
                if frame_name == 'global' and hasattr(frame_result, 'residue_rmsds') and frame_result.residue_rmsds:
                    result_summary['residue_rmsds'] = frame_result.residue_rmsds
                
                # Quality assessment
                rmsd = frame_result.overall_rmsd
                if rmsd < 2.0:
                    quality = "excellent"
                elif rmsd < 3.0:
                    quality = "good"
                elif rmsd < 4.0:
                    quality = "moderate"
                else:
                    quality = "poor"
                result_summary[f'{frame_name}_quality'] = quality
            
            batch_results.append(result_summary)
            detailed_results[pair_id] = alignment_results
            
            # Print summary
            global_rmsd = alignment_results['global'].overall_rmsd
            print(f"  ✓ Global RMSD: {global_rmsd:.3f} Å ({result_summary['global_quality']})")
            
            # Create individual structure report
            pair_dir = output_dir / f"{pair_id}_analysis"
            pair_dir.mkdir(exist_ok=True)
            
            # Export per-residue RMSD data for each frame
            from biostructbenchmark.core.alignment import export_residue_rmsd_csv, save_aligned_structures
            for frame_name, frame_result in alignment_results.items():
                if hasattr(frame_result, 'residue_rmsds') and frame_result.residue_rmsds:
                    csv_path = pair_dir / f"rmsd_{frame_name}.csv"
                    export_residue_rmsd_csv(frame_result.residue_rmsds, str(csv_path))
            
            # Save aligned structures used for RMSD calculations
            print(f"  📁 Saving aligned structures...")
            save_aligned_structures(alignment_results, pair_dir, pair_id, exp_structure, pred_structure)
            
            # Create individual pair summary
            pair_summary = {
                'pair_id': pair_id,
                'files': {
                    'experimental': str(exp_path),
                    'predicted': str(pred_path)
                },
                'alignment_summary': {
                    frame: {
                        'rmsd': float(result.overall_rmsd),
                        'atoms': int(getattr(result, 'aligned_atom_count', 0)),
                        'quality': result_summary[f'{frame}_quality']
                    }
                    for frame, result in alignment_results.items()
                },
                'overall_assessment': {
                    'primary_rmsd': float(global_rmsd),
                    'quality_category': result_summary['global_quality'],
                    'suitable_for_analysis': bool(global_rmsd < 4.0)
                }
            }
            
            with open(pair_dir / f"{pair_id}_summary.json", 'w') as f:
                json.dump(pair_summary, f, indent=2)
            
        except Exception as e:
            failed_analyses.append((pair_id, str(e)))
            print(f"  ✗ Error: {e}")
            continue
    
    print(f"\n=== BATCH ANALYSIS COMPLETE ===")
    print(f"Successful analyses: {len(batch_results)}")
    print(f"Failed analyses: {len(failed_analyses)}")
    
    # Generate batch summary files
    if batch_results:
        generate_batch_summary(batch_results, failed_analyses, output_dir)
        
        # Generate advanced analytics if we have enough data
        if len(batch_results) >= 3:
            generate_advanced_analytics(batch_results, detailed_results, output_dir)
        
        # Generate visualizations
        generate_batch_visualizations(batch_results, output_dir)
        
        return True
    else:
        print("No successful analyses to report")
        return False

def generate_batch_summary(batch_results: List[Dict], failed_analyses: List, output_dir: Path):
    """Generate comprehensive batch summary reports"""
    
    print(f"\n=== GENERATING BATCH SUMMARY ===")
    
    # Create comprehensive CSV summary
    batch_df = pd.DataFrame(batch_results)
    summary_csv = output_dir / "batch_analysis_summary.csv"
    batch_df.to_csv(summary_csv, index=False)
    print(f"✓ Generated {summary_csv.name}")
    
    # Calculate batch statistics
    global_rmsds = [r['global_rmsd'] for r in batch_results]
    dna_rmsds = [r.get('dna_centric_rmsd', np.nan) for r in batch_results if 'dna_centric_rmsd' in r]
    protein_rmsds = [r.get('protein_centric_rmsd', np.nan) for r in batch_results if 'protein_centric_rmsd' in r]
    
    # Remove NaN values
    dna_rmsds = [x for x in dna_rmsds if not np.isnan(x)]
    protein_rmsds = [x for x in protein_rmsds if not np.isnan(x)]
    
    batch_stats = {
        "analysis_metadata": {
            "analysis_date": datetime.now().isoformat(),
            "total_structure_pairs": len(batch_results) + len(failed_analyses),
            "successful_analyses": len(batch_results),
            "failed_analyses": len(failed_analyses),
            "success_rate": f"{len(batch_results)/(len(batch_results) + len(failed_analyses))*100:.1f}%"
        },
        "rmsd_statistics": {
            "global_alignment": {
                "mean": float(np.mean(global_rmsds)),
                "std": float(np.std(global_rmsds)),
                "min": float(np.min(global_rmsds)),
                "max": float(np.max(global_rmsds)),
                "median": float(np.median(global_rmsds)),
                "count": len(global_rmsds)
            }
        },
        "quality_distribution": {
            "excellent": sum(1 for r in batch_results if r['global_quality'] == 'excellent'),
            "good": sum(1 for r in batch_results if r['global_quality'] == 'good'),
            "moderate": sum(1 for r in batch_results if r['global_quality'] == 'moderate'),
            "poor": sum(1 for r in batch_results if r['global_quality'] == 'poor')
        }
    }
    
    # Add DNA and protein statistics if available
    if dna_rmsds:
        batch_stats["rmsd_statistics"]["dna_centric"] = {
            "mean": float(np.mean(dna_rmsds)),
            "std": float(np.std(dna_rmsds)),
            "min": float(np.min(dna_rmsds)),
            "max": float(np.max(dna_rmsds)),
            "count": len(dna_rmsds)
        }
    
    if protein_rmsds:
        batch_stats["rmsd_statistics"]["protein_centric"] = {
            "mean": float(np.mean(protein_rmsds)),
            "std": float(np.std(protein_rmsds)),
            "min": float(np.min(protein_rmsds)),
            "max": float(np.max(protein_rmsds)),
            "count": len(protein_rmsds)
        }
    
    # Add failure analysis
    if failed_analyses:
        batch_stats["failed_analyses"] = [
            {"pair_id": pair_id, "error": error} 
            for pair_id, error in failed_analyses
        ]
    
    # Add interpretation
    mean_rmsd = batch_stats["rmsd_statistics"]["global_alignment"]["mean"]
    if mean_rmsd < 2.0:
        overall_quality = "excellent"
        interpretation = "AlphaFold3 predictions show excellent agreement with experimental structures"
    elif mean_rmsd < 3.0:
        overall_quality = "good"
        interpretation = "AlphaFold3 predictions show good agreement with experimental structures"
    elif mean_rmsd < 4.0:
        overall_quality = "moderate"
        interpretation = "AlphaFold3 predictions show moderate agreement with experimental structures"
    else:
        overall_quality = "poor"
        interpretation = "AlphaFold3 predictions show significant deviations from experimental structures"
    
    batch_stats["interpretation"] = {
        "overall_quality": overall_quality,
        "summary": interpretation,
        "recommendation": "Continue analysis with individual structure investigation" if overall_quality in ["good", "excellent"] else "Investigate systematic prediction errors"
    }
    
    # Save batch statistics
    stats_json = output_dir / "batch_statistics.json"
    with open(stats_json, 'w') as f:
        json.dump(batch_stats, f, indent=2)
    print(f"✓ Generated {stats_json.name}")
    
    # Print summary to console
    print(f"\n=== BATCH ANALYSIS SUMMARY ===")
    print(f"Total Structures: {len(batch_results)}")
    print(f"Mean Global RMSD: {mean_rmsd:.3f} Å")
    print(f"Overall Quality: {overall_quality.title()}")
    print(f"Quality Distribution:")
    for quality, count in batch_stats["quality_distribution"].items():
        if count > 0:
            print(f"  - {quality.title()}: {count} structures")

def generate_advanced_analytics(batch_results: List[Dict], detailed_results: Dict, output_dir: Path):
    """Generate advanced analytics including PCA if sufficient data"""
    
    print(f"\n=== GENERATING ADVANCED ANALYTICS ===")
    
    try:
        from biostructbenchmark.analysis.pca import PCAAnalyzer
        from biostructbenchmark.visualization.pca_plots import PCAVisualizer
        
        # Prepare data for PCA analysis
        structure_data = {}
        for result in batch_results:
            pair_id = result['pair_id']
            structure_data[pair_id] = {
                'global_rmsd': result['global_rmsd'],
                'global_atoms': result['global_atoms']
            }
            
            # Add other metrics if available
            if 'dna_centric_rmsd' in result:
                structure_data[pair_id]['dna_centric_rmsd'] = result['dna_centric_rmsd']
            if 'protein_centric_rmsd' in result:
                structure_data[pair_id]['protein_centric_rmsd'] = result['protein_centric_rmsd']
        
        # Perform PCA analysis
        analyzer = PCAAnalyzer()
        pca_result = analyzer.perform_structure_pca(structure_data)
        outliers = analyzer.identify_structure_outliers(pca_result, list(structure_data.keys()))
        
        # Generate PCA visualizations
        visualizer = PCAVisualizer()
        pca_plots = visualizer.create_pca_summary_plots(
            pca_result, list(structure_data.keys()), outliers, 
            output_dir / "pca_analysis"
        )
        
        print(f"✓ Generated PCA analysis with {len(pca_plots)} plots")
        
        # Create PCA summary report
        pca_report = {
            "pca_analysis": {
                "total_structures": len(structure_data),
                "principal_components": len(pca_result.explained_variance_ratio),
                "variance_explained": {
                    f"PC{i+1}": float(ratio) 
                    for i, ratio in enumerate(pca_result.explained_variance_ratio)
                },
                "outliers_detected": len(outliers),
                "outlier_details": [
                    {
                        "structure_id": o.structure_id,
                        "outlier_type": o.outlier_type,
                        "outlier_score": float(o.outlier_score)
                    }
                    for o in outliers
                ]
            }
        }
        
        pca_report_path = output_dir / "pca_analysis_report.json"
        with open(pca_report_path, 'w') as f:
            json.dump(pca_report, f, indent=2)
        print(f"✓ Generated {pca_report_path.name}")
        
    except Exception as e:
        print(f"⚠ Advanced analytics generation failed: {e}")

def generate_batch_visualizations(batch_results: List[Dict], output_dir: Path):
    """Generate comprehensive batch visualization reports using unified residue analysis"""
    
    print(f"\n=== GENERATING VISUALIZATIONS ===")
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        from biostructbenchmark.visualization.residue_plots import create_residue_analysis, ResidueVisualizer
        
        viz_dir = output_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        # Extract all residue data for comprehensive analysis
        all_residue_data = []
        global_rmsds = []
        
        for result in batch_results:
            if 'residue_rmsds' in result and result['residue_rmsds']:
                all_residue_data.extend(result['residue_rmsds'])
            if 'global_rmsd' in result:
                global_rmsds.append(result['global_rmsd'])
        
        # Generate comprehensive residue analysis if we have residue data
        if all_residue_data:
            print(f"  → Generating comprehensive residue analysis for {len(all_residue_data)} residues...")
            try:
                residue_viz_paths = create_residue_analysis(all_residue_data, viz_dir)
                print(f"  ✓ Generated {len(residue_viz_paths)} detailed residue visualizations")
            except Exception as e:
                print(f"  ⚠ Detailed residue visualization failed: {e}")
        
        # Generate batch summary dashboard
        
        plt.figure(figsize=(12, 8))
        
        # Subplot 1: RMSD distribution
        plt.subplot(2, 2, 1)
        plt.hist(global_rmsds, bins=15, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(np.mean(global_rmsds), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(global_rmsds):.3f} Å')
        plt.xlabel('Global RMSD (Å)')
        plt.ylabel('Frequency')
        plt.title('RMSD Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Subplot 2: Quality pie chart
        plt.subplot(2, 2, 2)
        quality_counts = {}
        for result in batch_results:
            qual = result['global_quality']
            quality_counts[qual] = quality_counts.get(qual, 0) + 1
        
        colors = {'excellent': '#28a745', 'good': '#17a2b8', 'moderate': '#ffc107', 'poor': '#dc3545'}
        pie_colors = [colors.get(q, 'gray') for q in quality_counts.keys()]
        
        plt.pie(quality_counts.values(), labels=quality_counts.keys(), autopct='%1.1f%%',
                colors=pie_colors)
        plt.title('Quality Distribution')
        
        # Subplot 3: RMSD vs Structure ID
        plt.subplot(2, 2, 3)
        structure_ids = [r['pair_id'] for r in batch_results]
        plt.scatter(range(len(global_rmsds)), global_rmsds, alpha=0.7)
        plt.axhline(y=2.0, color='green', linestyle='--', alpha=0.7, label='Excellent threshold')
        plt.axhline(y=3.0, color='orange', linestyle='--', alpha=0.7, label='Good threshold')
        plt.xlabel('Structure Index')
        plt.ylabel('Global RMSD (Å)')
        plt.title('RMSD by Structure')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Subplot 4: Statistics table
        plt.subplot(2, 2, 4)
        plt.axis('off')
        
        stats_text = f"""
Batch Analysis Summary

Total Structures: {len(batch_results)}
Mean RMSD: {np.mean(global_rmsds):.3f} ± {np.std(global_rmsds):.3f} Å
Min RMSD: {np.min(global_rmsds):.3f} Å  
Max RMSD: {np.max(global_rmsds):.3f} Å
Median RMSD: {np.median(global_rmsds):.3f} Å

Quality Breakdown:
{chr(10).join([f"  {q.title()}: {c} structures" for q, c in quality_counts.items()])}
"""
        plt.text(0.1, 0.9, stats_text, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        dashboard_path = viz_dir / "batch_analysis_dashboard.png"
        plt.savefig(dashboard_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Generated {dashboard_path.name}")
        
    except Exception as e:
        print(f"⚠ Visualization generation failed: {e}")

def main():
    """Main batch analysis execution"""
    
    print("BioStructBenchmark Real Data Batch Analysis")
    print("=" * 50)
    
    # Define paths
    experimental_dir = Path("/Users/mesler/Documents/IonuI_StructureML/experimental")
    predicted_dir = Path("/Users/mesler/Documents/IonuI_StructureML/predicted_alphafold3")
    output_dir = PROJECT_ROOT / "tests" / "outputs" / "real_data_batch_analysis"
    
    # Verify directories exist
    if not experimental_dir.exists():
        print(f"Error: Experimental directory not found: {experimental_dir}")
        return False
    
    if not predicted_dir.exists():
        print(f"Error: Predicted directory not found: {predicted_dir}")
        return False
    
    print(f"Experimental structures: {experimental_dir}")
    print(f"Predicted structures: {predicted_dir}")
    print(f"Output directory: {output_dir}")
    
    # Find structure pairs
    structure_pairs = find_structure_pairs(experimental_dir, predicted_dir)
    
    if not structure_pairs:
        print("No matching structure pairs found!")
        return False
    
    # Run batch analysis
    success = run_batch_analysis(structure_pairs, output_dir)
    
    if success:
        print(f"\n🎉 Batch analysis complete!")
        print(f"📁 Results saved to: {output_dir}")
        print(f"📊 Generated comprehensive reports and visualizations")
    else:
        print(f"\n⚠ Batch analysis failed - check error messages above")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)