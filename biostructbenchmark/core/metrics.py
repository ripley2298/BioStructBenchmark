"""Error decomposition and analysis metrics"""

import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass
from biostructbenchmark.core.alignment import ResidueRMSD, AlignmentResult


@dataclass
class ErrorComponents:
    """Container for decomposed structural errors"""
    translation_error: float
    orientation_error: float  
    total_error: float
    translation_vector: np.ndarray
    rotation_angle: float


@dataclass
class StructureMetrics:
    """Comprehensive structure comparison metrics"""
    overall_rmsd: float
    protein_rmsd: float
    dna_rmsd: float
    error_components: ErrorComponents
    residue_error_stats: Dict[str, float]
    worst_residues: List[Tuple[str, float]]  # (residue_id, rmsd)
    best_residues: List[Tuple[str, float]]


def decompose_structural_error(alignment_result: AlignmentResult) -> ErrorComponents:
    """
    Decompose structural error into translation and orientation components.
    
    Args:
        alignment_result: Results from structure alignment
        
    Returns:
        ErrorComponents with translation and orientation errors
    """
    
    # Extract transformation components
    rotation_matrix = alignment_result.rotation_matrix
    translation_vector = alignment_result.translation_vector
    
    # Calculate translation error magnitude
    translation_error = np.linalg.norm(translation_vector)
    
    # Calculate rotation angle from rotation matrix
    # Using Rodrigues' rotation formula: trace(R) = 1 + 2*cos(θ)
    trace_r = np.trace(rotation_matrix)
    rotation_angle = np.arccos(np.clip((trace_r - 1) / 2, -1, 1))
    
    # Convert to degrees
    rotation_angle_deg = np.degrees(rotation_angle)
    
    # Estimate orientation error contribution
    # This is a simplified approximation - in practice you'd want more sophisticated analysis
    orientation_error = rotation_angle_deg * 0.1  # Rough conversion to Angstrom-like units
    
    return ErrorComponents(
        translation_error=translation_error,
        orientation_error=orientation_error,
        total_error=alignment_result.overall_rmsd,
        translation_vector=translation_vector,
        rotation_angle=rotation_angle_deg
    )


def calculate_residue_statistics(residue_rmsds: List[ResidueRMSD]) -> Dict[str, float]:
    """
    Calculate statistical summary of per-residue RMSDs.
    
    Args:
        residue_rmsds: List of per-residue RMSD values
        
    Returns:
        Dictionary with statistical measures
    """
    if not residue_rmsds:
        return {}
    
    rmsd_values = [r.rmsd for r in residue_rmsds]
    
    return {
        'mean_rmsd': np.mean(rmsd_values),
        'median_rmsd': np.median(rmsd_values),
        'std_rmsd': np.std(rmsd_values),
        'min_rmsd': np.min(rmsd_values),
        'max_rmsd': np.max(rmsd_values),
        'q25_rmsd': np.percentile(rmsd_values, 25),
        'q75_rmsd': np.percentile(rmsd_values, 75),
        'rmsd_range': np.max(rmsd_values) - np.min(rmsd_values)
    }


def identify_outlier_residues(residue_rmsds: List[ResidueRMSD], 
                            n_worst: int = 5, 
                            n_best: int = 5) -> Tuple[List[Tuple[str, float]], List[Tuple[str, float]]]:
    """
    Identify the best and worst residues by RMSD.
    
    Args:
        residue_rmsds: List of per-residue RMSD values
        n_worst: Number of worst residues to return
        n_best: Number of best residues to return
        
    Returns:
        Tuple of (worst_residues, best_residues) as (residue_id, rmsd) pairs
    """
    if not residue_rmsds:
        return [], []
    
    # Sort by RMSD
    sorted_residues = sorted(residue_rmsds, key=lambda x: x.rmsd)
    
    # Get best and worst
    best_residues = [(r.residue_id, r.rmsd) for r in sorted_residues[:n_best]]
    worst_residues = [(r.residue_id, r.rmsd) for r in sorted_residues[-n_worst:]]
    worst_residues.reverse()  # Worst first
    
    return worst_residues, best_residues


def calculate_molecule_specific_rmsd(residue_rmsds: List[ResidueRMSD]) -> Tuple[float, float]:
    """
    Calculate separate RMSD values for protein and DNA components.
    
    Args:
        residue_rmsds: List of per-residue RMSD values
        
    Returns:
        Tuple of (protein_rmsd, dna_rmsd)
    """
    protein_residues = [r for r in residue_rmsds if r.molecule_type == 'protein']
    dna_bases = [r for r in residue_rmsds if r.molecule_type == 'dna']
    
    protein_rmsd = 0.0
    if protein_residues:
        protein_rmsd = np.mean([r.rmsd for r in protein_residues])
    
    dna_rmsd = 0.0  
    if dna_bases:
        dna_rmsd = np.mean([r.rmsd for r in dna_bases])
    
    return protein_rmsd, dna_rmsd


def generate_comprehensive_metrics(alignment_result: AlignmentResult) -> StructureMetrics:
    """
    Generate comprehensive structural comparison metrics.
    
    Args:
        alignment_result: Results from structure alignment
        
    Returns:
        StructureMetrics with all calculated metrics
    """
    # Decompose errors
    error_components = decompose_structural_error(alignment_result)
    
    # Calculate residue statistics
    residue_stats = calculate_residue_statistics(alignment_result.residue_rmsds)
    
    # Identify outliers
    worst_residues, best_residues = identify_outlier_residues(alignment_result.residue_rmsds)
    
    # Calculate molecule-specific RMSDs
    protein_rmsd, dna_rmsd = calculate_molecule_specific_rmsd(alignment_result.residue_rmsds)
    
    return StructureMetrics(
        overall_rmsd=alignment_result.overall_rmsd,
        protein_rmsd=protein_rmsd,
        dna_rmsd=dna_rmsd,
        error_components=error_components,
        residue_error_stats=residue_stats,
        worst_residues=worst_residues,
        best_residues=best_residues
    )


def print_metrics_summary(metrics: StructureMetrics) -> None:
    """Print a formatted summary of structure metrics."""
    print("\n" + "="*50)
    print("STRUCTURE COMPARISON METRICS")
    print("="*50)
    
    print(f"\nOverall RMSD: {metrics.overall_rmsd:.3f} Å")
    print(f"Protein RMSD: {metrics.protein_rmsd:.3f} Å")
    print(f"DNA RMSD: {metrics.dna_rmsd:.3f} Å")
    
    print(f"\nError Decomposition:")
    print(f"  Translation error: {metrics.error_components.translation_error:.3f} Å")
    print(f"  Orientation error: {metrics.error_components.orientation_error:.3f}°")
    print(f"  Rotation angle: {metrics.error_components.rotation_angle:.1f}°")
    
    if metrics.residue_error_stats:
        stats = metrics.residue_error_stats
        print(f"\nResidue RMSD Statistics:")
        print(f"  Mean: {stats['mean_rmsd']:.3f} Å")
        print(f"  Median: {stats['median_rmsd']:.3f} Å")
        print(f"  Std Dev: {stats['std_rmsd']:.3f} Å")
        print(f"  Range: {stats['min_rmsd']:.3f} - {stats['max_rmsd']:.3f} Å")
    
    print(f"\nWorst Residues (by RMSD):")
    for residue_id, rmsd in metrics.worst_residues:
        print(f"  {residue_id}: {rmsd:.3f} Å")
    
    print(f"\nBest Residues (by RMSD):")
    for residue_id, rmsd in metrics.best_residues:
        print(f"  {residue_id}: {rmsd:.3f} Å")


def export_metrics_csv(metrics: StructureMetrics, output_path: str) -> None:
    """Export comprehensive metrics to CSV file."""
    import csv
    
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header and basic metrics
        writer.writerow(['Metric', 'Value', 'Unit'])
        writer.writerow(['Overall RMSD', f"{metrics.overall_rmsd:.3f}", 'Å'])
        writer.writerow(['Protein RMSD', f"{metrics.protein_rmsd:.3f}", 'Å'])  
        writer.writerow(['DNA RMSD', f"{metrics.dna_rmsd:.3f}", 'Å'])
        writer.writerow(['Translation Error', f"{metrics.error_components.translation_error:.3f}", 'Å'])
        writer.writerow(['Rotation Angle', f"{metrics.error_components.rotation_angle:.1f}", '°'])
        
        # Write residue statistics
        if metrics.residue_error_stats:
            writer.writerow([])
            writer.writerow(['Residue Statistics', '', ''])
            for stat, value in metrics.residue_error_stats.items():
                writer.writerow([stat, f"{value:.3f}", 'Å'])
        
        # Write outliers
        writer.writerow([])
        writer.writerow(['Worst Residues', '', ''])
        for residue_id, rmsd in metrics.worst_residues:
            writer.writerow([residue_id, f"{rmsd:.3f}", 'Å'])
        
        writer.writerow([])
        writer.writerow(['Best Residues', '', ''])
        for residue_id, rmsd in metrics.best_residues:
            writer.writerow([residue_id, f"{rmsd:.3f}", 'Å'])


def generate_batch_report(all_results: list, failed_pairs: list, output_dir):
    """Generate summary report for batch processing"""
    try:
        import json
        from pathlib import Path
        
        output_dir = Path(output_dir)
        
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