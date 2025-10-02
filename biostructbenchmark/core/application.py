"""
Core application logic for BioStructBenchmark
Contains the main analysis orchestration function
"""

from pathlib import Path
from typing import Dict, Any


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
        from biostructbenchmark.core.runners import run_multi_frame_alignment
        results['multi_frame'] = run_multi_frame_alignment(
            exp_path, pred_path, output_dir, args
        )
    
    # Single-frame RMSD (if not doing multi-frame)
    elif flags.get('rmsd') and not flags.get('multi_frame'):
        from biostructbenchmark.core.runners import run_basic_rmsd
        results['rmsd'] = run_basic_rmsd(
            exp_path, pred_path, output_dir, args
        )
    
    # B-factor analysis
    if flags.get('bfactor'):
        from biostructbenchmark.analysis.bfactor import run_bfactor_analysis
        results['bfactor'] = run_bfactor_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Consensus analysis
    if flags.get('consensus'):
        from biostructbenchmark.analysis.consensus import run_consensus_analysis
        results['consensus'] = run_consensus_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Mutation analysis
    if flags.get('mutations'):
        from biostructbenchmark.analysis.mutations import run_mutation_analysis
        results['mutations'] = run_mutation_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Hydrogen bond analysis
    if flags.get('hbond'):
        from biostructbenchmark.core.runners import run_hydrogen_bond_analysis
        results['hbond'] = run_hydrogen_bond_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # DSSR DNA structural analysis
    if flags.get('dssr'):
        from biostructbenchmark.core.runners import run_dssr_analysis
        results['dssr'] = run_dssr_analysis(
            exp_path, pred_path, output_dir, args
        )
    
    # Visualization
    if flags.get('visualize'):
        from biostructbenchmark.visualization.coordinator import generate_visualizations
        results['visualization'] = generate_visualizations(
            results, output_dir, args
        )
    
    return results