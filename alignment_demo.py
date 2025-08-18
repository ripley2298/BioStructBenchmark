#!/usr/bin/env python3
"""
Comprehensive Alignment Demo Script

This script demonstrates how to use the enhanced alignment functionality to generate:
1. PDB files of aligned structures
2. Per-residue/nucleotide RMSD data
3. Per-chain RMSD values  
4. Comprehensive alignment reports

Usage:
    python alignment_demo.py

This will align the p456_02 experimental and predicted structures and generate
all output files in the 'alignment_outputs' directory.
"""

from pathlib import Path
from biostructbenchmark.core.alignment import (
    align_structures_by_reference_frame,
    export_comprehensive_alignment_report
)
from biostructbenchmark.core.io import get_structure


def main():
    """Main demo function"""
    
    # Define input structure paths
    exp_path = Path("tests/data/proteins_pdb/p456_02_experimental.pdb")
    pred_path = Path("tests/data/proteins_cif/p456_02_predicted.cif")
    
    # Check if files exist
    if not exp_path.exists():
        print(f"Error: Experimental structure not found: {exp_path}")
        return
    if not pred_path.exists():
        print(f"Error: Predicted structure not found: {pred_path}")
        return
    
    print("=" * 60)
    print("COMPREHENSIVE STRUCTURE ALIGNMENT DEMO")
    print("=" * 60)
    
    # Load structures
    print(f"\nLoading structures...")
    print(f"  Experimental: {exp_path}")
    print(f"  Predicted:    {pred_path}")
    
    exp_structure = get_structure(exp_path)
    pred_structure = get_structure(pred_path)
    
    if exp_structure is None or pred_structure is None:
        print("Error: Failed to load structures")
        return
    
    print("✓ Structures loaded successfully")
    
    # Create output directory
    output_dir = Path("alignment_outputs")
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir.absolute()}")
    
    # Perform protein-only alignment
    print(f"\nPerforming protein-only alignment...")
    result = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='protein',
        align_subset='protein'
    )
    
    print(f"✓ Alignment completed")
    print(f"  Overall RMSD: {result.overall_rmsd:.3f} Å")
    print(f"  Aligned atoms: {result.aligned_atom_count}")
    print(f"  Chains found: {len(result.chain_rmsds)}")
    print(f"  Residues analyzed: {len(result.residue_rmsds)}")
    
    # Display per-chain RMSD summary
    print(f"\nPer-chain RMSD summary:")
    for chain in result.chain_rmsds:
        mol_types = ', '.join(chain.molecule_types)
        print(f"  Chain {chain.chain_id}: {chain.rmsd:.3f} Å ({mol_types})")
    
    # Export comprehensive report
    print(f"\nExporting comprehensive alignment report...")
    export_comprehensive_alignment_report(
        result=result,
        output_dir=output_dir,
        prefix="p456_protein_alignment"
    )
    
    # Also perform full structure alignment for comparison
    print(f"\nPerforming full structure alignment (including DNA/ligands)...")
    full_result = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='full',
        align_subset='full'
    )
    
    print(f"✓ Full structure alignment completed")
    print(f"  Overall RMSD: {full_result.overall_rmsd:.3f} Å")
    print(f"  Aligned atoms: {full_result.aligned_atom_count}")
    
    # Export full structure report
    export_comprehensive_alignment_report(
        result=full_result,
        output_dir=output_dir,
        prefix="p456_full_alignment"
    )
    
    print(f"\n" + "=" * 60)
    print("DEMO COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"\nGenerated outputs in: {output_dir.absolute()}")
    print(f"\nFiles created:")
    
    # List all created files
    output_files = sorted(output_dir.glob("*"))
    for file_path in output_files:
        print(f"  {file_path.name}")
    
    print(f"\nTo validate the alignment:")
    print(f"1. Open the aligned PDB files in ChimeraX, PyMOL, or similar")
    print(f"2. Check the CSV files for detailed RMSD data")
    print(f"3. Review the summary text files for comprehensive statistics")
    
    # Show some key statistics
    protein_residues = [r for r in result.residue_rmsds if r.molecule_type == 'protein']
    if protein_residues:
        protein_rmsds = [r.rmsd for r in protein_residues]
        under_2A = sum(1 for rmsd in protein_rmsds if rmsd < 2.0)
        under_5A = sum(1 for rmsd in protein_rmsds if rmsd < 5.0)
        
        print(f"\nProtein alignment quality:")
        print(f"  Mean RMSD: {sum(protein_rmsds)/len(protein_rmsds):.3f} Å")
        print(f"  Residues < 2.0 Å: {under_2A}/{len(protein_rmsds)} ({under_2A/len(protein_rmsds)*100:.1f}%)")
        print(f"  Residues < 5.0 Å: {under_5A}/{len(protein_rmsds)} ({under_5A/len(protein_rmsds)*100:.1f}%)")


if __name__ == "__main__":
    main()