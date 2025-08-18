#!/usr/bin/env python3
"""
Comprehensive alignment demo showing all three reference frame contexts:
1. Protein-aligned (protein reference frame)
2. DNA-aligned (DNA nucleotides reference frame) 
3. Center of mass average (full structure reference frame)
4. Kabsch algorithm with backbone-based reference frames

Each context shows RMSD values for all components.
"""

import sys
from pathlib import Path

# Add the local package to the path
sys.path.insert(0, str(Path(__file__).parent))

from biostructbenchmark.core.alignment import (
    align_structures_by_reference_frame, 
    align_structures_kabsch_backbone,
    export_comprehensive_alignment_report
)
from biostructbenchmark.core.io import get_structure

def print_alignment_summary(result, context_name):
    """Print detailed RMSD summary for an alignment result"""
    print(f"\n{context_name.upper()} ALIGNMENT RESULTS:")
    print("=" * 60)
    print(f"Reference Frame: {result.reference_frame}")
    print(f"Overall RMSD: {result.overall_rmsd:.3f} Å")
    print(f"Total Aligned Atoms: {result.aligned_atom_count}")
    
    print(f"\nPER-COMPONENT RMSD VALUES:")
    print("-" * 40)
    for chain_rmsd in result.chain_rmsds:
        chain_id = chain_rmsd.chain_id
        rmsd = chain_rmsd.rmsd
        atom_count = chain_rmsd.atom_count
        residue_count = chain_rmsd.residue_count
        molecule_type = chain_rmsd.molecule_types
        print(f"Chain {chain_id} ({molecule_type}): {rmsd:.3f} Å ({atom_count} atoms, {residue_count} residues)")

def main():
    """Main comprehensive alignment demo"""
    
    print("=" * 80)
    print("COMPREHENSIVE STRUCTURE ALIGNMENT - ALL REFERENCE FRAMES")
    print("=" * 80)
    
    # Define input structure paths
    exp_path = Path("tests/data/proteins_pdb/p456_02_experimental.pdb")
    pred_path = Path("tests/data/proteins_cif/p456_02_predicted.cif")
    
    if not exp_path.exists() or not pred_path.exists():
        print("Error: Structure files not found")
        return
    
    print("Loading structures...")
    print(f"  Experimental: {exp_path}")
    print(f"  Predicted:    {pred_path}")
    
    # Load structures
    exp_structure = get_structure(exp_path)
    pred_structure = get_structure(pred_path)
    
    print("✓ Structures loaded successfully")
    
    # Create output directory
    output_dir = Path("alignment_outputs_comprehensive")
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir.absolute()}")
    
    # =================================================================
    # 1. PROTEIN-ALIGNED CONTEXT
    # =================================================================
    print(f"\n1. PERFORMING PROTEIN-ALIGNED ANALYSIS...")
    
    # Protein-to-protein (protein reference, protein subset)
    protein_to_protein = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='protein',
        align_subset='protein'
    )
    
    # Protein-to-full (protein reference, all components)
    protein_to_full = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='protein',
        align_subset='full'
    )
    
    print_alignment_summary(protein_to_protein, "Protein-to-Protein")
    print_alignment_summary(protein_to_full, "Protein Reference with All Components")
    
    # Export protein-aligned results
    export_comprehensive_alignment_report(
        result=protein_to_protein,
        output_dir=output_dir,
        prefix="protein_aligned_protein_only"
    )
    
    export_comprehensive_alignment_report(
        result=protein_to_full,
        output_dir=output_dir,
        prefix="protein_aligned_full_structure"
    )
    
    # =================================================================
    # 2. DNA-ALIGNED CONTEXT  
    # =================================================================
    print(f"\n2. PERFORMING DNA-ALIGNED ANALYSIS...")
    
    # DNA-to-DNA (DNA reference, DNA subset)
    dna_to_dna = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='dna',
        align_subset='dna'
    )
    
    # DNA-to-full (DNA reference, all components)
    dna_to_full = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='dna',
        align_subset='full'
    )
    
    print_alignment_summary(dna_to_dna, "DNA-to-DNA")
    print_alignment_summary(dna_to_full, "DNA Reference with All Components")
    
    # Export DNA-aligned results
    export_comprehensive_alignment_report(
        result=dna_to_dna,
        output_dir=output_dir,
        prefix="dna_aligned_dna_only"
    )
    
    export_comprehensive_alignment_report(
        result=dna_to_full,
        output_dir=output_dir,
        prefix="dna_aligned_full_structure"
    )
    
    # =================================================================
    # 3. CENTER OF MASS AVERAGE CONTEXT
    # =================================================================
    print(f"\n3. PERFORMING CENTER OF MASS AVERAGE ANALYSIS...")
    
    # Full-to-full (all atoms reference, all atoms subset)
    full_to_full = align_structures_by_reference_frame(
        exp_structure, pred_structure,
        reference_frame='full',
        align_subset='full'
    )
    
    print_alignment_summary(full_to_full, "Center of Mass Average (Full-to-Full)")
    
    # Export center of mass results
    export_comprehensive_alignment_report(
        result=full_to_full,
        output_dir=output_dir,
        prefix="center_of_mass_full_structure"
    )
    
    # =================================================================
    # 4. KABSCH ALGORITHM WITH BACKBONE-BASED REFERENCE FRAMES
    # =================================================================
    print(f"\n4. PERFORMING KABSCH ALGORITHM ALIGNMENT...")
    
    # Kabsch protein backbone alignment
    kabsch_protein = align_structures_kabsch_backbone(
        exp_structure, pred_structure,
        reference_frame='protein',
        align_subset='full'
    )
    
    print_alignment_summary(kabsch_protein, "Kabsch Protein Backbone to Full Structure")
    
    # Export Kabsch protein results
    export_comprehensive_alignment_report(
        result=kabsch_protein,
        output_dir=output_dir,
        prefix="kabsch_protein_backbone_full_structure"
    )
    
    # Kabsch DNA backbone alignment
    kabsch_dna = align_structures_kabsch_backbone(
        exp_structure, pred_structure,
        reference_frame='dna',
        align_subset='full'
    )
    
    print_alignment_summary(kabsch_dna, "Kabsch DNA Backbone to Full Structure")
    
    # Export Kabsch DNA results
    export_comprehensive_alignment_report(
        result=kabsch_dna,
        output_dir=output_dir,
        prefix="kabsch_dna_backbone_full_structure"
    )
    
    # Kabsch full backbone alignment
    kabsch_full = align_structures_kabsch_backbone(
        exp_structure, pred_structure,
        reference_frame='full',
        align_subset='full'
    )
    
    print_alignment_summary(kabsch_full, "Kabsch Full Backbone to Full Structure")
    
    # Export Kabsch full results
    export_comprehensive_alignment_report(
        result=kabsch_full,
        output_dir=output_dir,
        prefix="kabsch_full_backbone_full_structure"
    )

    # =================================================================
    # SUMMARY COMPARISON
    # =================================================================
    print(f"\n" + "=" * 80)
    print("ALIGNMENT COMPARISON SUMMARY")
    print("=" * 80)
    
    alignments = [
        ("Protein-aligned (protein only)", protein_to_protein),
        ("Protein-aligned (all components)", protein_to_full),
        ("DNA-aligned (DNA only)", dna_to_dna),
        ("DNA-aligned (all components)", dna_to_full),
        ("Center of mass average", full_to_full),
        ("Kabsch protein backbone", kabsch_protein),
        ("Kabsch DNA backbone", kabsch_dna),
        ("Kabsch full backbone", kabsch_full)
    ]
    
    for name, result in alignments:
        print(f"\n{name}:")
        print(f"  Overall RMSD: {result.overall_rmsd:.3f} Å")
        for chain_rmsd in result.chain_rmsds:
            chain_id = chain_rmsd.chain_id
            rmsd = chain_rmsd.rmsd
            molecule_type = chain_rmsd.molecule_types
            print(f"    Chain {chain_id} ({molecule_type}): {rmsd:.3f} Å")
    
    print(f"\n" + "=" * 80)
    print("DEMO COMPLETED SUCCESSFULLY!")
    print("=" * 80)
    print(f"\nGenerated outputs in: {output_dir.absolute()}")
    
    # List all created files
    print(f"\nFiles created:")
    output_files = sorted(output_dir.glob("*"))
    for file_path in output_files:
        print(f"  {file_path.name}")
    
    print(f"\nSuperimposition files for visual validation:")
    superimposition_files = sorted(output_dir.glob("*superimposition.pdb"))
    for file_path in superimposition_files:
        print(f"  {file_path.name}")

if __name__ == "__main__":
    main()