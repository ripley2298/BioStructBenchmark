#!/usr/bin/env python3
"""
Test script for DNA P-atom structural alignment
"""

import sys
from pathlib import Path

# Add the local package to the path
sys.path.insert(0, str(Path(__file__).parent))

from biostructbenchmark.core.alignment import (
    align_structures_with_structural_superimposition,
    align_dna_structures_with_phosphates,
    create_dna_correspondence_map,
    export_comprehensive_alignment_report
)
from biostructbenchmark.core.io import get_structure

def test_dna_p_alignment():
    """Test DNA P-atom structural alignment"""
    
    print("=" * 80)
    print("DNA P-ATOM STRUCTURAL ALIGNMENT TEST")
    print("=" * 80)
    
    # Load structures
    exp_path = Path("tests/data/proteins_pdb/p456_02_experimental.pdb")
    pred_path = Path("tests/data/proteins_cif/p456_02_predicted.cif")
    
    if not exp_path.exists() or not pred_path.exists():
        print("Error: Structure files not found")
        return
    
    print(f"Loading structures...")
    print(f"  Experimental: {exp_path}")
    print(f"  Predicted:    {pred_path}")
    
    exp_structure = get_structure(exp_path)
    pred_structure = get_structure(pred_path)
    
    print("✓ Structures loaded successfully")
    
    # Create output directory
    output_dir = Path("dna_p_alignment_test")
    output_dir.mkdir(exist_ok=True)
    
    # Test 1: DNA correspondence mapping
    print(f"\n1. TESTING DNA CORRESPONDENCE MAPPING...")
    correspondence_map = create_dna_correspondence_map(exp_structure, pred_structure)
    print(f"Created correspondence map with {len(correspondence_map)} nucleotide pairs:")
    
    for i, (obs_id, pred_id) in enumerate(correspondence_map.items()):
        if i < 10:  # Show first 10 for brevity
            print(f"  {obs_id} ↔ {pred_id}")
        elif i == 10:
            print(f"  ... and {len(correspondence_map) - 10} more pairs")
            break
    
    # Test 2: P-atom structural alignment
    print(f"\n2. TESTING P-ATOM STRUCTURAL ALIGNMENT...")
    superimposer, p_rmsd, p_count, aligned_structure = align_dna_structures_with_phosphates(
        exp_structure, pred_structure
    )
    
    if superimposer:
        print(f"✓ P-atom alignment successful:")
        print(f"  P-atoms used: {p_count}")
        print(f"  P-atom RMSD: {p_rmsd:.4f} Å")
    else:
        print("✗ P-atom alignment failed")
        return
    
    # Test 3: Three reference frame alignments
    print(f"\n3. TESTING THREE REFERENCE FRAME ALIGNMENTS...")
    
    reference_frames = ['global', 'dna_centric', 'protein_centric']
    results = {}
    
    for ref_frame in reference_frames:
        print(f"\n--- {ref_frame.upper()} ALIGNMENT ---")
        result = align_structures_with_structural_superimposition(
            exp_structure, pred_structure,
            reference_frame=ref_frame
        )
        results[ref_frame] = result
        
        # Export results for each reference frame
        export_comprehensive_alignment_report(
            result=result,
            output_dir=output_dir,
            prefix=f"{ref_frame}_alignment"
        )
    
    # Test 4: Comparison of all three reference frames
    print(f"\n4. COMPARING ALL THREE REFERENCE FRAMES...")
    print(f"{'Reference Frame':<15} {'Overall RMSD':<12} {'Protein RMSD':<13} {'DNA RMSD':<10}")
    print("-" * 55)
    
    for ref_frame, result in results.items():
        # Calculate separate protein and DNA RMSDs from chain data
        protein_rmsd = "N/A"
        dna_rmsd = "N/A"
        
        for chain in result.chain_rmsds:
            if chain.molecule_types == 'protein':
                protein_rmsd = f"{chain.rmsd:.4f}"
            elif 'dna' in chain.molecule_types:
                dna_rmsd = f"{chain.rmsd:.4f}"
        
        print(f"{ref_frame:<15} {result.overall_rmsd:<12.4f} {protein_rmsd:<13} {dna_rmsd:<10}")
    
    # Show detailed per-chain comparison
    print(f"\nDetailed per-chain RMSD comparison:")
    print(f"{'Chain':<8} {'Global':<8} {'DNA-centric':<12} {'Protein-centric':<16}")
    print("-" * 50)
    
    all_chains = set()
    for result in results.values():
        all_chains.update(chain.chain_id for chain in result.chain_rmsds)
    
    for chain_id in sorted(all_chains):
        rmsds = {}
        for ref_frame, result in results.items():
            chain_rmsd = next((chain.rmsd for chain in result.chain_rmsds if chain.chain_id == chain_id), 0.0)
            rmsds[ref_frame] = chain_rmsd
        
        print(f"{chain_id:<8} {rmsds.get('global', 0.0):<8.4f} {rmsds.get('dna_centric', 0.0):<12.4f} {rmsds.get('protein_centric', 0.0):<16.4f}")
    
    print(f"\n" + "=" * 80)
    print("DNA P-ATOM ALIGNMENT TEST COMPLETED!")
    print("=" * 80)
    print(f"\nGenerated outputs in: {output_dir.absolute()}")
    
    # List all created files
    output_files = sorted(output_dir.glob("*"))
    if output_files:
        print(f"\nFiles created:")
        for file_path in output_files:
            print(f"  {file_path.name}")

if __name__ == "__main__":
    test_dna_p_alignment()