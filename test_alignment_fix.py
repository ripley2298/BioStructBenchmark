#!/usr/bin/env python3
"""
Quick test to verify alignment fix is working
"""

import sys
from pathlib import Path

# Add the local package to the path
sys.path.insert(0, str(Path(__file__).parent))

from biostructbenchmark.core.alignment import align_structures_by_reference_frame, export_comprehensive_alignment_report
from biostructbenchmark.core.io import get_structure

def test_alignment_fix():
    """Test that alignment now produces realistic RMSD values"""
    
    # Load test structures
    exp_path = Path("tests/data/proteins_pdb/p456_02_experimental.pdb")
    pred_path = Path("tests/data/proteins_cif/p456_02_predicted.cif")
    
    if not exp_path.exists() or not pred_path.exists():
        print("Error: Test structure files not found")
        return False
    
    print("Loading structures...")
    exp_structure = get_structure(exp_path)
    pred_structure = get_structure(pred_path)
    
    print("Performing protein-only alignment...")
    result = align_structures_by_reference_frame(
        observed=exp_structure,
        predicted=pred_structure,
        reference_frame='protein',
        align_subset='protein'
    )
    
    print(f"\nAlignment Results:")
    print(f"Overall RMSD: {result.overall_rmsd:.3f} Å")
    print(f"Aligned atoms: {result.aligned_atom_count}")
    print(f"Reference frame: {result.reference_frame}")
    
    # Check if we have the transformed structures
    has_structures = hasattr(result, 'aligned_predicted_structure') and hasattr(result, 'reference_structure')
    print(f"Has transformed structures: {has_structures}")
    
    # Check per-chain RMSDs
    if result.chain_rmsds:
        print(f"\nPer-chain RMSDs:")
        for chain in result.chain_rmsds:
            print(f"  Chain {chain.chain_id}: {chain.rmsd:.3f} Å ({chain.atom_count} atoms)")
    
    # Test export
    output_dir = Path("alignment_outputs_test")
    output_dir.mkdir(exist_ok=True)
    
    print(f"\nExporting results to {output_dir}...")
    export_comprehensive_alignment_report(
        result=result,
        output_dir=output_dir,
        prefix="test_alignment"
    )
    
    # Validate RMSD is realistic
    if result.overall_rmsd < 5.0:  # Should be much better now
        print(f"\n✅ SUCCESS: RMSD is realistic ({result.overall_rmsd:.3f} Å < 5.0 Å)")
        return True
    else:
        print(f"\n❌ FAILURE: RMSD is still too high ({result.overall_rmsd:.3f} Å)")
        return False

if __name__ == "__main__":
    success = test_alignment_fix()
    sys.exit(0 if success else 1)