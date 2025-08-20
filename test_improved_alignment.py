#!/usr/bin/env python3
"""
Quick test to verify gap-tolerant alignment improvements with p456_07
"""

from Bio.PDB import PDBParser, MMCIFParser
from biostructbenchmark.core.alignment import align_structures_three_frames

def test_p456_07_improvement():
    """Test the improved alignment on p456_07"""
    
    print("=== TESTING IMPROVED GAP-TOLERANT ALIGNMENT ===")
    
    # Load structures
    experimental_path = "/Users/mesler/Documents/IonuI_StructureML/experimental/p456_07_experimental.pdb"
    predicted_path = "/Users/mesler/Documents/IonuI_StructureML/predicted_alphafold3/p456_07_alphafold3.cif"
    
    parser_pdb = PDBParser(QUIET=True)
    parser_cif = MMCIFParser(QUIET=True)
    
    obs_structure = parser_pdb.get_structure("obs", experimental_path)
    pred_structure = parser_cif.get_structure("pred", predicted_path)
    
    print(f"✓ Loaded p456_07 structures")
    
    # Run three-frame alignment with new algorithm
    print("\n--- RUNNING THREE-FRAME ALIGNMENT ---")
    results = align_structures_three_frames(obs_structure, pred_structure)
    
    print(f"\n=== ALIGNMENT RESULTS SUMMARY ===")
    for frame, result in results.items():
        print(f"\n{frame.upper()} ALIGNMENT:")
        print(f"  Overall RMSD: {result.overall_rmsd:.4f} Å")
        print(f"  Aligned atoms: {result.aligned_atom_count}")
        print(f"  Total residues: {len(result.residue_rmsds)}")
        
        # Count by molecule type
        protein_residues = sum(1 for r in result.residue_rmsds if r.is_protein)
        dna_residues = sum(1 for r in result.residue_rmsds if r.is_dna)
        print(f"  Protein residues: {protein_residues}")
        print(f"  DNA residues: {dna_residues}")
        
        # Show quality
        print(f"  Quality assessment: {result.reference_frame}")

if __name__ == "__main__":
    test_p456_07_improvement()