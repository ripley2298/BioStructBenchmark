#!/usr/bin/env python3
"""
Debug script to investigate correspondence mapping algorithm failure
Focuses on p456_07 where ChimeraX can align but our algorithm fails
"""

from Bio.PDB import PDBParser, MMCIFParser
from biostructbenchmark.core.alignment import align_sequences, create_correspondence_map

def analyze_p456_07_alignment():
    """Debug the specific alignment failure for p456_07"""
    
    print("=== P456_07 ALIGNMENT DEBUGGING ===")
    
    # Load structures
    experimental_path = "/Users/mesler/Documents/IonuI_StructureML/experimental/p456_07_experimental.pdb"
    predicted_path = "/Users/mesler/Documents/IonuI_StructureML/predicted_alphafold3/p456_07_alphafold3.cif"
    
    parser_pdb = PDBParser(QUIET=True)
    parser_cif = MMCIFParser(QUIET=True)
    
    try:
        obs_structure = parser_pdb.get_structure("obs", experimental_path)
        pred_structure = parser_cif.get_structure("pred", predicted_path)
    except Exception as e:
        print(f"Structure loading failed: {e}")
        return
    
    print(f"✓ Loaded structures successfully")
    
    # Analyze protein chains (A should be the main protein chain)
    valid_aa = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
    
    for chain_id in ['A', 'B', 'C']:
        print(f"\n--- CHAIN {chain_id} ANALYSIS ---")
        
        # Get residues for this chain
        try:
            obs_chain = obs_structure[0][chain_id]
            pred_chain = pred_structure[0][chain_id]
        except KeyError:
            print(f"Chain {chain_id} not found in one or both structures")
            continue
        
        if chain_id == 'A':  # Focus on protein chain
            obs_residues = [r for r in obs_chain if r.get_resname().strip() in valid_aa]
            pred_residues = [r for r in pred_chain if r.get_resname().strip() in valid_aa]
        else:  # DNA chains
            dna_types = {'DA','DT','DG','DC'}
            obs_residues = [r for r in obs_chain if r.get_resname().strip() in dna_types]
            pred_residues = [r for r in pred_chain if r.get_resname().strip() in dna_types]
        
        print(f"Observed residues: {len(obs_residues)}")
        print(f"Predicted residues: {len(pred_residues)}")
        
        if not obs_residues or not pred_residues:
            print("No valid residues found")
            continue
        
        # Show first/last 5 residue positions and types
        print("\nFirst 5 observed residues:")
        for i, res in enumerate(obs_residues[:5]):
            print(f"  {i+1}: {res.get_id()[1]:4d} {res.get_resname().strip()}")
        
        print(f"... ({len(obs_residues)-10} more) ...")
        print("Last 5 observed residues:")
        for i, res in enumerate(obs_residues[-5:], len(obs_residues)-5):
            print(f"  {i+1}: {res.get_id()[1]:4d} {res.get_resname().strip()}")
        
        print("\nFirst 5 predicted residues:")
        for i, res in enumerate(pred_residues[:5]):
            print(f"  {i+1}: {res.get_id()[1]:4d} {res.get_resname().strip()}")
        
        print(f"... ({len(pred_residues)-10} more) ...")
        print("Last 5 predicted residues:")
        for i, res in enumerate(pred_residues[-5:], len(pred_residues)-5):
            print(f"  {i+1}: {res.get_id()[1]:4d} {res.get_resname().strip()}")
        
        # Test current alignment algorithm
        print("\n--- CURRENT ALIGNMENT ALGORITHM RESULT ---")
        aligned_pairs = align_sequences(obs_residues, pred_residues)
        print(f"Aligned pairs found: {len(aligned_pairs)}")
        
        if aligned_pairs:
            print("First 10 aligned pairs:")
            for i, (obs_res, pred_res) in enumerate(aligned_pairs[:10]):
                obs_pos = obs_res.get_id()[1]
                pred_pos = pred_res.get_id()[1]
                obs_name = obs_res.get_resname().strip()
                pred_name = pred_res.get_resname().strip()
                print(f"  {i+1:2d}: Obs {obs_pos:3d}:{obs_name} → Pred {pred_pos:3d}:{pred_name}")
            
            print(f"... ({len(aligned_pairs)-10} more pairs)")
            
            print("Last aligned pair:")
            obs_res, pred_res = aligned_pairs[-1]
            obs_pos = obs_res.get_id()[1]
            pred_pos = pred_res.get_id()[1]
            print(f"  Final: Obs {obs_pos:3d}:{obs_res.get_resname().strip()} → Pred {pred_pos:3d}:{pred_res.get_resname().strip()}")
        
        # Look at sequences around where alignment likely stops
        if chain_id == 'A' and len(aligned_pairs) > 0:
            print(f"\n--- ANALYZING WHY ALIGNMENT STOPS ---")
            
            # Get sequence strings
            obs_seq = [res.get_resname().strip() for res in obs_residues]
            pred_seq = [res.get_resname().strip() for res in pred_residues]
            
            # Find where the longest match ends
            last_aligned_obs_idx = obs_residues.index(aligned_pairs[-1][0])
            last_aligned_pred_idx = pred_residues.index(aligned_pairs[-1][1])
            
            print(f"Last aligned observed residue index: {last_aligned_obs_idx}")
            print(f"Last aligned predicted residue index: {last_aligned_pred_idx}")
            
            # Show context around where alignment stops
            start_context = max(0, last_aligned_obs_idx - 5)
            end_context = min(len(obs_residues), last_aligned_obs_idx + 10)
            
            print(f"\nObserved sequence context (around position {last_aligned_obs_idx}):")
            for i in range(start_context, end_context):
                marker = " → STOPS HERE" if i == last_aligned_obs_idx else ""
                pos = obs_residues[i].get_id()[1]
                seq = obs_seq[i]
                print(f"  Obs[{i:3d}]: Pos {pos:3d} = {seq}{marker}")
            
            pred_start_context = max(0, last_aligned_pred_idx - 5)
            pred_end_context = min(len(pred_residues), last_aligned_pred_idx + 10)
            
            print(f"\nPredicted sequence context (around position {last_aligned_pred_idx}):")
            for i in range(pred_start_context, pred_end_context):
                marker = " → STOPS HERE" if i == last_aligned_pred_idx else ""
                pos = pred_residues[i].get_id()[1]
                seq = pred_seq[i]
                print(f"  Pred[{i:3d}]: Pos {pos:3d} = {seq}{marker}")
            
            # Check what happens next
            if last_aligned_obs_idx + 1 < len(obs_residues) and last_aligned_pred_idx + 1 < len(pred_residues):
                next_obs = obs_seq[last_aligned_obs_idx + 1]
                next_pred = pred_seq[last_aligned_pred_idx + 1]
                print(f"\nNext residue after alignment stops:")
                print(f"  Observed: {next_obs}")
                print(f"  Predicted: {next_pred}")
                print(f"  Match: {'YES' if next_obs == next_pred else 'NO - THIS IS WHY ALIGNMENT STOPS'}")
    
if __name__ == "__main__":
    analyze_p456_07_alignment()