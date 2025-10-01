"""
biostructbenchmark/analysis/hbond.py
Hydrogen bond analysis for protein-DNA interactions using X3DNA-DSSR
"""

import subprocess
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from Bio.PDB import PDBParser


@dataclass
class HydrogenBond:
    """Container for hydrogen bond information"""
    donor_atom: str  # Format: chain:residue:atom
    acceptor_atom: str  # Format: chain:residue:atom
    donor_residue: str  # Format: chain:residue_name:position
    acceptor_residue: str  # Format: chain:residue_name:position
    donor_type: str  # 'protein' or 'nucleic'
    acceptor_type: str  # 'protein' or 'nucleic'
    distance: float  # Angstroms
    angle: Optional[float] = None  # Degrees
    interaction_type: str = "unknown"  # Type of interaction
    
    @property
    def bond_id(self) -> str:
        """Unique identifier for this hydrogen bond"""
        return f"{self.donor_atom}-->{self.acceptor_atom}"


@dataclass
class HBondComparison:
    """Container for hydrogen bond network comparison with interaction type sub-groupings"""
    experimental_bonds: List[HydrogenBond]
    predicted_bonds: List[HydrogenBond]
    
    # Category 1: Present in both structures (conserved interactions)
    common_bonds: List[Tuple[HydrogenBond, HydrogenBond]]  # (exp, pred) pairs
    common_protein_protein: List[Tuple[HydrogenBond, HydrogenBond]]
    common_protein_nucleic: List[Tuple[HydrogenBond, HydrogenBond]]
    common_nucleic_nucleic: List[Tuple[HydrogenBond, HydrogenBond]]
    
    # Category 2: Absent in predicted structure (missing interactions)
    experimental_only: List[HydrogenBond]
    missing_protein_protein: List[HydrogenBond]
    missing_protein_nucleic: List[HydrogenBond]
    missing_nucleic_nucleic: List[HydrogenBond]
    
    # Category 3: Additional in predicted structure (false positive interactions)
    predicted_only: List[HydrogenBond]
    false_positive_protein_protein: List[HydrogenBond]
    false_positive_protein_nucleic: List[HydrogenBond]
    false_positive_nucleic_nucleic: List[HydrogenBond]
    
    bond_distance_differences: Dict[str, float]  # bond_id -> distance difference


@dataclass
class HBondStatistics:
    """Summary statistics for hydrogen bond analysis"""
    total_experimental: int
    total_predicted: int
    total_common: int
    total_experimental_only: int
    total_predicted_only: int
    conservation_rate: float  # fraction of experimental bonds preserved
    prediction_accuracy: float  # fraction of predicted bonds that are correct
    mean_distance_difference: float  # for common bonds
    protein_to_protein_counts: Dict[str, int]  # experimental vs predicted counts
    protein_to_nucleic_counts: Dict[str, int]  # experimental vs predicted counts
    nucleic_to_nucleic_counts: Dict[str, int]  # experimental vs predicted counts


class HBondAnalyzer:
    """Analyze hydrogen bonds in protein-DNA complexes using X3DNA-DSSR"""
    
    def __init__(self, distance_tolerance: float = 0.5):
        """
        Initialize hydrogen bond analyzer
        
        Args:
            distance_tolerance: Distance tolerance for matching bonds (Angstroms)
        """
        self.distance_tolerance = distance_tolerance
        self.parser = PDBParser(QUIET=True)
    
    def extract_hbonds_with_dssr(self, structure_path: Path) -> List[HydrogenBond]:
        """
        Extract hydrogen bonds using X3DNA-DSSR
        
        Uses: x3dna-dssr -i=input_file --get-hbond --json
        
        Args:
            structure_path: Path to structure file
            
        Returns:
            List of HydrogenBond objects
        """
        hbonds = []
        
        try:
            # Run x3dna-dssr to get all hydrogen bonds
            cmd = f'x3dna-dssr -i="{structure_path}" --get-hbond --json'
            
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True,
                timeout=60
            )
            
            if result.returncode == 0 and result.stdout.strip():
                # DSSR output contains debug info followed by JSON
                # Find the JSON portion (starts with {"num_hbonds":...)
                stdout_lines = result.stdout.strip().split('\n')
                json_started = False
                
                for line in stdout_lines:
                    line = line.strip()
                    
                    # Look for the main JSON object containing hydrogen bonds
                    if line.startswith('{"num_hbonds":'):
                        try:
                            # Parse the main JSON containing all hydrogen bonds
                            dssr_data = json.loads(line)
                            hbond_list = dssr_data.get('hbonds', [])
                            
                            print(f"DSSR found {len(hbond_list)} total hydrogen bonds")
                            
                            for hbond_data in hbond_list:
                                # Filter for the interaction types we want
                                residue_pair = hbond_data.get('residue_pair', '')
                                if residue_pair in ['nt:aa', 'aa:aa', 'nt:nt']:
                                    hbond = self._parse_dssr_hbond(hbond_data)
                                    if hbond:
                                        hbonds.append(hbond)
                            break
                            
                        except json.JSONDecodeError as e:
                            print(f"JSON decode error: {e}")
                            continue
            else:
                print(f"X3DNA-DSSR failed with return code {result.returncode}")
                if result.stderr:
                    print(f"DSSR stderr: {result.stderr}")
                                
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Error running X3DNA-DSSR: {e}")
        
        print(f"Extracted {len(hbonds)} hydrogen bonds from {structure_path}")
        return hbonds
    
    def _parse_dssr_hbond(self, hbond_data: dict) -> Optional[HydrogenBond]:
        """Parse DSSR hydrogen bond data into HydrogenBond object"""
        try:
            # Extract information from DSSR format
            # Example: "atom1_id": "O@A.ALA16", "atom2_id": "OP2@C.DA16"
            atom1_id = hbond_data.get('atom1_id', '')
            atom2_id = hbond_data.get('atom2_id', '')
            distance = float(hbond_data.get('distance', 0.0))
            residue_pair = hbond_data.get('residue_pair', '')
            
            # Parse atom1 (format: "ATOM@CHAIN.RESIDUE")
            if '@' in atom1_id and '.' in atom1_id:
                atom1_name, chain_res1 = atom1_id.split('@')
                chain1, residue1 = chain_res1.split('.')
                # Extract residue name and number (e.g., "ALA16" -> "ALA", "16")
                residue1_name = ''.join(c for c in residue1 if c.isalpha())
                residue1_num = ''.join(c for c in residue1 if c.isdigit())
            else:
                return None
            
            # Parse atom2 (format: "ATOM@CHAIN.RESIDUE")
            if '@' in atom2_id and '.' in atom2_id:
                atom2_name, chain_res2 = atom2_id.split('@')
                chain2, residue2 = chain_res2.split('.')
                # Extract residue name and number
                residue2_name = ''.join(c for c in residue2 if c.isalpha())
                residue2_num = ''.join(c for c in residue2 if c.isdigit())
            else:
                return None
            
            # Determine molecule types
            nucleic_residues = {'DA', 'DT', 'DG', 'DC', 'A', 'T', 'G', 'C'}
            
            type1 = "nucleic" if residue1_name in nucleic_residues else "protein"
            type2 = "nucleic" if residue2_name in nucleic_residues else "protein"
            
            # For DSSR output, the order matters for determining donor/acceptor
            # Based on residue_pair, determine which is donor and which is acceptor
            if residue_pair == "nt:aa":
                # nucleotide to amino acid
                donor_atom = atom1_id
                acceptor_atom = atom2_id
                donor_residue = f"{chain1}:{residue1_name}:{residue1_num}"
                acceptor_residue = f"{chain2}:{residue2_name}:{residue2_num}"
                donor_type = type1
                acceptor_type = type2
            elif residue_pair == "aa:nt":
                # amino acid to nucleotide
                donor_atom = atom1_id
                acceptor_atom = atom2_id
                donor_residue = f"{chain1}:{residue1_name}:{residue1_num}"
                acceptor_residue = f"{chain2}:{residue2_name}:{residue2_num}"
                donor_type = type1
                acceptor_type = type2
            else:
                # For aa:aa and nt:nt, use atom1 as donor, atom2 as acceptor
                donor_atom = atom1_id
                acceptor_atom = atom2_id
                donor_residue = f"{chain1}:{residue1_name}:{residue1_num}"
                acceptor_residue = f"{chain2}:{residue2_name}:{residue2_num}"
                donor_type = type1
                acceptor_type = type2
            
            return HydrogenBond(
                donor_atom=donor_atom,
                acceptor_atom=acceptor_atom,
                donor_residue=donor_residue,
                acceptor_residue=acceptor_residue,
                donor_type=donor_type,
                acceptor_type=acceptor_type,
                distance=distance,
                angle=None,  # DSSR doesn't provide angle in your example
                interaction_type=residue_pair
            )
            
        except (KeyError, ValueError, TypeError) as e:
            print(f"Warning: Failed to parse DSSR hydrogen bond data: {e}")
            print(f"Data: {hbond_data}")
            return None
    
    def create_sequence_correspondence_from_structures(self, exp_structure_path: Path, 
                                                     pred_structure_path: Path) -> Dict:
        """
        Create sequence-based correspondence mapping directly from PDB structures
        
        This maps experimental PDB numbering to predicted PDB numbering using sequence alignment,
        which is essential for matching DSSR hydrogen bond output that uses original PDB numbering.
        
        Args:
            exp_structure_path: Path to experimental structure file
            pred_structure_path: Path to predicted structure file
            
        Returns:
            Dict mapping experimental residue identifiers to predicted identifiers
            Format: {exp_residue_key: pred_residue_key} where key = "chain:resname:resnum"
        """
        print(f"Creating sequence-based correspondence mapping from PDB structures...")
        
        try:
            # Load structures with appropriate parsers
            exp_structure = self.parser.get_structure("experimental", exp_structure_path)
            
            # Use specialized CIF parser for CIF files
            if str(pred_structure_path).endswith('.cif'):
                from Bio.PDB import MMCIFParser
                cif_parser = MMCIFParser(QUIET=True)
                pred_structure = cif_parser.get_structure("predicted", pred_structure_path)
            else:
                pred_structure = self.parser.get_structure("predicted", pred_structure_path)
            
            print(f"DEBUG: Experimental structure has models: {[m.get_id() for m in exp_structure]}")
            print(f"DEBUG: Predicted structure has models: {[m.get_id() for m in pred_structure]}")
            
            # Get first available models (handle different numbering between PDB and CIF)
            if len(list(exp_structure)) == 0:
                print("ERROR: No models in experimental structure")
                return {}
            if len(list(pred_structure)) == 0:
                print("ERROR: No models in predicted structure")
                return {}
                
            exp_model = list(exp_structure)[0]
            pred_model = list(pred_structure)[0]
            
            print(f"DEBUG: Loaded experimental structure with chains: {[c.get_id() for c in exp_model]}")
            print(f"DEBUG: Loaded predicted structure with chains: {[c.get_id() for c in pred_model]}")
            
            correspondence = {}
            
            # Process each chain in experimental structure
            for exp_chain in exp_model:
                exp_chain_id = exp_chain.get_id()
                exp_residues = list(exp_chain.get_residues())
                
                print(f"DEBUG: Experimental chain {exp_chain_id} has {len(exp_residues)} residues")
                
                if not exp_residues:
                    continue
                
                # Extract sequence from experimental chain
                exp_sequence = []
                exp_residue_map = {}  # position -> (residue, pdb_info)
                
                for i, res in enumerate(exp_residues):
                    resname = res.get_resname().strip()
                    exp_sequence.append(resname)
                    # Store original PDB info: (chain, resname, position, full_id)
                    exp_residue_map[i] = (res, exp_chain_id, resname, res.get_id()[1], res.get_id())
                
                # Find matching chain in predicted structure by sequence similarity
                best_match_chain = None
                best_similarity = 0
                
                for pred_chain in pred_model:
                    pred_chain_id = pred_chain.get_id()
                    pred_residues = list(pred_chain.get_residues())
                    if not pred_residues:
                        continue
                    
                    pred_sequence = [res.get_resname().strip() for res in pred_residues]
                    
                    # Calculate sequence similarity
                    similarity = self._calculate_sequence_similarity(exp_sequence, pred_sequence)
                    
                    print(f"DEBUG: Chain similarity {exp_chain_id} vs {pred_chain_id}: {similarity:.3f} "
                          f"(exp:{len(exp_sequence)}, pred:{len(pred_sequence)})")
                    
                    if similarity > best_similarity:  # Lowered threshold for debugging
                        best_similarity = similarity
                        best_match_chain = pred_chain
                
                if best_match_chain and best_similarity > 0.3:  # Lower threshold temporarily
                    pred_chain_id = best_match_chain.get_id()
                    pred_residues = list(best_match_chain.get_residues())
                    pred_sequence = [res.get_resname().strip() for res in pred_residues]
                    
                    # Perform sequence alignment
                    alignment = self._align_sequences_for_correspondence(exp_sequence, pred_sequence)
                    
                    print(f"DEBUG: Alignment produced {len(alignment)} correspondences")
                    
                    # Create correspondence mapping using original PDB identifiers
                    for exp_idx, pred_idx in alignment:
                        if exp_idx < len(exp_residue_map) and pred_idx < len(pred_residues):
                            exp_res, exp_cid, exp_rname, exp_rnum, exp_full_id = exp_residue_map[exp_idx]
                            pred_res = pred_residues[pred_idx]
                            pred_rname = pred_res.get_resname().strip()
                            pred_rnum = pred_res.get_id()[1]
                            
                            # Use format that matches DSSR output: "CHAIN.RESNAMERESNUM"
                            # Handle negative residue numbers (DNA often starts at negative positions)
                            exp_key = f"{exp_cid}.{exp_rname}{exp_rnum}"
                            pred_key = f"{pred_chain_id}.{pred_rname}{pred_rnum}"
                            
                            correspondence[exp_key] = pred_key
                            
                            if len(correspondence) <= 3:  # Debug first few
                                print(f"DEBUG: Mapped {exp_key} -> {pred_key}")
                    
                    print(f"Chain {exp_chain_id} -> {pred_chain_id}: {len(alignment)} residues aligned "
                          f"(similarity: {best_similarity:.2f})")
                else:
                    print(f"Warning: No suitable match found for experimental chain {exp_chain_id} "
                          f"(best similarity: {best_similarity:.3f})")
            
            print(f"Created sequence correspondence for {len(correspondence)} residues")
            return correspondence
            
        except Exception as e:
            print(f"Error creating sequence correspondence: {e}")
            import traceback
            traceback.print_exc()
            return {}
    
    def _calculate_sequence_similarity(self, seq1: List[str], seq2: List[str]) -> float:
        """Calculate sequence similarity between two residue lists using substring matching"""
        if not seq1 or not seq2:
            return 0.0

        # Method 1: Direct position comparison
        min_len = min(len(seq1), len(seq2))
        direct_matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        direct_similarity = direct_matches / max(len(seq1), len(seq2))

        # Method 2: Check for substring match (handles N/C-terminal extensions)
        substring_similarity = 0.0
        if len(seq1) <= len(seq2):
            # Check if seq1 is a substring of seq2
            for offset in range(len(seq2) - len(seq1) + 1):
                matches = sum(1 for i in range(len(seq1)) if seq1[i] == seq2[offset + i])
                substring_similarity = max(substring_similarity, matches / len(seq1))
        else:
            # Check if seq2 is a substring of seq1
            for offset in range(len(seq1) - len(seq2) + 1):
                matches = sum(1 for i in range(len(seq2)) if seq2[i] == seq1[offset + i])
                substring_similarity = max(substring_similarity, matches / len(seq2))

        # Use the better similarity score
        return max(direct_similarity, substring_similarity)
    
    def _align_sequences_for_correspondence(self, exp_seq: List[str], pred_seq: List[str]) -> List[Tuple[int, int]]:
        """
        Robust sequence alignment using BioPython PairwiseAligner

        Returns list of (exp_index, pred_index) pairs for aligned positions
        """
        from Bio.Align import PairwiseAligner

        # Amino acid and nucleotide 3-letter to 1-letter code mapping
        aa_codes = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            # DNA/RNA nucleotides
            'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T',
            'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U', 'T': 'T'
        }

        # Convert residue lists to single-letter codes for alignment
        exp_str = ''.join([aa_codes.get(r, 'X') for r in exp_seq])
        pred_str = ''.join([aa_codes.get(r, 'X') for r in pred_seq])

        # Create pairwise aligner
        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        # Perform alignment
        alignments = aligner.align(exp_str, pred_str)
        if not alignments:
            print("No alignment found, returning empty correspondence")
            return []

        best_alignment = alignments[0]

        # Map alignment back to residue pairs
        aligned_pairs = []

        # Get alignment coordinates to map back to residues
        for exp_block, pred_block in zip(best_alignment.aligned[0], best_alignment.aligned[1]):
            exp_start, exp_end = exp_block
            pred_start, pred_end = pred_block

            # Align residues in this block
            block_length = min(exp_end - exp_start, pred_end - pred_start)

            for i in range(block_length):
                exp_idx = exp_start + i
                pred_idx = pred_start + i
                if exp_idx < len(exp_seq) and pred_idx < len(pred_seq):
                    aligned_pairs.append((exp_idx, pred_idx))

        return aligned_pairs
    
    def compare_hydrogen_bonds_with_correspondence(self, experimental_hbonds: List[HydrogenBond],
                                                  predicted_hbonds: List[HydrogenBond],
                                                  correspondence_map: Dict) -> HBondComparison:
        """
        Compare hydrogen bond networks using sequence alignment-based correspondence mapping
        
        Uses "atom_pair" and "residue_pair" matching as specified by user:
        - Matching hydrogen bonds share same "atom_pair" (N:O, O:O, etc.)
        - Matching hydrogen bonds share same "residue_pair" (nt:aa, aa:aa, nt:nt)  
        - Sequence alignment resolves atom_id differences between experimental and predicted
        
        Args:
            experimental_hbonds: H-bonds from experimental structure
            predicted_hbonds: H-bonds from predicted structure  
            correspondence_map: Dict mapping experimental DSSR keys to predicted DSSR keys
                               Format: {"A.ALA16": "A.ALA16", "B.DG-2": "B.DG-2", ...}
            
        Returns:
            HBondComparison with proper sequence alignment-based matching and sub-groupings
        """
        print(f"Comparing hydrogen bonds using atom_pair/residue_pair matching with {len(correspondence_map)} residue correspondences")
        
        if correspondence_map:
            sample_items = list(correspondence_map.items())[:3]
            print(f"DEBUG: Sample correspondences: {sample_items}")
        
        # Initialize all categories and sub-groupings
        conserved_bonds = []
        conserved_protein_protein = []
        conserved_protein_nucleic = []
        conserved_nucleic_nucleic = []
        
        missing_in_predicted = []
        missing_protein_protein = []
        missing_protein_nucleic = []
        missing_nucleic_nucleic = []

        false_positive_in_predicted = list(predicted_hbonds)  # Start with all predicted bonds
        false_positive_protein_protein = []
        false_positive_protein_nucleic = []
        false_positive_nucleic_nucleic = []
        
        distance_differences = {}
        
        # Process each experimental hydrogen bond
        for i, exp_hb in enumerate(experimental_hbonds):
            if i < 3:  # Debug first 3 bonds
                print(f"DEBUG: Processing exp bond {i+1}: {exp_hb.donor_atom} -> {exp_hb.acceptor_atom}")
                print(f"       Residues: {exp_hb.donor_residue} -> {exp_hb.acceptor_residue}")
            
            # Extract DSSR-format keys from hydrogen bond atoms
            # Format: "N@A.ALA16" -> "A.ALA16"
            exp_donor_key = self._extract_dssr_residue_key(exp_hb.donor_atom)
            exp_acceptor_key = self._extract_dssr_residue_key(exp_hb.acceptor_atom)
            
            if exp_donor_key and exp_acceptor_key:
                # Look up corresponding predicted residues using sequence alignment
                pred_donor_key = correspondence_map.get(exp_donor_key)
                pred_acceptor_key = correspondence_map.get(exp_acceptor_key)
                
                if pred_donor_key and pred_acceptor_key:
                    # Find matching hydrogen bond in predicted structure
                    matched_bond = self._find_corresponding_hbond_by_dssr_keys(
                        exp_hb, predicted_hbonds, pred_donor_key, pred_acceptor_key)
                    
                    if matched_bond:
                        # Category 1: Present in both structures
                        bond_pair = (exp_hb, matched_bond)
                        conserved_bonds.append(bond_pair)
                        
                        # Sub-categorize by interaction type
                        interaction_type = self._get_interaction_type(exp_hb)
                        if interaction_type == "protein:protein":
                            conserved_protein_protein.append(bond_pair)
                        elif interaction_type == "protein:nucleic":
                            conserved_protein_nucleic.append(bond_pair)
                        elif interaction_type == "nucleic:nucleic":
                            conserved_nucleic_nucleic.append(bond_pair)
                        
                        distance_differences[exp_hb.bond_id] = matched_bond.distance - exp_hb.distance

                        # Remove the specific matched bond from false positive list
                        if matched_bond in false_positive_in_predicted:
                            false_positive_in_predicted.remove(matched_bond)

                        # Also remove any other predicted bonds between the same residue pair
                        # to prevent double-counting when structures have multiple H-bonds per residue pair
                        bonds_to_remove = []
                        for fp_bond in false_positive_in_predicted:
                            fp_donor_key = self._extract_dssr_residue_key(fp_bond.donor_atom)
                            fp_acceptor_key = self._extract_dssr_residue_key(fp_bond.acceptor_atom)
                            if (fp_donor_key == pred_donor_key and fp_acceptor_key == pred_acceptor_key):
                                bonds_to_remove.append(fp_bond)

                        for bond in bonds_to_remove:
                            false_positive_in_predicted.remove(bond)
                        
                        print(f"CONSERVED [{interaction_type}]: {exp_donor_key} -> {exp_acceptor_key} "
                              f"(Δd={matched_bond.distance - exp_hb.distance:.2f}Å)")
                    else:
                        # Category 2: Absent in predicted structure
                        self._categorize_missing_bond(exp_hb, missing_in_predicted, 
                                                    missing_protein_protein, missing_protein_nucleic, missing_nucleic_nucleic)
                else:
                    # No correspondence found (residues not aligned)
                    if i < 3:
                        print(f"DEBUG: No sequence correspondence for {exp_donor_key} or {exp_acceptor_key}")
                    self._categorize_missing_bond(exp_hb, missing_in_predicted,
                                                missing_protein_protein, missing_protein_nucleic, missing_nucleic_nucleic)
            else:
                # Malformed DSSR identifier
                if i < 3:
                    print(f"DEBUG: Could not extract DSSR keys from {exp_hb.donor_atom}, {exp_hb.acceptor_atom}")
                self._categorize_missing_bond(exp_hb, missing_in_predicted,
                                            missing_protein_protein, missing_protein_nucleic, missing_nucleic_nucleic)
        
        # Category 3: Sub-categorize false positive interactions
        for fp_hb in false_positive_in_predicted:
            interaction_type = self._get_interaction_type(fp_hb)
            if interaction_type == "protein:protein":
                false_positive_protein_protein.append(fp_hb)
            elif interaction_type == "protein:nucleic":
                false_positive_protein_nucleic.append(fp_hb)
            elif interaction_type == "nucleic:nucleic":
                false_positive_nucleic_nucleic.append(fp_hb)

            print(f"FALSE_POSITIVE [{interaction_type}]: {fp_hb.donor_residue} -> {fp_hb.acceptor_residue} "
                  f"(d={fp_hb.distance:.2f}Å)")
        
        # Print summary
        print(f"\n=== Hydrogen Bond Alignment Summary ===")
        print(f"Category 1 - Conserved: {len(conserved_bonds)} (PP:{len(conserved_protein_protein)}, "
              f"PN:{len(conserved_protein_nucleic)}, NN:{len(conserved_nucleic_nucleic)})")
        print(f"Category 2 - Missing: {len(missing_in_predicted)} (PP:{len(missing_protein_protein)}, "
              f"PN:{len(missing_protein_nucleic)}, NN:{len(missing_nucleic_nucleic)})")
        print(f"Category 3 - False Positives: {len(false_positive_in_predicted)} (PP:{len(false_positive_protein_protein)}, "
              f"PN:{len(false_positive_protein_nucleic)}, NN:{len(false_positive_nucleic_nucleic)})")
        
        return HBondComparison(
            experimental_bonds=experimental_hbonds,
            predicted_bonds=predicted_hbonds,
            common_bonds=conserved_bonds,
            common_protein_protein=conserved_protein_protein,
            common_protein_nucleic=conserved_protein_nucleic,
            common_nucleic_nucleic=conserved_nucleic_nucleic,
            experimental_only=missing_in_predicted,
            missing_protein_protein=missing_protein_protein,
            missing_protein_nucleic=missing_protein_nucleic,
            missing_nucleic_nucleic=missing_nucleic_nucleic,
            predicted_only=false_positive_in_predicted,
            false_positive_protein_protein=false_positive_protein_protein,
            false_positive_protein_nucleic=false_positive_protein_nucleic,
            false_positive_nucleic_nucleic=false_positive_nucleic_nucleic,
            bond_distance_differences=distance_differences
        )
    
    def _find_corresponding_hbond(self, exp_hb: HydrogenBond, predicted_hbonds: List[HydrogenBond],
                                 pred_donor_key: str, pred_acceptor_key: str) -> Optional[HydrogenBond]:
        """
        Find corresponding hydrogen bond using atom_pair and residue_pair matching
        
        Matching criteria based on user requirement:
        - Same "atom_pair" (e.g., N:O, O:O, etc.)  
        - Same "residue_pair" (e.g., nt:aa, aa:aa, nt:nt)
        - Use sequence alignment to resolve atom_id differences
        """
        # Extract atom pair from experimental bond
        exp_atom_pair = self._get_atom_pair(exp_hb)
        exp_residue_pair = exp_hb.interaction_type  # This is set from DSSR residue_pair
        
        best_match = None
        best_compatibility = 0
        
        for pred_hb in predicted_hbonds:
            pred_parts_donor = pred_hb.donor_residue.split(':')
            pred_parts_acceptor = pred_hb.acceptor_residue.split(':')
            
            if len(pred_parts_donor) >= 3 and len(pred_parts_acceptor) >= 3:
                pred_hb_donor_key = f"{pred_parts_donor[0]}:{pred_parts_donor[2]}"
                pred_hb_acceptor_key = f"{pred_parts_acceptor[0]}:{pred_parts_acceptor[2]}"
                
                # Check if this bond matches the sequence-aligned positions
                if (pred_hb_donor_key == pred_donor_key and pred_hb_acceptor_key == pred_acceptor_key):
                    
                    # Primary matching criteria: atom_pair and residue_pair
                    pred_atom_pair = self._get_atom_pair(pred_hb)
                    pred_residue_pair = pred_hb.interaction_type
                    
                    # Calculate compatibility score
                    compatibility = 0
                    
                    # Must have same residue_pair (nt:aa, aa:aa, nt:nt)
                    if pred_residue_pair == exp_residue_pair:
                        compatibility += 3  # High weight for residue pair match
                        
                        # Must have same atom_pair (N:O, O:O, etc.)
                        if pred_atom_pair == exp_atom_pair:
                            compatibility += 2  # High weight for atom pair match
                            
                            # Bonus for similar distance
                            if abs(pred_hb.distance - exp_hb.distance) <= self.distance_tolerance:
                                compatibility += 1
                        
                        # Store best match
                        if compatibility > best_compatibility:
                            best_compatibility = compatibility
                            best_match = pred_hb
        
        # Require minimum compatibility (residue_pair + atom_pair match)
        if best_compatibility >= 5:  # 3 + 2 = 5 minimum
            return best_match
        
        return None
    
    def _get_atom_pair(self, hbond: HydrogenBond) -> str:
        """
        Extract atom pair type from hydrogen bond (e.g., N:O, O:O)
        
        Parses DSSR atom IDs like "N@A.ALA16" to get atom type
        """
        try:
            # Extract atom names from DSSR format: "ATOM@CHAIN.RESIDUE"
            donor_atom_name = hbond.donor_atom.split('@')[0] if '@' in hbond.donor_atom else hbond.donor_atom
            acceptor_atom_name = hbond.acceptor_atom.split('@')[0] if '@' in hbond.acceptor_atom else hbond.acceptor_atom
            
            # Clean atom names (remove numbers, keep only element)
            donor_element = ''.join(c for c in donor_atom_name if c.isalpha())
            acceptor_element = ''.join(c for c in acceptor_atom_name if c.isalpha())
            
            return f"{donor_element}:{acceptor_element}"
            
        except (IndexError, AttributeError):
            return "unknown:unknown"
    
    def _extract_dssr_residue_key(self, atom_id: str) -> Optional[str]:
        """
        Extract DSSR residue key from atom ID
        
        Converts "N@A.ALA16" -> "A.ALA16"
        """
        try:
            if '@' in atom_id:
                return atom_id.split('@')[1]  # "A.ALA16"
            else:
                return None
        except (IndexError, AttributeError):
            return None
    
    def _find_corresponding_hbond_by_dssr_keys(self, exp_hb: HydrogenBond, predicted_hbonds: List[HydrogenBond],
                                              pred_donor_key: str, pred_acceptor_key: str) -> Optional[HydrogenBond]:
        """
        Find corresponding hydrogen bond using DSSR keys and atom_pair/residue_pair matching
        
        Args:
            exp_hb: Experimental hydrogen bond
            predicted_hbonds: List of predicted hydrogen bonds
            pred_donor_key: Predicted donor residue key (e.g., "A.ALA16")
            pred_acceptor_key: Predicted acceptor residue key (e.g., "B.DG-2")
        """
        exp_atom_pair = self._get_atom_pair(exp_hb)
        exp_residue_pair = exp_hb.interaction_type
        
        best_match = None
        best_compatibility = 0
        
        for pred_hb in predicted_hbonds:
            # Extract DSSR keys from predicted bond
            pred_hb_donor_key = self._extract_dssr_residue_key(pred_hb.donor_atom)
            pred_hb_acceptor_key = self._extract_dssr_residue_key(pred_hb.acceptor_atom)
            
            # Check if this bond involves the sequence-aligned residues
            if (pred_hb_donor_key == pred_donor_key and pred_hb_acceptor_key == pred_acceptor_key):
                
                # Primary matching criteria: atom_pair and residue_pair
                pred_atom_pair = self._get_atom_pair(pred_hb)
                pred_residue_pair = pred_hb.interaction_type
                
                # Calculate compatibility score
                compatibility = 0
                
                # Must have same residue_pair (nt:aa, aa:aa, nt:nt)
                if pred_residue_pair == exp_residue_pair:
                    compatibility += 3  # High weight for residue pair match
                    
                    # Must have same atom_pair (N:O, O:O, etc.)
                    if pred_atom_pair == exp_atom_pair:
                        compatibility += 2  # High weight for atom pair match
                        
                        # Bonus for similar distance
                        if abs(pred_hb.distance - exp_hb.distance) <= self.distance_tolerance:
                            compatibility += 1
                
                # Store best match
                if compatibility > best_compatibility:
                    best_compatibility = compatibility
                    best_match = pred_hb
        
        # Require minimum compatibility (residue_pair + atom_pair match)
        if best_compatibility >= 5:  # 3 + 2 = 5 minimum
            return best_match
        
        return None
    
    def _get_interaction_type(self, hbond: HydrogenBond) -> str:
        """Determine interaction type based on donor and acceptor molecule types"""
        donor_type = "nucleic" if hbond.donor_type == "dna" else hbond.donor_type
        acceptor_type = "nucleic" if hbond.acceptor_type == "dna" else hbond.acceptor_type
        
        if donor_type == "protein" and acceptor_type == "protein":
            return "protein:protein"
        elif donor_type == "nucleic" and acceptor_type == "nucleic":
            return "nucleic:nucleic"
        elif (donor_type == "protein" and acceptor_type == "nucleic") or \
             (donor_type == "nucleic" and acceptor_type == "protein"):
            return "protein:nucleic"
        else:
            return f"{donor_type}:{acceptor_type}"
    
    def _categorize_missing_bond(self, exp_hb: HydrogenBond, missing_list: List,
                               missing_pp: List, missing_pn: List, missing_nn: List):
        """Helper to categorize missing bonds by interaction type"""
        missing_list.append(exp_hb)
        interaction_type = self._get_interaction_type(exp_hb)
        
        if interaction_type == "protein:protein":
            missing_pp.append(exp_hb)
        elif interaction_type == "protein:nucleic":
            missing_pn.append(exp_hb)
        elif interaction_type == "nucleic:nucleic":
            missing_nn.append(exp_hb)
        
        print(f"MISSING [{interaction_type}]: {exp_hb.donor_residue} -> {exp_hb.acceptor_residue} "
              f"(d={exp_hb.distance:.2f}Å)")
    
    def calculate_statistics(self, comparison: HBondComparison) -> HBondStatistics:
        """Calculate summary statistics from hydrogen bond comparison"""
        total_exp = len(comparison.experimental_bonds)
        total_pred = len(comparison.predicted_bonds) 
        total_common = len(comparison.common_bonds)
        total_exp_only = len(comparison.experimental_only)
        total_pred_only = len(comparison.predicted_only)
        
        conservation_rate = total_common / total_exp if total_exp > 0 else 0.0
        prediction_accuracy = total_common / total_pred if total_pred > 0 else 0.0
        
        distance_diffs = list(comparison.bond_distance_differences.values())
        mean_distance_diff = np.mean(distance_diffs) if distance_diffs else 0.0
        
        # Count by interaction types
        pp_exp = len(comparison.common_protein_protein) + len(comparison.missing_protein_protein)
        pp_pred = len(comparison.common_protein_protein) + len(comparison.false_positive_protein_protein)

        pn_exp = len(comparison.common_protein_nucleic) + len(comparison.missing_protein_nucleic)
        pn_pred = len(comparison.common_protein_nucleic) + len(comparison.false_positive_protein_nucleic)

        nn_exp = len(comparison.common_nucleic_nucleic) + len(comparison.missing_nucleic_nucleic)
        nn_pred = len(comparison.common_nucleic_nucleic) + len(comparison.false_positive_nucleic_nucleic)
        
        return HBondStatistics(
            total_experimental=total_exp,
            total_predicted=total_pred,
            total_common=total_common,
            total_experimental_only=total_exp_only,
            total_predicted_only=total_pred_only,
            conservation_rate=conservation_rate,
            prediction_accuracy=prediction_accuracy,
            mean_distance_difference=mean_distance_diff,
            protein_to_protein_counts={'experimental': pp_exp, 'predicted': pp_pred},
            protein_to_nucleic_counts={'experimental': pn_exp, 'predicted': pn_pred},
            nucleic_to_nucleic_counts={'experimental': nn_exp, 'predicted': nn_pred}
        )
    
    def analyze_structures_with_correspondence(self, experimental_path: Path, predicted_path: Path, 
                                             correspondence_map: Optional[Dict] = None) -> Tuple[HBondComparison, HBondStatistics]:
        """
        Analyze hydrogen bond networks using X3DNA-DSSR with sequence-based correspondence mapping
        
        Args:
            experimental_path: Path to experimental structure
            predicted_path: Path to predicted structure  
            correspondence_map: Optional pre-computed correspondence map. If None, will create from sequences.
                               
        Returns:
            (HBondComparison, HBondStatistics) with X3DNA-DSSR analysis
        """
        try:
            # Use X3DNA-DSSR for hydrogen bond extraction
            exp_hbonds = self.extract_hbonds_with_dssr(experimental_path)
            pred_hbonds = self.extract_hbonds_with_dssr(predicted_path)
            
            print(f"Extracted {len(exp_hbonds)} experimental and {len(pred_hbonds)} predicted hydrogen bonds")
            
            # Create sequence-based correspondence if not provided
            if correspondence_map is None:
                print("Creating sequence-based correspondence mapping from PDB structures...")
                correspondence_map = self.create_sequence_correspondence_from_structures(
                    experimental_path, predicted_path)
            
            # Use correspondence-aware comparison with atom_pair/residue_pair matching
            comparison = self.compare_hydrogen_bonds_with_correspondence(
                exp_hbonds, pred_hbonds, correspondence_map)
            statistics = self.calculate_statistics(comparison)
            
            return comparison, statistics
            
        except Exception as e:
            print(f"Error: X3DNA-DSSR hydrogen bond analysis failed: {e}")
            # Return empty results rather than crashing
            empty_comparison = HBondComparison(
                experimental_bonds=[], predicted_bonds=[], common_bonds=[],
                common_protein_protein=[], common_protein_nucleic=[], common_nucleic_nucleic=[],
                experimental_only=[], missing_protein_protein=[], missing_protein_nucleic=[], missing_nucleic_nucleic=[],
                predicted_only=[], false_positive_protein_protein=[], false_positive_protein_nucleic=[], false_positive_nucleic_nucleic=[],
                bond_distance_differences={}
            )
            empty_stats = HBondStatistics(
                total_experimental=0, total_predicted=0, total_common=0, total_experimental_only=0, total_predicted_only=0,
                conservation_rate=0.0, prediction_accuracy=0.0, mean_distance_difference=0.0,
                protein_to_protein_counts={'experimental': 0, 'predicted': 0},
                protein_to_nucleic_counts={'experimental': 0, 'predicted': 0},
                nucleic_to_nucleic_counts={'experimental': 0, 'predicted': 0}
            )
            
            return empty_comparison, empty_stats
    
    def export_results(self, comparison: HBondComparison, statistics: HBondStatistics,
                      output_dir: Path, pair_id: str):
        """
        Export hydrogen bond analysis results
        
        Args:
            comparison: HBondComparison object
            statistics: HBondStatistics object
            output_dir: Output directory
            pair_id: Structure pair identifier
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        # Export detailed hydrogen bond lists
        self._export_hbond_details(comparison, output_dir / f"{pair_id}_hbond_details.csv")

        # Export comparison summary
        self._export_hbond_summary(comparison, statistics, output_dir / f"{pair_id}_hbond_summary.csv")

        # Export detailed text report
        self._export_hbond_report(comparison, statistics, output_dir / f"{pair_id}_detailed_report.txt", pair_id)
    
    def _export_hbond_details(self, comparison: HBondComparison, output_path: Path):
        """Export detailed hydrogen bond information"""
        data = []
        
        # Common bonds
        for exp_hb, pred_hb in comparison.common_bonds:
            data.append({
                'bond_type': 'conserved',
                'interaction_type': self._get_interaction_type(exp_hb),
                'donor_residue_exp': exp_hb.donor_residue,
                'acceptor_residue_exp': exp_hb.acceptor_residue,
                'distance_exp': exp_hb.distance,
                'donor_residue_pred': pred_hb.donor_residue,
                'acceptor_residue_pred': pred_hb.acceptor_residue,
                'distance_pred': pred_hb.distance,
                'distance_difference': round(pred_hb.distance - exp_hb.distance, 3)
            })
        
        # Missing bonds
        for hb in comparison.experimental_only:
            data.append({
                'bond_type': 'missing',
                'interaction_type': self._get_interaction_type(hb),
                'donor_residue_exp': hb.donor_residue,
                'acceptor_residue_exp': hb.acceptor_residue,
                'distance_exp': hb.distance,
                'donor_residue_pred': None,
                'acceptor_residue_pred': None,
                'distance_pred': None,
                'distance_difference': None
            })
        
        # False positive bonds
        for hb in comparison.predicted_only:
            data.append({
                'bond_type': 'false_positive',
                'interaction_type': self._get_interaction_type(hb),
                'donor_residue_exp': None,
                'acceptor_residue_exp': None,
                'distance_exp': None,
                'donor_residue_pred': hb.donor_residue,
                'acceptor_residue_pred': hb.acceptor_residue,
                'distance_pred': hb.distance,
                'distance_difference': None
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
    
    def _export_hbond_summary(self, comparison: HBondComparison, statistics: HBondStatistics, 
                             output_path: Path):
        """Export hydrogen bond comparison summary"""
        summary_data = [
            {'metric': 'total_experimental_bonds', 'value': statistics.total_experimental},
            {'metric': 'total_predicted_bonds', 'value': statistics.total_predicted},
            {'metric': 'conserved_bonds_total', 'value': statistics.total_common},
            {'metric': 'conserved_protein_protein', 'value': len(comparison.common_protein_protein)},
            {'metric': 'conserved_protein_nucleic', 'value': len(comparison.common_protein_nucleic)},
            {'metric': 'conserved_nucleic_nucleic', 'value': len(comparison.common_nucleic_nucleic)},
            {'metric': 'missing_bonds_total', 'value': statistics.total_experimental_only},
            {'metric': 'missing_protein_protein', 'value': len(comparison.missing_protein_protein)},
            {'metric': 'missing_protein_nucleic', 'value': len(comparison.missing_protein_nucleic)},
            {'metric': 'missing_nucleic_nucleic', 'value': len(comparison.missing_nucleic_nucleic)},
            {'metric': 'false_positive_bonds_total', 'value': statistics.total_predicted_only},
            {'metric': 'false_positive_protein_protein', 'value': len(comparison.false_positive_protein_protein)},
            {'metric': 'false_positive_protein_nucleic', 'value': len(comparison.false_positive_protein_nucleic)},
            {'metric': 'false_positive_nucleic_nucleic', 'value': len(comparison.false_positive_nucleic_nucleic)},
            {'metric': 'conservation_rate', 'value': statistics.conservation_rate},
            {'metric': 'prediction_accuracy', 'value': statistics.prediction_accuracy},
            {'metric': 'mean_distance_difference', 'value': statistics.mean_distance_difference}
        ]
        
        df = pd.DataFrame(summary_data)
        df.to_csv(output_path, index=False)

    def _export_hbond_report(self, comparison: HBondComparison, statistics: HBondStatistics,
                            output_path: Path, pair_id: str):
        """Export detailed text report of hydrogen bond analysis"""
        with open(output_path, 'w') as f:
            # Header
            f.write("HYDROGEN BOND ANALYSIS REPORT\n")
            f.write("EXPERIMENTAL vs PREDICTED STRUCTURE COMPARISON\n")
            f.write("=" * 70 + "\n\n")

            # Structure information
            f.write("STRUCTURE INFORMATION\n")
            f.write("-" * 25 + "\n")
            parts = pair_id.split("_vs_")
            if len(parts) == 2:
                f.write(f"Experimental: {parts[0]}\n")
                f.write(f"Predicted: {parts[1]}\n\n")
            else:
                f.write(f"Structure Pair: {pair_id}\n\n")

            # Analysis summary
            f.write("ANALYSIS SUMMARY\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total experimental hydrogen bonds: {statistics.total_experimental}\n")
            f.write(f"Total predicted hydrogen bonds: {statistics.total_predicted}\n")
            f.write(f"Conserved bonds: {statistics.total_common}\n")
            f.write(f"Missing bonds: {statistics.total_experimental_only}\n")
            f.write(f"False positive bonds: {statistics.total_predicted_only}\n")
            f.write(f"Conservation rate: {statistics.conservation_rate:.1%}\n")
            f.write(f"Prediction accuracy: {statistics.prediction_accuracy:.1%}\n")
            f.write(f"Mean distance difference: {statistics.mean_distance_difference:.3f} Å\n\n")

            # Detailed interaction analysis
            f.write("DETAILED INTERACTION ANALYSIS\n")
            f.write("-" * 35 + "\n")
            f.write("PROTEIN-PROTEIN INTERACTIONS:\n")
            f.write(f"  Conserved: {len(comparison.common_protein_protein)}\n")
            f.write(f"  Missing:   {len(comparison.missing_protein_protein)}\n")
            f.write(f"  False Pos: {len(comparison.false_positive_protein_protein)}\n\n")

            f.write("PROTEIN-NUCLEIC INTERACTIONS:\n")
            f.write(f"  Conserved: {len(comparison.common_protein_nucleic)}\n")
            f.write(f"  Missing:   {len(comparison.missing_protein_nucleic)}\n")
            f.write(f"  False Pos: {len(comparison.false_positive_protein_nucleic)}\n\n")

            f.write("NUCLEIC-NUCLEIC INTERACTIONS:\n")
            f.write(f"  Conserved: {len(comparison.common_nucleic_nucleic)}\n")
            f.write(f"  Missing:   {len(comparison.missing_nucleic_nucleic)}\n")
            f.write(f"  False Pos: {len(comparison.false_positive_nucleic_nucleic)}\n\n")

            # Conserved bonds (first 20)
            f.write("CONSERVED HYDROGEN BONDS (First 20)\n")
            f.write("-" * 40 + "\n")
            for i, (exp_hb, pred_hb) in enumerate(comparison.common_bonds[:20], 1):
                delta = pred_hb.distance - exp_hb.distance
                sign = "+" if delta >= 0 else ""
                f.write(f"{i:2d}. {exp_hb.donor_residue} -> {exp_hb.acceptor_residue}\n")
                f.write(f"    Exp: {exp_hb.distance:.2f} Å, Pred: {pred_hb.distance:.2f} Å (Δ={sign}{delta:.2f} Å)\n")
                f.write(f"    Type: {self._get_interaction_type(exp_hb)}\n")
                f.write(f"    Atoms: {exp_hb.donor_atom} -> {exp_hb.acceptor_atom}\n\n")

            # Missing bonds (first 10)
            f.write("MISSING HYDROGEN BONDS (First 10)\n")
            f.write("-" * 35 + "\n")
            for i, hb in enumerate(comparison.experimental_only[:10], 1):
                f.write(f"{i:2d}. {hb.donor_residue} -> {hb.acceptor_residue}\n")
                f.write(f"    Distance: {hb.distance:.2f} Å\n")
                f.write(f"    Type: {self._get_interaction_type(hb)}\n")
                f.write(f"    Atoms: {hb.donor_atom} -> {hb.acceptor_atom}\n\n")

            # False positive bonds (first 10)
            f.write("FALSE POSITIVE HYDROGEN BONDS (First 10)\n")
            f.write("-" * 32 + "\n")
            for i, hb in enumerate(comparison.predicted_only[:10], 1):
                f.write(f"{i:2d}. {hb.donor_residue} -> {hb.acceptor_residue}\n")
                f.write(f"    Distance: {hb.distance:.2f} Å\n")
                f.write(f"    Type: {self._get_interaction_type(hb)}\n")
                f.write(f"    Atoms: {hb.donor_atom} -> {hb.acceptor_atom}\n\n")