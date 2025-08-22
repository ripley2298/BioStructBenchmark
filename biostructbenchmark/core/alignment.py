#!/usr/bin/env python3
"""
Streamlined structure alignment module for protein-DNA complexes
Supports sequence-based correspondence and three reference frames
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from Bio.PDB import Structure as BioStructure, Superimposer
from Bio.Align import PairwiseAligner
from dataclasses import dataclass


# --- DATA STRUCTURES ---

@dataclass
class ResidueRMSD:
    """Container for structural unit RMSD data (residues for protein, bases for DNA)"""
    residue_id: str  # For compatibility: unit_id (residue for protein, base for DNA)
    residue_type: str  # For compatibility: unit_type (amino acid or nucleotide type)
    chain_id: str
    position: int
    rmsd: float
    atom_count: int
    molecule_type: str  # 'protein' or 'dna'
    dna_atoms: Optional[List['Bio.PDB.Atom.Atom']] = None  # DNA atoms for distance calculations
    
    @property
    def unit_id(self) -> str:
        """Proper nomenclature: unit_id (residue for protein, base for DNA)"""
        return self.residue_id
    
    @property
    def unit_type(self) -> str:
        """Proper nomenclature: unit_type (amino acid for protein, nucleotide for DNA)"""
        return self.residue_type
    
    @property
    def is_protein(self) -> bool:
        """True if this represents a protein residue"""
        return self.molecule_type == 'protein'
    
    @property
    def is_dna(self) -> bool:
        """True if this represents a DNA base/nucleotide"""
        return self.molecule_type == 'dna'
    
    def min_phosphate_distance(self, protein_atoms: Optional[List['Bio.PDB.Atom.Atom']] = None) -> float:
        """
        Returns min distance (Å) from protein residue to DNA phosphate atoms (P, O1P, O2P).
        
        Args:
            protein_atoms: List of protein atoms for distance calculation.
                          If None, assumes this is a protein residue and uses self.dna_atoms for DNA
        
        Returns:
            Minimum distance to phosphate atoms, or float('inf') if no DNA atoms available
        """
        if not self.dna_atoms:
            return float('inf')
        
        # Determine source and target atoms based on molecule type
        if self.molecule_type == 'protein':
            # This is a protein residue, calculate distance to DNA phosphates
            source_atoms = protein_atoms if protein_atoms else []
            target_atoms = [atom for atom in self.dna_atoms if any(p in atom.name for p in ['P', 'O1P', 'O2P'])]
        else:
            # This is a DNA residue, no phosphate distance calculation needed
            return float('inf')
        
        if not source_atoms or not target_atoms:
            return float('inf')
        
        # Calculate minimum distance using Biopython's atom distance calculation
        min_dist = float('inf')
        for protein_atom in source_atoms:
            for dna_atom in target_atoms:
                try:
                    # Use Biopython's distance calculation (coord vector subtraction + norm)
                    distance = np.linalg.norm(protein_atom.coord - dna_atom.coord)
                    min_dist = min(min_dist, distance)
                except (AttributeError, TypeError):
                    # Handle cases where coordinates are not available
                    continue
        
        return min_dist


@dataclass  
class ChainRMSD:
    chain_id: str
    rmsd: float
    atom_count: int
    residue_count: int
    molecule_types: str


@dataclass
class AlignmentResult:
    overall_rmsd: float
    residue_rmsds: List[ResidueRMSD]
    chain_rmsds: List[ChainRMSD]
    transformation_matrix: np.ndarray
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray
    aligned_atom_count: int
    reference_frame: str
    aligned_predicted_structure: BioStructure
    reference_structure: BioStructure


# --- SEQUENCE ALIGNMENT ---

def align_sequences(obs_residues: List, pred_residues: List) -> List[Tuple]:
    """
    Gap-tolerant sequence alignment for experimental structures with missing residues.
    Uses sequence alignment to properly handle gaps and missing residues.
    """
    if not obs_residues or not pred_residues:
        return []
    
    # Create position-to-residue maps for efficient lookup
    obs_pos_map = {res.get_id()[1]: res for res in obs_residues}
    pred_pos_map = {res.get_id()[1]: res for res in pred_residues}
    
    # Get position ranges
    obs_positions = sorted(obs_pos_map.keys())
    pred_positions = sorted(pred_pos_map.keys())
    
    print(f"Gap-tolerant alignment: Obs positions {obs_positions[0]}-{obs_positions[-1]} ({len(obs_positions)} residues)")
    print(f"                       Pred positions {pred_positions[0]}-{pred_positions[-1]} ({len(pred_positions)} residues)")
    
    # Detect gaps in experimental structure
    gaps = []
    for i in range(len(obs_positions)-1):
        if obs_positions[i+1] - obs_positions[i] > 1:
            gap_start = obs_positions[i] + 1
            gap_end = obs_positions[i+1] - 1
            gaps.append((gap_start, gap_end))
    
    if gaps:
        print(f"Detected gaps in experimental structure: {gaps}")
    
    # First try direct position matching for cases without complex gaps
    direct_matches = []
    
    for pos in obs_positions:
        if pos in pred_pos_map:
            obs_res = obs_pos_map[pos]
            pred_res = pred_pos_map[pos]
            # Verify sequence identity for the match
            if obs_res.get_resname().strip() == pred_res.get_resname().strip():
                direct_matches.append((obs_res, pred_res))
    
    print(f"Direct position matches: {len(direct_matches)}")
    
    # If direct matching worked well, use it
    if len(direct_matches) >= min(len(obs_residues), len(pred_residues)) * 0.8:
        print(f"Using direct position matching: {len(direct_matches)} pairs")
        return direct_matches
    
    # Strategy 2: Use sequence alignment to handle gaps properly
    print("Using sequence alignment to handle gaps...")
    aligned_pairs = align_sequences_with_gaps(obs_residues, pred_residues)
    
    print(f"Sequence alignment resulted in: {len(aligned_pairs)} aligned pairs")
    return aligned_pairs


def get_single_letter_code(resname: str) -> str:
    """
    Convert 3-letter amino acid code to single letter, handling DNA/RNA nucleotides.
    """
    code_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # DNA nucleotides
        'DA': 'A', 'DT': 'T', 'DG': 'G', 'DC': 'C',
        # RNA nucleotides
        'A': 'A', 'U': 'U', 'G': 'G', 'C': 'C'
    }
    return code_map.get(resname.strip(), 'X')


def align_sequences_with_gaps(obs_residues: List, pred_residues: List) -> List[Tuple]:
    """
    Use BioPython sequence alignment to properly handle gaps in experimental structures.
    """
    # Build sequence strings
    obs_seq_str = ""
    pred_seq_str = ""
    
    for res in obs_residues:
        obs_seq_str += get_single_letter_code(res.get_resname())
    
    for res in pred_residues:
        pred_seq_str += get_single_letter_code(res.get_resname())
    
    # Create pairwise aligner
    aligner = PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(obs_seq_str, pred_seq_str)
    if not alignments:
        print("No alignment found, falling back to direct position matching")
        return []
    
    best_alignment = alignments[0]  # Get best alignment
    
    print(f"Alignment score: {best_alignment.score}")
    
    # Map alignment back to residue pairs
    aligned_pairs = []
    obs_idx = 0
    pred_idx = 0
    
    # Get alignment coordinates to map back to residues
    for obs_block, pred_block in zip(best_alignment.aligned[0], best_alignment.aligned[1]):
        obs_start, obs_end = obs_block
        pred_start, pred_end = pred_block
        
        # Align residues in this block
        block_length = min(obs_end - obs_start, pred_end - pred_start)
        
        for i in range(block_length):
            if obs_start + i < len(obs_residues) and pred_start + i < len(pred_residues):
                obs_res = obs_residues[obs_start + i]
                pred_res = pred_residues[pred_start + i]
                aligned_pairs.append((obs_res, pred_res))
    
    return aligned_pairs


def find_alignment_segments(obs_seq: List[str], pred_seq: List[str], 
                          obs_residues: List, pred_residues: List) -> List[Dict]:
    """
    Find multiple alignment segments to handle gaps in experimental structures.
    Returns segments with their aligned pairs.
    """
    segments = []
    min_segment_length = 5  # Minimum consecutive matches to form a segment
    
    # Track which positions have been used to avoid double-counting
    used_obs = set()
    used_pred = set()
    
    # Find all possible alignment segments
    for obs_start in range(len(obs_seq)):
        if obs_start in used_obs:
            continue
            
        for pred_start in range(len(pred_seq)):
            if pred_start in used_pred:
                continue
                
            # Try to extend alignment from this starting point
            segment_pairs = []
            obs_pos, pred_pos = obs_start, pred_start
            consecutive_matches = 0
            
            while (obs_pos < len(obs_seq) and pred_pos < len(pred_seq) and 
                   obs_pos not in used_obs and pred_pos not in used_pred):
                
                if obs_seq[obs_pos] == pred_seq[pred_pos]:
                    # Match found
                    segment_pairs.append((obs_residues[obs_pos], pred_residues[pred_pos]))
                    consecutive_matches += 1
                    obs_pos += 1
                    pred_pos += 1
                else:
                    # Mismatch - try to handle gaps
                    if consecutive_matches >= min_segment_length:
                        # Save current segment if it's long enough
                        break
                    else:
                        # Try skipping gaps in either sequence
                        gap_handled = False
                        
                        # Try skipping in observed sequence (experimental gap)
                        if (obs_pos + 1 < len(obs_seq) and 
                            obs_seq[obs_pos + 1] == pred_seq[pred_pos]):
                            obs_pos += 1  # Skip gap in observed
                            gap_handled = True
                        
                        # Try skipping in predicted sequence  
                        elif (pred_pos + 1 < len(pred_seq) and 
                              obs_seq[obs_pos] == pred_seq[pred_pos + 1]):
                            pred_pos += 1  # Skip gap in predicted
                            gap_handled = True
                        
                        # Try skipping in both (rare but possible)
                        elif (obs_pos + 1 < len(obs_seq) and pred_pos + 1 < len(pred_seq) and
                              obs_seq[obs_pos + 1] == pred_seq[pred_pos + 1]):
                            obs_pos += 1
                            pred_pos += 1
                            gap_handled = True
                        
                        if not gap_handled:
                            # Can't handle this gap, break segment
                            break
            
            # Add segment if it's substantial
            if len(segment_pairs) >= min_segment_length:
                # Mark positions as used
                for i, (obs_res, pred_res) in enumerate(segment_pairs):
                    used_obs.add(obs_start + i)
                    used_pred.add(pred_start + i)
                
                segments.append({
                    'obs_start': obs_start,
                    'pred_start': pred_start,
                    'length': len(segment_pairs),
                    'pairs': segment_pairs
                })
                
                print(f"Found alignment segment: obs[{obs_start}:{obs_start + len(segment_pairs)}] → "
                      f"pred[{pred_start}:{pred_start + len(segment_pairs)}] ({len(segment_pairs)} residues)")
    
    # Sort segments by length (longest first) to prioritize quality alignments
    segments.sort(key=lambda x: x['length'], reverse=True)
    
    return segments


def create_correspondence_map(observed: BioStructure, predicted: BioStructure, 
                            mol_type: str) -> Dict:
    """Create sequence-based correspondence map for protein or DNA"""
    correspondence = {}
    
    # Define residue types
    if mol_type == 'protein':
        valid_types = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                      'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
    else:  # DNA
        valid_types = {'DA','DT','DG','DC'}
    
    # Extract residues by chain
    for chain in observed[0]:
        chain_id = chain.get_id()
        obs_residues = [r for r in chain if r.get_resname().strip() in valid_types]
        
        if chain_id in [c.get_id() for c in predicted[0]] and obs_residues:
            pred_chain = predicted[0][chain_id]
            pred_residues = [r for r in pred_chain if r.get_resname().strip() in valid_types]
            
            if pred_residues:
                aligned_pairs = align_sequences(obs_residues, pred_residues)
                for obs_res, pred_res in aligned_pairs:
                    obs_key = (chain_id, obs_res.get_id()[1])
                    pred_key = (chain_id, pred_res.get_id()[1])
                    correspondence[obs_key] = pred_key
    
    return correspondence


# --- STRUCTURE ALIGNMENT ---

def get_backbone_atoms(structure: BioStructure, correspondence: Dict, 
                      atom_name: str, is_observed: bool = True) -> List:
    """Extract backbone atoms (CA or P) using correspondence mapping"""
    atoms = []
    for (obs_chain, obs_pos), (pred_chain, pred_pos) in correspondence.items():
        try:
            if is_observed:
                # Use observed structure positions
                residue = structure[0][obs_chain][obs_pos]
            else:
                # Use predicted structure positions  
                residue = structure[0][pred_chain][pred_pos]
            
            if atom_name in residue:
                atoms.append(residue[atom_name])
        except KeyError:
            continue
    return atoms


def get_matched_backbone_atoms(obs_structure: BioStructure, pred_structure: BioStructure,
                              correspondence: Dict, atom_name: str) -> Tuple[List, List]:
    """
    Extract matched backbone atoms ensuring equal counts for superimposition.
    Only returns atoms where both observed and predicted residues have the required atom.
    """
    obs_atoms = []
    pred_atoms = []
    
    for (obs_chain, obs_pos), (pred_chain, pred_pos) in correspondence.items():
        try:
            obs_residue = obs_structure[0][obs_chain][obs_pos]
            pred_residue = pred_structure[0][pred_chain][pred_pos]
            
            # Only add if both residues have the required atom
            if atom_name in obs_residue and atom_name in pred_residue:
                obs_atoms.append(obs_residue[atom_name])
                pred_atoms.append(pred_residue[atom_name])
        except KeyError:
            continue
    
    return obs_atoms, pred_atoms


def align_structures_three_frames(observed: BioStructure, predicted: BioStructure) -> Dict:
    """Perform alignment in three reference frames: global, DNA-centric, protein-centric"""
    import copy
    
    results = {}
    
    # Get correspondence maps
    protein_corr = create_correspondence_map(observed, predicted, 'protein')
    dna_corr = create_correspondence_map(observed, predicted, 'dna')
    
    print(f"Correspondence: {len(protein_corr)} protein, {len(dna_corr)} DNA")
    
    for ref_frame in ['global', 'dna_centric', 'protein_centric']:
        print(f"\n--- {ref_frame.upper()} ALIGNMENT ---")
        
        # Create copy for transformation
        pred_copy = copy.deepcopy(predicted)
        
        # Get alignment atoms based on reference frame using matched atom extraction
        if ref_frame == 'global':
            # Use both CA and P atoms with matched extraction
            protein_obs, protein_pred = get_matched_backbone_atoms(observed, pred_copy, protein_corr, 'CA')
            dna_obs, dna_pred = get_matched_backbone_atoms(observed, pred_copy, dna_corr, 'P')
            obs_atoms = protein_obs + dna_obs
            pred_atoms = protein_pred + dna_pred
            
        elif ref_frame == 'dna_centric':
            # Use only P atoms with matched extraction
            obs_atoms, pred_atoms = get_matched_backbone_atoms(observed, pred_copy, dna_corr, 'P')
            
        else:  # protein_centric
            # Use only CA atoms with matched extraction
            obs_atoms, pred_atoms = get_matched_backbone_atoms(observed, pred_copy, protein_corr, 'CA')
        
        # Perform superimposition
        if len(obs_atoms) >= 3:
            superimposer = Superimposer()
            superimposer.set_atoms(obs_atoms, pred_atoms)
            superimposer.apply(pred_copy.get_atoms())
            
            backbone_rmsd = superimposer.rms
            print(f"Backbone RMSD: {backbone_rmsd:.4f} Å ({len(obs_atoms)} atoms)")
            
            # Calculate per-residue and per-chain RMSDs
            residue_rmsds = calculate_per_residue_rmsds(observed, pred_copy, 
                                                      protein_corr, dna_corr)
            chain_rmsds = calculate_per_chain_rmsds_from_residues(residue_rmsds)
            
            # Calculate overall RMSD using weighted average from residue RMSDs (correspondence-based)
            total_sq_diff = 0.0
            total_atoms = 0
            for rmsd in residue_rmsds:
                total_sq_diff += rmsd.atom_count * (rmsd.rmsd ** 2)
                total_atoms += rmsd.atom_count
            overall_rmsd = np.sqrt(total_sq_diff / total_atoms) if total_atoms > 0 else 0.0
            
            print(f"Overall RMSD: {overall_rmsd:.4f} Å")
            
            # Store result
            results[ref_frame] = AlignmentResult(
                overall_rmsd=overall_rmsd,
                residue_rmsds=residue_rmsds,
                chain_rmsds=chain_rmsds,
                transformation_matrix=np.eye(4),
                rotation_matrix=np.eye(3), 
                translation_vector=np.zeros(3),
                aligned_atom_count=len(obs_atoms),
                reference_frame=f"{ref_frame}_alignment",
                aligned_predicted_structure=pred_copy,
                reference_structure=observed
            )
        else:
            print(f"Insufficient atoms ({len(obs_atoms)}) for {ref_frame} alignment")
    
    return results


# --- RMSD CALCULATIONS ---

def calculate_residue_rmsd(obs_res, pred_res) -> Optional[ResidueRMSD]:
    """Calculate RMSD between two residues/nucleotides"""
    obs_atoms = {atom.get_name(): atom for atom in obs_res.get_atoms()}
    pred_atoms = {atom.get_name(): atom for atom in pred_res.get_atoms()}
    
    # Find common atoms
    common_names = list(obs_atoms.keys() & pred_atoms.keys())
    if len(common_names) < 2:
        return None
    
    # Calculate RMSD
    obs_coords = np.array([obs_atoms[name].get_coord() for name in common_names])
    pred_coords = np.array([pred_atoms[name].get_coord() for name in common_names])
    diff = obs_coords - pred_coords
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    # Create result
    chain_id = obs_res.get_parent().get_id()
    is_dna = obs_res.get_resname().strip() in {'DA','DT','DG','DC'}
    
    return ResidueRMSD(
        residue_id=f"{chain_id}:{obs_res.get_id()[1]}",
        residue_type=obs_res.get_resname(),
        chain_id=chain_id,
        position=obs_res.get_id()[1],
        rmsd=rmsd,
        atom_count=len(common_names),
        molecule_type='dna' if is_dna else 'protein'
    )


def calculate_per_residue_rmsds(observed: BioStructure, predicted: BioStructure,
                               protein_corr: Dict, dna_corr: Dict) -> List[ResidueRMSD]:
    """Calculate per-residue RMSDs using sequence correspondence"""
    residue_rmsds = []
    
    # Process protein residues
    for (obs_chain, obs_pos), (pred_chain, pred_pos) in protein_corr.items():
        try:
            obs_res = observed[0][obs_chain][obs_pos]
            pred_res = predicted[0][pred_chain][pred_pos]
            rmsd_result = calculate_residue_rmsd(obs_res, pred_res)
            if rmsd_result:
                residue_rmsds.append(rmsd_result)
        except KeyError:
            continue
    
    # Process DNA nucleotides
    for (obs_chain, obs_pos), (pred_chain, pred_pos) in dna_corr.items():
        try:
            obs_res = observed[0][obs_chain][obs_pos] 
            pred_res = predicted[0][pred_chain][pred_pos]
            rmsd_result = calculate_residue_rmsd(obs_res, pred_res)
            if rmsd_result:
                residue_rmsds.append(rmsd_result)
        except KeyError:
            continue
    
    return residue_rmsds


def calculate_per_chain_rmsds_from_residues(residue_rmsds: List[ResidueRMSD]) -> List[ChainRMSD]:
    """Calculate per-chain RMSDs by aggregating residue-level data (post-alignment)."""
    chain_data = {}
    
    for res_rmsd in residue_rmsds:
        chain_id = res_rmsd.chain_id
        if chain_id not in chain_data:
            chain_data[chain_id] = {
                'weighted_sq_sum': 0.0,
                'total_atoms': 0,
                'residue_count': 0,
                'molecule_type': res_rmsd.molecule_type
            }
        
        # Weighted sum for proper averaging
        chain_data[chain_id]['weighted_sq_sum'] += res_rmsd.atom_count * (res_rmsd.rmsd ** 2)
        chain_data[chain_id]['total_atoms'] += res_rmsd.atom_count
        chain_data[chain_id]['residue_count'] += 1
    
    # Calculate final chain RMSDs
    chain_rmsds = []
    for chain_id, data in chain_data.items():
        if data['total_atoms'] > 0:
            chain_rmsd = np.sqrt(data['weighted_sq_sum'] / data['total_atoms'])
            
            chain_rmsds.append(ChainRMSD(
                chain_id=chain_id,
                rmsd=chain_rmsd,
                atom_count=data['total_atoms'],
                residue_count=data['residue_count'],
                molecule_types=data['molecule_type']
            ))
    
    return chain_rmsds


def calculate_per_chain_rmsds(observed: BioStructure, predicted: BioStructure) -> List[ChainRMSD]:
    """Calculate per-chain RMSDs using raw coordinates (DEPRECATED - may give inflated values)."""
    import warnings
    warnings.warn("calculate_per_chain_rmsds uses raw coordinates and may give inflated RMSD values. "
                  "Use calculate_per_chain_rmsds_from_residues instead.", DeprecationWarning, stacklevel=2)
    
    chain_rmsds = []
    
    for obs_chain in observed[0]:
        chain_id = obs_chain.get_id()
        if chain_id in [c.get_id() for c in predicted[0]]:
            pred_chain = predicted[0][chain_id]
            
            # Get all atoms from both chains
            obs_atoms = list(obs_chain.get_atoms())
            pred_atoms = list(pred_chain.get_atoms())
            
            if obs_atoms and pred_atoms:
                min_len = min(len(obs_atoms), len(pred_atoms))
                obs_coords = np.array([atom.get_coord() for atom in obs_atoms[:min_len]])
                pred_coords = np.array([atom.get_coord() for atom in pred_atoms[:min_len]])
                
                diff = obs_coords - pred_coords
                rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
                
                # Determine molecule type
                residues = list(obs_chain.get_residues())
                if residues:
                    first_res = residues[0].get_resname().strip()
                    mol_type = 'dna' if first_res in {'DA','DT','DG','DC'} else 'protein'
                else:
                    mol_type = 'unknown'
                
                chain_rmsds.append(ChainRMSD(
                    chain_id=chain_id,
                    rmsd=rmsd,
                    atom_count=min_len,
                    residue_count=len(residues),
                    molecule_types=mol_type
                ))
    
    return chain_rmsds


# --- MAIN ALIGNMENT FUNCTION ---

def align_structures_with_structural_superimposition(observed: BioStructure,
                                                   predicted: BioStructure,
                                                   reference_frame: str = 'global') -> AlignmentResult:
    """
    Main alignment function supporting three reference frames:
    - global: Combined protein + DNA backbone alignment
    - dna_centric: DNA P-atom alignment for DNA quality assessment  
    - protein_centric: Protein CA-atom alignment for interface assessment
    """
    results = align_structures_three_frames(observed, predicted)
    
    if reference_frame in results:
        return results[reference_frame]
    else:
        raise ValueError(f"Reference frame '{reference_frame}' not found")


# Legacy Functions


def calculate_per_residue_rmsd_for_subset(observed: BioStructure, predicted: BioStructure,
                                        molecule_type: str = 'full') -> List[ResidueRMSD]:
    """Legacy function for backward compatibility"""
    protein_corr = create_correspondence_map(observed, predicted, 'protein')
    dna_corr = create_correspondence_map(observed, predicted, 'dna')
    return calculate_per_residue_rmsds(observed, predicted, protein_corr, dna_corr)


# --- LEGACY COMPATIBILITY FUNCTIONS ---

def align_dna_structures_with_phosphates(observed: BioStructure, predicted: BioStructure) -> Tuple:
    """Legacy function - returns DNA correspondence and P-atom RMSD"""
    dna_corr = create_correspondence_map(observed, predicted, 'dna')
    
    # Calculate P-atom RMSD without transformation
    p_coords_obs = []
    p_coords_pred = []
    
    for (obs_chain, obs_pos), (pred_chain, pred_pos) in dna_corr.items():
        try:
            obs_nuc = observed[0][obs_chain][obs_pos]
            pred_nuc = predicted[0][pred_chain][pred_pos] 
            if 'P' in obs_nuc and 'P' in pred_nuc:
                p_coords_obs.append(obs_nuc['P'].get_coord())
                p_coords_pred.append(pred_nuc['P'].get_coord())
        except KeyError:
            continue
    
    if p_coords_obs:
        p_coords_obs = np.array(p_coords_obs)
        p_coords_pred = np.array(p_coords_pred)
        diff = p_coords_obs - p_coords_pred
        p_rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        print(f"Found {len(p_coords_obs)} P atom pairs for DNA nucleotide correspondence")
        print(f"  P-atom RMSD (original coordinates): {p_rmsd:.4f} Å")
        return dna_corr, p_rmsd, len(p_coords_obs), predicted
    else:
        return dna_corr, None, 0, predicted


def create_dna_correspondence_map(observed: BioStructure, predicted: BioStructure) -> Dict:
    """Legacy function - returns DNA correspondence map"""
    return create_correspondence_map(observed, predicted, 'dna')


def export_comprehensive_alignment_report(result: AlignmentResult, output_dir, prefix: str):
    """Legacy function - simplified export"""
    from pathlib import Path
    import csv
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Export per-unit RMSD CSV with proper nomenclature
    csv_path = output_dir / f"{prefix}_per_unit_rmsd.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['unit_id', 'unit_type', 'chain_id', 'position', 
                        'rmsd', 'atom_count', 'molecule_type', 'reference_frame', 'unit_class'])
        for rmsd in result.residue_rmsds:
            unit_class = 'amino_acid' if rmsd.is_protein else 'nucleotide'
            writer.writerow([rmsd.unit_id, rmsd.unit_type, rmsd.chain_id, 
                           rmsd.position, rmsd.rmsd, rmsd.atom_count, 
                           rmsd.molecule_type, result.reference_frame, unit_class])
    
    # Export per-chain RMSD CSV
    chain_csv_path = output_dir / f"{prefix}_per_chain_rmsd.csv"
    with open(chain_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['chain_id', 'rmsd', 'atom_count', 'residue_count', 
                        'molecule_types', 'reference_frame'])
        for chain in result.chain_rmsds:
            writer.writerow([chain.chain_id, chain.rmsd, chain.atom_count,
                           chain.residue_count, chain.molecule_types, result.reference_frame])
    
    # Export summary
    summary_path = output_dir / f"{prefix}_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("============================================================\n")
        f.write("STRUCTURE ALIGNMENT COMPREHENSIVE REPORT\n") 
        f.write("============================================================\n\n")
        f.write(f"Reference Frame: {result.reference_frame}\n")
        f.write(f"Overall RMSD: {result.overall_rmsd:.3f} Å\n")
        f.write(f"Total Aligned Atoms: {result.aligned_atom_count}\n\n")
        f.write("PER-CHAIN RMSD VALUES:\n")
        f.write("----------------------------------------\n")
        for chain in result.chain_rmsds:
            f.write(f"Chain {chain.chain_id}: {chain.rmsd:.3f} Å "
                   f"({chain.atom_count} atoms, {chain.residue_count} residues, "
                   f"{chain.molecule_types})\n")
    
    print(f"Exported per-unit RMSD data (residues+bases) to: {csv_path}")
    print(f"Exported per-chain RMSD data to: {chain_csv_path}")
    print(f"Summary report: {summary_path}")


# --- API COMPATIBILITY FUNCTIONS ---

def compare_structures(observed: BioStructure, predicted: BioStructure, 
                      method: str = 'structural') -> AlignmentResult:
    """Main API function for structure comparison"""
    if method == 'structural':
        return align_structures_with_structural_superimposition(observed, predicted, 'global')
    else:
        # Fallback to global alignment
        return align_structures_with_structural_superimposition(observed, predicted, 'global')


def export_residue_rmsd_csv(residue_rmsds: List[ResidueRMSD], output_path: str):
    """Export structural unit RMSD data to CSV file with proper nomenclature"""
    import csv
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        # Use proper scientific nomenclature in headers
        writer.writerow(['unit_id', 'unit_type', 'chain_id', 'position',
                        'rmsd', 'atom_count', 'molecule_type', 'unit_class'])
        for rmsd in residue_rmsds:
            unit_class = 'amino_acid' if rmsd.is_protein else 'nucleotide'
            writer.writerow([rmsd.unit_id, rmsd.unit_type, rmsd.chain_id,
                           rmsd.position, rmsd.rmsd, rmsd.atom_count, rmsd.molecule_type, unit_class])


def save_aligned_structures(alignment_results: Dict, output_dir, pair_id: str,
                          experimental_structure: BioStructure, predicted_structure: BioStructure):
    """
    Save the aligned structures used for RMSD calculations as superimposed structure files.
    
    Args:
        alignment_results: Dictionary containing alignment results for each frame
        output_dir: Output directory path
        pair_id: Structure pair identifier  
        experimental_structure: Original experimental structure
        predicted_structure: Original predicted structure
    """
    from Bio.PDB import PDBIO
    from Bio.PDB.Structure import Structure as PDBStructure
    from pathlib import Path
    import copy
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    io = PDBIO()
    
    # Create superimposed structure files for each reference frame
    for frame_name, result in alignment_results.items():
        if hasattr(result, 'aligned_predicted_structure') and result.aligned_predicted_structure:
            # Create a combined structure with both experimental (reference) and predicted (aligned)
            combined_structure = PDBStructure('superimposed')
            
            # Add experimental structure as model 0 (reference - will show in one color)
            exp_model = copy.deepcopy(experimental_structure[0])
            exp_model.id = 0
            exp_model.serial_num = 0
            exp_model.full_id = (f"{pair_id}_experimental",)
            combined_structure.add(exp_model)
            
            # Add aligned predicted structure as model 1 (aligned - will show in different color)
            pred_model = copy.deepcopy(result.aligned_predicted_structure[0])
            pred_model.id = 1
            pred_model.serial_num = 1
            pred_model.full_id = (f"{pair_id}_predicted_aligned",)
            combined_structure.add(pred_model)
            
            # Save the superimposed structure with custom MODEL records
            superimposed_file = f"{pair_id}_superimposed_{frame_name}.pdb"
            output_file_path = output_dir / superimposed_file
            
            # Save with custom MODEL identification
            _save_superimposed_with_model_names(combined_structure, str(output_file_path), pair_id, io)
            
            print(f"    Saved superimposed structure: {superimposed_file}")
        else:
            print(f"    Warning: No aligned structure available for {frame_name}")

def _save_superimposed_with_model_names(combined_structure, file_path: str, pair_id: str, io):
    """
    Save the combined structure with custom MODEL records that include descriptive names.
    """
    import tempfile
    import os
    
    # Save to a temporary file first
    io.set_structure(combined_structure)
    
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False) as temp_file:
        temp_file_path = temp_file.name
        
    # Save the structure to temporary file
    io.save(temp_file_path)
        
    # Read the temporary file content
    with open(temp_file_path, 'r') as f:
        pdb_lines = f.readlines()
    
    # Remove temporary file
    os.unlink(temp_file_path)
    
    # Process and write the final file with custom MODEL records
    with open(file_path, 'w') as f:
        # Write header information
        f.write(f"HEADER    SUPERIMPOSED STRUCTURE                    {pair_id}_COMPARISON\n")
        f.write(f"REMARK   1 MODEL 1: {pair_id}_EXPERIMENTAL (reference structure)\n")
        f.write(f"REMARK   1 MODEL 2: {pair_id}_PREDICTED_ALIGNED (superimposed structure)\n")
        f.write("REMARK   1 Models are colored differently in molecular viewers\n")
        
        model_num = 1
        for line in pdb_lines:
            if line.startswith('MODEL'):
                if model_num == 1:
                    f.write(f"MODEL        1 {pair_id}_EXPERIMENTAL\n")
                elif model_num == 2:
                    f.write(f"MODEL        2 {pair_id}_PREDICTED_ALIGNED\n")
                model_num += 1
            else:
                f.write(line)

def _add_structure_identification_header(file_path: str, pair_id: str):
    """
    Add identification headers to the PDB file to distinguish between experimental and predicted structures.
    """
    # Read the original file
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Create header lines to identify the structures
    header_lines = [
        f"HEADER    SUPERIMPOSED STRUCTURE                    {pair_id}_COMPARISON\n",
        f"REMARK   1 MODEL 0: {pair_id}_EXPERIMENTAL (reference structure)\n", 
        f"REMARK   1 MODEL 1: {pair_id}_PREDICTED_ALIGNED (superimposed structure)\n",
        "REMARK   1 Models are colored differently in molecular viewers\n"
    ]
    
    # Write the file with headers
    with open(file_path, 'w') as f:
        # Write headers first
        f.writelines(header_lines)
        # Write original content
        f.writelines(lines)