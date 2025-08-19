#!/usr/bin/env python3
"""
Streamlined structure alignment module for protein-DNA complexes
Supports sequence-based correspondence and three reference frames
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from Bio.PDB import Structure as BioStructure, Superimposer
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
    """Align sequences using sliding window to find optimal correspondence"""
    obs_seq = [res.get_resname().strip() for res in obs_residues]
    pred_seq = [res.get_resname().strip() for res in pred_residues]
    
    # Find best consecutive match
    best_score = 0
    best_obs_start = 0 
    best_pred_start = 0
    
    for obs_start in range(len(obs_seq)):
        for pred_start in range(len(pred_seq)):
            score = 0
            obs_pos, pred_pos = obs_start, pred_start
            
            while (obs_pos < len(obs_seq) and pred_pos < len(pred_seq) and
                   obs_seq[obs_pos] == pred_seq[pred_pos]):
                score += 1
                obs_pos += 1
                pred_pos += 1
            
            if score > best_score:
                best_score = score
                best_obs_start = obs_start
                best_pred_start = pred_start
    
    # Return aligned pairs
    return [(obs_residues[best_obs_start + i], pred_residues[best_pred_start + i]) 
            for i in range(best_score)]


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
        
        # Get alignment atoms based on reference frame
        if ref_frame == 'global':
            # Use both CA and P atoms
            obs_atoms = (get_backbone_atoms(observed, protein_corr, 'CA', is_observed=True) + 
                        get_backbone_atoms(observed, dna_corr, 'P', is_observed=True))
            pred_atoms = (get_backbone_atoms(pred_copy, protein_corr, 'CA', is_observed=False) +
                         get_backbone_atoms(pred_copy, dna_corr, 'P', is_observed=False))
            
        elif ref_frame == 'dna_centric':
            # Use only P atoms
            obs_atoms = get_backbone_atoms(observed, dna_corr, 'P', is_observed=True)
            pred_atoms = get_backbone_atoms(pred_copy, dna_corr, 'P', is_observed=False)
            
        else:  # protein_centric
            # Use only CA atoms  
            obs_atoms = get_backbone_atoms(observed, protein_corr, 'CA', is_observed=True)
            pred_atoms = get_backbone_atoms(pred_copy, protein_corr, 'CA', is_observed=False)
        
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