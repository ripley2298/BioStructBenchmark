"""
biostructbenchmark/core/alignment.py
Handles structure superposition with multiple reference frames for comprehensive benchmarking
"""

from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union
import numpy as np
from dataclasses import dataclass
import warnings
import copy

from biostructbenchmark.core.io import get_structure
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB import Structure, PDBIO, Select
from Bio.PDB.Structure import Structure as BioStructure
try:
    from scipy.spatial.distance import cdist
except ImportError:
    def cdist(XA, XB, metric='euclidean'):
        """Simple fallback implementation of cdist for Euclidean distance"""
        import numpy as np
        XA = np.atleast_2d(XA)
        XB = np.atleast_2d(XB)
        distances = np.sqrt(np.sum((XA[:, np.newaxis] - XB[np.newaxis, :])**2, axis=2))
        return distances


@dataclass
class ResidueRMSD:
    """Container for per-residue RMSD data"""
    residue_id: str
    residue_type: str
    chain_id: str
    position: int
    rmsd: float
    atom_count: int
    molecule_type: str  # 'protein' or 'dna_nucleotides'


@dataclass
class ChainRMSD:
    """Container for per-chain RMSD data"""
    chain_id: str
    rmsd: float
    atom_count: int
    residue_count: int
    molecule_types: List[str]  # Types of molecules in this chain


@dataclass
class AlignmentResult:
    """Container for alignment results"""
    overall_rmsd: float
    residue_rmsds: List[ResidueRMSD]
    chain_rmsds: List[ChainRMSD]
    transformation_matrix: np.ndarray
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray
    aligned_atom_count: int
    reference_frame: str  # 'full', 'protein', 'dna_nucleotides'
    aligned_predicted_structure: BioStructure  # Transformed predicted structure
    reference_structure: BioStructure  # Original observed structure
    
    
@dataclass
class MultiFrameAlignmentResult:
    """Container for all three reference frame alignments"""
    full_structure: AlignmentResult
    dna_to_protein: AlignmentResult
    dna_to_dna: AlignmentResult
    
    def get_summary(self) -> Dict:
        """Generate summary statistics"""
        return {
            'full_structure_rmsd': self.full_structure.overall_rmsd,
            'dna_nucleotides_positioning_rmsd': self.dna_to_protein.overall_rmsd,
            'dna_nucleotides_standalone_rmsd': self.dna_to_dna.overall_rmsd,
            'full_atom_count': self.full_structure.aligned_atom_count,
            'dna_nucleotides_atom_count': self.dna_to_dna.aligned_atom_count
        }


class MoleculeSelector(Select):
    """Bio.PDB selector for filtering specific molecule types"""
    
    def __init__(self, molecule_type: str = 'all'):
        """
        Args:
            molecule_type: 'protein', 'dna_nucleotides', or 'all'
        """
        self.molecule_type = molecule_type
    
    def accept_residue(self, residue):
        """Filter residues based on molecule type"""
        if self.molecule_type == 'all':
            return True
        elif self.molecule_type == 'protein':
            return is_protein_residue(residue)
        elif self.molecule_type == 'dna':
            return is_dna_residue(residue)
        return False


def perform_multi_frame_alignment(observed_path: Path, 
                                 predicted_path: Path,
                                 output_dir: Optional[Path] = None) -> Optional[MultiFrameAlignmentResult]:
    """
    Perform all three reference frame alignments and optionally save aligned structures
    
    Args:
        observed_path: Path to experimental structure
        predicted_path: Path to predicted structure  
        output_dir: Optional directory to save aligned structures
        
    Returns:
        MultiFrameAlignmentResult containing all three alignments
    """
    observed = get_structure(observed_path)
    predicted = get_structure(predicted_path)
    
    if observed is None or predicted is None:
        return None
    
    try:
        # Create output directory if specified
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # Alignment 1: Full structure to experimental center of mass
        print("Performing Alignment 1: Full computational structure → experimental structure center of mass")
        full_alignment = align_structures_by_reference_frame(
            observed, predicted, 
            reference_frame='full',
            align_subset='full'
        )
        
        if output_dir:
            save_aligned_structure(
                predicted, 
                output_dir / "aligned_1_full_to_experimental.pdb",
                "Full structure aligned to experimental center of mass"
            )
        
        # Alignment 2: Computational DNA nucleotides to experimental protein center of mass
        print("Performing Alignment 2: Computational DNA nucleotides → experimental protein center of mass")
        dna_to_protein = align_structures_by_reference_frame(
            observed, predicted,
            reference_frame='protein',  # Use protein as reference
            align_subset='dna'  # But calculate RMSD for DNA nucleotides
        )
        
        if output_dir:
            save_aligned_structure(
                predicted,
                output_dir / "aligned_2_dna_to_protein_reference.pdb",
                "DNA nucleotides aligned using protein center of mass as reference"
            )
        
        # Alignment 3: Computational DNA nucleotides to experimental DNA nucleotides center of mass  
        print("Performing Alignment 3: Computational DNA nucleotides → experimental DNA nucleotides center of mass")
        dna_to_dna = align_structures_by_reference_frame(
            observed, predicted,
            reference_frame='dna',
            align_subset='dna'
        )
        
        if output_dir:
            save_aligned_structure(
                predicted,
                output_dir / "aligned_3_dna_to_dna.pdb", 
                "DNA nucleotides aligned to experimental DNA nucleotides center of mass"
            )
            
            # Also save the experimental structure for reference
            save_aligned_structure(
                observed,
                output_dir / "experimental_reference.pdb",
                "Experimental structure (reference)"
            )
        
        # Print summary
        print("\n=== Alignment Summary ===")
        print(f"1. Full structure RMSD: {full_alignment.overall_rmsd:.2f} Å ({full_alignment.aligned_atom_count} atoms)")
        print(f"2. DNA nucleotides positioning RMSD (relative to protein): {dna_to_protein.overall_rmsd:.2f} Å")
        print(f"3. DNA nucleotides standalone RMSD: {dna_to_dna.overall_rmsd:.2f} Å ({dna_to_dna.aligned_atom_count} atoms)")
        
        return MultiFrameAlignmentResult(
            full_structure=full_alignment,
            dna_to_protein=dna_to_protein,
            dna_to_dna=dna_to_dna
        )
        
    except Exception as e:
        print(f"Error in multi-frame alignment: {e}")
        return None


def align_structures_with_structural_superimposition(observed: BioStructure,
                                                   predicted: BioStructure,
                                                   reference_frame: str = 'global') -> AlignmentResult:
    """
    Align structures using three distinct reference frames:
    
    1. 'global': Minimize overall RMSD for both protein and DNA together
    2. 'dna_centric': DNA superimposed to assess DNA model quality independently  
    3. 'protein_centric': Protein superimposed to assess protein-DNA interface prediction
    
    Args:
        observed: Experimental structure (reference)
        predicted: Predicted structure (to be aligned)
        reference_frame: Frame to use for superposition ('global', 'dna_centric', 'protein_centric')
    
    Returns:
        AlignmentResult with structurally aligned coordinates and separate protein/DNA RMSDs
    """
    import copy
    from Bio.PDB import Superimposer
    
    # Create deep copy for transformation
    predicted_copy = copy.deepcopy(predicted)
    
    # Perform structural alignment based on reference frame
    if reference_frame == 'global':
        # Global alignment: minimize overall RMSD for both protein and DNA
        print("Performing global structural alignment (protein + DNA)...")
        
        # Get all backbone atoms (CA for protein, P for DNA)
        backbone_obs = []
        backbone_pred = []
        
        # Add protein CA atoms with sequence-based correspondence  
        protein_correspondence = create_protein_correspondence_map(observed, predicted_copy)
        for obs_full_id, pred_full_id in protein_correspondence.items():
            try:
                obs_chain_id, obs_pos = obs_full_id
                pred_chain_id, pred_pos = pred_full_id
                obs_residue = observed[0][obs_chain_id][obs_pos]
                pred_residue = predicted_copy[0][pred_chain_id][pred_pos]
                if 'CA' in obs_residue and 'CA' in pred_residue:
                    backbone_obs.append(obs_residue['CA'])
                    backbone_pred.append(pred_residue['CA'])
            except KeyError:
                continue
                
        # Add DNA P atoms with sequence correspondence
        correspondence_map, _, _, _ = align_dna_structures_with_phosphates(observed, predicted_copy)
        for obs_full_id, pred_full_id in correspondence_map.items():
            try:
                obs_chain_id, obs_pos = obs_full_id
                pred_chain_id, pred_pos = pred_full_id
                obs_nucleotide = observed[0][obs_chain_id][obs_pos]
                pred_nucleotide = predicted_copy[0][pred_chain_id][pred_pos]
                if 'P' in obs_nucleotide and 'P' in pred_nucleotide:
                    backbone_obs.append(obs_nucleotide['P'])
                    backbone_pred.append(pred_nucleotide['P'])
            except KeyError:
                continue
        
        if len(backbone_obs) >= 3:
            super_imposer = Superimposer()
            super_imposer.set_atoms(backbone_obs, backbone_pred)
            super_imposer.apply(predicted_copy.get_atoms())
            backbone_rmsd = super_imposer.rms
            print(f"Global alignment: {backbone_rmsd:.4f} Å ({len(backbone_obs)} backbone atoms)")
        else:
            raise ValueError("Insufficient backbone atoms for global alignment")
            
    elif reference_frame == 'dna_centric':
        # DNA-centric: superimpose DNA to assess DNA model quality
        print("Performing DNA-centric structural alignment...")
        
        # Get DNA P atoms with sequence correspondence
        dna_correspondence_map, _, _, _ = align_dna_structures_with_phosphates(observed, predicted_copy)
        p_atoms_obs = []
        p_atoms_pred = []
        
        for obs_full_id, pred_full_id in dna_correspondence_map.items():
            try:
                obs_chain_id, obs_pos = obs_full_id
                pred_chain_id, pred_pos = pred_full_id
                obs_nucleotide = observed[0][obs_chain_id][obs_pos]
                pred_nucleotide = predicted_copy[0][pred_chain_id][pred_pos]
                if 'P' in obs_nucleotide and 'P' in pred_nucleotide:
                    p_atoms_obs.append(obs_nucleotide['P'])
                    p_atoms_pred.append(pred_nucleotide['P'])
            except KeyError:
                continue
        
        if len(p_atoms_obs) >= 3:
            super_imposer = Superimposer()
            super_imposer.set_atoms(p_atoms_obs, p_atoms_pred)
            super_imposer.apply(predicted_copy.get_atoms())
            backbone_rmsd = super_imposer.rms
            print(f"DNA-centric alignment: {backbone_rmsd:.4f} Å ({len(p_atoms_obs)} P atoms)")
        else:
            raise ValueError("Insufficient DNA P atoms for DNA-centric alignment")
        
    elif reference_frame == 'protein_centric':
        # Protein-centric: superimpose protein using sequence-based correspondence
        print("Performing protein-centric structural alignment...")
        
        # Get protein CA atoms with sequence correspondence
        protein_correspondence = create_protein_correspondence_map(observed, predicted_copy)
        ca_atoms_obs = []
        ca_atoms_pred = []
        
        for obs_full_id, pred_full_id in protein_correspondence.items():
            try:
                obs_chain_id, obs_pos = obs_full_id
                pred_chain_id, pred_pos = pred_full_id
                obs_residue = observed[0][obs_chain_id][obs_pos]
                pred_residue = predicted_copy[0][pred_chain_id][pred_pos]
                if 'CA' in obs_residue and 'CA' in pred_residue:
                    ca_atoms_obs.append(obs_residue['CA'])
                    ca_atoms_pred.append(pred_residue['CA'])
            except KeyError:
                continue
        
        if len(ca_atoms_obs) >= 3:
            super_imposer = Superimposer()
            super_imposer.set_atoms(ca_atoms_obs, ca_atoms_pred)
            super_imposer.apply(predicted_copy.get_atoms())
            backbone_rmsd = super_imposer.rms
            print(f"Protein-centric alignment: {backbone_rmsd:.4f} Å ({len(ca_atoms_obs)} CA atoms)")
        else:
            raise ValueError("Insufficient CA atoms for protein-centric alignment")
    else:
        raise ValueError(f"Unknown reference frame: {reference_frame}")
    
    # Calculate separate RMSDs for protein and DNA using aligned structure
    protein_rmsd = calculate_rmsd_for_molecule_type(observed, predicted_copy, 'protein')
    dna_rmsd = calculate_rmsd_for_molecule_type(observed, predicted_copy, 'dna')
    
    # Calculate overall RMSD for full structure
    obs_all_atoms = get_atoms_by_molecule_type(observed, 'full')
    pred_all_atoms = get_atoms_by_molecule_type(predicted_copy, 'full')
    overall_rmsd = calculate_rmsd_with_matching_atoms(obs_all_atoms, pred_all_atoms)
    
    # Calculate per-residue RMSDs using aligned structure
    residue_rmsds = calculate_per_residue_rmsd_for_subset(
        observed, predicted_copy, 'full'
    )
    
    # Calculate per-chain RMSDs using aligned structure
    chain_rmsds = calculate_per_chain_rmsd(
        observed, predicted_copy, 'full'
    )
    
    print(f"RMSDs after {reference_frame} alignment:")
    print(f"  Overall: {overall_rmsd:.4f} Å")
    if protein_rmsd is not None:
        print(f"  Protein: {protein_rmsd:.4f} Å")
    if dna_rmsd is not None:
        print(f"  DNA: {dna_rmsd:.4f} Å")
    
    # Create transformation matrix placeholder
    transformation_matrix = np.eye(4)
    rotation_matrix = np.eye(3)
    translation_vector = np.zeros(3)
    
    return AlignmentResult(
        overall_rmsd=overall_rmsd,
        residue_rmsds=residue_rmsds,
        chain_rmsds=chain_rmsds,
        transformation_matrix=transformation_matrix,
        rotation_matrix=rotation_matrix,
        translation_vector=translation_vector,
        aligned_atom_count=len(obs_all_atoms),
        reference_frame=f"{reference_frame}_alignment",
        aligned_predicted_structure=predicted_copy,
        reference_structure=observed
    )


def align_structures_kabsch_backbone(observed: BioStructure,
                                    predicted: BioStructure,
                                    reference_frame: str = 'protein',
                                    align_subset: str = 'full') -> AlignmentResult:
    """
    Align structures using Kabsch algorithm with backbone-based reference frames
    
    Args:
        observed: Experimental structure (reference)
        predicted: Predicted structure (to be aligned)  
        reference_frame: Backbone atoms to use for superposition ('protein', 'dna', 'full')
        align_subset: Subset to calculate final RMSD for ('protein', 'dna', 'full')
    
    Returns:
        AlignmentResult with Kabsch-optimized transformation
    """
    import copy
    
    # Create a deep copy of the predicted structure for transformation
    predicted_copy = copy.deepcopy(predicted)
    
    # Get backbone atoms for reference frame alignment
    obs_backbone = get_backbone_atoms(observed, reference_frame)
    pred_backbone = get_backbone_atoms(predicted_copy, reference_frame)
    
    if not obs_backbone or not pred_backbone:
        raise ValueError(f"No backbone atoms found for reference frame: {reference_frame}")
    
    # Ensure equal number of backbone atoms for alignment
    min_len = min(len(obs_backbone), len(pred_backbone))
    obs_backbone = obs_backbone[:min_len]
    pred_backbone = pred_backbone[:min_len]
    
    # Extract coordinates for Kabsch algorithm
    obs_coords = np.array([atom.get_coord() for atom in obs_backbone])
    pred_coords = np.array([atom.get_coord() for atom in pred_backbone])
    
    # Apply Kabsch algorithm for optimal rotation and translation
    rotation_matrix, translation_vector, backbone_rmsd = kabsch_algorithm(obs_coords, pred_coords)
    
    # Apply transformation to entire predicted structure
    predicted_copy = apply_transformation_to_structure(predicted_copy, rotation_matrix, translation_vector)
    
    # Calculate RMSD for specified subset using transformed structure
    obs_calc_atoms = get_atoms_by_molecule_type(observed, align_subset)
    pred_calc_atoms = get_atoms_by_molecule_type(predicted_copy, align_subset)
    
    # Use identity-based matching for accurate RMSD calculation
    rmsd = calculate_rmsd_with_matching_atoms(obs_calc_atoms, pred_calc_atoms)
    
    # Calculate per-residue RMSDs for the aligned subset using transformed structure
    residue_rmsds = calculate_per_residue_rmsd_for_subset(
        observed, predicted_copy, align_subset
    )
    
    # Calculate per-chain RMSDs using transformed structure
    chain_rmsds = calculate_per_chain_rmsd(
        observed, predicted_copy, align_subset
    )
    
    # Create transformation matrix (4x4 homogeneous coordinates)
    transformation_matrix = np.eye(4)
    transformation_matrix[:3, :3] = rotation_matrix
    transformation_matrix[:3, 3] = translation_vector
    
    return AlignmentResult(
        overall_rmsd=rmsd,
        residue_rmsds=residue_rmsds,
        chain_rmsds=chain_rmsds,
        transformation_matrix=transformation_matrix,
        rotation_matrix=rotation_matrix,
        translation_vector=translation_vector,
        aligned_atom_count=len(obs_calc_atoms),
        reference_frame=f"{reference_frame}_backbone_to_{align_subset}",
        aligned_predicted_structure=predicted_copy,  # Kabsch-transformed structure
        reference_structure=observed  # Original observed structure
    )


def align_structures_by_reference_frame(observed: BioStructure,
                                       predicted: BioStructure,
                                       reference_frame: str = 'full',
                                       align_subset: str = 'full') -> AlignmentResult:
    """
    Align structures by reference frame and return comprehensive results including transformed structures
    
    Args:
        observed: Experimental structure (reference)
        predicted: Predicted structure (to be aligned)
        reference_frame: Frame to use for superposition ('full', 'protein', 'dna')
        align_subset: Subset to calculate RMSD for ('full', 'protein', 'dna')
    
    Returns:
        AlignmentResult with transformed structures for visualization
    """
    import copy
    
    # Create a deep copy of the predicted structure for transformation
    predicted_copy = copy.deepcopy(predicted)
    # ... [previous code] ...

    # Get atoms for alignment (using reference_frame)
    obs_ref_atoms = get_atoms_by_molecule_type(observed, reference_frame)
    pred_ref_atoms = get_atoms_by_molecule_type(predicted_copy, reference_frame)
    
    if not obs_ref_atoms or not pred_ref_atoms:
        raise ValueError(f"No atoms found for reference frame: {reference_frame}")
    
    # Calculate centers of mass (using ALL atoms in reference frame)
    obs_com = calculate_center_of_mass(obs_ref_atoms)
    pred_com = calculate_center_of_mass(pred_ref_atoms)
    
    # Align structures using ALL atoms in reference frame
    min_len = min(len(obs_ref_atoms), len(pred_ref_atoms))
    obs_ref_atoms = obs_ref_atoms[:min_len]
    pred_ref_atoms = pred_ref_atoms[:min_len]
    
    superimposer = Superimposer()
    superimposer.set_atoms(obs_ref_atoms, pred_ref_atoms)
    superimposer.apply(predicted_copy[0].get_atoms())  # Apply to entire predicted structure copy
    
    # Calculate RMSD using identity-based matching (same as per-residue calculation)
    # This ensures we compare the same atoms after superposition
    obs_calc_atoms = get_atoms_by_molecule_type(observed, align_subset)
    pred_calc_atoms = get_atoms_by_molecule_type(predicted_copy, align_subset)  # Use transformed copy
    
    # Use identity-based matching for accurate RMSD calculation
    rmsd = calculate_rmsd_with_matching_atoms(obs_calc_atoms, pred_calc_atoms)
    
    # Calculate per-residue RMSDs for the aligned subset using transformed structure
    residue_rmsds = calculate_per_residue_rmsd_for_subset(
        observed, predicted_copy, align_subset
    )
    
    # Calculate per-chain RMSDs using transformed structure
    chain_rmsds = calculate_per_chain_rmsd(
        observed, predicted_copy, align_subset
    )
    
    # Get transformation details
    rot_matrix = superimposer.rotran[0] if hasattr(superimposer, 'rotran') else np.eye(3)
    trans_vector = superimposer.rotran[1] if hasattr(superimposer, 'rotran') else np.zeros(3)
    
    return AlignmentResult(
        overall_rmsd=rmsd,
        residue_rmsds=residue_rmsds,
        chain_rmsds=chain_rmsds,
        transformation_matrix=np.vstack([
            np.hstack([rot_matrix, trans_vector.reshape(3, 1)]),
            [0, 0, 0, 1]
        ]),
        rotation_matrix=rot_matrix,
        translation_vector=trans_vector,
        aligned_atom_count=len(obs_calc_atoms),
        reference_frame=f"{reference_frame}_to_{align_subset}",
        aligned_predicted_structure=predicted_copy,  # Transformed structure
        reference_structure=observed  # Original observed structure
    )


def get_atoms_by_molecule_type(structure: BioStructure, 
                              molecule_type: str = 'full') -> List:
    """
    Get ALL atoms for specified molecule type (not just backbone), excluding waters and ligands
    
    Args:
        structure: Bio.PDB Structure object
        molecule_type: 'full', 'protein', or 'dna_nucleotides'
        
    Returns:
        List of ALL atoms (not just backbone), excluding waters, ions, and ligands
    """
    atoms = []
    excluded_residues = {
        'HOH', 'WAT', 'H2O',  # Water molecules
        'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'MN',  # Common ions
        'SO4', 'PO4', 'NO3',  # Common anions
        'EDO', 'PEG', 'GOL', 'DMS', 'ACE', 'NH2',  # Common ligands/buffers
        'UNK', 'UNL'  # Unknown residues
    }
    
    if molecule_type == 'full':
        # Get ALL atoms in the entire structure (excluding waters/ligands)
        for chain in structure[0]:
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name not in excluded_residues:
                    for atom in residue.get_atoms():
                        atoms.append(atom)
    elif molecule_type == 'protein':
        # Get ALL atoms from protein residues
        for chain in structure[0]:
            for residue in chain:
                if is_protein_residue(residue):
                    for atom in residue.get_atoms():
                        atoms.append(atom)
    elif molecule_type == 'dna':
        # Get ALL atoms from DNA nucleotides
        for chain in structure[0]:
            for residue in chain:
                if is_dna_residue(residue):
                    for atom in residue.get_atoms():
                        atoms.append(atom)
    return atoms



def calculate_center_of_mass(atoms: List) -> np.ndarray:
    """Calculate center of mass for a list of atoms"""
    if not atoms:
        return np.zeros(3)
    
    coords = np.array([atom.get_coord() for atom in atoms])
    return np.mean(coords, axis=0)


def get_backbone_atoms(structure: BioStructure, molecule_type: str = 'protein') -> List:
    """
    Extract backbone atoms for centroid-based reference frame calculation
    
    Args:
        structure: Bio.PDB Structure object
        molecule_type: 'protein', 'dna', or 'full'
        
    Returns:
        List of backbone atoms (CA for protein, P for DNA)
    """
    backbone_atoms = []
    excluded_residues = {
        'HOH', 'WAT', 'H2O',  # Water molecules
        'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'MN',  # Common ions
        'SO4', 'PO4', 'NO3',  # Common anions
        'EDO', 'PEG', 'GOL', 'DMS', 'ACE', 'NH2',  # Common ligands/buffers
        'UNK', 'UNL'  # Unknown residues
    }
    
    for chain in structure[0]:
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name in excluded_residues:
                continue
                
            if molecule_type == 'full':
                # Include both protein and DNA backbone atoms
                if is_protein_residue(residue):
                    # CA atom for protein
                    if 'CA' in residue:
                        backbone_atoms.append(residue['CA'])
                elif is_dna_residue(residue):
                    # P atom for DNA
                    if 'P' in residue:
                        backbone_atoms.append(residue['P'])
            elif molecule_type == 'protein' and is_protein_residue(residue):
                # CA atom for protein
                if 'CA' in residue:
                    backbone_atoms.append(residue['CA'])
            elif molecule_type == 'dna' and is_dna_residue(residue):
                # P atom for DNA nucleotides
                if 'P' in residue:
                    backbone_atoms.append(residue['P'])
    
    return backbone_atoms


def kabsch_algorithm(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Implement Kabsch algorithm for optimal rotation and translation
    
    Args:
        coords_ref: Reference coordinates (N x 3)
        coords_mobile: Mobile coordinates to align (N x 3)
        
    Returns:
        Tuple of (rotation_matrix, translation_vector, rmsd)
    """
    assert coords_ref.shape == coords_mobile.shape, "Coordinate arrays must have same shape"
    assert coords_ref.shape[1] == 3, "Coordinates must be 3D"
    
    # Center both coordinate sets
    centroid_ref = np.mean(coords_ref, axis=0)
    centroid_mobile = np.mean(coords_mobile, axis=0)
    
    coords_ref_centered = coords_ref - centroid_ref
    coords_mobile_centered = coords_mobile - centroid_mobile
    
    # Compute covariance matrix
    H = coords_mobile_centered.T @ coords_ref_centered
    
    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)
    
    # Calculate rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (determinant = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Calculate translation vector
    t = centroid_ref - R @ centroid_mobile
    
    # Apply transformation and calculate RMSD
    coords_mobile_aligned = (R @ coords_mobile.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((coords_ref - coords_mobile_aligned)**2, axis=1)))
    
    return R, t, rmsd


def apply_transformation_to_structure(structure: BioStructure, rotation_matrix: np.ndarray, 
                                    translation_vector: np.ndarray) -> BioStructure:
    """
    Apply rotation and translation transformation to all atoms in a structure
    
    Args:
        structure: Bio.PDB Structure object
        rotation_matrix: 3x3 rotation matrix
        translation_vector: 3x1 translation vector
        
    Returns:
        Transformed structure (modifies in place and returns copy)
    """
    structure_copy = copy.deepcopy(structure)
    
    for model in structure_copy:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Get current coordinates
                    coord = atom.get_coord()
                    # Apply rotation and translation
                    new_coord = (rotation_matrix @ coord) + translation_vector
                    # Set new coordinates
                    atom.set_coord(new_coord)
    
    return structure_copy


def calculate_rmsd(atoms1, atoms2) -> float:
    """
    Calculate RMSD between two sets of atoms or coordinates
    
    Args:
        atoms1: List of Bio.PDB atoms OR numpy array of coordinates
        atoms2: List of Bio.PDB atoms OR numpy array of coordinates
    
    Returns:
        RMSD value in Angstroms
    """
    import numpy as np
    import warnings
    
    # Check if we're dealing with numpy arrays or atom objects
    if isinstance(atoms1, np.ndarray):
        coords1 = atoms1
    else:
        coords1 = np.array([atom.get_coord() for atom in atoms1])
    
    if isinstance(atoms2, np.ndarray):
        coords2 = atoms2
    else:
        coords2 = np.array([atom.get_coord() for atom in atoms2])
    
    # Ensure same length
    if len(coords1) != len(coords2):
        warnings.warn(f"Coordinate count mismatch: {len(coords1)} vs {len(coords2)}")
        min_len = min(len(coords1), len(coords2))
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]
    
    # Calculate RMSD
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


# Should this be nuked later? This may become a vestigial artifact.
def align_dna_chain_sequences(obs_chain_residues, pred_chain_residues):
    """
    Perform sequence alignment for DNA nucleotide chains to handle numbering mismatches.
    Uses a sliding window approach to find the best local alignment between sequences
    of different lengths.
    
    Args:
        obs_chain_residues: List of residues from observed structure chain
        pred_chain_residues: List of residues from predicted structure chain
        
    Returns:
        List of tuples (obs_residue, pred_residue) for aligned pairs
    """
    aligned_pairs = []
    obs_len = len(obs_chain_residues)
    pred_len = len(pred_chain_residues)
    
    # Convert residues to nucleotide type strings for easier comparison
    obs_sequence = [res.get_resname().strip() for res in obs_chain_residues]
    pred_sequence = [res.get_resname().strip() for res in pred_chain_residues]
    
    print(f"Experimental sequence ({obs_len}): {' '.join(obs_sequence)}")
    print(f"Predicted sequence ({pred_len}): {' '.join(pred_sequence)}")
    
    # Find the best alignment using a sliding window approach
    # Try different starting positions for both sequences
    best_score = 0
    best_obs_start = 0
    best_pred_start = 0
    best_length = 0
    
    # Try all possible starting positions for observed sequence
    for obs_start in range(obs_len):
        # Try all possible starting positions for predicted sequence
        for pred_start in range(pred_len):
            # Calculate how many consecutive matches we can get from these starting positions
            score = 0
            length = 0
            obs_pos = obs_start
            pred_pos = pred_start
            
            while obs_pos < obs_len and pred_pos < pred_len:
                if obs_sequence[obs_pos] == pred_sequence[pred_pos]:
                    score += 1
                    length += 1
                    obs_pos += 1
                    pred_pos += 1
                else:
                    # Stop at first mismatch for this alignment
                    break
            
            # Also try allowing some mismatches (more flexible alignment)
            if score == 0:  # If no exact matches, try flexible alignment
                obs_pos = obs_start
                pred_pos = pred_start
                matches = 0
                total = 0
                
                while obs_pos < obs_len and pred_pos < pred_len and total < min(obs_len - obs_start, pred_len - pred_start):
                    if obs_sequence[obs_pos] == pred_sequence[pred_pos]:
                        matches += 1
                    total += 1
                    obs_pos += 1
                    pred_pos += 1
                
                # Use a score that considers both matches and total length
                if total > 0:
                    score = matches
                    length = total
            
            # Update best alignment if this one is better
            if score > best_score:
                best_score = score
                best_obs_start = obs_start
                best_pred_start = pred_start
                best_length = length
    
    print(f"Best alignment: obs_start={best_obs_start}, pred_start={best_pred_start}, score={best_score}, length={best_length}")
    
    # Create aligned pairs using the best alignment found
    if best_score > 0:
        for i in range(best_length):
            obs_idx = best_obs_start + i
            pred_idx = best_pred_start + i
            
            if obs_idx < obs_len and pred_idx < pred_len:
                obs_res = obs_chain_residues[obs_idx]
                pred_res = pred_chain_residues[pred_idx]
                aligned_pairs.append((obs_res, pred_res))
    
    return aligned_pairs


def create_dna_correspondence_map(observed: BioStructure, predicted: BioStructure) -> dict:
    """
    Create a mapping dictionary for DNA nucleotide correspondence between structures
    
    Args:
        observed: Experimental structure
        predicted: Predicted structure
        
    Returns:
        Dictionary mapping (chain_id, position) from observed to predicted
    """
    correspondence_map = {}
    
    # Excluded residues (waters, ions, ligands)
    excluded_residues = {
        'HOH', 'WAT', 'H2O', 'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'MN',
        'SO4', 'PO4', 'NO3', 'EDO', 'PEG', 'GOL', 'DMS', 'ACE', 'NH2',
        'UNK', 'UNL'
    }
    
    # Group DNA residues by chain
    obs_chains = {}
    pred_chains = {}
    
    for chain in observed[0]:
        chain_id = chain.get_id()
        obs_chains[chain_id] = []
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name not in excluded_residues and is_dna_residue(residue):
                obs_chains[chain_id].append(residue)
    
    for chain in predicted[0]:
        chain_id = chain.get_id()
        pred_chains[chain_id] = []
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name not in excluded_residues and is_dna_residue(residue):
                pred_chains[chain_id].append(residue)
    
    # Create correspondence using sequence alignment
    for chain_id in obs_chains:
        if chain_id in pred_chains and obs_chains[chain_id] and pred_chains[chain_id]:
            aligned_pairs = align_dna_chain_sequences(obs_chains[chain_id], pred_chains[chain_id])
            
            for obs_res, pred_res in aligned_pairs:
                obs_full_id = (obs_res.get_parent().get_id(), obs_res.get_id()[1])
                pred_full_id = (pred_res.get_parent().get_id(), pred_res.get_id()[1])
                correspondence_map[obs_full_id] = pred_full_id
    
    return correspondence_map


def align_dna_structures_with_phosphates(observed: BioStructure, predicted: BioStructure) -> tuple:
    """
    Calculate DNA nucleotide correspondence and P-atom RMSD without structural transformation
    
    This function establishes nucleotide correspondence using sequence alignment and calculates
    P-atom RMSD to assess DNA backbone alignment quality. It does NOT transform coordinates
    since DNA and protein occupy different spatial regions.
    
    Args:
        observed: Experimental structure (reference)
        predicted: Predicted structure (for comparison)
        
    Returns:
        Tuple of (correspondence_map, p_atom_rmsd, aligned_pairs_count, original_predicted)
    """
    import copy
    
    # Get DNA correspondence mapping using our sequence alignment
    correspondence_map = create_dna_correspondence_map(observed, predicted)
    
    # Extract phosphate atoms for aligned nucleotide pairs
    p_coords_obs = []
    p_coords_pred = []
    
    for obs_full_id, pred_full_id in correspondence_map.items():
        try:
            # Get nucleotides using full IDs
            obs_chain_id, obs_pos = obs_full_id
            pred_chain_id, pred_pos = pred_full_id
            
            obs_nucleotide = observed[0][obs_chain_id][obs_pos]
            pred_nucleotide = predicted[0][pred_chain_id][pred_pos]
            
            # Extract phosphate atoms (P) - DNA backbone equivalent of CA
            if 'P' in obs_nucleotide and 'P' in pred_nucleotide:
                p_coords_obs.append(obs_nucleotide['P'].get_coord())
                p_coords_pred.append(pred_nucleotide['P'].get_coord())
            
        except KeyError:
            # Skip if nucleotide not found (shouldn't happen with proper mapping)
            continue
    
    print(f"Found {len(p_coords_obs)} P atom pairs for DNA nucleotide correspondence")
    
    # Calculate P-atom RMSD without transformation (assessment only)
    if len(p_coords_obs) >= 3:
        p_coords_obs = np.array(p_coords_obs)
        p_coords_pred = np.array(p_coords_pred)
        
        # Calculate RMSD between P atoms in their original positions
        diff = p_coords_obs - p_coords_pred
        p_atom_rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        print(f"  P-atom RMSD (original coordinates): {p_atom_rmsd:.4f} Å")
        
        return correspondence_map, p_atom_rmsd, len(p_coords_obs), predicted
    else:
        print(f"Insufficient P atoms ({len(p_coords_obs)}) for DNA nucleotide analysis")
        return correspondence_map, None, 0, predicted


def align_protein_chain_sequences(obs_residues: List, pred_residues: List) -> List[tuple]:
    """
    Align protein sequences using sliding window approach to find optimal correspondence
    
    Args:
        obs_residues: List of observed protein residues
        pred_residues: List of predicted protein residues
        
    Returns:
        List of (observed_residue, predicted_residue) pairs for aligned positions
    """
    # Extract amino acid sequences
    obs_sequence = [res.get_resname().strip() for res in obs_residues]
    pred_sequence = [res.get_resname().strip() for res in pred_residues]
    
    obs_len = len(obs_sequence)
    pred_len = len(pred_sequence)
    
    print(f"Experimental protein sequence ({obs_len}): {' '.join(obs_sequence[:10])}{'...' if obs_len > 10 else ''}")
    print(f"Predicted protein sequence ({pred_len}): {' '.join(pred_sequence[:10])}{'...' if pred_len > 10 else ''}")
    
    # Find the best alignment using a sliding window approach
    best_score = 0
    best_obs_start = 0
    best_pred_start = 0
    best_length = 0
    
    for obs_start in range(obs_len):
        for pred_start in range(pred_len):
            # Calculate consecutive matches from these starting positions
            obs_pos = obs_start
            pred_pos = pred_start
            score = 0
            
            while (obs_pos < obs_len and pred_pos < pred_len and 
                   obs_sequence[obs_pos] == pred_sequence[pred_pos]):
                score += 1
                obs_pos += 1
                pred_pos += 1
            
            if score > best_score:
                best_score = score
                best_obs_start = obs_start
                best_pred_start = pred_start
                best_length = score
    
    print(f"Best alignment: obs_start={best_obs_start}, pred_start={best_pred_start}, score={best_score}, length={best_length}")
    
    # Create aligned pairs using the best alignment
    aligned_pairs = []
    for i in range(best_length):
        obs_idx = best_obs_start + i
        pred_idx = best_pred_start + i
        if obs_idx < obs_len and pred_idx < pred_len:
            aligned_pairs.append((obs_residues[obs_idx], pred_residues[pred_idx]))
    
    return aligned_pairs


def create_protein_correspondence_map(observed: BioStructure, predicted: BioStructure) -> dict:
    """
    Create correspondence map between protein residues using sequence alignment
    
    Args:
        observed: Experimental structure
        predicted: Predicted structure
        
    Returns:
        Dictionary mapping (chain_id, position) from observed to predicted
    """
    correspondence_map = {}
    
    # Protein residue types
    protein_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    # Group protein residues by chain
    obs_chains = {}
    pred_chains = {}
    
    for chain in observed[0]:
        chain_id = chain.get_id()
        obs_chains[chain_id] = []
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name in protein_residues:
                obs_chains[chain_id].append(residue)
    
    for chain in predicted[0]:
        chain_id = chain.get_id()
        pred_chains[chain_id] = []
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name in protein_residues:
                pred_chains[chain_id].append(residue)
    
    # Create correspondence using sequence alignment
    for chain_id in obs_chains:
        if chain_id in pred_chains and obs_chains[chain_id] and pred_chains[chain_id]:
            aligned_pairs = align_protein_chain_sequences(obs_chains[chain_id], pred_chains[chain_id])
            
            print(f"Protein chain {chain_id}: {len(aligned_pairs)} aligned pairs from {len(obs_chains[chain_id])} obs, {len(pred_chains[chain_id])} pred")
            
            for obs_res, pred_res in aligned_pairs:
                obs_full_id = (obs_res.get_parent().get_id(), obs_res.get_id()[1])
                pred_full_id = (pred_res.get_parent().get_id(), pred_res.get_id()[1])
                correspondence_map[obs_full_id] = pred_full_id
    
    return correspondence_map


def calculate_rmsd_for_molecule_type(observed: BioStructure, predicted: BioStructure, molecule_type: str) -> float:
    """
    Calculate RMSD for a specific molecule type (protein or DNA)
    
    Args:
        observed: Experimental structure (reference)
        predicted: Predicted structure (aligned)
        molecule_type: 'protein' or 'dna'
        
    Returns:
        RMSD value for the specified molecule type, or None if no atoms found
    """
    obs_atoms = get_atoms_by_molecule_type(observed, molecule_type)
    pred_atoms = get_atoms_by_molecule_type(predicted, molecule_type)
    
    if not obs_atoms or not pred_atoms:
        return None
        
    return calculate_rmsd_with_matching_atoms(obs_atoms, pred_atoms)


def calculate_nucleotide_rmsd(obs_res, pred_res):
    """
    Calculate RMSD between two individual residues
    
    Args:
        obs_res: Observed residue object
        pred_res: Predicted residue object
        
    Returns:
        ResidueRMSD object or None if calculation fails
    """
    # Get all atoms from both residues/nucleotides (internal: using 'residue' for both protein and DNA)
    obs_atoms_list = list(obs_res.get_atoms())
    pred_atoms_list = list(pred_res.get_atoms())
    
    # Skip if insufficient atoms (minimum 2 required for RMSD calculation)
    if len(obs_atoms_list) < 2 or len(pred_atoms_list) < 2:
        unit_name = "nucleotide" if is_dna_residue(obs_res) else "residue"
        print(f"Insufficient atoms in {unit_name} {obs_res.get_id()}: {len(obs_atoms_list)} vs {len(pred_atoms_list)}")
        return None
    
    # Find matching atom names between the two residues/nucleotides
    obs_atom_names = {atom.get_name(): atom for atom in obs_atoms_list}
    pred_atom_names = {atom.get_name(): atom for atom in pred_atoms_list}
    
    # Get common atom names that exist in both units
    common_atom_names = list(obs_atom_names.keys() & pred_atom_names.keys())
    
    # Skip if no matching atoms found
    if not common_atom_names:
        unit_name = "nucleotides" if is_dna_residue(obs_res) else "residues"
        print(f"No matching atoms found between {unit_name} {obs_res.get_id()}")
        return None
    
    # Extract coordinates for common atoms only (for accurate RMSD)
    obs_coords = np.array([obs_atom_names[name].get_coord() for name in common_atom_names])
    pred_coords = np.array([pred_atom_names[name].get_coord() for name in common_atom_names])
    
    # Calculate RMSD using the matching atoms
    diff = obs_coords - pred_coords
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    # Format unit ID for display (e.g., "A:123")
    chain_id = obs_res.get_parent().get_id()
    position = obs_res.get_id()[1]
    residue_id = f"{chain_id}:{position}"
    
    # Create and return the RMSD result
    residue_rmsd = ResidueRMSD(
        residue_id=residue_id,
        residue_type=obs_res.get_resname(),
        chain_id=chain_id,
        position=position,
        rmsd=rmsd,
        atom_count=len(common_atom_names),
        molecule_type='protein' if is_protein_residue(obs_res) else 'dna'
    )
    
    return residue_rmsd


def calculate_per_residue_rmsd_for_subset(observed: BioStructure,
                                         predicted: BioStructure,
                                         molecule_type: str = 'full') -> List[ResidueRMSD]:
    """
    Calculate per-residue RMSD for a specific molecule type.
    Uses sequence-based alignment for DNA nucleotides to handle numbering differences.
    
    Args:
        observed: Experimental structure
        predicted: Predicted structure  
        molecule_type: 'protein', 'dna_nucleotides', or 'full'
        
    Returns:
        List of ResidueRMSD objects with RMSD values for each matched residue
    """
    residue_rmsds = []
    
    # Get atoms from both structures for the specified molecule type
    obs_atoms = get_atoms_by_molecule_type(observed, molecule_type)
    pred_atoms = get_atoms_by_molecule_type(predicted, molecule_type)
    
    # If no atoms found for this molecule type, return empty list
    if not obs_atoms or not pred_atoms:
        return residue_rmsds
    
    # Excluded residues (waters, ions, ligands)
    excluded_residues = {
        'HOH', 'WAT', 'H2O',  # Water molecules
        'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'MN',  # Common ions
        'SO4', 'PO4', 'NO3',  # Common anions
        'EDO', 'PEG', 'GOL', 'DMS', 'ACE', 'NH2',  # Common ligands/buffers
        'UNK', 'UNL'  # Unknown residues
    }
    
    # For DNA nucleotides, use sequence-based alignment; for proteins, use position-based
    if molecule_type == 'dna' or (molecule_type == 'full' and any(is_dna_residue(res) for chain in observed[0] for res in chain)):
        # Group residues by chain for sequence alignment
        obs_chains = {}
        pred_chains = {}
        
        # Extract DNA residues from observed structure grouped by chain
        for chain in observed[0]:
            chain_id = chain.get_id()
            obs_chains[chain_id] = []
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name in excluded_residues:
                    continue
                if is_dna_residue(residue):
                    obs_chains[chain_id].append(residue)
        
        # Extract DNA residues from predicted structure grouped by chain
        for chain in predicted[0]:
            chain_id = chain.get_id()
            pred_chains[chain_id] = []
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name in excluded_residues:
                    continue
                if is_dna_residue(residue):
                    pred_chains[chain_id].append(residue)
        
        # Align DNA sequences for each matching chain
        for chain_id in obs_chains:
            if chain_id in pred_chains and obs_chains[chain_id] and pred_chains[chain_id]:
                # Perform sequence alignment for this chain
                aligned_pairs = align_dna_chain_sequences(obs_chains[chain_id], pred_chains[chain_id])
                
                print(f"DNA chain {chain_id}: {len(aligned_pairs)} aligned pairs from {len(obs_chains[chain_id])} obs, {len(pred_chains[chain_id])} pred")
                
                # Calculate RMSD for each aligned pair
                for obs_res, pred_res in aligned_pairs:
                    rmsd_result = calculate_nucleotide_rmsd(obs_res, pred_res)
                    if rmsd_result:
                        residue_rmsds.append(rmsd_result)
        
        # If molecule_type is 'full', also process protein residues with sequence-based matching
        if molecule_type == 'full':
            protein_rmsds = calculate_per_residue_rmsd_for_subset(observed, predicted, 'protein')
            residue_rmsds.extend(protein_rmsds)
    else:
        # Use sequence-based correspondence for proteins
        if molecule_type == 'protein':
            # Get protein correspondence mapping using sequence alignment
            protein_correspondence = create_protein_correspondence_map(observed, predicted)
            
            # Calculate RMSD for each aligned protein residue pair
            for obs_full_id, pred_full_id in protein_correspondence.items():
                try:
                    obs_chain_id, obs_pos = obs_full_id
                    pred_chain_id, pred_pos = pred_full_id
                    obs_res = observed[0][obs_chain_id][obs_pos]
                    pred_res = predicted[0][pred_chain_id][pred_pos]
                    
                    # Calculate RMSD for this matched pair
                    rmsd_result = calculate_nucleotide_rmsd(obs_res, pred_res)
                    if rmsd_result:
                        residue_rmsds.append(rmsd_result)
                except KeyError:
                    continue
            
            print(f"Protein sequence alignment: {len(residue_rmsds)} aligned residues")
        
        elif molecule_type == 'full':
            # For 'full' type, process both protein and DNA with respective sequence alignments
            # This case should have been handled above, but fallback to empty
            pass
        
        # Print summary of matching process
        print(f"Residue matching summary - Matched: {matched_count}/{total_count} residues")
    
    return residue_rmsds


class CleanStructureSelector(Select):
    """Bio.PDB selector for excluding waters, ions, and other unwanted molecules"""
    
    def __init__(self, molecule_type: str = 'protein'):
        """
        Args:
            molecule_type: 'protein', 'dna_nucleotides', 'protein_dna', or 'all'
        """
        self.molecule_type = molecule_type
        self.excluded_residues = {
            'HOH', 'WAT', 'H2O',  # Water molecules
            'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'MN',  # Common ions
            'SO4', 'PO4', 'NO3',  # Common anions
            'EDO', 'PEG', 'GOL', 'DMS', 'ACE', 'NH2',  # Common ligands/buffers
            'UNK', 'UNL'  # Unknown residues
        }
    
    def accept_residue(self, residue):
        """Filter residues based on molecule type, excluding waters and ligands"""
        res_name = residue.get_resname().strip()
        
        # Always exclude water and common ligands/ions
        if res_name in self.excluded_residues:
            return False
        
        # Filter by molecule type
        if self.molecule_type == 'all':
            return True
        elif self.molecule_type == 'protein':
            return is_protein_residue(residue)
        elif self.molecule_type == 'dna_nucleotides':
            return is_dna_residue(residue)
        elif self.molecule_type == 'protein_dna':
            return is_protein_residue(residue) or is_dna_residue(residue)
        
        return False


def save_aligned_structure(structure: BioStructure, 
                          output_path: Path,
                          description: str = "",
                          molecule_type: str = 'protein_dna') -> None:
    """
    Save aligned structure to PDB file with description, excluding waters and ligands
    
    Args:
        structure: Structure to save
        output_path: Output file path
        description: Description for logging
        molecule_type: Type of molecules to include ('protein', 'dna_nucleotides', 'protein_dna', 'all')
    """
    io = PDBIO()
    io.set_structure(structure)
    
    # Use selector to exclude waters and unwanted molecules
    selector = CleanStructureSelector(molecule_type)
    io.save(str(output_path), selector)
    
    print(f"Saved: {output_path}")
    if description:
        print(f"  Description: {description}")
    print(f"  Excluded: waters, ions, and ligands")


def is_protein_residue(residue) -> bool:
    """Check if residue is a standard amino acid"""
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
        'THR', 'TRP', 'TYR', 'VAL'
    }
    return residue.get_resname().strip() in standard_aa


def is_dna_residue(residue) -> bool:
    """Check if residue is a DNA nucleotide"""
    dna_bases = {'DA', 'DT', 'DG', 'DC', 'A', 'T', 'G', 'C',
                 'DI', 'DU', 'OA', 'OT', 'OM', 'OB', 'OP', 'CL'}
    # Include alternative representations and modified bases
    res_name = residue.get_resname().strip()
    return res_name.upper() in dna_bases or res_name.lower() in [base.lower() for base in dna_bases]


def get_common_atoms(res1, res2) -> List[str]:
    """Get list of common atom names between two residues (all atoms)"""
    atoms1 = {atom.get_name() for atom in res1}
    atoms2 = {atom.get_name() for atom in res2}
    return list(atoms1 & atoms2)


def calculate_rmsd_with_matching_atoms(atoms1, atoms2) -> float:
    """Calculate RMSD between two sets of atoms with matching names"""
    atom_names1 = {atom.get_name(): atom for atom in atoms1}
    atom_names2 = {atom.get_name(): atom for atom in atoms2}

    common_atom_names = list(atom_names1.keys() & atom_names2.keys())
    if not common_atom_names:
        return 0.0

    coords1 = np.array([atom_names1[name].get_coord() for name in common_atom_names])
    coords2 = np.array([atom_names2[name].get_coord() for name in common_atom_names])

    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))



def calculate_per_chain_rmsd(observed: BioStructure,
                            predicted: BioStructure,
                            molecule_type: str = 'full') -> List[ChainRMSD]:
    """
    Calculate per-chain RMSD values.
    
    Args:
        observed: Experimental structure
        predicted: Predicted structure  
        molecule_type: 'protein', 'dna_nucleotides', or 'full'
        
    Returns:
        List of ChainRMSD objects with RMSD values for each chain
    """
    chain_rmsds = []
    
    # Get all chain IDs from observed structure
    obs_chain_ids = set()
    for chain in observed[0]:
        obs_chain_ids.add(chain.get_id())
    
    # Calculate RMSD for each chain
    for chain_id in obs_chain_ids:
        obs_chain_atoms = []
        pred_chain_atoms = []
        residue_count = 0
        molecule_types = set()
        
        # Get atoms from observed chain
        try:
            obs_chain = observed[0][chain_id]
            for residue in obs_chain:
                if molecule_type == 'full' or \
                   (molecule_type == 'protein' and is_protein_residue(residue)) or \
                   (molecule_type == 'dna_nucleotides' and is_dna_residue(residue)):
                    residue_count += 1
                    if is_protein_residue(residue):
                        molecule_types.add('protein')
                    elif is_dna_residue(residue):
                        molecule_types.add('dna_nucleotides')
                    else:
                        molecule_types.add('other')
                    
                    for atom in residue.get_atoms():
                        obs_chain_atoms.append(atom)
        except KeyError:
            continue
        
        # Get atoms from predicted chain
        try:
            pred_chain = predicted[0][chain_id]
            for residue in pred_chain:
                if molecule_type == 'full' or \
                   (molecule_type == 'protein' and is_protein_residue(residue)) or \
                   (molecule_type == 'dna_nucleotides' and is_dna_residue(residue)):
                    
                    for atom in residue.get_atoms():
                        pred_chain_atoms.append(atom)
        except KeyError:
            continue
        
        # Calculate RMSD for this chain if atoms exist
        if obs_chain_atoms and pred_chain_atoms:
            chain_rmsd = calculate_rmsd_with_matching_atoms(obs_chain_atoms, pred_chain_atoms)
            
            chain_rmsds.append(ChainRMSD(
                chain_id=chain_id,
                rmsd=chain_rmsd,
                atom_count=min(len(obs_chain_atoms), len(pred_chain_atoms)),
                residue_count=residue_count,
                molecule_types=list(molecule_types)
            ))
    
    return chain_rmsds


def export_residue_rmsd_csv(residue_rmsds: List[ResidueRMSD], 
                           output_path: Path,
                           reference_frame: str = "") -> None:
    """Export per-residue RMSD data to CSV"""
    import csv

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['residue_id', 'residue_type', 'chain_id', 
                        'position', 'rmsd', 'atom_count', 'molecule_type',
                        'reference_frame'])
        
        for r in residue_rmsds:
            writer.writerow([r.residue_id, r.residue_type, r.chain_id,
                           r.position, r.rmsd, r.atom_count, r.molecule_type,
                           reference_frame])
    
    print(f"Exported per-residue RMSD data to: {output_path}")


def export_chain_rmsd_csv(chain_rmsds: List[ChainRMSD], 
                         output_path: Path,
                         reference_frame: str = "") -> None:
    """Export per-chain RMSD data to CSV"""
    import csv

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['chain_id', 'rmsd', 'atom_count', 'residue_count', 
                        'molecule_types', 'reference_frame'])
        
        for c in chain_rmsds:
            writer.writerow([c.chain_id, c.rmsd, c.atom_count, c.residue_count,
                           ';'.join(c.molecule_types), reference_frame])
    
    print(f"Exported per-chain RMSD data to: {output_path}")


def export_comprehensive_alignment_report(result: AlignmentResult,
                                         output_dir: Path,
                                         prefix: str = "alignment") -> None:
    """
    Export comprehensive alignment report with all outputs
    
    Args:
        result: AlignmentResult object containing transformed structures
        output_dir: Directory to save outputs
        prefix: Prefix for output filenames
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Export aligned structures as PDB files
    exp_output = output_dir / f"{prefix}_experimental.pdb"
    pred_output = output_dir / f"{prefix}_predicted_aligned.pdb"
    superimposition_output = output_dir / f"{prefix}_superimposition.pdb"
    
    # Use the transformed structures from the AlignmentResult
    save_aligned_structure(result.reference_structure, exp_output, "Experimental structure (reference)")
    save_aligned_structure(result.aligned_predicted_structure, pred_output, "Predicted structure (aligned)")
    
    # Export combined superimposition for visual validation
    save_superimposition_pdb(result, superimposition_output, f"Superimposition: {result.reference_frame}")
    
    # Export per-residue RMSD data
    residue_csv = output_dir / f"{prefix}_per_residue_rmsd.csv"
    export_residue_rmsd_csv(result.residue_rmsds, residue_csv, result.reference_frame)
    
    # Export per-chain RMSD data
    chain_csv = output_dir / f"{prefix}_per_chain_rmsd.csv"
    export_chain_rmsd_csv(result.chain_rmsds, chain_csv, result.reference_frame)
    
    # Export summary report
    summary_file = output_dir / f"{prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("STRUCTURE ALIGNMENT COMPREHENSIVE REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Reference Frame: {result.reference_frame}\n")
        f.write(f"Overall RMSD: {result.overall_rmsd:.3f} Å\n")
        f.write(f"Total Aligned Atoms: {result.aligned_atom_count}\n\n")
        
        f.write("PER-CHAIN RMSD VALUES:\n")
        f.write("-" * 40 + "\n")
        for chain in result.chain_rmsds:
            mol_types = ', '.join(chain.molecule_types)
            f.write(f"Chain {chain.chain_id}: {chain.rmsd:.3f} Å "
                   f"({chain.atom_count} atoms, {chain.residue_count} residues, {mol_types})\n")
        
        f.write(f"\nPER-RESIDUE RMSD STATISTICS:\n")
        f.write("-" * 40 + "\n")
        if result.residue_rmsds:
            rmsds = [r.rmsd for r in result.residue_rmsds]
            f.write(f"Number of residues: {len(result.residue_rmsds)}\n")
            f.write(f"Mean RMSD: {np.mean(rmsds):.3f} Å\n")
            f.write(f"Median RMSD: {np.median(rmsds):.3f} Å\n")
            f.write(f"Min RMSD: {np.min(rmsds):.3f} Å\n")
            f.write(f"Max RMSD: {np.max(rmsds):.3f} Å\n")
            f.write(f"Std Dev: {np.std(rmsds):.3f} Å\n")
            
            # Count by molecule type
            protein_count = sum(1 for r in result.residue_rmsds if r.molecule_type == 'protein')
            dna_count = sum(1 for r in result.residue_rmsds if r.molecule_type == 'dna_nucleotides')
            f.write(f"\nMolecule type breakdown:\n")
            f.write(f"  Protein residues: {protein_count}\n")
            f.write(f"  DNA nucleotides: {dna_count}\n")
            
            # RMSD quality assessment
            under_2A = sum(1 for r in result.residue_rmsds if r.rmsd < 2.0)
            under_5A = sum(1 for r in result.residue_rmsds if r.rmsd < 5.0)
            f.write(f"\nRMSD Quality Assessment:\n")
            f.write(f"  Residues with RMSD < 2.0 Å: {under_2A} ({under_2A/len(result.residue_rmsds)*100:.1f}%)\n")
            f.write(f"  Residues with RMSD < 5.0 Å: {under_5A} ({under_5A/len(result.residue_rmsds)*100:.1f}%)\n")
        
        f.write(f"\nOUTPUT FILES:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Experimental structure: {exp_output.name}\n")
        f.write(f"Aligned predicted structure: {pred_output.name}\n")
        f.write(f"Combined superimposition: {superimposition_output.name}\n")
        f.write(f"Per-residue RMSD data: {residue_csv.name}\n")
        f.write(f"Per-chain RMSD data: {chain_csv.name}\n")
        f.write(f"This summary report: {summary_file.name}\n")
    
    print(f"\nComprehensive alignment report exported to: {output_dir}")
    print(f"Summary report: {summary_file}")


def save_superimposition_pdb(result: AlignmentResult, 
                            output_path: Path,
                            description: str = "Structural superimposition") -> None:
    """
    Save a combined PDB file showing both experimental and aligned predicted structures
    for visual validation of the superimposition.
    
    Args:
        result: AlignmentResult containing both structures
        output_path: Path to save the combined PDB file
        description: Description for the PDB header
    """
    from Bio.PDB import PDBIO, Structure, Model
    import copy
    
    # Create a new structure to hold both models
    combined_structure = Structure.Structure("superimposition")
    
    # Add experimental structure as Model 0
    exp_model = copy.deepcopy(result.reference_structure[0])
    exp_model.id = 0
    exp_model.detach_parent()
    combined_structure.add(exp_model)
    
    # Add aligned predicted structure as Model 1
    pred_model = copy.deepcopy(result.aligned_predicted_structure[0])
    pred_model.id = 1  
    pred_model.detach_parent()
    combined_structure.add(pred_model)
    
    # Write the combined structure
    io = PDBIO()
    io.set_structure(combined_structure)
    
    # Use our CleanStructureSelector to include all biological molecules but exclude waters/ligands
    selector = CleanStructureSelector(molecule_type='all')
    io.save(str(output_path), select=selector)
    
    print(f"Saved superimposition: {output_path}")
    print(f"  Description: {description}")
    print(f"  Model 0: Experimental structure (reference)")
    print(f"  Model 1: Predicted structure (aligned)")
    print(f"  Included: protein and DNA nucleotides")
    print(f"  Excluded: waters, ions, and ligands")


def save_component_specific_superimpositions(result: AlignmentResult,
                                           output_dir: Path,
                                           prefix: str) -> None:
    """
    Save component-specific superimposition files for better visualization of
    different molecular components (protein, DNA chains, etc.)
    
    Args:
        result: AlignmentResult containing both structures
        output_dir: Directory to save the component-specific files
        prefix: Prefix for output filenames
    """
    from Bio.PDB import PDBIO, Structure, Model, Chain
    import copy
    
    # Get unique chain IDs and their molecule types from chain RMSDs
    chain_info = {}
    for chain_rmsd in result.chain_rmsds:
        chain_id = chain_rmsd.chain_id
        molecule_type = chain_rmsd.molecule_types
        chain_info[chain_id] = {
            'molecule_type': molecule_type,
            'rmsd': chain_rmsd.rmsd,
            'atom_count': chain_rmsd.atom_count
        }
    
    # Create component-specific superimposition files
    for chain_id, info in chain_info.items():
        molecule_type = info['molecule_type']
        rmsd = info['rmsd']
        
        # Create filename based on chain and molecule type
        if molecule_type == 'protein':
            component_file = output_dir / f"{prefix}_superimposition_protein_chain_{chain_id}.pdb"
        elif molecule_type == 'dna_nucleotides':
            component_file = output_dir / f"{prefix}_superimposition_dna_chain_{chain_id}.pdb"
        else:
            component_file = output_dir / f"{prefix}_superimposition_chain_{chain_id}.pdb"
        
        # Create a new structure with only this chain from both models
        component_structure = Structure.Structure(f"superimposition_chain_{chain_id}")
        
        # Add experimental chain as Model 0
        if chain_id in [chain.id for chain in result.reference_structure[0]]:
            exp_model = Model.Model(0)
            exp_chain = copy.deepcopy(result.reference_structure[0][chain_id])
            exp_chain.detach_parent()
            exp_model.add(exp_chain)
            component_structure.add(exp_model)
        
        # Add aligned predicted chain as Model 1
        if chain_id in [chain.id for chain in result.aligned_predicted_structure[0]]:
            pred_model = Model.Model(1)
            pred_chain = copy.deepcopy(result.aligned_predicted_structure[0][chain_id])
            pred_chain.detach_parent()
            pred_model.add(pred_chain)
            component_structure.add(pred_model)
        
        # Write the component-specific structure
        io = PDBIO()
        io.set_structure(component_structure)
        
        # Use our CleanStructureSelector to exclude waters/ligands
        selector = CleanStructureSelector()
        io.save(str(component_file), select=selector)
        
        print(f"Saved component superimposition: {component_file.name}")
        print(f"  Chain {chain_id} ({molecule_type}): RMSD = {rmsd:.3f} Å")
        print(f"  Model 0: Experimental chain {chain_id}")
        print(f"  Model 1: Predicted chain {chain_id} (aligned)")


# Backwards compatibility functions
def compare_structures(observed_path: Path, predicted_path: Path) -> Optional[AlignmentResult]:
    """Legacy function - performs full structure alignment"""
    observed = get_structure(observed_path)
    predicted = get_structure(predicted_path)
    
    if observed is None or predicted is None:
        return None
    
    return align_structures_by_reference_frame(observed, predicted, 'full', 'full')


def align_structures(observed: BioStructure, predicted: BioStructure) -> Dict:
    """Legacy function - performs basic alignment"""
    result = align_structures_by_reference_frame(observed, predicted, 'full', 'full')
    
    return {
        'rmsd': result.overall_rmsd,
        'transformation': result.transformation_matrix,
        'rotation': result.rotation_matrix,
        'translation': result.translation_vector,
        'atom_count': result.aligned_atom_count,
        'chain_rmsds': result.chain_rmsds,
        'residue_rmsds': result.residue_rmsds
    }
