"""
biostructbenchmark/core/alignment.py
Handles structure superposition with multiple reference frames for comprehensive benchmarking
"""

from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union
import numpy as np
from dataclasses import dataclass
import warnings

from biostructbenchmark.core.io import get_structure
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB import Structure, PDBIO, Select
from Bio.PDB.Structure import Structure as BioStructure


@dataclass
class ResidueRMSD:
    """Container for per-residue RMSD data"""
    residue_id: str
    residue_type: str
    chain_id: str
    position: int
    rmsd: float
    atom_count: int
    molecule_type: str  # 'protein' or 'dna'


@dataclass
class AlignmentResult:
    """Container for alignment results"""
    overall_rmsd: float
    residue_rmsds: List[ResidueRMSD]
    transformation_matrix: np.ndarray
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray
    aligned_atom_count: int
    reference_frame: str  # 'full', 'protein', 'dna'
    
    
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
            'dna_positioning_rmsd': self.dna_to_protein.overall_rmsd,
            'dna_standalone_rmsd': self.dna_to_dna.overall_rmsd,
            'full_atom_count': self.full_structure.aligned_atom_count,
            'dna_atom_count': self.dna_to_dna.aligned_atom_count
        }


class MoleculeSelector(Select):
    """Bio.PDB selector for filtering specific molecule types"""
    
    def __init__(self, molecule_type: str = 'all'):
        """
        Args:
            molecule_type: 'protein', 'dna', or 'all'
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
        
        # Alignment 2: Computational DNA to experimental protein center of mass
        print("Performing Alignment 2: Computational DNA → experimental protein center of mass")
        dna_to_protein = align_structures_by_reference_frame(
            observed, predicted,
            reference_frame='protein',  # Use protein as reference
            align_subset='dna'  # But calculate RMSD for DNA
        )
        
        if output_dir:
            save_aligned_structure(
                predicted,
                output_dir / "aligned_2_dna_to_protein_reference.pdb",
                "DNA aligned using protein center of mass as reference"
            )
        
        # Alignment 3: Computational DNA to experimental DNA center of mass  
        print("Performing Alignment 3: Computational DNA → experimental DNA center of mass")
        dna_to_dna = align_structures_by_reference_frame(
            observed, predicted,
            reference_frame='dna',
            align_subset='dna'
        )
        
        if output_dir:
            save_aligned_structure(
                predicted,
                output_dir / "aligned_3_dna_to_dna.pdb", 
                "DNA aligned to experimental DNA center of mass"
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
        print(f"2. DNA positioning RMSD (relative to protein): {dna_to_protein.overall_rmsd:.2f} Å")
        print(f"3. DNA standalone RMSD: {dna_to_dna.overall_rmsd:.2f} Å ({dna_to_dna.aligned_atom_count} atoms)")
        
        return MultiFrameAlignmentResult(
            full_structure=full_alignment,
            dna_to_protein=dna_to_protein,
            dna_to_dna=dna_to_dna
        )
        
    except Exception as e:
        print(f"Error in multi-frame alignment: {e}")
        return None


def align_structures_by_reference_frame(observed: BioStructure,
                                       predicted: BioStructure,
                                       reference_frame: str = 'full',
                                       align_subset: str = 'full') -> AlignmentResult:
    # ... [previous code] ...

    # Get atoms for alignment (using reference_frame)
    obs_ref_atoms = get_atoms_by_molecule_type(observed, reference_frame)
    pred_ref_atoms = get_atoms_by_molecule_type(predicted, reference_frame)
    
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
    superimposer.apply(predicted[0].get_atoms())  # Apply to entire predicted structure
    
    # Calculate RMSD using ALL atoms for the specified subset (no special case for 'full')
    obs_calc_atoms = get_atoms_by_molecule_type(observed, align_subset)
    pred_calc_atoms = get_atoms_by_molecule_type(predicted, align_subset)
    
    # Ensure we have matching atom counts
    min_calc = min(len(obs_calc_atoms), len(pred_calc_atoms))
    obs_calc_atoms = obs_calc_atoms[:min_calc]
    pred_calc_atoms = pred_calc_atoms[:min_calc]
    
    rmsd = calculate_rmsd(obs_calc_atoms, pred_calc_atoms)
    
    # Calculate per-residue RMSDs for the aligned subset
    residue_rmsds = calculate_per_residue_rmsd_for_subset(
        observed, predicted, align_subset
    )
    
    # Get transformation details
    rot_matrix = superimposer.rotran[0] if hasattr(superimposer, 'rotran') else np.eye(3)
    trans_vector = superimposer.rotran[1] if hasattr(superimposer, 'rotran') else np.zeros(3)
    
    return AlignmentResult(
        overall_rmsd=rmsd,
        residue_rmsds=residue_rmsds,
        transformation_matrix=np.vstack([
            np.hstack([rot_matrix, trans_vector.reshape(3, 1)]),
            [0, 0, 0, 1]
        ]),
        rotation_matrix=rot_matrix,
        translation_vector=trans_vector,
        aligned_atom_count=len(obs_calc_atoms),
        reference_frame=f"{reference_frame}_to_{align_subset}"
    )


def get_atoms_by_molecule_type(structure: BioStructure, 
                              molecule_type: str = 'full') -> List:
    """
    Get ALL atoms for specified molecule type (not just backbone)
    
    Args:
        structure: Bio.PDB Structure object
        molecule_type: 'full', 'protein', or 'dna'
        
    Returns:
        List of ALL atoms (not just backbone)
    """
    atoms = []
    if molecule_type == 'full':
        # Get ALL atoms in the entire structure
        for chain in structure[0]:
            for residue in chain:
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
        # Get ALL atoms from DNA residues
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
def calculate_per_residue_rmsd_for_subset(observed: BioStructure,
                                         predicted: BioStructure,
                                         molecule_type: str = 'full') -> List[ResidueRMSD]:
    """
    Calculate per-residue RMSD for a specific molecule type.
    
    Args:
        observed: Experimental structure
        predicted: Predicted structure  
        molecule_type: 'protein', 'dna', or 'full'
        
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
    
    # Create dictionaries mapping (resname, position, chain_id) -> residue
    # This enables identity-based matching instead of index-based matching
    obs_residues = []
    pred_residues = []
    
    # Extract all residues from observed structure that match the molecule type
    for chain in observed[0]:
        for residue in chain:
            if molecule_type == 'full' or \
               (molecule_type == 'protein' and is_protein_residue(residue)) or \
               (molecule_type == 'dna' and is_dna_residue(residue)):
                obs_residues.append(residue)
    
    # Extract all residues from predicted structure that match the molecule type
    for chain in predicted[0]:
        for residue in chain:
            if molecule_type == 'full' or \
               (molecule_type == 'protein' and is_protein_residue(residue)) or \
               (molecule_type == 'dna' and is_dna_residue(residue)):
                pred_residues.append(residue)
    
    # Create lookup dictionaries for identity-based matching
    # Key: (residue_name, residue_position, chain_id) -> residue object
    obs_res_dict = {}
    for res in obs_residues:
        key = (res.get_resname(), res.get_id()[1], res.get_parent().get_id())
        obs_res_dict[key] = res
    
    pred_res_dict = {}
    for res in pred_residues:
        key = (res.get_resname(), res.get_id()[1], res.get_parent().get_id())
        pred_res_dict[key] = res
    
    # Match residues by identity (name, position, chain) instead of index
    matched_count = 0
    total_count = len(obs_res_dict)
    
    # Iterate through all observed residues and find matches in predicted
    for obs_key, obs_res in obs_res_dict.items():
        if obs_key in pred_res_dict:
            pred_res = pred_res_dict[obs_key]
            matched_count += 1
            
            # Get all atoms from both residues
            obs_atoms_list = list(obs_res.get_atoms())
            pred_atoms_list = list(pred_res.get_atoms())
            
            # Skip if insufficient atoms (minimum 2 required for RMSD calculation)
            if len(obs_atoms_list) < 2 or len(pred_atoms_list) < 2:
                print(f"Insufficient atoms in residue {obs_res.get_id()}: {len(obs_atoms_list)} vs {len(pred_atoms_list)}")
                continue
            
            # Find matching atom names between the two residues
            obs_atom_names = {atom.get_name(): atom for atom in obs_atoms_list}
            pred_atom_names = {atom.get_name(): atom for atom in pred_atoms_list}
            
            # Get common atom names that exist in both residues
            common_atom_names = list(obs_atom_names.keys() & pred_atom_names.keys())
            
            # Skip if no matching atoms found
            if not common_atom_names:
                print(f"No matching atoms found between residues {obs_res.get_id()}")
                continue
            
            # Extract coordinates for common atoms only (for accurate RMSD)
            obs_coords = np.array([obs_atom_names[name].get_coord() for name in common_atom_names])
            pred_coords = np.array([pred_atom_names[name].get_coord() for name in common_atom_names])
            
            # Calculate RMSD using the matching atoms
            diff = obs_coords - pred_coords
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            # Format residue ID for display (e.g., "A:123")
            chain_id = obs_res.get_parent().get_id()
            position = obs_res.get_id()[1]
            residue_id = f"{chain_id}:{position}"
            
            # Create and store the residue RMSD result
            residue_rmsd = ResidueRMSD(
                residue_id=residue_id,
                residue_type=obs_res.get_resname(),
                chain_id=chain_id,
                position=position,
                rmsd=rmsd,
                atom_count=len(common_atom_names),
                molecule_type='protein' if is_protein_residue(obs_res) else 'dna'
            )
            residue_rmsds.append(residue_rmsd)
        else:
            # Optional: Print unmatched residues for debugging
            print(f"Unmatched residue in observed structure: {obs_key}")
    
    # Print summary of matching process
    print(f"Residue matching summary - Matched: {matched_count}/{total_count} residues")
    
    return residue_rmsds


def save_aligned_structure(structure: BioStructure, 
                          output_path: Path,
                          description: str = "") -> None:
    """Save aligned structure to PDB file with description"""
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))
    
    print(f"Saved: {output_path}")
    if description:
        print(f"  Description: {description}")


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
        'atom_count': result.aligned_atom_count
    }
