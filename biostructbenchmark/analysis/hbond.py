"""
biostructbenchmark/analysis/hbond.py
Hydrogen bond analysis for protein-DNA interactions
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import itertools

from Bio.PDB import Structure, Selection, NeighborSearch
from Bio.PDB.vectors import Vector
from biostructbenchmark.core.io import get_structure


@dataclass
class HydrogenBond:
    """Container for hydrogen bond information"""
    donor_atom: str  # Format: chain:residue:atom (e.g., A:SER:4:OG)
    acceptor_atom: str  # Format: chain:residue:atom (e.g., B:DG:15:O6)
    donor_residue: str  # Format: chain:residue_name:position (e.g., A:SER:4)
    acceptor_residue: str  # Format: chain:residue_name:position (e.g., B:DG:15)
    donor_type: str  # 'protein' or 'dna'
    acceptor_type: str  # 'protein' or 'dna'
    distance: float  # Angstroms
    angle: Optional[float] = None  # Degrees (D-H-A angle if hydrogen available)
    interaction_type: str = "protein_dna"  # Type of interaction
    
    @property
    def bond_id(self) -> str:
        """Unique identifier for this hydrogen bond"""
        return f"{self.donor_atom}-->{self.acceptor_atom}"
    
    @property
    def is_protein_dna(self) -> bool:
        """Check if this is a protein-DNA interaction"""
        return (self.donor_type == 'protein' and self.acceptor_type == 'dna') or \
               (self.donor_type == 'dna' and self.acceptor_type == 'protein')


@dataclass
class HBondComparison:
    """Container for hydrogen bond network comparison"""
    experimental_bonds: List[HydrogenBond]
    predicted_bonds: List[HydrogenBond]
    common_bonds: List[Tuple[HydrogenBond, HydrogenBond]]  # (exp, pred) pairs
    experimental_only: List[HydrogenBond]
    predicted_only: List[HydrogenBond]
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
    protein_to_dna_bonds: Dict[str, int]  # experimental vs predicted counts
    dna_to_protein_bonds: Dict[str, int]  # experimental vs predicted counts


class HBondAnalyzer:
    """Analyze hydrogen bonds in protein-DNA complexes"""
    
    def __init__(self, distance_cutoff: float = 3.5, angle_cutoff: float = 120.0):
        """
        Initialize hydrogen bond analyzer
        
        Args:
            distance_cutoff: Maximum distance for hydrogen bond (Angstroms)
            angle_cutoff: Minimum D-H-A angle for hydrogen bond (degrees)
        """
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff
        
        # Define hydrogen bond donors and acceptors
        self.protein_donors = {
            # Backbone
            'N': ['H'],  # Amide nitrogen
            # Side chains
            'OH': ['H'],  # Serine, Threonine, Tyrosine
            'SH': ['H'],  # Cysteine
            'NH': ['H'],  # Asparagine, Glutamine amide
            'NH2': ['H1', 'H2'],  # Asparagine, Glutamine amide
            'NH3': ['H1', 'H2', 'H3'],  # Lysine
            'NE': ['HE'],  # Arginine
            'NH1': ['HH11', 'HH12'],  # Arginine
            'NH2': ['HH21', 'HH22'],  # Arginine
            'NE2': ['HE2'],  # Histidine
            'ND1': ['HD1'],  # Histidine
            'NZ': ['HZ1', 'HZ2', 'HZ3']  # Lysine
        }
        
        self.protein_acceptors = {
            # Backbone
            'O': [],  # Carbonyl oxygen
            # Side chains
            'OD1': [], 'OD2': [],  # Aspartate
            'OE1': [], 'OE2': [],  # Glutamate
            'OG': [], 'OG1': [],  # Serine, Threonine
            'OH': [],  # Tyrosine
            'ND1': [], 'NE2': [],  # Histidine (can be acceptor)
            'SD': []  # Cysteine sulfur
        }
        
        self.dna_donors = {
            'N1': ['H1'],  # Guanine
            'N2': ['H21', 'H22'],  # Guanine
            'N4': ['H41', 'H42'],  # Cytosine
            'N6': ['H61', 'H62']   # Adenine
        }
        
        self.dna_acceptors = {
            'O2': [],  # Cytosine, Thymine
            'O4': [],  # Thymine
            'O6': [],  # Guanine
            'N1': [],  # Adenine (can be acceptor)
            'N3': [],  # Adenine, Cytosine
            'N7': [],  # Adenine, Guanine
            # Phosphate groups
            'O1P': [], 'O2P': [], 'OP1': [], 'OP2': [],
            # Sugar groups
            "O2'": [], "O3'": [], "O4'": [], "O5'": []
        }
    
    def _classify_molecule_type(self, residue) -> str:
        """Classify residue as protein or DNA"""
        dna_residues = {'DA', 'DT', 'DG', 'DC', 'A', 'T', 'G', 'C'}
        protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        resname = residue.get_resname().strip()
        
        if resname in dna_residues:
            return 'dna'
        elif resname in protein_residues:
            return 'protein'
        else:
            # Try to infer from atom names
            atom_names = [atom.get_name() for atom in residue.get_atoms()]
            if any(name.startswith(('P', 'O1P', 'O2P', "O5'", "C5'")) for name in atom_names):
                return 'dna'
            else:
                return 'protein'
    
    def _get_residue_identifier(self, residue) -> str:
        """Get residue identifier in format chain:residue_name:position"""
        chain_id = residue.get_parent().id
        res_name = residue.get_resname().strip()
        res_num = residue.id[1]
        return f"{chain_id}:{res_name}:{res_num}"
    
    def _get_atom_identifier(self, atom) -> str:
        """Get atom identifier in format chain:residue:position:atom"""
        residue = atom.get_parent()
        chain_id = residue.get_parent().id
        res_name = residue.get_resname().strip()
        res_num = residue.id[1]
        atom_name = atom.get_name()
        return f"{chain_id}:{res_name}:{res_num}:{atom_name}"
    
    def _is_hydrogen_bond_geometry(self, donor_atom, acceptor_atom, 
                                   hydrogen_atom=None) -> Tuple[bool, float, Optional[float]]:
        """
        Check if atoms satisfy hydrogen bond geometric criteria
        
        Returns:
            (is_hbond, distance, angle)
        """
        # Calculate distance
        distance = donor_atom - acceptor_atom
        
        if distance > self.distance_cutoff:
            return False, distance, None
        
        # If hydrogen is available, calculate D-H-A angle
        angle = None
        if hydrogen_atom is not None:
            try:
                # Vectors for angle calculation
                dh_vector = hydrogen_atom.get_vector() - donor_atom.get_vector()
                ha_vector = acceptor_atom.get_vector() - hydrogen_atom.get_vector()
                
                # Calculate angle in degrees
                cos_angle = (dh_vector * ha_vector) / (dh_vector.norm() * ha_vector.norm())
                cos_angle = max(-1.0, min(1.0, cos_angle))  # Clamp to [-1, 1]
                angle = np.degrees(np.arccos(cos_angle))
                
                if angle < self.angle_cutoff:
                    return False, distance, angle
                    
            except (ValueError, ZeroDivisionError):
                # If angle calculation fails, rely only on distance
                pass
        
        return True, distance, angle
    
    def find_hydrogen_bonds(self, structure) -> List[HydrogenBond]:
        """
        Find hydrogen bonds in a structure
        
        Args:
            structure: BioPython Structure object
            
        Returns:
            List of HydrogenBond objects
        """
        if not structure:
            return []
        
        hbonds = []
        
        # Get all atoms for neighbor search
        atoms = Selection.unfold_entities(structure, 'A')
        ns = NeighborSearch(atoms)
        
        # Classify all residues
        residue_types = {}
        protein_residues = []
        dna_residues = []
        
        for residue in Selection.unfold_entities(structure, 'R'):
            mol_type = self._classify_molecule_type(residue)
            residue_types[residue] = mol_type
            
            if mol_type == 'protein':
                protein_residues.append(residue)
            elif mol_type == 'dna':
                dna_residues.append(residue)
        
        # Find protein-DNA hydrogen bonds
        for protein_res in protein_residues:
            for dna_res in dna_residues:
                # Check protein donor -> DNA acceptor
                hbonds.extend(self._find_hbonds_between_residues(
                    protein_res, dna_res, 'protein', 'dna', ns))
                
                # Check DNA donor -> protein acceptor
                hbonds.extend(self._find_hbonds_between_residues(
                    dna_res, protein_res, 'dna', 'protein', ns))
        
        return hbonds
    
    def _find_hbonds_between_residues(self, donor_res, acceptor_res, 
                                      donor_type: str, acceptor_type: str,
                                      neighbor_search) -> List[HydrogenBond]:
        """Find hydrogen bonds between two specific residues"""
        hbonds = []
        
        # Get donor and acceptor sets based on molecule type
        if donor_type == 'protein':
            donor_atoms = self.protein_donors
        else:
            donor_atoms = self.dna_donors
            
        if acceptor_type == 'protein':
            acceptor_atoms = self.protein_acceptors
        else:
            acceptor_atoms = self.dna_acceptors
        
        # Check all potential donor-acceptor pairs
        for donor_atom in donor_res.get_atoms():
            donor_name = donor_atom.get_name()
            
            if donor_name not in donor_atoms:
                continue
                
            for acceptor_atom in acceptor_res.get_atoms():
                acceptor_name = acceptor_atom.get_name()
                
                if acceptor_name not in acceptor_atoms:
                    continue
                
                # Look for hydrogen atoms bonded to donor
                hydrogen_atom = None
                potential_hydrogens = donor_atoms[donor_name]
                
                for h_name in potential_hydrogens:
                    try:
                        h_atom = donor_res[h_name]
                        hydrogen_atom = h_atom
                        break
                    except KeyError:
                        continue
                
                # Check geometry
                is_hbond, distance, angle = self._is_hydrogen_bond_geometry(
                    donor_atom, acceptor_atom, hydrogen_atom)
                
                if is_hbond:
                    hbond = HydrogenBond(
                        donor_atom=self._get_atom_identifier(donor_atom),
                        acceptor_atom=self._get_atom_identifier(acceptor_atom),
                        donor_residue=self._get_residue_identifier(donor_res),
                        acceptor_residue=self._get_residue_identifier(acceptor_res),
                        donor_type=donor_type,
                        acceptor_type=acceptor_type,
                        distance=distance,
                        angle=angle,
                        interaction_type="protein_dna"
                    )
                    hbonds.append(hbond)
        
        return hbonds
    
    def compare_hydrogen_bonds(self, experimental_hbonds: List[HydrogenBond],
                              predicted_hbonds: List[HydrogenBond],
                              tolerance: float = 0.5) -> HBondComparison:
        """
        Compare hydrogen bond networks between experimental and predicted structures
        
        Args:
            experimental_hbonds: List of hydrogen bonds from experimental structure
            predicted_hbonds: List of hydrogen bonds from predicted structure
            tolerance: Distance tolerance for matching bonds (Angstroms)
            
        Returns:
            HBondComparison object
        """
        # Create sets for quick lookup
        exp_bond_map = {hb.bond_id: hb for hb in experimental_hbonds}
        pred_bond_map = {hb.bond_id: hb for hb in predicted_hbonds}
        
        # Find exact matches first
        common_bonds = []
        experimental_only = []
        predicted_only = []
        distance_differences = {}
        
        # Check for exact bond ID matches
        for bond_id, exp_hb in exp_bond_map.items():
            if bond_id in pred_bond_map:
                pred_hb = pred_bond_map[bond_id]
                common_bonds.append((exp_hb, pred_hb))
                distance_differences[bond_id] = pred_hb.distance - exp_hb.distance
            else:
                # Try to find similar bonds (same residues, different atoms)
                matched = False
                for pred_bond_id, pred_hb in pred_bond_map.items():
                    if (exp_hb.donor_residue == pred_hb.donor_residue and 
                        exp_hb.acceptor_residue == pred_hb.acceptor_residue):
                        # Similar bond found (same residues)
                        common_bonds.append((exp_hb, pred_hb))
                        distance_differences[f"{bond_id}~{pred_bond_id}"] = \
                            pred_hb.distance - exp_hb.distance
                        matched = True
                        # Remove from predicted map to avoid double matching
                        del pred_bond_map[pred_bond_id]
                        break
                
                if not matched:
                    experimental_only.append(exp_hb)
        
        # Remaining predicted bonds are predicted-only
        for bond_id, pred_hb in pred_bond_map.items():
            if not any(bond_id in pair[1].bond_id or 
                      f"{pair[0].bond_id}~{bond_id}" in distance_differences
                      for pair in common_bonds):
                predicted_only.append(pred_hb)
        
        return HBondComparison(
            experimental_bonds=experimental_hbonds,
            predicted_bonds=predicted_hbonds,
            common_bonds=common_bonds,
            experimental_only=experimental_only,
            predicted_only=predicted_only,
            bond_distance_differences=distance_differences
        )

    def compare_hydrogen_bonds_with_correspondence(self, experimental_hbonds: List[HydrogenBond],
                                                  predicted_hbonds: List[HydrogenBond],
                                                  correspondence_map: Dict[str, str],
                                                  tolerance: float = 0.5) -> HBondComparison:
        """
        Compare hydrogen bond networks using structural correspondence mapping
        
        This method fixes the critical issue where different residue numbering schemes
        between experimental and predicted structures cause false negative matches.
        Instead of direct residue ID matching, it uses the structural alignment
        correspondence to properly identify equivalent bonds.
        
        Args:
            experimental_hbonds: H-bonds from experimental structure
            predicted_hbonds: H-bonds from predicted structure  
            correspondence_map: Maps exp residue IDs to pred residue IDs
            tolerance: Distance tolerance for matching bonds (Angstroms)
            
        Returns:
            HBondComparison with proper correspondence-based matching
        """
        # Create reverse correspondence map (pred -> exp) for efficiency
        reverse_correspondence = {v: k for k, v in correspondence_map.items()}
        
        common_bonds = []
        experimental_only = []
        predicted_only = list(predicted_hbonds)  # Start with all predicted bonds
        distance_differences = {}
        
        # Process each experimental bond
        for exp_hb in experimental_hbonds:
            matched = False
            
            # Check if both donor and acceptor residues have correspondence
            exp_donor_id = exp_hb.donor_residue
            exp_acceptor_id = exp_hb.acceptor_residue
            
            if exp_donor_id in correspondence_map and exp_acceptor_id in correspondence_map:
                # Get corresponding predicted residue IDs
                pred_donor_id = correspondence_map[exp_donor_id]
                pred_acceptor_id = correspondence_map[exp_acceptor_id]
                
                # Look for matching predicted bonds with corresponding residues
                for i, pred_hb in enumerate(predicted_hbonds):
                    if (pred_hb.donor_residue == pred_donor_id and 
                        pred_hb.acceptor_residue == pred_acceptor_id):
                        
                        # Found a correspondence match - check if atoms are similar
                        exp_donor_atom = exp_hb.donor_atom.split(':')[-1]  # Get atom name
                        exp_acceptor_atom = exp_hb.acceptor_atom.split(':')[-1]
                        pred_donor_atom = pred_hb.donor_atom.split(':')[-1] 
                        pred_acceptor_atom = pred_hb.acceptor_atom.split(':')[-1]
                        
                        # Match if same atom types or within distance tolerance
                        if ((exp_donor_atom == pred_donor_atom and exp_acceptor_atom == pred_acceptor_atom) or
                            abs(pred_hb.distance - exp_hb.distance) <= tolerance):
                            
                            # This is a true match
                            common_bonds.append((exp_hb, pred_hb))
                            distance_differences[exp_hb.bond_id] = pred_hb.distance - exp_hb.distance
                            
                            # Remove from predicted_only list
                            if pred_hb in predicted_only:
                                predicted_only.remove(pred_hb)
                            
                            matched = True
                            break
            
            if not matched:
                # Try looser matching - same residue types in correspondence
                for i, pred_hb in enumerate(predicted_hbonds):
                    # Check if residue types match even if exact correspondence missing
                    exp_donor_type = exp_hb.donor_residue.split(':')[1]  # Get residue type
                    exp_acceptor_type = exp_hb.acceptor_residue.split(':')[1]
                    pred_donor_type = pred_hb.donor_residue.split(':')[1]
                    pred_acceptor_type = pred_hb.acceptor_residue.split(':')[1]
                    
                    if (exp_donor_type == pred_donor_type and 
                        exp_acceptor_type == pred_acceptor_type and
                        abs(pred_hb.distance - exp_hb.distance) <= tolerance):
                        
                        # Loose match based on residue types and distance
                        common_bonds.append((exp_hb, pred_hb))
                        distance_differences[f"{exp_hb.bond_id}~loose"] = pred_hb.distance - exp_hb.distance
                        
                        if pred_hb in predicted_only:
                            predicted_only.remove(pred_hb)
                        
                        matched = True
                        break
            
            if not matched:
                experimental_only.append(exp_hb)
        
        return HBondComparison(
            experimental_bonds=experimental_hbonds,
            predicted_bonds=predicted_hbonds,
            common_bonds=common_bonds,
            experimental_only=experimental_only,
            predicted_only=predicted_only,
            bond_distance_differences=distance_differences
        )
    
    def calculate_statistics(self, comparison: HBondComparison) -> HBondStatistics:
        """Calculate summary statistics from hydrogen bond comparison"""
        
        total_exp = len(comparison.experimental_bonds)
        total_pred = len(comparison.predicted_bonds) 
        total_common = len(comparison.common_bonds)
        total_exp_only = len(comparison.experimental_only)
        total_pred_only = len(comparison.predicted_only)
        
        # Calculate conservation rate and prediction accuracy
        conservation_rate = total_common / total_exp if total_exp > 0 else 0.0
        prediction_accuracy = total_common / total_pred if total_pred > 0 else 0.0
        
        # Calculate mean distance difference for common bonds
        distance_diffs = list(comparison.bond_distance_differences.values())
        mean_distance_diff = np.mean(distance_diffs) if distance_diffs else 0.0
        
        # Count directional bonds
        protein_to_dna_exp = sum(1 for hb in comparison.experimental_bonds 
                                if hb.donor_type == 'protein' and hb.acceptor_type == 'dna')
        protein_to_dna_pred = sum(1 for hb in comparison.predicted_bonds
                                 if hb.donor_type == 'protein' and hb.acceptor_type == 'dna')
        
        dna_to_protein_exp = sum(1 for hb in comparison.experimental_bonds
                                if hb.donor_type == 'dna' and hb.acceptor_type == 'protein')
        dna_to_protein_pred = sum(1 for hb in comparison.predicted_bonds
                                 if hb.donor_type == 'dna' and hb.acceptor_type == 'protein')
        
        return HBondStatistics(
            total_experimental=total_exp,
            total_predicted=total_pred,
            total_common=total_common,
            total_experimental_only=total_exp_only,
            total_predicted_only=total_pred_only,
            conservation_rate=conservation_rate,
            prediction_accuracy=prediction_accuracy,
            mean_distance_difference=mean_distance_diff,
            protein_to_dna_bonds={
                'experimental': protein_to_dna_exp,
                'predicted': protein_to_dna_pred
            },
            dna_to_protein_bonds={
                'experimental': dna_to_protein_exp,
                'predicted': dna_to_protein_pred
            }
        )
    
    def analyze_structures(self, experimental_path: Path, predicted_path: Path) -> Tuple[HBondComparison, HBondStatistics]:
        """
        Analyze hydrogen bond networks in experimental vs predicted structures
        
        Args:
            experimental_path: Path to experimental structure
            predicted_path: Path to predicted structure
            
        Returns:
            (HBondComparison, HBondStatistics)
        """
        # Load structures
        exp_structure = get_structure(experimental_path)
        pred_structure = get_structure(predicted_path)
        
        if not exp_structure or not pred_structure:
            raise ValueError(f"Could not load structures from {experimental_path} or {predicted_path}")
        
        # Find hydrogen bonds
        exp_hbonds = self.find_hydrogen_bonds(exp_structure)
        pred_hbonds = self.find_hydrogen_bonds(pred_structure)
        
        # Compare networks
        comparison = self.compare_hydrogen_bonds(exp_hbonds, pred_hbonds)
        statistics = self.calculate_statistics(comparison)
        
        return comparison, statistics

    def analyze_structures_with_correspondence(self, experimental_path: Path, predicted_path: Path, 
                                             correspondence_map: Dict[str, str]) -> Tuple[HBondComparison, HBondStatistics]:
        """
        Analyze hydrogen bond networks using structural correspondence mapping
        
        This method addresses the critical issue where residue numbering differences
        between experimental and predicted structures cause false negative matches
        in hydrogen bond comparison.
        
        Args:
            experimental_path: Path to experimental structure
            predicted_path: Path to predicted structure  
            correspondence_map: Dict mapping experimental residue IDs to predicted residue IDs
                               Format: {exp_chain:res_name:num -> pred_chain:res_name:num}
            
        Returns:
            (HBondComparison, HBondStatistics) with proper correspondence-based matching
        """
        # Load structures
        exp_structure = get_structure(experimental_path)
        pred_structure = get_structure(predicted_path)
        
        if not exp_structure or not pred_structure:
            raise ValueError(f"Could not load structures from {experimental_path} or {predicted_path}")
        
        # Find hydrogen bonds in both structures
        exp_hbonds = self.find_hydrogen_bonds(exp_structure)
        pred_hbonds = self.find_hydrogen_bonds(pred_structure)
        
        # Use correspondence-aware comparison
        comparison = self.compare_hydrogen_bonds_with_correspondence(
            exp_hbonds, pred_hbonds, correspondence_map)
        statistics = self.calculate_statistics(comparison)
        
        return comparison, statistics
    
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
        
        # Export statistics
        self._export_statistics(statistics, output_dir / f"{pair_id}_hbond_statistics.json")
    
    def _export_hbond_details(self, comparison: HBondComparison, output_path: Path):
        """Export detailed hydrogen bond information"""
        
        data = []
        
        # Common bonds
        for exp_hb, pred_hb in comparison.common_bonds:
            data.append({
                'bond_type': 'common',
                'donor_atom_exp': exp_hb.donor_atom,
                'acceptor_atom_exp': exp_hb.acceptor_atom,
                'donor_residue_exp': exp_hb.donor_residue,
                'acceptor_residue_exp': exp_hb.acceptor_residue,
                'distance_exp': exp_hb.distance,
                'angle_exp': exp_hb.angle,
                'donor_atom_pred': pred_hb.donor_atom,
                'acceptor_atom_pred': pred_hb.acceptor_atom,
                'donor_residue_pred': pred_hb.donor_residue,
                'acceptor_residue_pred': pred_hb.acceptor_residue,
                'distance_pred': pred_hb.distance,
                'angle_pred': pred_hb.angle,
                'distance_difference': pred_hb.distance - exp_hb.distance,
                'donor_type': exp_hb.donor_type,
                'acceptor_type': exp_hb.acceptor_type
            })
        
        # Experimental only bonds
        for hb in comparison.experimental_only:
            data.append({
                'bond_type': 'experimental_only',
                'donor_atom_exp': hb.donor_atom,
                'acceptor_atom_exp': hb.acceptor_atom,
                'donor_residue_exp': hb.donor_residue,
                'acceptor_residue_exp': hb.acceptor_residue,
                'distance_exp': hb.distance,
                'angle_exp': hb.angle,
                'donor_atom_pred': None,
                'acceptor_atom_pred': None,
                'donor_residue_pred': None,
                'acceptor_residue_pred': None,
                'distance_pred': None,
                'angle_pred': None,
                'distance_difference': None,
                'donor_type': hb.donor_type,
                'acceptor_type': hb.acceptor_type
            })
        
        # Predicted only bonds
        for hb in comparison.predicted_only:
            data.append({
                'bond_type': 'predicted_only',
                'donor_atom_exp': None,
                'acceptor_atom_exp': None,
                'donor_residue_exp': None,
                'acceptor_residue_exp': None,
                'distance_exp': None,
                'angle_exp': None,
                'donor_atom_pred': hb.donor_atom,
                'acceptor_atom_pred': hb.acceptor_atom,
                'donor_residue_pred': hb.donor_residue,
                'acceptor_residue_pred': hb.acceptor_residue,
                'distance_pred': hb.distance,
                'angle_pred': hb.angle,
                'distance_difference': None,
                'donor_type': hb.donor_type,
                'acceptor_type': hb.acceptor_type
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
    
    def _export_hbond_summary(self, comparison: HBondComparison, statistics: HBondStatistics, 
                             output_path: Path):
        """Export hydrogen bond comparison summary"""
        
        summary_data = [{
            'metric': 'total_experimental_bonds',
            'value': statistics.total_experimental
        }, {
            'metric': 'total_predicted_bonds', 
            'value': statistics.total_predicted
        }, {
            'metric': 'common_bonds',
            'value': statistics.total_common
        }, {
            'metric': 'experimental_only_bonds',
            'value': statistics.total_experimental_only
        }, {
            'metric': 'predicted_only_bonds',
            'value': statistics.total_predicted_only
        }, {
            'metric': 'conservation_rate',
            'value': statistics.conservation_rate
        }, {
            'metric': 'prediction_accuracy',
            'value': statistics.prediction_accuracy
        }, {
            'metric': 'mean_distance_difference',
            'value': statistics.mean_distance_difference
        }, {
            'metric': 'protein_to_dna_experimental',
            'value': statistics.protein_to_dna_bonds['experimental']
        }, {
            'metric': 'protein_to_dna_predicted',
            'value': statistics.protein_to_dna_bonds['predicted']
        }, {
            'metric': 'dna_to_protein_experimental',
            'value': statistics.dna_to_protein_bonds['experimental']
        }, {
            'metric': 'dna_to_protein_predicted',
            'value': statistics.dna_to_protein_bonds['predicted']
        }]
        
        df = pd.DataFrame(summary_data)
        df.to_csv(output_path, index=False)
    
    def _export_statistics(self, statistics: HBondStatistics, output_path: Path):
        """Export statistics as JSON"""
        import json
        import numpy as np
        
        def convert_numpy_types(obj):
            """Convert numpy types to native Python types for JSON serialization"""
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_numpy_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy_types(item) for item in obj]
            else:
                return obj
        
        stats_dict = {
            'total_experimental': statistics.total_experimental,
            'total_predicted': statistics.total_predicted,
            'total_common': statistics.total_common,
            'total_experimental_only': statistics.total_experimental_only,
            'total_predicted_only': statistics.total_predicted_only,
            'conservation_rate': statistics.conservation_rate,
            'prediction_accuracy': statistics.prediction_accuracy,
            'mean_distance_difference': statistics.mean_distance_difference,
            'protein_to_dna_bonds': statistics.protein_to_dna_bonds,
            'dna_to_protein_bonds': statistics.dna_to_protein_bonds
        }
        
        # Convert numpy types to native Python types
        stats_dict = convert_numpy_types(stats_dict)
        
        with open(output_path, 'w') as f:
            json.dump(stats_dict, f, indent=2)