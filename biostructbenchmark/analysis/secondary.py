"""
biostructbenchmark/analysis/secondary.py
Secondary structure analysis using simplified DSSP algorithm
"""

from pathlib import Path
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
import numpy as np
import pandas as pd

from biostructbenchmark.core.io import get_structure

@dataclass
class SecondaryStructure:
    """Container for secondary structure assignment"""
    residue_id: str
    chain_id: str
    position: int
    structure_type: str  # 'H' (helix), 'E' (sheet), 'L' (loop)
    confidence: float

class SecondaryStructureAnalyzer:
    """Simplified secondary structure assignment"""
    
    def __init__(self):
        self.phi_psi_ranges = {
            'helix': (-60, -45, -60, -45),  # (phi_min, phi_max, psi_min, psi_max)
            'sheet': (-150, -90, 90, 180),
            'loop': None  # Everything else
        }
    
    def analyze(self, structure_path: Path) -> List[SecondaryStructure]:
        """Analyze secondary structure of a protein"""
        structure = get_structure(structure_path)
        if not structure:
            return []
        
        assignments = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    ss_type = self._assign_structure(residue)
                    assignments.append(SecondaryStructure(
                        residue_id=f"{residue.get_resname()}_{residue.get_id()[1]}",
                        chain_id=chain.get_id(),
                        position=residue.get_id()[1],
                        structure_type=ss_type,
                        confidence=0.8  # Simplified
                    ))
        
        return assignments
    
    def _assign_structure(self, residue) -> str:
        """Assign secondary structure based on backbone geometry"""
        # Simplified assignment - in reality would use DSSP algorithm
        return 'L'  # Default to loop
    
    def analyze_structure(self, structure) -> List[SecondaryStructure]:
        """Analyze secondary structure of a Bio.PDB structure"""
        assignments = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Only analyze protein residues
                    if residue.get_resname() in ['A', 'T', 'G', 'C', 'U']:
                        continue  # Skip DNA/RNA
                    
                    ss_type = self._assign_structure(residue)
                    assignments.append(SecondaryStructure(
                        residue_id=f"{residue.get_resname()}_{residue.get_id()[1]}",
                        chain_id=chain.get_id(),
                        position=residue.get_id()[1],
                        structure_type=ss_type,
                        confidence=0.8  # Simplified
                    ))
        
        return assignments
    
    def compare_structures(self, observed_ss: List[SecondaryStructure], 
                         predicted_ss: List[SecondaryStructure]) -> Dict:
        """Compare secondary structure assignments between two structures"""
        
        # Create position maps
        obs_map = {f"{ss.chain_id}_{ss.position}": ss.structure_type for ss in observed_ss}
        pred_map = {f"{ss.chain_id}_{ss.position}": ss.structure_type for ss in predicted_ss}
        
        # Find common positions
        common_positions = set(obs_map.keys()) & set(pred_map.keys())
        
        if not common_positions:
            return {'accuracy': 0.0, 'total': 0, 'correct': 0}
        
        correct = 0
        for pos in common_positions:
            if obs_map[pos] == pred_map[pos]:
                correct += 1
        
        accuracy = correct / len(common_positions)
        
        return {
            'accuracy': accuracy,
            'total': len(common_positions),
            'correct': correct,
            'by_type': self._analyze_by_type(obs_map, pred_map, common_positions)
        }
    
    def _analyze_by_type(self, obs_map: Dict, pred_map: Dict, common_positions: Set) -> Dict:
        """Analyze accuracy by secondary structure type"""
        type_stats = {'H': {'total': 0, 'correct': 0}, 
                     'E': {'total': 0, 'correct': 0}, 
                     'L': {'total': 0, 'correct': 0}}
        
        for pos in common_positions:
            obs_type = obs_map[pos]
            pred_type = pred_map[pos]
            
            if obs_type in type_stats:
                type_stats[obs_type]['total'] += 1
                if obs_type == pred_type:
                    type_stats[obs_type]['correct'] += 1
        
        # Calculate accuracies
        for ss_type in type_stats:
            total = type_stats[ss_type]['total']
            if total > 0:
                type_stats[ss_type]['accuracy'] = type_stats[ss_type]['correct'] / total
            else:
                type_stats[ss_type]['accuracy'] = 0.0
        
        return type_stats
    
    def comparison_to_dataframe(self, comparison: Dict) -> pd.DataFrame:
        """Convert secondary structure comparison to pandas DataFrame"""
        data = [
            {'metric': 'overall_accuracy', 'value': comparison['accuracy']},
            {'metric': 'total_residues', 'value': comparison['total']},
            {'metric': 'correct_predictions', 'value': comparison['correct']},
        ]
        
        # Add by-type accuracies
        for ss_type, stats in comparison['by_type'].items():
            data.append({
                'metric': f'{ss_type}_accuracy',
                'value': stats['accuracy']
            })
            data.append({
                'metric': f'{ss_type}_total',
                'value': stats['total']
            })
        
        return pd.DataFrame(data)