# Create the file with the consensus analyzer code
"""
biostructbenchmark/analysis/consensus.py
Consensus error mapping across multiple structure comparisons
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from biostructbenchmark.core.alignment import ResidueRMSD


@dataclass
class ConsensusError:
    """Container for consensus error at a position"""
    position: str  # e.g., "A_10_G" for chain A position 10 with base G (DNA) or "A_10" (protein)
    mean_rmsd: float
    std_rmsd: float
    frequency: int  # Number of structures with this position
    confidence: float  # Statistical confidence of the error
    residue_types: List[str]  # Residue types seen at this position
    dna_interface_residues: List[str] = field(default_factory=list)  # DNA-binding protein residues (e.g., ["ARG", "LYS"])


class ConsensusAnalyzer:
    """Identify consistently mispredicted regions across structures"""
    
    def __init__(self, rmsd_threshold: float = 3.0):
        """
        Initialize consensus analyzer
        
        Args:
            rmsd_threshold: RMSD threshold for considering a position as mispredicted
        """
        self.rmsd_threshold = rmsd_threshold
        self.position_errors = defaultdict(list)
    
    def add_structure_comparison(self, residue_rmsds: List[ResidueRMSD], 
                               structure_id: str) -> None:
        """
        Add RMSD data from one structure comparison
        
        Args:
            residue_rmsds: Per-residue RMSD data from alignment
            structure_id: Identifier for this structure pair
        """
        for rmsd_data in residue_rmsds:
            # Include DNA base identity in position format: "A_10_G" for DNA base G at position 10
            if rmsd_data.is_dna:
                position_key = f"{rmsd_data.chain_id}_{rmsd_data.position}_{rmsd_data.residue_type}"
            else:
                position_key = f"{rmsd_data.chain_id}_{rmsd_data.position}"
            
            self.position_errors[position_key].append({
                'rmsd': rmsd_data.rmsd,
                'residue_type': rmsd_data.residue_type,
                'structure_id': structure_id,
                'molecule_type': rmsd_data.molecule_type,
                'rmsd_data': rmsd_data  # Store full data for phosphate distance calculations
            })
    
    def calculate_consensus(self, min_structures: int = 3) -> List[ConsensusError]:
        """
        Calculate consensus errors across all added structures
        
        Args:
            min_structures: Minimum number of structures required for consensus
            
        Returns:
            List of consensus errors sorted by mean RMSD
        """
        consensus_errors = []
        
        for position, errors in self.position_errors.items():
            if len(errors) < min_structures:
                continue
            
            rmsds = [e['rmsd'] for e in errors]
            residue_types = list(set(e['residue_type'] for e in errors))
            
            mean_rmsd = np.mean(rmsds)
            std_rmsd = np.std(rmsds)
            
            # Calculate DNA-specific metrics for interface analysis
            dna_interface_residues = []
            electrostatic_conf = 1.0
            
            # Check if this involves protein residues near DNA (for interface analysis)
            protein_errors = [e for e in errors if e.get('molecule_type') == 'protein']
            if protein_errors:
                # Collect interface residue types
                interface_types = [e['residue_type'] for e in protein_errors]
                dna_interface_residues = list(set(interface_types))
                
                # Calculate electrostatic confidence for protein-DNA interfaces
                phosphate_distances = []
                for error in protein_errors:
                    rmsd_data = error.get('rmsd_data')
                    if rmsd_data and hasattr(rmsd_data, 'min_phosphate_distance'):
                        try:
                            dist = rmsd_data.min_phosphate_distance()
                            if dist != float('inf'):
                                phosphate_distances.append(dist)
                        except:
                            continue
                
                if phosphate_distances:
                    mean_phosphate_dist = np.mean(phosphate_distances)
                    # Critical threshold: electrostatic interaction disrupted if > 3.5 Å
                    electrostatic_conf = 1.0 if mean_phosphate_dist < 3.5 else 0.0
            
            # Calculate enhanced confidence with electrostatic component
            if mean_rmsd > 0:
                base_confidence = (1 - std_rmsd/mean_rmsd) * min(1.0, mean_rmsd/self.rmsd_threshold)
                confidence = base_confidence * electrostatic_conf
            else:
                confidence = 0.0
            
            consensus_error = ConsensusError(
                position=position,
                mean_rmsd=mean_rmsd,
                std_rmsd=std_rmsd,
                frequency=len(errors),
                confidence=max(0, min(1, confidence)),
                residue_types=residue_types,
                dna_interface_residues=dna_interface_residues
            )
            
            consensus_errors.append(consensus_error)
        
        # Sort by mean RMSD (worst first)
        consensus_errors.sort(key=lambda x: x.mean_rmsd, reverse=True)
        
        return consensus_errors
    
    def identify_problem_regions(self, consensus_errors: List[ConsensusError],
                                rmsd_cutoff: float = None) -> Dict[str, List[ConsensusError]]:
        """
        Identify problematic regions based on consensus errors
        
        Args:
            consensus_errors: List of consensus errors
            rmsd_cutoff: RMSD threshold (uses self.rmsd_threshold if None)
            
        Returns:
            Dictionary categorizing errors by severity
        """
        if rmsd_cutoff is None:
            rmsd_cutoff = self.rmsd_threshold
        
        regions = {
            'critical': [],  # Very high RMSD, high confidence
            'problematic': [],  # High RMSD, medium confidence
            'moderate': [],  # Moderate RMSD
            'acceptable': []  # Low RMSD
        }
        
        for error in consensus_errors:
            if error.mean_rmsd > rmsd_cutoff * 1.5 and error.confidence > 0.7:
                regions['critical'].append(error)
            elif error.mean_rmsd > rmsd_cutoff and error.confidence > 0.5:
                regions['problematic'].append(error)
            elif error.mean_rmsd > rmsd_cutoff * 0.5:
                regions['moderate'].append(error)
            else:
                regions['acceptable'].append(error)
        
        return regions
    
    def export_consensus_map(self, consensus_errors: List[ConsensusError],
                            output_path: Path) -> None:
        """
        Export consensus error map to CSV
        
        Args:
            consensus_errors: List of consensus errors
            output_path: Path for output CSV file
        """
        data = []
        for error in consensus_errors:
            data.append({
                'position': error.position,
                'mean_rmsd': f"{error.mean_rmsd:.3f}",
                'std_rmsd': f"{error.std_rmsd:.3f}",
                'frequency': error.frequency,
                'confidence': f"{error.confidence:.3f}",
                'residue_types': ','.join(error.residue_types)
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        print(f"Consensus error map saved to: {output_path}")
    
    def generate_summary_report(self, consensus_errors: List[ConsensusError]) -> Dict:
        """
        Generate summary statistics of consensus errors
        
        Args:
            consensus_errors: List of consensus errors
            
        Returns:
            Dictionary with summary statistics
        """
        if not consensus_errors:
            return {'total_positions': 0}
        
        rmsds = [e.mean_rmsd for e in consensus_errors]
        high_error_positions = [e for e in consensus_errors if e.mean_rmsd > self.rmsd_threshold]
        
        return {
            'total_positions': len(consensus_errors),
            'mean_rmsd_overall': np.mean(rmsds),
            'std_rmsd_overall': np.std(rmsds),
            'max_rmsd': max(rmsds),
            'min_rmsd': min(rmsds),
            'high_error_count': len(high_error_positions),
            'high_error_percentage': 100 * len(high_error_positions) / len(consensus_errors),
            'worst_positions': [
                {'position': e.position, 'rmsd': e.mean_rmsd}
                for e in consensus_errors[:10]
            ]
        }
    
    def identify_consensus_errors(self, rmsd_datasets: List[List]) -> List[ConsensusError]:
        """
        Identify consensus errors from multiple RMSD datasets
        
        Args:
            rmsd_datasets: List of ResidueRMSD lists from multiple structure comparisons
            
        Returns:
            List of consensus errors
        """
        # Reset analyzer state
        self.position_errors = defaultdict(list)
        
        # Add each dataset
        for i, dataset in enumerate(rmsd_datasets):
            self.add_structure_comparison(dataset, f"structure_{i}")
        
        # Calculate and return consensus
        return self.calculate_consensus(min_structures=1)  # Lower threshold for testing
    
    def to_dataframe(self, consensus_errors: List[ConsensusError]) -> pd.DataFrame:
        """Convert consensus errors to pandas DataFrame"""
        data = []
        for error in consensus_errors:
            data.append({
                'position': error.position,
                'mean_rmsd': error.mean_rmsd,
                'std_rmsd': error.std_rmsd,
                'frequency': error.frequency,
                'confidence': error.confidence,
                'residue_types': ','.join(error.residue_types),
                'dna_interface_residues': ','.join(error.dna_interface_residues) if error.dna_interface_residues else ''
            })
        return pd.DataFrame(data)
    
    def analyze_dna_specific_errors(self, consensus_errors: List[ConsensusError]) -> Dict[str, any]:
        """
        Analyze DNA-specific systemic errors in protein-DNA complex modeling
        
        Args:
            consensus_errors: List of consensus errors
            
        Returns:
            Dictionary with DNA-specific error analysis
        """
        dna_errors = [e for e in consensus_errors if '_' in e.position and len(e.position.split('_')) > 2]
        protein_interface_errors = [e for e in consensus_errors if e.dna_interface_residues]
        
        # Base-specific error patterns (e.g., "G" vs "A" contexts)
        base_specific_errors = {}
        for error in dna_errors:
            try:
                base = error.position.split('_')[-1]  # Extract base from "A_10_G" format
                if base not in base_specific_errors:
                    base_specific_errors[base] = []
                base_specific_errors[base].append(error.mean_rmsd)
            except:
                continue
        
        # Calculate base-specific statistics
        base_stats = {}
        for base, rmsds in base_specific_errors.items():
            if rmsds:
                base_stats[base] = {
                    'mean_rmsd': np.mean(rmsds),
                    'count': len(rmsds),
                    'max_rmsd': max(rmsds)
                }
        
        # Interface residue analysis
        interface_residue_counts = {}
        total_interface_positions = 0
        arg_lys_positions = 0
        
        for error in protein_interface_errors:
            total_interface_positions += 1
            for residue in error.dna_interface_residues:
                interface_residue_counts[residue] = interface_residue_counts.get(residue, 0) + 1
                if residue in ['ARG', 'LYS']:
                    arg_lys_positions += 1
        
        # Calculate critical interface metrics
        arg_lys_percentage = (arg_lys_positions / total_interface_positions * 100) if total_interface_positions > 0 else 0
        
        # Identify critical phosphate distance errors (electrostatic disruption)
        critical_electrostatic_errors = [e for e in consensus_errors if e.confidence == 0.0 and e.dna_interface_residues]
        
        return {
            'base_specific_patterns': base_stats,
            'interface_residue_distribution': interface_residue_counts,
            'arg_lys_interface_percentage': arg_lys_percentage,
            'total_dna_positions': len(dna_errors),
            'total_interface_positions': total_interface_positions,
            'critical_electrostatic_errors': len(critical_electrostatic_errors),
            'systemic_issues': {
                'insufficient_arg_lys': arg_lys_percentage < 30.0,
                'base_bias_detected': len(base_stats) > 0 and max(base_stats.values(), key=lambda x: x['mean_rmsd'])['mean_rmsd'] > 2.5 * min(base_stats.values(), key=lambda x: x['mean_rmsd'])['mean_rmsd'] if len(base_stats) > 1 else False,
                'electrostatic_disruption': len(critical_electrostatic_errors) > 0
            }
        }