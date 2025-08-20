# Create the mutations analyzer file
"""
biostructbenchmark/analysis/mutations.py
Mutation detection and impact analysis
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass

from biostructbenchmark.core.io import get_structure
from biostructbenchmark.core.alignment import ResidueRMSD


@dataclass
class Mutation:
    """Container for mutation information"""
    chain_id: str
    position: int
    wild_type: str  # Original residue
    mutant: str  # Mutated residue
    local_rmsd: float  # RMSD at mutation site
    neighbor_rmsd: float  # Average RMSD of neighboring residues
    impact_score: float  # Calculated impact on structure
    mutation_type: str  # 'conservative', 'non-conservative', 'radical'


class MutationAnalyzer:
    """Detect and analyze mutations between structures"""
    
    def __init__(self):
        """Initialize mutation analyzer"""
        self.aa_properties = self._load_aa_properties()
    
    def _load_aa_properties(self) -> Dict[str, Dict]:
        """Load amino acid properties for classification"""
        return {
            # Hydrophobic
            'ALA': {'hydrophobic': True, 'charge': 0, 'size': 'small'},
            'VAL': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            'ILE': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            'LEU': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            'MET': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            'PHE': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'TRP': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'PRO': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            
            # Polar
            'SER': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'THR': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'CYS': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'TYR': {'hydrophobic': False, 'charge': 0, 'size': 'large'},
            'ASN': {'hydrophobic': False, 'charge': 0, 'size': 'medium'},
            'GLN': {'hydrophobic': False, 'charge': 0, 'size': 'medium'},
            'GLY': {'hydrophobic': False, 'charge': 0, 'size': 'tiny'},
            
            # Charged
            'ASP': {'hydrophobic': False, 'charge': -1, 'size': 'medium'},
            'GLU': {'hydrophobic': False, 'charge': -1, 'size': 'medium'},
            'LYS': {'hydrophobic': False, 'charge': 1, 'size': 'large'},
            'ARG': {'hydrophobic': False, 'charge': 1, 'size': 'large'},
            'HIS': {'hydrophobic': False, 'charge': 1, 'size': 'medium'},
        }
    
    def detect_mutations(self, observed_path: Path, predicted_path: Path) -> List[Mutation]:
        """
        Detect mutations between two structures
        
        Args:
            observed_path: Path to observed/reference structure
            predicted_path: Path to predicted/mutated structure
            
        Returns:
            List of detected mutations
        """
        obs_struct = get_structure(observed_path)
        pred_struct = get_structure(predicted_path)
        
        if not obs_struct or not pred_struct:
            return []
        
        mutations = []
        
        # Compare residues in each chain
        for obs_chain, pred_chain in zip(obs_struct[0], pred_struct[0]):
            if obs_chain.get_id() != pred_chain.get_id():
                continue
            
            obs_residues = {r.get_id()[1]: r for r in obs_chain}
            pred_residues = {r.get_id()[1]: r for r in pred_chain}
            
            # Find common positions
            common_positions = set(obs_residues.keys()) & set(pred_residues.keys())
            
            for position in common_positions:
                obs_res = obs_residues[position]
                pred_res = pred_residues[position]
                
                obs_name = obs_res.get_resname().strip()
                pred_name = pred_res.get_resname().strip()
                
                # Check if mutation occurred
                if obs_name != pred_name and obs_name in self.aa_properties and pred_name in self.aa_properties:
                    mutation = Mutation(
                        chain_id=obs_chain.get_id(),
                        position=position,
                        wild_type=obs_name,
                        mutant=pred_name,
                        local_rmsd=0.0,  # Will be calculated
                        neighbor_rmsd=0.0,  # Will be calculated
                        impact_score=0.0,  # Will be calculated
                        mutation_type=self._classify_mutation(obs_name, pred_name)
                    )
                    mutations.append(mutation)
        
        return mutations
    
    def _classify_mutation(self, wild_type: str, mutant: str) -> str:
        """
        Classify mutation type based on amino acid properties
        
        Args:
            wild_type: Original amino acid
            mutant: Mutated amino acid
            
        Returns:
            Mutation classification
        """
        if wild_type not in self.aa_properties or mutant not in self.aa_properties:
            return 'unknown'
        
        wt_props = self.aa_properties[wild_type]
        mut_props = self.aa_properties[mutant]
        
        # Check property changes
        hydrophobic_change = wt_props['hydrophobic'] != mut_props['hydrophobic']
        charge_change = wt_props['charge'] != mut_props['charge']
        size_change = wt_props['size'] != mut_props['size']
        
        if charge_change:
            return 'radical'  # Charge change is most disruptive
        elif hydrophobic_change:
            return 'non-conservative'  # Hydrophobicity change
        elif size_change:
            return 'non-conservative'  # Size change
        else:
            return 'conservative'  # Similar properties
    
    def identify_mutations(self, observed_structure, predicted_structure) -> List[Mutation]:
        """
        Identify mutations between two Bio.PDB structures
        
        Args:
            observed_structure: Bio.PDB structure (observed/reference)
            predicted_structure: Bio.PDB structure (predicted/mutated)
            
        Returns:
            List of detected mutations
        """
        mutations = []
        
        # Compare residues in each chain
        for obs_model in observed_structure:
            for pred_model in predicted_structure:
                obs_chains = {chain.get_id(): chain for chain in obs_model}
                pred_chains = {chain.get_id(): chain for chain in pred_model}
                
                common_chains = set(obs_chains.keys()) & set(pred_chains.keys())
                
                for chain_id in common_chains:
                    obs_chain = obs_chains[chain_id]
                    pred_chain = pred_chains[chain_id]
                    
                    obs_residues = {r.get_id()[1]: r for r in obs_chain}
                    pred_residues = {r.get_id()[1]: r for r in pred_chain}
                    
                    # Find common positions
                    common_positions = set(obs_residues.keys()) & set(pred_residues.keys())
                    
                    for position in common_positions:
                        obs_res = obs_residues[position]
                        pred_res = pred_residues[position]
                        
                        obs_name = obs_res.get_resname().strip()
                        pred_name = pred_res.get_resname().strip()
                        
                        # Check if mutation occurred (only for protein residues)
                        if (obs_name != pred_name and 
                            obs_name in self.aa_properties and 
                            pred_name in self.aa_properties):
                            
                            mutation = Mutation(
                                chain_id=chain_id,
                                position=position,
                                wild_type=obs_name,
                                mutant=pred_name,
                                local_rmsd=0.0,  # Will be calculated
                                neighbor_rmsd=0.0,  # Will be calculated
                                impact_score=0.0,  # Will be calculated
                                mutation_type=self._classify_mutation(obs_name, pred_name)
                            )
                            mutations.append(mutation)
                break  # Only process first model
        
        return mutations
    
    def analyze_mutation_impact(self, mutations: List[Mutation], 
                               residue_rmsds: List[ResidueRMSD]) -> Dict:
        """
        Analyze the impact of mutations on structural accuracy
        
        Args:
            mutations: List of detected mutations
            residue_rmsds: Per-residue RMSD data from alignment
            
        Returns:
            Dictionary with mutation impact analysis
        """
        if not mutations or not residue_rmsds:
            return {'total_mutations': 0, 'avg_impact': 0.0}
        
        # Create RMSD lookup
        rmsd_map = {}
        for rmsd_data in residue_rmsds:
            key = f"{rmsd_data.chain_id}_{rmsd_data.position}"
            rmsd_map[key] = rmsd_data.rmsd
        
        # Analyze each mutation
        mutation_impacts = []
        for mutation in mutations:
            mut_key = f"{mutation.chain_id}_{mutation.position}"
            
            if mut_key in rmsd_map:
                local_rmsd = rmsd_map[mut_key]
                
                # Calculate neighbor RMSD (positions within ±2)
                neighbor_rmsds = []
                for offset in [-2, -1, 1, 2]:
                    neighbor_key = f"{mutation.chain_id}_{mutation.position + offset}"
                    if neighbor_key in rmsd_map:
                        neighbor_rmsds.append(rmsd_map[neighbor_key])
                
                neighbor_rmsd = np.mean(neighbor_rmsds) if neighbor_rmsds else 0.0
                
                # Calculate impact score based on RMSD and mutation type
                type_multiplier = {
                    'radical': 3.0,
                    'non-conservative': 2.0,
                    'conservative': 1.0
                }.get(mutation.mutation_type, 1.0)
                
                impact_score = local_rmsd * type_multiplier
                
                # Update mutation with calculated values
                mutation.local_rmsd = local_rmsd
                mutation.neighbor_rmsd = neighbor_rmsd
                mutation.impact_score = impact_score
                
                mutation_impacts.append(impact_score)
        
        # Calculate summary statistics
        avg_impact = np.mean(mutation_impacts) if mutation_impacts else 0.0
        max_impact = max(mutation_impacts) if mutation_impacts else 0.0
        
        # Categorize mutations by impact
        high_impact = [m for m in mutations if m.impact_score > 5.0]
        medium_impact = [m for m in mutations if 2.0 < m.impact_score <= 5.0]
        low_impact = [m for m in mutations if m.impact_score <= 2.0]
        
        return {
            'total_mutations': len(mutations),
            'avg_impact': avg_impact,
            'max_impact': max_impact,
            'high_impact_count': len(high_impact),
            'medium_impact_count': len(medium_impact),
            'low_impact_count': len(low_impact),
            'by_type': self._analyze_by_mutation_type(mutations)
        }
    
    def _analyze_by_mutation_type(self, mutations: List[Mutation]) -> Dict:
        """Analyze mutations by type"""
        type_counts = {'radical': 0, 'non-conservative': 0, 'conservative': 0}
        type_impacts = {'radical': [], 'non-conservative': [], 'conservative': []}
        
        for mutation in mutations:
            mut_type = mutation.mutation_type
            if mut_type in type_counts:
                type_counts[mut_type] += 1
                type_impacts[mut_type].append(mutation.impact_score)
        
        # Calculate average impacts
        type_analysis = {}
        for mut_type in type_counts:
            impacts = type_impacts[mut_type]
            type_analysis[mut_type] = {
                'count': type_counts[mut_type],
                'avg_impact': np.mean(impacts) if impacts else 0.0,
                'max_impact': max(impacts) if impacts else 0.0
            }
        
        return type_analysis