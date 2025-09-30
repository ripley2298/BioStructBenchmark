#!/usr/bin/env python3
"""
Systematic Modeling Error Analysis for Protein Mutants via Residue-Level Geometric PCA

A pipeline for identifying systematic modeling errors in protein mutants through
residue-level geometric PCA (NOT global RMSD). Uses backbone dihedrals and
x3DNA-DSSR DNA geometry features for comprehensive structural analysis.

Author: Senior Structural Bioinformatician
Dependencies: Biopython, x3DNA-DSSR, scikit-learn, numpy, pandas
"""

import numpy as np
import pandas as pd
import json
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# Biopython imports
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_dihedral


@dataclass
class MutantAnalysisResult:
    """Container for mutant analysis results"""
    mutant_id: str
    global_rmsd: float  # Supplementary context only
    backbone_features: np.ndarray
    dna_features: np.ndarray
    combined_features: np.ndarray
    pca_scores: np.ndarray
    error_hotspots: List[str]  # Residue IDs with high PC loadings
    systematic_errors: Dict[str, float]  # PC loadings by feature type


@dataclass 
class ResidueFeatures:
    """Container for residue-level structural features"""
    residue_id: str
    chain_id: str
    position: int
    phi: Optional[float]
    psi: Optional[float]
    dna_distance: Optional[float]  # Distance to nearest DNA if present
    interface_score: Optional[float]  # DNA interface involvement


class MutantPCAAnalyzer:
    """
    Systematic modeling error analysis for protein mutants using residue-level geometric PCA.
    
    This analyzer identifies systematic modeling errors through:
    1. Backbone dihedral angle analysis (phi/psi)
    2. DNA interface geometry features (x3DNA-DSSR)
    3. Principal Component Analysis for error pattern detection
    4. Residue-level error hotspot identification
    """
    
    def __init__(self, wildtype_pdb: str, mutant_pdb_dict: Dict[str, str]):
        """
        Initialize the mutant PCA analyzer.
        
        Args:
            wildtype_pdb: Path to wild-type reference structure (PDB)
            mutant_pdb_dict: Dictionary mapping mutant IDs to PDB file paths
                            e.g., {"mutant_A": "mutant_A.pdb", "mutant_B": "mutant_B.pdb"}
        """
        self.wildtype_pdb = Path(wildtype_pdb)
        self.mutant_pdb_dict = {k: Path(v) for k, v in mutant_pdb_dict.items()}
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        
        # Analysis parameters
        self.n_components = 3  # Principal components to analyze
        self.error_threshold = 2.0  # Standard deviations for error hotspots
        
        # Storage for analysis results
        self.wildtype_structure = None
        self.mutant_structures = {}
        self.common_residues = []
        self.global_rmsds = {}
        self.feature_matrix = None
        self.pca_model = None
        self.scaler = None
        
        print(f"🧬 Initializing Mutant PCA Analyzer")
        print(f"   Wild-type: {self.wildtype_pdb}")
        print(f"   Mutants: {list(self.mutant_pdb_dict.keys())}")

    def _get_parser_for_file(self, file_path: Path):
        """Get appropriate parser based on file extension"""
        if file_path.suffix.lower() in ['.cif', '.mmcif']:
            return self.cif_parser
        else:
            return self.pdb_parser

    def load_and_align_structures(self) -> None:
        """
        Step 1: Load structures and align mutants to wild-type reference.
        
        - Loads wild-type and all mutant structures
        - Identifies common CA atoms across all structures
        - Performs structural alignment of mutants to wild-type
        - Calculates global RMSD for supplementary context (NOT for analysis)
        """
        print(f"\n📂 Step 1: Loading and aligning structures...")
        
        # Load wild-type structure
        wt_parser = self._get_parser_for_file(self.wildtype_pdb)
        self.wildtype_structure = wt_parser.get_structure("WT", self.wildtype_pdb)
        wt_residues = self._get_protein_residues(self.wildtype_structure)
        wt_residue_ids = {(res.get_parent().id, res.id[1]): res for res in wt_residues}
        
        print(f"   ✓ Loaded wild-type: {len(wt_residues)} protein residues")
        
        # Load mutant structures and find common residues
        self.mutant_structures = {}
        common_residue_ids = set(wt_residue_ids.keys())
        
        for mutant_id, pdb_path in self.mutant_pdb_dict.items():
            mutant_parser = self._get_parser_for_file(pdb_path)
            structure = mutant_parser.get_structure(mutant_id, pdb_path)
            self.mutant_structures[mutant_id] = structure
            
            mutant_residues = self._get_protein_residues(structure)
            mutant_residue_ids = {(res.get_parent().id, res.id[1]): res for res in mutant_residues}
            
            # Update common residues (intersection)
            common_residue_ids &= set(mutant_residue_ids.keys())
            print(f"   ✓ Loaded {mutant_id}: {len(mutant_residues)} protein residues")
        
        # Store common residues for alignment
        self.common_residues = sorted(common_residue_ids)
        print(f"   ✓ Common residues for alignment: {len(self.common_residues)}")
        
        # Perform structural alignment and calculate global RMSD
        self._align_structures_and_calculate_rmsd()
        
        # Report global RMSD (supplementary context only)
        rmsd_list = [f"{mut_id}={rmsd:.2f}Å" for mut_id, rmsd in self.global_rmsds.items()]
        print(f"\n📊 RMSD: {', '.join(rmsd_list)}")
        print(f"   (Global RMSD reported for context - NOT used for error analysis)")

    def _get_protein_residues(self, structure: Structure) -> List[Residue]:
        """Extract protein residues (standard amino acids) from structure"""
        protein_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residue (not heteroatom)
                        if residue.get_resname() in self._standard_amino_acids():
                            protein_residues.append(residue)
        return protein_residues
    
    def _standard_amino_acids(self) -> set:
        """Standard 20 amino acid codes"""
        return {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
    
    def _align_structures_and_calculate_rmsd(self) -> None:
        """Align mutant structures to wild-type and calculate global RMSD"""
        # Get wild-type CA atoms for common residues
        wt_ca_atoms = []
        for chain_id, res_id in self.common_residues:
            try:
                residue = self.wildtype_structure[0][chain_id][res_id]
                if 'CA' in residue:
                    wt_ca_atoms.append(residue['CA'])
            except KeyError:
                continue
        
        print(f"   → Aligning {len(wt_ca_atoms)} CA atoms...")
        
        # Align each mutant to wild-type
        for mutant_id, structure in self.mutant_structures.items():
            # Get corresponding CA atoms from mutant
            mutant_ca_atoms = []
            for chain_id, res_id in self.common_residues:
                try:
                    residue = structure[0][chain_id][res_id]
                    if 'CA' in residue:
                        mutant_ca_atoms.append(residue['CA'])
                except KeyError:
                    continue
            
            # Ensure same number of atoms
            min_atoms = min(len(wt_ca_atoms), len(mutant_ca_atoms))
            if min_atoms < 10:  # Minimum atoms for meaningful alignment
                warnings.warn(f"Too few common CA atoms ({min_atoms}) for {mutant_id}")
                continue
            
            # Perform superimposition
            superimposer = Superimposer()
            superimposer.set_atoms(wt_ca_atoms[:min_atoms], mutant_ca_atoms[:min_atoms])
            
            # Apply transformation to entire mutant structure
            superimposer.apply(structure.get_atoms())
            
            # Calculate and store global RMSD
            self.global_rmsds[mutant_id] = superimposer.rms
            
        print(f"   ✓ Structural alignment complete")

    def extract_residue_features(self) -> None:
        """
        Step 2: Extract residue-level features for PCA analysis.
        
        Extracts:
        - Backbone dihedral angles (phi/psi) for all residues
        - DNA interface features using x3DNA-DSSR (if DNA present)
        - Residue-level structural descriptors
        """
        print(f"\n🔬 Step 2: Extracting residue-level features...")
        
        all_features = []
        
        # Process wild-type first (reference)
        wt_features = self._extract_structure_features(
            self.wildtype_structure, "WT", is_reference=True
        )
        all_features.extend(wt_features)
        
        # Process all mutants
        for mutant_id, structure in self.mutant_structures.items():
            mutant_features = self._extract_structure_features(
                structure, mutant_id, is_reference=False
            )
            all_features.extend(mutant_features)
        
        # Convert to feature matrix
        self.feature_matrix = self._build_feature_matrix(all_features)
        print(f"   ✓ Feature matrix: {self.feature_matrix.shape}")

    def _extract_structure_features(self, structure: Structure, struct_id: str, 
                                  is_reference: bool = False) -> List[Dict]:
        """Extract comprehensive residue-level features from structure"""
        features = []
        
        # Check for DNA chains
        has_dna = self._check_dna_presence(structure)
        dna_features = {}
        
        if has_dna:
            print(f"   → Extracting DNA features for {struct_id}...")
            dna_features = self._extract_dssr_features(structure, struct_id)
        
        # Extract backbone dihedral features
        protein_residues = self._get_protein_residues(structure)
        print(f"   → Extracting backbone dihedrals for {struct_id}: {len(protein_residues)} residues")
        
        for residue in protein_residues:
            chain_id = residue.get_parent().id
            res_id = residue.id[1]
            residue_key = f"{chain_id}_{res_id}"
            
            # Skip if not in common residues
            if (chain_id, res_id) not in self.common_residues:
                continue
            
            feature_dict = {
                'structure_id': struct_id,
                'residue_id': residue_key,
                'chain_id': chain_id,
                'position': res_id,
                'resname': residue.get_resname(),
                'is_reference': is_reference
            }
            
            # Extract backbone dihedrals
            phi, psi = self._calculate_backbone_dihedrals(residue)
            feature_dict['phi'] = phi
            feature_dict['psi'] = psi
            
            # Add DNA interface features if available
            if has_dna and residue_key in dna_features:
                feature_dict.update(dna_features[residue_key])
            else:
                # Default DNA features for non-interface residues
                feature_dict.update({
                    'minor_groove_distance': np.nan,
                    'major_groove_distance': np.nan,
                    'dna_interface_score': 0.0,
                    'base_pair_rise': np.nan,
                    'helical_twist': np.nan
                })
            
            features.append(feature_dict)
        
        return features

    def _check_dna_presence(self, structure: Structure) -> bool:
        """Check if structure contains DNA chains"""
        dna_chains = []
        for model in structure:
            for chain in model:
                # Check for DNA residues (A, T, G, C, DA, DT, DG, DC)
                dna_residues = {'A', 'T', 'G', 'C', 'DA', 'DT', 'DG', 'DC'}
                for residue in chain:
                    if residue.get_resname().strip() in dna_residues:
                        dna_chains.append(chain.id)
                        break
        
        has_dna = len(dna_chains) > 0
        if has_dna:
            print(f"   ✓ DNA detected in chains: {', '.join(set(dna_chains))}")
        
        return has_dna

    def _extract_dssr_features(self, structure: Structure, struct_id: str) -> Dict[str, Dict]:
        """Extract DNA geometry features using x3DNA-DSSR"""
        dna_features = {}
        
        try:
            # Save temporary PDB file for DSSR
            temp_pdb = f"temp_{struct_id}.pdb"
            io = PDBIO()
            io.set_structure(structure)
            io.save(temp_pdb)
            
            # Run x3DNA-DSSR
            dssr_output = f"dssr_{struct_id}.json"
            cmd = f"x3dna-dssr --input={temp_pdb} --output={dssr_output} --json"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and Path(dssr_output).exists():
                # Parse DSSR output
                with open(dssr_output, 'r') as f:
                    dssr_data = json.load(f)
                
                dna_features = self._parse_dssr_output(dssr_data, structure)
                print(f"     ✓ DSSR analysis complete: {len(dna_features)} interface residues")
                
                # Cleanup temporary files
                Path(temp_pdb).unlink(missing_ok=True)
                Path(dssr_output).unlink(missing_ok=True)
            else:
                print(f"     ⚠ DSSR failed for {struct_id}: {result.stderr}")
                
        except Exception as e:
            print(f"     ⚠ DSSR extraction failed for {struct_id}: {e}")
        
        return dna_features

    def _parse_dssr_output(self, dssr_data: Dict, structure: Structure) -> Dict[str, Dict]:
        """Parse x3DNA-DSSR output and map to protein residues"""
        interface_features = {}
        
        # Extract base pair parameters
        bp_params = {}
        if 'basePairs' in dssr_data:
            for bp in dssr_data['basePairs']:
                if 'bp' in bp and 'rise' in bp:
                    bp_id = bp.get('index', len(bp_params))
                    bp_params[bp_id] = {
                        'rise': bp['rise'],
                        'twist': bp.get('twist', np.nan),
                        'major_groove': bp.get('majorGroove', np.nan),
                        'minor_groove': bp.get('minorGroove', np.nan)
                    }
        
        # Extract DNA-protein interface information
        protein_residues = self._get_protein_residues(structure)
        
        for residue in protein_residues:
            chain_id = residue.get_parent().id
            res_id = residue.id[1]
            residue_key = f"{chain_id}_{res_id}"
            
            # Calculate distance to nearest DNA
            min_dna_distance = self._calculate_dna_distance(residue, structure)
            
            # Determine interface involvement (< 5Å cutoff)
            is_interface = min_dna_distance < 5.0 if min_dna_distance is not None else False
            
            interface_features[residue_key] = {
                'minor_groove_distance': min_dna_distance,
                'major_groove_distance': min_dna_distance,  # Simplified
                'dna_interface_score': 1.0 if is_interface else 0.0,
                'base_pair_rise': np.mean([bp['rise'] for bp in bp_params.values()]) if bp_params else np.nan,
                'helical_twist': np.mean([bp['twist'] for bp in bp_params.values() if not np.isnan(bp['twist'])]) if bp_params else np.nan
            }
        
        return interface_features

    def _calculate_dna_distance(self, protein_residue: Residue, structure: Structure) -> Optional[float]:
        """Calculate minimum distance from protein residue to DNA"""
        min_distance = None
        
        if 'CA' not in protein_residue:
            return None
        
        protein_ca = protein_residue['CA']
        
        # Find DNA residues and calculate distances
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in {'A', 'T', 'G', 'C', 'DA', 'DT', 'DG', 'DC'}:
                        # Get phosphate atom for distance calculation
                        if 'P' in residue:
                            distance = protein_ca - residue['P']
                            if min_distance is None or distance < min_distance:
                                min_distance = distance
        
        return min_distance

    def _calculate_backbone_dihedrals(self, residue: Residue) -> Tuple[Optional[float], Optional[float]]:
        """Calculate phi and psi dihedral angles for a residue"""
        try:
            # Get previous and next residues for dihedral calculation
            chain = residue.get_parent()
            residues = list(chain)
            current_idx = residues.index(residue)
            
            phi = psi = None
            
            # Calculate phi (C-1, N, CA, C)
            if current_idx > 0:
                prev_residue = residues[current_idx - 1]
                if all(atom in prev_residue for atom in ['C']) and \
                   all(atom in residue for atom in ['N', 'CA', 'C']):
                    phi = calc_dihedral(
                        prev_residue['C'].get_vector(),
                        residue['N'].get_vector(),
                        residue['CA'].get_vector(),
                        residue['C'].get_vector()
                    )
                    phi = np.degrees(phi)
            
            # Calculate psi (N, CA, C, N+1)
            if current_idx < len(residues) - 1:
                next_residue = residues[current_idx + 1]
                if all(atom in residue for atom in ['N', 'CA', 'C']) and \
                   all(atom in next_residue for atom in ['N']):
                    psi = calc_dihedral(
                        residue['N'].get_vector(),
                        residue['CA'].get_vector(),
                        residue['C'].get_vector(),
                        next_residue['N'].get_vector()
                    )
                    psi = np.degrees(psi)
            
        except (KeyError, ValueError, IndexError):
            phi = psi = None
        
        return phi, psi

    def _build_feature_matrix(self, all_features: List[Dict]) -> np.ndarray:
        """Build numerical feature matrix for PCA"""
        df = pd.DataFrame(all_features)
        
        # Select numerical features for PCA
        feature_columns = [
            'phi', 'psi', 
            'minor_groove_distance', 'major_groove_distance',
            'dna_interface_score', 'base_pair_rise', 'helical_twist'
        ]
        
        # Handle missing values
        feature_data = df[feature_columns].fillna(0.0)
        
        # Store metadata for later use
        self.feature_metadata = df[['structure_id', 'residue_id', 'chain_id', 'position', 'is_reference']].copy()
        
        return feature_data.values

    def perform_pca_analysis(self) -> None:
        """
        Step 3: Perform Principal Component Analysis on residue features.
        
        - Scales features using StandardScaler
        - Applies PCA to identify systematic error patterns
        - Identifies error hotspots from PC loadings
        """
        print(f"\n📈 Step 3: Performing PCA analysis...")
        
        # Scale features (critical for PCA with mixed units)
        self.scaler = StandardScaler()
        scaled_features = self.scaler.fit_transform(self.feature_matrix)
        
        # Apply PCA
        self.pca_model = PCA(n_components=self.n_components)
        pca_scores = self.pca_model.fit_transform(scaled_features)
        
        # Store PCA results
        self.pca_scores = pca_scores
        
        # Print PCA summary
        explained_variance = self.pca_model.explained_variance_ratio_
        print(f"   ✓ PCA complete:")
        for i, var in enumerate(explained_variance):
            print(f"     PC{i+1}: {var:.1%} variance explained")
        print(f"     Total: {sum(explained_variance):.1%} variance explained")
        
        # Analyze feature contributions
        self._analyze_feature_contributions()
        
        # Identify error hotspots
        self._identify_error_hotspots()

    def _analyze_feature_contributions(self) -> None:
        """Analyze which features contribute most to each principal component"""
        feature_names = [
            'phi', 'psi', 
            'minor_groove_distance', 'major_groove_distance',
            'dna_interface_score', 'base_pair_rise', 'helical_twist'
        ]
        
        print(f"\n   🔍 Feature contributions to principal components:")
        
        for pc in range(self.n_components):
            loadings = self.pca_model.components_[pc]
            
            # Sort features by absolute loading
            feature_importance = list(zip(feature_names, loadings))
            feature_importance.sort(key=lambda x: abs(x[1]), reverse=True)
            
            print(f"\n     PC{pc+1} (explains {self.pca_model.explained_variance_ratio_[pc]:.1%}):")
            for feat_name, loading in feature_importance[:3]:  # Top 3 features
                direction = "↑" if loading > 0 else "↓"
                print(f"       {direction} {feat_name}: {loading:.3f}")

    def _identify_error_hotspots(self) -> None:
        """Identify residue-level error hotspots from PCA scores"""
        print(f"\n   🎯 Identifying error hotspots...")
        
        # Calculate composite error score from first two PCs
        composite_scores = np.sqrt(self.pca_scores[:, 0]**2 + self.pca_scores[:, 1]**2)
        
        # Identify outliers (> 2 standard deviations)
        threshold = np.mean(composite_scores) + self.error_threshold * np.std(composite_scores)
        hotspot_indices = np.where(composite_scores > threshold)[0]
        
        # Map back to residues
        hotspot_residues = []
        for idx in hotspot_indices:
            residue_info = self.feature_metadata.iloc[idx]
            if not residue_info['is_reference']:  # Only consider mutant residues
                hotspot_residues.append({
                    'structure_id': residue_info['structure_id'],
                    'residue_id': residue_info['residue_id'],
                    'error_score': composite_scores[idx],
                    'pc1_score': self.pca_scores[idx, 0],
                    'pc2_score': self.pca_scores[idx, 1]
                })
        
        # Sort by error score
        hotspot_residues.sort(key=lambda x: x['error_score'], reverse=True)
        
        print(f"     ✓ Identified {len(hotspot_residues)} error hotspots")
        
        # Report top hotspots
        if hotspot_residues:
            print(f"     Top error hotspots:")
            for i, hotspot in enumerate(hotspot_residues[:5]):  # Top 5
                print(f"       {i+1}. {hotspot['structure_id']}:{hotspot['residue_id']} "
                      f"(score: {hotspot['error_score']:.2f})")
        
        self.error_hotspots = hotspot_residues

    def generate_analysis_report(self, output_dir: str = "mutant_pca_analysis") -> None:
        """
        Generate comprehensive analysis report with visualizations.
        
        Args:
            output_dir: Directory to save analysis outputs
        """
        print(f"\n📊 Generating analysis report...")
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Save PCA results
        self._save_pca_results(output_path)
        
        # Save error hotspots
        self._save_error_hotspots(output_path)
        
        # Generate visualizations
        self._generate_visualizations(output_path)
        
        print(f"   ✓ Analysis complete! Results saved to: {output_path.absolute()}")

    def _save_pca_results(self, output_path: Path) -> None:
        """Save PCA results to CSV files"""
        # PCA scores with metadata
        pca_df = self.feature_metadata.copy()
        for i in range(self.n_components):
            pca_df[f'PC{i+1}'] = self.pca_scores[:, i]
        
        pca_df.to_csv(output_path / "pca_scores.csv", index=False)
        
        # Feature loadings
        feature_names = [
            'phi', 'psi', 
            'minor_groove_distance', 'major_groove_distance',
            'dna_interface_score', 'base_pair_rise', 'helical_twist'
        ]
        
        loadings_df = pd.DataFrame(
            self.pca_model.components_.T,
            columns=[f'PC{i+1}' for i in range(self.n_components)],
            index=feature_names
        )
        loadings_df.to_csv(output_path / "feature_loadings.csv")
        
        # Explained variance
        variance_df = pd.DataFrame({
            'Principal_Component': [f'PC{i+1}' for i in range(self.n_components)],
            'Explained_Variance_Ratio': self.pca_model.explained_variance_ratio_,
            'Cumulative_Variance': np.cumsum(self.pca_model.explained_variance_ratio_)
        })
        variance_df.to_csv(output_path / "explained_variance.csv", index=False)

    def _save_error_hotspots(self, output_path: Path) -> None:
        """Save error hotspots analysis"""
        if self.error_hotspots:
            hotspots_df = pd.DataFrame(self.error_hotspots)
            hotspots_df.to_csv(output_path / "error_hotspots.csv", index=False)
        
        # Global RMSD summary
        rmsd_df = pd.DataFrame([
            {'mutant_id': mut_id, 'global_rmsd_CA': rmsd}
            for mut_id, rmsd in self.global_rmsds.items()
        ])
        rmsd_df.to_csv(output_path / "global_rmsd_summary.csv", index=False)

    def _generate_visualizations(self, output_path: Path) -> None:
        """Generate PCA visualizations"""
        try:
            import matplotlib.pyplot as plt
            try:
                import seaborn as sns
                HAS_SEABORN = True
                # Set style
                plt.style.use('default')
                sns.set_palette("husl")
            except ImportError:
                HAS_SEABORN = False
                plt.style.use('default')
            
            # PCA biplot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # PC1 vs PC2 scatter plot
            mutant_mask = ~self.feature_metadata['is_reference']
            wt_mask = self.feature_metadata['is_reference']
            
            # Plot wild-type points
            ax1.scatter(self.pca_scores[wt_mask, 0], self.pca_scores[wt_mask, 1], 
                       c='gray', alpha=0.6, s=20, label='Wild-type')
            
            # Plot mutant points by structure
            for mutant_id in self.mutant_pdb_dict.keys():
                mutant_specific_mask = self.feature_metadata['structure_id'] == mutant_id
                ax1.scatter(self.pca_scores[mutant_specific_mask, 0], 
                           self.pca_scores[mutant_specific_mask, 1],
                           s=30, alpha=0.7, label=mutant_id)
            
            ax1.set_xlabel(f'PC1 ({self.pca_model.explained_variance_ratio_[0]:.1%})')
            ax1.set_ylabel(f'PC2 ({self.pca_model.explained_variance_ratio_[1]:.1%})')
            ax1.set_title('PCA: Residue-Level Geometric Features')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Explained variance plot
            pcs = range(1, self.n_components + 1)
            ax2.bar(pcs, self.pca_model.explained_variance_ratio_, alpha=0.7)
            ax2.set_xlabel('Principal Component')
            ax2.set_ylabel('Explained Variance Ratio')
            ax2.set_title('PCA Explained Variance')
            ax2.set_xticks(pcs)
            
            plt.tight_layout()
            plt.savefig(output_path / "pca_analysis.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            # Feature importance heatmap
            feature_names = [
                'phi', 'psi', 'minor_groove_dist', 'major_groove_dist',
                'dna_interface', 'bp_rise', 'helical_twist'
            ]
            
            fig, ax = plt.subplots(figsize=(8, 6))
            if HAS_SEABORN:
                sns.heatmap(self.pca_model.components_, 
                           xticklabels=feature_names,
                           yticklabels=[f'PC{i+1}' for i in range(self.n_components)],
                           annot=True, cmap='RdBu_r', center=0, ax=ax)
            else:
                # Manual heatmap without seaborn
                im = ax.imshow(self.pca_model.components_, aspect='auto', cmap='RdBu_r')
                ax.set_xticks(range(len(feature_names)))
                ax.set_xticklabels(feature_names, rotation=45, ha='right')
                ax.set_yticks(range(self.n_components))
                ax.set_yticklabels([f'PC{i+1}' for i in range(self.n_components)])
                plt.colorbar(im, ax=ax)
            ax.set_title('Feature Loadings in Principal Components')
            plt.tight_layout()
            plt.savefig(output_path / "feature_loadings_heatmap.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"     ✓ Visualizations saved")
            
        except ImportError:
            print(f"     ⚠ Matplotlib/Seaborn not available - skipping visualizations")


def main():
    """
    Example usage of the Mutant PCA Analyzer
    """
    # Example configuration
    wildtype_pdb = "wildtype.pdb"
    mutant_pdb_dict = {
        "mutant_A": "mutant_A.pdb",
        "mutant_B": "mutant_B.pdb", 
        "mutant_C": "mutant_C.pdb"
    }
    
    # Initialize analyzer
    analyzer = MutantPCAAnalyzer(wildtype_pdb, mutant_pdb_dict)
    
    try:
        # Run analysis pipeline
        analyzer.load_and_align_structures()
        analyzer.extract_residue_features()
        analyzer.perform_pca_analysis()
        analyzer.generate_analysis_report()
        
        print(f"\n🎉 Mutant PCA analysis complete!")
        print(f"   → Systematic modeling errors identified through residue-level geometric PCA")
        print(f"   → Global RMSD reported for context (not used for error mapping)")
        print(f"   → Results saved with comprehensive visualizations")
        
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()