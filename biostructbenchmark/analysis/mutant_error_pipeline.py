#!/usr/bin/env python3
"""
Focused Mutant Error Analysis Pipeline
Identifies systematic modeling errors in protein mutants using residue-level geometric PCA
with specific focus on DNA interface validation.

Critical Implementation:
- Uses Bio.PDB.Superimposer for CA alignment (NOT Bio.PDB.RMSD)
- Handles missing dihedrals with np.nan (no row dropping)
- DNA features only if DNA chains present
- PCA with 2 components, validation against x3DNA-DSSR interface
"""

import numpy as np
import pandas as pd
import json
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import matplotlib.pyplot as plt
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Biopython imports
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_dihedral


class MutantErrorPipeline:
    """
    Focused pipeline for systematic modeling error analysis in protein mutants.
    
    Implements residue-level geometric PCA following strict bioinformatics protocols:
    1. CA-based structural alignment to wild-type
    2. Backbone dihedral extraction (phi/psi)
    3. x3DNA-DSSR DNA interface features (if DNA present)
    4. PCA analysis with validation against interface residues
    """
    
    def __init__(self, wildtype_pdb: str, mutant_pdb_dict: Dict[str, str]):
        """
        Initialize the focused mutant error pipeline.
        
        Args:
            wildtype_pdb: Path to wild-type reference PDB
            mutant_pdb_dict: {mutant_id: pdb_path} dictionary
        """
        self.wildtype_pdb = Path(wildtype_pdb)
        self.mutant_pdb_dict = {k: Path(v) for k, v in mutant_pdb_dict.items()}
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        
        # Analysis storage
        self.wildtype_structure = None
        self.mutant_structures = {}
        self.common_residues = []
        self.global_rmsds = {}
        self.feature_matrix = None
        self.residue_ids = []
        self.pca_model = None
        self.scaler = None
        self.dna_interface_residues = set()
        
        print(f"🧬 Mutant Error Analysis Pipeline")
        print(f"   Wild-type: {self.wildtype_pdb}")
        print(f"   Mutants: {list(self.mutant_pdb_dict.keys())}")

    def _get_parser_for_file(self, file_path: Path):
        """Get appropriate parser based on file extension"""
        if file_path.suffix.lower() in ['.cif', '.mmcif']:
            return self.cif_parser
        else:
            return self.pdb_parser

    def step1_load_and_align(self) -> None:
        """
        Step 1: Load structures and align mutants to wild-type via CA atoms.
        
        - Identifies common residues across all structures
        - Performs CA-based superimposition using Bio.PDB.Superimposer
        - Reports global RMSD for supplementary context only
        """
        print(f"\n📂 Step 1: Loading structures and CA alignment...")
        
        # Load wild-type
        wt_parser = self._get_parser_for_file(self.wildtype_pdb)
        self.wildtype_structure = wt_parser.get_structure("WT", self.wildtype_pdb)
        wt_residues = self._get_protein_residues(self.wildtype_structure)
        wt_residue_map = {(res.get_parent().id, res.id[1]): res for res in wt_residues}
        
        print(f"   ✓ Wild-type loaded: {len(wt_residues)} protein residues")
        
        # Load mutants and find common residues
        common_residue_keys = set(wt_residue_map.keys())
        
        for mutant_id, pdb_path in self.mutant_pdb_dict.items():
            mutant_parser = self._get_parser_for_file(Path(pdb_path))
            structure = mutant_parser.get_structure(mutant_id, pdb_path)
            self.mutant_structures[mutant_id] = structure
            
            mut_residues = self._get_protein_residues(structure)
            mut_residue_map = {(res.get_parent().id, res.id[1]): res for res in mut_residues}
            
            # Intersection for common residues
            common_residue_keys &= set(mut_residue_map.keys())
            print(f"   ✓ {mutant_id} loaded: {len(mut_residues)} protein residues")
        
        self.common_residues = sorted(common_residue_keys)
        print(f"   ✓ Common residues for alignment: {len(self.common_residues)}")
        
        # Perform CA-based alignment
        self._align_mutants_to_wildtype()
        
        # Report global RMSD (supplementary context only)
        rmsd_values = [f"{mid}={rmsd:.2f}Å" for mid, rmsd in self.global_rmsds.items()]
        print(f"\n📊 RMSD: {', '.join(rmsd_values)}")
        print(f"   (Global RMSD for context - NOT used for error mapping)")

    def _get_protein_residues(self, structure: Structure) -> List[Residue]:
        """Extract standard protein residues"""
        standard_aa = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        protein_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if (residue.id[0] == ' ' and 
                        residue.get_resname() in standard_aa):
                        protein_residues.append(residue)
        
        return protein_residues

    def _align_mutants_to_wildtype(self) -> None:
        """Align mutant structures to wild-type using CA atoms"""
        # Extract WT CA atoms for common residues
        wt_ca_atoms = []
        for chain_id, res_id in self.common_residues:
            try:
                residue = self.wildtype_structure[0][chain_id][res_id]
                if 'CA' in residue:
                    wt_ca_atoms.append(residue['CA'])
            except KeyError:
                continue
        
        print(f"   → CA alignment with {len(wt_ca_atoms)} atoms")
        
        # Align each mutant using Bio.PDB.Superimposer
        for mutant_id, structure in self.mutant_structures.items():
            mutant_ca_atoms = []
            for chain_id, res_id in self.common_residues:
                try:
                    residue = structure[0][chain_id][res_id]
                    if 'CA' in residue:
                        mutant_ca_atoms.append(residue['CA'])
                except KeyError:
                    continue
            
            # Ensure equal number of atoms
            min_atoms = min(len(wt_ca_atoms), len(mutant_ca_atoms))
            if min_atoms < 10:
                warnings.warn(f"Insufficient CA atoms ({min_atoms}) for {mutant_id}")
                continue
            
            # Perform superimposition
            superimposer = Superimposer()
            superimposer.set_atoms(wt_ca_atoms[:min_atoms], mutant_ca_atoms[:min_atoms])
            
            # Apply transformation to entire structure
            superimposer.apply(structure.get_atoms())
            
            # Store global RMSD
            self.global_rmsds[mutant_id] = superimposer.rms
        
        print(f"   ✓ Structural alignment complete")

    def step2_extract_features(self) -> None:
        """
        Step 2: Extract residue-level features.
        
        Features per residue:
        - phi, psi backbone dihedrals
        - minor_groove_width, base_pair_rise (if DNA present)
        
        Handles missing dihedrals with np.nan (no row dropping).
        """
        print(f"\n🔬 Step 2: Extracting residue-level features...")
        
        # Check DNA presence across all structures
        has_dna = any(self._check_dna_presence(struct) 
                     for struct in [self.wildtype_structure] + list(self.mutant_structures.values()))
        
        if has_dna:
            print(f"   ✓ DNA detected - will extract interface features")
            self._extract_dssr_interface_residues()
        else:
            print(f"   ℹ No DNA chains - skipping DNA features")
        
        # Extract features for all structures
        all_features = []
        
        # Process each structure
        for struct_id in ["WT"] + list(self.mutant_pdb_dict.keys()):
            if struct_id == "WT":
                structure = self.wildtype_structure
            else:
                structure = self.mutant_structures[struct_id]
            
            struct_features = self._extract_structure_features(structure, struct_id, has_dna)
            all_features.append(struct_features)
        
        # Build feature matrix
        self._build_feature_matrix(all_features, has_dna)
        print(f"   ✓ Feature matrix built: {self.feature_matrix.shape}")

    def _check_dna_presence(self, structure: Structure) -> bool:
        """Check if structure contains DNA chains"""
        dna_residues = {'A', 'T', 'G', 'C', 'DA', 'DT', 'DG', 'DC'}
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname().strip() in dna_residues:
                        return True
        return False

    def _extract_dssr_interface_residues(self) -> None:
        """Extract DNA interface residues using x3DNA-DSSR"""
        print(f"   → Running x3DNA-DSSR analysis...")
        
        # Use wild-type for interface definition
        temp_pdb = "temp_wt_dssr.pdb"
        io = PDBIO()
        io.set_structure(self.wildtype_structure)
        io.save(temp_pdb)
        
        try:
            # Run DSSR
            dssr_output = "dna_features.json"
            cmd = f"x3dna-dssr --input={temp_pdb} --output={dssr_output} --json"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and Path(dssr_output).exists():
                with open(dssr_output, 'r') as f:
                    dssr_data = json.load(f)
                
                # Parse interface residues
                if 'interface_residues' in dssr_data:
                    for res_info in dssr_data['interface_residues']:
                        if 'residue_id' in res_info:
                            self.dna_interface_residues.add(res_info['residue_id'])
                
                print(f"     ✓ DSSR identified {len(self.dna_interface_residues)} interface residues")
                
                # Cleanup
                Path(temp_pdb).unlink(missing_ok=True)
                Path(dssr_output).unlink(missing_ok=True)
            else:
                print(f"     ⚠ DSSR failed: {result.stderr}")
                
        except Exception as e:
            print(f"     ⚠ DSSR extraction failed: {e}")
        finally:
            # Cleanup temp files
            Path(temp_pdb).unlink(missing_ok=True)
            Path("dna_features.json").unlink(missing_ok=True)

    def _extract_structure_features(self, structure: Structure, struct_id: str, has_dna: bool) -> List[float]:
        """Extract features for a single structure"""
        features = []
        
        for chain_id, res_id in self.common_residues:
            try:
                residue = structure[0][chain_id][res_id]
                
                # Calculate backbone dihedrals
                phi, psi = self._calculate_dihedrals(residue)
                
                # DNA features if present
                if has_dna:
                    minor_groove, base_rise = self._get_dna_features(residue, structure)
                    residue_features = [phi, psi, minor_groove, base_rise]
                else:
                    residue_features = [phi, psi]
                
                features.extend(residue_features)
                
            except KeyError:
                # Handle missing residues with np.nan
                if has_dna:
                    features.extend([np.nan, np.nan, np.nan, np.nan])
                else:
                    features.extend([np.nan, np.nan])
        
        return features

    def _calculate_dihedrals(self, residue: Residue) -> Tuple[float, float]:
        """Calculate phi and psi dihedrals, return np.nan if missing"""
        try:
            chain = residue.get_parent()
            residues = list(chain)
            current_idx = residues.index(residue)
            
            phi = psi = np.nan
            
            # Calculate phi
            if current_idx > 0:
                prev_res = residues[current_idx - 1]
                if (all(atom in prev_res for atom in ['C']) and 
                    all(atom in residue for atom in ['N', 'CA', 'C'])):
                    phi = calc_dihedral(
                        prev_res['C'].get_vector(),
                        residue['N'].get_vector(),
                        residue['CA'].get_vector(),
                        residue['C'].get_vector()
                    )
                    phi = np.degrees(phi)
            
            # Calculate psi
            if current_idx < len(residues) - 1:
                next_res = residues[current_idx + 1]
                if (all(atom in residue for atom in ['N', 'CA', 'C']) and 
                    'N' in next_res):
                    psi = calc_dihedral(
                        residue['N'].get_vector(),
                        residue['CA'].get_vector(),
                        residue['C'].get_vector(),
                        next_res['N'].get_vector()
                    )
                    psi = np.degrees(psi)
            
        except (KeyError, ValueError, IndexError):
            phi = psi = np.nan
        
        return phi, psi

    def _get_dna_features(self, residue: Residue, structure: Structure) -> Tuple[float, float]:
        """Get DNA interface features for residue"""
        # Simplified DNA features - in real implementation would use DSSR output
        # For now, calculate distance to DNA as proxy
        
        min_dna_distance = self._calculate_min_dna_distance(residue, structure)
        
        # Placeholder DNA features (would come from DSSR JSON)
        minor_groove_width = 12.0 if min_dna_distance and min_dna_distance < 5.0 else np.nan
        base_pair_rise = 3.4 if min_dna_distance and min_dna_distance < 5.0 else np.nan
        
        return minor_groove_width, base_pair_rise

    def _calculate_min_dna_distance(self, protein_residue: Residue, structure: Structure) -> Optional[float]:
        """Calculate minimum distance to DNA"""
        if 'CA' not in protein_residue:
            return None
        
        protein_ca = protein_residue['CA']
        min_dist = None
        
        dna_residues = {'A', 'T', 'G', 'C', 'DA', 'DT', 'DG', 'DC'}
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname().strip() in dna_residues and 'P' in residue:
                        dist = protein_ca - residue['P']
                        if min_dist is None or dist < min_dist:
                            min_dist = dist
        
        return min_dist

    def _build_feature_matrix(self, all_features: List[List[float]], has_dna: bool) -> None:
        """Build feature matrix from extracted features"""
        n_structures = len(all_features)
        n_residues = len(self.common_residues)
        n_features_per_residue = 4 if has_dna else 2
        
        # Reshape to (structures, residues * features_per_residue)
        self.feature_matrix = np.array(all_features)
        
        # Create residue IDs for reference
        self.residue_ids = []
        for chain_id, res_id in self.common_residues:
            self.residue_ids.append(f"{chain_id}_{res_id}")
        
        print(f"     → Matrix shape: {self.feature_matrix.shape}")
        print(f"     → Features per residue: {n_features_per_residue}")

    def step3_pca_analysis(self) -> None:
        """
        Step 3: Build feature matrix and perform PCA.
        
        - Scales features using StandardScaler
        - Runs PCA with 2 components
        - Saves PCA loadings as error_hotspots.csv
        """
        print(f"\n📈 Step 3: PCA Analysis...")
        
        # Handle NaN values by filling with feature means
        feature_matrix_clean = np.copy(self.feature_matrix)
        
        # Fill NaN with column means (handling empty columns)
        for col in range(feature_matrix_clean.shape[1]):
            col_data = feature_matrix_clean[:, col]
            if np.all(np.isnan(col_data)):
                # If entire column is NaN, fill with 0
                feature_matrix_clean[:, col] = 0.0
            else:
                col_mean = np.nanmean(col_data)
                if np.isnan(col_mean):
                    col_mean = 0.0
                feature_matrix_clean[:, col] = np.where(np.isnan(col_data), col_mean, col_data)
        
        # Final check for any remaining NaNs
        if np.any(np.isnan(feature_matrix_clean)):
            print(f"   ⚠ Replacing remaining NaN values with 0")
            feature_matrix_clean = np.nan_to_num(feature_matrix_clean, nan=0.0)
        
        # Scale features
        self.scaler = StandardScaler()
        scaled_features = self.scaler.fit_transform(feature_matrix_clean)
        
        # PCA with 2 components
        self.pca_model = PCA(n_components=2)
        pca_scores = self.pca_model.fit_transform(scaled_features)
        
        # Report explained variance
        var_explained = self.pca_model.explained_variance_ratio_
        print(f"   ✓ PCA complete:")
        print(f"     PC1: {var_explained[0]:.1%} variance explained")
        print(f"     PC2: {var_explained[1]:.1%} variance explained")
        print(f"     Total: {sum(var_explained):.1%} variance explained")
        
        # Save PCA loadings as error hotspots
        self._save_error_hotspots()
        
        # Store for visualization
        self.pca_scores = pca_scores

    def _save_error_hotspots(self) -> None:
        """Save PCA loadings as error_hotspots.csv"""
        # Calculate residue-level loadings
        n_features_per_residue = 4 if len(self.pca_model.components_[0]) // len(self.residue_ids) == 4 else 2
        
        residue_loadings = []
        
        for i, residue_id in enumerate(self.residue_ids):
            # Average loadings across features for this residue
            start_idx = i * n_features_per_residue
            end_idx = start_idx + n_features_per_residue
            
            pc1_loading = np.mean(np.abs(self.pca_model.components_[0][start_idx:end_idx]))
            pc2_loading = np.mean(np.abs(self.pca_model.components_[1][start_idx:end_idx]))
            
            residue_loadings.append({
                'residue_id': residue_id,
                'PC1_loading': pc1_loading,
                'PC2_loading': pc2_loading,
                'PC1_abs_loading': abs(pc1_loading)
            })
        
        # Sort by PC1 absolute loading
        residue_loadings.sort(key=lambda x: x['PC1_abs_loading'], reverse=True)
        
        # Save to CSV
        df = pd.DataFrame(residue_loadings)
        df.to_csv("error_hotspots.csv", index=False)
        
        print(f"   ✓ Error hotspots saved to error_hotspots.csv")
        
        # Store top hotspots for validation
        self.top_pc1_hotspots = [item['residue_id'] for item in residue_loadings[:5]]

    def step4_validate_dna_interface(self) -> None:
        """
        Step 4: Validate with x3DNA-DSSR interface overlap.
        
        Checks if top 5 PC1 hotspots overlap with DNA interface residues.
        """
        print(f"\n🔍 Step 4: DNA Interface Validation...")
        
        if not self.dna_interface_residues:
            print(f"   ℹ No DNA interface residues identified - skipping validation")
            return
        
        # Check overlap between top PC1 hotspots and DNA interface
        interface_overlap = []
        
        for hotspot_id in self.top_pc1_hotspots:
            # Convert residue format for comparison
            if hotspot_id in self.dna_interface_residues:
                interface_overlap.append(hotspot_id)
        
        # Print validation result
        if interface_overlap:
            print(f"   ✓ ERROR HOTSPOTS IN DNA INTERFACE: {interface_overlap}")
        else:
            print(f"   ℹ No overlap between top PC1 hotspots and DNA interface")
        
        print(f"     Top 5 PC1 hotspots: {self.top_pc1_hotspots}")
        print(f"     DNA interface residues: {list(self.dna_interface_residues)[:10]}...")

    def generate_visualizations(self) -> None:
        """Generate required visualizations"""
        print(f"\n📊 Generating visualizations...")
        
        # PCA scatter plot
        self._plot_pca_scatter()
        
        # Top 5 residues bar plot
        self._plot_top_residues()
        
        print(f"   ✓ Visualizations saved")

    def _plot_pca_scatter(self) -> None:
        """PCA scatter plot with mutant labels"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot points with labels
        mutant_ids = ["WT"] + list(self.mutant_pdb_dict.keys())
        colors = plt.cm.Set1(np.linspace(0, 1, len(mutant_ids)))
        
        for i, (mutant_id, color) in enumerate(zip(mutant_ids, colors)):
            ax.scatter(self.pca_scores[i, 0], self.pca_scores[i, 1], 
                      c=[color], s=100, label=mutant_id, alpha=0.8)
            ax.annotate(mutant_id, (self.pca_scores[i, 0], self.pca_scores[i, 1]),
                       xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        # Labels with variance explained
        var_explained = self.pca_model.explained_variance_ratio_
        ax.set_xlabel(f'PC1 ({var_explained[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({var_explained[1]:.1%} variance)')
        ax.set_title('PCA: Systematic Modeling Error Analysis')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig("pca_scatter_plot.png", dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_top_residues(self) -> None:
        """Bar plot of top 5 residues by PC1 absolute loading"""
        # Read error hotspots
        df = pd.read_csv("error_hotspots.csv")
        top_5 = df.head(5)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        bars = ax.bar(range(len(top_5)), top_5['PC1_abs_loading'], 
                     color='steelblue', alpha=0.7)
        
        ax.set_xlabel('Residue')
        ax.set_ylabel('PC1 Absolute Loading')
        ax.set_title('Top 5 Error Hotspots by PC1 Loading')
        ax.set_xticks(range(len(top_5)))
        ax.set_xticklabels(top_5['residue_id'], rotation=45)
        
        # Add value labels on bars
        for bar, value in zip(bars, top_5['PC1_abs_loading']):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                   f'{value:.3f}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        plt.savefig("top_residues_pc1_loading.png", dpi=300, bbox_inches='tight')
        plt.close()

    def run_pipeline(self) -> None:
        """Run the complete analysis pipeline"""
        try:
            self.step1_load_and_align()
            self.step2_extract_features()
            self.step3_pca_analysis()
            self.step4_validate_dna_interface()
            self.generate_visualizations()
            
            print(f"\n🎉 Mutant Error Analysis Complete!")
            print(f"   → Outputs: error_hotspots.csv, pca_scatter_plot.png, top_residues_pc1_loading.png")
            
        except Exception as e:
            print(f"\n❌ Pipeline failed: {e}")
            raise


def main():
    """Example usage"""
    # Configuration
    wildtype_pdb = "wildtype.pdb"
    mutant_pdb_dict = {
        "mutant_A": "mutant_A.pdb",
        "mutant_B": "mutant_B.pdb",
        "mutant_C": "mutant_C.pdb"
    }
    
    # Run pipeline
    pipeline = MutantErrorPipeline(wildtype_pdb, mutant_pdb_dict)
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()