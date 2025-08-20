"""
biostructbenchmark/analysis/pca.py
Principal Component Analysis for structural error pattern identification
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from scipy import stats

from biostructbenchmark.core.alignment import ResidueRMSD


@dataclass
class PCAResult:
    """Container for PCA analysis results"""
    explained_variance_ratio: np.ndarray
    principal_components: np.ndarray
    transformed_data: np.ndarray
    feature_names: List[str]
    loadings: np.ndarray
    eigenvalues: np.ndarray
    cumulative_variance: np.ndarray


@dataclass
class StructureOutlier:
    """Container for structure-level outlier information"""
    structure_id: str
    outlier_score: float
    distance_to_centroid: float
    principal_component_scores: np.ndarray
    outlier_type: str  # 'extreme', 'moderate', 'normal'
    contributing_features: List[str]


@dataclass
class ResidueOutlier:
    """Container for residue-level outlier information"""
    position: str  # e.g., "A_146"
    residue_type: str
    outlier_score: float
    rmsd_value: float
    pc_contributions: Dict[int, float]  # PC number -> contribution
    outlier_severity: str  # 'critical', 'significant', 'moderate'


class PCAAnalyzer:
    """Principal Component Analysis for structural error patterns"""
    
    def __init__(self, n_components: Optional[int] = None, 
                 outlier_threshold: float = 2.0):
        """
        Initialize PCA analyzer
        
        Args:
            n_components: Number of principal components (None for auto)
            outlier_threshold: Z-score threshold for outlier detection
        """
        self.n_components = n_components
        self.outlier_threshold = outlier_threshold
        self.scaler = StandardScaler()
        self.pca = None
        self.structure_data = {}
        self.residue_data = {}
        
    def add_structure_comparison(self, structure_id: str, 
                               residue_rmsds: List[ResidueRMSD],
                               metadata: Optional[Dict] = None) -> None:
        """
        Add structure comparison data for PCA analysis
        
        Args:
            structure_id: Identifier for the structure comparison
            residue_rmsds: Per-residue RMSD data
            metadata: Additional structure metadata
        """
        # Extract features for this structure
        features = self._extract_structure_features(residue_rmsds)
        
        self.structure_data[structure_id] = {
            'features': features,
            'residue_rmsds': residue_rmsds,
            'metadata': metadata or {}
        }
        
        # Store residue-level data
        for rmsd_data in residue_rmsds:
            position_key = f"{rmsd_data.chain_id}_{rmsd_data.position}"
            if position_key not in self.residue_data:
                self.residue_data[position_key] = {
                    'residue_type': rmsd_data.residue_type,
                    'rmsds': [],
                    'structure_ids': []
                }
            
            self.residue_data[position_key]['rmsds'].append(rmsd_data.rmsd)
            self.residue_data[position_key]['structure_ids'].append(structure_id)
    
    def _extract_structure_features(self, residue_rmsds: List[ResidueRMSD]) -> Dict[str, float]:
        """Extract statistical features from residue RMSD data"""
        rmsds = [r.rmsd for r in residue_rmsds]
        
        if not rmsds:
            return {}
        
        # Basic statistics
        features = {
            'mean_rmsd': np.mean(rmsds),
            'median_rmsd': np.median(rmsds),
            'std_rmsd': np.std(rmsds),
            'max_rmsd': np.max(rmsds),
            'min_rmsd': np.min(rmsds),
            'q75_rmsd': np.percentile(rmsds, 75),
            'q25_rmsd': np.percentile(rmsds, 25),
            'iqr_rmsd': np.percentile(rmsds, 75) - np.percentile(rmsds, 25),
            'skewness': stats.skew(rmsds),
            'kurtosis': stats.kurtosis(rmsds),
        }
        
        # Fraction of residues in different RMSD ranges
        total_residues = len(rmsds)
        features.update({
            'frac_low_rmsd': sum(1 for r in rmsds if r < 1.0) / total_residues,
            'frac_medium_rmsd': sum(1 for r in rmsds if 1.0 <= r < 3.0) / total_residues,
            'frac_high_rmsd': sum(1 for r in rmsds if r >= 3.0) / total_residues,
        })
        
        # Separate protein and DNA statistics if available
        protein_rmsds = [r.rmsd for r in residue_rmsds if r.molecule_type == 'protein']
        dna_rmsds = [r.rmsd for r in residue_rmsds if r.molecule_type == 'dna']
        
        if protein_rmsds:
            features.update({
                'protein_mean_rmsd': np.mean(protein_rmsds),
                'protein_max_rmsd': np.max(protein_rmsds),
                'protein_std_rmsd': np.std(protein_rmsds),
            })
        
        if dna_rmsds:
            features.update({
                'dna_mean_rmsd': np.mean(dna_rmsds),
                'dna_max_rmsd': np.max(dna_rmsds),
                'dna_std_rmsd': np.std(dna_rmsds),
            })
        
        return features
    
    def perform_structure_pca(self, variance_threshold: float = 0.95) -> PCAResult:
        """
        Perform PCA on structure-level features
        
        Args:
            variance_threshold: Cumulative variance threshold for component selection
            
        Returns:
            PCA analysis results
        """
        if len(self.structure_data) < 3:
            raise ValueError("Need at least 3 structures for meaningful PCA")
        
        # Prepare feature matrix
        structure_ids = list(self.structure_data.keys())
        feature_names = list(next(iter(self.structure_data.values()))['features'].keys())
        
        # Create feature matrix
        feature_matrix = []
        for struct_id in structure_ids:
            features = self.structure_data[struct_id]['features']
            feature_vector = [features.get(fname, 0.0) for fname in feature_names]
            feature_matrix.append(feature_vector)
        
        feature_matrix = np.array(feature_matrix)
        
        # Handle missing values
        feature_matrix = np.nan_to_num(feature_matrix, nan=0.0)
        
        # Standardize features
        scaled_features = self.scaler.fit_transform(feature_matrix)
        
        # Determine number of components
        if self.n_components is None:
            # Use variance threshold to determine components
            pca_temp = PCA()
            pca_temp.fit(scaled_features)
            cumvar = np.cumsum(pca_temp.explained_variance_ratio_)
            n_comp = np.argmax(cumvar >= variance_threshold) + 1
            n_comp = min(n_comp, len(feature_names), len(structure_ids))
        else:
            n_comp = min(self.n_components, len(feature_names), len(structure_ids))
        
        # Perform PCA
        self.pca = PCA(n_components=n_comp)
        transformed_data = self.pca.fit_transform(scaled_features)
        
        # Calculate loadings (correlations between original features and PCs)
        loadings = self.pca.components_.T * np.sqrt(self.pca.explained_variance_)
        
        return PCAResult(
            explained_variance_ratio=self.pca.explained_variance_ratio_,
            principal_components=self.pca.components_,
            transformed_data=transformed_data,
            feature_names=feature_names,
            loadings=loadings,
            eigenvalues=self.pca.explained_variance_,
            cumulative_variance=np.cumsum(self.pca.explained_variance_ratio_)
        )
    
    def identify_structure_outliers(self, pca_result: PCAResult, 
                                  method: str = 'zscore') -> List[StructureOutlier]:
        """
        Identify outlier structures based on PCA results
        
        Args:
            pca_result: Results from PCA analysis
            method: Outlier detection method ('zscore', 'isolation', 'dbscan')
            
        Returns:
            List of structure outliers
        """
        structure_ids = list(self.structure_data.keys())
        outliers = []
        
        if method == 'zscore':
            # Use distance from centroid in PC space
            centroid = np.mean(pca_result.transformed_data, axis=0)
            
            for i, struct_id in enumerate(structure_ids):
                pc_scores = pca_result.transformed_data[i]
                distance = np.linalg.norm(pc_scores - centroid)
                
                # Calculate z-score of distance
                distances = [np.linalg.norm(pca_result.transformed_data[j] - centroid) 
                           for j in range(len(structure_ids))]
                z_score = (distance - np.mean(distances)) / np.std(distances)
                
                # Determine outlier type
                if abs(z_score) > self.outlier_threshold * 1.5:
                    outlier_type = 'extreme'
                elif abs(z_score) > self.outlier_threshold:
                    outlier_type = 'moderate'
                else:
                    outlier_type = 'normal'
                
                # Identify contributing features
                contributing_features = self._identify_contributing_features(
                    pc_scores, pca_result, top_n=3
                )
                
                outlier = StructureOutlier(
                    structure_id=struct_id,
                    outlier_score=abs(z_score),
                    distance_to_centroid=distance,
                    principal_component_scores=pc_scores,
                    outlier_type=outlier_type,
                    contributing_features=contributing_features
                )
                outliers.append(outlier)
        
        # Sort by outlier score (highest first)
        outliers.sort(key=lambda x: x.outlier_score, reverse=True)
        
        return outliers
    
    def _identify_contributing_features(self, pc_scores: np.ndarray, 
                                      pca_result: PCAResult, 
                                      top_n: int = 3) -> List[str]:
        """Identify features that contribute most to outlier status"""
        contributions = []
        
        for i, pc_score in enumerate(pc_scores):
            if i >= len(pca_result.loadings[0]):
                break
                
            pc_loadings = pca_result.loadings[:, i]
            weighted_contributions = pc_loadings * pc_score
            
            for j, contrib in enumerate(weighted_contributions):
                if j < len(pca_result.feature_names):
                    contributions.append((
                        pca_result.feature_names[j], 
                        abs(contrib)
                    ))
        
        # Sort by contribution magnitude and return top features
        contributions.sort(key=lambda x: x[1], reverse=True)
        return [feature for feature, _ in contributions[:top_n]]
    
    def perform_residue_pca(self, min_structures: int = 3) -> Tuple[PCAResult, Dict]:
        """
        Perform PCA on residue-level error patterns
        
        Args:
            min_structures: Minimum number of structures required for a residue
            
        Returns:
            PCA results and residue position mapping
        """
        # Filter residues that appear in enough structures
        valid_positions = {
            pos: data for pos, data in self.residue_data.items()
            if len(data['rmsds']) >= min_structures
        }
        
        if len(valid_positions) < 10:
            raise ValueError("Need at least 10 residue positions for PCA")
        
        # Create RMSD matrix (positions x structures)
        position_names = list(valid_positions.keys())
        structure_ids = list(self.structure_data.keys())
        
        rmsd_matrix = np.full((len(position_names), len(structure_ids)), np.nan)
        
        for i, pos in enumerate(position_names):
            pos_data = valid_positions[pos]
            for j, rmsd in enumerate(pos_data['rmsds']):
                struct_idx = structure_ids.index(pos_data['structure_ids'][j])
                rmsd_matrix[i, struct_idx] = rmsd
        
        # Handle missing values (interpolate or use mean)
        for i in range(rmsd_matrix.shape[0]):
            row = rmsd_matrix[i, :]
            valid_mask = ~np.isnan(row)
            if np.any(valid_mask):
                mean_val = np.mean(row[valid_mask])
                rmsd_matrix[i, np.isnan(row)] = mean_val
        
        # Standardize by position (each residue position)
        scaled_matrix = self.scaler.fit_transform(rmsd_matrix)
        
        # Perform PCA
        n_comp = min(len(structure_ids), len(position_names), 10)
        pca = PCA(n_components=n_comp)
        transformed_data = pca.fit_transform(scaled_matrix)
        
        # Calculate loadings
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        
        result = PCAResult(
            explained_variance_ratio=pca.explained_variance_ratio_,
            principal_components=pca.components_,
            transformed_data=transformed_data,
            feature_names=structure_ids,  # In this case, structures are features
            loadings=loadings,
            eigenvalues=pca.explained_variance_,
            cumulative_variance=np.cumsum(pca.explained_variance_ratio_)
        )
        
        return result, {i: pos for i, pos in enumerate(position_names)}
    
    def identify_residue_outliers(self, residue_pca_result: PCAResult,
                                position_mapping: Dict,
                                rmsd_threshold: float = 3.0) -> List[ResidueOutlier]:
        """
        Identify residue positions with unusual error patterns
        
        Args:
            residue_pca_result: PCA results for residues
            position_mapping: Mapping from index to position name
            rmsd_threshold: RMSD threshold for high errors
            
        Returns:
            List of residue outliers
        """
        outliers = []
        
        # Calculate distance from centroid in PC space
        centroid = np.mean(residue_pca_result.transformed_data, axis=0)
        
        for i, pc_scores in enumerate(residue_pca_result.transformed_data):
            position = position_mapping[i]
            residue_type = self.residue_data[position]['residue_type']
            mean_rmsd = np.mean(self.residue_data[position]['rmsds'])
            
            # Calculate outlier score
            distance = np.linalg.norm(pc_scores - centroid)
            
            # Calculate PC contributions
            pc_contributions = {}
            for j, score in enumerate(pc_scores):
                if j < len(residue_pca_result.explained_variance_ratio):
                    weight = residue_pca_result.explained_variance_ratio[j]
                    pc_contributions[j] = abs(score) * weight
            
            # Determine severity
            if mean_rmsd > rmsd_threshold and distance > 2.0:
                severity = 'critical'
            elif mean_rmsd > rmsd_threshold or distance > 1.5:
                severity = 'significant'
            else:
                severity = 'moderate'
            
            outlier = ResidueOutlier(
                position=position,
                residue_type=residue_type,
                outlier_score=distance,
                rmsd_value=mean_rmsd,
                pc_contributions=pc_contributions,
                outlier_severity=severity
            )
            outliers.append(outlier)
        
        # Sort by outlier score
        outliers.sort(key=lambda x: x.outlier_score, reverse=True)
        
        return outliers
    
    
    def export_results(self, pca_result: PCAResult, outliers: List,
                      output_dir: Path, prefix: str = "pca") -> Dict[str, Path]:
        """
        Export PCA analysis results to files
        
        Args:
            pca_result: PCA analysis results
            outliers: List of identified outliers
            output_dir: Directory for output files
            prefix: Prefix for output filenames
            
        Returns:
            Dictionary mapping result type to file path
        """
        output_dir.mkdir(exist_ok=True)
        exported_files = {}
        
        # 1. Export PCA summary
        summary_data = {
            'component': [f'PC{i+1}' for i in range(len(pca_result.explained_variance_ratio))],
            'explained_variance_ratio': pca_result.explained_variance_ratio,
            'eigenvalue': pca_result.eigenvalues,
            'cumulative_variance': pca_result.cumulative_variance
        }
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / f"{prefix}_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        exported_files['summary'] = summary_file
        
        # 2. Export transformed data (PC scores)
        if hasattr(self, 'structure_data') and len(self.structure_data) > 0:
            # Structure-level PCA
            structure_ids = list(self.structure_data.keys())
            scores_data = {
                'structure_id': structure_ids
            }
            
            for i in range(pca_result.transformed_data.shape[1]):
                scores_data[f'PC{i+1}'] = pca_result.transformed_data[:, i]
        else:
            # Residue-level PCA - use indices as IDs
            scores_data = {
                'position_index': list(range(pca_result.transformed_data.shape[0]))
            }
            
            for i in range(pca_result.transformed_data.shape[1]):
                scores_data[f'PC{i+1}'] = pca_result.transformed_data[:, i]
        
        scores_df = pd.DataFrame(scores_data)
        scores_file = output_dir / f"{prefix}_scores.csv"
        scores_df.to_csv(scores_file, index=False)
        exported_files['scores'] = scores_file
        
        # 3. Export loadings
        loadings_data = {
            'feature': pca_result.feature_names
        }
        
        for i in range(pca_result.loadings.shape[1]):
            loadings_data[f'PC{i+1}_loading'] = pca_result.loadings[:, i]
        
        loadings_df = pd.DataFrame(loadings_data)
        loadings_file = output_dir / f"{prefix}_loadings.csv"
        loadings_df.to_csv(loadings_file, index=False)
        exported_files['loadings'] = loadings_file
        
        # 4. Export outliers
        if outliers:
            outlier_data = []
            for outlier in outliers:
                if hasattr(outlier, 'structure_id'):  # Structure outlier
                    outlier_data.append({
                        'id': outlier.structure_id,
                        'type': 'structure',
                        'outlier_score': outlier.outlier_score,
                        'outlier_type': outlier.outlier_type,
                        'distance_to_centroid': outlier.distance_to_centroid,
                        'contributing_features': ';'.join(outlier.contributing_features)
                    })
                else:  # Residue outlier
                    outlier_data.append({
                        'id': outlier.position,
                        'type': 'residue',
                        'outlier_score': outlier.outlier_score,
                        'outlier_type': outlier.outlier_severity,
                        'rmsd_value': outlier.rmsd_value,
                        'residue_type': outlier.residue_type
                    })
            
            outliers_df = pd.DataFrame(outlier_data)
            outliers_file = output_dir / f"{prefix}_outliers.csv"
            outliers_df.to_csv(outliers_file, index=False)
            exported_files['outliers'] = outliers_file
        
        return exported_files
    
    def generate_analysis_report(self, pca_result: PCAResult, 
                               structure_outliers: List[StructureOutlier],
                               output_path: Path) -> None:
        """
        Generate a comprehensive analysis report
        
        Args:
            pca_result: PCA analysis results
            structure_outliers: List of structure outliers
            output_path: Path for the report file
        """
        with open(output_path, 'w') as f:
            f.write("PCA Analysis Report\\n")
            f.write("=" * 50 + "\\n\\n")
            
            # Summary statistics
            f.write("ANALYSIS SUMMARY\\n")
            f.write("-" * 20 + "\\n")
            f.write(f"Number of structures analyzed: {len(self.structure_data)}\\n")
            f.write(f"Number of features: {len(pca_result.feature_names)}\\n")
            f.write(f"Number of principal components: {len(pca_result.explained_variance_ratio)}\\n")
            f.write(f"Total variance explained: {pca_result.cumulative_variance[-1]:.1%}\\n\\n")
            
            # Principal components
            f.write("PRINCIPAL COMPONENTS\\n")
            f.write("-" * 20 + "\\n")
            for i, (var_ratio, eigenval) in enumerate(zip(pca_result.explained_variance_ratio, 
                                                        pca_result.eigenvalues)):
                f.write(f"PC{i+1}: {var_ratio:.1%} variance (eigenvalue: {eigenval:.3f})\\n")
            f.write("\\n")
            
            # Outliers
            f.write("STRUCTURE OUTLIERS\\n")
            f.write("-" * 20 + "\\n")
            
            extreme_outliers = [o for o in structure_outliers if o.outlier_type == 'extreme']
            moderate_outliers = [o for o in structure_outliers if o.outlier_type == 'moderate']
            
            f.write(f"Extreme outliers: {len(extreme_outliers)}\\n")
            f.write(f"Moderate outliers: {len(moderate_outliers)}\\n\\n")
            
            if extreme_outliers:
                f.write("TOP EXTREME OUTLIERS:\\n")
                for outlier in extreme_outliers[:5]:
                    f.write(f"  {outlier.structure_id}: score={outlier.outlier_score:.2f}, ")
                    f.write(f"features={', '.join(outlier.contributing_features)}\\n")
                f.write("\\n")
            
            # Feature importance
            f.write("FEATURE IMPORTANCE (PC1)\\n")
            f.write("-" * 20 + "\\n")
            if len(pca_result.loadings) > 0:
                pc1_loadings = pca_result.loadings[:, 0]
                feature_importance = list(zip(pca_result.feature_names, abs(pc1_loadings)))
                feature_importance.sort(key=lambda x: x[1], reverse=True)
                
                for feature, importance in feature_importance[:10]:
                    f.write(f"  {feature}: {importance:.3f}\\n")