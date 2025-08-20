"""
biostructbenchmark/visualization/pca_plots.py
Visualization methods for PCA analysis results
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Optional

from biostructbenchmark.analysis.pca import PCAResult, StructureOutlier, ResidueOutlier


class PCAVisualizer:
    """Create visualization plots for PCA analysis results"""
    
    def __init__(self, style: str = 'default', palette: str = 'husl'):
        """
        Initialize PCA visualizer
        
        Args:
            style: Matplotlib style to use
            palette: Seaborn color palette
        """
        plt.style.use(style)
        sns.set_palette(palette)
        self.figure_size = (10, 8)
        self.dpi = 300
    
    def create_scree_plot(self, pca_result: PCAResult, output_path: Path,
                         title: str = "PCA Scree Plot") -> Path:
        """
        Create scree plot showing explained variance
        
        Args:
            pca_result: PCA analysis results
            output_path: Path for saving the plot
            title: Plot title
            
        Returns:
            Path to saved plot
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        n_components = len(pca_result.explained_variance_ratio)
        x = range(1, n_components + 1)
        
        # Bar plot for individual variance
        ax.bar(x, pca_result.explained_variance_ratio * 100, 
               alpha=0.7, label='Individual', color='skyblue')
        
        # Line plot for cumulative variance
        ax.plot(x, pca_result.cumulative_variance * 100, 
                'ro-', label='Cumulative', linewidth=2, markersize=6)
        
        ax.set_xlabel('Principal Component')
        ax.set_ylabel('Explained Variance (%)')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add text annotations for first few components
        for i in range(min(3, len(x))):
            ax.annotate(f'{pca_result.explained_variance_ratio[i]:.1%}',
                       (x[i], pca_result.explained_variance_ratio[i] * 100),
                       textcoords="offset points", xytext=(0,10), ha='center')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def create_pc_scatter_plot(self, pca_result: PCAResult, 
                              structure_ids: List[str],
                              outliers: Optional[List[StructureOutlier]] = None,
                              output_path: Path = None,
                              pc_x: int = 1, pc_y: int = 2,
                              title: str = None) -> Path:
        """
        Create scatter plot of principal components
        
        Args:
            pca_result: PCA analysis results
            structure_ids: List of structure identifiers
            outliers: Optional list of structure outliers for highlighting
            output_path: Path for saving the plot
            pc_x: Principal component for x-axis (1-indexed)
            pc_y: Principal component for y-axis (1-indexed)
            title: Plot title
            
        Returns:
            Path to saved plot
        """
        if len(pca_result.transformed_data[0]) < max(pc_x, pc_y):
            raise ValueError(f"Not enough components for PC{pc_x} vs PC{pc_y}")
        
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        pc_x_idx = pc_x - 1
        pc_y_idx = pc_y - 1
        
        x_data = pca_result.transformed_data[:, pc_x_idx]
        y_data = pca_result.transformed_data[:, pc_y_idx]
        
        # Color points based on outlier status if provided
        colors = ['blue'] * len(structure_ids)
        sizes = [60] * len(structure_ids)
        
        if outliers:
            outlier_dict = {o.structure_id: o for o in outliers}
            for i, struct_id in enumerate(structure_ids):
                if struct_id in outlier_dict:
                    outlier = outlier_dict[struct_id]
                    if outlier.outlier_type == 'extreme':
                        colors[i] = 'red'
                        sizes[i] = 100
                    elif outlier.outlier_type == 'moderate':
                        colors[i] = 'orange'
                        sizes[i] = 80
        
        scatter = ax.scatter(x_data, y_data, c=colors, s=sizes, alpha=0.7, edgecolors='black')
        
        # Add structure labels
        for i, struct_id in enumerate(structure_ids):
            ax.annotate(struct_id, (x_data[i], y_data[i]), 
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
        
        ax.set_xlabel(f'PC{pc_x} ({pca_result.explained_variance_ratio[pc_x_idx]:.1%} variance)')
        ax.set_ylabel(f'PC{pc_y} ({pca_result.explained_variance_ratio[pc_y_idx]:.1%} variance)')
        
        if title is None:
            title = f'PCA: PC{pc_x} vs PC{pc_y}'
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        # Add legend for outliers if present
        if outliers:
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='blue', label='Normal'),
                Patch(facecolor='orange', label='Moderate Outlier'),
                Patch(facecolor='red', label='Extreme Outlier')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def create_loadings_plot(self, pca_result: PCAResult, output_path: Path,
                           pc_x: int = 1, pc_y: int = 2,
                           max_features: int = 20,
                           title: str = None) -> Path:
        """
        Create loadings plot showing feature contributions
        
        Args:
            pca_result: PCA analysis results
            output_path: Path for saving the plot
            pc_x: Principal component for x-axis (1-indexed)
            pc_y: Principal component for y-axis (1-indexed)
            max_features: Maximum number of features to show
            title: Plot title
            
        Returns:
            Path to saved plot
        """
        if len(pca_result.feature_names) > max_features:
            # Select top contributing features
            pc_x_idx = pc_x - 1
            total_contrib = np.abs(pca_result.loadings[:, pc_x_idx])
            if len(pca_result.loadings[0]) > pc_y - 1:
                total_contrib += np.abs(pca_result.loadings[:, pc_y - 1])
            
            top_indices = np.argsort(total_contrib)[-max_features:]
            selected_features = [pca_result.feature_names[i] for i in top_indices]
            selected_loadings = pca_result.loadings[top_indices, :]
        else:
            selected_features = pca_result.feature_names
            selected_loadings = pca_result.loadings
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        pc_x_idx = pc_x - 1
        pc_y_idx = pc_y - 1
        
        if selected_loadings.shape[1] > pc_x_idx and selected_loadings.shape[1] > pc_y_idx:
            pc_x_loadings = selected_loadings[:, pc_x_idx]
            pc_y_loadings = selected_loadings[:, pc_y_idx]
            
            # Draw arrows and labels
            for i, feature in enumerate(selected_features):
                ax.arrow(0, 0, pc_x_loadings[i], pc_y_loadings[i],
                        head_width=0.02, head_length=0.02, fc='red', ec='red',
                        alpha=0.7, linewidth=1.5)
                
                # Position labels to avoid overlap
                label_x = pc_x_loadings[i] * 1.1
                label_y = pc_y_loadings[i] * 1.1
                
                ax.text(label_x, label_y, feature,
                       fontsize=9, ha='center', va='center',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
            
            ax.set_xlabel(f'PC{pc_x} ({pca_result.explained_variance_ratio[pc_x_idx]:.1%})')
            ax.set_ylabel(f'PC{pc_y} ({pca_result.explained_variance_ratio[pc_y_idx]:.1%})')
            
            if title is None:
                title = f'PCA Loadings: PC{pc_x} vs PC{pc_y}'
            ax.set_title(title)
            
            ax.grid(True, alpha=0.3)
            ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
            
            # Set equal aspect ratio and reasonable limits
            max_loading = max(np.max(np.abs(pc_x_loadings)), np.max(np.abs(pc_y_loadings)))
            limit = min(1.2, max_loading * 1.2)
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)
            ax.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def create_biplot(self, pca_result: PCAResult, structure_ids: List[str],
                     outliers: Optional[List[StructureOutlier]] = None,
                     output_path: Path = None,
                     pc_x: int = 1, pc_y: int = 2,
                     max_features: int = 10) -> Path:
        """
        Create biplot combining scores and loadings
        
        Args:
            pca_result: PCA analysis results
            structure_ids: List of structure identifiers
            outliers: Optional list of structure outliers
            output_path: Path for saving the plot
            pc_x: Principal component for x-axis (1-indexed)
            pc_y: Principal component for y-axis (1-indexed)
            max_features: Maximum number of feature arrows to show
            
        Returns:
            Path to saved plot
        """
        fig, ax = plt.subplots(figsize=(14, 10))
        
        pc_x_idx = pc_x - 1
        pc_y_idx = pc_y - 1
        
        # Plot structure scores
        x_scores = pca_result.transformed_data[:, pc_x_idx]
        y_scores = pca_result.transformed_data[:, pc_y_idx]
        
        # Color by outlier status
        colors = ['blue'] * len(structure_ids)
        sizes = [60] * len(structure_ids)
        
        if outliers:
            outlier_dict = {o.structure_id: o for o in outliers}
            for i, struct_id in enumerate(structure_ids):
                if struct_id in outlier_dict:
                    outlier = outlier_dict[struct_id]
                    if outlier.outlier_type == 'extreme':
                        colors[i] = 'red'
                        sizes[i] = 100
                    elif outlier.outlier_type == 'moderate':
                        colors[i] = 'orange'
                        sizes[i] = 80
        
        scatter = ax.scatter(x_scores, y_scores, c=colors, s=sizes, 
                           alpha=0.7, edgecolors='black', label='Structures')
        
        # Add structure labels
        for i, struct_id in enumerate(structure_ids):
            ax.annotate(struct_id, (x_scores[i], y_scores[i]), 
                       xytext=(3, 3), textcoords='offset points',
                       fontsize=7, alpha=0.8)
        
        # Scale factor for loadings
        scale_factor = max(np.max(np.abs(x_scores)), np.max(np.abs(y_scores))) / \
                      max(np.max(np.abs(pca_result.loadings[:, pc_x_idx])), 
                          np.max(np.abs(pca_result.loadings[:, pc_y_idx])))
        
        # Select top contributing features
        if len(pca_result.feature_names) > max_features:
            total_contrib = (np.abs(pca_result.loadings[:, pc_x_idx]) + 
                           np.abs(pca_result.loadings[:, pc_y_idx]))
            top_indices = np.argsort(total_contrib)[-max_features:]
        else:
            top_indices = range(len(pca_result.feature_names))
        
        # Plot feature loadings as arrows
        for i in top_indices:
            arrow_x = pca_result.loadings[i, pc_x_idx] * scale_factor * 0.8
            arrow_y = pca_result.loadings[i, pc_y_idx] * scale_factor * 0.8
            
            ax.arrow(0, 0, arrow_x, arrow_y,
                    head_width=max(x_scores)/50, head_length=max(x_scores)/50,
                    fc='darkred', ec='darkred', alpha=0.6, linewidth=2)
            
            ax.text(arrow_x * 1.1, arrow_y * 1.1, pca_result.feature_names[i],
                   fontsize=8, ha='center', va='center', weight='bold',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))
        
        ax.set_xlabel(f'PC{pc_x} ({pca_result.explained_variance_ratio[pc_x_idx]:.1%} variance)')
        ax.set_ylabel(f'PC{pc_y} ({pca_result.explained_variance_ratio[pc_y_idx]:.1%} variance)')
        ax.set_title(f'PCA Biplot: PC{pc_x} vs PC{pc_y}')
        ax.grid(True, alpha=0.3)
        
        # Add legend
        if outliers:
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='blue', label='Normal Structure'),
                Patch(facecolor='orange', label='Moderate Outlier'),
                Patch(facecolor='red', label='Extreme Outlier'),
                Patch(facecolor='darkred', label='Feature Loading')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def create_outlier_heatmap(self, outliers: List[StructureOutlier],
                              pca_result: PCAResult,
                              output_path: Path,
                              top_n: int = 20) -> Path:
        """
        Create heatmap of outlier scores across principal components
        
        Args:
            outliers: List of structure outliers
            pca_result: PCA analysis results
            output_path: Path for saving the plot
            top_n: Number of top outliers to show
            
        Returns:
            Path to saved plot
        """
        # Select top outliers
        top_outliers = sorted(outliers, key=lambda x: x.outlier_score, reverse=True)[:top_n]
        
        # Create matrix of PC scores
        pc_scores_matrix = []
        outlier_labels = []
        
        for outlier in top_outliers:
            pc_scores_matrix.append(outlier.principal_component_scores)
            outlier_labels.append(f"{outlier.structure_id} ({outlier.outlier_score:.2f})")
        
        pc_scores_matrix = np.array(pc_scores_matrix)
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, max(6, len(top_outliers) * 0.4)))
        
        sns.heatmap(pc_scores_matrix, 
                   xticklabels=[f'PC{i+1}' for i in range(pc_scores_matrix.shape[1])],
                   yticklabels=outlier_labels,
                   cmap='RdBu_r', center=0, 
                   cbar_kws={'label': 'PC Score'},
                   ax=ax)
        
        ax.set_title('Outlier Principal Component Scores')
        ax.set_xlabel('Principal Component')
        ax.set_ylabel('Structure (Outlier Score)')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def create_pca_summary_plots(self, pca_result: PCAResult, 
                               structure_ids: List[str],
                               outliers: List[StructureOutlier],
                               output_dir: Path,
                               prefix: str = "pca") -> Dict[str, Path]:
        """
        Create a comprehensive set of PCA visualization plots
        
        Args:
            pca_result: PCA analysis results
            structure_ids: List of structure identifiers
            outliers: List of structure outliers
            output_dir: Directory for output plots
            prefix: Prefix for plot filenames
            
        Returns:
            Dictionary mapping plot type to file path
        """
        output_dir.mkdir(exist_ok=True)
        plot_files = {}
        
        # 1. Scree plot
        scree_path = output_dir / f"{prefix}_scree_plot.png"
        plot_files['scree'] = self.create_scree_plot(pca_result, scree_path)
        
        # 2. PC1 vs PC2 scatter
        if len(pca_result.transformed_data[0]) >= 2:
            scatter_path = output_dir / f"{prefix}_pc1_pc2_scatter.png"
            plot_files['scatter'] = self.create_pc_scatter_plot(
                pca_result, structure_ids, outliers, scatter_path
            )
        
        # 3. Loadings plot
        if len(pca_result.feature_names) > 0:
            loadings_path = output_dir / f"{prefix}_loadings.png"
            plot_files['loadings'] = self.create_loadings_plot(pca_result, loadings_path)
        
        # 4. Biplot
        if len(pca_result.transformed_data[0]) >= 2:
            biplot_path = output_dir / f"{prefix}_biplot.png"
            plot_files['biplot'] = self.create_biplot(
                pca_result, structure_ids, outliers, biplot_path
            )
        
        # 5. Outlier heatmap
        if len(outliers) > 1:
            heatmap_path = output_dir / f"{prefix}_outlier_heatmap.png"
            plot_files['outlier_heatmap'] = self.create_outlier_heatmap(
                outliers, pca_result, heatmap_path
            )
        
        return plot_files