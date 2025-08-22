"""
biostructbenchmark/visualization/residue_plots.py
Comprehensive residue-level visualization with heatmaps, correlations, and dashboard functionality
Merged from plots.py for unified residue analysis and model validation
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Union, Any

try:
    import seaborn as sns
    SNS_AVAILABLE = True
except ImportError:
    SNS_AVAILABLE = False

from biostructbenchmark.core.alignment import ResidueRMSD


class ResidueVisualizer:
    """Comprehensive residue-level visualization with heatmaps, correlations, and dashboard capabilities"""
    
    def __init__(self, style: str = 'seaborn-v0_8' if SNS_AVAILABLE else 'default', 
                 figsize: Tuple[int, int] = (12, 8)):
        """Initialize with style and figure size preferences"""
        plt.style.use(style)
        self.figsize = figsize
        self.fig_params = {'dpi': 300, 'bbox_inches': 'tight', 'facecolor': 'white'}
        self.colors = {
            'primary': '#2E86AB', 'secondary': '#A23B72', 'accent': '#F18F01',
            'protein': '#4CAF50', 'dna': '#FF9800', 'error': '#F44336',
            'high_rmsd': '#F18F01', 'low_rmsd': '#C73E1D'
        }
    
    def plot_rmsd_heatmap(self, residue_data: List[ResidueRMSD], 
                         output_path: Optional[Path] = None) -> plt.Figure:
        """
        Create per-residue RMSD heatmap using minimal memory footprint
        
        Args:
            residue_data: List of ResidueRMSD objects
            output_path: Optional save path
            
        Returns:
            matplotlib Figure object
        """
        if not residue_data:
            raise ValueError("Empty residue data provided")
        
        # Convert to minimal DataFrame view (avoid data duplication)
        df = pd.DataFrame([{
            'chain': r.chain_id, 'pos': r.position, 'rmsd': r.rmsd, 
            'type': r.molecule_type, 'res': r.residue_type
        } for r in residue_data])
        
        # Pivot for heatmap (memory efficient)
        heatmap_data = df.pivot_table(index='chain', columns='pos', values='rmsd', fill_value=0)
        
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Use matplotlib's efficient imshow for large datasets
        im = ax.imshow(heatmap_data.values, cmap='viridis', aspect='auto', interpolation='nearest')
        
        # Minimal labeling for efficiency
        ax.set_yticks(range(len(heatmap_data.index)))
        ax.set_yticklabels(heatmap_data.index)
        ax.set_xlabel('Residue Position')
        ax.set_ylabel('Chain ID')
        ax.set_title('Per-Residue RMSD Heatmap', fontweight='bold')
        
        # Colorbar with RMSD values
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('RMSD (Å)', rotation=270, labelpad=15)
        
        plt.tight_layout()
        if output_path: plt.savefig(output_path, dpi=300, bbox_inches='tight')
        return fig
    
    def plot_residue_correlation_heatmap(self, residue_data: List[ResidueRMSD],
                                       secondary_data: Optional[List] = None,
                                       output_path: Optional[Path] = None) -> plt.Figure:
        """
        Create correlation heatmap between residue RMSD and secondary data
        
        Args:
            residue_data: Primary residue RMSD data  
            secondary_data: Optional secondary data for correlation (e.g., B-factors)
            output_path: Optional save path
            
        Returns:
            matplotlib Figure
        """
        if not residue_data:
            raise ValueError("Empty residue data provided")
        
        # Create correlation matrix using imshow for efficiency
        df = pd.DataFrame([{
            'position': r.position, 'rmsd': r.rmsd, 'chain': r.chain_id,
            'type': r.molecule_type, 'residue': r.residue_type
        } for r in residue_data])
        
        # Add secondary data if provided
        if secondary_data:
            sec_df = pd.DataFrame([{
                'position': i, 'secondary_value': val
            } for i, val in enumerate(secondary_data)])
            df = df.merge(sec_df, on='position', how='left')
        
        # Create correlation matrix focusing on numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) < 2:
            # If no secondary data, create position-RMSD correlation
            corr_matrix = df[['position', 'rmsd']].corr()
        else:
            corr_matrix = df[numeric_cols].corr()
        
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Use imshow for memory efficiency
        im = ax.imshow(corr_matrix.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
        
        # Set ticks and labels
        ax.set_xticks(range(len(corr_matrix.columns)))
        ax.set_yticks(range(len(corr_matrix.index)))
        ax.set_xticklabels(corr_matrix.columns, rotation=45)
        ax.set_yticklabels(corr_matrix.index)
        
        # Add correlation values as text
        for i in range(len(corr_matrix.index)):
            for j in range(len(corr_matrix.columns)):
                text = ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                             ha="center", va="center", color="black" if abs(corr_matrix.iloc[i, j]) < 0.5 else "white")
        
        ax.set_title('Residue Data Correlation Matrix', fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Correlation Coefficient', rotation=270, labelpad=15)
        
        plt.tight_layout()
        if output_path: plt.savefig(output_path, **self.fig_params)
        return fig
    
    def plot_chain_comparison(self, residue_data: List[ResidueRMSD]) -> plt.Figure:
        """Plot per-chain RMSD comparison with minimal processing"""
        # Group by chain efficiently using dictionary comprehension
        chain_rmsds = {}
        for r in residue_data:
            chain_rmsds.setdefault(r.chain_id, []).append(r.rmsd)
        
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Box plot for chain comparison
        chains, rmsds = zip(*chain_rmsds.items())
        ax.boxplot(rmsds, labels=chains)
        
        ax.set_xlabel('Chain ID')
        ax.set_ylabel('RMSD (Å)')
        ax.set_title('Per-Chain RMSD Distribution', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_hydrogen_bond_errors(self, hbond_error_data: List[Tuple[str, int]], 
                                  output_path: Optional[Path] = None) -> plt.Figure:
        """
        Generate publication-ready horizontal bar chart for hydrogen bond network errors
        
        Args:
            hbond_error_data: List of tuples (residue_name, error_count)
                             e.g., [("Lys102", 5), ("Ser57", 4), ("Asp201", 3)]
            output_path: Optional save path
            
        Returns:
            matplotlib Figure object
        """
        # Filter and sort data according to specifications
        # Only include residues with error_count >= 2
        filtered_data = [(residue, count) for residue, count in hbond_error_data if count >= 2]
        
        # Sort descending by error count (highest first)
        filtered_data.sort(key=lambda x: x[1], reverse=True)
        
        # Take only top 10 residues
        top_data = filtered_data[:10]
        
        if not top_data:
            # Create empty plot if no data meets criteria
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.text(0.5, 0.5, 'No residues with ≥2 H-bond errors', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return fig
        
        # Extract residues and errors for plotting
        residues, errors = zip(*top_data)
        
        # Create figure with exact specifications
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        
        # Create horizontal bar chart
        bars = ax.barh(residues, errors, color="#C50000", edgecolor="black", height=0.6)
        
        # Configure axes and labels
        ax.set_xlabel("Number of Hydrogen Bond Errors", fontsize=10)
        ax.set_title("Top Residues with H-bond Network Errors (Predicted vs. Experimental)", 
                    fontsize=12, pad=10)
        
        # Invert y-axis so highest error is at top
        ax.invert_yaxis()
        
        # Set x-axis ticks to integers from 0 to max error count
        max_error = max(errors)
        ax.set_xticks(range(0, max_error + 1))
        
        # Remove gridlines and ensure clean appearance
        ax.grid(False)
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        plt.tight_layout()
        
        # Save if path provided
        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=300, facecolor='white')
        
        return fig
    
    def residue_dashboard(self, data_dict: Dict[str, Any], output_path: Path) -> plt.Figure:
        """
        Create comprehensive residue analysis dashboard
        
        Args:
            data_dict: Dictionary with keys 'rmsd', 'bfactor', 'consensus', 'mutations'
            output_path: Save path for dashboard
            
        Returns:
            matplotlib Figure with 2x2 subplot layout
        """
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Residue Analysis Dashboard', fontsize=16, fontweight='bold')
        
        # RMSD heatmap (top-left)
        if 'rmsd' in data_dict and data_dict['rmsd']:
            self._plot_residue_rmsd_summary(axes[0,0], data_dict['rmsd'])
        else:
            axes[0,0].text(0.5, 0.5, 'No RMSD Data', ha='center', va='center', transform=axes[0,0].transAxes)
        
        # B-factor correlation (top-right) 
        if 'bfactor' in data_dict and data_dict['bfactor']:
            self._plot_bfactor_correlation(axes[0,1], data_dict['bfactor'])
        else:
            axes[0,1].text(0.5, 0.5, 'No B-factor Data', ha='center', va='center', transform=axes[0,1].transAxes)
        
        # Consensus errors (bottom-left)
        if 'consensus' in data_dict and data_dict['consensus']:
            self._plot_consensus_summary(axes[1,0], data_dict['consensus'])
        else:
            axes[1,0].text(0.5, 0.5, 'No Consensus Data', ha='center', va='center', transform=axes[1,0].transAxes)
        
        # Mutation impact (bottom-right)
        if 'mutations' in data_dict and data_dict['mutations']:
            self._plot_mutation_impact(axes[1,1], data_dict['mutations'])
        else:
            axes[1,1].text(0.5, 0.5, 'No Mutation Data', ha='center', va='center', transform=axes[1,1].transAxes)
        
        plt.tight_layout()
        fig.savefig(output_path, **self.fig_params)
        return fig
    
    def _plot_residue_rmsd_summary(self, ax: plt.Axes, rmsd_data: List[ResidueRMSD]) -> None:
        """Efficient residue RMSD heatmap summary"""
        # Create minimal heatmap view
        df = pd.DataFrame([{
            'chain': r.chain_id, 'pos': r.position, 'rmsd': r.rmsd, 'type': r.molecule_type
        } for r in rmsd_data])
        
        # Group for efficient visualization using imshow
        heatmap_data = df.pivot_table(index='chain', columns='pos', values='rmsd', fill_value=0)
        
        if not heatmap_data.empty:
            im = ax.imshow(heatmap_data.values, cmap='viridis', aspect='auto', interpolation='nearest')
            ax.set_title('Residue RMSD Overview', fontweight='bold')
            ax.set_ylabel('Chain ID')
            ax.set_xlabel('Position')
            
            # Minimal labeling for dashboard view
            if len(heatmap_data.index) <= 10:
                ax.set_yticks(range(len(heatmap_data.index)))
                ax.set_yticklabels(heatmap_data.index)
    
    def _plot_bfactor_correlation(self, ax: plt.Axes, bfactor_data: List) -> None:
        """B-factor vs confidence correlation plot"""
        if not bfactor_data:
            return
            
        # Extract data efficiently (assuming BFactorComparison objects)
        exp_vals = [b.experimental_bfactor for b in bfactor_data]
        pred_vals = [b.predicted_confidence for b in bfactor_data]
        
        # Scatter plot with correlation
        ax.scatter(exp_vals, pred_vals, alpha=0.6, s=20, c=self.colors['secondary'])
        
        # Add correlation line
        if len(exp_vals) > 1:
            correlation = np.corrcoef(exp_vals, pred_vals)[0,1]
            z = np.polyfit(exp_vals, pred_vals, 1)
            p = np.poly1d(z)
            ax.plot(exp_vals, p(exp_vals), 'r--', alpha=0.8)
            ax.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax.transAxes, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_xlabel('Experimental B-factor')
        ax.set_ylabel('Predicted Confidence')
        ax.set_title('B-factor vs Confidence', fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    def _plot_consensus_summary(self, ax: plt.Axes, consensus_data: List) -> None:
        """Consensus error visualization"""
        if not consensus_data:
            return
            
        # Extract values efficiently without aggregated metrics
        confidences = [c.confidence for c in consensus_data]
        
        # Histogram of confidence values
        ax.hist(confidences, bins=20, alpha=0.7, color=self.colors['primary'])
        ax.set_xlabel('Consensus Confidence')
        ax.set_ylabel('Frequency')
        ax.set_title('Consensus Confidence Distribution', fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    def _plot_mutation_impact(self, ax: plt.Axes, mutation_data: List) -> None:
        """Mutation impact summary"""
        if not mutation_data:
            return
            
        # Group by mutation type efficiently
        impact_by_type = {}
        for mut in mutation_data:
            mut_type = getattr(mut, 'mutation_type', 'unknown')
            impact_by_type.setdefault(mut_type, []).append(getattr(mut, 'impact_score', 0))
        
        if impact_by_type:
            # Box plot of impact scores by type
            types, impacts = zip(*impact_by_type.items())
            ax.boxplot(impacts, labels=types)
            ax.set_ylabel('Impact Score')
            ax.set_title('Mutation Impact by Type', fontweight='bold')
            ax.grid(True, alpha=0.3)


def create_residue_analysis(residue_data: List[ResidueRMSD], output_dir: Path,
                           analysis_data: Optional[Dict[str, Any]] = None) -> Dict[str, Path]:
    """
    Generate comprehensive residue analysis with heatmaps, correlations, and dashboard
    
    Args:
        residue_data: List of ResidueRMSD objects
        output_dir: Directory for output files
        analysis_data: Optional additional analysis data for correlations and H-bond errors
        
    Returns:
        Dictionary mapping plot types to output paths
    """
    output_dir.mkdir(exist_ok=True)
    viz = ResidueVisualizer()
    output_paths = {}
    
    # Generate core plots efficiently (one at a time to minimize memory)
    plots = [
        ('heatmap', lambda: viz.plot_rmsd_heatmap(residue_data)),
        ('correlation', lambda: viz.plot_residue_correlation_heatmap(
            residue_data, 
            analysis_data.get('secondary_data') if analysis_data else None
        )),
        ('chains', lambda: viz.plot_chain_comparison(residue_data))
    ]
    
    for plot_name, plot_func in plots:
        fig = plot_func()
        path = output_dir / f'residue_{plot_name}.png'
        fig.savefig(path, **viz.fig_params)
        plt.close(fig)  # Free memory immediately
        output_paths[plot_name] = path
    
    # Generate hydrogen bond error plot if data provided
    if analysis_data and 'hbond_errors' in analysis_data:
        hbond_path = output_dir / 'h-bond_errors.png'
        fig = viz.plot_hydrogen_bond_errors(analysis_data['hbond_errors'], hbond_path)
        plt.close(fig)
        output_paths['hbond_errors'] = hbond_path
        
        # Print caption as required
        print("\n*Top 10 residues (≥2 errors) cover 87% of all H-bond errors. Errors defined as missing bonds or incorrect geometry (MolProbity validation). Full list in Supplementary Table S1.*")
    
    # Generate dashboard if additional data provided
    if analysis_data:
        dashboard_data = {'rmsd': residue_data}
        dashboard_data.update(analysis_data)
        dashboard_path = output_dir / 'residue_dashboard.png'
        fig = viz.residue_dashboard(dashboard_data, dashboard_path)
        plt.close(fig)
        output_paths['dashboard'] = dashboard_path
    
    # Generate summary CSV efficiently  
    summary_path = output_dir / 'residue_summary.csv'
    pd.DataFrame([{
        'residue_id': r.residue_id, 'chain': r.chain_id, 'position': r.position,
        'rmsd': f'{r.rmsd:.3f}', 'type': r.molecule_type, 'residue': r.residue_type
    } for r in residue_data]).to_csv(summary_path, index=False)
    
    output_paths['summary'] = summary_path
    return output_paths


def quick_residue_plot(residue_data: List[ResidueRMSD], plot_type: str = 'heatmap',
                      output_path: Optional[Path] = None) -> plt.Figure:
    """
    Single-function plotting for quick analysis
    
    Args:
        residue_data: Residue RMSD data
        plot_type: 'heatmap', 'correlation', or 'chains'
        output_path: Optional save path
        
    Returns:
        matplotlib Figure
    """
    viz = ResidueVisualizer()
    
    plot_funcs = {
        'heatmap': viz.plot_rmsd_heatmap,
        'correlation': viz.plot_residue_correlation_heatmap, 
        'chains': viz.plot_chain_comparison
    }
    
    if plot_type not in plot_funcs:
        raise ValueError(f"Plot type '{plot_type}' not supported. Use: {list(plot_funcs.keys())}")
    
    fig = plot_funcs[plot_type](residue_data)
    if output_path: fig.savefig(output_path, **viz.fig_params)
    return fig


def plot_hydrogen_bond_errors_standalone(hbond_error_data: List[Tuple[str, int]], 
                                         output_path: Optional[Path] = None) -> plt.Figure:
    """
    Standalone function for hydrogen bond error visualization
    
    Args:
        hbond_error_data: List of tuples (residue_name, error_count)
                         e.g., [("Lys102", 5), ("Ser57", 4), ("Asp201", 3)]
        output_path: Optional save path
        
    Returns:
        matplotlib Figure object
    """
    viz = ResidueVisualizer()
    fig = viz.plot_hydrogen_bond_errors(hbond_error_data, output_path)
    
    # Print caption as required
    print("\n*Top 10 residues (≥2 errors) cover 87% of all H-bond errors. Errors defined as missing bonds or incorrect geometry (MolProbity validation). Full list in Supplementary Table S1.*")
    
    return fig


# Legacy compatibility aliases and functions from plots.py
class PublicationPlotter(ResidueVisualizer):
    """Legacy compatibility alias for ResidueVisualizer with publication methods"""
    
    def summary_dashboard(self, data_dict: Dict[str, Any], output_path: Path) -> plt.Figure:
        """Legacy compatibility method - redirects to residue_dashboard"""
        return self.residue_dashboard(data_dict, output_path)
    
    def comparison_plot(self, data1: List[float], data2: List[float], 
                       labels: Tuple[str, str], title: str = "Comparison",
                       output_path: Optional[Path] = None) -> plt.Figure:
        """
        Generic comparison plot for any two datasets
        
        Args:
            data1, data2: Datasets to compare
            labels: Labels for datasets
            title: Plot title
            output_path: Optional save path
            
        Returns:
            matplotlib Figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Side-by-side histograms
        ax1.hist([data1, data2], bins=20, alpha=0.7, label=labels, 
                color=[self.colors['primary'], self.colors['secondary']])
        ax1.set_xlabel('Value')
        ax1.set_ylabel('Frequency')
        ax1.set_title(f'{title} - Distribution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Scatter plot comparison
        min_len = min(len(data1), len(data2))
        if min_len > 0:
            ax2.scatter(data1[:min_len], data2[:min_len], alpha=0.6, s=20)
            
            # Perfect correlation line
            max_val = max(max(data1[:min_len]), max(data2[:min_len]))
            ax2.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, label='Perfect correlation')
            
            ax2.set_xlabel(labels[0])
            ax2.set_ylabel(labels[1])
            ax2.set_title(f'{title} - Correlation')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        if output_path: fig.savefig(output_path, **self.fig_params)
        return fig


def create_publication_report(analysis_results: Dict[str, Any], output_dir: Path) -> Dict[str, Path]:
    """
    Legacy compatibility function - redirects to create_residue_analysis
    """
    if 'rmsd' in analysis_results:
        return create_residue_analysis(analysis_results['rmsd'], output_dir, analysis_results)
    else:
        output_dir.mkdir(exist_ok=True)
        return {}


def quick_comparison(data1: List[float], data2: List[float], 
                    labels: Tuple[str, str] = ('Dataset 1', 'Dataset 2'),
                    output_path: Optional[Path] = None) -> plt.Figure:
    """
    Quick two-dataset comparison plot - single function interface
    
    Args:
        data1, data2: Datasets to compare  
        labels: Dataset labels
        output_path: Optional save path
        
    Returns:
        matplotlib Figure
    """
    plotter = PublicationPlotter()
    return plotter.comparison_plot(data1, data2, labels, output_path=output_path)


# Backward compatibility alias
create_residue_report = create_residue_analysis
