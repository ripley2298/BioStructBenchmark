"""
biostructbenchmark/visualization/curves_plots.py

Visualization functions for CURVES+ analysis results
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap

from ..analysis.curves import CurvesParameters, HydrogenBond


class CurvesVisualizer:
    """Visualization tools for CURVES+ analysis results"""
    
    def __init__(self, style: str = 'seaborn-v0_8', figsize: Tuple[int, int] = (12, 8)):
        """
        Initialize visualizer
        
        Args:
            style: Matplotlib style
            figsize: Default figure size
        """
        plt.style.use(style)
        self.figsize = figsize
        self.setup_custom_colors()
    
    def setup_custom_colors(self):
        """Set up custom color schemes for DNA visualization"""
        # Base-specific colors
        self.base_colors = {
            'A': '#FF6B6B',  # Red
            'T': '#4ECDC4',  # Teal  
            'G': '#45B7D1',  # Blue
            'C': '#FFA07A',  # Orange
            'DA': '#FF6B6B', 'DT': '#4ECDC4', 'DG': '#45B7D1', 'DC': '#FFA07A'
        }
        
        # Parameter-specific colormaps
        self.param_cmaps = {
            'twist': 'RdYlBu_r',
            'roll': 'RdBu_r', 
            'tilt': 'PRGn',
            'rise': 'viridis',
            'slide': 'plasma',
            'shift': 'coolwarm'
        }
    
    def plot_parameter_comparison(self, comparison_df: pd.DataFrame, 
                                parameter: str, output_path: Optional[Path] = None) -> plt.Figure:
        """
        Plot comparison of a specific geometric parameter
        
        Args:
            comparison_df: DataFrame from CurvesAnalyzer.compare_geometries()
            parameter: Parameter to plot (e.g., 'twist_diff', 'roll_diff')
            output_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Bar plot of differences
        steps = comparison_df['step'].values
        differences = comparison_df[f'{parameter}_diff'].values
        
        # Color bars by magnitude
        colors = ['red' if x > 0 else 'blue' for x in differences]
        
        bars = ax1.bar(range(len(steps)), differences, color=colors, alpha=0.7)
        ax1.set_xlabel('Base Pair Step')
        ax1.set_ylabel(f'{parameter.replace("_", " ").title()} Difference (Predicted - Experimental)')
        ax1.set_title(f'{parameter.replace("_", " ").title()} Parameter Comparison')
        ax1.set_xticks(range(len(steps)))
        ax1.set_xticklabels(steps, rotation=45, ha='right')
        ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax1.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, differences)):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + np.sign(height)*0.1,
                    f'{value:.1f}', ha='center', va='bottom' if height > 0 else 'top',
                    fontsize=8)
        
        # Distribution plot
        ax2.hist(differences, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.set_xlabel(f'{parameter.replace("_", " ").title()} Difference')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution of Differences')
        ax2.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='No difference')
        ax2.axvline(x=np.mean(differences), color='orange', linestyle='-', 
                   label=f'Mean: {np.mean(differences):.2f}')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_parameter_heatmap(self, comparison_df: pd.DataFrame, 
                             output_path: Optional[Path] = None) -> plt.Figure:
        """
        Create heatmap of all parameter differences
        
        Args:
            comparison_df: DataFrame from CurvesAnalyzer.compare_geometries()
            output_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        # Select difference columns
        diff_cols = [col for col in comparison_df.columns if col.endswith('_diff')]
        
        # Create matrix for heatmap
        heatmap_data = comparison_df.set_index('step')[diff_cols]
        
        # Clean up column names for display
        clean_names = [col.replace('_diff', '').replace('_', ' ').title() for col in diff_cols]
        heatmap_data.columns = clean_names
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Use diverging colormap centered at zero
        vmax = max(abs(heatmap_data.min().min()), abs(heatmap_data.max().max()))
        
        sns.heatmap(heatmap_data.T, 
                   cmap='RdBu_r', 
                   center=0, 
                   vmin=-vmax, 
                   vmax=vmax,
                   annot=True, 
                   fmt='.1f',
                   cbar_kws={'label': 'Parameter Difference (Predicted - Experimental)'},
                   ax=ax)
        
        ax.set_title('DNA Geometric Parameter Differences Heatmap', fontsize=14, fontweight='bold')
        ax.set_xlabel('Base Pair Step', fontsize=12)
        ax.set_ylabel('Geometric Parameters', fontsize=12)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_groove_geometry(self, experimental_params: Dict[str, CurvesParameters],
                           predicted_params: Dict[str, CurvesParameters],
                           output_path: Optional[Path] = None) -> plt.Figure:
        """
        Plot major and minor groove geometry comparison
        
        Args:
            experimental_params: Experimental CURVES+ parameters
            predicted_params: Predicted CURVES+ parameters  
            output_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Collect groove data
        steps = []
        exp_major_width, pred_major_width = [], []
        exp_minor_width, pred_minor_width = [], []
        exp_major_depth, pred_major_depth = [], []
        exp_minor_depth, pred_minor_depth = [], []
        
        common_steps = set(experimental_params.keys()) & set(predicted_params.keys())
        
        for step in sorted(common_steps):
            exp_param = experimental_params[step]
            pred_param = predicted_params[step]
            
            if (exp_param.major_groove_width is not None and 
                pred_param.major_groove_width is not None):
                steps.append(step)
                exp_major_width.append(exp_param.major_groove_width)
                pred_major_width.append(pred_param.major_groove_width)
                exp_minor_width.append(exp_param.minor_groove_width or 0)
                pred_minor_width.append(pred_param.minor_groove_width or 0)
                exp_major_depth.append(exp_param.major_groove_depth or 0)
                pred_major_depth.append(pred_param.major_groove_depth or 0)
                exp_minor_depth.append(exp_param.minor_groove_depth or 0)
                pred_minor_depth.append(pred_param.minor_groove_depth or 0)
        
        x_pos = np.arange(len(steps))
        width = 0.35
        
        # Major groove width
        ax1.bar(x_pos - width/2, exp_major_width, width, label='Experimental', 
                color='skyblue', alpha=0.8)
        ax1.bar(x_pos + width/2, pred_major_width, width, label='Predicted', 
                color='lightcoral', alpha=0.8)
        ax1.set_xlabel('Base Pair Step')
        ax1.set_ylabel('Major Groove Width (Å)')
        ax1.set_title('Major Groove Width Comparison')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(steps, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Minor groove width
        ax2.bar(x_pos - width/2, exp_minor_width, width, label='Experimental', 
                color='skyblue', alpha=0.8)
        ax2.bar(x_pos + width/2, pred_minor_width, width, label='Predicted', 
                color='lightcoral', alpha=0.8)
        ax2.set_xlabel('Base Pair Step')
        ax2.set_ylabel('Minor Groove Width (Å)')
        ax2.set_title('Minor Groove Width Comparison')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(steps, rotation=45, ha='right')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Correlation plots
        if exp_major_width and pred_major_width:
            ax3.scatter(exp_major_width, pred_major_width, alpha=0.7, color='blue')
            # Add perfect correlation line
            min_val, max_val = min(min(exp_major_width), min(pred_major_width)), \
                             max(max(exp_major_width), max(pred_major_width))
            ax3.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7, label='Perfect correlation')
            ax3.set_xlabel('Experimental Major Groove Width (Å)')
            ax3.set_ylabel('Predicted Major Groove Width (Å)')
            ax3.set_title('Major Groove Width Correlation')
            
            # Calculate and display correlation
            corr = np.corrcoef(exp_major_width, pred_major_width)[0, 1]
            ax3.text(0.05, 0.95, f'R = {corr:.3f}', transform=ax3.transAxes, 
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            ax3.legend()
            ax3.grid(True, alpha=0.3)
        
        if exp_minor_width and pred_minor_width:
            ax4.scatter(exp_minor_width, pred_minor_width, alpha=0.7, color='orange')
            min_val, max_val = min(min(exp_minor_width), min(pred_minor_width)), \
                             max(max(exp_minor_width), max(pred_minor_width))
            ax4.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7, label='Perfect correlation')
            ax4.set_xlabel('Experimental Minor Groove Width (Å)')
            ax4.set_ylabel('Predicted Minor Groove Width (Å)')
            ax4.set_title('Minor Groove Width Correlation')
            
            corr = np.corrcoef(exp_minor_width, pred_minor_width)[0, 1]
            ax4.text(0.05, 0.95, f'R = {corr:.3f}', transform=ax4.transAxes,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_hydrogen_bond_network(self, hbonds: List[HydrogenBond], 
                                 structure_name: str = "Structure",
                                 output_path: Optional[Path] = None) -> plt.Figure:
        """
        Create network visualization of hydrogen bonds
        
        Args:
            hbonds: List of HydrogenBond objects
            structure_name: Name for plot title
            output_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        if not hbonds:
            fig, ax = plt.subplots(figsize=self.figsize)
            ax.text(0.5, 0.5, 'No hydrogen bonds detected', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title(f'{structure_name} - Hydrogen Bond Network')
            return fig
        
        # Create network graph
        G = nx.Graph()
        
        # Add nodes (residues)
        residues = set()
        for hb in hbonds:
            residues.add(hb.donor_residue)
            residues.add(hb.acceptor_residue)
        
        for residue in residues:
            # Determine if protein or DNA
            res_type = 'protein' if any(aa in residue for aa in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']) else 'DNA'
            G.add_node(residue, type=res_type)
        
        # Add edges (hydrogen bonds)
        for hb in hbonds:
            G.add_edge(hb.donor_residue, hb.acceptor_residue, 
                      weight=1.0/hb.distance,  # Stronger bonds have higher weight
                      distance=hb.distance,
                      strength=hb.strength,
                      atoms=f"{hb.donor_atom}-{hb.acceptor_atom}")
        
        # Create layout
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Use spring layout for better visualization
        pos = nx.spring_layout(G, k=3, iterations=50)
        
        # Separate protein and DNA nodes
        protein_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'protein']
        dna_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'DNA']
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, 
                              node_color='lightblue', node_size=500, 
                              alpha=0.8, label='Protein', ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=dna_nodes, 
                              node_color='lightcoral', node_size=500, 
                              alpha=0.8, label='DNA', ax=ax)
        
        # Draw edges with different styles for bond strength
        strong_edges = [(u, v) for u, v, d in G.edges(data=True) if d['strength'] == 'strong']
        medium_edges = [(u, v) for u, v, d in G.edges(data=True) if d['strength'] == 'medium']
        weak_edges = [(u, v) for u, v, d in G.edges(data=True) if d['strength'] == 'weak']
        
        nx.draw_networkx_edges(G, pos, edgelist=strong_edges, width=3, 
                              alpha=0.8, edge_color='red', label='Strong H-bonds', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=medium_edges, width=2, 
                              alpha=0.6, edge_color='orange', label='Medium H-bonds', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=weak_edges, width=1, 
                              alpha=0.4, edge_color='gray', label='Weak H-bonds', ax=ax)
        
        # Add labels
        nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)
        
        ax.set_title(f'{structure_name} - Hydrogen Bond Network\n'
                    f'({len(hbonds)} H-bonds, {len(protein_nodes)} protein residues, {len(dna_nodes)} DNA residues)', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper right')
        ax.axis('off')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_hbond_statistics(self, exp_hbonds: List[HydrogenBond], 
                            pred_hbonds: List[HydrogenBond],
                            output_path: Optional[Path] = None) -> plt.Figure:
        """
        Compare hydrogen bond statistics between experimental and predicted structures
        
        Args:
            exp_hbonds: Experimental hydrogen bonds
            pred_hbonds: Predicted hydrogen bonds
            output_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Distance distributions
        exp_distances = [hb.distance for hb in exp_hbonds]
        pred_distances = [hb.distance for hb in pred_hbonds]
        
        ax1.hist(exp_distances, bins=20, alpha=0.7, label='Experimental', color='blue')
        ax1.hist(pred_distances, bins=20, alpha=0.7, label='Predicted', color='red')
        ax1.set_xlabel('H-bond Distance (Å)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('H-bond Distance Distribution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Strength distribution
        exp_strengths = [hb.strength for hb in exp_hbonds]
        pred_strengths = [hb.strength for hb in pred_hbonds]
        
        strength_categories = ['strong', 'medium', 'weak']
        exp_counts = [exp_strengths.count(s) for s in strength_categories]
        pred_counts = [pred_strengths.count(s) for s in strength_categories]
        
        x = np.arange(len(strength_categories))
        width = 0.35
        
        ax2.bar(x - width/2, exp_counts, width, label='Experimental', alpha=0.8, color='blue')
        ax2.bar(x + width/2, pred_counts, width, label='Predicted', alpha=0.8, color='red')
        ax2.set_xlabel('H-bond Strength')
        ax2.set_ylabel('Count')
        ax2.set_title('H-bond Strength Distribution')
        ax2.set_xticks(x)
        ax2.set_xticklabels(strength_categories)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Residue type analysis
        def get_residue_types(hbonds):
            protein_count = sum(1 for hb in hbonds if any(aa in hb.donor_residue for aa in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']) or
                                                      any(aa in hb.acceptor_residue for aa in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']))
            return protein_count
        
        exp_protein_hbonds = get_residue_types(exp_hbonds)
        pred_protein_hbonds = get_residue_types(pred_hbonds)
        
        categories = ['Total H-bonds', 'Protein-DNA H-bonds']
        exp_values = [len(exp_hbonds), exp_protein_hbonds]
        pred_values = [len(pred_hbonds), pred_protein_hbonds]
        
        x = np.arange(len(categories))
        
        ax3.bar(x - width/2, exp_values, width, label='Experimental', alpha=0.8, color='blue')
        ax3.bar(x + width/2, pred_values, width, label='Predicted', alpha=0.8, color='red')
        ax3.set_xlabel('H-bond Type')
        ax3.set_ylabel('Count')
        ax3.set_title('H-bond Count Comparison')
        ax3.set_xticks(x)
        ax3.set_xticklabels(categories)
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Summary statistics table
        ax4.axis('tight')
        ax4.axis('off')
        
        summary_data = [
            ['Metric', 'Experimental', 'Predicted', 'Difference'],
            ['Total H-bonds', len(exp_hbonds), len(pred_hbonds), len(pred_hbonds) - len(exp_hbonds)],
            ['Mean distance (Å)', f'{np.mean(exp_distances):.2f}' if exp_distances else 'N/A', 
             f'{np.mean(pred_distances):.2f}' if pred_distances else 'N/A', 
             f'{np.mean(pred_distances) - np.mean(exp_distances):.2f}' if exp_distances and pred_distances else 'N/A'],
            ['Strong bonds', exp_counts[0], pred_counts[0], pred_counts[0] - exp_counts[0]],
            ['Medium bonds', exp_counts[1], pred_counts[1], pred_counts[1] - exp_counts[1]],
            ['Weak bonds', exp_counts[2], pred_counts[2], pred_counts[2] - exp_counts[2]]
        ]
        
        table = ax4.table(cellText=summary_data[1:], colLabels=summary_data[0],
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        ax4.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def create_pymol_script(self, hbonds: List[HydrogenBond], 
                          structure_path: Path, 
                          output_path: Path) -> None:
        """
        Generate PyMOL script to visualize hydrogen bonds
        
        Args:
            hbonds: List of hydrogen bonds
            structure_path: Path to structure file
            output_path: Path to save PyMOL script
        """
        script_lines = [
            "# PyMOL script for hydrogen bond visualization",
            f"# Generated for structure: {structure_path.name}",
            "",
            f"load {structure_path.absolute()}",
            "",
            "# Basic structure visualization",
            "hide everything",
            "show cartoon, polymer",
            "show sticks, nucleic",
            "color skyblue, polymer and chain A",
            "color lightpink, nucleic",
            "",
            "# Color DNA bases",
            "color red, resn DA+A",
            "color blue, resn DT+T", 
            "color green, resn DG+G",
            "color orange, resn DC+C",
            "",
            "# Hydrogen bond visualization",
            "distance hbonds, none, none",
            ""
        ]
        
        for i, hb in enumerate(hbonds):
            # Parse residue information
            donor_parts = hb.donor_residue.split('_')
            acceptor_parts = hb.acceptor_residue.split('_')
            
            if len(donor_parts) >= 2 and len(acceptor_parts) >= 2:
                donor_res = donor_parts[1]
                acceptor_res = acceptor_parts[1]
                
                # Create distance measurement
                script_lines.append(
                    f"distance hbond_{i}, "
                    f"resi {donor_res} and name {hb.donor_atom}, "
                    f"resi {acceptor_res} and name {hb.acceptor_atom}"
                )
        
        script_lines.extend([
            "",
            "# Styling for hydrogen bonds",
            "hide labels, hbonds",
            "color yellow, hbonds",
            "",
            "# Set view",
            "orient",
            "zoom",
            "",
            "# Save image",
            f"png {output_path.parent / 'hbond_visualization.png'}, dpi=300"
        ])
        
        with open(output_path, 'w') as f:
            f.write('\n'.join(script_lines))


def create_curves_report(analysis_dir: Path, output_dir: Path) -> None:
    """
    Generate comprehensive CURVES+ analysis report
    
    Args:
        analysis_dir: Directory containing CURVES+ analysis results
        output_dir: Directory for report output
    """
    visualizer = CurvesVisualizer()
    output_dir.mkdir(exist_ok=True)
    
    # Load comparison data
    comparison_file = analysis_dir / 'all_geometry_comparisons.csv'
    if not comparison_file.exists():
        print(f"No comparison data found at {comparison_file}")
        return
    
    comparison_df = pd.read_csv(comparison_file)
    
    # Create parameter-specific plots
    parameters = ['twist', 'roll', 'tilt', 'rise', 'slide', 'shift',
                 'shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']
    
    for param in parameters:
        param_col = f'{param}_diff'
        if param_col in comparison_df.columns:
            fig = visualizer.plot_parameter_comparison(
                comparison_df, param, 
                output_dir / f'{param}_comparison.png'
            )
            plt.close(fig)
    
    # Create heatmap
    heatmap_fig = visualizer.plot_parameter_heatmap(
        comparison_df, 
        output_dir / 'parameter_heatmap.png'
    )
    plt.close(heatmap_fig)
    
    # Process individual structure hydrogen bond data
    structure_dirs = [d for d in analysis_dir.iterdir() if d.is_dir()]
    
    for struct_dir in structure_dirs:
        exp_hbonds_file = struct_dir / 'experimental_hbonds.csv'
        pred_hbonds_file = struct_dir / 'predicted_hbonds.csv'
        
        if exp_hbonds_file.exists() and pred_hbonds_file.exists():
            # Load hydrogen bond data
            exp_hbonds_df = pd.read_csv(exp_hbonds_file)
            pred_hbonds_df = pd.read_csv(pred_hbonds_file)
            
            # Convert to HydrogenBond objects
            exp_hbonds = [HydrogenBond(**row) for _, row in exp_hbonds_df.iterrows()]
            pred_hbonds = [HydrogenBond(**row) for _, row in pred_hbonds_df.iterrows()]
            
            # Create visualizations
            struct_output = output_dir / struct_dir.name
            struct_output.mkdir(exist_ok=True)
            
            # Hydrogen bond statistics
            stats_fig = visualizer.plot_hbond_statistics(
                exp_hbonds, pred_hbonds,
                struct_output / 'hbond_statistics.png'
            )
            plt.close(stats_fig)
            
            # Network visualizations
            exp_network_fig = visualizer.plot_hydrogen_bond_network(
                exp_hbonds, f"{struct_dir.name} (Experimental)",
                struct_output / 'experimental_hbond_network.png'
            )
            plt.close(exp_network_fig)
            
            pred_network_fig = visualizer.plot_hydrogen_bond_network(
                pred_hbonds, f"{struct_dir.name} (Predicted)",
                struct_output / 'predicted_hbond_network.png'
            )
            plt.close(pred_network_fig)
    
    print(f"CURVES+ analysis report generated in {output_dir}")


if __name__ == "__main__":
    # Example usage
    from pathlib import Path
    
    # Generate example report
    analysis_dir = Path("results/curves_analysis")
    output_dir = Path("results/curves_visualizations")
    
    if analysis_dir.exists():
        create_curves_report(analysis_dir, output_dir)
    else:
        print(f"Analysis directory {analysis_dir} not found")
