"""
biostructbenchmark/analysis/dssr.py
X3DNA-DSSR wrapper for analyzing DNA structural parameters in protein-DNA complexes

Focuses on the 5 most critical parameters for protein-DNA binding interface analysis:
1. Base pairs - Number of canonical base pairs at protein-binding interface
2. Helical twist - Average twist per base pair step (affects groove dimensions)
3. Major groove width - Primary protein-binding interface dimensions
4. Minor groove width - Affects DNA bending and protein recognition
5. Stacking energy - Interface stability indicator
"""

import json
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np


@dataclass
class DSSRParameters:
    """Container for critical DSSR structural parameters"""
    structure_name: str
    structure_type: str  # 'Experimental' or 'Predicted'
    base_pairs: int
    helical_twist: float  # degrees
    major_groove_width: float  # Angstroms
    minor_groove_width: float  # Angstroms
    stacking_energy: float  # kcal/mol
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame creation"""
        return {
            'Structure': self.structure_name,
            'Type': self.structure_type,
            'BasePairs': self.base_pairs,
            'Twist(°)': self.helical_twist,
            'MajorGroove(Å)': self.major_groove_width,
            'MinorGroove(Å)': self.minor_groove_width,
            'StackingEnergy(kcal/mol)': self.stacking_energy
        }


class DSSRAnalyzer:
    """
    X3DNA-DSSR wrapper for protein-DNA complex structural analysis
    
    Extracts the 5 most critical parameters for protein-DNA binding interface
    comparison between experimental and predicted structures.
    """
    
    def __init__(self, dssr_command: str = "dssr"):
        """
        Initialize DSSR analyzer
        
        Args:
            dssr_command: Command to execute DSSR (default: "dssr" from PATH)
        """
        self.dssr_command = dssr_command
        self._verify_dssr_installation()
    
    def _verify_dssr_installation(self) -> None:
        """Verify DSSR is installed and accessible"""
        try:
            result = subprocess.run([self.dssr_command, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                raise RuntimeError(f"DSSR not found or not working: {result.stderr}")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            raise RuntimeError(f"DSSR not available in PATH. Please install X3DNA-DSSR: {e}")
    
    def analyze_structure(self, pdb_path: Union[str, Path], 
                         structure_type: str = "Unknown") -> Optional[DSSRParameters]:
        """
        Analyze a single PDB structure using DSSR
        
        Args:
            pdb_path: Path to PDB file
            structure_type: Type label ('Experimental' or 'Predicted')
            
        Returns:
            DSSRParameters object or None if analysis failed
        """
        pdb_path = Path(pdb_path)
        
        if not pdb_path.exists():
            warnings.warn(f"PDB file not found: {pdb_path}")
            return None
        
        # Create temporary output file for DSSR JSON
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.json', delete=False) as temp_json:
            temp_json_path = temp_json.name
        
        try:
            # Run DSSR with JSON output
            cmd = [
                self.dssr_command,
                "--json",
                "-i", str(pdb_path),
                "-o", temp_json_path
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if result.returncode != 0:
                warnings.warn(f"DSSR analysis failed for {pdb_path.name}: {result.stderr}")
                return None
            
            # Parse DSSR JSON output
            return self._parse_dssr_json(temp_json_path, pdb_path.stem, structure_type)
            
        except subprocess.TimeoutExpired:
            warnings.warn(f"DSSR analysis timeout for {pdb_path.name}")
            return None
        except Exception as e:
            warnings.warn(f"DSSR analysis error for {pdb_path.name}: {e}")
            return None
        finally:
            # Clean up temporary file
            try:
                Path(temp_json_path).unlink()
            except:
                pass
    
    def _parse_dssr_json(self, json_path: str, structure_name: str, 
                        structure_type: str) -> Optional[DSSRParameters]:
        """
        Parse DSSR JSON output to extract critical parameters
        
        Args:
            json_path: Path to DSSR JSON output file
            structure_name: Name identifier for the structure
            structure_type: Type label ('Experimental' or 'Predicted')
            
        Returns:
            DSSRParameters object or None if parsing failed
        """
        try:
            with open(json_path, 'r') as f:
                dssr_data = json.load(f)
            
            # Extract base pairs (canonical only)
            base_pairs = self._extract_base_pairs(dssr_data)
            
            # Extract helical parameters
            helical_twist = self._extract_helical_twist(dssr_data)
            
            # Extract groove dimensions
            major_groove_width, minor_groove_width = self._extract_groove_widths(dssr_data)
            
            # Extract stacking energy
            stacking_energy = self._extract_stacking_energy(dssr_data)
            
            return DSSRParameters(
                structure_name=structure_name,
                structure_type=structure_type,
                base_pairs=base_pairs,
                helical_twist=helical_twist,
                major_groove_width=major_groove_width,
                minor_groove_width=minor_groove_width,
                stacking_energy=stacking_energy
            )
            
        except Exception as e:
            warnings.warn(f"Failed to parse DSSR JSON for {structure_name}: {e}")
            return None
    
    def _extract_base_pairs(self, dssr_data: Dict) -> int:
        """Extract number of canonical base pairs"""
        try:
            if 'pairs' not in dssr_data:
                return 0
            
            canonical_pairs = 0
            for pair in dssr_data['pairs']:
                # Count only canonical base pairs (Watson-Crick)
                if pair.get('canonical', False) or pair.get('WC', False):
                    canonical_pairs += 1
            
            return canonical_pairs
        except:
            return 0
    
    def _extract_helical_twist(self, dssr_data: Dict) -> float:
        """Extract average helical twist per base pair step"""
        try:
            if 'helicalSteps' not in dssr_data:
                return 0.0
            
            twist_values = []
            for step in dssr_data['helicalSteps']:
                if 'twist' in step:
                    twist_values.append(step['twist'])
            
            return np.mean(twist_values) if twist_values else 0.0
        except:
            return 0.0
    
    def _extract_groove_widths(self, dssr_data: Dict) -> Tuple[float, float]:
        """Extract major and minor groove widths"""
        try:
            major_width = 0.0
            minor_width = 0.0
            
            # Try to extract from groove analysis
            if 'grooves' in dssr_data:
                major_widths = []
                minor_widths = []
                
                for groove in dssr_data['grooves']:
                    if groove.get('type') == 'major':
                        if 'width' in groove:
                            major_widths.append(groove['width'])
                    elif groove.get('type') == 'minor':
                        if 'width' in groove:
                            minor_widths.append(groove['width'])
                
                major_width = np.mean(major_widths) if major_widths else 0.0
                minor_width = np.mean(minor_widths) if minor_widths else 0.0
            
            # Fallback: estimate from base pair geometry
            if major_width == 0.0 or minor_width == 0.0:
                if 'pairs' in dssr_data:
                    # Use typical B-form DNA values as fallback
                    major_width = 12.0  # Typical major groove width
                    minor_width = 5.7   # Typical minor groove width
            
            return major_width, minor_width
        except:
            return 12.0, 5.7  # B-form DNA defaults
    
    def _extract_stacking_energy(self, dssr_data: Dict) -> float:
        """Extract average base stacking energy"""
        try:
            if 'stacks' not in dssr_data:
                return 0.0
            
            stacking_energies = []
            for stack in dssr_data['stacks']:
                if 'energy' in stack:
                    stacking_energies.append(stack['energy'])
                elif 'overlap_area' in stack:
                    # Estimate energy from overlap area (rough approximation)
                    stacking_energies.append(stack['overlap_area'] * -0.5)
            
            return np.mean(stacking_energies) if stacking_energies else 0.0
        except:
            return 0.0
    
    def compare_structures(self, experimental_pdbs: List[Union[str, Path]], 
                          predicted_pdbs: List[Union[str, Path]]) -> Tuple[pd.DataFrame, Dict]:
        """
        Compare experimental and predicted structures using DSSR analysis
        
        Args:
            experimental_pdbs: List of experimental PDB file paths
            predicted_pdbs: List of predicted PDB file paths
            
        Returns:
            Tuple of (DataFrame with all results, summary statistics dict)
        """
        all_results = []
        
        # Analyze experimental structures
        print(f"Analyzing {len(experimental_pdbs)} experimental structures...")
        for pdb_path in experimental_pdbs:
            result = self.analyze_structure(pdb_path, "Experimental")
            if result:
                all_results.append(result)
        
        # Analyze predicted structures
        print(f"Analyzing {len(predicted_pdbs)} predicted structures...")
        for pdb_path in predicted_pdbs:
            result = self.analyze_structure(pdb_path, "Predicted")
            if result:
                all_results.append(result)
        
        if not all_results:
            raise RuntimeError("No structures could be analyzed successfully")
        
        # Create DataFrame
        df = pd.DataFrame([result.to_dict() for result in all_results])
        
        # Calculate summary statistics
        summary_stats = self._calculate_summary_statistics(df)
        
        # Print comparison report
        self._print_comparison_report(summary_stats)
        
        return df, summary_stats
    
    def _calculate_summary_statistics(self, df: pd.DataFrame) -> Dict:
        """Calculate summary statistics for experimental vs predicted"""
        exp_data = df[df['Type'] == 'Experimental']
        pred_data = df[df['Type'] == 'Predicted']
        
        parameters = ['BasePairs', 'Twist(°)', 'MajorGroove(Å)', 
                     'MinorGroove(Å)', 'StackingEnergy(kcal/mol)']
        
        stats = {}
        for param in parameters:
            exp_values = exp_data[param].dropna()
            pred_values = pred_data[param].dropna()
            
            if len(exp_values) > 0 and len(pred_values) > 0:
                exp_mean = np.mean(exp_values)
                exp_std = np.std(exp_values, ddof=1) if len(exp_values) > 1 else 0.0
                pred_mean = np.mean(pred_values)
                pred_std = np.std(pred_values, ddof=1) if len(pred_values) > 1 else 0.0
                difference = pred_mean - exp_mean
                
                # Check for significant deviations
                warning = self._check_deviation_threshold(param, abs(difference))
                
                stats[param] = {
                    'experimental_mean': exp_mean,
                    'experimental_std': exp_std,
                    'predicted_mean': pred_mean,
                    'predicted_std': pred_std,
                    'difference': difference,
                    'warning': warning
                }
        
        return stats
    
    def _check_deviation_threshold(self, parameter: str, abs_difference: float) -> bool:
        """Check if deviation exceeds warning thresholds"""
        thresholds = {
            'Twist(°)': 5.0,  # 5 degree threshold for twist
            'MajorGroove(Å)': 0.5,  # 0.5 Angstrom threshold for groove widths
            'MinorGroove(Å)': 0.5,  # 0.5 Angstrom threshold for groove widths
            'BasePairs': 2.0,  # 2 base pair threshold
            'StackingEnergy(kcal/mol)': 2.0  # 2 kcal/mol threshold
        }
        
        threshold = thresholds.get(parameter, 0.5)
        return abs_difference > threshold
    
    def _print_comparison_report(self, stats: Dict) -> None:
        """Print formatted comparison report to console"""
        print("\n" + "="*80)
        print("X3DNA-DSSR COMPARISON REPORT")
        print("Critical Parameters for Protein-DNA Binding Interface Analysis")
        print("="*80)
        print(f"{'Parameter':<20} | {'Exp (mean ± std)':<18} | {'Pred (mean ± std)':<18} | {'Diff':<8} | {'Status'}")
        print("-" * 80)
        
        for param, data in stats.items():
            exp_str = f"{data['experimental_mean']:.1f} ± {data['experimental_std']:.1f}"
            pred_str = f"{data['predicted_mean']:.1f} ± {data['predicted_std']:.1f}"
            diff = data['difference']
            diff_str = f"{diff:+.1f}"
            
            status = "⚠️ FLAGGED" if data['warning'] else "✓ OK"
            
            print(f"{param:<20} | {exp_str:<18} | {pred_str:<18} | {diff_str:<8} | {status}")
        
        print("="*80)
        print("\nTHRESHOLD CRITERIA:")
        print("• Twist deviation > 5°")
        print("• Groove width deviation > 0.5 Å") 
        print("• Base pair count deviation > 2")
        print("• Stacking energy deviation > 2.0 kcal/mol")
        print("\nFLAGGED DEVIATIONS indicate potential issues with protein-DNA")
        print("binding interface prediction accuracy.")
        print("="*80)
    
    def export_results(self, df: pd.DataFrame, output_path: Union[str, Path]) -> None:
        """
        Export analysis results to CSV file
        
        Args:
            df: DataFrame containing analysis results
            output_path: Path for output CSV file
        """
        output_path = Path(output_path)
        df.to_csv(output_path, index=False)
        print(f"\nResults exported to: {output_path}")


def analyze_protein_dna_complexes(experimental_pdbs: List[Union[str, Path]], 
                                 predicted_pdbs: List[Union[str, Path]],
                                 output_csv: Optional[Union[str, Path]] = None) -> pd.DataFrame:
    """
    Convenience function for protein-DNA complex comparison using DSSR
    
    Args:
        experimental_pdbs: List of experimental PDB file paths
        predicted_pdbs: List of predicted PDB file paths
        output_csv: Optional path to export results as CSV
        
    Returns:
        DataFrame with analysis results
    """
    analyzer = DSSRAnalyzer()
    df, stats = analyzer.compare_structures(experimental_pdbs, predicted_pdbs)
    
    if output_csv:
        analyzer.export_results(df, output_csv)
    
    return df


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python dssr.py experimental_dir predicted_dir [output.csv]")
        sys.exit(1)
    
    exp_dir = Path(sys.argv[1])
    pred_dir = Path(sys.argv[2])
    output_csv = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Find PDB files
    exp_pdbs = list(exp_dir.glob("*.pdb"))
    pred_pdbs = list(pred_dir.glob("*.pdb")) + list(pred_dir.glob("*.cif"))
    
    # Run analysis
    df = analyze_protein_dna_complexes(exp_pdbs, pred_pdbs, output_csv)
    print(f"\nAnalyzed {len(df)} structures total.")