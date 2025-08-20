"""
biostructbenchmark/analysis/dssr.py

DSSR (3DNA suite) integration for nucleic acid geometry analysis
Replaces CURVES+ with the more accessible and widely-used DSSR tool
"""

import os
import subprocess
import tempfile
import shutil
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
import pandas as pd
import numpy as np
from Bio.PDB import Structure

from biostructbenchmark.core.io import get_structure


@dataclass
class BasePairParameters:
    """Container for DSSR base pair parameters"""
    # Base pair identifiers
    bp_name: str  # e.g., "A.1_T.22"
    strand1: str  # First nucleotide
    strand2: str  # Second nucleotide
    
    # Geometric parameters (DSSR provides these in JSON output)
    shear: Optional[float] = None
    stretch: Optional[float] = None
    stagger: Optional[float] = None
    buckle: Optional[float] = None
    propeller: Optional[float] = None
    opening: Optional[float] = None
    
    # Additional DSSR-specific parameters
    overlap_area: Optional[float] = None
    h_bonds: Optional[int] = None
    bp_type: Optional[str] = None  # Watson-Crick, Wobble, etc.


@dataclass
class BaseStepParameters:
    """Container for DSSR base step parameters"""
    # Step identifiers
    step_name: str  # e.g., "AA/TT"
    bp1: str  # First base pair
    bp2: str  # Second base pair
    
    # Step parameters
    shift: Optional[float] = None
    slide: Optional[float] = None
    rise: Optional[float] = None
    tilt: Optional[float] = None
    roll: Optional[float] = None
    twist: Optional[float] = None
    
    # Additional step geometry
    helical_rise: Optional[float] = None
    helical_twist: Optional[float] = None
    x_displacement: Optional[float] = None
    y_displacement: Optional[float] = None
    z_displacement: Optional[float] = None


@dataclass
class GrooveParameters:
    """Container for DSSR groove geometry parameters"""
    # Position identifier
    position: str  # Position in sequence
    
    # Groove widths and depths
    major_groove_width: Optional[float] = None
    minor_groove_width: Optional[float] = None
    major_groove_depth: Optional[float] = None
    minor_groove_depth: Optional[float] = None
    
    # Additional groove properties
    major_groove_distance: Optional[float] = None
    minor_groove_distance: Optional[float] = None


@dataclass
class NucleotideGeometry:
    """Container for individual nucleotide geometry"""
    # Nucleotide identifier
    nucleotide_id: str  # e.g., "A.1"
    chain: str
    position: int
    base_type: str  # A, T, G, C, U
    
    # Sugar pucker parameters
    pucker_amplitude: Optional[float] = None
    pucker_phase: Optional[float] = None
    pucker_type: Optional[str] = None  # C2'-endo, C3'-endo, etc.
    
    # Glycosidic bond
    chi_angle: Optional[float] = None
    
    # Backbone torsions
    alpha: Optional[float] = None
    beta: Optional[float] = None
    gamma: Optional[float] = None
    delta: Optional[float] = None
    epsilon: Optional[float] = None
    zeta: Optional[float] = None


@dataclass
class ProteinNAContact:
    """Container for protein-nucleic acid contacts"""
    # Contact identifiers
    protein_residue: str  # e.g., "ARG.45"
    nucleotide: str  # e.g., "A.10"
    
    # Contact details
    contact_type: str  # hydrogen_bond, van_der_waals, stacking
    distance: float
    
    # Atoms involved
    protein_atom: str
    nucleotide_atom: str
    
    # Optional parameters
    angle: Optional[float] = None
    strength: Optional[str] = None  # strong, medium, weak


class DSSRAnalyzer:
    """Main class for DSSR nucleic acid structure analysis"""
    
    def __init__(self, dssr_executable: Optional[str] = None):
        """
        Initialize DSSR analyzer
        
        Args:
            dssr_executable: Path to DSSR executable. If None, searches PATH
        """
        self.dssr_exe = self._find_dssr_executable(dssr_executable)
        if not self.dssr_exe:
            raise RuntimeError(
                "DSSR executable not found. Please install DSSR from https://x3dna.org/ "
                "and ensure it's in PATH or provide explicit path"
            )
    
    def _find_dssr_executable(self, dssr_path: Optional[str]) -> Optional[str]:
        """Find DSSR executable in common locations"""
        if dssr_path and os.path.isfile(dssr_path):
            return dssr_path
        
        # Check common locations
        common_paths = [
            'dssr',  # Standard name
            'x3dna-dssr',
            '/usr/local/bin/dssr',
            '/opt/3dna/bin/dssr',
            os.path.expanduser('~/3dna/bin/dssr'),
            os.path.expanduser('~/dssr/dssr')
        ]
        
        for path in common_paths:
            if shutil.which(path):
                return path
        
        return None
    
    def analyze_structure(self, structure_path: Path, 
                         output_dir: Optional[Path] = None,
                         options: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Run DSSR analysis on a nucleic acid structure
        
        Args:
            structure_path: Path to PDB/CIF file
            output_dir: Directory for output files (temporary if None)
            options: Additional DSSR options
            
        Returns:
            Dictionary containing parsed DSSR analysis results
        """
        # Create temporary directory if needed
        temp_dir = output_dir or Path(tempfile.mkdtemp(prefix='dssr_'))
        temp_dir.mkdir(exist_ok=True)
        
        try:
            # Run DSSR analysis
            dssr_output = self._run_dssr(structure_path, temp_dir, options or {})
            
            # Parse DSSR JSON output
            results = self._parse_dssr_output(dssr_output)
            
            return results
            
        finally:
            # Clean up temporary files if we created them
            if output_dir is None:
                shutil.rmtree(temp_dir)
    
    def _run_dssr(self, structure_path: Path, output_dir: Path, 
                  options: Dict[str, Any]) -> Path:
        """Execute DSSR on nucleic acid structure"""
        
        # DSSR output files
        json_output = output_dir / f"{structure_path.stem}_dssr.json"
        
        # Build DSSR command
        cmd = [
            self.dssr_exe,
            '--input', str(structure_path),
            '--output', str(json_output),
            '--json',  # JSON output format
            '--more',  # Include additional analysis
        ]
        
        # Add optional parameters
        if options.get('no_image', True):
            cmd.append('--no-image')
        
        if options.get('auxiliary', False):
            cmd.append('--auxiliary')
        
        if options.get('hydrogen_bonds', True):
            cmd.append('--hbonds')
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=output_dir,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"DSSR failed: {result.stderr}")
            
            if not json_output.exists():
                raise RuntimeError("DSSR did not generate expected JSON output")
            
            return json_output
            
        except subprocess.TimeoutExpired:
            raise RuntimeError("DSSR execution timed out")
        except FileNotFoundError:
            raise RuntimeError(f"DSSR executable not found: {self.dssr_exe}")
    
    def _parse_dssr_output(self, json_file: Path) -> Dict[str, Any]:
        """Parse DSSR JSON output into structured data"""
        
        with open(json_file, 'r') as f:
            dssr_data = json.load(f)
        
        results = {
            'base_pairs': [],
            'base_steps': [],
            'nucleotides': [],
            'protein_na_contacts': [],
            'groove_parameters': [],
            'summary': {}
        }
        
        # Parse base pairs
        if 'pairs' in dssr_data:
            for bp_data in dssr_data['pairs']:
                bp = self._parse_base_pair(bp_data)
                results['base_pairs'].append(bp)
        
        # Parse base steps
        if 'steps' in dssr_data:
            for step_data in dssr_data['steps']:
                step = self._parse_base_step(step_data)
                results['base_steps'].append(step)
        
        # Parse individual nucleotides
        if 'nucleotides' in dssr_data:
            for nt_data in dssr_data['nucleotides']:
                nt = self._parse_nucleotide(nt_data)
                results['nucleotides'].append(nt)
        
        # Parse protein-NA contacts (if present)
        if 'protein_na_contacts' in dssr_data:
            for contact_data in dssr_data['protein_na_contacts']:
                contact = self._parse_protein_na_contact(contact_data)
                results['protein_na_contacts'].append(contact)
        
        # Parse groove parameters
        if 'groove_widths' in dssr_data:
            for groove_data in dssr_data['groove_widths']:
                groove = self._parse_groove_parameters(groove_data)
                results['groove_parameters'].append(groove)
        
        # Extract summary information
        results['summary'] = self._extract_summary(dssr_data)
        
        return results
    
    def _parse_base_pair(self, bp_data: Dict) -> BasePairParameters:
        """Parse base pair data from DSSR JSON"""
        return BasePairParameters(
            bp_name=bp_data.get('name', ''),
            strand1=bp_data.get('nt1', ''),
            strand2=bp_data.get('nt2', ''),
            shear=bp_data.get('shear'),
            stretch=bp_data.get('stretch'),
            stagger=bp_data.get('stagger'),
            buckle=bp_data.get('buckle'),
            propeller=bp_data.get('propeller'),
            opening=bp_data.get('opening'),
            overlap_area=bp_data.get('overlap_area'),
            h_bonds=bp_data.get('h_bonds'),
            bp_type=bp_data.get('type')
        )
    
    def _parse_base_step(self, step_data: Dict) -> BaseStepParameters:
        """Parse base step data from DSSR JSON"""
        return BaseStepParameters(
            step_name=step_data.get('name', ''),
            bp1=step_data.get('bp1', ''),
            bp2=step_data.get('bp2', ''),
            shift=step_data.get('shift'),
            slide=step_data.get('slide'),
            rise=step_data.get('rise'),
            tilt=step_data.get('tilt'),
            roll=step_data.get('roll'),
            twist=step_data.get('twist'),
            helical_rise=step_data.get('helical_rise'),
            helical_twist=step_data.get('helical_twist'),
            x_displacement=step_data.get('x_displacement'),
            y_displacement=step_data.get('y_displacement'),
            z_displacement=step_data.get('z_displacement')
        )
    
    def _parse_nucleotide(self, nt_data: Dict) -> NucleotideGeometry:
        """Parse individual nucleotide geometry from DSSR JSON"""
        return NucleotideGeometry(
            nucleotide_id=nt_data.get('id', ''),
            chain=nt_data.get('chain', ''),
            position=nt_data.get('position', 0),
            base_type=nt_data.get('base', ''),
            pucker_amplitude=nt_data.get('pucker_amplitude'),
            pucker_phase=nt_data.get('pucker_phase'),
            pucker_type=nt_data.get('pucker_type'),
            chi_angle=nt_data.get('chi'),
            alpha=nt_data.get('alpha'),
            beta=nt_data.get('beta'),
            gamma=nt_data.get('gamma'),
            delta=nt_data.get('delta'),
            epsilon=nt_data.get('epsilon'),
            zeta=nt_data.get('zeta')
        )
    
    def _parse_protein_na_contact(self, contact_data: Dict) -> ProteinNAContact:
        """Parse protein-nucleic acid contact from DSSR JSON"""
        return ProteinNAContact(
            protein_residue=contact_data.get('protein_residue', ''),
            nucleotide=contact_data.get('nucleotide', ''),
            contact_type=contact_data.get('contact_type', ''),
            distance=contact_data.get('distance', 0.0),
            angle=contact_data.get('angle'),
            protein_atom=contact_data.get('protein_atom', ''),
            nucleotide_atom=contact_data.get('nucleotide_atom', ''),
            strength=contact_data.get('strength')
        )
    
    def _parse_groove_parameters(self, groove_data: Dict) -> GrooveParameters:
        """Parse groove parameters from DSSR JSON"""
        return GrooveParameters(
            position=groove_data.get('position', ''),
            major_groove_width=groove_data.get('major_groove_width'),
            minor_groove_width=groove_data.get('minor_groove_width'),
            major_groove_depth=groove_data.get('major_groove_depth'),
            minor_groove_depth=groove_data.get('minor_groove_depth'),
            major_groove_distance=groove_data.get('major_groove_distance'),
            minor_groove_distance=groove_data.get('minor_groove_distance')
        )
    
    def _extract_summary(self, dssr_data: Dict) -> Dict:
        """Extract summary statistics from DSSR analysis"""
        summary = {}
        
        # Basic counts
        summary['num_base_pairs'] = len(dssr_data.get('pairs', []))
        summary['num_base_steps'] = len(dssr_data.get('steps', []))
        summary['num_nucleotides'] = len(dssr_data.get('nucleotides', []))
        summary['num_chains'] = len(dssr_data.get('chains', []))
        
        # Structure type
        summary['structure_type'] = dssr_data.get('structure_type', 'unknown')
        
        # Sequence information
        if 'chains' in dssr_data:
            sequences = []
            for chain in dssr_data['chains']:
                if 'sequence' in chain:
                    sequences.append(chain['sequence'])
            summary['sequences'] = sequences
        
        return summary
    
    def compare_structures(self, experimental_results: Dict, 
                          predicted_results: Dict) -> pd.DataFrame:
        """
        Compare DSSR analysis results between experimental and predicted structures
        
        Args:
            experimental_results: DSSR results from experimental structure
            predicted_results: DSSR results from predicted structure
            
        Returns:
            DataFrame with parameter differences
        """
        comparisons = []
        
        # Compare base pair parameters
        exp_bps = {bp.bp_name: bp for bp in experimental_results['base_pairs']}
        pred_bps = {bp.bp_name: bp for bp in predicted_results['base_pairs']}
        
        common_bps = set(exp_bps.keys()) & set(pred_bps.keys())
        
        for bp_name in common_bps:
            exp_bp = exp_bps[bp_name]
            pred_bp = pred_bps[bp_name]
            
            comparison = {
                'type': 'base_pair',
                'name': bp_name,
                'shear_diff': self._safe_subtract(pred_bp.shear, exp_bp.shear),
                'stretch_diff': self._safe_subtract(pred_bp.stretch, exp_bp.stretch),
                'stagger_diff': self._safe_subtract(pred_bp.stagger, exp_bp.stagger),
                'buckle_diff': self._safe_subtract(pred_bp.buckle, exp_bp.buckle),
                'propeller_diff': self._safe_subtract(pred_bp.propeller, exp_bp.propeller),
                'opening_diff': self._safe_subtract(pred_bp.opening, exp_bp.opening),
            }
            comparisons.append(comparison)
        
        # Compare base step parameters
        exp_steps = {step.step_name: step for step in experimental_results['base_steps']}
        pred_steps = {step.step_name: step for step in predicted_results['base_steps']}
        
        common_steps = set(exp_steps.keys()) & set(pred_steps.keys())
        
        for step_name in common_steps:
            exp_step = exp_steps[step_name]
            pred_step = pred_steps[step_name]
            
            comparison = {
                'type': 'base_step',
                'name': step_name,
                'shift_diff': self._safe_subtract(pred_step.shift, exp_step.shift),
                'slide_diff': self._safe_subtract(pred_step.slide, exp_step.slide),
                'rise_diff': self._safe_subtract(pred_step.rise, exp_step.rise),
                'tilt_diff': self._safe_subtract(pred_step.tilt, exp_step.tilt),
                'roll_diff': self._safe_subtract(pred_step.roll, exp_step.roll),
                'twist_diff': self._safe_subtract(pred_step.twist, exp_step.twist),
            }
            comparisons.append(comparison)
        
        return pd.DataFrame(comparisons)
    
    def _safe_subtract(self, a: Optional[float], b: Optional[float]) -> Optional[float]:
        """Safely subtract two values that might be None"""
        if a is not None and b is not None:
            return a - b
        return None
    
    def analyze_protein_dna_interface(self, structure_path: Path) -> Dict[str, Any]:
        """
        Specialized analysis for protein-DNA interfaces
        
        Args:
            structure_path: Path to protein-DNA complex structure
            
        Returns:
            Dictionary with interface analysis results
        """
        # Run DSSR with additional options for protein-DNA analysis
        options = {
            'auxiliary': True,
            'hydrogen_bonds': True,
            'no_image': True
        }
        
        results = self.analyze_structure(structure_path, options=options)
        
        # Add interface-specific analysis
        interface_analysis = {
            'dssr_results': results,
            'interface_statistics': self._calculate_interface_stats(results),
            'binding_sites': self._identify_binding_sites(results),
            'geometric_distortions': self._analyze_geometric_distortions(results)
        }
        
        return interface_analysis
    
    def _calculate_interface_stats(self, results: Dict) -> Dict:
        """Calculate statistics for protein-DNA interface"""
        stats = {}
        
        contacts = results.get('protein_na_contacts', [])
        stats['total_contacts'] = len(contacts)
        
        # Count by contact type
        contact_types = {}
        for contact in contacts:
            contact_type = contact.contact_type
            contact_types[contact_type] = contact_types.get(contact_type, 0) + 1
        
        stats['contact_types'] = contact_types
        
        # Distance statistics
        distances = [c.distance for c in contacts if c.distance is not None]
        if distances:
            stats['mean_contact_distance'] = np.mean(distances)
            stats['std_contact_distance'] = np.std(distances)
            stats['min_contact_distance'] = np.min(distances)
            stats['max_contact_distance'] = np.max(distances)
        
        return stats
    
    def _identify_binding_sites(self, results: Dict) -> List[Dict]:
        """Identify major protein-DNA binding sites"""
        # This would implement sophisticated binding site identification
        # For now, return a placeholder structure
        binding_sites = []
        
        contacts = results.get('protein_na_contacts', [])
        
        # Group contacts by protein residue
        residue_contacts = {}
        for contact in contacts:
            res = contact.protein_residue
            if res not in residue_contacts:
                residue_contacts[res] = []
            residue_contacts[res].append(contact)
        
        # Identify residues with multiple contacts as potential binding sites
        for residue, res_contacts in residue_contacts.items():
            if len(res_contacts) >= 2:  # Threshold for binding site
                binding_site = {
                    'protein_residue': residue,
                    'num_contacts': len(res_contacts),
                    'contact_types': [c.contact_type for c in res_contacts],
                    'nucleotides_contacted': [c.nucleotide for c in res_contacts]
                }
                binding_sites.append(binding_site)
        
        return binding_sites
    
    def _analyze_geometric_distortions(self, results: Dict) -> Dict:
        """Analyze geometric distortions in protein-bound DNA"""
        distortions = {}
        
        # Analyze base pair parameters for unusual values
        base_pairs = results.get('base_pairs', [])
        
        # Define normal ranges (these would be refined based on literature)
        normal_ranges = {
            'shear': (-2.0, 2.0),
            'stretch': (-1.0, 1.0),
            'stagger': (-1.0, 1.0),
            'buckle': (-30.0, 30.0),
            'propeller': (-30.0, 30.0),
            'opening': (-10.0, 10.0)
        }
        
        distorted_bps = []
        for bp in base_pairs:
            bp_distortions = {}
            for param, (min_val, max_val) in normal_ranges.items():
                value = getattr(bp, param, None)
                if value is not None and (value < min_val or value > max_val):
                    bp_distortions[param] = value
            
            if bp_distortions:
                distorted_bps.append({
                    'base_pair': bp.bp_name,
                    'distortions': bp_distortions
                })
        
        distortions['distorted_base_pairs'] = distorted_bps
        distortions['num_distorted'] = len(distorted_bps)
        
        return distortions
    
    def export_results(self, results: Dict, output_path: Path) -> None:
        """
        Export DSSR analysis results to files
        
        Args:
            results: DSSR analysis results dictionary
            output_path: Base path for output files
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Export base pairs
        if results['base_pairs']:
            bp_data = []
            for bp in results['base_pairs']:
                bp_dict = {
                    'bp_name': bp.bp_name,
                    'strand1': bp.strand1,
                    'strand2': bp.strand2,
                    'shear': bp.shear,
                    'stretch': bp.stretch,
                    'stagger': bp.stagger,
                    'buckle': bp.buckle,
                    'propeller': bp.propeller,
                    'opening': bp.opening,
                    'overlap_area': bp.overlap_area,
                    'h_bonds': bp.h_bonds,
                    'bp_type': bp.bp_type
                }
                bp_data.append(bp_dict)
            
            bp_df = pd.DataFrame(bp_data)
            bp_df.to_csv(f"{output_path}_base_pairs.csv", index=False)
        
        # Export base steps
        if results['base_steps']:
            step_data = []
            for step in results['base_steps']:
                step_dict = {
                    'step_name': step.step_name,
                    'bp1': step.bp1,
                    'bp2': step.bp2,
                    'shift': step.shift,
                    'slide': step.slide,
                    'rise': step.rise,
                    'tilt': step.tilt,
                    'roll': step.roll,
                    'twist': step.twist,
                    'helical_rise': step.helical_rise,
                    'helical_twist': step.helical_twist
                }
                step_data.append(step_dict)
            
            step_df = pd.DataFrame(step_data)
            step_df.to_csv(f"{output_path}_base_steps.csv", index=False)
        
        # Export protein-NA contacts
        if results['protein_na_contacts']:
            contact_data = []
            for contact in results['protein_na_contacts']:
                contact_dict = {
                    'protein_residue': contact.protein_residue,
                    'nucleotide': contact.nucleotide,
                    'contact_type': contact.contact_type,
                    'distance': contact.distance,
                    'angle': contact.angle,
                    'protein_atom': contact.protein_atom,
                    'nucleotide_atom': contact.nucleotide_atom,
                    'strength': contact.strength
                }
                contact_data.append(contact_dict)
            
            contact_df = pd.DataFrame(contact_data)
            contact_df.to_csv(f"{output_path}_protein_na_contacts.csv", index=False)
        
        # Export summary as JSON
        with open(f"{output_path}_summary.json", 'w') as f:
            json.dump(results['summary'], f, indent=2)


def analyze_dna_batch(experimental_dir: Path, predicted_dir: Path,
                     output_dir: Path) -> None:
    """
    Batch analysis of DNA geometry using DSSR for structure pairs
    
    Args:
        experimental_dir: Directory containing experimental structures
        predicted_dir: Directory containing predicted structures  
        output_dir: Directory for results
    """
    analyzer = DSSRAnalyzer()
    output_dir.mkdir(exist_ok=True)
    
    # Find matching structure pairs
    exp_files = list(experimental_dir.glob("*.pdb")) + list(experimental_dir.glob("*.cif"))
    
    all_comparisons = []
    
    for exp_file in exp_files:
        # Find corresponding predicted structure
        pred_file = None
        for ext in ['.pdb', '.cif']:
            potential_pred = predicted_dir / f"{exp_file.stem}{ext}"
            if potential_pred.exists():
                pred_file = potential_pred
                break
        
        if not pred_file:
            print(f"Warning: No predicted structure found for {exp_file.name}")
            continue
        
        print(f"Analyzing: {exp_file.name} vs {pred_file.name}")
        
        try:
            # Analyze both structures
            exp_results = analyzer.analyze_structure(exp_file)
            pred_results = analyzer.analyze_structure(pred_file)
            
            # Compare parameters
            comparison_df = analyzer.compare_structures(exp_results, pred_results)
            comparison_df['structure_pair'] = exp_file.stem
            
            all_comparisons.append(comparison_df)
            
            # Save individual results
            structure_output = output_dir / exp_file.stem
            structure_output.mkdir(exist_ok=True)
            
            comparison_df.to_csv(structure_output / 'dssr_comparison.csv', index=False)
            
            # Export full DSSR results
            analyzer.export_results(exp_results, structure_output / 'experimental_dssr')
            analyzer.export_results(pred_results, structure_output / 'predicted_dssr')
            
        except Exception as e:
            print(f"Error processing {exp_file.name}: {e}")
            continue
    
    # Combine all comparisons
    if all_comparisons:
        combined_df = pd.concat(all_comparisons, ignore_index=True)
        combined_df.to_csv(output_dir / 'all_dssr_comparisons.csv', index=False)
        
        # Generate summary statistics
        numeric_cols = combined_df.select_dtypes(include=[np.number]).columns
        summary_stats = combined_df[numeric_cols].describe()
        summary_stats.to_csv(output_dir / 'dssr_summary_statistics.csv')
        
        print(f"DSSR analysis complete. Results saved to {output_dir}")
    else:
        print("No successful analyses completed")


if __name__ == "__main__":
    # Example usage (will work once DSSR is installed)
    try:
        analyzer = DSSRAnalyzer()
        
        # Test on a single structure
        test_structure = Path("test_structure.pdb")
        if test_structure.exists():
            results = analyzer.analyze_structure(test_structure)
            
            print(f"Found {len(results['base_pairs'])} base pairs")
            print(f"Found {len(results['base_steps'])} base steps")
            print(f"Found {len(results['protein_na_contacts'])} protein-NA contacts")
        
    except RuntimeError as e:
        print(f"DSSR not available: {e}")
        print("Please install DSSR from https://x3dna.org/ to use this module")