#!/usr/bin/env python3
"""
Integration test for core components: io.py, alignment.py, metrics.py
Tests the complete pipeline functionality and saves outputs for inspection.
"""

import sys
import os
from pathlib import Path
import traceback
from datetime import datetime

# Add the package to Python path
sys.path.insert(0, str(Path(__file__).parent))

def test_core_integration():
    """Test all core components working together"""
    print("="*60)
    print("BioStructBenchmark Core Components Integration Test")
    print("="*60)
    print(f"Test started: {datetime.now()}")
    
    # Create output directory
    output_dir = Path("test_core_outputs")
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir.absolute()}")
    
    try:
        # Import core modules
        print("\n1. Testing imports...")
        from biostructbenchmark.core import io, alignment, metrics
        from biostructbenchmark.core.io import get_structure, validate_file
        from biostructbenchmark.core.alignment import compare_structures, align_structures_three_frames
        from biostructbenchmark.core.metrics import generate_comprehensive_metrics
        print("✓ All imports successful")
        
        # Test data files - use similar protein structures
        test_data_dir = Path("tests/data")
        observed_file = test_data_dir / "proteins_pdb" / "p456_02_experimental.pdb"
        predicted_file = test_data_dir / "proteins_cif" / "p456_02_predicted.cif"  # Use matching structure
        
        print(f"\n2. Testing io.py functionality...")
        print(f"   Observed file: {observed_file}")
        print(f"   Predicted file: {predicted_file}")
        
        # Test file validation
        if not observed_file.exists():
            print(f"   ERROR: {observed_file} not found")
            return False
        if not predicted_file.exists():
            print(f"   ERROR: {predicted_file} not found")
            return False
            
        print("   Testing file validation...")
        obs_valid = validate_file(observed_file)
        pred_valid = validate_file(predicted_file)
        print(f"   ✓ Observed file valid: {obs_valid}")
        print(f"   ✓ Predicted file valid: {pred_valid}")
        
        # Test structure loading
        print("   Testing structure loading...")
        observed_structure = get_structure(observed_file)
        predicted_structure = get_structure(predicted_file)
        
        if observed_structure is None:
            print(f"   ERROR: Failed to load {observed_file}")
            return False
        if predicted_structure is None:
            print(f"   ERROR: Failed to load {predicted_file}")
            return False
            
        print(f"   ✓ Loaded observed structure: {len(list(observed_structure[0]))} chains")
        print(f"   ✓ Loaded predicted structure: {len(list(predicted_structure[0]))} chains")
        
        print(f"\n3. Testing alignment.py functionality...")
        print("   Running structure comparison...")
        
        # Test main comparison function - use structural superimposition method directly
        from biostructbenchmark.core.alignment import align_structures_with_structural_superimposition
        try:
            alignment_result = align_structures_with_structural_superimposition(observed_structure, predicted_structure, 'global')
        except ValueError as e:
            print(f"   Note: Global alignment failed ({e}), trying protein_centric...")
            try:
                alignment_result = align_structures_with_structural_superimposition(observed_structure, predicted_structure, 'protein_centric')
            except ValueError as e2:
                print(f"   Note: Protein_centric alignment failed ({e2}), using compare_structures...")
                alignment_result = compare_structures(observed_structure, predicted_structure)
        
        if alignment_result is None:
            print("   ERROR: Structure comparison failed")
            return False
            
        print(f"   ✓ Overall RMSD: {alignment_result.overall_rmsd:.3f} Å")
        print(f"   ✓ Aligned atoms: {alignment_result.aligned_atom_count}")
        print(f"   ✓ Residue RMSDs calculated: {len(alignment_result.residue_rmsds)}")
        print(f"   ✓ Chain RMSDs calculated: {len(alignment_result.chain_rmsds)}")
        
        # Test three-frame alignment
        print("   Testing three-frame alignment...")
        frame_results = align_structures_three_frames(observed_structure, predicted_structure)
        print(f"   ✓ Generated {len(frame_results)} reference frames")
        for frame, result in frame_results.items():
            print(f"     - {frame}: {result.overall_rmsd:.3f} Å")
        
        print(f"\n4. Testing metrics.py functionality...")
        print("   Generating comprehensive metrics...")
        
        # Test metrics generation
        structure_metrics = generate_comprehensive_metrics(alignment_result)
        
        print(f"   ✓ Overall RMSD: {structure_metrics.overall_rmsd:.3f} Å")
        print(f"   ✓ Protein RMSD: {structure_metrics.protein_rmsd:.3f} Å") 
        print(f"   ✓ DNA RMSD: {structure_metrics.dna_rmsd:.3f} Å")
        print(f"   ✓ Translation error: {structure_metrics.error_components.translation_error:.3f} Å")
        print(f"   ✓ Rotation angle: {structure_metrics.error_components.rotation_angle:.1f}°")
        print(f"   ✓ Best residues identified: {len(structure_metrics.best_residues)}")
        print(f"   ✓ Worst residues identified: {len(structure_metrics.worst_residues)}")
        
        print(f"\n5. Testing export functionality...")
        
        # Export alignment results
        from biostructbenchmark.core.alignment import export_comprehensive_alignment_report, export_residue_rmsd_csv
        
        print("   Exporting alignment results...")
        export_comprehensive_alignment_report(alignment_result, output_dir, "integration_test")
        
        # Export residue RMSD CSV
        residue_csv_path = output_dir / "integration_test_residue_rmsds.csv"
        export_residue_rmsd_csv(alignment_result.residue_rmsds, str(residue_csv_path))
        
        # Export metrics
        from biostructbenchmark.core.metrics import export_metrics_csv, print_metrics_summary
        
        print("   Exporting metrics...")
        metrics_csv_path = output_dir / "integration_test_metrics.csv"
        export_metrics_csv(structure_metrics, str(metrics_csv_path))
        
        # Print summary to console
        print("\n   Metrics Summary:")
        print_metrics_summary(structure_metrics)
        
        print(f"\n6. Verifying output files...")
        output_files = list(output_dir.glob("*"))
        print(f"   Generated {len(output_files)} output files:")
        for file in sorted(output_files):
            size = file.stat().st_size if file.is_file() else 0
            print(f"     - {file.name} ({size} bytes)")
        
        print(f"\n" + "="*60)
        print("INTEGRATION TEST SUCCESSFUL!")
        print("="*60)
        print(f"All outputs saved to: {output_dir.absolute()}")
        print(f"Test completed: {datetime.now()}")
        
        return True
        
    except Exception as e:
        print(f"\nERROR during integration test:")
        print(f"Exception: {e}")
        print("Traceback:")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_core_integration()
    sys.exit(0 if success else 1)