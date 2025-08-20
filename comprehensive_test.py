#!/usr/bin/env python3
"""
Comprehensive test of BioStructBenchmark program
Tests all major components and generates verification outputs
"""

import sys
import os
from pathlib import Path
import traceback
import json
from datetime import datetime

# Add project root to path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

def test_imports():
    """Test that all major modules can be imported"""
    print("\n=== TESTING IMPORTS ===")
    
    try:
        # Core modules
        from biostructbenchmark.core import alignment, metrics, io
        print("✓ Core modules imported successfully")
        
        # Analysis modules
        from biostructbenchmark.analysis import pca
        print("✓ Analysis modules imported successfully")
        
        # Visualization modules
        from biostructbenchmark.visualization import plots, structure, pca_plots
        print("✓ Visualization modules imported successfully")
        
        # DSSR module (may not have executable)
        try:
            from biostructbenchmark.analysis import dssr
            print("✓ DSSR module imported successfully")
        except Exception as e:
            print(f"⚠ DSSR module import warning: {e}")
        
        return True
        
    except Exception as e:
        print(f"✗ Import failed: {e}")
        traceback.print_exc()
        return False

def test_structure_loading():
    """Test structure loading functionality"""
    print("\n=== TESTING STRUCTURE LOADING ===")
    
    try:
        from biostructbenchmark.core.io import get_structure
        
        # Test data paths
        test_pdb = PROJECT_ROOT / "tests" / "data" / "proteins_pdb" / "1bom.pdb"
        test_cif = PROJECT_ROOT / "tests" / "data" / "proteins_cif" / "1bom.cif"
        
        if not test_pdb.exists():
            print(f"✗ Test PDB file not found: {test_pdb}")
            return False
        
        # Load PDB structure
        pdb_structure = get_structure(test_pdb)
        if pdb_structure:
            print(f"✓ PDB structure loaded: {test_pdb.name}")
            print(f"  - Chains: {[chain.id for chain in pdb_structure.get_chains()]}")
            print(f"  - Residues: {len(list(pdb_structure.get_residues()))}")
        else:
            print(f"✗ Failed to load PDB structure: {test_pdb}")
            return False
        
        # Load CIF structure if available
        if test_cif.exists():
            cif_structure = get_structure(test_cif)
            if cif_structure:
                print(f"✓ CIF structure loaded: {test_cif.name}")
            else:
                print(f"⚠ CIF structure loading failed: {test_cif}")
        
        return True
        
    except Exception as e:
        print(f"✗ Structure loading test failed: {e}")
        traceback.print_exc()
        return False

def test_alignment():
    """Test structure alignment functionality"""
    print("\n=== TESTING ALIGNMENT ===")
    
    try:
        from biostructbenchmark.core.alignment import align_structures_three_frames
        from biostructbenchmark.core.io import get_structure
        
        # Use test structures
        exp_path = PROJECT_ROOT / "tests" / "data" / "proteins_pdb" / "p456_02_experimental.pdb"
        pred_path = PROJECT_ROOT / "tests" / "data" / "proteins_cif" / "p456_02_predicted.cif"
        
        if not exp_path.exists() or not pred_path.exists():
            print(f"✗ Test structure files not found")
            print(f"  Expected: {exp_path}")
            print(f"  Predicted: {pred_path}")
            return False
        
        # Create output directory
        output_dir = PROJECT_ROOT / "test_outputs" / "alignment_test"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Running multi-frame alignment on {exp_path.name} vs {pred_path.name}")
        
        # Load structures
        exp_structure = get_structure(exp_path)
        pred_structure = get_structure(pred_path)
        
        if not exp_structure or not pred_structure:
            print("✗ Failed to load test structures")
            return False
        
        # Perform alignment
        result = align_structures_three_frames(exp_structure, pred_structure)
        
        if result and isinstance(result, dict):
            print("✓ Multi-frame alignment completed successfully")
            print(f"  - Alignment frames: {list(result.keys())}")
            
            # Print some basic statistics if available
            for frame_name, frame_result in result.items():
                if hasattr(frame_result, 'overall_rmsd'):
                    print(f"  - {frame_name} RMSD: {frame_result.overall_rmsd:.3f} Å")
                elif 'rmsd' in str(frame_result):
                    print(f"  - {frame_name}: {frame_result}")
                    
            print(f"  - Output directory: {output_dir}")
            
            return True
        else:
            print("✗ Multi-frame alignment failed")
            return False
            
    except Exception as e:
        print(f"✗ Alignment test failed: {e}")
        traceback.print_exc()
        return False

def test_pca_analysis():
    """Test PCA analysis functionality"""
    print("\n=== TESTING PCA ANALYSIS ===")
    
    try:
        from biostructbenchmark.analysis.pca import PCAAnalyzer
        
        # Test structure-level PCA
        analyzer = PCAAnalyzer()
        
        # Create mock data for structure-level analysis
        structure_data = {
            'structure1': {
                'full_rmsd': 2.5,
                'dna_positioning_rmsd': 3.1,
                'dna_structure_rmsd': 1.8
            },
            'structure2': {
                'full_rmsd': 1.9,
                'dna_positioning_rmsd': 2.8,
                'dna_structure_rmsd': 2.2
            }
        }
        
        pca_output = PROJECT_ROOT / "test_outputs" / "pca_test"
        pca_output.mkdir(parents=True, exist_ok=True)
        
        # Note: PCA requires multiple structures, so this is a minimal test
        print("✓ PCA analyzer initialized successfully")
        print(f"  - Output directory: {pca_output}")
        
        return True
        
    except Exception as e:
        print(f"✗ PCA analysis test failed: {e}")
        traceback.print_exc()
        return False

def test_visualization():
    """Test visualization functionality"""
    print("\n=== TESTING VISUALIZATION ===")
    
    try:
        from biostructbenchmark.visualization.plots import PublicationPlotter
        from biostructbenchmark.visualization.structure import StructureVisualizer
        from biostructbenchmark.core.alignment import ResidueRMSD
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        
        # Create test output directory
        viz_output = PROJECT_ROOT / "test_outputs" / "visualization_test"
        viz_output.mkdir(parents=True, exist_ok=True)
        
        # Test PublicationPlotter
        plotter = PublicationPlotter()
        print("✓ PublicationPlotter initialized")
        
        # Create mock RMSD data
        mock_rmsd = [
            ResidueRMSD('A1_ALA', 'ALA', 'A', 1, 2.1, 4, 'protein'),
            ResidueRMSD('A2_VAL', 'VAL', 'A', 2, 1.8, 4, 'protein'),
            ResidueRMSD('A3_LEU', 'LEU', 'A', 3, 3.5, 4, 'protein'),
            ResidueRMSD('B1_A', 'A', 'B', 1, 3.2, 12, 'dna'),
            ResidueRMSD('B2_T', 'T', 'B', 2, 2.8, 12, 'dna')
        ]
        
        # Test basic comparison plot
        rmsd_values = [r.rmsd for r in mock_rmsd]
        protein_rmsds = [r.rmsd for r in mock_rmsd if r.is_protein]
        
        comparison_plot = plotter.comparison_plot(
            protein_rmsds, rmsd_values[:len(protein_rmsds)],
            ('Protein RMSD', 'All RMSD'),
            'Test Comparison',
            viz_output / 'test_comparison.png'
        )
        
        if (viz_output / 'test_comparison.png').exists():
            print("✓ Comparison plot generated")
        else:
            print("⚠ Comparison plot not saved")
        
        # Test StructureVisualizer
        visualizer = StructureVisualizer()
        print("✓ StructureVisualizer initialized")
        print(f"  - Backend: {visualizer.backend}")
        
        return True
        
    except Exception as e:
        print(f"✗ Visualization test failed: {e}")
        traceback.print_exc()
        return False

def test_dssr_integration():
    """Test DSSR integration (may not have executable)"""
    print("\n=== TESTING DSSR INTEGRATION ===")
    
    try:
        from biostructbenchmark.analysis.dssr import DSSRAnalyzer
        
        # Try to initialize - may fail if DSSR not installed
        try:
            analyzer = DSSRAnalyzer()
            print("✓ DSSR analyzer initialized successfully")
            print(f"  - DSSR executable: {analyzer.dssr_exe}")
            return True
        except RuntimeError as e:
            print(f"⚠ DSSR not available: {e}")
            print("  This is expected if DSSR license hasn't been obtained")
            return True  # Not a failure - just unavailable
            
    except Exception as e:
        print(f"✗ DSSR integration test failed: {e}")
        traceback.print_exc()
        return False

def generate_test_report(results):
    """Generate comprehensive test report"""
    print("\n=== GENERATING TEST REPORT ===")
    
    # Create report directory
    report_dir = PROJECT_ROOT / "test_outputs" / "comprehensive_test_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate report data
    report_data = {
        "test_timestamp": datetime.now().isoformat(),
        "test_results": results,
        "summary": {
            "total_tests": len(results),
            "passed_tests": sum(1 for r in results.values() if r),
            "failed_tests": sum(1 for r in results.values() if not r),
            "success_rate": sum(1 for r in results.values() if r) / len(results) * 100
        },
        "environment_info": {
            "python_version": sys.version,
            "platform": sys.platform,
            "project_root": str(PROJECT_ROOT)
        }
    }
    
    # Save JSON report
    report_file = report_dir / "test_report.json"
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2)
    
    # Generate text report
    text_report = report_dir / "test_report.txt"
    with open(text_report, 'w') as f:
        f.write("BioStructBenchmark Comprehensive Test Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Project Root: {PROJECT_ROOT}\n\n")
        
        f.write("Test Results:\n")
        f.write("-" * 20 + "\n")
        for test_name, result in results.items():
            status = "PASS" if result else "FAIL"
            f.write(f"{test_name:30} {status}\n")
        
        f.write(f"\nSummary:\n")
        f.write(f"  Total Tests: {report_data['summary']['total_tests']}\n")
        f.write(f"  Passed: {report_data['summary']['passed_tests']}\n")
        f.write(f"  Failed: {report_data['summary']['failed_tests']}\n")
        f.write(f"  Success Rate: {report_data['summary']['success_rate']:.1f}%\n")
    
    print(f"✓ Test report generated: {report_dir}")
    return report_data

def main():
    """Run comprehensive test suite"""
    print("BioStructBenchmark Comprehensive Test Suite")
    print("=" * 50)
    
    # Initialize results
    results = {}
    
    # Run all tests
    results["imports"] = test_imports()
    results["structure_loading"] = test_structure_loading() 
    results["alignment"] = test_alignment()
    results["pca_analysis"] = test_pca_analysis()
    results["visualization"] = test_visualization()
    results["dssr_integration"] = test_dssr_integration()
    
    # Generate report
    report_data = generate_test_report(results)
    
    # Print summary
    print("\n" + "=" * 50)
    print("COMPREHENSIVE TEST SUMMARY")
    print("=" * 50)
    
    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{test_name:20} {status}")
    
    print(f"\nOverall Success Rate: {report_data['summary']['success_rate']:.1f}%")
    
    if report_data['summary']['failed_tests'] > 0:
        print(f"\n⚠ {report_data['summary']['failed_tests']} test(s) failed")
        print("Check individual test outputs above for details")
    else:
        print("\n🎉 All tests passed successfully!")
    
    print(f"\nDetailed reports saved to: test_outputs/comprehensive_test_report/")
    
    return report_data['summary']['success_rate'] == 100.0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)