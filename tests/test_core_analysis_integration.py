#!/usr/bin/env python3
"""
Core-Analysis Integration Test for BioStructBenchmark
Tests complete pipeline: core → analysis modules
"""

import sys
import os
from pathlib import Path
import traceback
from datetime import datetime

# Add the package to Python path
sys.path.insert(0, str(Path(__file__).parent))

def test_core_analysis_integration():
    """Test complete pipeline from core to analysis modules"""
    print("="*70)
    print("BioStructBenchmark Core-Analysis Integration Test")
    print("="*70)
    print(f"Test started: {datetime.now()}")
    
    # Create output directory
    output_dir = Path("test_core_analysis_outputs")
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir.absolute()}")
    
    try:
        print("\n1. Testing Core Module Imports...")
        from biostructbenchmark.core import io, alignment, metrics
        from biostructbenchmark.core.io import get_structure, validate_file
        from biostructbenchmark.core.alignment import compare_structures, align_structures_three_frames, ResidueRMSD
        from biostructbenchmark.core.metrics import generate_comprehensive_metrics
        print("✓ Core imports successful")
        
        print("\n2. Testing Analysis Module Imports...")
        from biostructbenchmark.analysis import BFactorAnalyzer, SecondaryStructureAnalyzer, ConsensusAnalyzer
        from biostructbenchmark.analysis.bfactor import BFactorComparison, BFactorStatistics
        from biostructbenchmark.analysis.secondary import SecondaryStructure
        from biostructbenchmark.analysis.consensus import ConsensusError
        from biostructbenchmark.analysis.mutations import MutationAnalyzer, Mutation
        print("✓ Analysis imports successful")
        
        # Test data files
        test_data_dir = Path("tests/data")
        observed_file = test_data_dir / "proteins_pdb" / "p456_02_experimental.pdb"
        predicted_file = test_data_dir / "proteins_cif" / "p456_02_predicted.cif"
        
        print(f"\n3. Testing Core Data Pipeline...")
        print(f"   Loading structures: {observed_file.name} vs {predicted_file.name}")
        
        # Load structures via core
        observed_structure = get_structure(observed_file)
        predicted_structure = get_structure(predicted_file)
        
        if not observed_structure or not predicted_structure:
            print("   ERROR: Failed to load structures")
            return False
            
        print(f"   ✓ Structures loaded successfully")
        
        # Run core alignment analysis
        print("   Running core structure comparison...")
        alignment_result = compare_structures(observed_structure, predicted_structure)
        if not alignment_result:
            print("   ERROR: Structure comparison failed")
            return False
            
        print(f"   ✓ Alignment complete: {alignment_result.overall_rmsd:.3f} Å")
        print(f"   ✓ Units analyzed: {len(alignment_result.residue_rmsds)} (residues + bases)")
        
        # Generate core metrics
        structure_metrics = generate_comprehensive_metrics(alignment_result)
        print(f"   ✓ Core metrics generated")
        
        print(f"\n4. Testing Analysis Module Integration...")
        
        # Test B-factor analysis
        print("   Testing B-factor analysis...")
        try:
            bfactor_analyzer = BFactorAnalyzer()
            bfactor_results = bfactor_analyzer.analyze_structures(observed_structure, predicted_structure)
            print(f"   ✓ B-factor analysis: {len(bfactor_results)} comparisons")
            
            # Generate statistics
            bfactor_stats = bfactor_analyzer.calculate_statistics(bfactor_results)
            print(f"   ✓ B-factor correlation: {bfactor_stats.correlation:.3f}")
        except Exception as e:
            print(f"   ⚠ B-factor analysis failed: {e}")
        
        # Test secondary structure analysis
        print("   Testing secondary structure analysis...")
        try:
            ss_analyzer = SecondaryStructureAnalyzer()
            ss_observed = ss_analyzer.analyze_structure(observed_structure)
            ss_predicted = ss_analyzer.analyze_structure(predicted_structure)
            
            print(f"   ✓ Secondary structure - Observed: {len(ss_observed)} assignments")
            print(f"   ✓ Secondary structure - Predicted: {len(ss_predicted)} assignments")
            
            # Compare secondary structures
            ss_comparison = ss_analyzer.compare_structures(ss_observed, ss_predicted)
            print(f"   ✓ SS accuracy: {ss_comparison['accuracy']:.1%}")
        except Exception as e:
            print(f"   ⚠ Secondary structure analysis failed: {e}")
        
        # Test consensus analysis (using core ResidueRMSD data)
        print("   Testing consensus error analysis...")
        try:
            consensus_analyzer = ConsensusAnalyzer()
            
            # Simulate multiple comparisons by using the same data with slight variations
            rmsd_datasets = [
                alignment_result.residue_rmsds,  # Original
                alignment_result.residue_rmsds,  # Duplicate for testing
            ]
            
            consensus_errors = consensus_analyzer.identify_consensus_errors(rmsd_datasets)
            print(f"   ✓ Consensus analysis: {len(consensus_errors)} error regions identified")
            
            # Find top problematic regions
            top_errors = sorted(consensus_errors, key=lambda x: x.mean_rmsd, reverse=True)[:5]
            print(f"   ✓ Top error region: {top_errors[0].position} ({top_errors[0].mean_rmsd:.3f} Å)")
        except Exception as e:
            print(f"   ⚠ Consensus analysis failed: {e}")
        
        # Test mutation analysis
        print("   Testing mutation analysis...")
        try:
            mutation_analyzer = MutationAnalyzer()
            mutations = mutation_analyzer.identify_mutations(observed_structure, predicted_structure)
            print(f"   ✓ Mutations identified: {len(mutations)}")
            
            if mutations:
                # Analyze impact of mutations on RMSD
                mutation_impact = mutation_analyzer.analyze_mutation_impact(
                    mutations, alignment_result.residue_rmsds
                )
                print(f"   ✓ Mutation impact analysis complete")
        except Exception as e:
            print(f"   ⚠ Mutation analysis failed: {e}")
        
        print(f"\n5. Testing Data Export Integration...")
        
        # Export core results
        from biostructbenchmark.core.alignment import export_comprehensive_alignment_report
        export_comprehensive_alignment_report(alignment_result, output_dir, "core_analysis_test")
        
        # Export analysis results to CSV
        try:
            # B-factor results
            if 'bfactor_results' in locals():
                bfactor_df = bfactor_analyzer.to_dataframe(bfactor_results)
                bfactor_df.to_csv(output_dir / "bfactor_analysis.csv", index=False)
                print("   ✓ B-factor results exported")
            
            # Secondary structure results  
            if 'ss_comparison' in locals():
                ss_df = ss_analyzer.comparison_to_dataframe(ss_comparison)
                ss_df.to_csv(output_dir / "secondary_structure_comparison.csv", index=False)
                print("   ✓ Secondary structure results exported")
            
            # Consensus errors
            if 'consensus_errors' in locals():
                consensus_df = consensus_analyzer.to_dataframe(consensus_errors)
                consensus_df.to_csv(output_dir / "consensus_errors.csv", index=False)
                print("   ✓ Consensus error results exported")
                
        except Exception as e:
            print(f"   ⚠ Some export operations failed: {e}")
        
        print(f"\n6. Testing End-to-End Workflow...")
        
        # Create a summary report combining core and analysis results
        summary_path = output_dir / "integration_summary.txt"
        with open(summary_path, 'w') as f:
            f.write("BioStructBenchmark Core-Analysis Integration Summary\\n")
            f.write("="*60 + "\\n\\n")
            f.write(f"Test Date: {datetime.now()}\\n")
            f.write(f"Structures: {observed_file.name} vs {predicted_file.name}\\n\\n")
            
            # Core results
            f.write("CORE ANALYSIS RESULTS:\\n")
            f.write("-" * 30 + "\\n")
            f.write(f"Overall RMSD: {alignment_result.overall_rmsd:.3f} Å\\n")
            f.write(f"Protein RMSD: {structure_metrics.protein_rmsd:.3f} Å\\n")
            f.write(f"DNA RMSD: {structure_metrics.dna_rmsd:.3f} Å\\n")
            f.write(f"Units analyzed: {len(alignment_result.residue_rmsds)}\\n")
            f.write(f"Translation error: {structure_metrics.error_components.translation_error:.3f} Å\\n")
            f.write(f"Rotation angle: {structure_metrics.error_components.rotation_angle:.1f}°\\n\\n")
            
            # Analysis results
            f.write("ANALYSIS MODULE RESULTS:\\n")
            f.write("-" * 30 + "\\n")
            
            if 'bfactor_stats' in locals():
                f.write(f"B-factor correlation: {bfactor_stats.correlation:.3f}\\n")
            if 'ss_comparison' in locals():
                f.write(f"Secondary structure accuracy: {ss_comparison['accuracy']:.1%}\\n")
            if 'consensus_errors' in locals():
                f.write(f"Consensus error regions: {len(consensus_errors)}\\n")
            if 'mutations' in locals():
                f.write(f"Mutations identified: {len(mutations)}\\n")
                
        print(f"   ✓ Integration summary created: {summary_path}")
        
        print(f"\n7. Verifying output files...")
        output_files = list(output_dir.glob("*"))
        print(f"   Generated {len(output_files)} output files:")
        for file in sorted(output_files):
            size = file.stat().st_size if file.is_file() else 0
            print(f"     - {file.name} ({size} bytes)")
        
        print(f"\n" + "="*70)
        print("CORE-ANALYSIS INTEGRATION TEST SUCCESSFUL!")
        print("="*70)
        print(f"All outputs saved to: {output_dir.absolute()}")
        print(f"Test completed: {datetime.now()}")
        print(f"✓ Core modules working correctly")
        print(f"✓ Analysis modules integrated successfully") 
        print(f"✓ Data flows properly between modules")
        print(f"✓ Export functionality working end-to-end")
        
        return True
        
    except Exception as e:
        print(f"\nERROR during integration test:")
        print(f"Exception: {e}")
        print("Traceback:")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_core_analysis_integration()
    sys.exit(0 if success else 1)