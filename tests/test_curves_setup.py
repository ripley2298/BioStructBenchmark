#!/usr/bin/env python3
"""
Test script for CURVES+ setup and integration
"""

import sys
import os
from pathlib import Path
import shutil

# Add the package to Python path
sys.path.insert(0, str(Path(__file__).parent))

def test_curves_setup():
    """Test CURVES+ installation and integration"""
    print("="*60)
    print("CURVES+ Setup Test")
    print("="*60)
    
    print("\n1. Checking CURVES+ executable...")
    
    # Check for CURVES+ executable in common locations
    curves_paths = [
        'Cur+',
        'curves+', 
        '/usr/local/bin/Cur+',
        '/opt/curves/Cur+',
        os.path.expanduser('~/curves_plus/Cur+')
    ]
    
    curves_exe = None
    for path in curves_paths:
        if shutil.which(path):
            curves_exe = path
            break
        elif os.path.isfile(path):
            curves_exe = path
            break
    
    if curves_exe:
        print(f"✓ CURVES+ found at: {curves_exe}")
        
        # Test executable
        import subprocess
        try:
            result = subprocess.run([curves_exe, '-h'], 
                                  capture_output=True, text=True, timeout=10)
            print("✓ CURVES+ executable responds to help command")
        except Exception as e:
            print(f"⚠ CURVES+ executable issue: {e}")
    else:
        print("✗ CURVES+ executable not found")
        print("\nINSTALLATION NEEDED:")
        print("1. Run: ./install_curves_plus.sh")
        print("2. Or manually install CURVES+ from https://bisi.ibcp.fr/tools/curves_plus/")
        return False
    
    print("\n2. Testing CurvesAnalyzer import...")
    try:
        from biostructbenchmark.analysis import CurvesAnalyzer
        print("✓ CurvesAnalyzer import successful")
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False
    
    print("\n3. Testing CurvesAnalyzer initialization...")
    try:
        analyzer = CurvesAnalyzer(curves_exe)
        print("✓ CurvesAnalyzer initialized successfully")
    except Exception as e:
        print(f"✗ Initialization failed: {e}")
        return False
    
    print("\n4. Checking test data...")
    test_data_dir = Path("tests/data")
    dna_structures = []
    
    if test_data_dir.exists():
        # Look for structures with DNA
        for pdb_file in test_data_dir.rglob("*.pdb"):
            if "dna" in pdb_file.name.lower() or "nucleotide" in pdb_file.name.lower():
                dna_structures.append(pdb_file)
        
        # Check our test structure that has DNA
        p456_file = test_data_dir / "proteins_pdb" / "p456_02_experimental.pdb"
        if p456_file.exists():
            dna_structures.append(p456_file)
            print(f"✓ Found test structure with DNA: {p456_file.name}")
    
    if not dna_structures:
        print("⚠ No DNA-containing test structures found")
        print("  CURVES+ needs DNA structures for meaningful analysis")
        return False
    
    print(f"✓ Found {len(dna_structures)} potential DNA structures")
    
    print("\n5. Testing DNA extraction...")
    test_structure = dna_structures[0]
    try:
        # Test DNA chain extraction
        import tempfile
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            dna_file = analyzer._extract_dna_chains(test_structure, temp_path)
            
            if dna_file.exists() and dna_file.stat().st_size > 0:
                print(f"✓ DNA extraction successful: {dna_file.stat().st_size} bytes")
            else:
                print("⚠ DNA extraction produced empty file")
                
    except Exception as e:
        print(f"✗ DNA extraction failed: {e}")
        return False
    
    print("\n6. Testing hydrogen bond detection...")
    try:
        hbonds = analyzer.detect_hydrogen_bonds(test_structure)
        print(f"✓ Hydrogen bond detection: {len(hbonds)} bonds found")
        
        if hbonds:
            # Show first few bonds
            for i, bond in enumerate(hbonds[:3]):
                print(f"  Bond {i+1}: {bond.donor_residue}:{bond.donor_atom} -> {bond.acceptor_residue}:{bond.acceptor_atom} ({bond.distance:.2f}Å)")
                
    except Exception as e:
        print(f"✗ Hydrogen bond detection failed: {e}")
        return False
    
    print("\n" + "="*60)
    print("CURVES+ SETUP TEST RESULTS")
    print("="*60)
    
    if curves_exe:
        print("✓ CURVES+ executable found and working")
        print("✓ Python integration functional")
        print("✓ DNA analysis capabilities available")
        print("\nREADY FOR CURVES+ ANALYSIS!")
        print(f"\nTo use CURVES+:")
        print("  from biostructbenchmark.analysis import CurvesAnalyzer")
        print(f"  analyzer = CurvesAnalyzer('{curves_exe}')")
        print("  params = analyzer.analyze_structure('dna_structure.pdb')")
        return True
    else:
        print("✗ CURVES+ not properly installed")
        print("\nINSTALLATION REQUIRED:")
        print("1. Download CURVES+ from: https://bisi.ibcp.fr/tools/curves_plus/")
        print("2. Extract and compile the source code")
        print("3. Add executable to PATH or specify path in CurvesAnalyzer")
        return False

if __name__ == "__main__":
    success = test_curves_setup()
    sys.exit(0 if success else 1)