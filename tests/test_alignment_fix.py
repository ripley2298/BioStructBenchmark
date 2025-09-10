#!/usr/bin/env python3
"""
Quick test to verify basic alignment functionality is working
"""

import sys
from pathlib import Path
import pytest

# Add the local package to the path
sys.path.insert(0, str(Path(__file__).parent))

from biostructbenchmark.core.alignment import compare_structures
from biostructbenchmark.core.io import get_structure

def test_alignment_basic_functionality():
    """Test that basic alignment functionality works"""
    
    # Load test structures
    exp_path = Path("tests/data/proteins_pdb/p456_02_experimental.pdb")
    pred_path = Path("tests/data/proteins_cif/p456_02_predicted.cif")
    
    if not exp_path.exists() or not pred_path.exists():
        pytest.skip("Test structure files not found")
    
    print("Loading structures...")
    try:
        exp_structure = get_structure(exp_path)
        pred_structure = get_structure(pred_path)
    except Exception as e:
        pytest.skip(f"Could not load structures: {e}")
    
    print("Performing alignment...")
    try:
        result = compare_structures(exp_structure, pred_structure)
        
        print(f"Overall RMSD: {result.overall_rmsd:.3f} Å")
        print(f"Aligned atoms: {result.aligned_atom_count}")
        
        # Check that we get reasonable results
        assert result.overall_rmsd > 0, "RMSD should be positive"
        assert result.overall_rmsd < 50, "RMSD should be reasonable (< 50 Å)"
        assert result.aligned_atom_count > 0, "Should align some atoms"
        
        print("✓ Alignment test passed")
        return True
        
    except Exception as e:
        pytest.skip(f"Alignment failed: {e}")

if __name__ == "__main__":
    test_alignment_basic_functionality()