import pytest
import numpy as np
from pathlib import Path
from biostructbenchmark.core.alignment import calculate_rmsd, align_structures

class TestAlignment:
    def test_calculate_rmsd(self):
        """Test basic RMSD calculation"""
        coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        coords2 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        
        rmsd = calculate_rmsd(coords1, coords2)
        assert abs(rmsd) < 1e-10  # Should be zero for identical coords
    
    def test_rmsd_translation(self):
        """Test RMSD with translation"""
        coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        coords2 = np.array([[1, 1, 1], [2, 1, 1], [1, 2, 1]])  # Translated by (1,1,1)
        
        rmsd = calculate_rmsd(coords1, coords2)
        expected = np.sqrt(3)  # sqrt(1^2 + 1^2 + 1^2)
        assert abs(rmsd - expected) < 1e-10
