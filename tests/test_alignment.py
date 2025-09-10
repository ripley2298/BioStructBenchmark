"""
tests/test_alignment.py
Fixed test suite for alignment module
"""

import pytest
import numpy as np
from pathlib import Path
from biostructbenchmark.core.alignment import export_residue_rmsd_csv, ResidueRMSD


class TestAlignment:
   
    # For troubleshooting issues with csv export in alignment.py, may be removed after issue is addressed.
    def test_export_residue_rmsd_csv(self, tmp_path):
        """Test CSV export functionality"""
        
        residue_rmsds = [
            ResidueRMSD("ALA_10", "ALA", "A", 10, 2.5, 5, "protein"),
            ResidueRMSD("DG_5", "DG", "B", 5, 3.2, 10, "dna")
        ]
        
        output_file = tmp_path / "test_rmsd.csv"
        export_residue_rmsd_csv(residue_rmsds, output_file)
        
        assert output_file.exists()
        
        # TEMPORARY DEBUG CODE - ONLY FOR DEBUGGING
        import csv
        print("=== DEBUG: CSV Content ===")
        with open(output_file, 'r') as f:
            content = f.read()
            print(content)
        print("==========================")
