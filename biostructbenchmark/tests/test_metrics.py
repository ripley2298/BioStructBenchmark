import pytest
from biostructbenchmark.core.metrics import decompose_error, calculate_statistics

class TestMetrics:
    def test_error_decomposition(self):
        """Test error decomposition into translation/orientation"""
        # Mock test - implement based on actual metrics.py
        result = decompose_error([1, 2, 3], [1, 2, 3])
        assert 'translation' in result
        assert 'orientation' in result
    
    def test_statistics_calculation(self):
        """Test basic statistics calculation"""
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        stats = calculate_statistics(data)
        
        assert 'mean' in stats
        assert 'std' in stats
        assert abs(stats['mean'] - 3.0) < 1e-10
