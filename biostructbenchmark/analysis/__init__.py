"""
biostructbenchmark.analysis
Analysis modules for structure comparison
"""

# Import analysis classes
from .bfactor import BFactorAnalyzer, BFactorComparison, BFactorStatistics
from .secondary import SecondaryStructureAnalyzer, SecondaryStructure
from .consensus import ConsensusAnalyzer, ConsensusError
from .mutations import MutationAnalyzer, Mutation
from .pca import PCAAnalyzer, PCAResult, StructureOutlier, ResidueOutlier

# Conditional import for CURVES+ (may not have executable)
try:
    from .curves import CurvesAnalyzer, CurvesParameters, HydrogenBond
    _curves_available = True
except (ImportError, RuntimeError):
    _curves_available = False
    CurvesAnalyzer = None
    CurvesParameters = None
    HydrogenBond = None

__all__ = [
    # B-factor analysis
    'BFactorAnalyzer',
    'BFactorComparison', 
    'BFactorStatistics',
    
    # Secondary structure
    'SecondaryStructureAnalyzer',
    'SecondaryStructure',
    
    # Consensus errors
    'ConsensusAnalyzer',
    'ConsensusError',
    
    # Mutations
    'MutationAnalyzer',
    'Mutation',
    
    # PCA analysis
    'PCAAnalyzer',
    'PCAResult',
    'StructureOutlier',
    'ResidueOutlier',
]

# Add CURVES+ if available
if _curves_available:
    __all__.extend([
        'CurvesAnalyzer',
        'CurvesParameters',
        'HydrogenBond'
    ])