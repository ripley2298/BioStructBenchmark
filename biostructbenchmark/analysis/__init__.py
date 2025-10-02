"""
biostructbenchmark.analysis
Analysis modules for structure comparison
"""

# Import analysis classes
from .bfactor import BFactorAnalyzer, BFactorComparison, BFactorStatistics
from .secondary import SecondaryStructureAnalyzer, SecondaryStructure
from .consensus import ConsensusAnalyzer, ConsensusError
from .mutations import MutationAnalyzer, Mutation
from .hbond import HBondAnalyzer, HydrogenBond, HBondComparison, HBondStatistics


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

    # Hydrogen bond analysis
    'HBondAnalyzer',
    'HydrogenBond',
    'HBondComparison',
    'HBondStatistics',
]

