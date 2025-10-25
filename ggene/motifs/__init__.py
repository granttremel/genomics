"""Motif detection and analysis module for genomic sequences.

This module provides classes and functions for detecting various types of 
sequence motifs including patterns, repeats, structural features, and 
machine learning-based motifs.
"""

# Core motif classes
from .motif import BaseMotif, MotifDetector, IndexedMotifDetector

# Specific motif types
from .pattern import PatternMotif, string_to_re, repeat_to_re
from .repeat import RepeatMotif
from .metric import MetricMotif, gcpattern
from .mlmotif import MLMotif
from .structural import StructuralMotif
from .pwm import PWM

# Utilities
from .motifio import MotifIO

# Pre-defined motifs from pattern.py
from .pattern import (
    splice_donor_motif,
    splice_branch_motif, 
    splice_acceptor_motif
)

__all__ = [
    # Core classes
    'BaseMotif',
    'MotifDetector', 
    'IndexedMotifDetector',
    
    # Motif types
    'PatternMotif',
    'RepeatMotif',
    'MetricMotif',
    'MLMotif',
    'StructuralMotif',
    'PWM',
    
    # Utilities
    'MotifIO',
    'string_to_re',
    'repeat_to_re',
    
    # Pre-defined motifs
    'gcpattern',
    'splice_donor_motif',
    'splice_branch_motif',
    'splice_acceptor_motif',
]
