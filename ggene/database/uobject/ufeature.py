
from typing import Dict, Any

from ggene.database.uobject.uobject import UObject


class UFeature(UObject):
    """
    Common representation for all genomic features.

    Inherits namespace-style attribute access from BaseFeature.
    Adds genomic-specific functionality (overlaps, length, etc.).

    All attributes (core and extended) can be accessed directly.
    Missing attributes return None rather than raising AttributeError.
    """

    # Core fields that become instance attributes
    _core_fields = frozenset({
        'chrom', 'start', 'end', 'feature_type', 'source',
        'score', 'strand', 'name', 'id'
    })

    # Default values for core fields
    _defaults = {
        'chrom': "",
        'start': -1,
        'end': -1,
        'feature_type': "",
        'source': "",
        'score': None,
        'strand': "",
        'name': "",
        'id': "",
    }

    __slots__ = ('chrom', 'start', 'end', 'feature_type', 'source',
                 'score', 'strand', 'name', 'id', 'attributes', '_canonical')

    @property
    def chr(self):
        """Alias for chrom"""
        return self.chrom

    @property
    def length(self):
        """Length of the feature in base pairs"""
        return self.end - self.start

    @property
    def loc(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    def overlaps(self, start: int, end: int) -> bool:
        """Check if feature overlaps with given range."""
        return not (self.end < start or self.start > end)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for JSON serialization.

        Base class automatically includes 'chr' and 'length' properties.
        We just add the custom aliases for backwards compatibility.
        """
        result = super().to_dict()

        # Add aliases for backwards compatibility
        result['type'] = self.feature_type  # Alias for feature_type
        result['gene_id'] = self.id          # Alias for id

        return result

    def short_repr(self):
        """Ultra-compact representation"""
        loc = f"{self.chrom}:{self.start}-{self.end}"
        parts = [f"UFeature({self.feature_type!r}, {loc}"]

        if self.strand:
            parts.append(f", strand={self.strand!r}")
        if self.name:
            parts.append(f", name={self.name!r}")

        # Show non-empty attributes count
        n_attrs = len(self.attributes)
        if n_attrs > 0:
            parts.append(f", +{n_attrs} attrs")

        parts.append(")")
        return "".join(parts)
    
    def long_repr(self):
        
        outstrs = []
        outstrs.append(f"{self.name}, {self.feature_type}, {self.id}")
        outstrs.append(f"  {self.loc}")
        dd = self.to_dict()
        atts = dd.pop("attributes",{})
        for k, v in dd.items():
            outstrs.append(f"  {k}: {v}")
        
        outstrs.append("  Attributes:")
        for k, v in atts.items():
            outstrs.append(f"    {k}: {v}")
        
        return "\n".join(outstrs)    
    
    def __repr__(self) -> str:
        """Compact repr showing key fields and non-empty attributes."""
        # Core info
        loc = f"{self.chrom}:{self.start}-{self.end}"
        parts = [f"UFeature({self.name}, {self.feature_type}, {loc} ({self.length}bp)"]

        if self.strand:
            parts.append(f", strand={self.strand!r}")
        if self.name:
            parts.append(f", name={self.name!r}")

        for att, val in self.attributes.items():
            if att.startswith("_"):
                continue
            if val:
                parts.append(f", {att}={val}")

        parts.append(")")
        return "".join(parts)

    def __eq__(self, other):
        """Genomic feature equality based on location and type"""
        return (other.start == self.start and other.end == self.end and
                other.feature_type == self.feature_type and other.strand == self.strand)

    def __hash__(self):
        """Hash for genomic features"""
        return hash((self.feature_type, self.start, self.end, self.strand))

    def __lt__(self, other):
        """For heap-based merging and sorting by genomic position"""
        return (self.chrom, self.start, self.end) < (other.chrom, other.start, other.end)


# Alias for backwards compatibility
UFeat = UFeature
