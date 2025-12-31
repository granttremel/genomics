

from typing import Iterator, Dict, Any, List, Optional, Tuple, Union, Callable, TypeVar, TYPE_CHECKING

if TYPE_CHECKING:
    from ggene.database.annotations import ColumnSpec, DerivedSpec


class UFeature:
    """Common representation for all genomic features.

    Supports namespace-style attribute access: feature.gene_biotype
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

    def __init__(self, data: Optional[Dict[str, Any]] = None,
                 canonical: Optional[set] = None, **kwargs):
        """Initialize from a flat dict. Core fields are extracted, rest go to attributes.

        Args:
            data: Flat dict of all fields (core + extended). Core fields are
                  extracted to instance attributes, others go to self.attributes.
            canonical: Optional set of "expected" attribute names for this feature type.
                       Used to distinguish typos from intentionally-checked-but-empty.
            **kwargs: Alternative to data dict, for backwards compatibility.
        """
        # Merge data and kwargs
        all_data = {**(data or {}), **kwargs}

        # Initialize attributes dict first
        object.__setattr__(self, 'attributes', {})
        object.__setattr__(self, '_canonical', canonical or set())

        # Set core fields with defaults
        for field in self._core_fields:
            value = all_data.pop(field, self._defaults.get(field))
            object.__setattr__(self, field, value)

        # Everything else goes to attributes (skip empty values)
        for key, value in all_data.items():
            if key != 'attributes' and value not in (None, "", []):
                self.attributes[key] = value

        # Merge any pre-existing attributes dict
        if 'attributes' in (data or {}):
            for k, v in data['attributes'].items():
                if v not in (None, "", []):
                    self.attributes[k] = v
    
    
    @classmethod
    def from_parsed(cls, data: Dict[str, Any],
                    columns: Optional[List['ColumnSpec']] = None,
                    derived: Optional[List['DerivedSpec']] = None, raw_line:str = "") -> 'UFeature':
        """Create UFeature from parsed column data with derived attributes.

        This is the main factory method used by TabularStream.
        """
        # Track canonical attribute names from columns
        canonical = set()
        if columns:
            canonical = {col.name for col in columns if hasattr(col, "name")}

        # Create feature with flat data
        feature = cls(data, canonical=canonical)

        # Compute derived attributes (they can access all fields via namespace)
        if derived:
            for spec in derived:
                result = spec.evaluate(feature)
                
                if result is not None:
                    # Derived can override core fields too
                    if spec.name in cls._core_fields:
                        prev = getattr(feature, spec.name)
                        object.__setattr__(feature, spec.name, result)
                        curr = getattr(feature, spec.name)
                    else:
                        feature.attributes[spec.name] = result
                    canonical.add(spec.name)

        # feature.attributes["_raw_line"] = raw_line

        return feature

    @property
    def chr(self):
        return self.chrom

    @property
    def length(self):
        return self.end - self.start

    def __getattr__(self, name: str) -> Any:
        """Return attribute value, or None if not present.

        This allows namespace-style access without AttributeError for missing attrs.
        """
        # Avoid recursion on internal attributes
        if name.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")

        # Try attributes dict, return None if not found
        try:
            attrs = object.__getattribute__(self, 'attributes')
            return attrs.get(name)  # Returns None if missing
        except AttributeError:
            return None

    def __setattr__(self, name: str, value: Any) -> None:
        """Set attribute, routing non-core fields to attributes dict."""
        if name.startswith('_') or name in self._core_fields or name == 'attributes':
            object.__setattr__(self, name, value)
        else:
            # Route to attributes dict
            attrs = object.__getattribute__(self, 'attributes')
            if value not in (None, "", []):
                attrs[name] = value
            elif name in attrs:
                del attrs[name]  # Remove empty values

    def __lt__(self, other):
        """For heap-based merging."""
        return (self.chrom, self.start, self.end) < (other.chrom, other.start, other.end)

    def overlaps(self, start: int, end: int) -> bool:
        """Check if feature overlaps with given range."""
        return not (self.end < start or self.start > end)

    def get(self, att, default=None):
        """Get attribute by name, checking core fields then attributes dict."""
        if att in self._core_fields:
            return getattr(self, att, default)
        return self.attributes.get(att, default)

    def __getitem__(self, att):
        if isinstance(att, str):
            return self.get(att)
        else:
            raise IndexError

    def __eq__(self, other):
        return (other.start == self.start and other.end == self.end and
                other.feature_type == self.feature_type and other.strand == self.strand)

    def __hash__(self):
        return hash((self.feature_type, self.start, self.end, self.strand))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'type': self.feature_type,
            'source': self.source,
            'score': self.score,
            'strand': self.strand,
            'name': self.name,
            'id': self.id,
            'chr': self.chrom,
            'gene_id': self.id,
            'attributes': self.attributes
        }

    def has_attr(self, name: str) -> bool:
        """Check if attribute exists and has a non-empty value."""
        if name in self._core_fields:
            val = getattr(self, name)
            return val not in (None, "", -1)
        return name in self.attributes

    def is_canonical(self, name: str) -> bool:
        """Check if attribute name was defined in the stream's column specs."""
        return name in self._canonical

    def short_repr(self):
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

    def __repr__(self) -> str:
        """Compact repr showing key fields and non-empty attributes."""
        # Core info
        loc = f"{self.chrom}:{self.start}-{self.end}"
        parts = [f"UFeature({self.feature_type}, {loc}"]

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

UFeat = UFeature

