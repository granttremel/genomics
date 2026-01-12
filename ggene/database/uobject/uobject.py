"""
Base class for feature-like objects with namespace attribute access

Provides:
- Namespace-style attribute access (feature.attr returns None instead of AttributeError)
- Separation of core fields (instance attributes) vs extended attributes (dict)
- Factory method from parsed data
- Extensible repr/to_dict
"""

from typing import Dict, Any, List, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:
    from ggene.database.annotations import ColumnSpec, DerivedSpec


class UObject:
    """
    Base class for feature-like objects with flexible attribute access.

    Subclasses should define:
        _core_fields: frozenset of field names to store as instance attributes
        _defaults: dict of default values for core fields
        __slots__: tuple including core fields + 'attributes' + '_canonical'

    Features:
        - Missing attributes return None instead of raising AttributeError
        - Core fields stored as instance attributes (fast, typed)
        - Extended fields stored in self.attributes dict (flexible)
        - Factory method from_parsed() for column specs
    """

    # Subclasses MUST override these
    _core_fields: frozenset = frozenset()
    _defaults: Dict[str, Any] = {}

    # Subclasses SHOULD override __slots__ to include their core fields
    __slots__ = ('attributes', '_canonical')

    def __init__(self, data: Optional[Dict[str, Any]] = None,
                 canonical: Optional[Set[str]] = None, **kwargs):
        """
        Initialize from a flat dict.

        Core fields are extracted to instance attributes, rest go to self.attributes.

        Args:
            data: Flat dict of all fields (core + extended)
            canonical: Optional set of "expected" attribute names for this feature type
            **kwargs: Alternative to data dict
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
            if key != 'attributes' and value is not None:
                self.attributes[key] = value

        # Merge any pre-existing attributes dict
        if 'attributes' in (data or {}):
            for k, v in data['attributes'].items():
                if v is not None:
                    self.attributes[k] = v

    @classmethod
    def from_parsed(cls, data: Dict[str, Any],
                    columns: Optional[List['ColumnSpec']] = None,
                    derived: Optional[List['DerivedSpec']] = None):
        """
        Create feature from parsed column data with derived attributes.

        This is the main factory method used by TabularStream.

        Args:
            data: Parsed column data
            columns: Column specifications (for tracking canonical names)
            derived: Derived attribute specifications

        Returns:
            Feature instance
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
                        object.__setattr__(feature, spec.name, result)
                    else:
                        feature.attributes[spec.name] = result
                    canonical.add(spec.name)

        return feature

    def purge(self, atts):
        """
        Remove and return specified attributes.

        Args:
            atts: List of attribute names or specs

        Returns:
            Dict of removed attribute values
        """
        res = {}

        for att in atts:
            if hasattr(att, 'name'):
                name = att.name
            else:
                name = att

            if name in self.attributes:
                res[name] = self.attributes[name]
                del self.attributes[name]

        return res

    def __getattr__(self, name: str) -> Any:
        """
        Return attribute value, or None if not present.

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

    def get(self, att: str, default=None):
        """Get attribute by name, checking core fields then attributes dict."""
        if att in self._core_fields:
            return getattr(self, att, default)
        return self.attributes.get(att, default)

    def __getitem__(self, att: str):
        """Dict-style access to attributes."""
        if isinstance(att, str):
            return self.get(att)
        else:
            raise IndexError(f"Invalid index type: {type(att)}")

    def has_attr(self, name: str) -> bool:
        """Check if attribute exists and has a non-empty value."""
        if name in self._core_fields:
            val = getattr(self, name, None)
            return val not in (None, "", -1)
        return name in self.attributes

    def is_canonical(self, name: str) -> bool:
        """Check if attribute name was defined in the stream's column specs."""
        return name in self._canonical

    def get_core_fields(self) -> Dict[str, Any]:
        """Get dict of all core field values."""
        return {field: getattr(self, field) for field in self._core_fields}

    def get_all_fields(self) -> Dict[str, Any]:
        """Get dict of all fields (core + attributes)."""
        result = self.get_core_fields()
        result.update(self.attributes)
        return result

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for JSON serialization.

        Automatically includes all properties defined in the class.
        Subclasses can override to customize output.

        Default implementation returns core fields + properties + attributes.
        """
        result = self.get_core_fields()

        # Add all properties
        for name in dir(type(self)):
            # Skip private attributes and special methods
            if name.startswith('_'):
                continue

            # Get the class attribute (descriptor)
            try:
                attr = getattr(type(self), name)
            except AttributeError:
                continue

            # Check if it's a property
            if isinstance(attr, property):
                # Skip if already in core fields (avoid duplicates)
                if name in self._core_fields:
                    continue

                # Call the property to get its value
                try:
                    value = getattr(self, name)
                    # Only include non-None values
                    if value is not None:
                        result[name] = value
                except Exception:
                    # Skip properties that raise errors
                    pass

        # Add attributes dict
        result['attributes'] = self.attributes.copy()

        return result

    def __repr__(self) -> str:
        """
        Compact repr showing class name and key fields.

        Subclasses can override to customize representation.
        """
        class_name = type(self).__name__

        # Show core fields
        core_parts = []
        # for field in sorted(self._core_fields):
        for field in self._core_fields:
            value = getattr(self, field, None)
            if value not in (None, "", -1):
                if isinstance(value, str):
                    core_parts.append(f"{field}={value!r}")
                else:
                    core_parts.append(f"{field}={value}")

        att_parts = []
        for field in ['description']:
            value = getattr(self, field, None)
            if value not in (None, "", -1):
                if isinstance(value, str):
                    att_parts.append(f"{field}={value!r}")
                else:
                    att_parts.append(f"{field}={value}")

        # Show attributes count
        n_attrs = len(self.attributes)
        if n_attrs > 0:
            core_parts.append(f"+{n_attrs} attrs")
        
        return f"{class_name}({', '.join(core_parts)}, {', '.join(att_parts)})"

    def short_repr(self) -> str:
        """
        Ultra-compact representation.

        Subclasses should override for custom short format.
        """
        class_name = type(self).__name__
        n_attrs = len(self.attributes)
        return f"{class_name}(+{n_attrs} attrs)"



    def __eq__(self, other):
        """
        Equality comparison.

        Default compares all core fields. Subclasses should override
        if different equality semantics are needed.
        """
        if not isinstance(other, self.__class__):
            return False

        for field in self._core_fields:
            if getattr(self, field) != getattr(other, field):
                return False

        return True

    def __hash__(self):
        """
        Hash based on core fields.

        Subclasses should override if different hash semantics are needed.
        """
        return hash(tuple(getattr(self, field) for field in sorted(self._core_fields)))

    def __lt__(self, other):
        """
        Less-than comparison for sorting.

        Default compares core fields in definition order.
        Subclasses should override for custom sort order.
        """
        for field in self._core_fields:
            self_val = getattr(self, field)
            other_val = getattr(other, field)
            if self_val < other_val:
                return True
            elif self_val > other_val:
                return False
        return False
