# Genome Browser Refactoring Guide

## Overview

The genome browser has been refactored from a monolithic 2,115-line file into a modular architecture with clear separation of concerns. This guide explains the new structure and how to use it.

## New Module Structure

```
ggene/
├── browser/           # Browser orchestration (future)
├── display/           # All display/rendering logic
│   ├── colors.py      # ANSI color codes and color management
│   ├── formatters.py  # Label formatting and compression
│   ├── renderer.py    # Main display coordinator
│   ├── sequence_display.py    # DNA/RNA sequence rendering
│   ├── feature_display.py     # Feature bar rendering
│   └── amino_acid_display.py  # Amino acid translation display
├── processing/        # Data processing logic
│   ├── coordinate_mapper.py   # Coordinate transformations
│   ├── sequence_processor.py  # Sequence analysis and translation
│   └── feature_processor.py   # Feature extraction and filtering
└── visualization/     # Advanced visualization (future)
```

## Module Descriptions

### Display Modules

#### `display/colors.py`
- **Purpose**: Centralized color management for terminal display
- **Key Classes**: `Colors`
- **Extracted From**: Lines 28-65 of genome_browser.py
- **Usage**:
```python
from ggene.display.colors import Colors
print(f"{Colors.SNP}A→T{Colors.RESET}")
```

#### `display/formatters.py`
- **Purpose**: Format feature labels with intelligent compression
- **Key Classes**: `LabelFormatter`
- **Extracted From**: Lines 1302-1440 of genome_browser.py
- **Methods**:
  - `format_compressed_label()`: Compress multiple feature labels
  - `format_single_feature_label()`: Format single feature

#### `display/sequence_display.py`
- **Purpose**: Render DNA/RNA sequences with variant highlighting
- **Key Classes**: `SequenceDisplay`
- **Extracted From**: Lines 762-1045 of genome_browser.py
- **Features**:
  - Variant detection and highlighting
  - Codon annotation (start/stop)
  - Motif underlining
  - Position ruler generation

#### `display/feature_display.py`
- **Purpose**: Render genomic features as ASCII art bars
- **Key Classes**: `FeatureDisplay`
- **Extracted From**: Lines 1047-1300 of genome_browser.py
- **Features**:
  - Feature grouping by overlap
  - Row assignment for non-overlapping display
  - Customizable bar styles per feature type

#### `display/amino_acid_display.py`
- **Purpose**: Display amino acid translations and changes
- **Key Classes**: `AminoAcidDisplay`
- **Extracted From**: Lines 1442-1629 of genome_browser.py
- **Features**:
  - CDS region translation
  - Amino acid change detection
  - Conservative vs non-conservative change coloring

#### `display/renderer.py`
- **Purpose**: Coordinate all display components
- **Key Classes**: `DisplayRenderer`
- **New Addition**: Orchestrates all display modules
- **Features**:
  - Unified rendering pipeline
  - Section interleaving
  - Header/footer generation

### Processing Modules

#### `processing/coordinate_mapper.py`
- **Purpose**: Handle coordinate system transformations
- **Key Classes**: `CoordinateMapper`, `SequenceRenderer`
- **Extracted From**: Lines 217-228 of genome_browser.py
- **Features**:
  - Window to genomic coordinate mapping
  - Strand-aware transformations
  - Sequence rendering (reverse complement, RNA)

#### `processing/sequence_processor.py`
- **Purpose**: Sequence analysis and manipulation
- **Key Classes**: `SequenceProcessor`
- **Extracted From**: Lines 1605-1629 and various locations
- **Features**:
  - Variant detection
  - Codon finding
  - DNA/RNA translation
  - ORF detection
  - Motif searching
  - GC content calculation

#### `processing/feature_processor.py`
- **Purpose**: Feature extraction and organization
- **Key Classes**: `FeatureProcessor`
- **Extracted From**: Lines 599-617 and various locations
- **Features**:
  - Feature aggregation from multiple sources
  - Feature filtering by type
  - Grouping by overlap
  - Hierarchy management

## Integration with Existing Code

### Gradual Migration Approach

The refactored modules are designed to work alongside the existing code. You can migrate gradually:

1. **Phase 1**: Use new modules for new features
```python
from ggene.display.colors import Colors
from ggene.display.formatters import LabelFormatter

# Use in existing code
formatter = LabelFormatter(window_size)
label = formatter.format_single_feature_label(feature)
```

2. **Phase 2**: Replace individual methods
```python
# Instead of self._display_sequences()
from ggene.display.sequence_display import SequenceDisplay
seq_display = SequenceDisplay(self.state.window_size)
lines, n = seq_display.render_sequences(ref_seq, alt_seq, features, self.state)
```

3. **Phase 3**: Use the unified renderer
```python
from ggene.display.renderer import DisplayRenderer

class InteractiveGenomeBrowser:
    def __init__(self, genome_manager):
        self.renderer = DisplayRenderer(genome_manager)

    def display(self):
        all_lines = self.renderer.render_full_view(window, self.state)
        output = self.renderer.render_interleaved(all_lines)
        print(output)
```

## Adding New Features

### Example: Adding Regex Pattern Highlighting

```python
# In processing/sequence_processor.py
def find_regex_patterns(self, seq: str, pattern: str) -> List[int]:
    """Find all matches of a regex pattern in sequence."""
    import re
    matches = []
    for match in re.finditer(pattern, seq, re.IGNORECASE):
        matches.append(match.start())
    return matches

# In display/sequence_display.py
def highlight_patterns(self, seq: str, pattern_positions: List[int]) -> str:
    """Highlight regex pattern matches in sequence."""
    # Add highlighting logic
    pass
```

### Example: Adding New Data Visualizations

```python
# Create visualization/data_plots.py
class DataPlotter:
    def render_gc_content(self, seq: str, window_size: int) -> List[str]:
        """Render GC content as ASCII plot."""
        gc_values = self.calculate_gc_content(seq, window_size)
        return self.ascii_plot(gc_values)

    def render_correlation(self, seq1: str, seq2: str) -> List[str]:
        """Render sequence correlation plot."""
        # Implementation here
        pass
```

## Benefits of the New Structure

1. **Testability**: Each module can be unit tested independently
2. **Maintainability**: Changes to display don't affect processing logic
3. **Reusability**: Components can be used in other tools (GUI version, etc.)
4. **Extensibility**: Easy to add new display modes or data sources
5. **Performance**: Can optimize individual components
6. **Collaboration**: Multiple developers can work on different modules

## Migration Checklist

- [x] Extract Colors class
- [x] Extract formatters
- [x] Extract coordinate mapping
- [x] Extract sequence processing
- [x] Extract feature processing
- [x] Extract sequence display
- [x] Extract feature display
- [x] Extract amino acid display
- [x] Create unified renderer
- [ ] Update main browser to use renderer
- [ ] Add unit tests for each module
- [ ] Add integration with draw.py utilities
- [ ] Create navigation module
- [ ] Create state management module
- [ ] Add regex pattern highlighting
- [ ] Add additional data visualizations

## Example: Using the New Modules Together

```python
from ggene.display.renderer import DisplayRenderer
from ggene.processing.feature_processor import FeatureProcessor
from ggene.genome_iterator_v2 import UnifiedGenomeIterator

# Initialize components
genome_manager = GenomeManager()
renderer = DisplayRenderer(genome_manager)
feature_processor = FeatureProcessor()

# Get window data
iterator = UnifiedGenomeIterator(genome_manager, chrom=1, position=1000000)
window = iterator.get_window_at(position)

# Process features
features = feature_processor.get_all_features(window, state, genome_manager)

# Render display
display_sections = renderer.render_full_view(window, state)
output = renderer.render_interleaved(display_sections)
print(output)
```

## Testing

Each module should have corresponding tests:

```python
# tests/test_display/test_formatters.py
import pytest
from ggene.display.formatters import LabelFormatter

def test_format_single_gene():
    formatter = LabelFormatter()
    feature = {'feature_type': 'gene', 'name': 'BRCA1'}
    label = formatter.format_single_feature_label(feature)
    assert label == 'BRCA1'

def test_compress_multiple_transcripts():
    formatter = LabelFormatter(window_size=100)
    features = [
        {'feature_type': 'transcript', 'info': {'transcript_name': 'NM_001'}},
        {'feature_type': 'transcript', 'info': {'transcript_name': 'NM_002'}},
        {'feature_type': 'transcript', 'info': {'transcript_name': 'NM_003'}},
    ]
    label = formatter.format_compressed_label('transcript', features, 0)
    assert '001' in label or 'NM' in label
```

## Future Enhancements

1. **Browser Module**: Create dedicated browser orchestration module
2. **Visualization Module**: Integrate draw.py utilities
3. **Plugin System**: Allow external modules to add display components
4. **Configuration System**: Centralized configuration management
5. **Performance Monitoring**: Add profiling hooks to identify bottlenecks
6. **Caching Layer**: Add intelligent caching for expensive computations

## Notes for Contributors

- Keep modules focused on a single responsibility
- Use type hints for all public methods
- Add docstrings explaining purpose and usage
- Write unit tests for new functionality
- Follow existing naming conventions
- Update this guide when adding new modules