"""
Feature display module for genome browser.

This module handles the rendering of genomic features as visual bars and labels
in the terminal interface.
"""

from typing import List, Dict, Tuple, Optional, Any
from collections import defaultdict
from ggene.display.colors import Colors
from ggene.display.formatters import LabelFormatter
from ggene.processing.coordinate_mapper import CoordinateMapper
from ggene.processing.feature_processor import FeatureProcessor
import logging

logger = logging.getLogger(__name__)


class FeatureDisplay:
    """Renders genomic features as ASCII art bars with labels."""

    def __init__(self, render_state, window_size: int = 240):
        """Initialize the feature display renderer.

        Args:
            window_size: Size of the display window
        """
        self.rstate = render_state
        self.window_size = window_size
        self.coord_mapper = CoordinateMapper(window_size)
        self.label_formatter = LabelFormatter(window_size)
        self.feature_processor = FeatureProcessor()
        self.max_display_rows = 20  # Maximum rows for feature display

    def render_features(self, features: List, variant_features: List = None,
                        state=None) -> List[str]:
        """Render features as visual bars with labels.

        Args:
            features: List of genomic features
            variant_features: List of variant features (optional)
            state: Browser state object

        Returns:
            List of display lines
        """
        if not features:
            return []

        # Update settings based on state
        if state:
            self.window_size = state.window_size
            self.coord_mapper.update_settings(
                window_size=state.window_size,
                show_reverse_strand=state.show_reverse_strand
            )
            self.label_formatter.window_size = state.window_size

        # Filter features by selected types if state provides them
        if state and hasattr(state, 'selected_feature_types'):
            features = self.feature_processor.filter_by_types(
                features, state.selected_feature_types
            )

        # Group features by extent
        window_start = state.position if state else 0
        grouped_features = self._group_features_by_extent(features, window_start)

        # Assign features to display rows
        row_assignments = self._assign_features_to_rows(grouped_features)

        # Generate display lines
        lines = []
        for row_idx in range(min(len(row_assignments), self.max_display_rows)):
            line = self._render_feature_row(row_assignments[row_idx], window_start, state)
            if line:
                lines.append(line)

        return lines

    # would rather order by feature hierarchy.. 

    def _group_features_by_extent(self, features: List, window_start: int) -> List[Dict]:
        """Group overlapping features by their extent.

        Args:
            features: List of features
            window_start: Start position of window

        Returns:
            List of feature groups with metadata
        """
        # Group features by type and overlapping extent
        type_groups = defaultdict(list)

        for feature in features:
            ftype = self._get_feature_type(feature)
            type_groups[ftype].append(feature)

        # Process each type group
        grouped = []
        for ftype, type_features in type_groups.items():
            # Sort by start position
            sorted_features = sorted(type_features, key=lambda f: f.start)

            # Group overlapping features
            current_group = []
            group_start = -1
            group_end = -1

            for feature in sorted_features:
                feat_start = feature.start
                feat_end = feature.end

                if not current_group or feat_start <= group_end + 1:
                    # Add to current group
                    current_group.append(feature)
                    if group_start == -1:
                        group_start = feat_start
                    group_end = max(group_end, feat_end)
                else:
                    # Save current group and start new one
                    if current_group:
                        grouped.append({
                            'type': ftype,
                            'start': group_start,
                            'end': group_end,
                            'features': current_group
                        })
                    current_group = [feature]
                    group_start = feat_start
                    group_end = feat_end

            # Add final group
            if current_group:
                grouped.append({
                    'type': ftype,
                    'start': group_start,
                    'end': group_end,
                    'features': current_group
                })

        return grouped

    def _assign_features_to_rows(self, feature_groups: List[Dict]) -> List[List[Dict]]:
        """Assign feature groups to non-overlapping display rows.

        Args:
            feature_groups: List of feature groups

        Returns:
            List of rows, each containing non-overlapping feature groups
        """
        # Sort by priority (hierarchy) and start position
        sorted_groups = sorted(
            feature_groups,
            key=lambda g: (
                self.feature_processor.get_feature_hierarchy_level(g['type']),
                g['start']
            )
        )

        # Assign to rows
        rows = []
        for group in sorted_groups:
            # Find first row where this group fits
            placed = False
            for row in rows:
                # Check if group overlaps with any feature in this row
                overlaps = False
                for existing in row:
                    if not (group['end'] < existing['start'] or
                            group['start'] > existing['end']):
                        overlaps = True
                        break

                if not overlaps:
                    row.append(group)
                    placed = True
                    break

            # If no existing row works, create new row
            if not placed:
                rows.append([group])

        return rows

    def _render_feature_row(self, feature_groups: List[Dict], window_start: int,
                             state) -> Optional[str]:
        """Render a single row of features.

        Args:
            feature_groups: List of feature groups in this row
            window_start: Start position of window
            state: Browser state

        Returns:
            Rendered line string or None
        """
        if not feature_groups:
            return None
        
        # Create display buffer
        display_buffer = [' '] * self.window_size

        # Render each feature group
        for group in feature_groups:
            # Calculate display positions
            rel_start = max(0, group['start'] - window_start)
            rel_end = min(self.window_size - 1, group['end'] - window_start)

            # Apply strand-aware transformation
            if state and state.show_reverse_strand:
                disp_start, disp_end = self.coord_mapper.render_pos(rel_start, rel_end)
            else:
                disp_start, disp_end = rel_start, rel_end

            # Ensure valid range
            disp_start = max(0, min(disp_start, self.window_size - 1))
            disp_end = max(0, min(disp_end, self.window_size - 1))

            if disp_start > disp_end:
                disp_start, disp_end = disp_end, disp_start

            # Get color for this feature type
            color = self._get_feature_color(group['type'])

            # Generate label
            label = self._generate_feature_label(group, disp_end - disp_start + 1)

            # Render feature bar
            self._render_feature_bar(display_buffer, disp_start, disp_end,
                                      label, color, group['type'])

        # Convert buffer to string
        return f"{self.rstate._get_margin()}{''.join(display_buffer)}"

    def _render_feature_bar(self, buffer: List[str], start: int, end: int,
                             label: str, color: str, ftype: str):
        """Render a feature bar into the display buffer.

        Args:
            buffer: Display buffer (modified in place)
            start: Start position in buffer
            end: End position in buffer
            label: Feature label
            color: ANSI color code
            ftype: Feature type
        """
        # Calculate available space
        width = end - start + 1

        if width < 1:
            return

        left_char, fill_char, right_char = list(self.rstate.feature_disp)

        if width > 5:
            # Render with label
            buffer[start] = f"{color}{left_char}{Colors.RESET}"

            # Center the label
            label_start = start + (width - len(label)) // 2
            for i, char in enumerate(label):
                if label_start + i < end:
                    buffer[label_start + i] = f"{color}{Colors.BOLD}{char}{Colors.RESET}"

            # Fill remaining space
            for i in range(start + 1, label_start):
                if buffer[i] == ' ':
                    buffer[i] = f"{color}{fill_char}{Colors.RESET}"
            for i in range(label_start + len(label), end):
                if buffer[i] == ' ':
                    buffer[i] = f"{color}{fill_char}{Colors.RESET}"

            buffer[end] = f"{color}{right_char}{Colors.RESET}"

    def _generate_feature_label(self, feature_group: Dict, available_width: int) -> str:
        """Generate a label for a feature group.

        Args:
            feature_group: Feature group dictionary
            available_width: Available width for label

        Returns:
            Label string
        """
        features = feature_group['features']
        ftype = feature_group['type']

        # Use label formatter for compression
        if len(features) > 1:
            label = self.label_formatter.format_compressed_label(
                ftype, features, 0  # window_start not needed for relative labels
            )
        else:
            label = self.label_formatter.format_single_feature_label(features[0])

        # Truncate if too long
        if len(label) > available_width - 2:
            if available_width > 5:
                label = label[:available_width-5] + "..."
            else:
                label = ""

        return label

    def _get_feature_color(self, feature_type: str) -> str:
        """Get ANSI color code for a feature type.

        Args:
            feature_type: Type of feature

        Returns:
            ANSI color code string
        """
        return Colors.get_feature_color(feature_type)

    def _get_feature_type(self, feature) -> str:
        """Get the type of a feature.

        Args:
            feature: Feature object

        Returns:
            Feature type string
        """
        if hasattr(feature, 'feature_type'):
            return feature.feature_type
        elif hasattr(feature, 'type'):
            return feature.type
        elif isinstance(feature, dict):
            return feature.get('feature_type', feature.get('type', 'unknown'))
        return 'unknown'