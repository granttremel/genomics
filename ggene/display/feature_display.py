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
        
        detail = self._render_feature_detail(state, None, features, ["motif"])
        lines.extend(detail)
        
        return lines

    # would rather order by feature hierarchy.. 

    def _group_features_by_extent(self, features: List, window_start: int) -> List[Dict]:
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
        return Colors.get_feature_color(feature_type)

    def _get_feature_type(self, feature) -> str:
        if hasattr(feature, 'feature_type'):
            return feature.feature_type
        elif hasattr(feature, 'type'):
            return feature.type
        elif isinstance(feature, dict):
            return feature.get('feature_type', feature.get('type', 'unknown'))
        return 'unknown'
    
    def _render_feature_detail(self, state, window, features, selected_features):
        
        detail_lines = []
        
        
        if "motif" in selected_features:
            motifs = self.feature_processor.filter_by_types(features, ["motif"])
            mlines = self._render_motif_detail(state, motifs)
            detail_lines.extend(mlines)
        return detail_lines
    
    def _render_anno_detail(self, ):
        pass
    
    def _render_variant_detail(self, ):
        
        pass    
    def _render_motif_detail(self, state, motifs):
        lines = []
        header = self.rstate._get_margin("Motifs:")
        lines.append(header)
        
        # Group motifs by type
        motif_groups = {}
        for motif in motifs:
            motif_info = motif.get('info', {}) if isinstance(motif, dict) else {}
            motif_name = motif_info.get('name', motif.get('name', 'unknown'))
            strand = motif.get('strand', '+')
            
            if motif_name not in motif_groups:
                motif_groups[motif_name] = []
            motif_groups[motif_name].append(motif)
        
        # Display each motif type
        for motif_type, motif_list in sorted(motif_groups.items()):
            # Choose color for motif type
            marg = self.rstate._get_margin()
            lines.append(f"{Colors.MOTIF}{motif_type}{Colors.RESET}")
            
            for motif in motif_list[:5]:  # Show max 5 instances per type
                start = motif.get('start')
                end = motif.get('end')
                strand = motif.get('strand', '+')
                score = motif.get('score', motif.get('info', {}).get('score', 0))
                seq = motif.get('sequence', motif.get('info', {}).get('sequence', ''))
                
                # Determine which sequence (ref or personal)
                in_window_start = max(0, start - state.position)
                in_window_end = min(state.window_size, end - state.position)
                
                if in_window_start < in_window_end:
                    strand_symbol = '→' if strand == '+' else '←'
                    lines.append(f"{marg}{start:,} {strand_symbol} {seq} (score: {score:.2f})")
            
            if len(motif_list) > 5:
                lines.append(f"    ... and {len(motif_list) - 5} more")
        return lines