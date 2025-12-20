"""
Feature processing utilities for genome browser.

This module handles feature extraction, filtering, grouping, and processing
for genomic annotations.
"""

from typing import Dict, List, Optional, Set, Tuple, Any
from collections import defaultdict
from ggene.database.annotations import UFeature


class FeatureProcessor:
    """Processes and organizes genomic features for display."""

    def __init__(self):
        """Initialize the feature processor."""
        self.feature_hierarchy = {
            'gene': 0,
            'transcript': 1,
            'exon': 2,
            'CDS': 2,
            'five_prime_utr': 3,
            'three_prime_utr': 3,
            'intron': 3,
            'variant': 4,
            'motif': 5,
            'other': 6
        }

    def get_all_features(self, window, state, genome_manager,
                         existing_features: List = None) -> List:
        """Aggregate features from multiple sources.

        Args:
            window: Current window object with features
            state: Current browser state
            genome_manager: GenomeManager instance
            existing_features: Existing features to append to

        Returns:
            List of all features in the window
        """
        if existing_features is None:
            features = []
        else:
            features = list(existing_features)  # Make a copy

        # Add motifs as features
        if window.motifs:
            for motif in window.motifs:
                # Convert motif to feature format
                motif_copy = dict(motif)
                if 'sequence' in motif_copy:
                    motif_copy.pop('sequence')
                motif_copy['feature_type'] = motif_copy.pop('type', 'motif')
                motif_copy['source'] = 'motif'

                # Create UFeature
                motif_feat = UFeature(chrom=state.chrom, **motif_copy)
                features.append(motif_feat)

        # Add annotations from genome manager
        if hasattr(genome_manager, 'annotations'):
            unified_features = genome_manager.get_all_annotations(
                str(state.chrom),
                state.position,
                state.position + state.window_size - 1,
                include_motifs=True
            )
            features.extend(unified_features)

        return features

    def filter_by_types(self, features: List, feature_types: Set[str]) -> List:
        """Filter features by type.

        Args:
            features: List of features
            feature_types: Set of feature types to include

        Returns:
            Filtered list of features
        """
        return [f for f in features if self._get_feature_type(f) in feature_types]

    def group_features_by_type(self, features: List) -> Dict[str, List]:
        """Group features by their type.

        Args:
            features: List of features

        Returns:
            Dictionary mapping feature types to lists of features
        """
        grouped = defaultdict(list)
        for feature in features:
            ftype = self._get_feature_type(feature)
            grouped[ftype].append(feature)
        return dict(grouped)

    def group_features_by_extent(self, features: List, window_start: int,
                                  window_size: int) -> List[Tuple[str, int, int, List]]:
        """Group features by their spatial extent for display.

        Args:
            features: List of features
            window_start: Start position of the window
            window_size: Size of the display window

        Returns:
            List of (feature_type, start, end, features) tuples
        """
        # Group by type first
        by_type = self.group_features_by_type(features)

        # Then group by overlapping extent
        grouped_extents = []

        for ftype, type_features in by_type.items():
            # Sort features by start position
            sorted_features = sorted(type_features, key=lambda f: f.start)

            # Group overlapping features
            groups = []
            current_group = []
            current_end = -1

            for feature in sorted_features:
                # Convert to window coordinates
                rel_start = max(0, feature.start - window_start)
                rel_end = min(window_size - 1, feature.end - window_start)

                if not current_group or feature.start <= current_end + 1:
                    # Add to current group
                    current_group.append(feature)
                    current_end = max(current_end, feature.end)
                else:
                    # Start new group
                    if current_group:
                        group_start = max(0, current_group[0].start - window_start)
                        group_end = min(window_size - 1, current_end - window_start)
                        groups.append((ftype, group_start, group_end, current_group))

                    current_group = [feature]
                    current_end = feature.end

            # Add final group
            if current_group:
                group_start = max(0, current_group[0].start - window_start)
                group_end = min(window_size - 1, current_end - window_start)
                groups.append((ftype, group_start, group_end, current_group))

            grouped_extents.extend(groups)

        return grouped_extents

    def get_feature_hierarchy_level(self, feature_type: str) -> int:
        """Get the hierarchy level for a feature type.

        Args:
            feature_type: Type of the feature

        Returns:
            Hierarchy level (lower is higher priority)
        """
        return self.feature_hierarchy.get(feature_type.lower(), 6)

    def filter_overlapping_features(self, features: List) -> List:
        """Filter overlapping features based on hierarchy.

        Args:
            features: List of features

        Returns:
            List of non-overlapping features
        """
        # Sort by hierarchy level and start position
        sorted_features = sorted(
            features,
            key=lambda f: (self.get_feature_hierarchy_level(self._get_feature_type(f)),
                           f.start)
        )

        # Keep non-overlapping features
        kept_features = []
        covered_ranges = []

        for feature in sorted_features:
            # Check if this feature overlaps with any kept feature
            overlaps = False
            for start, end in covered_ranges:
                if not (feature.end < start or feature.start > end):
                    overlaps = True
                    break

            if not overlaps:
                kept_features.append(feature)
                covered_ranges.append((feature.start, feature.end))

        return kept_features

    def extract_cds_features(self, features: List) -> List:
        """Extract CDS features from transcript features.

        Args:
            features: List of features

        Returns:
            List of CDS features
        """
        cds_features = []

        for feature in features:
            ftype = self._get_feature_type(feature)

            if ftype == 'CDS':
                cds_features.append(feature)
            elif ftype == 'transcript' and hasattr(feature, 'cds_start'):
                # Create CDS feature from transcript
                cds_feat = UFeature(
                    chrom=feature.chrom,
                    start=feature.cds_start,
                    end=feature.cds_end,
                    feature_type='CDS',
                    source=feature.source,
                    attributes={'parent': feature.name}
                )
                cds_features.append(cds_feat)

        return cds_features

    def find_features_at_position(self, features: List, position: int) -> List:
        """Find all features that overlap a specific position.

        Args:
            features: List of features
            position: Genomic position

        Returns:
            List of overlapping features
        """
        return [f for f in features if f.start <= position <= f.end]

    def get_feature_info(self, feature) -> Dict[str, Any]:
        """Extract displayable information from a feature.

        Args:
            feature: Feature object

        Returns:
            Dictionary of feature information
        """
        info = {}

        # Get basic info
        info['type'] = self._get_feature_type(feature)
        info['start'] = feature.start
        info['end'] = feature.end

        # Get name
        if hasattr(feature, 'name'):
            info['name'] = feature.name
        elif hasattr(feature, 'attributes'):
            attrs = feature.attributes if isinstance(feature.attributes, dict) else {}
            info['name'] = attrs.get('gene_name', attrs.get('transcript_name', ''))

        # Get additional attributes
        if hasattr(feature, 'attributes') and isinstance(feature.attributes, dict):
            info.update(feature.attributes)

        return info

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
        else:
            return 'unknown'

    def cache_current_gene(self, features: List) -> Optional[Dict]:
        """Find and cache the current gene from features.

        Args:
            features: List of features

        Returns:
            Current gene information or None
        """
        genes = [f for f in features if self._get_feature_type(f) == 'gene']

        if genes:
            # Return the first gene (could be enhanced to pick the "best" gene)
            gene = genes[0]
            return self.get_feature_info(gene)

        return None

    def get_transcripts_for_gene(self, features: List, gene_name: str) -> List:
        """Get all transcripts for a specific gene.

        Args:
            features: List of features
            gene_name: Name of the gene

        Returns:
            List of transcript features
        """
        transcripts = []

        for feature in features:
            ftype = self._get_feature_type(feature)
            if ftype == 'transcript':
                # Check if transcript belongs to this gene
                info = self.get_feature_info(feature)
                if info.get('gene_name') == gene_name or info.get('parent') == gene_name:
                    transcripts.append(feature)

        return transcripts