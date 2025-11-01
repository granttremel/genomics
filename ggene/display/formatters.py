"""
Label formatting utilities for genome browser display.

This module provides functions for formatting feature labels with various
compression strategies for dense genomic regions.
"""


class LabelFormatter:
    """Formats labels for genomic features with intelligent compression."""

    def __init__(self, window_size: int = 240):
        """Initialize the label formatter.

        Args:
            window_size: Size of the display window
        """
        self.window_size = window_size

    def format_compressed_label(self, ftype: str, features, window_start: int) -> str:
        """Format a compressed label for multiple features with advanced compression.

        Args:
            ftype: Feature type
            features: List of features with same extent
            window_start: Start position of the current window

        Returns:
            Formatted compressed label string
        """
        if len(features) <= 5:
            return self.format_single_feature_label(features[0])

        # Check if all features span the full window - compress heavily if so
        window_end = window_start + self.window_size - 1
        all_span_window = all(
            f.start <= window_start and f.end >= window_end
            for f in features
        )

        if all_span_window and len(features) >= 5:
            # Heavy compression for fully spanning features
            if ftype == 'transcript':
                # Extract transcript suffixes and compress
                suffixes = []
                for f in features:
                    info = f.get('info', {})
                    name = info.get('transcript_name', info.get('transcript_id', '?'))
                    if '-' in name:
                        suffixes.append(name.split('-')[-1])
                    else:
                        suffixes.append(name)

                # Try to find patterns like 001, 002, 003...
                numeric_suffixes = [s for s in suffixes if s.isdigit()]
                if len(numeric_suffixes) >= 3:
                    sorted_nums = sorted(numeric_suffixes, key=int)
                    if len(sorted_nums) >= 5:
                        return f"{sorted_nums[0]}-{sorted_nums[-1]}({len(features)})"

                # Fallback to first few + count
                unique_suffixes = sorted(set(suffixes))[:3]
                return f"{'|'.join(unique_suffixes)}+{len(features)-len(unique_suffixes)}"

            elif ftype == 'exon':
                nums = [str(f.get('info', {}).get('exon_number', '?')) for f in features]
                numeric_nums = [n for n in nums if n.isdigit()]
                if len(numeric_nums) >= 3:
                    sorted_nums = sorted(numeric_nums, key=int)
                    return f"ex{sorted_nums[0]}-{sorted_nums[-1]}"
                return f"exons({len(features)})"

            else:
                return f"{ftype}({len(features)})"

        # Normal compression for smaller sets or partial spans
        if ftype == 'gene':
            names = [f.get('info', {}).get('gene_name', '?') for f in features]
            unique_names = sorted(set(names))
            if len(unique_names) > 3:
                return f"{unique_names[0]}/+{len(unique_names)-1}"
            return '/'.join(unique_names)[:25]

        elif ftype == 'transcript':
            names = []
            for f in features:
                info = f.get('info', {})
                name = info.get('transcript_name', info.get('transcript_id', '?'))
                if '-' in name:
                    names.append(name.split('-')[-1])
                else:
                    names.append(name)
            unique_names = sorted(set(names))
            if len(unique_names) > 4:
                return f"{unique_names[0]}/+{len(unique_names)-1}"
            return '/'.join(unique_names)[:25]

        elif ftype == 'exon':
            nums = [str(f.get('info', {}).get('exon_number', '?')) for f in features]
            unique_nums = sorted(set(nums), key=lambda x: int(x) if x.isdigit() else 0)
            if len(unique_nums) > 4:
                return f"ex{unique_nums[0]}-{unique_nums[-1]}"
            return f"ex{'/'.join(unique_nums)}"

        elif ftype == 'CDS':
            return f"CDS({len(features)})" if len(features) > 1 else "CDS"
        else:
            return f"{ftype}({len(features)})" if len(features) > 1 else ftype[:10]

    def format_single_feature_label(self, feature) -> str:
        """Format a label for a single feature.

        Args:
            feature: Feature object or dictionary

        Returns:
            Formatted label string
        """
        ftype = getattr(feature, 'feature_type', feature.get('feature_type', 'unknown'))

        # Handle both attribute-style and dictionary-style access
        if hasattr(feature, 'attributes'):
            info = feature.attributes
        else:
            info = feature.get('attributes', feature.get('info', {}))

        if ftype == 'gene':
            if hasattr(feature, 'name') and feature.name:
                return feature.name
            return info.get('gene_name', 'gene')
        elif ftype == 'transcript':
            name = info.get('transcript_name', info.get('transcript_id', 'transcript'))
            # Shorten long transcript names
            if len(name) > 20:
                if '-' in name:
                    return name.split('-')[-1]
                return name[:20]
            return name
        elif ftype == 'exon':
            exon_num = info.get('exon_number', '?')
            return f"exon {exon_num}"
        elif ftype == 'CDS':
            return "CDS"
        elif ftype == 'variant':
            ref = info.get('ref', '?')
            alt_list = info.get('alt', ['?'])
            alt = alt_list[0] if isinstance(alt_list, list) else alt_list

            # Apply sequence rendering if needed (RNA mode)
            if hasattr(self, 'render_seq_func'):
                ref = self.render_seq_func(ref)
                alt = self.render_seq_func(alt)

            if len(ref) == len(alt):
                return f"{ref}→{alt}"
            else:
                return f"VAR({len(ref)}→{len(alt)})"
        elif ftype == 'five_prime_utr':
            return "5'UTR"
        elif ftype == 'three_prime_utr':
            return "3'UTR"
        elif ftype == 'intron':
            return "intron"
        else:
            return ftype[:10]

    def set_render_seq_func(self, render_func):
        """Set a custom sequence rendering function for variant labels.

        Args:
            render_func: Function to transform sequences (e.g., DNA to RNA)
        """
        self.render_seq_func = render_func