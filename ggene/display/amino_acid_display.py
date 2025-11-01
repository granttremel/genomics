"""
Amino acid display module for genome browser.

This module handles the rendering of amino acid translations and changes
between reference and personal sequences.
"""

from typing import List, Tuple, Optional, Dict
from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA
from ggene.display.colors import Colors
from ggene.processing.sequence_processor import SequenceProcessor
import logging

logger = logging.getLogger(__name__)


class AminoAcidDisplay:
    """Renders amino acid translations and changes."""

    def __init__(self):
        """Initialize the amino acid display renderer."""
        self.seq_processor = SequenceProcessor()
        self.three_letter_codes = {
            'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
            'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
            'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
            'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
            '*': 'Ter', 'X': 'Unk'
        }

    def render_amino_acid_changes(self, ref_seq: str, personal_seq: str,
                                   features: List = None, state=None) -> List[str]:
        """Display amino acid changes between reference and personal sequences.

        Only shows changes in CDS regions where both sequences can be translated.

        Args:
            ref_seq: Reference DNA/RNA sequence
            personal_seq: Personal DNA/RNA sequence
            features: List of genomic features
            state: Browser state object

        Returns:
            List of display lines
        """
        if not ref_seq or not personal_seq:
            return []

        lines = []

        # Find CDS features
        cds_features = self._extract_cds_features(features)

        if not cds_features:
            # No CDS features, try to translate the entire sequence in frame 0
            lines.extend(self._translate_full_sequence(ref_seq, personal_seq))
        else:
            # Translate each CDS region
            for cds in cds_features:
                cds_lines = self._translate_cds_region(
                    ref_seq, personal_seq, cds, state
                )
                lines.extend(cds_lines)

        return lines

    def _extract_cds_features(self, features: List) -> List:
        """Extract CDS features from feature list.

        Args:
            features: List of features

        Returns:
            List of CDS features
        """
        if not features:
            return []

        cds_features = []
        for feature in features:
            ftype = self._get_feature_type(feature)
            if ftype == 'CDS':
                cds_features.append(feature)

        return cds_features

    def _translate_full_sequence(self, ref_seq: str, personal_seq: str) -> List[str]:
        """Translate full sequences in all three frames.

        Args:
            ref_seq: Reference sequence
            personal_seq: Personal sequence

        Returns:
            List of display lines
        """
        lines = []

        # Remove gaps for translation
        ref_no_gaps = ref_seq.replace('-', '')
        pers_no_gaps = personal_seq.replace('-', '')

        # Translate in all three frames
        for frame in range(3):
            if len(ref_no_gaps) > frame + 2 and len(pers_no_gaps) > frame + 2:
                # Translate sequences
                ref_aa = self._translate_sequence(ref_no_gaps[frame:], frame)
                pers_aa = self._translate_sequence(pers_no_gaps[frame:], frame)

                # Find changes
                changes = self._find_amino_acid_changes(ref_aa, pers_aa)

                if changes:
                    lines.append(f"Frame {frame}: {self._format_changes(changes)}")

        return lines

    def _translate_cds_region(self, ref_seq: str, personal_seq: str,
                               cds_feature, state) -> List[str]:
        """Translate a specific CDS region.

        Args:
            ref_seq: Reference sequence
            personal_seq: Personal sequence
            cds_feature: CDS feature object
            state: Browser state

        Returns:
            List of display lines
        """
        lines = []

        # Get CDS boundaries relative to window
        window_start = state.position if state else 0
        cds_start = max(0, cds_feature.start - window_start)
        cds_end = min(len(ref_seq), cds_feature.end - window_start + 1)

        if cds_start >= cds_end:
            return lines

        # Extract CDS sequences
        ref_cds = ref_seq[cds_start:cds_end].replace('-', '')
        pers_cds = personal_seq[cds_start:cds_end].replace('-', '')

        # Get frame from feature if available
        frame = getattr(cds_feature, 'frame', 0)

        # Translate
        ref_aa = self._translate_sequence(ref_cds, frame)
        pers_aa = self._translate_sequence(pers_cds, frame)

        # Find changes
        changes = self._find_amino_acid_changes(ref_aa, pers_aa)

        if changes:
            # Get CDS name if available
            cds_name = getattr(cds_feature, 'name', 'CDS')
            lines.append(f"{cds_name}: {self._format_changes(changes)}")

        return lines

    def _translate_sequence(self, dna_seq: str, frame: int = 0) -> str:
        """Translate DNA sequence to amino acids.

        Args:
            dna_seq: DNA sequence (gaps should be removed)
            frame: Reading frame (0, 1, or 2)

        Returns:
            Amino acid sequence
        """
        # Use sequence processor for translation
        if frame > 0 and frame < len(dna_seq):
            dna_seq = dna_seq[frame:]

        return self.seq_processor.translate_dna(dna_seq)

    def _find_amino_acid_changes(self, ref_aa: str, pers_aa: str) -> List[Tuple[int, str, str]]:
        """Find amino acid changes between two sequences.

        Args:
            ref_aa: Reference amino acid sequence
            pers_aa: Personal amino acid sequence

        Returns:
            List of (position, ref_aa, pers_aa) tuples
        """
        changes = []
        min_len = min(len(ref_aa), len(pers_aa))

        for i in range(min_len):
            if ref_aa[i] != pers_aa[i]:
                # Skip unknown amino acids
                if ref_aa[i] != 'X' and pers_aa[i] != 'X':
                    changes.append((i + 1, ref_aa[i], pers_aa[i]))

        return changes

    def _format_changes(self, changes: List[Tuple[int, str, str]],
                        max_changes: int = 10) -> str:
        """Format amino acid changes for display.

        Args:
            changes: List of (position, ref_aa, pers_aa) tuples
            max_changes: Maximum number of changes to display

        Returns:
            Formatted string of changes
        """
        if not changes:
            return "No changes"

        formatted = []
        for pos, ref_aa, pers_aa in changes[:max_changes]:
            # Get three-letter codes
            ref_3 = self.three_letter_codes.get(ref_aa, ref_aa)
            pers_3 = self.three_letter_codes.get(pers_aa, pers_aa)

            # Color based on change type
            if ref_aa == '*' or pers_aa == '*':
                # Stop codon change
                color = Colors.DELETION
            elif self._is_conservative_change(ref_aa, pers_aa):
                # Conservative change
                color = Colors.CDS
            else:
                # Non-conservative change
                color = Colors.SNP

            formatted.append(f"{color}{ref_3}{pos}{pers_3}{Colors.RESET}")

        result = " ".join(formatted)
        if len(changes) > max_changes:
            result += f" ... (+{len(changes) - max_changes} more)"

        return result

    def _is_conservative_change(self, aa1: str, aa2: str) -> bool:
        """Check if amino acid change is conservative.

        Args:
            aa1: First amino acid
            aa2: Second amino acid

        Returns:
            True if change is conservative
        """
        # Define amino acid groups for conservative changes
        groups = [
            {'A', 'V', 'L', 'I', 'M'},  # Hydrophobic aliphatic
            {'F', 'W', 'Y'},             # Aromatic
            {'S', 'T', 'C'},             # Polar uncharged small
            {'N', 'Q'},                  # Polar uncharged large
            {'D', 'E'},                  # Acidic
            {'K', 'R', 'H'},             # Basic
            {'G'},                       # Glycine
            {'P'}                        # Proline
        ]

        for group in groups:
            if aa1 in group and aa2 in group:
                return True

        return False

    def render_amino_acid_tracks(self, ref_seq: str, personal_seq: str,
                                  features: List = None) -> List[str]:
        """Render amino acid sequences as aligned tracks.

        Args:
            ref_seq: Reference DNA sequence
            personal_seq: Personal DNA sequence
            features: List of features (for CDS identification)

        Returns:
            List of display lines showing AA tracks
        """
        lines = []

        # Remove gaps
        ref_no_gaps = ref_seq.replace('-', '')
        pers_no_gaps = personal_seq.replace('-', '')

        # Translate in three frames
        for frame in range(3):
            if len(ref_no_gaps) > frame + 2:
                ref_aa = self._translate_sequence(ref_no_gaps[frame:], 0)
                pers_aa = self._translate_sequence(pers_no_gaps[frame:], 0)

                # Format as tracks with position markers
                ref_track = self._format_aa_track(ref_aa, "Ref", frame)
                pers_track = self._format_aa_track(pers_aa, "Alt", frame)

                if ref_track or pers_track:
                    lines.append(f"Frame {frame}:")
                    if ref_track:
                        lines.append(ref_track)
                    if pers_track:
                        lines.append(pers_track)

        return lines

    def _format_aa_track(self, aa_seq: str, label: str, frame: int) -> str:
        """Format amino acid sequence as a display track.

        Args:
            aa_seq: Amino acid sequence
            label: Track label
            frame: Reading frame

        Returns:
            Formatted track string
        """
        if not aa_seq:
            return ""

        # Use spacing to align with nucleotide display (3 nt per AA)
        spaced_aa = "  ".join(aa_seq)  # Two spaces between each AA

        return f"  {label}F{frame}: {spaced_aa}"

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