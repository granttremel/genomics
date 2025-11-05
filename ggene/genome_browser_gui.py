#!/usr/bin/env python3
"""
PyQt6-based interactive genome browser for visualizing personal variants and features.
"""
import sys
import os
import subprocess
import logging
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
from dataclasses import dataclass
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QTreeWidget, QTreeWidgetItem, QTextEdit, QLineEdit,
    QLabel, QScrollArea, QToolBar, QStatusBar, QDockWidget,
    QPushButton, QComboBox, QSpinBox, QCheckBox, QGroupBox,
    QMessageBox, QInputDialog, QGraphicsView, QGraphicsScene
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, QRectF, QPointF
from PyQt6.QtGui import (
    QFont, QFontMetrics, QColor, QPalette, QTextCharFormat,
    QTextCursor, QKeySequence, QShortcut, QAction, QPen,
    QBrush, QPainter, QPolygonF
)
from PyQt6.QtWidgets import QGraphicsTextItem, QGraphicsRectItem, QGraphicsPolygonItem

from PyQt6 import QtWidgets
from PyQt6 import QtCore
# QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) 

if TYPE_CHECKING:
    from .genomemanager import GenomeManager

from ggene import CODON_TABLE_DNA, to_rna, complement, reverse_complement
from .genome_browser_v2 import BrowserState, Colors as TermColors
from .genome_iterator import GenomeIterator

logger = logging.getLogger(__name__)


class GenomeColors:
    """Color scheme for the genome browser."""
    # Variant colors
    SNP = QColor(220, 50, 47)  # Red
    INSERTION = QColor(133, 153, 0)  # Green
    DELETION = QColor(181, 137, 0)  # Yellow/Orange
    
    # Feature colors
    GENE = QColor(38, 139, 210)  # Blue
    TRANSCRIPT = QColor(211, 54, 130)  # Magenta
    EXON = QColor(42, 161, 152)  # Cyan
    CDS = QColor(181, 137, 0)  # Yellow
    UTR = QColor(88, 110, 117)  # Gray
    INTRON = QColor(150, 150, 150)  # Gray for introns
    START_CODON = QColor(0, 200, 0)  # Bright green
    STOP_CODON = QColor(200, 0, 0)  # Bright red
    
    # Base colors (more muted for better readability)
    A = QColor(100, 140, 100)  # Muted green
    T = QColor(100, 140, 180)  # Muted blue
    G = QColor(180, 160, 100)  # Muted gold
    C = QColor(180, 100, 100)  # Muted crimson
    N = QColor(128, 128, 128)  # Gray
    GAP = QColor(200, 200, 200)  # Light gray
    
    # Sequence text color (light grey as requested)
    SEQUENCE_TEXT = QColor(150, 150, 150)  # Light grey
    
    # UI colors (dark mode)
    BACKGROUND = QColor(40, 44, 52)  # Dark background
    FOREGROUND = QColor(220, 220, 220)  # Light text
    HIGHLIGHT = QColor(70, 70, 70)  # Dark highlight
    PANEL_BACKGROUND = QColor(50, 54, 62)  # Panel background
    BORDER = QColor(80, 80, 80)  # Border color
    
    @classmethod
    def get_base_color(cls, base: str) -> QColor:
        """Get color for a DNA/RNA base."""
        base = base.upper()
        if base == 'A':
            return cls.A
        elif base == 'T' or base == 'U':
            return cls.T
        elif base == 'G':
            return cls.G
        elif base == 'C':
            return cls.C
        elif base == '-':
            return cls.GAP
        else:
            return cls.N
    
    @classmethod
    def get_variant_color(cls, ref: str, alt: str) -> QColor:
        """Get color based on variant type."""
        if len(ref) == len(alt):
            return cls.SNP
        elif len(ref) > len(alt):
            return cls.DELETION
        else:
            return cls.INSERTION


class FeatureTrackWidget(QWidget):
    """Custom widget for displaying feature tracks with banners."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.features = []
        self.window_size = 120
        self.position = 1000000
        self.show_reverse_strand = False
        self.setMinimumHeight(400)
        self.setMaximumHeight(400)
        
        # Calculate label offset to align with sequence viewer
        font = QFont("Courier New", 10)
        font.setStyleHint(QFont.StyleHint.Monospace)
        metrics = QFontMetrics(font)
        self.label_offset = max(
            metrics.horizontalAdvance("Position: "),
            metrics.horizontalAdvance("Ref:      "),
            metrics.horizontalAdvance("Personal: ")
        )
        
        # Feature hierarchy for display order
        self.feature_hierarchy = [
            'gene', 'transcript', 'exon', 'intron', 
            'CDS', 'five_prime_utr', 'three_prime_utr',
            'start_codon', 'stop_codon'
        ]
    
    def set_features(self, features: List[Dict], position: int, window_size: int, reverse_strand: bool = False):
        """Update the features to display."""
        self.features = features
        self.position = position
        self.window_size = window_size
        self.show_reverse_strand = reverse_strand
        self.update()
    
    def paintEvent(self, event):
        """Paint the feature tracks."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Background (dark mode)
        painter.fillRect(self.rect(), GenomeColors.BACKGROUND)
        
        if not self.features:
            painter.setPen(QPen(GenomeColors.FOREGROUND, 1))
            painter.drawText(self.label_offset, 30, "No features in this region")
            return
        
        # Calculate dimensions with label offset
        available_width = self.width() - self.label_offset - 10
        base_width = available_width // self.window_size if self.window_size > 0 else 1
        track_height = 20
        y_offset = 30
        
        # Draw position ruler
        painter.setPen(QPen(GenomeColors.BORDER, 1))
        painter.drawLine(self.label_offset, 20, self.label_offset + available_width, 20)
        
        # Draw tick marks
        for i in range(0, self.window_size, 10):
            x = self.label_offset + int(i * base_width)
            painter.drawLine(x, 18, x, 22)
            if i % 50 == 0:
                painter.setPen(QPen(GenomeColors.FOREGROUND, 1))
                # Note: This shows approximate positions, actual may vary with indels
                painter.drawText(x - 15, 15, str(self.position + i))
                painter.setPen(QPen(GenomeColors.BORDER, 1))
        
        # Group features by type
        feature_groups = {}
        for feature in self.features:
            ftype = feature.get('feature', 'unknown')
            if ftype == 'variant':
                continue  # Handle variants separately
            
            if ftype not in feature_groups:
                feature_groups[ftype] = []
            feature_groups[ftype].append(feature)
        
        # Add any missing feature types that are present
        for ftype in feature_groups.keys():
            if ftype not in self.feature_hierarchy:
                self.feature_hierarchy.append(ftype)
        
        # Draw features in hierarchical order
        for ftype in self.feature_hierarchy:
            if ftype not in feature_groups:
                continue
            
            color = self.get_feature_color(ftype)
            
            # Group by overlapping features to stack them
            tracks = []
            for feature in sorted(feature_groups[ftype], key=lambda f: f.get('start', 0)):
                feat_start = feature.get('start', 0)
                feat_end = feature.get('end', 0)
                
                # Find a track where this feature fits without overlap
                placed = False
                for track in tracks:
                    can_place = True
                    for existing in track:
                        exist_start = existing.get('start', 0)
                        exist_end = existing.get('end', 0)
                        if not (feat_end < exist_start or feat_start > exist_end):
                            can_place = False
                            break
                    if can_place:
                        track.append(feature)
                        placed = True
                        break
                
                if not placed:
                    tracks.append([feature])
            
            # Draw feature type label on the left (only once per type)
            if tracks:
                first_y = y_offset
                # Draw type label
                painter.setPen(QPen(GenomeColors.FOREGROUND, 1))
                font = QFont("Arial", 10, QFont.Weight.Bold)
                painter.setFont(font)
                
                # Format feature type name
                type_label = ftype.replace('_', ' ').title()
                if type_label == "Five Prime Utr":
                    type_label = "5' UTR"
                elif type_label == "Three Prime Utr":
                    type_label = "3' UTR"
                elif type_label == "Cds":
                    type_label = "CDS"
                
                # Draw label aligned to the right of the label area
                label_rect = painter.fontMetrics().boundingRect(type_label)
                painter.drawText(self.label_offset - label_rect.width() - 10, 
                               first_y + track_height // 2 + 5, type_label)
            
            # Draw each track
            for track_idx, track in enumerate(tracks):
                y = y_offset + track_idx * (track_height + 2)
                
                for feature in track:
                    feat_start = feature.get('start', 0)
                    feat_end = feature.get('end', 0)
                    info = feature.get('info', {})
                    
                    # Calculate display positions
                    rel_start = max(0, feat_start - self.position)
                    rel_end = min(self.window_size, feat_end - self.position + 1)
                    
                    if rel_end <= 0 or rel_start >= self.window_size:
                        continue
                    
                    x_start = self.label_offset + rel_start * base_width
                    x_end = self.label_offset + rel_end * base_width
                    
                    # Draw feature bar
                    painter.fillRect(x_start, y, x_end - x_start, track_height, color)
                    painter.setPen(QPen(color.darker(150), 1))
                    painter.drawRect(x_start, y, x_end - x_start, track_height)
                    
                    # Draw feature label
                    label = self.get_feature_label(ftype, info, feature)
                    painter.setPen(QPen(Qt.GlobalColor.white if color.lightness() < 150 else Qt.GlobalColor.black, 1))
                    font = QFont("Arial", 9)
                    painter.setFont(font)
                    
                    # Center label in feature bar
                    label_width = painter.fontMetrics().horizontalAdvance(label)
                    label_x = x_start + (x_end - x_start - label_width) / 2
                    if label_x > x_start and label_x + label_width < x_end:
                        painter.drawText(int(label_x), y + 14, label)
            
            y_offset += len(tracks) * (track_height + 2) + 5
        
        # Draw variants as flags
        self.draw_variants(painter, available_width, base_width)
    
    def draw_variants(self, painter: QPainter, available_width: float, base_width: float):
        """Draw variant flags."""
        variant_features = [f for f in self.features if f.get('feature') == 'variant']
        if not variant_features:
            return
        
        y_base = self.height() - 40
        
        for var in variant_features:
            var_start = var.get('start', 0)
            ref = var.get('ref', '')
            alt = var.get('alt', '')
            # qual = var.get('qual', 0)  # Commented to fix unused variable warning
            
            # Calculate position
            rel_pos = var_start - self.position
            if 0 <= rel_pos < self.window_size:
                x = self.label_offset + rel_pos * base_width
                
                # Determine variant type and color
                if len(ref) == len(alt):
                    color = GenomeColors.SNP
                    symbol = "●"
                elif len(ref) > len(alt):
                    color = GenomeColors.DELETION
                    symbol = "▼"
                else:
                    color = GenomeColors.INSERTION
                    symbol = "▲"
                
                # Draw flag pole
                painter.setPen(QPen(color, 2))
                painter.drawLine(x, y_base, x, y_base - 20)
                
                # Draw symbol
                painter.setPen(QPen(color, 3))
                font = QFont("Arial", 12)
                painter.setFont(font)
                painter.drawText(int(x - 5), y_base - 22, symbol)
                
                # Draw label
                label = f"{ref}→{alt}"
                if len(label) > 10:
                    label = f"VAR"
                font = QFont("Arial", 8)
                painter.setFont(font)
                painter.drawText(int(x - 10), y_base + 10, label)
    
    def get_feature_color(self, ftype: str) -> QColor:
        """Get color for a feature type."""
        color_map = {
            'gene': GenomeColors.GENE,
            'transcript': GenomeColors.TRANSCRIPT,
            'exon': GenomeColors.EXON,
            'CDS': GenomeColors.CDS,
            'five_prime_utr': GenomeColors.UTR,
            'three_prime_utr': GenomeColors.UTR,
            'intron': GenomeColors.INTRON,
            'start_codon': GenomeColors.START_CODON,
            'stop_codon': GenomeColors.STOP_CODON
        }
        return color_map.get(ftype, QColor(128, 128, 128))
    
    def get_feature_label(self, ftype: str, info: Dict, feature: Dict) -> str:
        """Get label for a feature."""
        if ftype == 'gene':
            return info.get('gene_name', 'gene')
        elif ftype == 'transcript':
            name = info.get('transcript_name', info.get('transcript_id', 'transcript'))
            if '-' in name and len(name) > 15:
                return name.split('-')[-1]
            return name[:15]
        elif ftype == 'exon':
            return f"E{info.get('exon_number', '?')}"
        elif ftype == 'CDS':
            return "CDS"
        elif ftype == 'five_prime_utr':
            return "5'UTR"
        elif ftype == 'three_prime_utr':
            return "3'UTR"
        elif ftype == 'intron':
            return "intron"
        elif ftype == 'start_codon':
            return "START"
        elif ftype == 'stop_codon':
            return "STOP"
        else:
            return ftype[:10]


class SequenceViewerWidget(QWidget):
    """Central widget for displaying DNA sequences with variants."""
    
    positionChanged = pyqtSignal(int)  # Emitted when position changes
    
    def __init__(self, genome_manager: 'GenomeManager', parent=None):
        super().__init__(parent)
        self.gm = genome_manager
        self.state = BrowserState(
            chrom='1',
            position=1000000,
            window_size=120,
            stride=40
        )
        
        self._current_features = []
        self._gene_cache = ""
        self._actual_end_position = 1000120  # Track actual genomic end position
        self.setup_ui()
        
    def setup_ui(self):
        """Set up the sequence viewer UI."""
        layout = QVBoxLayout(self)
        
        # Control bar
        control_layout = QHBoxLayout()
        
        # Calculate label width for proper alignment
        font = QFont("Courier New", 10)
        font.setStyleHint(QFont.StyleHint.Monospace)
        metrics = QFontMetrics(font)
        self.label_width = max(
            metrics.horizontalAdvance("Position: "),
            metrics.horizontalAdvance("Ref:      "),
            metrics.horizontalAdvance("Personal: ")
        )
        
        # Chromosome selector
        control_layout.addWidget(QLabel("Chr:"))
        self.chrom_combo = QComboBox()
        chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
        self.chrom_combo.addItems(chroms)
        self.chrom_combo.currentTextChanged.connect(self.on_chromosome_changed)
        control_layout.addWidget(self.chrom_combo)
        
        # Position input
        control_layout.addWidget(QLabel("Position:"))
        self.position_input = QLineEdit()
        self.position_input.setText(str(self.state.position))
        self.position_input.returnPressed.connect(self.on_position_changed)
        control_layout.addWidget(self.position_input)
        
        # Window size
        control_layout.addWidget(QLabel("Window:"))
        self.window_spin = QSpinBox()
        self.window_spin.setRange(40, 500)
        self.window_spin.setValue(self.state.window_size)
        self.window_spin.setSingleStep(20)
        self.window_spin.valueChanged.connect(self.on_window_changed)
        control_layout.addWidget(self.window_spin)
        
        # Display options
        self.show_features_cb = QCheckBox("Features")
        self.show_features_cb.setChecked(True)
        self.show_features_cb.toggled.connect(self.refresh_display)
        control_layout.addWidget(self.show_features_cb)
        
        self.show_amino_cb = QCheckBox("Amino Acids")
        self.show_amino_cb.toggled.connect(self.refresh_display)
        control_layout.addWidget(self.show_amino_cb)
        
        self.show_rna_cb = QCheckBox("RNA")
        self.show_rna_cb.toggled.connect(self.on_rna_toggled)
        control_layout.addWidget(self.show_rna_cb)
        
        self.reverse_strand_cb = QCheckBox("(-) Strand")
        self.reverse_strand_cb.toggled.connect(self.on_strand_toggled)
        control_layout.addWidget(self.reverse_strand_cb)
        
        control_layout.addStretch()
        layout.addLayout(control_layout)
        
        # Feature track display (new banner-style features)
        self.feature_track = FeatureTrackWidget()
        layout.addWidget(self.feature_track)
        
        # Sequence display area
        self.sequence_display = QTextEdit()
        self.sequence_display.setReadOnly(True)
        font = QFont("Courier New", 10)
        font.setStyleHint(QFont.StyleHint.Monospace)
        self.sequence_display.setFont(font)
        self.sequence_display.setMinimumHeight(200)
        self.sequence_display.setMaximumHeight(300)
        
        # Calculate character width for auto-sizing
        metrics = QFontMetrics(font)
        self.char_width = metrics.horizontalAdvance('A')
        
        layout.addWidget(self.sequence_display)
        
        # Feature details area (text listing)
        self.feature_display = QTextEdit()
        self.feature_display.setReadOnly(True)
        self.feature_display.setFont(font)
        self.feature_display.setMaximumHeight(150)
        layout.addWidget(self.feature_display)
        
        # Initial display
        self.refresh_display()
    
    def resizeEvent(self, event):
        """Handle window resize to adjust window size."""
        super().resizeEvent(event)
        self.adjust_window_size()
    
    def adjust_window_size(self):
        """Adjust window size based on available width."""
        if not hasattr(self, 'char_width') or not hasattr(self, 'label_width'):
            return
        
        # Get available width
        available_width = self.sequence_display.width() - 40  # Account for margins
        
        # Calculate how many characters can fit
        usable_width = available_width - self.label_width
        max_chars = max(40, int(usable_width // self.char_width))
        
        # Update window size if significantly different
        current_window = self.state.window_size
        if abs(max_chars - current_window) > 10:
            self.state.window_size = min(max_chars, 500)  # Cap at 500
            self.window_spin.setValue(self.state.window_size)
            self.refresh_display()
    
    def on_chromosome_changed(self, chrom: str):
        """Handle chromosome change."""
        self.state.chrom = chrom
        self.state.position = 1000000  # Reset to default position
        self.position_input.setText(str(self.state.position))
        self.refresh_display()
        self.positionChanged.emit(self.state.position)
    
    def on_position_changed(self):
        """Handle position change from input."""
        try:
            text = self.position_input.text().replace(',', '')
            
            # Try as number first
            try:
                pos = int(text)
                self.state.position = max(1, pos)
            except ValueError:
                # Try as gene name
                gene_name = text.upper()
                chrom, gene_info = self.gm.gene_map.find_gene(gene_name)
                
                if chrom and gene_info:
                    gene_data = gene_info[0]
                    self.state.chrom = str(chrom)
                    self.state.position = gene_data['start']
                    
                    # Update UI
                    self.chrom_combo.setCurrentText(str(chrom))
                    self.position_input.setText(str(self.state.position))
                else:
                    QMessageBox.warning(self, "Gene Not Found", 
                                       f"Gene '{gene_name}' not found")
                    return
            
            self.refresh_display()
            self.positionChanged.emit(self.state.position)
            
        except Exception as e:
            QMessageBox.warning(self, "Invalid Position", str(e))
    
    def on_window_changed(self, value: int):
        """Handle window size change."""
        self.state.window_size = value
        self.refresh_display()
    
    def on_rna_toggled(self, checked: bool):
        """Handle RNA mode toggle."""
        self.state.show_rna = checked
        self.refresh_display()
    
    def on_strand_toggled(self, checked: bool):
        """Handle strand toggle."""
        self.state.show_reverse_strand = checked
        self.refresh_display()
    
    def move_forward(self, large: bool = False):
        """Move forward in the genome using actual genomic end position."""
        move_distance = self.state.stride * (10 if large else 1)
        if self.state.show_reverse_strand:
            self.state.position = max(1, self.state.position - move_distance)
        else:
            # Use the actual end position to calculate next window start
            self.state.position = self._actual_end_position + 1
            if large:
                self.state.position += move_distance - self.state.stride
        
        self.position_input.setText(str(self.state.position))
        self.refresh_display()
        self.positionChanged.emit(self.state.position)
    
    def move_backward(self, large: bool = False):
        """Move backward in the genome."""
        move_distance = self.state.stride * (10 if large else 1)
        if self.state.show_reverse_strand:
            # Use actual end position for reverse strand
            self.state.position = self._actual_end_position + 1
            if large:
                self.state.position += move_distance - self.state.stride
        else:
            self.state.position = max(1, self.state.position - move_distance)
        
        self.position_input.setText(str(self.state.position))
        self.refresh_display()
        self.positionChanged.emit(self.state.position)
    
    def refresh_display(self):
        """Refresh the sequence display."""
        try:
            # Start with requested window size
            end_position = self.state.position + self.state.window_size - 1
            
            # Get sequences using GenomeIterator
            iterator = GenomeIterator(
                self.gm,
                self.state.chrom,
                self.state.position,
                end_position,
                window_size=self.state.window_size,
                integrate_variants=True,
                track_features=True
            )
            
            ref_seq, personal_seq, features = next(iterator)
            self._current_features = features
            self._cache_current_gene(features)
            
            # Adjust sequences to ensure window_size characters
            # If personal sequence is longer due to insertions, truncate
            # If shorter due to deletions, extend the window
            if len(personal_seq) > self.state.window_size:
                # Truncate to window size
                personal_seq = personal_seq[:self.state.window_size]
                ref_seq = ref_seq[:self.state.window_size]
                self._actual_end_position = end_position
            elif len(personal_seq) < self.state.window_size:
                # Need to fetch more bases to fill window
                needed = self.state.window_size - len(personal_seq)
                extended_end = end_position + needed
                
                # Try to get extended window
                try:
                    iterator2 = GenomeIterator(
                        self.gm,
                        self.state.chrom,
                        self.state.position,
                        extended_end,
                        window_size=self.state.window_size,
                        integrate_variants=True,
                        track_features=True
                    )
                    ref_seq, personal_seq, features = next(iterator2)
                    self._current_features = features
                    
                    # Now truncate to exact window size
                    if len(personal_seq) > self.state.window_size:
                        personal_seq = personal_seq[:self.state.window_size]
                        ref_seq = ref_seq[:len(personal_seq)]
                    
                    # Calculate actual genomic end position
                    self._actual_end_position = extended_end
                except:
                    # If we can't extend, use what we have
                    self._actual_end_position = end_position
            else:
                self._actual_end_position = end_position
            
            # Apply transformations
            if self.state.show_reverse_strand:
                ref_seq = reverse_complement(ref_seq, rna=False)
                personal_seq = reverse_complement(personal_seq, rna=False)
            
            if self.state.show_rna:
                ref_seq = to_rna(ref_seq)
                personal_seq = to_rna(personal_seq)
            
            # Display sequences
            self.display_sequences(ref_seq, personal_seq, features)
            
            # Update feature track widget with actual display window
            if self.show_features_cb.isChecked() and features:
                self.feature_track.set_features(features, self.state.position, 
                                               len(personal_seq),  # Use actual display length
                                               self.state.show_reverse_strand)
                self.feature_track.setVisible(True)
                self.display_features(features)
            else:
                self.feature_track.setVisible(False)
                self.feature_display.clear()
            
            # Update status with actual window info
            strand = '-' if self.state.show_reverse_strand else '+'
            molecule = 'RNA' if self.state.show_rna else 'DNA'
            status_text = (f"Chr {self.state.chrom} | "
                          f"Position {self.state.position:,}-{self._actual_end_position:,} | "
                          f"Display {len(personal_seq)}bp | "
                          f"({strand}) strand | {molecule}")
            
            if self.parent() and hasattr(self.parent(), 'statusBar'):
                self.parent().statusBar().showMessage(status_text)
                
        except StopIteration:
            self.sequence_display.setPlainText("No sequence data available at this position")
            self.feature_display.clear()
        except Exception as e:
            logger.error(f"Error refreshing display: {e}")
            self.sequence_display.setPlainText(f"Error: {e}")
    
    def display_sequences(self, ref_seq: str, personal_seq: str, features: List[Dict]):
        """Display the sequences with coloring."""
        if not ref_seq or not personal_seq:
            return
        
        # Clear display
        self.sequence_display.clear()
        cursor = self.sequence_display.textCursor()
        
        # Set default text format to light grey
        default_format = QTextCharFormat()
        default_format.setForeground(GenomeColors.SEQUENCE_TEXT)
        
        # Calculate actual genomic positions accounting for indels
        # Track the actual genomic position for each displayed character
        display_length = len(personal_seq)
        genomic_positions = []
        current_pos = self.state.position
        
        # Build position mapping (simplified - assumes personal seq is display seq)
        for i in range(display_length):
            genomic_positions.append(current_pos)
            # Advance position unless we're in a gap
            if i < len(ref_seq) and ref_seq[i] != '-':
                current_pos += 1
        
        # Position ruler with actual genomic positions
        cursor.setCharFormat(default_format)
        ruler_text = "Position: "
        for i in range(0, display_length, 10):
            if i < len(genomic_positions):
                ruler_text += f"{genomic_positions[i]:<10}"
        cursor.insertText(ruler_text[:display_length + 10] + "\n")
        
        ruler_line = "          "
        for i in range(display_length):
            if i % 10 == 0:
                ruler_line += "|"
            elif i % 5 == 0:
                ruler_line += "+"
            else:
                ruler_line += "."
        cursor.insertText(ruler_line + "\n\n")
        
        # Find start and stop codons in features for highlighting
        start_positions = set()
        stop_positions = set()
        if features:
            window_start = self.state.position
            for feature in features:
                ftype = feature.get('feature')
                if ftype == 'start_codon':
                    start = max(0, feature.get('start', 0) - window_start)
                    end = min(len(ref_seq), feature.get('end', 0) - window_start + 1)
                    for pos in range(start, end):
                        start_positions.add(pos)
                elif ftype == 'stop_codon':
                    start = max(0, feature.get('start', 0) - window_start)
                    end = min(len(ref_seq), feature.get('end', 0) - window_start + 1)
                    for pos in range(start, end):
                        stop_positions.add(pos)
        
        # Reference sequence (light grey)
        cursor.insertText("Ref:      ")
        for i, base in enumerate(ref_seq):
            char_format = QTextCharFormat()
            
            # Check for special codons
            if i in start_positions:
                char_format.setForeground(GenomeColors.START_CODON)
                char_format.setFontWeight(QFont.Weight.Bold)
            elif i in stop_positions:
                char_format.setForeground(GenomeColors.STOP_CODON)
                char_format.setFontWeight(QFont.Weight.Bold)
            else:
                char_format.setForeground(GenomeColors.SEQUENCE_TEXT)
            
            cursor.insertText(base, char_format)
        cursor.insertText("\n")
        
        # Personal sequence with variant highlighting
        cursor.insertText("Personal: ")
        for i, (ref_base, pers_base) in enumerate(zip(ref_seq, personal_seq)):
            char_format = QTextCharFormat()
            
            if ref_base != pers_base:
                # Variant - use bright colors with dark background
                char_format.setForeground(GenomeColors.get_variant_color(ref_base, pers_base))
                char_format.setFontWeight(QFont.Weight.Bold)
                char_format.setBackground(QColor(60, 60, 40))  # Dark yellow background
            elif i in start_positions:
                char_format.setForeground(GenomeColors.START_CODON)
                char_format.setFontWeight(QFont.Weight.Bold)
            elif i in stop_positions:
                char_format.setForeground(GenomeColors.STOP_CODON)
                char_format.setFontWeight(QFont.Weight.Bold)
            else:
                # Normal base - light grey
                char_format.setForeground(GenomeColors.SEQUENCE_TEXT)
            
            cursor.insertText(pers_base, char_format)
        cursor.insertText("\n")
        
        # Variant indicator line with symbols
        cursor.insertText("          ")
        for i in range(min(len(ref_seq), len(personal_seq))):
            if ref_seq[i] != personal_seq[i]:
                char_format = QTextCharFormat()
                char_format.setForeground(GenomeColors.get_variant_color(ref_seq[i], personal_seq[i]))
                char_format.setFontWeight(QFont.Weight.Bold)
                
                # Use different symbols for different variant types
                if ref_seq[i] == '-':
                    cursor.insertText("▲", char_format)  # Insertion
                elif personal_seq[i] == '-':
                    cursor.insertText("▼", char_format)  # Deletion
                else:
                    cursor.insertText("●", char_format)  # SNP
            else:
                cursor.insertText(" ")
        cursor.insertText("\n")
        
        # Amino acid translation if requested
        if self.show_amino_cb.isChecked():
            self.display_amino_acids(ref_seq, personal_seq, features, cursor)
    
    def display_amino_acids(self, ref_seq: str, personal_seq: str, 
                           features: List[Dict], cursor: QTextCursor):
        """Display amino acid translations for CDS regions."""
        # Find CDS features
        cds_features = [f for f in features if f.get('feature') == 'CDS']
        if not cds_features:
            return
        
        cursor.insertText("\n")
        window_start = self.state.position
        
        for cds in cds_features:
            cds_start = cds.get('start', 0)
            cds_end = cds.get('end', 0)
            strand = cds.get('strand', '+')
            
            # Calculate overlap with window
            overlap_start = max(cds_start, window_start)
            overlap_end = min(cds_end, window_start + len(ref_seq) - 1)
            
            if overlap_start > overlap_end:
                continue
            
            # Extract CDS portion
            window_cds_start = overlap_start - window_start
            window_cds_end = overlap_end - window_start + 1
            
            ref_coding = ref_seq[window_cds_start:window_cds_end]
            pers_coding = personal_seq[window_cds_start:window_cds_end]
            
            if strand == '-':
                ref_coding = reverse_complement(ref_coding)
                pers_coding = reverse_complement(pers_coding)
            
            # Translate
            ref_aa = self.translate_dna(ref_coding)
            pers_aa = self.translate_dna(pers_coding)
            
            # Display
            info = cds.get('info', {})
            gene_name = info.get('gene_name', 'Unknown')
            cursor.insertText(f"CDS ({gene_name}):\n")
            
            cursor.insertText("Ref AA:   " + " " * window_cds_start)
            for aa in ref_aa:
                cursor.insertText(f"{aa}  ")
            cursor.insertText("\n")
            
            cursor.insertText("Pers AA:  " + " " * window_cds_start)
            for ref_a, pers_a in zip(ref_aa, pers_aa):
                if ref_a != pers_a:
                    char_format = QTextCharFormat()
                    char_format.setForeground(GenomeColors.SNP)
                    char_format.setFontWeight(QFont.Weight.Bold)
                    cursor.insertText(f"{pers_a}  ", char_format)
                else:
                    cursor.insertText(f"{pers_a}  ")
            cursor.insertText("\n")
            break  # Show only first CDS
    
    def display_features(self, features: List[Dict]):
        """Display features in the feature display area."""
        self.feature_display.clear()
        cursor = self.feature_display.textCursor()
        
        cursor.insertText("Features:\n")
        
        # Group features by type
        feature_groups = {}
        for feature in features:
            ftype = feature.get('feature', 'unknown')
            if ftype not in feature_groups:
                feature_groups[ftype] = []
            feature_groups[ftype].append(feature)
        
        # Display by type
        for ftype, type_features in feature_groups.items():
            color = self.get_feature_color(ftype)
            
            char_format = QTextCharFormat()
            char_format.setForeground(color)
            
            for feature in type_features:
                info = feature.get('info', {})
                start = feature.get('start', 0)
                end = feature.get('end', 0)
                
                # Format feature info
                if ftype == 'gene':
                    name = info.get('gene_name', 'unknown')
                    text = f"  {ftype}: {name} ({start:,}-{end:,})"
                elif ftype == 'transcript':
                    name = info.get('transcript_name', info.get('transcript_id', 'unknown'))
                    text = f"  {ftype}: {name}"
                elif ftype == 'exon':
                    num = info.get('exon_number', '?')
                    text = f"  {ftype} {num}: {start:,}-{end:,}"
                elif ftype == 'variant':
                    ref = feature.get('ref', '')
                    alt = feature.get('alt', '')
                    qual = feature.get('qual', 0)
                    text = f"  {ftype}: {start:,} {ref}→{alt} (Q={qual:.1f})"
                else:
                    text = f"  {ftype}: {start:,}-{end:,}"
                
                cursor.insertText(text + "\n", char_format)
    
    def get_feature_color(self, ftype: str) -> QColor:
        """Get color for a feature type."""
        color_map = {
            'gene': GenomeColors.GENE,
            'transcript': GenomeColors.TRANSCRIPT,
            'exon': GenomeColors.EXON,
            'CDS': GenomeColors.CDS,
            'five_prime_utr': GenomeColors.UTR,
            'three_prime_utr': GenomeColors.UTR,
            'intron': GenomeColors.INTRON,
            'start_codon': GenomeColors.START_CODON,
            'stop_codon': GenomeColors.STOP_CODON,
            'variant': GenomeColors.SNP
        }
        return color_map.get(ftype, QColor(128, 128, 128))
    
    def translate_dna(self, dna_seq: str) -> str:
        """Translate DNA sequence to amino acids."""
        dna_seq = dna_seq.replace('-', '')
        amino_acids = []
        
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3].upper()
            if 'N' in codon or '-' in codon:
                amino_acids.append('X')
            else:
                amino_acids.append(CODON_TABLE_DNA.get(codon, 'X'))
        
        return ''.join(amino_acids)
    
    def _cache_current_gene(self, features):
        """Cache the current gene name."""
        for f in features:
            if f.get('feature') == 'gene':
                info = f.get('info', {})
                if 'gene_name' in info:
                    self._gene_cache = info['gene_name']
                    return
        self._gene_cache = ""


class CommandLineWidget(QWidget):
    """Command line widget at the bottom of the browser."""
    
    commandExecuted = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()
        self.command_history = []
        self.history_index = -1
    
    def setup_ui(self):
        """Set up the command line UI."""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.prompt_label = QLabel("Command:")
        layout.addWidget(self.prompt_label)
        
        self.command_input = QLineEdit()
        self.command_input.returnPressed.connect(self.execute_command)
        layout.addWidget(self.command_input)
        
        # Set up command history navigation
        self.command_input.installEventFilter(self)
    
    def eventFilter(self, obj, event):
        """Handle up/down arrow keys for command history."""
        if obj == self.command_input and event.type() == event.Type.KeyPress:
            key = event.key()
            if key == Qt.Key.Key_Up:
                self.navigate_history(-1)
                return True
            elif key == Qt.Key.Key_Down:
                self.navigate_history(1)
                return True
        return super().eventFilter(obj, event)
    
    def navigate_history(self, direction: int):
        """Navigate command history."""
        if not self.command_history:
            return
        
        self.history_index = max(0, min(len(self.command_history) - 1, 
                                       self.history_index + direction))
        self.command_input.setText(self.command_history[self.history_index])
    
    def execute_command(self):
        """Execute the entered command."""
        command = self.command_input.text().strip()
        if not command:
            return
        
        # Add to history
        self.command_history.append(command)
        self.history_index = len(self.command_history)
        
        # Clear input
        self.command_input.clear()
        
        # Emit signal
        self.commandExecuted.emit(command)


class FeatureTreeWidget(QTreeWidget):
    """Tree widget for displaying genomic features."""
    
    featureSelected = pyqtSignal(str, str, int)  # chrom, feature_type, position
    
    def __init__(self, genome_manager: 'GenomeManager', parent=None):
        super().__init__(parent)
        self.gm = genome_manager
        self.setup_ui()
    
    def setup_ui(self):
        """Set up the feature tree UI."""
        self.setHeaderLabels(["Feature", "Count"])
        self.itemDoubleClicked.connect(self.on_item_double_clicked)
    
    def update_features(self, chrom: str, position: int):
        """Update the feature tree for the current position."""
        self.clear()
        
        # Get features at current position
        features = self.gm.get_features_at_position(chrom, position)
        
        # Group by type
        feature_groups = {}
        for feature in features:
            ftype = feature.get('feature', 'unknown')
            if ftype not in feature_groups:
                feature_groups[ftype] = []
            feature_groups[ftype].append(feature)
        
        # Create tree items
        for ftype, type_features in feature_groups.items():
            type_item = QTreeWidgetItem(self, [ftype, str(len(type_features))])
            
            for feature in type_features:
                info = feature.get('info', {})
                name = self.get_feature_name(ftype, info)
                start = feature.get('start', 0)
                end = feature.get('end', 0)
                
                feature_item = QTreeWidgetItem(type_item, 
                    [f"{name} ({start:,}-{end:,})", ""])
                feature_item.setData(0, Qt.ItemDataRole.UserRole, 
                                   (chrom, ftype, start))
    
    def get_feature_name(self, ftype: str, info: Dict) -> str:
        """Get a display name for a feature."""
        if ftype == 'gene':
            return info.get('gene_name', 'unknown')
        elif ftype == 'transcript':
            return info.get('transcript_name', info.get('transcript_id', 'unknown'))
        elif ftype == 'exon':
            return f"Exon {info.get('exon_number', '?')}"
        else:
            return ftype
    
    def on_item_double_clicked(self, item: QTreeWidgetItem, column: int):
        """Handle double-click on a feature."""
        data = item.data(0, Qt.ItemDataRole.UserRole)
        if data:
            chrom, ftype, position = data
            self.featureSelected.emit(chrom, ftype, position)


class ChromosomeMapWidget(QGraphicsView):
    """Widget for displaying a chromosome map with features."""
    
    positionClicked = pyqtSignal(int)  # Emitted when user clicks on a position
    geneClicked = pyqtSignal(str, int, int)  # Emitted when clicking a gene (name, start, end)
    
    def __init__(self, genome_manager: 'GenomeManager', parent=None):
        super().__init__(parent)
        self.gm = genome_manager
        self.current_chrom = '1'
        self.current_position = 1000000
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.genes_data = {}  # Cache for chromosome gene data
        self.gene_items = []  # Track gene graphic items for interaction
        self.setup_ui()
        self.load_chromosome_genes()
    
    def setup_ui(self):
        """Set up the chromosome map UI."""
        self.setMinimumWidth(200)
        self.setMaximumWidth(300)
        self.setMouseTracking(True)  # Enable hover events
        
    def update_chromosome(self, chrom: str, position: int):
        """Update the chromosome map display."""
        if self.current_chrom != chrom:
            self.current_chrom = chrom
            self.load_chromosome_genes()
        self.current_position = position
        self.draw_chromosome()
    
    def load_chromosome_genes(self):
        """Load gene data from JSON file for current chromosome."""
        import json
        import os
        
        # Construct JSON file path
        json_path = os.path.join('data', f'chrom{self.current_chrom}_map.json')
        
        try:
            if os.path.exists(json_path):
                with open(json_path, 'r') as f:
                    all_genes = json.load(f)
                    # Filter for protein_coding genes only
                    self.genes_data = {
                        name: info for name, info in all_genes.items()
                        if info.get('biotype') == 'protein_coding'
                    }
                    logger.info(f"Loaded {len(self.genes_data)} protein-coding genes for chromosome {self.current_chrom}")
            else:
                self.genes_data = {}
                logger.warning(f"No gene map file found at {json_path}")
        except Exception as e:
            logger.error(f"Error loading chromosome gene data: {e}")
            self.genes_data = {}
    
    def draw_chromosome(self):
        """Draw the chromosome with major features."""
        self.scene.clear()
        self.gene_items = []
        
        # Get chromosome length (approximate)
        chrom_lengths = {
            '1': 249250621, '2': 242193529, '3': 198295559,
            '4': 190214555, '5': 181538259, '6': 170805979,
            '7': 159345973, '8': 145138636, '9': 138394717,
            '10': 133797422, '11': 135086622, '12': 133275309,
            '13': 114364328, '14': 107043718, '15': 101991189,
            '16': 90338345, '17': 83257441, '18': 80373285,
            '19': 58617616, '20': 64444167, '21': 46709983,
            '22': 50818468, 'X': 156040895, 'Y': 57227415,
            'MT': 16569
        }
        
        chrom_length = chrom_lengths.get(self.current_chrom, 250000000)
        
        # Draw chromosome outline (dark mode)
        width = 150
        height = self.height() - 60
        
        chrom_rect = self.scene.addRect(25, 20, width, height, 
                          QPen(GenomeColors.BORDER), 
                          QBrush(GenomeColors.PANEL_BACKGROUND))
        
        # Draw genes as notches on the chromosome
        if self.genes_data:
            for gene_name, gene_info in self.genes_data.items():
                gene_start = gene_info.get('start', 0)
                gene_end = gene_info.get('end', 0)
                strand = gene_info.get('strand', '+')
                
                # Calculate Y position on chromosome
                gene_y = 20 + (gene_start / chrom_length) * height
                
                # Draw gene notch (left side for +, right side for -)
                if strand == '+':
                    notch_x = 20
                    notch_width = 5
                else:
                    notch_x = 175
                    notch_width = 5
                
                # Create clickable gene notch
                gene_rect = self.scene.addRect(notch_x, gene_y, notch_width, 2,
                                              QPen(GenomeColors.GENE),
                                              QBrush(GenomeColors.GENE))
                
                # Store gene info with the item
                gene_rect.setData(0, {'name': gene_name, 'info': gene_info})
                gene_rect.setToolTip(f"{gene_name}\n{gene_start:,}-{gene_end:,}\nStrand: {strand}")
                gene_rect.setCursor(Qt.CursorShape.PointingHandCursor)
                self.gene_items.append(gene_rect)
        
        # Draw current position indicator (on top)
        pos_y = 20 + (self.current_position / chrom_length) * height
        position_indicator = self.scene.addRect(10, pos_y - 2, width + 30, 4,
                          QPen(GenomeColors.SNP),
                          QBrush(GenomeColors.SNP))
        position_indicator.setZValue(10)  # Ensure it's on top
        
        # Add chromosome label
        label = self.scene.addText(f"Chr {self.current_chrom}")
        label.setDefaultTextColor(GenomeColors.FOREGROUND)
        label.setPos(70, 0)
        
        # Add position label
        pos_label = self.scene.addText(f"{self.current_position:,}")
        pos_label.setDefaultTextColor(GenomeColors.FOREGROUND)
        pos_label.setPos(50, pos_y + 5)
        
        # Add gene count label
        if self.genes_data:
            count_label = self.scene.addText(f"{len(self.genes_data)} genes")
            count_label.setDefaultTextColor(GenomeColors.FOREGROUND)
            count_label.setPos(60, height + 25)
        
        # Fit view
        self.fitInView(self.scene.itemsBoundingRect(), 
                       Qt.AspectRatioMode.KeepAspectRatio)
    
    def mousePressEvent(self, event):
        """Handle mouse clicks on the chromosome map."""
        scene_pos = self.mapToScene(event.pos())
        
        # Check if clicking on a gene item
        item = self.scene.itemAt(scene_pos, self.transform())
        if item and item in self.gene_items:
            gene_data = item.data(0)
            if gene_data:
                gene_name = gene_data['name']
                gene_info = gene_data['info']
                gene_start = gene_info.get('start', 0)
                gene_end = gene_info.get('end', 0)
                self.geneClicked.emit(gene_name, gene_start, gene_end)
                return
        
        # Otherwise, convert y position to genomic position
        if 20 <= scene_pos.y() <= self.height() - 60:
            relative_y = (scene_pos.y() - 20) / (self.height() - 60)
            
            chrom_lengths = {
                '1': 249250621, '2': 242193529, '3': 198295559,
                '4': 190214555, '5': 181538259, '6': 170805979,
                '7': 159345973, '8': 145138636, '9': 138394717,
                '10': 133797422, '11': 135086622, '12': 133275309,
                '13': 114364328, '14': 107043718, '15': 101991189,
                '16': 90338345, '17': 83257441, '18': 80373285,
                '19': 58617616, '20': 64444167, '21': 46709983,
                '22': 50818468, 'X': 156040895, 'Y': 57227415,
                'MT': 16569
            }
            chrom_length = chrom_lengths.get(self.current_chrom, 250000000)
            
            position = int(relative_y * chrom_length)
            self.positionClicked.emit(position)
        
        super().mousePressEvent(event)


class GenomeBrowserMainWindow(QMainWindow):
    """Main window for the genome browser application."""
    
    def __init__(self, genome_manager: 'GenomeManager'):
        super().__init__()
        self.gm = genome_manager
        self.setup_ui()
        self.setup_shortcuts()
    
    def setup_ui(self):
        """Set up the main window UI."""
        self.setWindowTitle("Genome Browser")
        self.setGeometry(100, 100, 1400, 900)
        
        # Central widget with main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel: Feature tree
        left_dock = QDockWidget("Features", self)
        self.feature_tree = FeatureTreeWidget(self.gm)
        self.feature_tree.featureSelected.connect(self.on_feature_selected)
        left_dock.setWidget(self.feature_tree)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, left_dock)
        
        # Center: Main sequence viewer
        center_layout = QVBoxLayout()
        self.sequence_viewer = SequenceViewerWidget(self.gm)
        self.sequence_viewer.positionChanged.connect(self.on_position_changed)
        center_layout.addWidget(self.sequence_viewer)
        
        # Command line at bottom
        self.command_line = CommandLineWidget()
        self.command_line.commandExecuted.connect(self.execute_command)
        center_layout.addWidget(self.command_line)
        
        main_layout.addLayout(center_layout)
        
        # Right panel: Chromosome map
        right_dock = QDockWidget("Chromosome Map", self)
        self.chrom_map = ChromosomeMapWidget(self.gm)
        self.chrom_map.positionClicked.connect(self.on_map_position_clicked)
        self.chrom_map.geneClicked.connect(self.on_gene_clicked)
        right_dock.setWidget(self.chrom_map)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, right_dock)
        
        # Status bar
        self.statusBar().showMessage("Ready")
        
        # Toolbar
        self.create_toolbar()
        
        # Initial update
        self.on_position_changed(self.sequence_viewer.state.position)
    
    def create_toolbar(self):
        """Create the main toolbar."""
        toolbar = self.addToolBar("Navigation")
        
        # Navigation actions
        backward_action = QAction("← Back", self)
        backward_action.triggered.connect(lambda: self.sequence_viewer.move_backward())
        toolbar.addAction(backward_action)
        
        forward_action = QAction("Forward →", self)
        forward_action.triggered.connect(lambda: self.sequence_viewer.move_forward())
        toolbar.addAction(forward_action)
        
        toolbar.addSeparator()
        
        # Jump actions
        jump_back_action = QAction("⇐ Jump Back", self)
        jump_back_action.triggered.connect(lambda: self.sequence_viewer.move_backward(True))
        toolbar.addAction(jump_back_action)
        
        jump_forward_action = QAction("Jump Forward ⇒", self)
        jump_forward_action.triggered.connect(lambda: self.sequence_viewer.move_forward(True))
        toolbar.addAction(jump_forward_action)
        
        toolbar.addSeparator()
        
        # Gene card action
        genecard_action = QAction("Open GeneCard", self)
        genecard_action.triggered.connect(self.open_genecard)
        toolbar.addAction(genecard_action)
    
    def setup_shortcuts(self):
        """Set up keyboard shortcuts."""
        # Navigation
        QShortcut(QKeySequence(Qt.Key.Key_Right), self, 
                 lambda: self.sequence_viewer.move_forward())
        QShortcut(QKeySequence(Qt.Key.Key_Left), self,
                 lambda: self.sequence_viewer.move_backward())
        QShortcut(QKeySequence(Qt.Key.Key_Up), self,
                 lambda: self.sequence_viewer.move_forward(True))
        QShortcut(QKeySequence(Qt.Key.Key_Down), self,
                 lambda: self.sequence_viewer.move_backward(True))
        
        # Features
        QShortcut(QKeySequence('f'), self,
                 lambda: self.sequence_viewer.show_features_cb.toggle())
        QShortcut(QKeySequence('a'), self,
                 lambda: self.sequence_viewer.show_amino_cb.toggle())
        QShortcut(QKeySequence('r'), self,
                 lambda: self.sequence_viewer.show_rna_cb.toggle())
        
        # Jump to position
        QShortcut(QKeySequence('g'), self, self.goto_position)
        
        # Help
        QShortcut(QKeySequence('?'), self, self.show_help)
    
    def on_position_changed(self, position: int):
        """Handle position change from sequence viewer."""
        chrom = self.sequence_viewer.state.chrom
        
        # Update feature tree
        self.feature_tree.update_features(chrom, position)
        
        # Update chromosome map
        self.chrom_map.update_chromosome(chrom, position)
    
    def on_feature_selected(self, chrom: str, ftype: str, position: int):
        """Handle feature selection from tree."""
        self.sequence_viewer.state.chrom = chrom
        self.sequence_viewer.state.position = position
        self.sequence_viewer.chrom_combo.setCurrentText(chrom)
        self.sequence_viewer.position_input.setText(str(position))
        self.sequence_viewer.refresh_display()
    
    def on_map_position_clicked(self, position: int):
        """Handle position click from chromosome map."""
        self.sequence_viewer.state.position = position
        self.sequence_viewer.position_input.setText(str(position))
        self.sequence_viewer.refresh_display()
    
    def on_gene_clicked(self, gene_name: str, start: int, end: int):
        """Handle gene click from chromosome map."""
        # Jump to gene start position
        self.sequence_viewer.state.position = start
        self.sequence_viewer.position_input.setText(str(start))
        self.sequence_viewer.refresh_display()
        
        # Show message in status bar
        self.statusBar().showMessage(f"Jumped to {gene_name} ({start:,}-{end:,})")
    
    def execute_command(self, command: str):
        """Execute a command from the command line."""
        parts = command.split()
        if not parts:
            return
        
        cmd = parts[0].lower()
        
        if cmd == 'goto' and len(parts) > 1:
            try:
                pos = int(parts[1].replace(',', ''))
                self.sequence_viewer.state.position = pos
                self.sequence_viewer.position_input.setText(str(pos))
                self.sequence_viewer.refresh_display()
            except ValueError:
                self.statusBar().showMessage(f"Invalid position: {parts[1]}")
        
        elif cmd == 'chr' and len(parts) > 1:
            chrom = parts[1].upper()
            if chrom in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']:
                self.sequence_viewer.state.chrom = chrom
                self.sequence_viewer.chrom_combo.setCurrentText(chrom)
                self.sequence_viewer.refresh_display()
            else:
                self.statusBar().showMessage(f"Invalid chromosome: {chrom}")
        
        elif cmd == 'find' and len(parts) > 1:
            gene_name = parts[1].upper()
            chrom, gene_info = self.gm.gene_map.find_gene(gene_name)
            if chrom and gene_info:
                gene_data = gene_info[0]
                self.sequence_viewer.state.chrom = str(chrom)
                self.sequence_viewer.state.position = gene_data['start']
                self.sequence_viewer.chrom_combo.setCurrentText(str(chrom))
                self.sequence_viewer.position_input.setText(str(gene_data['start']))
                self.sequence_viewer.refresh_display()
                self.statusBar().showMessage(f"Found {gene_name}")
            else:
                self.statusBar().showMessage(f"Gene not found: {gene_name}")
        
        elif cmd == 'help' or cmd == '?':
            self.show_help()
        
        else:
            self.statusBar().showMessage(f"Unknown command: {cmd}")
    
    def goto_position(self):
        """Show dialog to jump to position."""
        text, ok = QInputDialog.getText(self, "Go to Position", 
                                        "Enter position or gene name:")
        if ok and text:
            self.sequence_viewer.position_input.setText(text)
            self.sequence_viewer.on_position_changed()
    
    def open_genecard(self):
        """Open GeneCard for current gene."""
        gene_name = self.sequence_viewer._gene_cache
        if gene_name:
            url = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_name}"
            subprocess.run(["xdg-open", url])
        else:
            self.statusBar().showMessage("No gene at current position")
    
    def show_help(self):
        """Show help dialog."""
        help_text = """
        <h3>Genome Browser Help</h3>
        
        <h4>Navigation:</h4>
        <ul>
        <li><b>→ / Left Arrow:</b> Move forward</li>
        <li><b>← / Right Arrow:</b> Move backward</li>
        <li><b>↑ / Up Arrow:</b> Jump forward (10x)</li>
        <li><b>↓ / Down Arrow:</b> Jump backward (10x)</li>
        <li><b>g:</b> Go to position or gene</li>
        </ul>
        
        <h4>Display Options:</h4>
        <ul>
        <li><b>f:</b> Toggle features</li>
        <li><b>a:</b> Toggle amino acids</li>
        <li><b>r:</b> Toggle RNA mode</li>
        </ul>
        
        <h4>Commands:</h4>
        <ul>
        <li><b>goto [position]:</b> Jump to position</li>
        <li><b>chr [chromosome]:</b> Switch chromosome</li>
        <li><b>find [gene]:</b> Find gene</li>
        <li><b>help:</b> Show this help</li>
        </ul>
        """
        
        QMessageBox.information(self, "Help", help_text)


def launch_genome_browser(genome_manager: 'GenomeManager'):
    """Launch the PyQt genome browser application."""
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Set dark theme palette
    palette = QPalette()
    
    # Window colors
    palette.setColor(QPalette.ColorRole.Window, GenomeColors.BACKGROUND)
    palette.setColor(QPalette.ColorRole.WindowText, GenomeColors.FOREGROUND)
    
    # Base colors (for input fields, etc.)
    palette.setColor(QPalette.ColorRole.Base, GenomeColors.PANEL_BACKGROUND)
    palette.setColor(QPalette.ColorRole.AlternateBase, GenomeColors.BACKGROUND)
    
    # Text colors
    palette.setColor(QPalette.ColorRole.Text, GenomeColors.FOREGROUND)
    palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 255, 255))
    
    # Button colors
    palette.setColor(QPalette.ColorRole.Button, GenomeColors.PANEL_BACKGROUND)
    palette.setColor(QPalette.ColorRole.ButtonText, GenomeColors.FOREGROUND)
    
    # Highlight colors
    palette.setColor(QPalette.ColorRole.Highlight, GenomeColors.HIGHLIGHT)
    palette.setColor(QPalette.ColorRole.HighlightedText, GenomeColors.FOREGROUND)
    
    # Disabled colors
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, QColor(120, 120, 120))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(120, 120, 120))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(120, 120, 120))
    
    app.setPalette(palette)
    
    # Create and show main window
    window = GenomeBrowserMainWindow(genome_manager)
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    # Test launch
    print("This module should be imported and used with a GenomeManager instance")
    print("Use: launch_genome_browser(genome_manager)")