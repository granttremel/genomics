#!/usr/bin/env python3
"""
Interactive genome browser for visualizing personal variants and features.
"""
import subprocess
import sys
import tty
import termios
import os
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
from dataclasses import dataclass
import logging
import numpy as np

from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, COMPLEMENT_MAP, COMPLEMENT_MAP_RNA, to_rna, complement, reverse_complement
from ggene.features import Gene, Feature
from ggene import seqs, draw
from ggene.genome_iterator_v2 import UnifiedGenomeIteratorfrom ggene import draw, seqs
from ggene.unified_stream import UnifiedFeature
from ggene.genome_iterator_v2 import UnifiedGenomeIterator

if TYPE_CHECKING:
    from .genomemanager import GenomeManager

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

# ANSI color codes
class Colors:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'
    
    # Variant colors
    SNP = '\033[91m'        # Red for SNPs
    INSERTION = '\033[92m'   # Green for insertions
    DELETION = '\033[93m'    # Yellow for deletions
    
    # Feature colors
    GENE = '\033[94m'        # Blue
    TRANSCRIPT = '\033[95m'  # Magenta
    EXON = '\033[96m'        # Cyan
    CDS = '\033[93m'         # Yellow
    UTR = '\033[90m'         # Gray
    
    # Motif colors (for underlines)
    MOTIF_SPLICE = '\033[91m'    # Red for splice sites
    MOTIF_TATA = '\033[95m'      # Magenta for TATA box
    MOTIF_CPG = '\033[92m'       # Green for CpG islands
    MOTIF_POLYA = '\033[94m'     # Blue for polyA signals
    MOTIF_DEFAULT = '\033[96m'   # Cyan for other motifs
    
    # Navigation
    POSITION = '\033[97m'    # White
    
    @classmethod
    def variant_color(cls, ref: str, alt: str) -> str:
        """Get color based on variant type."""
        if len(ref) == len(alt):
            return cls.SNP
        elif len(ref) > len(alt):
            return cls.DELETION
        else:
            return cls.INSERTION


@dataclass
class BrowserState:
    """Current state of the genome browser."""
    chrom: Union[str, int]
    position: int
    window_size: int
    stride: int
    show_features: bool = True
    show_quality: bool = False
    show_amino_acids: bool = False
    show_reverse_strand: bool = False
    show_rna: bool = False
    feature_types: Optional[List[str]] = None
    show_second:bool = False
    
    def copy(self):
        return type(self)(**self.__dict__)
    
class MotifAnnotationLayer:
    """Live motif detection for genome browser."""

    def __init__(self):
        self.fast_cache = {}  # Cache recent windows
        self.background_model = None  # For significance testing

    def annotate_window(self, seq, pos, features):
        """Return motifs as feature-like objects."""
        motifs = []

        # Quick patterns (every window)
        for pattern_motif in self.pattern_motifs:
            if hit := pattern_motif.find_in_window(seq):
                motifs.append(self._to_feature(hit, pos))

        # Statistical (sample occasionally)  
        if pos % 100 == 0:
            gc = self.gc_content(seq)
            if gc > 0.7:  # CpG island threshold
                motifs.append(Feature(type='cpg_island'))

        return motifs

class InteractiveGenomeBrowser:
    """Interactive terminal-based genome browser."""
    
    genecard="https://www.genecards.org/cgi-bin/carddisp.pl?gene={name}"
    wiki="https://en.wikipedia.org/wiki/{name}"
    
    def __init__(self, genome_manager: 'GenomeManager', debug = False):
        """Initialize the genome browser.
        
        Args:
            genome_manager: Initialized GenomeManager instance
        """
        self.gm = genome_manager
        self.state = BrowserState(1, int(10e6), 128, 20)
        self.state2 = BrowserState(1, int(10e6), 128, 20)
        self.history = []
        self.iterator = None
        self.debug = debug
        if self.debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
        
        self.iterator = UnifiedGenomeIterator(
            self.gm,
            self.state.chrom,
            self.state.position,
            window_size = self.state.window_size,
            stride = self.state.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
        )
        self.iterator2 = UnifiedGenomeIterator(
            self.gm,
            self.state2.chrom,
            self.state2.position,
            window_size = self.state2.window_size,
            stride = self.state2.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
        )
        
        self._ref_cache1 = self._ref_cache2 = self._alt_cache1 = self._alt_cache2 = ""
        self._current_features={}
        self._gene_cache=""
        
        
        self.data = ["autocorrelation"]
        
        # Feature display order (based on hierarchy)
        self.feature_hierarchy = [
            'gene', 'transcript', 'exon', 'intron', 
            'CDS', 'five_prime_utr', 'three_prime_utr',
            'start_codon', 'stop_codon', 
            'tf_binding', 'motif', 'splice_donor', 'splice_acceptor',
            'repeat', 'chip_peak', 'variant'
        ]
        
        self.arrows = draw.get_arrow("default", 5, "-") # L, R
    
    def start(self, chrom: Union[str, int], position: int = 1,
             window_size: int = 240, stride: int = 40, show_second = False):
        """Start the interactive browser.
        
        Args:
            chrom: Chromosome to browse
            position: Starting position
            window_size: Size of sequence window to display
            stride: How many bases to move on arrow key
        """
        self.state = BrowserState(
            chrom=chrom,
            position=position,
            window_size=window_size,
            stride=stride,
            show_second = show_second,
        )
        if show_second:
            self.state2 = self.state.copy()
        self.iterator = UnifiedGenomeIterator(
            self.gm,
            self.state.chrom,
            self.state.position,
            window_size = self.state.window_size,
            stride = self.state.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
            
        )
        # Clear screen
        if not self.debug:
            os.system('clear' if os.name == 'posix' else 'cls')
        
        print(f"{Colors.BOLD}Interactive Genome Browser{Colors.RESET}")
        print(f"Navigate with arrow keys, 'q' to quit, 'h' for help")
        print("-" * 80)
        
        try:
            self._run_browser()
        except KeyboardInterrupt:
            self._cleanup()
    
    def _render_pos(self,start, end):
        if self.state.show_reverse_strand:
            return self.state.window_size-end, self.state.window_size-start
        else:
            return start,end
        
    def _render_seq(self, seq):
        if self.state.show_reverse_strand:
            seq = reverse_complement(seq)
        if self.state.show_rna:
            seq = to_rna(seq)
        return seq
        
    def _run_browser(self):
        """Main browser loop."""
        while True:
            # Display current view
            self._update_iterator()
            if self.state.show_second:
                self._display_current_view_sec()
            else:
                self._display_current_view2()
            # Get user input
            key = self._get_key()
            
            # Handle navigation
            if key == 'q':
                break
            elif key == '\x1b[C':  # Right arrow or 'l'
                self._move_forward()
            elif key == '\x1b[D':  # Left arrow or 'h'
                self._move_backward()
            elif key == '\x1b[A':  # Up arrow or 'k'
                self._move_forward(large=True)
            elif key == '\x1b[B':  # Down arrow or 'j'
                self._move_backward(large=True)
            elif key == 'shift_right':  # Right arrow or 'l'
                self._move_forward(second=True)
            elif key == 'shift_left':  # Left arrow or 'h'
                self._move_backward(second=True)
            elif key == 'shift_up':  # Up arrow or 'k'
                self._move_forward(second=True,large=True)
            elif key == 'shift_down':  # Down arrow or 'j'
                self._move_backward(second=True,large=True)
            elif key == 'ctrl_right':  
                self._jump_to_next_gene()
            elif key == 'ctrl_left':  
                self._jump_to_prev_gene()
            elif key == ' ' or key == '\n' or key == '\r':  # Space or Enter
                self._move_forward()
            elif key == 'g':  # Go to position
                self._goto_position()
            elif key == 'G':  # Go to position
                self._goto_position(second=True)
            elif key == 'n':
                self._next_feature()
            elif key == 'p':
                self._prev_feature()
            elif key == 'c':
                self._switch_chromosome()
            elif key == 'C':
                self._switch_chromosome(second=True)
            elif key == 'f':  # Toggle features
                self.state.show_features = not self.state.show_features
            elif key == 'a':  # Toggle amino acid display
                self.state.show_amino_acids = not self.state.show_amino_acids
            elif key == 'r':  # Toggle RNA mode
                self.state.show_rna = not self.state.show_rna
            elif key == '-':  # Toggle reverse strand
                self.state.show_reverse_strand = not self.state.show_reverse_strand
            elif key == 'w':  # Change window size
                self._change_window_size()
            elif key == 's':  # Change stride
                self._change_stride()
            elif key == '?':  # Help
                self._show_help()
            elif key == 'i':  # Show position info
                self._show_position_info()
            elif key == 'l':
                self.gm.assemble_gene(self._gene_cache, self.state.chrom)
            elif key == 'o':
                self._open_gene_card()
            elif key == 'k':
                self._open_wiki()
            elif key == 't':
                self._translate_transcript()
                # self.gm.ribo.translate_transcript()
            elif key == 'ctrl_right':  # Ctrl+Right: Next gene
                self._jump_to_next_gene()
            elif key == 'ctrl_left':   # Ctrl+Left: Previous gene  
                self._jump_to_prev_gene()
            elif key == 'ctrl_up':     # Ctrl+Up: Save state
                self._save_state()
            elif key == 'ctrl_down':   # Ctrl+Down: Load state
                self._load_state()
            elif key == 'd':
                if not self.state.show_second:
                    self.initialize_second_sequence()
                else:
                    self.state.show_second = False
                    self.state2.show_second = False
        
        self._cleanup()
    
    def initialize_second_sequence(self):
        print("initializing second sequence:")
        c=input("chromosome: ")
        pos=input("position: ")
        if not pos:
            pos = self.state.position
        else:
            pos=int(pos)
        
        self.state2 = self.state.copy()
        if c in list(map(str,range(23))) + ['MT','X','Y']:
            self.state2.chrom=c
        else:
            self.state2.chrom = self.state.chrom
        self.state2.position=pos
        
        self.iterator2 = UnifiedGenomeIterator(
            self.gm,
            self.state2.chrom,
            self.state2.position,
            window_size = self.state2.window_size,
            stride = self.state2.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
        )
        
        self.state.show_second= self.state2.show_second = True
    
    # def update_iterator(self):
        
    #     iters = [self.iterator]
    #     states = [self.state]
    #     if self.state.show_second:
    #         iters.append(self.iterator2)
    #         states.append(self.state2)
        
    #     for iterator, state in zip(iters, states):
    #         iterator.start = state.position
    #         iterator.end = state.position + state.window_size
    #         iterator.window_size = state.window_size
    #         iterator.stride = state.stride
    #         iterator._preload_buffer()
        
    #         logger.debug(f"iterator updated: {repr(iterator)}")
    
    def _update_iterator(self):
        
        self.iterator = UnifiedGenomeIterator(
            self.gm,
            self.state.chrom,
            self.state.position,
            window_size = self.state.window_size,
            stride = self.state.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
        )
        self.iterator2 = UnifiedGenomeIterator(
            self.gm,
            self.state2.chrom,
            self.state2.position,
            window_size = self.state2.window_size,
            stride = self.state2.stride,
            integrate_variants=True,
            track_features = True,
            feature_types=["gene","exon","CDS","motif"]
        )
    
    def _display_current_view2(self):
        
        if not self.debug:
            print('\033[2J\033[H')  # Clear screen and move to top
        
        # self.update_iterator()
        
        lines_used = 0
        header, n = self.get_header_lines(second = False, suppress = True)
        lines_used += n + 1
        all_lines, n = self._get_display_lines(suppress = True)
        lines_used += n
        
        dlines = []
        if self._ref_cache1:
            dlines = self.get_data_lines(self._ref_cache1, self._ref_cache1, bit_depth = 16, threshs = [0.45])
            dlines2 = self.get_data_lines(self._ref_cache1, self._ref_cache1, bit_depth = 8, threshs = [0.45], scale = 8)
            dlines.append("\n")
            dlines.extend(dlines2)
            lines_used += len(dlines)
        # Navigation hints at fixed position (line 24)
        footer = self.get_footer_lines(lines_used, footer_pos = 24, suppress= True)
        
        print("\n".join(header))
        print("-" * 80)
        
        seqlines1 = all_lines.get("seq",[])
        ruler = seqlines1.get("ruler",[])
        seq1 = seqlines1.get("seq",[])
        ext1 = seqlines1.get("extras",[])
        
        print('\n'.join(ruler))
        print('\n'.join(seq1))
        print('\n'.join(ext1))
        
        print("\n".join(dlines))
        
        print("\n".join(all_lines.get("features", [])))
        
        print("\n".join(footer))
    
    def _display_current_view_sec(self):
        
        if not self.debug:
            print('\033[2J\033[H')  # Clear screen and move to top
        
        lines_used = 0
        header, n1 = self.get_header_lines(suppress = True)
        header2, n2 = self.get_header_lines(second=True, suppress = True)
        lines_used += n1+n2+1
        
        all_lines, n = self._get_display_lines(suppress = True)
        lines_used += n
        
        all_lines2, n = self._get_display_lines(second=True, suppress=True)
        lines_used += n
        
        dlines = []
        if self._ref_cache1 and self._ref_cache2:
            dlines = self.get_data_lines(self._ref_cache1, self._ref_cache2, threshs = [0.45])
            lines_used += len(dlines)
            
            dlines = self.get_data_lines(self._ref_cache1, self._ref_cache2, threshs = [0.45], scale = 8)
            lines_used += len(dlines)
        # Navigation hints at fixed position (line 24)
        footer = self.get_footer_lines(lines_used, footer_pos = 32, suppress= True)
        
        print("\n".join(header))
        print("\n".join(header2))
        print("-" * 80)
        
        seqlines1 = all_lines.get("seq",[])
        seqlines2 = all_lines2.get("seq",[])
        ruler = seqlines1.get("ruler",[])
        seq1 = seqlines1.get("seq",[])
        seq2 = seqlines2.get("seq",[])
        ext1 = seqlines1.get("extras",[])
        ext2 = seqlines2.get("extras",[])
        
        print('\n'.join(ruler))
        print('\n'.join(seq1))
        print('\n'.join(ext1))
        print('\n'.join(seq2))
        print('\n'.join(ext2))
        
        print("\n".join(dlines))
        
        if all_lines.get("features"):
            print("\n".join(all_lines.get("features",[])))
            print("\n".join(all_lines2.get("features",[])))
        
        print("\n".join(footer))
        
    def get_data_lines(self, seq1:str, seq2, show_stats = True, threshs = None, bit_depth = 32, scale = None):
        
        vars = [i for i, s in enumerate(seq1) if s=='-']
        seq1_nv = seq1.replace("-","")
        seq2_nv = seq1.replace("-","")
        
        bg, fg = draw.get_color_scheme("test")
        rcfg = fg - 12
        if scale is None:
            scale = len(seq1)//2
            fill = 0.25
            maxval = None
        else:
            fill = 0.25
            maxval = None
            # maxval = 2*fill
                
        conv_res, rc_conv_res = seqs.convolve(seq1_nv, seq2_nv, fill = fill, scale = scale)
        
        data1 = draw.scalar_to_text_nb(conv_res, minval = 0, maxval =maxval, fg_color = fg, bg_color = bg, bit_depth = bit_depth)
        data2 = draw.scalar_to_text_nb(rc_conv_res, minval = 0, maxval = maxval, fg_color = rcfg, bg_color = bg, flip = True, bit_depth = bit_depth)
        
        dn1 = "Corr"
        dn2 = "RCCorr"
        pre2 = f"{Colors.BOLD}{dn1}:    {Colors.RESET}"
        rcpre1 = f"{Colors.BOLD}{dn2}:  {Colors.RESET}"
        blank = "         "
        prefix1 = [blank] * (len(data1)-1) + [pre2]
        prefix2 = [rcpre1] + [blank] * (len(data2)-1) 
        
        outlines = []
        for p, d in zip(prefix1, data1):
            dv = list(d)
            for iv in vars[::-1]:
                dv.insert(iv, '-')
            dvstr = "".join(dv)
            outlines.append(p+dvstr)
        for p, d in zip(prefix2, data2):
            dv = list(d)
            for iv in vars[::-1]:
                dv.insert(iv, '-')
            dvstr = "".join(dv)
            outlines.append(p+dvstr)
        
        if show_stats:
            outlines.append("{:<10}{:<10}{:<10}{:<10}{:<10}".format(*["","Mean", "SD", "Min", "Max"]))
            outlines.append(f"{dn1:<10}{np.mean(conv_res):<10.2g}{np.std(conv_res):<10.2g}{min(conv_res):<10.2g}{max(conv_res):<10.2g}")
            outlines.append(f"{dn2:<10}{np.mean(rc_conv_res):<10.2g}{np.std(rc_conv_res):<10.2g}{min(rc_conv_res):<10.2g}{max(rc_conv_res):<10.2g}")
        
        if threshs:
            tlines = []
            
            tres1 = self.get_threshold_lines(conv_res, dn1, seq1, threshs)
            tlines.extend(tres1)
            tres2 = self.get_threshold_lines(rc_conv_res, dn2, seq2, threshs)
            tlines.extend(tres2)
            
            if len(tlines) > 0:
                outlines.append(f"{Colors.BOLD}Threshs: {Colors.RESET}")
                outlines.extend(tlines)
            
        return outlines
    
    def get_threshold_lines(self, data, data_name, seq, threshs):
        tlines = []
        t = threshs[0]
        tpos = np.argmax(data)
        tval = data[tpos]
        if abs(tval) > t:
            minidx = max(0, tpos - 8)
            maxidx = max(0, tpos + 8 + 1)
            tseq = seq[minidx:maxidx]
            tlines.append(f"{data_name} met threshold {t:0.2g} with {tval:0.2f} at {tpos}, seq {tseq}")
        return tlines
    
    def get_header_lines(self, second = False, suppress =False):
        
        if second:
            state = self.state2
        else:
            state = self.state
        
        n =  0
        hlines = []
        # Header with mode indicators
        mode_indicators = []
        if state.show_reverse_strand:
            mode_indicators.append("(-) strand")
        else:
            mode_indicators.append("(+) strand")
        if state.show_rna:
            mode_indicators.append("RNA")
        else:
            mode_indicators.append("DNA")
        n += 2
        
        if not suppress:
            print(f"{Colors.BOLD}Chromosome {state.chrom} | "
                f"Position {state.position:,} | "
                f"Window {state.window_size}bp | "
                f"{' | '.join(mode_indicators)}{Colors.RESET}")
            # print("-" * 80)
        hlines.append(f"{Colors.BOLD}Chromosome {state.chrom} | "
              f"Position {state.position:,} | "
              f"Window {state.window_size}bp | "
              f"{' | '.join(mode_indicators)}{Colors.RESET}")
        n += 1
        return hlines, n
    
    def get_footer_lines(self, lines_used, footer_pos=24, suppress = False):
        
        lines = []
        remaining_lines = footer_pos - lines_used
        if remaining_lines > 0:
            if not suppress:
                print("\n" * remaining_lines)
            lines.append("\n" * remaining_lines)
        
        if not suppress:
            print("-" * 80)
            print(f"{Colors.DIM}[←/→: move | Ctrl+←/→: gene | ↑/↓: jump | g: goto | Ctrl+↑: save | Ctrl+↓: load | ?: help | q: quit]{Colors.RESET}")
        lines.append("-" * 80)
        lines.append(f"{Colors.DIM}[←/→: move | Ctrl+←/→: gene | ↑/↓: jump | g: goto | Ctrl+↑: save | Ctrl+↓: load | ?: help | q: quit]{Colors.RESET}")
        return lines
    
    def _display_current_view(self):
        """Display the current genomic view."""
        # Clear previous display and ensure consistent height
        
        if not self.debug:
            print('\033[2J\033[H')  # Clear screen and move to top
        
        # Track line count to maintain 24-line height
        lines_used = 0
        
        # Header with mode indicators
        mode_indicators = []
        if self.state.show_reverse_strand:
            mode_indicators.append("(-) strand")
        else:
            mode_indicators.append("(+) strand")
        if self.state.show_rna:
            mode_indicators.append("RNA")
        else:
            mode_indicators.append("DNA")
        
        print(f"{Colors.BOLD}Chromosome {self.state.chrom} | "
              f"Position {self.state.position:,} | "
              f"Window {self.state.window_size}bp | "
              f"{' | '.join(mode_indicators)}{Colors.RESET}")
        print("-" * 80)
        
        try:
            window = self.iterator.get_window_at(self.state.position)
            ref_seq, personal_seq, features = window.ref_seq, window.alt_seq, window.features
            variant_features = window.variant_features if window.variant_features else []
            variant_deltas = window.variant_deltas  # For coordinate tracking
            
            # Add detected motifs to features list
            if window.motifs:
                logger.debug(f"Adding {len(window.motifs)} detected motifs to features")
                # Convert motif dicts to feature-like objects
                for motif in window.motifs:
                    features.append({
                        'feature': 'motif',
                        'feature_type': 'motif',
                        'start': motif['start'],
                        'end': motif['end'],
                        'strand': '+',
                        'info': {
                            'name': motif['name'],
                            'score': motif.get('score', 0),
                            'sequence': motif.get('sequence', '')
                        },
                        'source': 'motif_detector'
                    })
            
            # Get additional annotations from unified system
            if hasattr(self.gm, 'annotations'):
                unified_features = self.gm.get_all_annotations(
                    str(self.state.chrom),
                    self.state.position,
                    self.state.position + self.state.window_size - 1,
                    include_motifs=True
                )
                features = unified_features
            
            self._cache_current_gene(features)
            
            # Apply strand and RNA transformations
            if self.state.show_reverse_strand:
                ref_seq = reverse_complement(ref_seq, rna=False)
                personal_seq = reverse_complement(personal_seq, rna=False)
                # Reverse features positions for display
                # (features stay the same, just displayed differently)
            
            if self.state.show_rna:
                ref_seq = to_rna(ref_seq)
                personal_seq = to_rna(personal_seq)
            
            # Display sequences with codon highlighting
            self._display_sequences(ref_seq, personal_seq, features)
            
            # scalar data
            data_results = self.calculate_data(self.data, personal_seq)
            self._display_data(self.data, data_results, window.variant_features)
            
            # Display amino acid translation if in CDS
            if self.state.show_amino_acids:
                self._display_amino_acid_changes(ref_seq, personal_seq, features)
            
            # Display features
            if self.state.show_features and features:
                self._display_features(features, variant_features)
            
            # Count actual lines used so far
            lines_used = 15  # Approximate count from header, sequences, features, etc.
            
            # Navigation hints at fixed position (line 24)
            remaining_lines = 24 - lines_used
            if remaining_lines > 0:
                print("\n" * remaining_lines, end="")
            
            print("-" * 80)
            print(f"{Colors.DIM}[←/→: move | Ctrl+←/→: gene | ↑/↓: jump | g: goto | Ctrl+↑: save | Ctrl+↓: load | ?: help | q: quit]{Colors.RESET}")
            
        except StopIteration:
            print(f"{Colors.DIM}No sequence data available at this position{Colors.RESET}")
    
    def get_all_features(self, window, state,  features =[]):
        
        if window.motifs:
            for motif in window.motifs:
                motif.pop("sequence")
                motif["feature_type"] = motif.pop("type")
                motif["source"] = "motif"
                motif_feat = UnifiedFeature(chrom = state.chrom, **motif)
                features.append(motif_feat)
                
        if hasattr(self.gm, 'annotations'):
            unified_features = self.gm.get_all_annotations(
                str(state.chrom),
                state.position,
                state.position + state.window_size - 1,
                include_motifs=True
            )
            features.extend(unified_features)
        return features
    
    def _get_display_lines(self, second = False, suppress = False):
        
        nused = 0
        all_lines = {}
        
        if second:
            state = self.state2
            iterator = self.iterator2
        else:
            state = self.state
            iterator = self.iterator
        
        window = iterator.get_window_at(state.position)
        ref_seq, personal_seq, features = window.ref_seq, window.alt_seq, window.features
        
        logger.debug(f"sequence from window with second {second}, pos {state.position}")
        logger.debug(f"iterator data: idx {iterator.chrom}:{iterator.start}-{iterator.end}, windowsz {iterator.window_size}, stride {iterator.stride}")
        
        variant_features = window.variant_features if window.variant_features else []
        variant_deltas = window.variant_deltas  # For coordinate tracking
        
        features = self.get_all_features(window, state, features)
        
        self._cache_current_gene(features, second=second)
        
        if ref_seq:
            logger.debug(f"ref seq {ref_seq[:32]}")
        else:
            logger.debug("no ref seq?")
        if personal_seq:
            logger.debug(f"alt seq {personal_seq[:32]}")
        else:
            logger.debug("no alt seq?")
        
        if state.show_rna:
            ref_seq = to_rna(ref_seq)
            personal_seq = to_rna(personal_seq)
        if state.show_reverse_strand:
            ref_seq = reverse_complement(ref_seq, rna=False)
            personal_seq = reverse_complement(personal_seq, rna=False)
        
        if second:
            self._ref_cache2 = ref_seq
            self._alt_cache2 = personal_seq
        else:
            self._ref_cache1 = ref_seq
            self._alt_cache1 = personal_seq
        
        # Display sequences with codon highlighting
        slines, n = self._display_sequences(ref_seq, personal_seq, features, second = second, suppress = suppress)
        if n < 1:
            raise Exception("no ref or personal seq to display")
        
        all_lines["seq"] = slines
        nused += n
        
        if not state.show_second:
            if state.show_amino_acids:
                alines = self._display_amino_acid_changes(ref_seq, personal_seq, features, second = second, suppress = suppress)
                all_lines["aa"]=alines
                nused += len(alines)
        
        if state.show_features and features:
            flines = self._display_features(features, variant_features, second = second, suppress = suppress)
            all_lines["features"] = flines
            nused += len(flines)
        return all_lines, nused
        
    def _cache_current_features(self, features):
        self._current_features={}
        for f in features:
            ftype = f['feature']
            if not ftype in self._current_features:
                self._current_features=[]
            self._current_features.append(f)
    
    def _cache_current_gene(self, features, second = False):
        if features:
            for f in features:
                if f.feature_type=="gene":
                    if 'gene_name' in f.attributes:
                        if second:
                            self._gene_cache2=f.get('gene_name')
                        else:
                            self._gene_cache=f.get('gene_name')
                        return
                    elif 'info' in f.attributes:
                        if second:
                            self._gene_cache2 = f.attributes.get('info').get('gene_name',"?")
                        else:
                            self._gene_cache = f.attributes.get('info').get('gene_name',"?")
                        return 
                
        if second:
            self._gene_cache2 = ""
        else:
            self._gene_cache=""
        
    
    def _open_gene_card(self):
        
        if not self._gene_cache:
            return
        
        subprocess.run(["firefox",self.genecard.format(name=self._gene_cache)])
    
    def _open_wiki(self):
        
        if not self._gene_cache:
            return
        
        subprocess.run(["firefox",self.wiki.format(name=self._gene_cache)])
    
    def calculate_data(self, data, seq):
        
        data_res = []
        if "autocorrelation" in data:
            res, rcres = seqs.convolve(seq, seq, fill = 0.25)
            data_res.append([res, rcres])
        
        return data_res
    
    def _display_data(self, data, data_results, variants, bit_depth = 8):
        bc, fc = draw.get_color_scheme("test")
        
        for d, res in zip(data, data_results):
            
            if not isinstance(res, list):
                res = [res]
            
                
            data_disp = f"{Colors.BOLD}{d[:7]}: {Colors.RESET}"
            for i, eachres in enumerate(res):
                
                data_render = draw.scalar_to_text_16b(eachres, fg_color = fc, bg_color = bc, flip = i%2)
                for r in data_render:
                    data_render_algn = self.iterator.annotations.sequence_stream.apply_variants_to_sequence(r, variants, self.state.position, substitution = " ")
                    data_disp += "".join(data_render_algn)
                    print(data_disp)
                    data_disp = "         "
                print()
            
    
    def _display_sequences(self, ref_seq: str, personal_seq: str, features= None, second = False, suppress = False):
        """Display reference and personal sequences with variant highlighting.
        
        Note: Sequences are already aligned with gaps from GenomeIterator.
        
        Args:
            ref_seq: Aligned reference sequence (may contain gaps '-')
            personal_seq: Aligned personal sequence (may contain gaps '-')
            features: Features in current window for codon highlighting
        """
        if second:
            state = self.state2
        else:
            state = self.state
        
        currstrand = '-' if state.show_reverse_strand else '+'
        # Extract motifs from features for underline display
        motifs = [f for f in features if f.get('feature_type') == 'motif' or f.get('feature') == 'motif'] if features else []
        
        lines = {}
        n = 0
        
        if motifs:
            logger.debug(f"displaying {len(motifs)} detected motifs")
            logger.debug(f"motifs: {motifs}")
        if not ref_seq or not personal_seq:
            return {}, 0
        
        # Sequences are already aligned from GenomeIterator
        aligned_ref = ref_seq
        aligned_pers = personal_seq
        
        # data_results = []
        # if self.data:
        #     data_results = self.calculate_data(self.data, aligned_pers)
        
        # Detect variant positions (including gaps)
        variant_positions = []
        for i in range(min(len(aligned_ref), len(aligned_pers))):
            # A position is a variant if bases differ or if there's a gap
            is_variant = (aligned_ref[i] != aligned_pers[i] or 
                         aligned_ref[i] == '-' or aligned_pers[i] == '-')
            variant_positions.append(is_variant)
        
        # Find start and stop codons in features
        start_codons = []
        stop_codons = []
        disp_motifs = []
        if features:
            window_start = state.position
            for feature in features:
                ftype = feature.feature_type if isinstance(feature, UnifiedFeature) else feature.get("feature_type")
                if ftype == 'start_codon':
                    # Calculate raw positions relative to window
                    raw_start = max(0, feature.start - window_start)
                    raw_end = min(len(aligned_ref), feature.end - window_start + 1)
                    
                    # Apply strand-aware rendering
                    if state.show_reverse_strand:
                        # Reflect positions across window center
                        display_start = state.window_size - raw_end
                        display_end = state.window_size - raw_start
                    else:
                        display_start = raw_start
                        display_end = raw_end
                    
                    # Ensure positions are within bounds
                    display_start = max(0, min(display_start, len(aligned_ref)))
                    display_end = max(0, min(display_end, len(aligned_ref)))
                    
                    if display_start < display_end:
                        start_codons.append((display_start, display_end))
                        
                elif ftype == 'stop_codon':
                    # Calculate raw positions relative to window
                    raw_start = max(0, feature.start - window_start)
                    raw_end = min(len(aligned_ref), feature.end - window_start + 1)
                    
                    # Apply strand-aware rendering
                    if state.show_reverse_strand:
                        # Reflect positions across window center
                        display_start = state.window_size - raw_end
                        display_end = state.window_size - raw_start
                    else:
                        display_start = raw_start
                        display_end = raw_end
                    
                    # Ensure positions are within bounds
                    display_start = max(0, min(display_start, len(aligned_ref)))
                    display_end = max(0, min(display_end, len(aligned_ref)))
                    
                    if display_start < display_end:
                        stop_codons.append((display_start, display_end))
                elif ftype == 'motif':
                    # Calculate raw positions relative to window
                    raw_start = max(0, feature.start - window_start)
                    raw_end = min(len(aligned_ref), feature.end - window_start + 1)
                    
                    # Apply strand-aware rendering
                    if state.show_reverse_strand:
                        # Reflect positions across window center
                        display_start = state.window_size - raw_end
                        display_end = state.window_size - raw_start
                    else:
                        display_start = raw_start
                        display_end = raw_end
                    
                    # Ensure positions are within bounds
                    display_start = max(0, min(display_start, len(aligned_ref)))
                    display_end = max(0, min(display_end, len(aligned_ref)))
                    
                    if display_start < display_end:
                        disp_motifs.append((display_start, display_end))
        
        # Calculate position labels
        start_pos = state.position
        pos_labels = []
        for i in range(0, len(ref_seq), 10):
            pos_labels.append(f"{start_pos + i:,}")
        
        rlines = []
        # Position ruler (adjusted for alignment)
        if not suppress:
            print(f"\n{Colors.POSITION}Position:{Colors.RESET}")
        rlines.append(f"\n{Colors.POSITION}Position:{Colors.RESET}")
        strand = '-' if state.show_reverse_strand else '+'
        ruler = ""
        pos_counter = 0  # Track actual genomic positions despite gaps
        for i in range(len(aligned_ref)):
            if aligned_ref[i] != '-':  # Only count non-gap positions
                if pos_counter % 10 == 0:
                    ruler += "|"
                elif pos_counter % 5 == 0:
                    ruler += strand
                else:
                    ruler += "."
                pos_counter += 1
            else:
                ruler += " "  # Space for gaps
        if not suppress:
            print(f"         {ruler}")
        rlines.append(f"         {ruler}")
        
        # Position numbers (every 10 bases, adjusted for gaps)
        pos_line = "         "
        pos_counter = 0
        for i in range(len(aligned_ref)):
            if aligned_ref[i] != '-':
                if pos_counter % 10 == 0:
                    label = f"{start_pos + pos_counter:,}"
                    pos_line = pos_line[:len(pos_line)] + label
                pos_counter += 1
            pos_line += " " if i < len(aligned_ref) - 1 else ""
        if not suppress:
            print(pos_line[:len(aligned_ref) + 9])
        rlines.append(pos_line[:len(aligned_ref) + 9])
        lines["ruler"] = rlines
        n += len(rlines)
        
        slines = []
        # Reference sequence with codon highlighting
        ref_display = f"\n{Colors.BOLD}Ref:    {Colors.RESET} "
        for i, base in enumerate(aligned_ref):
            # Check if in start/stop codon
            in_start = any(start <= i < end for start, end in start_codons)
            in_stop = any(start <= i < end for start, end in stop_codons)
            in_motif = any(start <= i < end for start, end in disp_motifs)
            
            if base == '-':
                ref_display += f"{Colors.DIM}-{Colors.RESET}"
            elif in_start:
                ref_display += f"{Colors.INSERTION}{base}{Colors.RESET}"
            elif in_stop:
                ref_display += f"{Colors.SNP}{base}{Colors.RESET}"
            elif in_motif:
                ref_display += f"{Colors.EXON}{base}{Colors.RESET}"
            else:
                ref_display += base
        if not suppress:
            print(ref_display)
        slines.append(ref_display)
        
        # Personal sequence with variant and codon coloring
        personal_display = f"{Colors.BOLD}Persnl: {Colors.RESET} "
        
        for i in range(len(aligned_pers)):
            ref_base = aligned_ref[i] if i < len(aligned_ref) else ''
            pers_base = aligned_pers[i]
            
            # Check if in start/stop codon
            in_start = any(start <= i < end for start, end in start_codons)
            in_stop = any(start <= i < end for start, end in stop_codons)
            in_motif = any(start <= i < end for start, end in disp_motifs)
            
            if pers_base == '-':
                personal_display += f"{Colors.DIM}-{Colors.RESET}"
            elif variant_positions[i] if i < len(variant_positions) else False:
                # Variant detected
                color = Colors.variant_color(ref_base, pers_base)
                personal_display += f"{color}{pers_base}{Colors.RESET}"
            elif in_start:
                personal_display += f"{Colors.INSERTION}{pers_base}{Colors.RESET}"
            elif in_stop:
                personal_display += f"{Colors.SNP}{pers_base}{Colors.RESET}"
            elif in_motif:
                personal_display += f"{Colors.EXON}{pers_base}{Colors.RESET}"
            else:
                personal_display += pers_base
        
        if not suppress:
            print(personal_display)
        slines.append(personal_display)
        lines["seq"] = slines
        n+= len(slines)
        

        
        marker_line = ["         "]
        offset = len(marker_line)
        la, ra = self.arrows
        stem = '-'
        extras  = []
        # Motif underline with strand indicators
        if motifs or True:
            
            # marker_line = "         "
            motif_colors = "         "
            
            for i in range(len(aligned_ref)):
                # Check if position is in any motif
                motif_found = None
                for motif in motifs:
                    window_start = state.position
                    motif_start = max(0, motif['start'] - window_start)
                    motif_end = min(len(aligned_ref), motif['end'] - window_start)
                    
                    if motif_start <= i < motif_end:
                        motif_found = motif
                        break
                
                if motif_found:
                    # Choose color based on motif type
                    motif_name = motif_found.get('info', {}).get('name', motif_found.get('name', ''))
                    strand = motif_found.get('strand', '+')
                    
                    if 'splice' in motif_name.lower():
                        color = Colors.MOTIF_SPLICE
                    elif 'tata' in motif_name.lower():
                        color = Colors.MOTIF_TATA
                    elif 'cpg' in motif_name.lower():
                        color = Colors.MOTIF_CPG
                    elif 'polya' in motif_name.lower():
                        color = Colors.MOTIF_POLYA
                    else:
                        color = Colors.MOTIF_DEFAULT
                    
                    
                    if i == motif_start:
                        marker_line.append(f"{color}{ra if strand == '+' else la}{Colors.RESET}")
                    else:
                        marker_line.append(f"{color}{stem}{Colors.RESET}")
                    
                else:
                    marker_line.append(" ")
            
            if not suppress:
                print(marker_line)
            extras.append(marker_line)
        
        # Variant indicator line
        for i in range(offset, len(variant_positions)):
            if variant_positions[i]:
                marker_line[i] = f"{Colors.DIM}*{Colors.RESET}"
            # else:
            #     marker_line += " "
        
        # if any(variant_positions):
        #     if not suppress:
        #         print(f"{Colors.DIM}{variant_line}{Colors.RESET}")
        #     extras.append(f"{Colors.DIM}{variant_line}{Colors.RESET}")
        # lines["extras"] = extras
        # n+=len(extras)
        
        return lines, n
    
    
    
    def _display_features(self, features: List[Dict[str, Any]], variant_features, second=False, suppress = False):
        """Display feature annotations as labeled bars with compression.
        
        Args:
            features: List of features in the current window
            variant_features: List of variant feature objects
        """
        if not features:
            logger.debug('no features?')
            return []
        
        if second:
            state = self.state2
        else:
            state = self.state
        
        feature_types = ["gene","exon","CDS","motif"]
        
        lines = []
        
        # Separate motifs from other features
        motifs = []
        other_features = []
        for feature in features:
            ftype = feature.feature_type if hasattr(feature, 'feature_type') else feature.get('feature_type', feature.get('feature'))
            if not ftype in feature_types:
                continue
            if ftype == 'motif':
                motifs.append(feature)
            else:
                other_features.append(feature)
        
        if not suppress:
            print(f"\n{Colors.BOLD}Features:{Colors.RESET}")
        lines.append(f"\n{Colors.BOLD}Features:{Colors.RESET}")
        
        # Group features by type and extent
        feature_groups = {}
        for feature in other_features:
            ftype = feature.feature_type if hasattr(feature, 'feature_type') else feature.get('feature_type', feature.get('feature'))
            # Skip start/stop codons as they're shown differently
            if ftype in ['start_codon', 'stop_codon']:
                continue
            if ftype not in feature_groups:
                feature_groups[ftype] = {}
            
            # Group by extent (start and end positions)
            feat_start = feature.start if hasattr(feature, 'start') else feature.get('start')
            feat_end = feature.end if hasattr(feature, 'end') else feature.get('end')
            extent_key = (feat_start, feat_end)
            
            if extent_key not in feature_groups[ftype]:
                feature_groups[ftype][extent_key] = []
            feature_groups[ftype][extent_key].append(feature)
        
        # Display features in hierarchical order
        window_start = state.position
        
        currstrand = '-' if state.show_reverse_strand else '+'
        for ftype in self.feature_hierarchy:
            if ftype not in feature_groups:
                continue
            
            color = self._get_feature_color(ftype)
            
            # Display each unique extent
            for extent_key, extent_features in feature_groups[ftype].items():
                feat_start, feat_end = extent_key
                
                # Calculate display positions using render_pos for strand awareness
                # First get the feature position relative to window
                rel_start = feat_start - window_start
                rel_end = feat_end - window_start
                
                # Apply strand-aware rendering
                if state.show_reverse_strand:
                    # Reverse the positions for display
                    display_start, display_end = self._render_pos(rel_start, rel_end + 1)
                else:
                    display_start = max(0, rel_start)
                    display_end = min(state.window_size, rel_end + 1)
                
                # Only display if feature overlaps with window
                # Check if the feature overlaps with the visible window at all
                if feat_end >= window_start and feat_start <= window_start + state.window_size - 1:
                    # Create feature bar with type label
                    type_label = f"{ftype[:7]:>7}"
                    bar = f"{type_label}: "  # 9 chars total
                    
                    for i in range(state.window_size):
                        # Check if this position is within the feature
                        actual_pos = window_start + i
                        
                        # For reverse strand, we need to check positions differently
                        if state.show_reverse_strand:
                            # On reverse strand, features appear flipped
                            render_i = state.window_size - 1 - i
                            if feat_start <= window_start + render_i <= feat_end:
                                # Determine what character to use
                                if render_i == rel_start and feat_start >= window_start:
                                    bar += "|"  # Feature starts in window
                                elif render_i == rel_end and feat_end <= window_start + state.window_size - 1:
                                    bar += "|"  # Feature ends in window
                                else:
                                    bar += "-"
                            else:
                                bar += " "
                        else:
                            # Normal forward strand display
                            if feat_start <= actual_pos <= feat_end:
                                # Determine what character to use
                                if i == display_start and feat_start >= window_start:
                                    bar += "|"  # Feature starts in window
                                elif i == display_end - 1 and feat_end <= window_start + state.window_size - 1:
                                    bar += "|"  # Feature ends in window
                                else:
                                    bar += "-"
                            else:
                                bar += " "
                    
                    # Create compressed label for multiple features
                    label = self._format_compressed_label(ftype, extent_features)
                    
                    if ftype=='gene':
                        feat = feature_groups[ftype][extent_key][0]
                        if not feat.strand==currstrand:
                            labelstr = label + f'({feat.strand})'
                            label=labelstr
                            
                    
                    # Calculate label position (keep label visible and centered)
                    if state.show_reverse_strand:
                        # For reverse strand, adjust label position
                        visible_start = max(0, min(display_start, display_end))
                        visible_end = min(state.window_size, max(display_start, display_end))
                    else:
                        visible_start = max(0, display_start)
                        visible_end = min(state.window_size, display_end)
                    label_center = visible_start + (visible_end - visible_start) // 2
                    label_pos = 9 + label_center - len(label) // 2
                    
                    # Ensure label stays within the visible part of the feature
                    label_pos = max(9 + visible_start, label_pos)
                    label_pos = min(9 + visible_end - len(label), label_pos)
                    
                    # Insert label into bar if it fits
                    if label_pos > 9 and label_pos + len(label) < len(bar):
                        bar = bar[:label_pos] + label + bar[label_pos + len(label):]
                    if not suppress:
                        print(f"{color}{bar}{Colors.RESET}")
                    lines.append(f"{color}{bar}{Colors.RESET}")
                    
        # lastly display variants separately with details
        # Combine variant features from both sources
        all_variant_features = []
        
        # Add variants from features list (UnifiedFeature objects)
        all_variant_features.extend([f for f in features if f.feature_type == 'variant'])
        
        # Add variants from variant_features list if provided
        if variant_features:
            all_variant_features.extend(variant_features)
        
        # Remove duplicates based on position
        seen_positions = set()
        unique_variants = []
        for var in all_variant_features:
            if not var.feature_type in feature_types:
                continue
            if var.start not in seen_positions:
                unique_variants.append(var)
                seen_positions.add(var.start)
        
        if unique_variants:
            if not suppress:
                print(f"  {Colors.SNP}Variants:{Colors.RESET}")
            lines.append(f"  {Colors.SNP}Variants:{Colors.RESET}")
            for var in unique_variants:
                var_start = var.start
                var_end = var.end
                
                # Get ref and alt from attributes
                attrs = var.attributes if hasattr(var, 'attributes') else var.get('attributes', {})
                ref = attrs.get('ref', '')
                alt = attrs.get('alt', [''])[0] if isinstance(attrs.get('alt'), list) else attrs.get('alt', '')
                qual = attrs.get('qual', var.score if hasattr(var, 'score') else 0)
                gt = attrs.get("genotype","")
                zyg = attrs.get("zygosity")
                
                ref = self._render_seq(ref)
                alt = self._render_seq(alt)
                
                # Only show if overlaps with window
                if var_end >= state.position and var_start <= state.position + state.window_size - 1:
                    var_type = 'SNP' if len(ref) == len(alt) else 'INDEL'
                    print(f"    {var_start}: {ref}→{gt} ({var_type}, Q={qual:.1f}, {zyg})")
        
        # Display motifs section
        if motifs:
            if not suppress:
                print(f"\n{Colors.BOLD}Motifs:{Colors.RESET}")
            lines.append(f"\n{Colors.BOLD}Motifs:{Colors.RESET}")
            
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
                if 'splice' in motif_type.lower():
                    color = Colors.MOTIF_SPLICE
                elif 'tata' in motif_type.lower():
                    color = Colors.MOTIF_TATA
                elif 'cpg' in motif_type.lower():
                    color = Colors.MOTIF_CPG
                elif 'polya' in motif_type.lower():
                    color = Colors.MOTIF_POLYA
                else:
                    color = Colors.MOTIF_DEFAULT
                
                if not suppress:
                    print(f"  {color}{motif_type}:{Colors.RESET}")
                lines.append(f"  {color}{motif_type}:{Colors.RESET}")
                
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
                        print(f"    {start:,} {strand_symbol} {seq} (score: {score:.2f})")
                
                if len(motif_list) > 5:
                    if not suppress:
                        print(f"    ... and {len(motif_list) - 5} more")
                    lines.append(f"    ... and {len(motif_list) - 5} more")
        return lines
    
    def _format_compressed_label(self, ftype: str, features) -> str:
        """Format a compressed label for multiple features with advanced compression.
        
        Args:
            ftype: Feature type
            features: List of features with same extent
            
        Returns:
            Formatted compressed label string
        """
        if len(features) <= 5:
            return self._format_single_feature_label(features[0])
        
        # Check if all features span the full window - compress heavily if so
        window_start = self.state.position
        window_end = window_start + self.state.window_size - 1
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
    
    def _format_single_feature_label(self, feature) -> str:
        """Format a label for a single feature.
        
        Args:
            feature: Feature dictionary
            
        Returns:
            Formatted label string
        """
        ftype = feature.feature_type
        info = feature.get('attributes', {})
        
        if ftype == 'gene':
            # return info.get('gene_name', 'gene')
            return feature.name if feature.name else "?"
        elif ftype == 'transcript':
            name = info.get('transcript_name', info.get('transcript_id', 'transcript'))
            # name=feature.name
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
            info = feature.attributes
            ref = info.get('ref', '?')
            alt = info.get('alt', '?')[0]
            ref = self._render_seq(ref)
            alt = self._render_seq(alt)
            
            
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
    
    def _display_amino_acid_changes(self, ref_seq: str, personal_seq: str, features, second = False):
        """Display amino acid changes for variants in CDS regions.
        
        Note: Sequences may contain alignment gaps ('-') which are removed before translation.
        
        Args:
            ref_seq: Reference DNA sequence (may contain gaps from alignment)
            personal_seq: Personal DNA sequence (may contain gaps from alignment)
            features: Features in current window
        """
        if second:
            state = self.state2
        else:
            state = self.state
        
        # Check if we're in a CDS region
        cds_features = [f for f in features if f.feature_type == 'CDS']
        if not cds_features:
            return
        
        window_start = self.state.position
        window_end = window_start + len(ref_seq) - 1
        
        # For each CDS, determine reading frame and translate
        for cds in cds_features:
            cds_start = cds.start
            cds_end = cds.end
            strand = cds.strand
            
            # Determine the overlap between CDS and window
            overlap_start = max(cds_start, window_start)
            overlap_end = min(cds_end, window_end)
            
            if overlap_start > overlap_end:
                continue  # No overlap
            
            # Calculate positions in the window
            window_cds_start = overlap_start - window_start
            window_cds_end = overlap_end - window_start + 1
            
            # Calculate frame offset based on position in CDS
            cds_offset = overlap_start - cds_start
            frame_offset = cds_offset % 3
            
            # Adjust to codon boundaries within the CDS region only
            adj_start = window_cds_start + (3 - frame_offset) % 3
            adj_end = window_cds_end - ((window_cds_end - window_cds_start - frame_offset) % 3)
            
            if adj_end <= adj_start:
                continue
            
            # Extract only the CDS portion and remove gaps for translation
            ref_coding_with_gaps = ref_seq[adj_start:adj_end]
            pers_coding_with_gaps = personal_seq[adj_start:adj_end]
            
            # Remove gaps for translation (gaps are from alignment display)
            ref_coding = ref_coding_with_gaps.replace('-', '')
            pers_coding = pers_coding_with_gaps.replace('-', '')
            
            # Reverse complement if negative strand
            if strand == '-':
                ref_coding = reverse_complement(ref_coding)
                pers_coding = reverse_complement(pers_coding)
            
            # Translate to amino acids
            ref_aa = self._translate_dna(ref_coding)
            pers_aa = self._translate_dna(pers_coding)
            
            # Display amino acid sequences
            info = cds.get('info', {})
            gene_name = info.get('gene_name', 'Unknown')
            print(f"\n{Colors.BOLD}AA Translation (CDS - {gene_name}):{Colors.RESET}")
            
            # Calculate proper amino acid positioning for current strand
            if state.show_reverse_strand:
                # For reverse strand, amino acids should be positioned from right
                # adj_start and adj_end are already calculated relative to the window
                display_adj_start = state.window_size - adj_end
                display_adj_end = state.window_size - adj_start
            else:
                # For forward strand, use positions as calculated
                display_adj_start = adj_start
                display_adj_end = adj_end
            
            # Ensure display positions are within bounds
            display_adj_start = max(0, min(display_adj_start, state.window_size))
            display_adj_end = max(0, min(display_adj_end, state.window_size))
            
            # Calculate leading spaces to align with sequence display (no 3x multiplication!)
            leading_spaces = display_adj_start
            
            # Position indicators for codons (align with actual sequence positions)
            aa_ruler = " " * 9  # Match sequence labels
            aa_ruler += " " * leading_spaces
            for i in range(len(ref_aa)):
                aa_ruler += f"{i%10}  "
            print(aa_ruler)
            
            # Reference amino acids (properly aligned)
            ref_aa_display = f"Ref AA:  {' ' * leading_spaces}"
            for aa in ref_aa:
                ref_aa_display += f"{aa}  "
            print(f"{Colors.DIM}{ref_aa_display}{Colors.RESET}")
            
            # Personal amino acids with change highlighting
            pers_aa_display = f"Pers AA: {' ' * leading_spaces}"
            for i, (ref_a, pers_a) in enumerate(zip(ref_aa, pers_aa)):
                if ref_a != pers_a:
                    # Missense mutation
                    if pers_a == '*':
                        # Nonsense mutation (stop gained)
                        pers_aa_display += f"{Colors.SNP}{pers_a}{Colors.RESET}  "
                    elif ref_a == '*':
                        # Stop lost
                        pers_aa_display += f"{Colors.INSERTION}{pers_a}{Colors.RESET}  "
                    else:
                        # Regular missense
                        pers_aa_display += f"{Colors.DELETION}{pers_a}{Colors.RESET}  "
                else:
                    pers_aa_display += f"{pers_a}  "
            print(pers_aa_display)
            
            # Show CDS boundaries with bars (strand-aware)
            boundary_line = " " * 9
            
            # Calculate display positions for CDS boundaries
            if state.show_reverse_strand:
                display_cds_start = state.window_size - window_cds_end
                display_cds_end = state.window_size - window_cds_start
            else:
                display_cds_start = window_cds_start
                display_cds_end = window_cds_end
                
            for i in range(state.window_size):
                if i >= display_cds_start and i < display_cds_end:
                    if i == display_cds_start:
                        boundary_line += "|"
                    elif i == display_cds_end - 1:
                        boundary_line += "|"
                    else:
                        boundary_line += "-"
                else:
                    boundary_line += " "
            print(f"{Colors.DIM}{boundary_line} CDS region{Colors.RESET}")
            
            # Show mutation types
            mutations = []
            for i, (ref_a, pers_a) in enumerate(zip(ref_aa, pers_aa)):
                if ref_a != pers_a:
                    # Calculate actual position in the CDS
                    aa_position_in_cds = (overlap_start - cds_start + frame_offset) // 3 + i + 1
                    if pers_a == '*':
                        mutations.append(f"p.{ref_a}{aa_position_in_cds}* (nonsense)")
                    elif ref_a == '*':
                        mutations.append(f"p.*{aa_position_in_cds}{pers_a} (stop lost)")
                    else:
                        mutations.append(f"p.{ref_a}{aa_position_in_cds}{pers_a}")
            
            if mutations:
                print(f"{Colors.DIM}Mutations: {', '.join(mutations[:5])}{Colors.RESET}")
            
            break  # Only show first CDS for clarity
    
    def _translate_dna(self, dna_seq: str) -> str:
        """Translate DNA sequence to amino acids.
        
        Args:
            dna_seq: DNA sequence (should be multiple of 3, gaps removed)
            
        Returns:
            Amino acid sequence
        """
        # Remove any remaining gaps just in case
        dna_seq = dna_seq.replace('-', '')
        
        # Standard genetic code
        codon_table = CODON_TABLE_DNA
        
        amino_acids = []
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3].upper()
            if 'N' in codon or '-' in codon:
                amino_acids.append('X')  # Unknown
            else:
                amino_acids.append(codon_table.get(codon, 'X'))
        
        return ''.join(amino_acids)
    
    def _translate_transcript(self):
        
        feats = self.gm.get_features_at_position(self.state.chrom, self.state.position, include_types=['transcript'])
        i=0
        for f in feats:
            print(i, f['info']['transcript_id'])
            i+=1
        
        txi = input("transcript to transcribe")
        txi = int(txi)
        
        if not txi < i:
            print('invalid transcript index')
            return
        
        txid = feats[txi]['info']['transcript_id']
        seq = self.gm.transcribe_gene(self._gene_cache, txi)
        res = seq['result']
        
        print(f"Warnings: {', '.join(res.warnings)}")
        print(f"Protein sequence for gene {self._gene_cache}, {txid}:")
        for i in range(len(res.protein_sequence)//80 + 1):
            print(res.protein_sequence[i*80:i*80+80]+'...')
    
    def _get_feature_color(self, ftype: str) -> str:
        """Get color for a feature type.
        
        Args:
            ftype: Feature type
            
        Returns:
            Color code
        """
        color_map = {
            'gene': Colors.GENE,
            'transcript': Colors.TRANSCRIPT,
            'exon': Colors.EXON,
            'CDS': Colors.CDS,
            'five_prime_utr': Colors.UTR,
            'three_prime_utr': Colors.UTR,
            'start_codon': Colors.SNP,
            'stop_codon': Colors.SNP,
            'variant': Colors.DELETION,
            'motif': Colors.EXON,  # Green for motifs
            'tf_binding': Colors.INSERTION,  # Green for TF binding
            'splice_donor': Colors.SNP,  # Red for splice sites
            'splice_acceptor': Colors.SNP,
            'repeat': Colors.DIM,  # Gray for repeats
            'chip_peak': Colors.EXON  # Cyan for ChIP peaks
        }
        return color_map.get(ftype, Colors.DIM)
    
    def _move_forward(self, large: bool = False, second = False):
        """Move forward in the genome (considers strand direction)."""
        if second:
            state = self.state2
        else:
            state = self.state
        
        move_distance = state.stride * (10 if large else 1)
        if state.show_reverse_strand:
            # On reverse strand, forward means lower positions
            state.position = max(1, state.position - move_distance)
        else:
            # On forward strand, forward means higher positions
            state.position += move_distance
    
    def _move_backward(self, large: bool = False, second = False):
        """Move backward in the genome (considers strand direction)."""
        if second:
            state = self.state2
        else:
            state = self.state
        move_distance = state.stride * (10 if large else 1)
        if state.show_reverse_strand:
            # On reverse strand, backward means higher positions
            state.position += move_distance
        else:
            # On forward strand, backward means lower positions
            state.position = max(1, state.position - move_distance)
    
    def _switch_chromosome(self, second = False):
                
        if second:
            state = self.state2
        else:
            state = self.state
        c=input("new chromosome: ")
        pos=input("new position (default=1,000,000): ")
        if not pos:
            pos = 1000000
        else:
            pos=int(pos)
        
        if c in list(map(str,range(23))) + ['MT','X','Y']:
            state.chrom=c
            state.position=pos
            
        # self._display_current_view()
        
    def _goto_position(self, second = False):
        """Jump to a specific position."""
        
        if second:
            state = self.state2
        else:
            state = self.state
        print("\n" + "-" * 40)
        try:
            pos_input = input("Enter position (or gene name): ").strip()
            
            # Try to parse as number
            try:
                new_pos = int(pos_input.replace(',', ''))
                state.position = max(1, new_pos)
            except ValueError:
                # Try as gene name
                gene_name = pos_input.upper()
                chrom, gene_info = self.gm.gene_map.find_gene(gene_name)
                
                if chrom and gene_info:
                    gene_data = gene_info[0]
                    state.chrom = chrom
                    state.position = gene_data['start']
                    print(f"Jumping to {gene_name} on chromosome {chrom}")
                else:
                    print(f"Gene '{gene_name}' not found")
                    input("Press Enter to continue...")
        except Exception as e:
            print(f"Error: {e}")
            input("Press Enter to continue...")
    
    def _next_feature(self):
        """Jump to the next feature of specified type."""
        print("\n" + "-" * 40)
        feat_type = input("Next feature type (default: gene): ").strip()
        if not feat_type:
            feat_type = "gene"
        
        # gene_end = self.gm
        
        result = self.gm.find_next_feature(
            self.state.chrom, 
            self.state.position,
            feat_type
        )
        
        if result:
            self.state.position = result['start']
            print(f"Jumped to {result['type']}: {result['name']} at position {result['start']:,}")
            print(f"Distance: {result['distance']:,} bp")
        else:
            print(f"No {feat_type} found downstream")
            input("Press Enter to continue...")
    
    def _prev_feature(self):
        """Jump to the previous feature of specified type."""
        print("\n" + "-" * 40)
        feat_type = input("Previous feature type (default: gene): ").strip()
        if not feat_type:
            feat_type = "gene"
        
        result = self.gm.find_prev_feature(
            self.state.chrom,
            self.state.position,
            feat_type
        )
        
        if result:
            self.state.position = result['start']
            print(f"Jumped to {result['type']}: {result['name']} at position {result['start']:,}")
            print(f"Distance: {result['distance']:,} bp")
        else:
            print(f"No {feat_type} found upstream")
            input("Press Enter to continue...")
    
    def _change_window_size(self):
        """Change the window size."""
        print("\n" + "-" * 40)
        try:
            new_size = int(input(f"Enter new window size (current: {self.state.window_size}): "))
            if 10 <= new_size <= 400:
                self.state.window_size = new_size
                self.state2.window_size = new_size
            else:
                print("Window size must be between 10 and 200")
                input("Press Enter to continue...")
        except ValueError:
            print("Invalid input")
            input("Press Enter to continue...")
    
    def _change_stride(self):
        """Change the stride."""
        print("\n" + "-" * 40)
        try:
            new_stride = int(input(f"Enter new stride (current: {self.state.stride}): "))
            if 1 <= new_stride <= 1000:
                self.state.stride = new_stride
                self.state2.stride = new_stride
            else:
                print("Stride must be between 1 and 1000")
                input("Press Enter to continue...")
        except ValueError:
            print("Invalid input")
            input("Press Enter to continue...")
    
    def _show_position_info(self):
        """Show detailed information about current position."""
        print("\n" + "=" * 60)
        print(f"Position Information")
        print("=" * 60)
        
        # Get all features at current position
        features = self.gm.get_features_at_position(
            self.state.chrom, 
            self.state.position
        )
        
        print(f"Chromosome: {self.state.chrom}")
        print(f"Position: {self.state.position:,}")
        print(f"Features at this position: {len(features)}")
        
        if features:
            print("\nFeature details:")
            for feature in features:
                ftype = feature.feature_type
                info = feature.get('info', {})
                print(f"  - {ftype}:")
                
                # Show relevant info
                if 'gene_name' in info:
                    print(f"      Gene: {info['gene_name']}")
                if 'transcript_name' in info:
                    print(f"      Transcript: {info['transcript_name']}")
                if 'exon_number' in info:
                    print(f"      Exon number: {info['exon_number']}")
                if 'strand' in feature:
                    print(f"      Strand: {feature['strand']}")
        
        # Check for variants at this position
        try:
            variants = []
            for var in self.gm.vcf(self.gm._make_index(
                self.state.chrom, 
                self.state.position, 
                self.state.position
            )):
                variants.append(var)
            
            if variants:
                print(f"\nVariants at this position: {len(variants)}")
                for var in variants:
                    print(f"  - {var.REF} -> {var.ALT[0] if var.ALT else 'N/A'}")
                    print(f"      Quality: {var.QUAL:.1f}")
        except:
            pass
        
        input("\nPress Enter to continue...")
    
    def _show_help(self):
        """Show help information."""
        print("\n" + "=" * 60)
        print("Genome Browser Help")
        print("=" * 60)
        print("Navigation:")
        print("  → / l       : Move forward by stride")
        print("  ← / h       : Move backward by stride")
        print("  ↑ / k       : Jump forward (10x stride)")
        print("  ↓ / j       : Jump backward (10x stride)")
        print("  Ctrl+→      : Jump to next gene (fast)")
        print("  Ctrl+←      : Jump to previous gene (fast)")
        print("  Space/Enter : Move forward by stride")
        print("  g           : Go to position or gene")
        print("  n           : Jump to next feature (interactive)")
        print("  p           : Jump to previous feature (interactive)")
        print("\nState Management:")
        print("  Ctrl+↑      : Save current state")
        print("  Ctrl+↓      : Load saved state")
        print("\nDisplay:")
        print("  f           : Toggle feature display")
        print("  a           : Toggle amino acid translation")
        print("  r           : Toggle RNA mode (T→U)")
        print("  -           : Toggle reverse strand (reflects codons & AAs)")
        print("  w           : Change window size")
        print("  s           : Change stride")
        print("  i           : Show position info")
        print("\nColors:")
        print(f"  {Colors.SNP}Red{Colors.RESET}        : SNPs / Stop codons")
        print(f"  {Colors.INSERTION}Green{Colors.RESET}      : Insertions / Start codons")
        print(f"  {Colors.DELETION}Yellow{Colors.RESET}     : Deletions / Missense mutations")
        print(f"  {Colors.GENE}Blue{Colors.RESET}       : Genes")
        print(f"  {Colors.EXON}Cyan{Colors.RESET}       : Exons")
        print("\nFeature Labels:")
        print("  Left side shows feature types with smart compression")
        print("  Multiple features: ex1/2/3, GENE1/+2, 001-015(25)")
        print("\nOther:")
        print("  ?           : Show this help")
        print("  q           : Quit")
        print("=" * 60)
        input("Press Enter to continue...")
    
    def _get_key(self):
        """Get a single keypress from the user."""
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            key = sys.stdin.read(1)
            
            logger.debug(f"raw key: {repr(key)}")
            # Check for escape sequences (arrow keys and ctrl combinations)
            if key == '\x1b':
                next_chars = sys.stdin.read(2)
                key += next_chars
                
                prefix = suffix = ""
                # Check for Ctrl+Arrow combinations
                if next_chars == '[1':
                    # Read the next character to identify Ctrl+Arrow
                    extra = sys.stdin.read(2)
                    if extra == ';5':
                        prefix = "ctrl_"
                    elif extra == ';2':
                        prefix = "shift_"
                    logger.debug(f"prefix {prefix} from next_chars {repr(next_chars)} and extra {repr(extra)}")
                    
                    final = sys.stdin.read(1)
                    if final == 'C':  # Ctrl+Right
                        suffix = 'right'
                    elif final == 'D':  # Ctrl+Right
                        suffix = 'left'
                    elif final == 'A':  # Ctrl+Right
                        suffix = 'up'
                    elif final == 'B':  # Ctrl+Right
                        suffix = 'down'
                    elif final in "c":
                        suffix = final
                    
                    logger.debug(f"suffix {suffix} from final {repr(final)}")
                    if prefix and suffix:
                        return prefix+suffix
                        
                    key += extra
            
            return key
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    
    def _jump_to_next_gene(self):
        """Quick jump to next gene using Ctrl+Right."""
        result = self.gm.find_next_feature(
            self.state.chrom, 
            self.state.position,
            'gene'
        )
        
        if result:
            self.state.position = result['start']
    
    def _jump_to_prev_gene(self):
        """Quick jump to previous gene using Ctrl+Left."""
        result = self.gm.find_prev_feature(
            self.state.chrom,
            self.state.position,
            'gene'
        )
        
        if result:
            self.state.position = result['start']
    
    def _save_state(self):
        """Save current browser state."""
        import json
        
        state_data = {
            'chrom': self.state.chrom,
            'position': self.state.position,
            'window_size': self.state.window_size,
            'stride': self.state.stride,
            'show_features': self.state.show_features,
            'show_amino_acids': self.state.show_amino_acids,
            'show_reverse_strand': self.state.show_reverse_strand,
            'show_rna': self.state.show_rna,
            'timestamp': __import__('time').time()
        }
        
        # Save to .genome_browser_state.json in working directory
        state_file = '.genome_browser_state.json'
        try:
            with open(state_file, 'w') as f:
                json.dump(state_data, f, indent=2)
            # Brief status message without disrupting display
            print(f"\\033[s\\033[25;1H{Colors.DIM}State saved{Colors.RESET}\\033[u", end="", flush=True)
            __import__('time').sleep(0.5)
        except Exception as e:
            print(f"\\033[s\\033[25;1H{Colors.SNP}Save failed: {e}{Colors.RESET}\\033[u", end="", flush=True)
            __import__('time').sleep(1)
    
    def _load_state(self):
        """Load previously saved browser state."""
        import json
        import os
        
        state_file = '.genome_browser_state.json'
        if not os.path.exists(state_file):
            print(f"\\033[s\\033[25;1H{Colors.DIM}No saved state found{Colors.RESET}\\033[u", end="", flush=True)
            __import__('time').sleep(0.5)
            return
        
        try:
            with open(state_file, 'r') as f:
                state_data = json.load(f)
            
            # Restore state
            self.state.chrom = state_data['chrom']
            self.state.position = state_data['position']
            self.state.window_size = state_data['window_size']
            self.state.stride = state_data['stride']
            self.state.show_features = state_data['show_features']
            self.state.show_amino_acids = state_data['show_amino_acids']
            self.state.show_reverse_strand = state_data['show_reverse_strand']
            self.state.show_rna = state_data['show_rna']
            
            print(f"\\033[s\\033[25;1H{Colors.DIM}State loaded{Colors.RESET}\\033[u", end="", flush=True)
            __import__('time').sleep(0.5)
        except Exception as e:
            print(f"\\033[s\\033[25;1H{Colors.SNP}Load failed: {e}{Colors.RESET}\\033[u", end="", flush=True)
            __import__('time').sleep(1)
    
    def _cleanup(self):
        """Clean up before exiting."""
        print("\n" + "-" * 40)
        print("Exiting genome browser")


def browse_genome(genome_manager: 'GenomeManager', 
                  chrom: Union[str, int] = 1,
                  position: int = 1000000,
                  window_size: int = 240,
                  debug=False,
                  use_gui: bool = False):
    """Convenience function to start the genome browser.
    
    Args:
        genome_manager: Initialized GenomeManager
        chrom: Starting chromosome
        position: Starting position
        window_size: Initial window size
        use_gui: If True, launch the PyQt GUI version instead of terminal version
    """
    if debug:
        logger.setLevel("DEBUG")
    
    if use_gui:
        try:
            from .genome_browser_gui import launch_genome_browser
            launch_genome_browser(genome_manager)
        except ImportError as e:
            print(f"GUI mode requires PyQt6: {e}")
            print("Falling back to terminal mode...")
            input()
            browser = InteractiveGenomeBrowser(genome_manager)
            browser.start(chrom, position, window_size)
    else:
        browser = InteractiveGenomeBrowser(genome_manager, debug = debug)
        # browser.start(chrom, position, window_size, show_second = True)
        browser.start(chrom, position, window_size)