
from typing import Tuple, Optional
import argparse
import random
import logging

from dataclasses import dataclass

from ggene.seqs import lambdas

from ggene.browser.base_browser import BaseBrowserState, BaseBrowser
from ggene.browser.genome_browser import BrowserState
from ggene.browser.commands import CommandRegistry, KeybindingManager
from ggene.browser.bindings.defaults import bind_defaults, elicit_input
from ggene.browser.bindings.info import bind_info_commands
from ggene.browser.bindings.save import bind_save_commands
from ggene.browser.bindings.nav import bind_nav_commands

from ggene.database.genome_iterator import UGenomeIterator
from ggene.database.genome_manager import GenomeManager
from ggene.database.annotations import get_experiment_stream
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.renderer import ArtistRenderer
from ggene.display.layout import Label, FullWidthRow, Alignment, VAlignment
from ggene.display.colors import FColors

logger = logging.getLogger(__name__)
# logger.setLevel("WARNING")
logger.setLevel("DEBUG")

@dataclass
class ExpBrowserState(BaseBrowserState):    
    # feature_types:Optional[Tuple[str]] = ("gene", "exon", "CDS", "three_prime_utr", "five_prime_utr", "variant", "motif", "dfam_hit", "repeat")
    feature_types:Optional[Tuple[str]] = ("gene", "exon", "CDS", "motif", "dfam_hit", "pseudogene", "lncRNA", "ncRNA", "tf_binding")
    detect_motifs:bool = False
    _detect_motifs_once:bool = False
    
    mm_height:int = 8
    seq_height:int = 5
    sc_height:int = 7
    line_height:int = 20
    text_height:int = 24
    log_height:int = 0
    
    def get_height(self, artist):
        
        at = type(artist).__name__
        if at == "MapArtist":
            return self.mm_height
        elif at == "SeqArtist":
            return self.seq_height
        elif at == "ScalarArtist":
            return self.sc_height
        elif at == "LineArtist":
            return self.line_height
        elif at == "TextArtist":
            return self.text_height
        elif at == "LogArtist":
            return self.log_height
    
    def to_dict(self):
        
        out_dict = self.__dict__
        
        if self._detect_motifs_once:
            out_dict["detect_motifs"] = True
            self._detect_motifs_once = False
        else:
            pass
        
        return out_dict
    
    
class ExpBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, exp_path, gm = None, **kwargs):
        
        if not gm:
            gm = GenomeManager()
        
        min_exp = kwargs.pop("min_exp", 0.5)
        filter = lambda f: f.mean_exp > min_exp
        
        self.exp_stream = get_experiment_stream(exp_path, kwargs.get("source_name","exp"), kwargs.get("feature_type", "exp_read"), filter=filter, gm=gm)
        
        logger.debug(f"{self.exp_stream.feature_type}")
        
        state = ExpBrowserState(
            chrom = kwargs.get("chrom", "1"), 
            position = kwargs.get("position", 10e6), 
            window_size = kwargs.get("window_size", 240), 
            stride = kwargs.get("stride", int(kwargs.get("window_size")/8)), 
        )
        
        state.update(**kwargs)
        state.feature_types = tuple([self.exp_stream.feature_type]) + state.feature_types
        
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        self.logger = logger
        self._log_artist = None
        
        self.logger.debug(f"feature_types: {self.state.feature_types}")
        
        dw = kwargs.get("display_width", 256)
        artists = self.build_artists(display_width = dw)
        rndr = self.build_renderer_full(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
        
        self.log_to_display("finished with init!", "TestBrowser.intialize", skip_update = True)
    
    def start(self, **kwargs):
        
        for f in self.exp_stream.stream("2", start = 1, end = 101e6):
            logger.debug(f)
            
        super().start(**kwargs)
        pass
    
    def log_to_display(self, message, caller, skip_update = False):
        
        if self._log_artist:
            self._log_artist.log(message, caller)
            if not skip_update:
                self.refresh_display()
        
    
    def update(self):
        super().update()
        
        for f in self.window.features:
            logger.debug(f)
        
        logger.debug(self.iterator.feature_types)
        
        pass
    
    def register_parameters(self):
        
        self.params.register("display.window_size", self.state.window_size)
        self.params.register("display.stride", self.state.stride)
    
    def register_commands(self):
        super().register_commands()
    
    def register_default_commands(self):
        
        binds = {}
        binds.update(bind_info_commands(self, self.registry))
        binds.update(bind_save_commands(self, self.registry))
        binds.update(bind_nav_commands(self, self.registry))
        binds.update(bind_defaults(self, self.registry))
        for k, n in binds.items():
            self.keybindings.bind(k, n)

    def build_artists(self, display_width = 240):
        
        artists = {}
        
        mm_height = self.state.mm_height
        seq_height = self.state.seq_height
        sc_height = self.state.sc_height
        line_height = self.state.line_height
        text_height = self.state.text_height
        log_height = self.state.log_height
        
        seq_feat = lambda seq, feats: sum(1 for f in feats if f.feature_type == self.exp_stream.feature_type and self.exp_stream.validate_feature(f))
        
        scale = 0.01
        ps = [
            [0.05, "genes", seq_feat],
            [0.01, "genes", seq_feat],
            [0.002, "genes", seq_feat]
        ]
        
        num_split = 3
        tmmgp = []
        for i in range(len(ps)):
            scale, qt1, qt2 = ps[i]
            seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = [qt1,qt2], needs_feats = [self.exp_stream.feature_type])
            arth = MapArtist(
                f"map_{qt1}-{qt2}",
                MapArtistParams(
                    display_width = display_width//num_split, 
                    display_height = mm_height, 
                    quantity = qt1, 
                    quantity2 = qt2,
                    scale = scale,
                    show_ruler = True,
                    show_fella = i==0
                    ), 
                    sequence_generator = seq_gen, feature_generator = feat_gen, 
                    top_label = f"{qt1}/{qt2} overmap ({scale}x)"
                )
            tmmgp.append(arth)
        
        artists["top_minimap_group"] = tmmgp
        
        sqa = SeqArtist("seq",SeqArtistParams(display_width = display_width, display_height =seq_height, display_strategy = "crop"), top_label = "Sequences:")
        artists["sequence_artist"] = sqa
        
        scalar_qts = ["correlation"]
        qt = scalar_qts[0]
        sca = ScalarArtist(f"scalar_{qt}",ScalarArtistParams(display_width = display_width, display_height = sc_height,  quantity = qt), top_label = qt)
        artists["scalar_artist"] = sca
        
        fa = LineArtist("feature_lines", LineArtistParams(
            display_width = display_width, 
            display_height = line_height, 
            show_ruler = True,
            # use_global_features = False,
            # feature_types = ('gene','transcript','exon','CDS')
        ), top_label = "Genomic Features")
        artists["feature_artist"] = fa
        
        txta = TextArtist("feature_text",TextArtistParams(
            display_width = display_width,
            display_height = text_height,
            # use_global_features = False, 
            # feature_types = ('gene','transcript','exon','CDS')
        ), top_label = "Features")
        artists["text_artist"] = txta
        
        loga = LogArtist("logs", LogArtistParams(display_width=display_width,display_height=log_height))
        artists["log_artist"] = loga
        self._log_artist = loga
        
        return artists

    def build_renderer(self, artists, display_height, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height)
        
        for artist in artists:
            rndr.add_full_width_row(artist, height = 8, top_label = artist.top_label)
        
        return rndr
    
    def build_renderer_full(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        minimap_artists = artists.get("top_minimap_group", [])
        nmm = num_map = len(minimap_artists)
        
        for i in range(num_map//nmm):
            row1 = rndr.add_row(height = self.state.mm_height, valign = VAlignment.CENTER)
            for ra in minimap_artists:
                row1.add_artist(
                    ra,
                    width = display_width // nmm,
                    top_label = ra.top_label,
                )
        
        """
        - curses (stdlib) - low-level but cross-platform terminal control
        - blessed / blessings - nicer curses wrapper
        - rich - excellent for progress bars, tables, styled output
        - textual - full TUI framework (by the rich author)
        """
        
        # for artist in artists:
        seqa = artists.get("sequence_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = height)
        
        seqa = artists.get("scalar_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = height)
        
        seqa = artists.get("feature_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER)
        
        seqa = artists.get("text_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, fixed_height = True)
        
        # loga = artists.get("log_artist")
        # height = self.state.get_height(loga)
        # rndr.add_full_width_row(loga, height = height, top_label = "Logs", fixed_height = False)
        
        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
    