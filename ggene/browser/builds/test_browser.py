
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
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.renderer import ArtistRenderer
from ggene.display.layout import Label, FullWidthRow, Alignment, VAlignment
from ggene.display.colors import FColors

logger = logging.getLogger(__name__)
# logger.setLevel("WARNING")
logger.setLevel("DEBUG")

@dataclass
class TestBrowserState(BaseBrowserState):    
    # feature_types:Optional[Tuple[str]] = ("gene", "exon", "CDS", "three_prime_utr", "five_prime_utr", "variant", "motif", "dfam_hit", "repeat")
    feature_types:Optional[Tuple[str]] = ("gene", "exon", "CDS", "motif", "dfam_hit", "pseudogene", "lncRNA", "ncRNA", "tf_binding")
    _detect_motifs_once:bool = False
    
    mm_height:int = 8
    seq_height:int = 5
    sc_height:int = 7
    line_height:int = 20
    text_height:int = 16
    
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
    
    def to_dict(self):
        
        out_dict = self.__dict__
        
        if self._detect_motifs_once:
            out_dict["detect_motifs"] = True
            self._detect_motifs_once = False
        else:
            out_dict["detect_motifs"] = False
        
        return out_dict
    
    
class TestBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, **kwargs):
        
        gm = GenomeManager()
        
        state = TestBrowserState(
            chrom = kwargs.get("chrom", "1"), 
            position = kwargs.get("position", 10e6), 
            window_size = kwargs.get("window_size", 240), 
            stride = kwargs.get("stride", int(kwargs.get("window_size")/8)), 
        )
        
        state.update(**kwargs)
        # print()
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        self.logger = logger
        
        dw = kwargs.get("display_width", 256)
        artists = self.build_artists(display_width = dw)
        rndr = self.build_renderer_mini(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
    
    def update(self):
        super().update()
        # logger.debug("features in iterator window:")
        
        # feats = {}
        # for feat in self.window.features:
        #     if not feat.feature_type in feats:
        #         feats[feat.feature_type] = []
            
        #     feats[feat.feature_type].append(feat)
        
        # for ft in feats:
        #     print(f"features of type {ft}: {len(feats[ft])}")
        
        # logger.debug(f"detect motifs: {self.iterator.detect_motifs}")
    
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
        
        scale = 0.01
        ps = [
            [scale, "genes", "pseudo"],
            [scale, "protein_coding", "lncrna"],
            [scale, "polya", "cpg"]
        ]
        
        num_split = 3
        tmmgp = []
        # for scale, qt1, qt2 in ps:
        for i in range(len(ps)):
            scale, qt1, qt2 = ps[i]
            seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = [qt1,qt2])
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
            # orientation = "in",
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
        
        return artists

    def build_renderer(self, artists, display_height, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height)
        
        for artist in artists:
            rndr.add_full_width_row(artist, height = 8, top_label = artist.top_label)
        
        return rndr
    
    def build_renderer_mini(self, artists, display_height, display_width, full_display_width):
        
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
        
        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
    