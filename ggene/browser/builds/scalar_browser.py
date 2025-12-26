
from typing import Tuple, Optional
import random

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

@dataclass
class ScalarBrowserState(BrowserState, BaseBrowserState):
    mm_height = 8
    feature_height = 24
    
class ScalarBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, **kwargs):
        
        gm = GenomeManager()
        
        state = ScalarBrowserState(kwargs.get("chrom", "1"), kwargs.get("position", 10e6), kwargs.get("window_size", 240), kwargs.get("stride", 20), kwargs.get("display_height", 32))
        state.update(**kwargs)
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        
        dw = kwargs.get("display_width", 240)
        artists = self.build_artists(display_width = dw)
        rndr = self.build_renderer_mini(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
    
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
        
        seq_L2 = lambda seq, feats: lambdas._seq_te(seq, feats, "L2")
        seq_L3 = lambda seq, feats: lambdas._seq_te(seq, feats, "L3")
        
        ps = [
            # [-1, "features", "genes"],
            [0.05, "genes", "gc"],
            [0.05, "genes", "gc"],
            # [-1, "gc", "ag"],
            # [0.05, "polya", "cpg"],
            # [0.05, "protein_coding", "lncrna"],
            [0.05, "Alu", "Line1"],
            [0.05, "AluJ", "AluY"],
            [0.05, "AluS", "cpg"],
            [0.05, "L1M", "L1P"],
            [0.05, "SVA", "LTR"],
            [0.05, "MIR", "L2"],
            # [0.05, "TTCTT", "TTGA"],
            # [0.05, "TATAAT", "TATATC"],
            # [0.05, "AATTAAG", "ATGT"],
            # [0.05, "GAGGAT", "CCCTGC"],
            
            
            # [0.05, "alu", "line1"],
            # [0.05, "motifs", "simple_repeats"],
            # [0.05, seq_L2, seq_L3]
        ]
        
        seq_specs = []
        for _, ss1, ss2 in ps:
            seq_specs.extend([ss1, ss2])
        
        seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = seq_specs)
        num_split = 3
        tmmgp = []
        for scale, qt1, qt2 in ps:
            arth = MapArtist(
                f"map_{qt1}-{qt2}",
                MapArtistParams(
                    display_width = display_width//num_split, 
                    display_height = self.state.mm_height, 
                    quantity = qt1, 
                    quantity2 = qt2,
                    scale = scale,
                    show_ruler = True), 
                    sequence_generator = seq_gen, feature_generator = feat_gen, 
                    top_label = f"{qt1}/{qt2} overmap ({scale}x)"
                )
            tmmgp.append(arth)
        
        artists["minimap_group"] = tmmgp

        # scalar_qts = ["correlation"]
        # qt = scalar_qts[0]
        # sca = ScalarArtist(f"scalar_{qt}",ScalarArtistParams(display_width = display_width, display_height = sc_height,  quantity = qt), top_label = qt)
        # artists["scalar_artist"] = sca
        
        fa = LineArtist("feature_lines", LineArtistParams(
            display_width = display_width, 
            display_height = self.state.feature_height, 
            show_ruler = True,
            use_global_features = False,
            feature_types = ('gene','exon', 'dfam_hit')), top_label = "Genomic Features")
        artists["feature_artist"] = fa
        
        return artists
    

    def build_renderer(self, artists, display_height, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height)
        
        for artist in artists:
            rndr.add_full_width_row(artist, height = 8, top_label = artist.top_label)
        
        return rndr
    
    def build_renderer_mini(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        minimap_artists = artists.get("minimap_group", [])
        num_mini = len(minimap_artists)
        mini_per_row = 2
        
        if True:
            for i in range(num_mini//mini_per_row):
                row1 = rndr.add_row(height = self.state.mm_height, valign = VAlignment.CENTER)
                mas = [minimap_artists.pop(0) for n in range(mini_per_row)]
                for ra in mas:
                    row1.add_artist(
                        ra,
                        width = display_width // mini_per_row,
                        top_label = ra.top_label,
                    )
                    # print(f"added artist {ra.name} to mini map")
        
        # seqa = artists.get("scalar_artist")
        # height = self.state.get_height(seqa)
        # rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = height)
        
        seqa = artists.get("feature_artist")
        rndr.add_full_width_row(seqa, height = self.state.feature_height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = False)

        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
