
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
class JustLinesBrowserState(BrowserState, BaseBrowserState):
    # mm_height = 8
    feature_height = 40
    
class JustLinesBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, **kwargs):
        
        self.gm = GenomeManager()
        
        state = JustLinesBrowserState(kwargs.get("chrom", "1"), kwargs.get("position", 10e6), kwargs.get("window_size", 240), kwargs.get("stride", 20), kwargs.get("display_height", 40))
        state.update(**kwargs)
        iterator = UGenomeIterator(self.gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(state = state, iterator = iterator, **kwargs)
        
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
        
        # seqa = artists.get("scalar_artist")
        # height = self.state.get_height(seqa)
        # rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = height)
        
        fa = artists.get("feature_artist")
        rndr.add_full_width_row(fa, height = self.state.feature_height, top_label = fa.top_label, valign = VAlignment.CENTER, fixed_height = True)

        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
