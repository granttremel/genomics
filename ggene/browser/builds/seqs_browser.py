
from typing import Tuple, Optional
import argparse
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
class SeqsBrowserState(BaseBrowserState):
    seq_height:int = 32
    sc_height:int = 8
    feats_height:int = 5
    text_height:int = 8
    feature_types:Optional[Tuple[str]] = ("gene", "exon", "CDS", "motif", "dfam_hit")
    
    def get_height(self, artist):
        
        at = type(artist).__name__
        if at == "MapArtist":
            return self.mm_height
        elif at == "SeqArtist":
            return self.seq_height
        elif at == "ScalarArtist":
            return self.sc_height
        elif at == "LineArtist":
            return self.feats_height
        elif at == "TextArtist":
            return self.text_height
    
    
class SeqsBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, gm = None, **kwargs):
        
        if not gm:
            gm = GenomeManager()
        
        state = SeqsBrowserState(kwargs.get("chrom", "1"), kwargs.get("position", 10e6), kwargs.get("window_size", 240), kwargs.get("stride", 20), kwargs.get("display_height", 32))
        state.update(**kwargs)
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        
        dw = kwargs.get("display_width", 240)
        artists = self.build_artists(qt1 = "gc", qt2 = "genes", display_width = dw)
        rndr = self.build_renderer(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
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

    def build_artists(self, qt1 = "gc", qt2 = "genes", display_width = 240):
        
        artists = {}
        
        
        seq_height = self.state.seq_height
        sc_height = self.state.sc_height
        text_height = self.state.text_height
        
        sqa = SeqArtist("seq",SeqArtistParams(display_width = display_width, display_height =seq_height, display_strategy = "fold"), top_label = "Sequences:")
        artists["sequence_artist"] = sqa
        
        scalar_qts = ["correlation"]
        qt = scalar_qts[0]
        sca = ScalarArtist(f"scalar_{qt}",ScalarArtistParams(display_width = display_width, display_height = sc_height,  quantity = qt), top_label = qt)
        artists["scalar_artist"] = sca
        
        fa = LineArtist("feature_lines", LineArtistParams(
            display_width = display_width, 
            display_height = self.state.feats_height, 
            show_ruler = True,
            use_global_features = False,
            feature_types = ('gene','exon')), top_label = "Genes")
        artists["feature_artist"] = fa
        
        txta = TextArtist("feature_text",TextArtistParams(
            display_width = display_width,
            display_height = text_height,
            use_global_features = False, 
            feature_types = ('gene','transcript','exon','CDS')), top_label = "Features")
        artists["text_artist"] = txta
        
        return artists
    
    def build_renderer(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        seqa = artists.get("sequence_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = True)
        
        sca = artists.get("scalar_artist")
        height = self.state.get_height(sca)
        rndr.add_full_width_row(sca, height = height, top_label = sca.top_label, valign = VAlignment.CENTER, fixed_height = False)
        
        fa = artists.get("feature_artist")
        height = self.state.get_height(fa)
        rndr.add_full_width_row(fa, height = height, top_label = fa.top_label, valign = VAlignment.CENTER, fixed_height = True)
        
        txta = artists.get("text_artist")
        height = self.state.get_height(txta)
        rndr.add_full_width_row(txta, height = height, top_label = txta.top_label)
        
        
        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
    

def main():
    parser = argparse.ArgumentParser(description='Interactive Genome Browser')
    parser.add_argument('--chrom', '-c', type=str, default='1',
                       help='Chromosome to browse (default: 1)')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--window', '-w', type=int, default=240, dest = "window_size",
                       help='Starting position (default: 1000000)')
    parser.add_argument('--stride', '-s', type=int, default=20,
                       help='Window size in base pairs (default: 80)')
    parser.add_argument('--debug', '-d', action="store_true",
                       help='debug', dest = "debug")
    parser.add_argument('--random', '-r', action="store_true",
                       help='random', dest = "random")
    
    args = parser.parse_args()
    
    if args.random:
        args.chrom, args.position = GenomeManager.get_random_location(margin = 10e6)
    
    # window_size = 256
    max_disp = 256
    display_width = max_disp
    display_height = 32
    
    brws = SeqsBrowser(window_size = args.window_size, display_width = display_width, display_height = display_height, **args.__dict__)
    brws.start(display_width = display_width, **args.__dict__)
    
    write_view(brws._rendered_view, "with_align")
    
    brws.to_yaml("./data/browser_data/test_browser_1.yaml")
    

if __name__=="__main__":
    main()

