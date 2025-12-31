
from typing import Tuple, Optional
import random

from dataclasses import dataclass

from ggene.seqs import lambdas

from ggene.browser.base_browser import BaseBrowserState, BaseBrowser
from ggene.browser.genome_browser import BrowserState
from ggene.browser.commands import CommandRegistry, KeybindingManager
from ggene.browser.bindings.defaults import bind_defaults, elicit_input, elicit_and_change
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
class InterBrowserState(BaseBrowserState):
    mm_height = 8
    seq_height = 6
    inter_height = 40
    feature_height = 8
    
class InterBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, gm = None, **kwargs):
        
        if not gm:
            gm = GenomeManager()
        
        dw = kwargs.pop("display_width", 240)    
        state = InterBrowserState(
            chrom=kwargs.pop("chrom", "1"), 
            position=kwargs.pop("position", 10e6), 
            window_size=kwargs.pop("window_size", dw), 
            stride=kwargs.pop("stride", 20), 
            display_height=kwargs.pop("display_height", 32),
            feature_types = kwargs.pop("feature_types", tuple())
        )
        
        state.update(**kwargs)
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride, feature_types = state.feature_types)
        
        debug = kwargs.pop("debug", False)
        
        super().__init__(gm, state = state, iterator = iterator, debug=debug, **kwargs)
        
        kernel = kwargs.pop("kernel", "ATAT"*4)
        artists = self.build_artists(kernel, display_width = dw, **kwargs)
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
        
        self.change_focus = self.registry.register("change_focus", "general",takes_input = False, call_update = True)(self.change_focus)
        self.keybindings.bind("F", "change_focus")
        
        
        pass
    
    def register_default_commands(self):
        
        binds = {}
        binds.update(bind_info_commands(self, self.registry))
        binds.update(bind_save_commands(self, self.registry))
        binds.update(bind_nav_commands(self, self.registry))
        binds.update(bind_defaults(self, self.registry))
        for k, n in binds.items():
            self.keybindings.bind(k, n)

    def change_focus(self, brws, state, window):
        
        res = elicit_input("Change focus (all, SW, RY, MK):")
        intera = self.artists.get("inter_artist")
        intera.params.focus = res
        return

    def build_artists(self, kernel, display_width = 240, **kwargs):
        
        artists = {}
        
        ps = [
            [0.01, "genes", "ncRNA"],
            [0.01, "gc", "polya"],
            [0.01, "Line1", "Alu"],
        ]
        
        seq_specs = []
        for _, ss1, ss2 in ps:
            seq_specs.extend([ss1, ss2])
        
        seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = seq_specs, needs_feats = ["gene","ncRNA","lncRNA","dfam_hit"])
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
        
        seqa = SeqArtist(
            "seq_artist",
            SeqArtistParams(display_width = display_width, display_height = self.state.seq_height),
            top_label = "Sequence"
        )
        artists['seq_artist'] = seqa
        
        ia = InterArtist(
            "inter_artist", 
            InterArtistParams(display_width=display_width, display_height = self.state.inter_height, **kwargs),
            kernel = kernel, 
            top_label = "Interaction"
        )
        artists["inter_artist"] = ia

        fa = LineArtist("feature_lines", LineArtistParams(
                display_width = display_width, 
                display_height = self.state.feature_height, 
                show_ruler = True,
                feature_types = ("gene","transcript","exon","CDS","simple_repeat","dfam_hit","clinical_variant","pseudogene","lncRNA","ncRNA"),
                # feature_types = ('gene','exon', 'dfam_hit')
            ), top_label = "Genomic Features")
        artists["feature_artist"] = fa
        
        self.artists = artists
        
        return self.artists
    

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
        
        for i in range(num_mini//mini_per_row):
            row1 = rndr.add_row(height = self.state.mm_height, valign = VAlignment.CENTER)
            mas = [minimap_artists.pop(0) for n in range(mini_per_row)]
            for ra in mas:
                row1.add_artist(
                    ra,
                    width = int(display_width / mini_per_row),
                    top_label = ra.top_label,
                )
        
        seqa = artists.get("seq_artist")
        rndr.add_full_width_row(seqa, height=self.state.seq_height, top_label = seqa.top_label)
        
        ia = artists.get("inter_artist")
        rndr.add_full_width_row(ia, height=self.state.inter_height, top_label = ia.top_label)
        
        seqa = artists.get("feature_artist")
        rndr.add_full_width_row(seqa, height = self.state.feature_height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = False)

        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
