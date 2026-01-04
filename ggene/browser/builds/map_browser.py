
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
class MapBrowserState(BaseBrowserState):
    display_height = 48
    feature_types = ('gene','lncRNA')
    map_height = 9
    line_height = 8
    maps_per_row = 3
    map_scales = [-1, 0.1, 0.02, 0.005, 0.001]
    map_spec_pairs = [
        ("polyy","polyr"),
        ("genes","lncRNA"),
        ("SVA","MIR"),
        ]
    
class MapBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, gm = None, **kwargs):
        
        if not gm:
            gm = GenomeManager()
            
        state = MapBrowserState(kwargs.get("chrom", "1"), kwargs.get("position", 10e6), kwargs.get("window_size", 240), kwargs.get("stride", 20), kwargs.get("display_height", 32), feature_types = ('gene','lncRNA'))
        state.update(**kwargs)
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride, feature_types = state.feature_types)
        
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        
        self.logger.debug(f"state on init: {self.state}")
        
        dw = kwargs.get("display_width", 240)
        artists = self.build_artists(display_width = dw)
        rndr = self.build_renderer(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride, feature_types = self.state.feature_types)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
    
    def register_parameters(self):
        
        self.params.register("display.window_size", self.state.window_size)
        self.params.register("display.stride", self.state.stride)
    
    def update(self):
        
        self.logger.debug(f"iterator: {self.iterator.track_features}, {self.iterator.feature_types}")
        self.logger.debug(f"state:  {self.state.feature_types}")
        
        super().update()
        
        self.logger.debug(f"iterator: {self.iterator.track_features}, {self.iterator.feature_types}")
        
        pass
    
    def register_commands(self):
        super().register_commands()
        
        # have to go through and update (re-init? ) the artists...
        
        @self.registry.register("change_qt","Change a displayed quantity", call_update = True)
        def change_qt(brws, state, window):
            
            avails = lambdas.lambda_map.keys()
            print(f"Change quantity: available options: {", ".join(avails)}")
            
            for nm in range(state.maps_per_row):
                qtpair = state.map_spec_pairs[nm]
                qttypes = ['upper','lower']
                
                new_pair = []
                for i in range(len(qtpair)):
                    res = elicit_input(f"Change quantity for column {nm}, {qttypes[i]} (currently {qtpair[i]})")
                    
                    if res:
                        new_qt = res.strip()
                        if not new_qt in lambdas.lambda_map:
                            print(f"Input {new_qt} not valid. no change")
                            new_pair.append(qtpair[i])
                        else:
                            new_pair.append(new_qt)
                    else:
                        new_pair.append(qtpair[i])
                
                state.map_spec_pairs[nm] = qtpair
        self.change_qt = change_qt
        
        self.keybindings.bind("h","change_qt")
        
            
            
    
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
        
        seq_specs = []
        for ss1, ss2 in self.state.map_spec_pairs:
            seq_specs.extend([ss1, ss2])
        
        seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = seq_specs)
        tmmgp = []
        
        for scale in self.state.map_scales:
            for qt1, qt2 in self.state.map_spec_pairs:
                arth = MapArtist(
                    f"map_{scale}_{qt1}-{qt2}",
                    MapArtistParams(
                        display_width = int(display_width/self.state.maps_per_row), 
                        display_height = self.state.map_height, 
                        quantity = qt1, 
                        quantity2 = qt2,
                        scale = scale,
                        show_ruler = True), 
                        sequence_generator = seq_gen, feature_generator = feat_gen, 
                        top_label = f"{qt1}/{qt2} overmap ({scale}x)"
                    )
                tmmgp.append(arth)
        
        artists["map_artists"] = tmmgp
        
        fa = LineArtist("feature_lines", LineArtistParams(
            display_width = display_width, 
            display_height = self.state.line_height, 
            show_ruler = True,
            use_global_features = True,
            feature_types = ('gene','lncRNA')
            ), top_label = "Genomic Features")
        artists["feature_artist"] = fa
        
        return artists

    def build_renderer(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        minimap_artists = artists.get("map_artists", [])
        
        num_mini = len(minimap_artists)
        
        for i in range(num_mini//self.state.maps_per_row):
            row1 = rndr.add_row(height = self.state.map_height, valign = VAlignment.CENTER)
            mas = [minimap_artists.pop(0) for n in range(self.state.maps_per_row)]
            for ra in mas:
                row1.add_artist(
                    ra,
                    width = int(display_width/self.state.maps_per_row), 
                    top_label = ra.top_label,
                )
        
        # seqa = artists.get("scalar_artist")
        # height = self.state.get_height(seqa)
        # rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = height)
        
        seqa = artists.get("feature_artist")
        rndr.add_full_width_row(seqa, height = self.state.line_height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = True)

        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
