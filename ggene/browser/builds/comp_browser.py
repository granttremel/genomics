
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

from ggene.database.genome_iterator import UGenomeIterator, BaseWindow, GenomeWindow
from ggene.database.genome_manager import GenomeManager
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.renderer import ArtistRenderer
from ggene.display.layout import Label, FullWidthRow, Alignment, VAlignment
from ggene.display.colors import FColors

class DualWindow(GenomeWindow):
    
    def __init__(self, window2, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.window2 = window2
    
    @classmethod
    def from_windows(cls, window, window2):
        
        dw = DualWindow(
            window2,
            window.chrom,
            window.start_ref,
            window.end_ref,
            window.start_alt,
            window.end_alt,
            window.start_display,
            window.end_display,
            window.ref_seq,
            window.alt_seq,
            window.features,
            window.variant_deltas,
            window.variant_features,
            window.motifs,
        )
        
        return dw

@dataclass
class CompBrowserState(BaseBrowserState):
    
    seqb_delta: int = 0
    
    mm_height = 8
    seq_height = 6
    comp_height = 16
    feature_height = 8
    txt_height = 8
    
class CompBrowser(BaseBrowser):
    
    footer_text = "[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]"
    
    def __init__(self, gm = None, **kwargs):
        
        if not gm:
            gm = GenomeManager()
        
        dw = kwargs.pop("display_width", 240)    
        window_size=kwargs.pop("window_size", dw)
        state = CompBrowserState(
            chrom=kwargs.pop("chrom", "1"), 
            position=kwargs.pop("position", 10e6), 
            window_size=window_size, 
            stride=window_size//4, 
            display_height=kwargs.pop("display_height", 32),
            feature_types = kwargs.pop("feature_types", ("gene","exon","CDS","simple_repeat","dfam_hit","pseudogene","lncRNA","ncRNA", "miRNA"))
        )
        state.update(**kwargs)
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride, feature_types = state.feature_types)
        
        debug = kwargs.pop("debug", False)
        
        super().__init__(gm, state = state, iterator = iterator, debug=debug, **kwargs)
        
        artists = self.build_artists(display_width = dw, **kwargs)
        rndr = self.build_renderer(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.state.stride = self.state.window_size//4
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
    
    def update(self):
        
        update_dict = self.state.to_dict()
        self.iterator.update(**update_dict)
        window = self.iterator.get_window()
        
        window2 = self.iterator.get_window_at(update_dict["position"] + update_dict["seqb_delta"])
        
        self.window = DualWindow.from_windows(window, window2)
        
    
    def register_commands(self):
        super().register_commands()
        
        def move(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', move_delta):
            state.position += move_delta
            state.seqb_delta -= move_delta
            return state
        
        @self.registry.register("move_forward", "Move sequence window forward")
        def move_forward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            # self.logger.debug("call move_forward")
            return move(brws, state, window, delta)
        self.move_forward = move_forward

        @self.registry.register("move_backward", "Move sequence window backward")
        def move_backward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            return move(brws, state, window, -delta)
        self.move_backward = move_backward
        
        def move_second(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', move_delta):
            state.seqb_delta += move_delta
            return state
        
        @self.registry.register("move_sec_forward", "Move second sequence window forward")
        def move_sec_forward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            return move_second(brws, state, window, delta)
        self.move_sec_forward = move_sec_forward

        @self.registry.register("move_sec_backward", "Move second sequence window backward")
        def move_sec_backward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            return move_second(brws, state, window, -delta)
        self.move_sec_backward = move_sec_backward
        
        @self.registry.register("move_both_forward", "Move both sequence windows forward")
        def move_both_forward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            state = move(brws, state, window, delta)
            return move_second(brws, state, window, delta)
        self.move_both_forward = move_both_forward
        
        @self.registry.register("move_both_backward", "Move both sequence windows forward")
        def move_both_backward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            delta = state.stride * 8 if large else state.stride
            state = move(brws, state, window, -delta)
            return move_second(brws, state, window, -delta)
        self.move_both_backward = move_both_backward
        
        @self.registry.register("set_active_gene", "Set active gene")
        def set_active_gene(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
            
            max_num = 10
            genes = sorted([f for f in window.features if f.feature_type in ['gene', 'miRNA','lncRNA']], key = lambda f:-f.length)
            
            for i, gene in enumerate(genes):
                print(f"{i} {gene.name} {gene}")
                if i >= max_num:
                    break
            
            ind = elicit_input("Select gene to set as active gene:", cast_to_type = int, max_val = len(genes), default = 0)
            
            ag = genes[ind]
            
            ca = brws.artists.get("comp_artist")
            if ca:
                ca.set_active_gene(ag.name, active_gene = ag)
                print(f"set active gene to {ag.name}")
        self.set_active_gene = set_active_gene
        
        self.keybindings.bind("right_arrow", "move_forward")
        self.keybindings.bind("left_arrow", "move_backward")
        self.keybindings.bind("up_arrow", "move_sec_backward")
        self.keybindings.bind("down_arrow", "move_sec_forward")
        self.keybindings.bind("ctrl_right", "move_both_forward")
        self.keybindings.bind("ctrl_left", "move_both_backward")
        self.keybindings.bind("a", "set_active_gene")
        
    def register_default_commands(self):
        
        binds = {}
        binds.update(bind_info_commands(self, self.registry))
        binds.update(bind_save_commands(self, self.registry))
        binds.update(bind_nav_commands(self, self.registry))
        binds.update(bind_defaults(self, self.registry))
        
        for k, n in binds.items():
            if k in ['up','down']:
                continue
            self.keybindings.bind(k, n)

    def build_artists(self, display_width = 240, **kwargs):
        
        artists = {}
        
        ps = [
            [0.01, "genes", "ncRNA"],
            [0.01, "cpg", "polya"],
            [0.01, "Line1", "Alu"],
        ]
        
        seq_specs = []
        for _, ss1, ss2 in ps:
            seq_specs.extend([ss1, ss2])
        
        seq_gen, feat_gen = MapArtist.get_generators(self.gm, seq_specs = seq_specs, use_features = ["gene","ncRNA","lncRNA","dfam_hit"])
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
            SeqArtistParams(display_width = display_width, display_height = self.state.seq_height, display_strategy = "crop", max_seq_rows = 6),
            top_label = "Sequence"
        )
        artists['seq_artist'] = seqa
        
        ca = CompArtist("comp_artist",
                CompArtistParams(display_width = display_width, display_height = self.state.comp_height,
                    seq_name = "ref",
                    score_mode = "corrs",
                    heatmap_minval = 0.3,
                    heatmap_maxval = 0.75,
                    mini_chunksz = 16,
                    topk = 3,
                    q = 2,
                    corr_scale = None,
                    color_scheme = "autumn_canopy",
                    half_block = True,
                    show_minimap = True,
                )
            )
        artists["comp_artist"] = ca
        
        fa = LineArtist("feature_lines", LineArtistParams(
                display_width = display_width, 
                display_height = self.state.feature_height, 
                show_ruler = True
            ), top_label = "Genomic Features (seq a)")
        artists["feature_artist"] = fa
        
        fab = LineArtist("feature_lines_b", LineArtistParams(
                display_width = display_width, 
                display_height = self.state.feature_height, 
                show_ruler = True,
                use_second = True,
            ), top_label = "Genomic Features (seq b)")
        artists["feature_artist_b"] = fab
        
        txta = TextArtist("feature_text",TextArtistParams(
            display_width = display_width,
            display_height = self.state.txt_height,
            # use_global_features = False, 
            # feature_types = ('gene','transcript','exon','CDS')
        ), top_label = "Features")
        artists["text_artist"] = txta
        
        self.artists = artists
        
        return self.artists
    
    def build_renderer(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        # minimap_artists = artists.get("minimap_group", [])
        # num_mini = len(minimap_artists)
        # mini_per_row = 2
        
        # for i in range(num_mini//mini_per_row):
        #     row1 = rndr.add_row(height = self.state.mm_height, valign = VAlignment.CENTER)
        #     mas = [minimap_artists.pop(0) for n in range(mini_per_row)]
        #     for ra in mas:
        #         row1.add_artist(
        #             ra,
        #             width = int(display_width / mini_per_row),
        #             top_label = ra.top_label,
        #         )
        
        seqa = artists.get("seq_artist")
        rndr.add_full_width_row(seqa, height=self.state.seq_height, top_label = seqa.top_label)
        
        ca = artists.get("comp_artist")
        rndr.add_full_width_row(ca, height=self.state.comp_height, top_label = ca.top_label)
        
        fa = artists.get("feature_artist")
        rndr.add_full_width_row(fa, height = self.state.feature_height, top_label = fa.top_label, valign = VAlignment.CENTER, fixed_height = False)
        
        fab = artists.get("feature_artist_b")
        rndr.add_full_width_row(fab, height = self.state.feature_height, top_label = fab.top_label, valign = VAlignment.CENTER, fixed_height = False)
        
        txta = artists.get("text_artist")
        rndr.add_full_width_row(txta, height = self.state.txt_height, top_label = fab.top_label, valign = VAlignment.CENTER, fixed_height = True)

        return rndr
        
def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")
