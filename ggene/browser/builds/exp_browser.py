
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
from ggene.database.databases import *
from ggene.seqs.lambdas import make_custom_lambda
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.renderer import ArtistRenderer
from ggene.display.layout import Label, FullWidthRow, Alignment, VAlignment
from ggene.display.colors import FColors

logger = logging.getLogger(__name__)
logger.setLevel("WARNING")
# logger.setLevel("DEBUG")

@dataclass
class ExpBrowserState(BaseBrowserState):
    sources:Optional[Tuple[str]] = ('genes','dfam','repeats')
    feature_types:Optional[Tuple[str]] = ("gene", "exon", "pseudogene", "lncRNA", "miRNA", "ncRNA")
    active_streams:Optional[Tuple[str]] = tuple()
    _active_streams:Optional[Tuple[str]] = tuple()
    detect_motifs:bool = False
    _detect_motifs_once:bool = False
    
    mm_height:int = 8
    expmap_height:int = 6
    seq_height:int = 5
    sc_height:int = 7
    line_height:int = 16
    text_height:int = 16
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
    
    def __init__(self, exp_keys, gm = None, **kwargs):
        self.state: ExpBrowserState # type: ignore[assignment]
        
        if not gm:
            gm = GenomeManager()
        
        self.exp_keys = exp_keys
        self.exp_ftypes = []
        self.exp_streams = {}
        self.exp_filter = kwargs.get("exp_filter", None)
        
        wsz = kwargs.pop("window_size", 256)
        state = ExpBrowserState(
            chrom = kwargs.pop("chrom", "1"), 
            position = kwargs.pop("position", 10e6), 
            window_size = wsz, 
            stride = kwargs.pop("stride", round(wsz/8)), 
        )
        state.update(**kwargs)
        
        iterator = UGenomeIterator(gm, state.chrom, state.position, window_size = state.window_size, stride = state.stride)
        
        super().__init__(gm, state = state, iterator = iterator, **kwargs)
        self.logger = logger
        self._log_artist = None
        
        self.initialize_streams(kwargs.get("exp_keys", self.exp_keys), exp_suffix = kwargs.get("exp_suffix",""), exp_parent_dir = kwargs.get("exp_parent_dir", ""))
        
        self.logger.debug(f"feature_types: {self.state.feature_types}")
        
        dw = kwargs.get("display_width", 256)
        artists = self.build_artists(display_width = dw)
        rndr = self.build_renderer(artists, display_height = 32, display_width = dw, full_display_width = dw)
        
        self.set_renderer(rndr)
    
    def initialize_streams(self, exp_keys, exp_suffix = "", exp_parent_dir = ""):
        
        srcs = []
        fts = []
        
        db = find_and_load_database(gm = self.gm, feature_type = "auto", has_suffix = exp_suffix, parent_dir = exp_parent_dir)
        fts.append(db.feature_type)
        self.exp_streams[db.source_name] = db
        srcs.append(db.source_name)
        for k in exp_keys:
            # db = find_and_load_database(gm = self.gm, feature_type = "auto", in_name = k)
            db = find_and_load_database(gm = self.gm, feature_type = "auto", in_name = k, has_suffix = exp_suffix, parent_dir = exp_parent_dir)
            if not db:
                logger.debug(f"failed to load database with key {k}")
                continue
            fts.append(db.feature_type)
            self.exp_streams[db.source_name] = db
            srcs.append(db.source_name)
        self.exp_ftypes = tuple(fts)
        self.state.feature_types = tuple(set(self.exp_ftypes + self.state.feature_types))
        # self.state.sources = tuple(list(self.state.sources))
        self.state.active_streams = tuple(srcs)
        
    def initialize(self, **kwargs):
        
        self.state.update(**kwargs)
        self.iterator.update(chrom = self.state.chrom, position=self.state.position, window_size = self.state.window_size, stride = self.state.stride, sources = self.state.sources, feature_types = self.state.feature_types)
        self.renderer.update(**kwargs)
        self.renderer.update_all(**kwargs)
        
        self.log_to_display("finished with init!", "TestBrowser.intialize", skip_update = True)
    
    def start(self, **kwargs):
            
        super().start(**kwargs)
    
    def log_to_display(self, message, caller, skip_update = False):
        
        if self._log_artist:
            self._log_artist.log(message, caller)
            if not skip_update:
                self.refresh_display()
        
    
    def update(self):
        
        # self.logger.debug(f"update: {self.state.active_streams}")
        # self.logger.debug(f"update: {self.state.__dict__}")
        
        update_dict = self.state.to_dict()
        update_dict["sources"] = tuple(set(update_dict.get("active_streams", tuple()) + update_dict["sources"]))
        
        logger.debug(f"all sources: {update_dict.get("sources")}")
        logger.debug(f"all feature types: {update_dict.get("feature_types")}")
        # logger.debug(f"update: sources {self.state.sources}")
        # logger.debug(f"update: active_streams {self.state.active_streams}")
        # logger.debug(f"update: _active_streams {self.state._active_streams}")
        
        self.iterator.update(**update_dict)
        self.window = self.iterator.get_window(filter = None)
        
        logger.debug(f"iterator sources: {self.iterator.sources}")
        logger.debug(f"iterator feature types: {self.iterator.feature_types}")
    
    def register_parameters(self):
        
        self.params.register("display.window_size", self.state.window_size)
        self.params.register("display.stride", self.state.stride)
    
    def register_commands(self):
        super().register_commands()
        
        
        @self.registry.register("toggle_stream", "Toggle an experimental data stream", call_update = True)
        def toggle_stream(brws, state, window):
            
            sns = list(self.exp_streams.keys())
            
            for i, source_name in enumerate(sns):
                print(i, source_name, source_name in state.active_streams)
            
            msg = f"Select stream to toggle."
            res = elicit_input(msg, int)
            
            if res in ["", None]:
                return
            
            else:
                sn_toggle = sns[res]
                state.active_streams = tuple(sn for sn in sns if not sn==sn_toggle)
            
            # logger.debug(f"toggle_stream {state.sources} {state.active_streams}")
                
            return state
        self.toggle_stream = toggle_stream
        
        @self.registry.register("focus_on_stream", "Focus on an experimental data stream", call_update = True)
        def focus_on_stream(brws, state, window):
            
            if state._active_streams:
                state.active_streams = state._active_streams
                state._active_streams = tuple()
                return state
            
            sns = list(self.exp_streams.keys())
            
            for i, source_name in enumerate(sns):
                print(i, source_name, source_name in state.active_streams)
            
            msg = f"Select stream to bring into focus."
            res = elicit_input(msg, int)
            
            logger.debug(f"focus_on_stream raw res {res}")
            
            if res in ["", None]:
                return
            
            else:
                sn_toggle = sns[res]
                logger.debug(f"input {res} selected {sn_toggle}")
                state._active_streams = state.active_streams
                state.active_streams = (sn_toggle,)
                
            # logger.debug(f"focus_on_stream {state.sources} {state.active_streams}")
            
            return state
        self.focus_on_stream = focus_on_stream
        
        self.keybindings.bind("x", "focus_on_stream")
        self.keybindings.bind("X", "toggle_stream")
        
    
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
        
        # seq_feat = lambda seq, feats: sum(1 for f in feats if f.feature_type == self.exp_stream.feature_type and self.exp_stream.validate_feature(f))
        
        scale = 0.01
        ps = [
            [0.05, "genes", "ncRNA"],
            [0.01, "lncRNA", "miRNA"],
            [0.002, "Alu", "Line1"]
        ]
        
        num_split = 3
        tmmgp = []
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
                        show_fella = False
                    ), 
                    sequence_generator = seq_gen, feature_generator = feat_gen, 
                    top_label = f"{qt1}/{qt2} overmap ({scale}x)"
                )
            tmmgp.append(arth)
        
        artists["top_minimap_group"] = tmmgp
        
        seq_spec = "exp_feats"
        if self.exp_filter:
            filt = self.exp_filter
        else:
            filt = lambda f:True
        seq_lam = lambda seq, feats: sum(1 for f in feats if f.feature_type in self.exp_ftypes and filt(f))
        source_names = list(self.exp_streams.keys())
        seq_fn = make_custom_lambda(seq_spec, seq_lam, feature_types = self.exp_ftypes, source_names = source_names)
        sg, fg = MapArtist.get_generators(self.gm, [seq_fn], use_sources = True)
        
        mm_scales = [-1, 0.1, 0.01]
        expmas = []
        
        for scale in mm_scales:
            expma = MapArtist(
                "exp_map",
                params = MapArtistParams(
                    display_width = round(display_width/len(mm_scales)),
                    display_height = self.state.expmap_height,
                    quantity = seq_lam,
                    scale = scale,
                    show_ruler = True,
                    show_fella = True,
                ),
                sequence_generator = sg, feature_generator = fg,
                top_label = f"Experiment Map ({scale})"
            )
            expmas.append(expma)
        artists["experiment_map"] = expmas
        
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
            emphasized_features=self.exp_ftypes,
        ), top_label = "Genomic Features")
        artists["feature_artist"] = fa
        
        cols = next(iter(self.exp_streams.values())).get_headers()
        def emphasized_formatter(f):
            
            parts = []
            for att in cols:
                if att == "name":
                    continue
                
                v = getattr(f, att)
                if not v:
                    continue
                
                if isinstance(v, float):
                    vstr = format(v, "0.3f")
                else:
                    vstr = str(v)
                
                parts.append(f"{att}={vstr}")
            
            return f"{f.name}: {", ".join(parts)}"
        
        txta = TextArtist("feature_text",TextArtistParams(
            display_width = display_width,
            display_height = text_height,
            emphasized_features=self.exp_ftypes,
            emphasized_formatter=emphasized_formatter,
        ), top_label = "Features")
        artists["text_artist"] = txta
        
        loga = LogArtist("logs", LogArtistParams(display_width=display_width, display_height=log_height))
        artists["log_artist"] = loga
        self._log_artist = loga
        
        return artists
    
    def build_renderer(self, artists, display_height, display_width, full_display_width):
        
        rndr = ArtistRenderer(full_display_width, display_height, row_marker = "─")
        
        minimap_artists = artists.get("top_minimap_group", [])
        nmm = num_map = len(minimap_artists)
        
        row1 = rndr.add_row(height = self.state.mm_height, valign = VAlignment.CENTER)
        for i in range(round(num_map/nmm)):
            for ra in minimap_artists:
                row1.add_artist(
                    ra,
                    width = round(display_width / nmm),
                    top_label = ra.top_label,
                )
        
        """
        - curses (stdlib) - low-level but cross-platform terminal control
        - blessed / blessings - nicer curses wrapper
        - rich - excellent for progress bars, tables, styled output
        - textual - full TUI framework (by the rich author)
        """
        
        expmas = artists.get("experiment_map")
        
        row = rndr.add_row(height = self.state.expmap_height, valign = VAlignment.CENTER)
        for expma in expmas:
            row.add_artist(expma, width = round(display_width / len(expmas)), top_label = expma.top_label)
        
        seqa = artists.get("sequence_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = True)
        
        seqa = artists.get("scalar_artist")
        height = self.state.get_height(seqa)
        rndr.add_full_width_row(seqa, height = height, top_label = seqa.top_label, valign = VAlignment.CENTER, fixed_height = True)
        
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
    