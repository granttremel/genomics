
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class TextArtistParams(BaseArtistParams):
    display_width:int = 256
    display_height:int = 8
    use_global_features:bool = True
    feature_types:Tuple[str] = ("gene","exon","CDS","motif","dfam_hit","repeat", "variant")

class TextArtist(BaseArtist):
    
    _lmargin = 2
    _rmargin = 2
    
    def __init__(self, name, params:TextArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
    
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        if self.params.use_global_features and state.feature_types:
            fts = state.feature_types
        else:
            fts = self.params.feature_types
        
        txtlines = self.get_text_lines(window.features, fts)
        
        if not txtlines:
            txtlines = ["no features in here"]
        
        return txtlines
    
    
    def get_text_lines(self, features, feature_types):
        
        txtlines = []
        
        
        features_filt = self.collect_features(features, feature_types, unique_startends= True )
        features_ft = {f.feature_type:[] for f in features_filt}
        
        for f in features_filt:
            features_ft[f.feature_type].append(f)
        
        max_feat_per_type = 5
        
        for ftype in feature_types:
            
            col = FColors.get_feature_type_color(ftype)
            feats = features_ft.get(ftype, [])
            
            if not feats:
                continue
            fi = 0
            txtlines.append(f"{FColors.BOLD}{col}{ftype}:{FColors.RESET}")
            for f in feats:
                fline = self.describe_feature(f)
                flines = self.wrap_line(fline, self.params.display_width, tabs = 1)
                txtlines.extend(flines)
                fi += 1
                if fi > max_feat_per_type:
                    break
        
        return txtlines
    
    def wrap_line(self, line:str, max_width, tabs = 0):
        
        if FColors.visible_len(line) < max_width:
            return [line]
        
        outlines = []
        
        lparts = list(line.split(" "))
        
        llen = 2*tabs
        split_lines = []
        while lparts:
            lp = lparts.pop(0)
            if llen + len(lp)+1 > max_width:
                outlines.append(" ".join(split_lines))
                split_lines = ["  "*(tabs + 1)]
                llen = len(split_lines[0])
                
            split_lines.append(lp)
            llen += len(lp) + 1
        
        if split_lines:
            outlines.append(" ".join(split_lines))
        
        return outlines
            
    def describe_feature(self, feat, tabs = 1):
        
        atts = ["feature_type","strand", "id"] + list(feat.attributes.keys())
        tabstr = "  "*tabs
        
        fname = feat.get("name","?")
        parts = []
        
        for att in atts:
            if not hasattr(feat, att):
                v = feat.attributes.get(att, "??")
            else:
                v = getattr(feat, att, "?")
            
            if isinstance(v, float):
                vstr = format(v, "0.1f")
            else:
                vstr = str(v)
            
            parts.append(f"{att}={vstr}")
        
        return f"{tabstr}{fname}: "+", ".join(parts)
    