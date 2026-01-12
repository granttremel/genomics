
from typing import List, Tuple, Dict, Any, Callable, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES
from ggene.draw.color import RESET, BOLD
from ggene.draw.color import Color
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
    emphasized_features: Tuple[str] = tuple()
    emphasized_formatter: Optional[Callable] = None
    highlight_color: int = 148
    biotypes:Tuple[str] = ("protein_coding","lncRNA","processed_pseudogene")

class TextArtist(BaseArtist):
    
    _lmargin = 2
    _rmargin = 2
    
    def __init__(self, name, params:TextArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
    
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        logger.debug(f"emphs {type(self.params.emphasized_features)}, params {type(self.params.feature_types)}, state {type(state.feature_types)}")
        
        if self.params.use_global_features and state.feature_types:
            fts = self.params.emphasized_features + state.feature_types
        else:
            fts = self.params.emphasized_features + self.params.feature_types
        all_fts = []
        for ft in fts:
            if not ft in all_fts:
                all_fts.append(ft)
        
        
        features = window.features + window.motifs
        txtlines = self.get_text_lines(features, all_fts)
        
        if not txtlines:
            txtlines = ["no features in here"]
        
        return txtlines
    
    
    def get_text_lines(self, features, feature_types):
        
        txtlines = []
        
        # pred = lambda f:f.get("gene_biotype") in self.params.biotypes or f.get("transcript_biotype") in self.params.biotypes
        
        num_feats = len(features)
        num_hidden = num_feats
        
        features_filt = self.collect_features(features, feature_types, unique_startends= True)
        features_ft = {f.feature_type:[] for f in features_filt}
        
        for f in features_filt:
            features_ft[f.feature_type].append(f)
        
        max_feat_per_type = 5
        
        for ftype in feature_types:
            
            if ftype in self.params.emphasized_features:
                col = Color.from_8bit(self.params.highlight_color).to_ansi()
            else:
                col = FColors.get_feature_type_color(ftype)
            feats = features_ft.get(ftype, [])

            if not feats:
                continue

            # if ftype == "variant":
            #     print(feats[0])

            txtlines.append(f"{BOLD}{col}{ftype}:{RESET}")
            if len(feats) > max_feat_per_type:
                flines = self.summarize_features(ftype, feats, tabs= 1)
                txtlines.extend(flines)
            else:
                for f in feats:
                    fline = self.describe_feature(f)
                    flines = self.wrap_line(fline, self.params.display_width, tabs = 1)
                    txtlines.extend(flines)
                
            num_hidden -= len(feats)
        
        txtlines.append(f"({num_hidden} features hidden)")
        
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
    
    def summarize_features(self, ftype, feats, tabs = 1):
        
        atts = self.get_feature_atts(ftype)
        if ftype in self.params.emphasized_features:
            return self.summarize_emphasized(feats, tabs = tabs)
        elif ftype in ["gene","transcript","exon","CDS"]:
            return self.summarize_genomic_features(feats, atts, tabs=tabs)
        else:
            return self.summarize_simple_features(feats, atts[0], tabs=tabs)

    def summarize_genomic_features(self, feats, atts, tabs = 1):
        
        tabstr = "  "*tabs
        other_feats = []
        genefeats = {}
        for f in feats:
            gn = f.get("gene_name")
            if gn:
                if not gn in genefeats:
                    genefeats[gn] = []
                genefeats[gn].append(f)
            else:
                other_feats.append(f)
        
        outstrs = []
        
        for gene in genefeats:
            fs = genefeats[gene]
            
            parts = []
            for att in atts:
                vals = set()
                for f in fs:
                    v = f.get(att)
                    if v:
                        vals.add(str(v))
                
                parts.append(f"{att}=[{",".join(vals)}]")
            
            outstrs.append(f"{tabstr}{gene}: {", ".join(parts)}")
        
        return outstrs
    
    def summarize_simple_features(self, feats, summary_att, tabs):
        
        tabstr = "  "*tabs
        cts = {}
        
        for f in feats:
            v = f.get(summary_att)
            if not v in cts:
                cts[v] = 0
            cts[v] += 1
        
        parts = []
        for name, ct in cts.items():
            parts.append(f"{ct}x {name}")
        
        return [tabstr+", ".join(parts)]
    
    def summarize_emphasized(self, feats, tabs = 1):
        
        fstrs = []
        
        if self.params.emphasized_formatter:
            for feat in feats:
                fstr = self.describe_emphasized(feat)
                fstrs.append(fstr)
        else:
            fstrs = self.summarize_genomic_features(feats, tabs=tabs)
        
        return fstrs
    
    def describe_emphasized(self, feat, tabs = 1):
        
        if self.params.emphasized_formatter:
            fstr = "  "*tabs + self.params.emphasized_formatter(feat)

        return fstr
        
    def describe_feature(self, feat, tabs = 1):
        
        if feat.feature_type in self.params.emphasized_features:
            return self.describe_emphasized(feat, tabs=tabs)
        
        atts = self.get_feature_atts(feat.feature_type)
        name_att = atts.pop(0)
        fname = feat.get(name_att, "?")
        
        tabstr = "  "*tabs
        
        parts = []
        
        for att in atts:
            if not hasattr(feat, att):
                v = feat.attributes.get(att, "??")
            else:
                v = getattr(feat, att, "?")
            
            if isinstance(v, float):
                vstr = format(v, "0.2f")
            else:
                vstr = str(v)
            
            parts.append(f"{att}={vstr}")
        
        return f"{tabstr}{fname}: "+", ".join(parts)
    
    def get_feature_atts(self, feature_type):
        
        if feature_type == "gene":
            return ["name","start","end","length","strand","gene_name", "description", "gene_source", "gene_biotype"]
        elif feature_type == "transcript":
            return ["name","start","end","length","transcript_name","transcript_source","transcript_biotype"]
        elif feature_type == "exon":
            return ["name","start","end","length","transcript_name","exon_number","transcript_source","transcript_biotype"]
        elif feature_type == "CDS":
            return ["name","start","end","length","id"]
        elif feature_type == "variant":
            return ["name","start","end","length","ref","genotype","qual"]
        elif feature_type == "repeat":
            return ["name","start","end","length","type","motif"]
        elif feature_type == "dfam_hit":
            return ["name", "start","end","length","bits","family_acc","e-value","bias","kimura_div","start_hmm","end_hmm"]
        elif feature_type == "tf_binding":
            return ["name","start","end","length","id","match_seq","score"]
        else:
            return ["name","start","end","length","id", "attributes"]
    