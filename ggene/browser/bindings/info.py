
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
import subprocess

from . import register_bindings
from ggene.browser.commands import CommandRegistry, Command
from .defaults import elicit_input, elicit_and_change, toggle_param, modify_param

if TYPE_CHECKING:
    from ggene.browser.base_browser import BaseBrowser,BaseBrowserState
    from ggene.database.genome_iterator import BaseWindow, BaseIterator, UGenomeIterator, GenomeWindow
    from ggene.display.artists.base import BaseArtist


genecard_frm="https://www.genecards.org/cgi-bin/carddisp.pl?gene={name}"
wiki_frm="https://en.wikipedia.org/wiki/{name}"
ensembl_frm = "http://www.ensembl.org/Homo_sapiens/geneview?gene={gene_id};db=core"

def open_genecard(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    gene = None
    for feat in window.features:
        if feat.feature_type == "gene":
            
            if feat.name:
                gene = feat
                break
    
    if gene is None:
        # logger.debug("no gene cache, returning")
        return
    # atts = gene.get("attributes",{})
    # if not atts:
    #     atts = gene.get("info",{})
    
    # name = atts.get("name","")
    
    if not gene.name:
        # logger.debug(f"no name on {gene}")
        return 
    
    # logger.debug(f"opening gene card with gene name {name}")
    subprocess.run(["firefox",genecard_frm.format(name=gene.name)])

    
def open_wiki(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    gene = None
    for feat in window.features:
        if feat.feature_type == "gene":
            
            if feat.name:
                gene = feat
                break
    
    if gene is None:
        # logger.debug("no gene cache, returning")
        return
    # atts = gene.get("attributes",{})
    # if not atts:
    #     atts = gene.get("info",{})
    
    # name = atts.get("name","")
    
    if not gene.name:
        # logger.debug(f"no name on {gene}")
        return 
        
    # logger.debug(f"opening wiki with gene name {name}")
    subprocess.run(["firefox",wiki_frm.format(name=gene.name)])
    
        
def open_ensembl(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
        
        gene = None
        for feat in window.features:
            if feat.feature_type == "gene":
                
                if feat.name and "gene_id" in feat.attributes:
                    gene = feat
                    # break
        name = feat.name
        gene_id = feat.attributes.get("gene_id","")
        
        subprocess.run(["firefox",ensembl_frm.format(gene_id=gene_id)])
        

info_bindings = {
    "r":{
        "key":"q",
        "name":"open_genecard",
        "method":open_genecard.__name__,
        "description":"Open genecard of gene in view",
        "_bound_method": open_genecard
    },
    "w":{
        "key":"w",
        "name":"open_wiki",
        "method":open_wiki.__name__,
        "description":"Open wiki page of gene in view",
        "_bound_method": open_wiki
    },
    "e":{
        "key":"e",
        "name":"open_ensembl",
        "method":open_ensembl.__name__,
        "description":"Open ensemble of gene in view",
        "_bound_method": open_ensembl
    },
}

def bind_info_commands(obj, registry:CommandRegistry):
    binding = register_bindings(obj, registry, info_bindings)
    return binding