
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
import random

from ggene.display.colors import FColors
from . import register_bindings
from ggene.browser.commands import CommandRegistry, Command
from .defaults import elicit_input, elicit_and_change, toggle_param, modify_param

if TYPE_CHECKING:
    from ggene.browser.base_browser import BaseBrowserState, BaseBrowser
    from ggene.database.genome_iterator import BaseWindow, BaseIterator, UGenomeIterator, GenomeWindow
    from ggene.display.artists.base import BaseArtist


genecard_frm="https://www.genecards.org/cgi-bin/carddisp.pl?gene={name}"
wiki_frm="https://en.wikipedia.org/wiki/{name}"
ensembl_frm = "http://www.ensembl.org/Homo_sapiens/geneview?gene={gene_id};db=core"

def save_view(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    view = brws._rendered_view
    
    if not view:
        return
    
    ftag = elicit_input("tag for filename?")
    if not ftag:
        ftag = str(random.randint(0, 9999))
    
    with open(f"./data/browser_views/view_{ftag}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")


save_bindings = {
    "ctrl_d":{
        "key":"ctrl_d",
        "name":"save_view",
        "method":save_view.__name__,
        "description":"Save view to file",
        "_bound_method": save_view
    },
}

def bind_save_commands(obj, registry:CommandRegistry):
    binding = register_bindings(obj, registry, save_bindings)
    return binding