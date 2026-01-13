
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING

from . import register_bindings
from ggene.browser.commands import CommandRegistry, Command
from ggene.browser.bindings.defaults import elicit_and_change, elicit_input

if TYPE_CHECKING:
    from ggene.browser.base_browser import BaseBrowserState, BaseBrowser
    from ggene.database.genome_iterator import BaseWindow, BaseIterator, UGenomeIterator, GenomeWindow
    from ggene.display.artists.base import BaseArtist

def parse_position(pos_str, default = 1000000):
    
    if not pos_str:
        return default
    
    fact = 1
    if pos_str[-1] == "M":
        fact = int(1e6)
        pos_str = pos_str[:-1]
    elif pos_str[-1] == "k":
        fact = int(1e3)
        pos_str = pos_str[:-1]
    
    return int(fact*float(pos_str.replace(',', '')))

def move(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', delta):
    state.position += delta
    state.seqb_delta -= delta

def move_second(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', delta):
    state.seqb_delta += delta

def move_forward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
    delta = state.stride * 8 if large else state.stride
    return move(brws, state, window, delta)

def move_backward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
    delta = state.stride * 8 if large else state.stride
    return move(brws, state, window, -delta)

def move_sec_forward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
    delta = state.stride * 8 if large else state.stride
    return move_second(brws, state, window, delta)

def move_sec_backward(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', large = False):
    delta = state.stride * 8 if large else state.stride
    return move_second(brws, state, window, -delta)

def goto_position(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    pos = state.position
    
    param_name = 'position'
    message = "Enter position"
    cast_to_type = str
    min_val = None
    max_val = None
    
    pos_str = elicit_input(message, cast_to_type, min_val, max_val)
    
    try:
        pos = parse_position(pos_str) - state.window_size//2
        modify_param(brws, state, window, param_name, pos)
    except:
        if pos_str == 'r':
            chrom, pos = brws.gm.get_random_location()
        
        state.chrom = chrom
        state.position = pos - state.window_size//2
    
    return pos
    
def switch_chromosome(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
        
    param_name = 'chrom'
    message = "Enter chromosome"
    cast_to_type = str
    min_val = None
    max_val = None
    
    return elicit_and_change(brws, state, window, param_name, message, cast_to_type, min_val, max_val)

def change_window_size(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    param_name = 'window_size'
    message = "Enter new window size"
    cast_to_type = int
    min_val = None
    max_val = None
    
    elicit_and_change(brws, state, window, param_name, message, cast_to_type, min_val, max_val)
    
    state.stride = round(state.window_size/8)
    
def change_stride(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    param_name = 'stride'
    message = "Enter new stride"
    cast_to_type = int
    min_val = None
    max_val = None
    return elicit_and_change(brws, state, window, param_name, message, cast_to_type, min_val, max_val)

def change_features(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    param_name = 'feature_types'
    message = f"Enter feature type to toggle. current: {", ".join(state.feature_types)}"
    cast_to_type = str
    
    toggle_feat = elicit_input(message, cast_to_type)
    if toggle_feat:
        return toggle_param_collection(brws, state, window, param_name, toggle_feat)
    else:
        print(f"unable to parse input {toggle_feat}")
        

def focus_on_feature(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    if hasattr(state, "_feature_types"):
        state.feature_types = state._feature_types
        del state._feature_types
        return
    else:
        featmap = {"g":"gene","t":"transcript","e":"exon","c":"CDS","v":"variant","r":"repeat","d":"dfam_hit","m":"motif"}
        
        message = f"Enter feature type to focus on (or comma separated list). g=gene, e=exon, c=cds, v=variant, r=repeat, d=dfam_hit. current: {", ".join(state.feature_types)}"
        
        inp_focus_feat = elicit_input(message)
        
        focus_feats = inp_focus_feat.split(",")
        full_feats = []
        
        for focus_feat in focus_feats:
            if focus_feat in featmap:
                focus_feat = featmap[focus_feat]
        
            if focus_feat in featmap.values():
                full_feats.append(focus_feat)
        
        if full_feats:
            state._feature_types = state.feature_types
            state.feature_types = tuple(full_feats)

def find_motifs(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    state._detect_motifs_once = True
    brws.logger.debug(f"updating motifs once.. {state.__dict__}")
    
    return state

def search_for_feature(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    chrom = elicit_input("Chromosome?", str)
    start = elicit_input("Starting position?", str)
    end = elicit_input("Ending position?", str)
    query = elicit_input("Input search query as k1=v1 k2=v2 k3=v3 ...")

    start = parse_position(start)
    end = parse_position(end, None)
    
    conds = {}
    
    condstrs = query.strip().split(" ")
    for cond in condstrs:
        kv = cond.split("=")
        if len(kv) == 2:
            k,v = kv
            k = k.strip(",").strip()
            v = v.strip(",").strip()
            
            conds[k] = v
    f = None
    for f in brws.gm.search.feature_search(chrom, start, end, conds):
        
        if not f:
            break
        
        brws.log_to_display(f"feature identified:{f}", "commands.search_for_feature")
        res = elicit_input("go to feature?", cast_to_type = bool)
    
        if res:
            break
        else:
            start = f.end + 1

    if f:
        state.chrom = f.chrom
        state.position = f.start - state.window_size//2
    else:
        brws.log_to_display("failed to identify a suitable feature.", "commands.search_for_feature")
    
    return state

def modify_param(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', param_name, param_value):
    
    if hasattr(state, param_name):
        v = getattr(state, param_name)
        
        try:
            setattr(state, param_name, param_value)
        except:
            print(f"failed to set param {param_name} from {v} to {param_value}")

def toggle_param(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', param_name):
    
    if hasattr(state, param_name):
        v = getattr(state, param_name)
        
        if isinstance(v, bool):
            setattr(state, param_name, not v)

    else:
        print(f"unable to toggle param {param_name}")

def toggle_param_collection(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', param_name, param_item):
    
    if hasattr(state, param_name):
        v = list(getattr(state, param_name))
        
        if param_item in v:
            v.remove(param_item)
        else:
            v.append(param_item)
        
        setattr(state, param_name, tuple(v))

    else:
        print(f"unable to toggle param {param_name}")

def toggle_panel(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    panel_names = []
    panels = {}
    
    i = 0
    for cell in brws.renderer.layout.cells:
        
        if cell.panel:
            aname = cell.panel.artist.name
            panel_names.append(aname)
            panels[aname] = cell.panel
            print(f"{i} {aname}: enabled = {cell.panel.enabled}")
            i += 1
    
    res = elicit_input("Select panel to toggle", str)
    
    try:
        ind = int(res)
        panel = panels[panel_names[ind]]
    except:
        panel = panels.get(res)
    
    panel.enabled ^= True

def save_state(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    brws.add_to_history()
    brws.save_state()
    
def save_to_history(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    brws.add_to_history()

def echo_key(brws, state, window):
    if hasattr(state, "_last_key"):
        print(f"last keypress: {state._last_key}")
    else:
        print("last keypress unknown")

default_bindings = {
    "right_arrow":{
        "key":"right_arrow",
        "name":"move_forward",
        "method":move_forward.__name__,
        "description":"Move forward by stride",
        "_bound_method": move_forward
    },

    "g":{
        "key":"g",
        "name":"goto",
        "method":goto_position.__name__,
        "description":"Go to position",
        "_bound_method": goto_position
    },
    "c":{
        "key":"c",
        "name":"switch_chromosome",
        "method":switch_chromosome.__name__,
        "description":"Switch to chromosome",
        "_bound_method": switch_chromosome,
        "call_update":True,
    },
    "w":{
        "key":"w",
        "name":"change_window_size",
        "method":change_window_size.__name__,
        "description":"Change window size",
        "_bound_method": change_window_size
    },
    "s":{
        "key":"s",
        "name":"change_stride",
        "method":change_stride.__name__,
        "description":"Change stride",
        "_bound_method": change_stride
    },
    "ctrl_s":{
        "key":"ctrl_s",
        "name":"save_state",
        "method":save_state.__name__,
        "description":"Save browser state",
        "_bound_method": save_state
    },
    "S":{
        "key":"S",
        "name":"search_for_feature",
        "method":search_for_feature.__name__,
        "description":"Search for feature",
        "_bound_method": search_for_feature
    },
    "h":{
        "key":"h",
        "name":"save_to_history",
        "method":save_to_history.__name__,
        "description":"Take snapshot of current view and store in history",
        "_bound_method": save_to_history
    },
    "F":{
        "key":"F",
        "name":"change_features",
        "method":change_features.__name__,
        "description":"Change feature types to display",
        "_bound_method":change_features
    },
    "f":{
        "key":"f",
        "name":"focus_on_feature",
        "method":focus_on_feature.__name__,
        "description":"Focus on a single feature type",
        "_bound_method":focus_on_feature
    },
    "m":{
        "key":"m",
        "name":"find_motifs",
        "method":find_motifs.__name__,
        "description":"Identify motifs (once)",
        "call_update":True,
        "_bound_method":find_motifs
    }
    # "m":{
    #     "key":"m",
    #     "name":"toggle_panel",
    #     "method":toggle_panel.__name__,
    #     "description":"toggle display state of panel",
    #     "_bound_method":toggle_panel
    # }
}

def bind_defaults(obj, registry:CommandRegistry):
    binding = register_bindings(obj, registry, default_bindings)
    return binding
    

def assign_default_registry(obj):
    
    cmdreg = CommandRegistry()
    
    for _, cmd_data in default_bindings.items():
        
        key = cmd_data.pop("key", "")
        cmd = Command(**cmd_data)
        
        if not hasattr(obj, cmd.method):
            
            setattr(obj, cmd.method, cmd._bound_method)
            cmdreg.register(getattr(obj, cmd.method), cmd.description, category = "default", takes_input = False)
            
    return cmdreg
