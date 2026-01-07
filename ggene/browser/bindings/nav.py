

from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
import subprocess

from . import register_bindings
from ggene.browser.commands import CommandRegistry, Command
from .defaults import elicit_input, elicit_and_change, toggle_param, modify_param, parse_position

if TYPE_CHECKING:
    from ggene.browser.base_browser import BaseBrowser,BaseBrowserState
    from ggene.database.genome_iterator import BaseWindow, BaseIterator, UGenomeIterator, GenomeWindow
    from ggene.display.artists.base import BaseArtist

def _window_zoom(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', z):
    
    old_pos = state.position + state.window_size // 2
    new_window_size = int(state.window_size / z)
    
    if new_window_size < brws.renderer.params.display_width:
        new_window_size = brws.renderer.params.display_width
        z = state.window_size / new_window_size
    
    # state.stride = int(state.stride / z)
    state.stride = round(state.window_size / 8)
    state.window_size = new_window_size
    state.position = old_pos - state.window_size // 2

default_window_zoom = 3

def window_zoom(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    z = elicit_input("Change zoom (10 to zoom in 10x, 0.1 or -10 to zoom out 10x)", cast_to_type = float)
    
    if z < 0:
        z = 1/abs(z)
    
    _window_zoom(brws, state, window, z)

def window_zoom_in(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    z = default_window_zoom
    return _window_zoom(brws, state, window, z)

def window_zoom_out(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    z = 1/default_window_zoom
    return _window_zoom(brws, state, window, z)

default_map_zoom = 3

def _map_zoom(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow', z):
    for artist in brws.renderer.iter_artists(artist_type = "MapArtist"):
        if artist.params.scale < 0:
            continue
        scale = artist.params.scale
        new_scale = scale / z
        artist.update(scale = new_scale)

def map_zoom(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    z = elicit_input("Change minimap zoom (10 to zoom out 10x, 0.1 or -10 to zoom in 10x)", cast_to_type = float)
    
    if z < 0:
        z = 1/abs(z)
    
    _map_zoom(brws, state, window, z)
    
def map_zoom_in(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    z = default_map_zoom
    return _map_zoom(brws, state, window, z)

def map_zoom_out(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    z = 1/default_map_zoom
    return _map_zoom(brws, state, window, z)

def jump(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    res = elicit_input("Size of jump? (e.g. -20M)")
    
    pos = parse_position(res)
    if pos is not None:
        state.position += pos
    
    return state
    
    

def next_feature(brws:'BaseBrowser', state:'BaseBrowserState', window:'BaseWindow'):
    
    brws.log_to_display("Called command next_feature", "commands.next_feature")
    f = elicit_input("Feature type?", str, default = 'gene')
    
    cond = {"feature_type": f}
    
    fgen = brws.gm.search.feature_search(state.chrom, state.position, None, cond)
    
    f = next(iter(fgen))
    state.position = f.start - state.window_size//2
    
    brws.log_to_display(f"moving to feature {f}", "commands.next_feature")
    
    return state
    

nav_bindings = {
    "z":{
        "key":"z",
        "name":"window_zoom",
        "method":window_zoom.__name__,
        "description":"Zoom in or out of window by given factor",
        "_bound_method": window_zoom
    },
    
    "=":{
        "key":"=",
        "name":"window_zoom_in",
        "method":window_zoom_in.__name__,
        "description":"Zoom in to window by factor (default=5)",
        "_bound_method": window_zoom_in
    },
    
    "-":{
        "key":"-",
        "name":"window_zoom_out",
        "method":window_zoom_out.__name__,
        "description":"Zoom out of window by factor (default=5)",
        "_bound_method": window_zoom_out
    },
    "Z":{
        "key":"Z",
        "name":"map_zoom",
        "method":map_zoom.__name__,
        "description":"Zoom in or out of window by given factor",
        "_bound_method": map_zoom
    },
    
    "shift_=":{
        "key":"shift_=",
        "name":"map_zoom_in",
        "method":map_zoom_in.__name__,
        "description":"Zoom in to window by factor (default=5)",
        "_bound_method": map_zoom_in
    },
    
    "shift_-":{
        "key":"shift_-",
        "name":"map_zoom_out",
        "method":map_zoom_out.__name__,
        "description":"Zoom out of window by factor (default=5)",
        "_bound_method": map_zoom_out
    },
    "j":{
        "key":"j",
        "name":"jump",
        "method":jump.__name__,
        "description":"Jump by a given distance (e.g. -20M)",
        "_bound_method": jump
    },
    "n":{
        "key":"n",
        "name":"next_feature",
        "method":next_feature.__name__,
        "description":"Move to next feature of desired type",
        "_bound_method": next_feature
    },
}

def bind_nav_commands(obj, registry:CommandRegistry):
    binding = register_bindings(obj, registry, nav_bindings)
    return binding