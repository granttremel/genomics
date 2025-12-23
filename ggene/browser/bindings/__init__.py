
from typing import Dict, Any

# from .defaults import default_bindings, bind_defaults, toggle_param, modify_param, elicit_input, elicit_and_change
# from .info import info_bindings, bind_info_commands
from ggene.browser.commands import CommandRegistry, Command

def register_bindings(obj, registry:CommandRegistry, bindings:Dict[str,Any]):
    
    binding = {}
    
    for _, cmd_data in bindings.items():
        key = cmd_data.pop("key", "")
        cmd = Command(**cmd_data)
        if not hasattr(obj, cmd.method):
            setattr(obj, cmd.method, cmd._bound_method)
            registry.commands[cmd.name] = cmd
            binding[key] = cmd.name
            
    return binding

__all__=["register_bindings", 
        #  "default_bindings", "bind_defaults", "info_bindings", "bind_info_commands", "toggle_param", "modify_param", "elicit_input", "elicit_and_change"
        ]
