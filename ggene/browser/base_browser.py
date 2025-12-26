
import subprocess
import sys
import tty
import termios
import os
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
from dataclasses import dataclass
import logging
import numpy as np
import time
from pathlib import Path
from datetime import datetime
import yaml
import uuid

from ggene import GenomeManager
from ggene.browser.commands import CommandRegistry, KeybindingManager
from ggene.browser.bindings.defaults import bind_defaults
from ggene.browser.config_menu import ParameterRegistry, ConfigMenu
from ggene.browser import utils
from ggene.draw.colors import Colors

if TYPE_CHECKING:
    from ggene.database.genome_iterator import BaseWindow, BaseIterator, UGenomeIterator, GenomeWindow
    from ggene.display.artists.base import BaseArtist
    from ggene.display.renderer import ArtistRenderer

logger = logging.getLogger(__name__)

@dataclass
class BaseBrowserState:

    def copy(self):
        return type(self)(**self.__dict__)
    
    def update(self, **kwargs):
        for k,v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)

class BaseBrowser:
    """Base class for interactive genome browsers with extensible keybindings."""

    footer_text = "press a key to do stuff!"

    def __init__(self, genome_manager, *args, **kwargs):

        self.debug = kwargs.get("debug", False)
        
        if self.debug:
            logger.setLevel("DEBUG")
        else:
            logger.setLevel("WARNING")
        self.session_id = str(uuid.uuid4())[:6]
        
        self.gm:GenomeManager = genome_manager
        
        load_state = kwargs.get("load",False)
        if load_state:
            self.load_last_state()
        else:
            self.state:BaseBrowserState = kwargs.get("state", BaseBrowserState())
            self.window:BaseWindow = None

            self.iterator:'BaseIterator' = kwargs.get("iterator")
        
        self.renderer: Union['BaseArtist','ArtistRenderer'] = kwargs.get("renderer")
        
        # Command system
        self.registry = CommandRegistry()
        self.keybindings = KeybindingManager(self.registry)

        # Parameter configuration system
        self.params = ParameterRegistry()
        self.register_parameters()  # Subclasses override this

        self.register_default_commands()
        self.register_commands()  # Subclasses override this
        self.setup_default_keybindings()

        # Load custom keybindings if provided
        keybinding_file = kwargs.get("keybinding_file")
        if keybinding_file and Path(keybinding_file).exists():
            self.keybindings.load_from_file(keybinding_file)

        self._rendered_view:List[str] = []

        # History tracking
        self.history:List[Dict[str, Any]] = []
    
    def set_state(self, state:BaseBrowserState):
        self.state = state
        
    def set_iterator(self, iterator:'BaseIterator'):
        self.iterator = iterator
    
    def set_renderer(self, renderer:'ArtistRenderer'):
        self.renderer = renderer
    
    def start(self, **kwargs):

        self.initialize(**kwargs)

        try:
            self.run_browser()
        except KeyboardInterrupt:
            self.cleanup()
    
    def run_browser(self):
        """Main browser loop."""

        while True:
            self.update()
            self.render()
            self.display()

            key = self.get_input()

            # Handle input and check for quit
            result = self.handle_input(key)
            if result == "quit":
                self.add_to_history()
                break

        self.cleanup()
    
    def initialize(self, **kwargs):
        """ override in subclass """
        
        # re-initialize state
        # initialize iterator
        # initialize artist(s)
        
    
    def update(self):
        
        # update iterator
        self.iterator.update(**self.state.__dict__)
        self.window = self.iterator.get_window()
        
        # render
        # self.render()
        
    
    def render(self):
        self._rendered_view = self.renderer.render(self.state, self.window, footer_text = self.footer_text)
    
    def display(self):
        
        if not self.debug:
            os.system('clear' if os.name == 'posix' else 'cls')
        
        for line in self._rendered_view:
            print(line + Colors.RESET)
        
    
    def get_input(self):
        k = utils.get_user_keypress()
        self.state._last_key = k
        return k
        
    def handle_input(self, key: str):
        """Handle user input by executing bound command.

        Args:
            key: Key pressed by user

        Returns:
            Result from command execution, or None
        """
        command_name = self.keybindings.get_command_for_key(key)

        if not command_name:
            logger.debug(f"No command bound to key: {repr(key)}")
            return None

        command = self.registry.get_command(command_name)
        if not command:
            logger.warning(f"Command not found in registry: {command_name}")
            return None
        
        self.registry.get_bound_method(self, command)
        
        # Get the method from the browser instance
        method_name = command.method
        if not hasattr(self, method_name):
            logger.error(f"Browser has no method: {method_name}")
            return None

        method = getattr(self, method_name)

        # Execute the command
        try:
            # Commands receive state and window, and should return updated state
            result = method(self, self.state, self.window)

            # if command.call_update:
            #     self.update()
                # self.render()
                # self.display()

            # If command returns a new state, update it
            if isinstance(result, BaseBrowserState):
                self.state = result

            return result

        except Exception as e:
            logger.error(f"Error executing command {command_name}: {e}")
            if self.debug:
                raise
            return None

    def register_parameters(self):
        """Register configurable parameters. Override in subclass."""
        # Example: self.params.register("display.window_size", self.state.window_size)
        pass

    def register_default_commands(self):

        binds = bind_defaults(self, self.registry)

        for k, n in binds.items():
            self.keybindings.bind(k, n)


    def register_commands(self):
        """Register available commands. Override in subclass to add commands."""
        # Base commands that all browsers should have
        @self.registry.register("quit", "Quit the browser", category="general")
        def quit(browser, state, window):
            # Special handling - will break the loop
            return "quit"

        self.quit = quit

        @self.registry.register("show_help", "Show help", category="general")
        def show_help_cmd(browser, state, window):
            self.show_help()
            return state

        self._show_help_cmd = show_help_cmd

        @self.registry.register("open_config_menu", "Open configuration menu", category="general")
        def open_config_menu(browser, state, window):
            self.open_config_menu()
            return state

        self._open_config_menu = open_config_menu

    def setup_default_keybindings(self):
        """Setup default keybindings. Override to customize."""
        # Bind basic keys
        self.keybindings.bind('q', 'quit')
        self.keybindings.bind('?', 'show_help')
        self.keybindings.bind('p', 'open_config_menu')  # Ctrl+P for preferences

    def show_help(self):
        """Display help screen with keybindings."""
        help_text = self.keybindings.generate_help_text()
        print(help_text)
        input("\nPress Enter to continue...")

    def open_config_menu(self, state, window):
        """Open interactive configuration menu."""
        menu = ConfigMenu(self.params)
        menu.show()

        # After menu closes, reinitialize components with new parameters
        self.reload_config()

    def reload_config(self):
        """Reload components after configuration changes.

        Override in subclass to reinitialize components with new parameter values.
        For example:
            - Reinitialize iterator with new window size
            - Reinitialize artist with new color scheme
            - Update state with new values
        """
        pass

    def cleanup(self):
        """Cleanup on browser exit. Override to add cleanup logic."""
        pass

    # ===== YAML Serialization =====

    def to_dict(self) -> dict:
        """
        Recursively serialize browser and all nested objects to a dictionary.

        Returns:
            Dictionary representation of the browser structure.
        """
        return {
            'browser_type': type(self).__name__,
            'debug': self.debug,
            'state': self._serialize_object(self.state),
            'iterator': self._serialize_object(self.iterator),
            'renderer': self._serialize_object(self.renderer),
            'commands': self._serialize_commands(),
            'keybindings': self._serialize_keybindings(),
            'parameters': self._serialize_parameters(),
        }

    def _serialize_object(self, obj, depth: int = 0, max_depth: int = 5) -> dict:
        """
        Recursively serialize an object to a dictionary.

        Args:
            obj: Object to serialize
            depth: Current recursion depth
            max_depth: Maximum recursion depth

        Returns:
            Dictionary representation of the object
        """
        if obj is None:
            return None

        if depth > max_depth:
            return {'_type': type(obj).__name__, '_truncated': True}

        # Handle primitive types
        if isinstance(obj, (str, int, float, bool)):
            return obj

        # Handle lists
        if isinstance(obj, (list, tuple)):
            return [self._serialize_object(item, depth + 1, max_depth) for item in obj]

        # Handle dicts
        if isinstance(obj, dict):
            return {k: self._serialize_object(v, depth + 1, max_depth) for k, v in obj.items()}

        # Handle numpy arrays
        if hasattr(obj, 'tolist'):
            return {'_type': 'ndarray', 'shape': list(obj.shape), 'dtype': str(obj.dtype)}

        # Handle dataclasses and objects with __dict__
        result = {'_type': type(obj).__name__}

        # Get serializable attributes
        if hasattr(obj, '__dataclass_fields__'):
            # Dataclass
            for field_name in obj.__dataclass_fields__:
                value = getattr(obj, field_name, None)
                result[field_name] = self._serialize_object(value, depth + 1, max_depth)
        elif hasattr(obj, '__dict__'):
            # Regular object
            for attr, value in obj.__dict__.items():
                # Skip private attributes and methods
                if attr.startswith('_') and not attr.startswith('__'):
                    continue
                # Skip callable attributes
                if callable(value):
                    continue
                result[attr] = self._serialize_object(value, depth + 1, max_depth)

        # Handle special cases for common types
        if hasattr(obj, 'params') and not 'params' in result:
            result['params'] = self._serialize_object(obj.params, depth + 1, max_depth)

        if hasattr(obj, 'artists') and not 'artists' in result:
            result['artists'] = {
                name: self._serialize_object(artist, depth + 1, max_depth)
                for name, artist in obj.artists.items()
            }

        return result

    def _serialize_commands(self) -> dict:
        """Serialize registered commands."""
        commands = {}
        for name, cmd in self.registry.commands.items():
            commands[name] = {
                'method': cmd.method,
                'description': cmd.description,
                'category': cmd.category,
                'call_update': cmd.call_update,
            }
        return commands

    def _serialize_keybindings(self) -> dict:
        """Serialize keybindings."""
        return dict(self.keybindings.bindings)

    def _serialize_parameters(self) -> dict:
        """Serialize registered parameters."""
        params = {}
        for path, param in self.params.params.items():
            params[path] = {
                'name': param.name,
                'path': param.path,
                'current_value': self._serialize_object(param.current_value, max_depth=2),
                'param_type': param.param_type.__name__ if param.param_type else None,
                'description': param.description,
            }
        return params

    def to_yaml(self, file_path: str = None) -> str:
        """
        Serialize browser to YAML format.

        Args:
            file_path: Optional path to save YAML file. If None, returns string.

        Returns:
            YAML string representation
        """
        import yaml

        data = self.to_dict()
        yaml_str = yaml.dump(data, default_flow_style=False, sort_keys=False, allow_unicode=True)

        if file_path:
            with open(file_path, 'w') as f:
                f.write(yaml_str)

        return yaml_str

    @classmethod
    def from_yaml(cls, file_path: str = None, yaml_str: str = None) -> dict:
        """
        Load browser configuration from YAML.

        Note: This returns the configuration dictionary. Full deserialization
        requires subclass-specific logic to reconstruct objects.

        Args:
            file_path: Path to YAML file
            yaml_str: YAML string (alternative to file_path)

        Returns:
            Configuration dictionary
        """
        import yaml

        if file_path:
            with open(file_path, 'r') as f:
                yaml_str = f.read()

        return yaml.safe_load(yaml_str)

    def save_config(self, file_path: str):
        """
        Save browser configuration to YAML file.

        Args:
            file_path: Path to save configuration
        """
        self.to_yaml(file_path)
        logger.info(f"Browser configuration saved to {file_path}")

    def describe(self) -> str:
        """
        Return a human-readable description of the browser structure.

        Returns:
            Multi-line string describing the browser
        """
        lines = [
            f"Browser: {type(self).__name__}",
            f"  Debug: {self.debug}",
            "",
            "State:",
        ]

        if self.state:
            for attr in dir(self.state):
                if not attr.startswith('_'):
                    value = getattr(self.state, attr)
                    if not callable(value):
                        lines.append(f"  {attr}: {value}")

        lines.append("")
        lines.append("Iterator:")
        if self.iterator:
            lines.append(f"  Type: {type(self.iterator).__name__}")

        lines.append("")
        lines.append("Renderer:")
        if self.renderer:
            lines.append(f"  Type: {type(self.renderer).__name__}")
            if hasattr(self.renderer, 'artists'):
                lines.append("  Artists:")
                for name, artist in self.renderer.artists.items():
                    lines.append(f"    {name}: {type(artist).__name__}")
                    if hasattr(artist, 'params'):
                        lines.append(f"      params: {artist.params}")

        lines.append("")
        lines.append(f"Commands: {len(self.registry._commands)}")
        lines.append(f"Keybindings: {len(self.keybindings._bindings)}")
        lines.append(f"Parameters: {len(self.params._params)}")

        return "\n".join(lines)

    # ===== History and State Management =====

    def take_snapshot(self, include_view: bool = True) -> Dict[str, Any]:
        """
        Take a snapshot of the current browser state.

        Args:
            include_view: Whether to include the rendered view in the snapshot

        Returns:
            Dictionary containing the snapshot data
        """
        snapshot = {
            'timestamp': datetime.now().isoformat(),
            'state': self._serialize_object(self.state, max_depth=3),
            'window': self._serialize_object(self.window, max_depth=3),
        }

        if include_view:
            from ggene.display.colors import FColors
            snapshot['view'] = [FColors.scrub_codes(line) for line in self._rendered_view]

        return snapshot

    def add_to_history(self, include_view: bool = True):
        """
        Add current state to the history.

        Args:
            include_view: Whether to include the rendered view in the snapshot
        """
        snapshot = self.take_snapshot(include_view=include_view)
        self.history.append(snapshot)

    def save_state(self, directory: str = "./data/browser/states") -> str:
        """
        Save the current browser state to a YAML file with timestamp.

        Args:
            directory: Directory to save state files

        Returns:
            Path to the saved file
        """
        Path(directory).mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{type(self).__name__}_state_{timestamp}.yaml"
        filepath = Path(directory) / filename

        state_data = {
            'timestamp': datetime.now().isoformat(),
            'session_id':self.session_id,
            'browser_type': type(self).__name__,
            'state': self._serialize_object(self.state, max_depth=3),
            'iterator': {
                '_type': type(self.iterator).__name__,
                'chrom': getattr(self.iterator, 'chrom', None),
                'position': getattr(self.iterator, 'position', None),
                'window_size': getattr(self.iterator, 'window_size', None),
                'stride': getattr(self.iterator, 'stride', None),
            },
            'features':self.window.features,
        }

        with open(filepath, 'w') as f:
            yaml.dump(state_data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

        logger.info(f"Browser state saved to {filepath}")
        return str(filepath)

    def load_state(self, filepath: str):
        """
        Load browser state from a YAML file and update the browser.

        Args:
            filepath: Path to the saved state file
        """
        with open(filepath, 'r') as f:
            state_data = yaml.safe_load(f)

        if 'state' in state_data:
            state_dict = state_data['state']
            for key, value in state_dict.items():
                if key != '_type' and hasattr(self.state, key):
                    setattr(self.state, key, value)

        if 'iterator' in state_data:
            iter_dict = state_data['iterator']
            if hasattr(self.iterator, 'chrom') and 'chrom' in iter_dict:
                self.iterator.chrom = iter_dict['chrom']
            if hasattr(self.iterator, 'position') and 'position' in iter_dict:
                self.iterator.position = iter_dict['position']
            if hasattr(self.iterator, 'window_size') and 'window_size' in iter_dict:
                self.iterator.window_size = iter_dict['window_size']
            if hasattr(self.iterator, 'stride') and 'stride' in iter_dict:
                self.iterator.stride = iter_dict['stride']

        self.update()

        logger.info(f"Browser state loaded from {filepath}")

    def load_last_state(self):
        
        state_dir = Path("./data/browser/states")
        self_class = type(self).__name__
        max_ts = 0
        max_state_file = ""
        
        for stf in state_dir.iterdir():
            
            parts = stf.stem.split("_")
            
            if not parts[0] == self_class:
                continue
            
            ts = int(parts[2]) * 1e7 + int(parts[3])
            
            if ts > max_ts:
                max_ts = ts
                max_state_file = str(stf)
            
            if max_state_file:
                return self.load_state(max_state_file)
            else:
                return None

    def save_history(self, directory: str = "./data/browser/histories") -> str:
        """
        Save the entire browser history to a YAML file with timestamp.
        Includes all snapshots with their views for debugging/inspection.

        Args:
            directory: Directory to save history files

        Returns:
            Path to the saved file
        """
        Path(directory).mkdir(parents=True, exist_ok=True)

        # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"browser_history_{self.session_id}.yaml"
        filepath = Path(directory) / filename

        history_data = {
            'saved_at': datetime.now().isoformat(),
            'browser_type': type(self).__name__,
            'num_snapshots': len(self.history),
            'snapshots': self.history
        }

        with open(filepath, 'w') as f:
            yaml.dump(history_data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

        logger.info(f"Browser history ({len(self.history)} snapshots) saved to {filepath}")
        return str(filepath)

    def clear_history(self):
        """Clear the browser history."""
        self.history.clear()
        logger.info("Browser history cleared")
