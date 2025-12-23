"""
Command system for genome browser keybindings.

Provides a decorator-based system for registering commands that can be:
- Bound to keys via configuration
- Discovered automatically
- Serialized to/from config files
- Documented with help text
"""

from typing import Dict, Callable, Optional, Union, Any, TYPE_CHECKING
from dataclasses import dataclass
import yaml
import json
from pathlib import Path
import re

if TYPE_CHECKING:
    from ggene.browser.base_browser import BaseBrowserState, BaseWindow


# Terminal escape sequences mapped to readable names
KEY_CODES = {
    # Arrow keys
    '\x1b[A': 'up_arrow',
    '\x1b[B': 'down_arrow',
    '\x1b[C': 'right_arrow',
    '\x1b[D': 'left_arrow',

    # Shift + Arrow keys
    '\x1b[1;2A': 'shift_up_arrow',
    '\x1b[1;2B': 'shift_down_arrow',
    '\x1b[1;2C': 'shift_right_arrow',
    '\x1b[1;2D': 'shift_left_arrow',

    # Ctrl + Arrow keys
    '\x1b[1;5A': 'ctrl_up_arrow',
    '\x1b[1;5B': 'ctrl_down_arrow',
    '\x1b[1;5C': 'ctrl_right_arrow',
    '\x1b[1;5D': 'ctrl_left_arrow',

    # Alt + Arrow keys
    '\x1b[1;3A': 'alt_up_arrow',
    '\x1b[1;3B': 'alt_down_arrow',
    '\x1b[1;3C': 'alt_right_arrow',
    '\x1b[1;3D': 'alt_left_arrow',

    # Function keys
    '\x1bOP': 'f1',
    '\x1bOQ': 'f2',
    '\x1bOR': 'f3',
    '\x1bOS': 'f4',
    '\x1b[15~': 'f5',
    '\x1b[17~': 'f6',
    '\x1b[18~': 'f7',
    '\x1b[19~': 'f8',
    '\x1b[20~': 'f9',
    '\x1b[21~': 'f10',
    '\x1b[23~': 'f11',
    '\x1b[24~': 'f12',

    # Other special keys
    '\x1b[H': 'home',
    '\x1b[F': 'end',
    '\x1b[2~': 'insert',
    '\x1b[3~': 'delete',
    '\x1b[5~': 'page_up',
    '\x1b[6~': 'page_down',
    '\x7f': 'backspace',
    '\x1b': 'escape',
    '\t': 'tab',
    '\n': 'enter',
    '\r': 'enter',
    ' ': 'space',
}

# Reverse mapping for converting readable names back to codes
CODE_NAMES = {v: k for k, v in KEY_CODES.items()}


def normalize_key(key: str) -> str:
    """Convert a key to its canonical form (escape sequence).

    Args:
        key: Either an escape sequence or readable name like 'right_arrow'

    Returns:
        The escape sequence for the key
    """
    # If it's already an escape sequence, return as-is
    if key in KEY_CODES:
        return key

    # Try to find it in the reverse mapping
    if key in CODE_NAMES:
        return CODE_NAMES[key]

    # Handle ctrl+letter combinations (e.g., 'ctrl_s' -> '\x13')
    ctrl_match = re.match(r'^ctrl_([a-z])$', key.lower())
    if ctrl_match:
        letter = ctrl_match.group(1)
        # Ctrl+letter is represented as ASCII value - 96 (or ord(letter) - ord('a') + 1)
        return chr(ord(letter) - ord('a') + 1)

    # Handle alt+letter combinations (e.g., 'alt_s' -> '\x1bs')
    alt_match = re.match(r'^alt_([a-z])$', key.lower())
    if alt_match:
        letter = alt_match.group(1)
        return f'\x1b{letter}'

    # If no conversion found, return as-is (might be a single character)
    return key


def get_readable_key(key: str) -> str:
    """Convert an escape sequence to a readable name.

    Args:
        key: An escape sequence

    Returns:
        Readable name for the key, or the key itself if no mapping exists
    """
    return KEY_CODES.get(key, key)


@dataclass
class Command:
    """Represents a browser command that can be bound to keys."""
    name: str
    method: str  # Method name to call on the browser
    description: str
    category: str = "general"
    takes_input: bool = False  # Whether command prompts for user input
    call_update:bool = True
    _bound_method: Callable = None

    def to_dict(self):
        """Serialize command to dict (for config file)."""
        return {
            'method': self.method,
            'description': self.description,
            'category': self.category,
            'takes_input': self.takes_input
        }


class CommandRegistry:
    """Registry of available browser commands."""

    def __init__(self):
        self.commands: Dict[str, Command] = {}

    def register(self, name: str, description: str, category: str = "general",
                 takes_input: bool = False, call_update = True):
        """Decorator to register a command.

        Usage:
            @registry.register("move_forward", "Move forward in genome", category="navigation")
            def _move_forward(self, state, window):
                state.position += state.stride
                return state
        """
        def decorator(method: Callable):
            command = Command(
                name=name,
                method=method.__name__,
                description=description,
                category=category,
                takes_input=takes_input
            )
            self.commands[name] = command
            return method
        return decorator

    def get_command(self, name: str) -> Optional[Command]:
        """Get command by name."""
        return self.commands.get(name)

    def get_bound_method(self, obj, command: Command):
        
        if command._bound_method:
            return command._bound_method
        else:
            if hasattr(obj, command.method):
                return getattr(obj, command.method)
        
        return None

    def get_by_category(self, category: str) -> Dict[str, Command]:
        """Get all commands in a category."""
        return {
            name: cmd for name, cmd in self.commands.items()
            if cmd.category == category
        }

    def list_commands(self) -> Dict[str, Command]:
        """List all registered commands."""
        return self.commands.copy()


class KeybindingManager:
    """Manages keybindings and their mapping to commands."""

    def __init__(self, registry: CommandRegistry):
        self.registry = registry
        self.bindings: Dict[str, str] = {}  # normalized_key -> command_name

    def bind(self, key: str, command_name: str):
        """Bind a key to a command.

        Args:
            key: Key (e.g., 'q', 'right_arrow', 'ctrl_s', '\x1b[C')
            command_name: Name of registered command
        """
        if command_name not in self.registry.commands:
            raise ValueError(f"Unknown command: {command_name}")

        # Normalize the key to its escape sequence form
        normalized_key = normalize_key(key)
        self.bindings[normalized_key] = command_name


        
    def bind_multiple(self, keys: list[str], command_name: str):
        """Bind multiple keys to the same command."""
        for key in keys:
            self.bind(key, command_name)

    def unbind(self, key: str):
        """Remove a keybinding."""
        normalized_key = normalize_key(key)
        self.bindings.pop(normalized_key, None)

    def get_command_for_key(self, key: str) -> Optional[str]:
        """Get command name for a key.

        Args:
            key: Raw key input from terminal (escape sequence)

        Returns:
            Command name if bound, None otherwise
        """
        # Key from terminal is already in escape sequence form
        return self.bindings.get(key)

    def get_keys_for_command(self, command_name: str) -> list[str]:
        """Get all keys bound to a command."""
        return [key for key, cmd in self.bindings.items() if cmd == command_name]

    def load_from_file(self, filepath: str):
        """Load keybindings from YAML or JSON file.

        Expected format (can use readable names or escape sequences):
            bindings:
              q: quit
              right_arrow: move_forward
              ctrl_s: save
              shift_up_arrow: move_forward_large
        """
        path = Path(filepath)

        if path.suffix in ['.yaml', '.yml']:
            with open(path) as f:
                config = yaml.safe_load(f)
        elif path.suffix == '.json':
            with open(path) as f:
                config = json.load(f)
        else:
            raise ValueError(f"Unsupported file format: {path.suffix}")

        # Load bindings (keys will be automatically normalized)
        for key, command_name in config.get('bindings', {}).items():
            self.bind(key, command_name)

    def save_to_file(self, filepath: str, use_readable_names: bool = True):
        """Save keybindings to YAML or JSON file.

        Args:
            filepath: Path to save to
            use_readable_names: If True, convert escape sequences to readable names
        """
        path = Path(filepath)

        # Convert bindings to readable format if requested
        if use_readable_names:
            readable_bindings = {
                get_readable_key(k): v for k, v in self.bindings.items()
            }
        else:
            readable_bindings = self.bindings.copy()

        config = {'bindings': readable_bindings}

        if path.suffix in ['.yaml', '.yml']:
            with open(path, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)
        elif path.suffix == '.json':
            with open(path, 'w') as f:
                json.dump(config, f, indent=2)
        else:
            raise ValueError(f"Unsupported file format: {path.suffix}")

    def generate_help_text(self) -> str:
        """Generate help text showing all keybindings organized by category."""
        lines = ["Keybindings:\n"]

        # Group commands by category
        categories = {}
        for key, cmd_name in self.bindings.items():
            command = self.registry.get_command(cmd_name)
            if command:
                if command.category not in categories:
                    categories[command.category] = []
                categories[command.category].append((key, command))

        # Format output
        for category, items in sorted(categories.items()):
            lines.append(f"\n{category.upper()}:")
            for key, command in sorted(items, key=lambda x: x[1].name):
                # Format key nicely
                key_display = self._format_key_display(key)
                lines.append(f"  {key_display:20} - {command.description}")

        return "\n".join(lines)

    def get_footer_text(self):
        
        parts = []
        
        for key, cmd_name in self.bindings.items():
            cmd = self.registry.get_command(cmd_name)
            # if cmd and cmd.category == "default":
            if cmd:
                rkey = self._format_key_display(key)
                parts.append(f"{rkey}: {cmd.description}")
        
        return f"Navigation: {" | ".join(parts)}"

    def _format_key_display(self, key: str) -> str:
        """Format key for display in help text."""
        readable = get_readable_key(key)

        # Add some nice formatting for common keys
        display_map = {
            'up_arrow': '↑',
            'down_arrow': '↓',
            'right_arrow': '→',
            'left_arrow': '←',
            'space': 'Space',
            'enter': 'Enter',
            'tab': 'Tab',
            'backspace': 'Backspace',
            'escape': 'Esc',
        }

        return display_map.get(readable, readable)


# Example default configuration (now with readable key names!)
DEFAULT_KEYBINDINGS = """
bindings:
  # General
  q: quit
  '?': show_help

  # Navigation - Arrow keys
  right_arrow: move_forward
  left_arrow: move_backward
  up_arrow: move_forward_large
  down_arrow: move_backward_large
  space: move_forward

  # Navigation - Vim style
  l: move_forward
  h: move_backward
  k: move_forward_large
  j: move_backward_large

  # Navigation - With modifiers
  shift_right_arrow: move_forward_second
  shift_left_arrow: move_backward_second
  shift_up_arrow: move_forward_large_second
  shift_down_arrow: move_backward_large_second
  ctrl_right_arrow: next_gene
  ctrl_left_arrow: prev_gene

  # Chromosome navigation
  c: switch_chromosome
  g: goto_position
  n: next_feature
  p: prev_feature

  # Display toggles
  f: toggle_features
  a: toggle_amino_acids
  r: toggle_rna
  '-': toggle_reverse_strand

  # Window controls
  w: change_window_size
  s: change_stride
  ctrl_s: save_keybindings

  # Info
  i: show_position_info
  ctrl_h: set_highlight
"""
