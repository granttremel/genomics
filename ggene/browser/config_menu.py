"""
Interactive configuration menu system for genome browser.

Allows runtime parameter editing with a navigable menu interface.
"""

from typing import Any, Callable, Dict, List, Optional, Union, get_type_hints
from dataclasses import dataclass, field
import sys
import tty
import termios
import os
from collections import defaultdict


@dataclass
class ConfigParam:
    """Metadata for a configurable parameter."""
    name: str
    path: str  # Full path like "renderer.artist.color_scheme"
    current_value: Any
    param_type: type
    description: str = ""
    category: str = "general"
    min_value: Optional[Any] = None
    max_value: Optional[Any] = None
    choices: Optional[List[Any]] = None  # For enum-like params
    on_change: Optional[Callable] = None  # Callback when value changes

    def validate(self, value: Any) -> tuple[bool, str]:
        """Validate a new value.

        Returns:
            (is_valid, error_message)
        """
        # Type check
        try:
            if self.param_type == bool:
                # Accept various boolean representations
                if isinstance(value, str):
                    value = value.lower() in ('true', 'yes', '1', 'on')
                value = bool(value)
            elif self.param_type in (int, float, str):
                value = self.param_type(value)
        except (ValueError, TypeError) as e:
            return False, f"Cannot convert to {self.param_type.__name__}: {e}"

        # Range check
        if self.min_value is not None and value < self.min_value:
            return False, f"Value must be >= {self.min_value}"
        if self.max_value is not None and value > self.max_value:
            return False, f"Value must be <= {self.max_value}"

        # Choices check
        if self.choices is not None and value not in self.choices:
            return False, f"Must be one of: {', '.join(map(str, self.choices))}"

        return True, ""

    def format_value(self) -> str:
        """Format current value for display."""
        if isinstance(self.current_value, bool):
            return "✓" if self.current_value else "✗"
        return str(self.current_value)


class ParameterRegistry:
    """Registry of all configurable parameters in the browser."""

    def __init__(self):
        self.params: Dict[str, ConfigParam] = {}
        self._categories: Dict[str, List[str]] = defaultdict(list)

    def register(
        self,
        path: str,
        current_value: Any,
        param_type: type = None,
        description: str = "",
        category: str = "general",
        min_value: Optional[Any] = None,
        max_value: Optional[Any] = None,
        choices: Optional[List[Any]] = None,
        on_change: Optional[Callable] = None
    ):
        """Register a parameter.

        Args:
            path: Dot-separated path (e.g., "renderer.artist.color_scheme")
            current_value: Current value
            param_type: Type of parameter (inferred if not provided)
            description: Human-readable description
            category: Category for organization
            min_value: Minimum allowed value (for numeric types)
            max_value: Maximum allowed value (for numeric types)
            choices: List of valid choices (for enum-like params)
            on_change: Callback function(new_value) called when value changes
        """
        if param_type is None:
            param_type = type(current_value)

        name = path.split('.')[-1]

        param = ConfigParam(
            name=name,
            path=path,
            current_value=current_value,
            param_type=param_type,
            description=description,
            category=category,
            min_value=min_value,
            max_value=max_value,
            choices=choices,
            on_change=on_change
        )

        self.params[path] = param
        self._categories[category].append(path)

    def get(self, path: str) -> Optional[ConfigParam]:
        """Get parameter by path."""
        return self.params.get(path)

    def set_value(self, path: str, value: Any) -> tuple[bool, str]:
        """Set parameter value with validation.

        Returns:
            (success, message)
        """
        param = self.params.get(path)
        if not param:
            return False, f"Unknown parameter: {path}"

        is_valid, error = param.validate(value)
        if not is_valid:
            return False, error

        # Convert to proper type
        if param.param_type == bool and isinstance(value, str):
            value = value.lower() in ('true', 'yes', '1', 'on')
        else:
            value = param.param_type(value)

        param.current_value = value

        # Call on_change callback if provided
        if param.on_change:
            param.on_change(value)

        return True, "Value updated"

    def get_by_category(self, category: str) -> List[ConfigParam]:
        """Get all parameters in a category."""
        return [self.params[path] for path in self._categories.get(category, [])]

    def get_categories(self) -> List[str]:
        """Get list of all categories."""
        return sorted(self._categories.keys())

    def get_tree_structure(self) -> Dict[str, Any]:
        """Get parameters organized as a tree structure.

        Returns nested dict representing the object tree.
        """
        tree = {}
        for path in self.params:
            parts = path.split('.')
            current = tree
            for part in parts[:-1]:
                if part not in current:
                    current[part] = {}
                current = current[part]
            current[parts[-1]] = self.params[path]
        return tree


class ConfigMenu:
    """Interactive menu for editing configuration parameters."""

    def __init__(self, registry: ParameterRegistry):
        self.registry = registry
        self.selected_index = 0
        self.current_category = None
        self.mode = "category"  # "category", "param_list", or "edit"
        self.editing_param = None

    def show(self):
        """Display and run the interactive menu."""
        # Save terminal settings
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)

        try:
            tty.setraw(fd)
            self._run_menu()
        finally:
            # Restore terminal settings
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)

    def _run_menu(self):
        """Main menu loop."""
        while True:
            self._display()
            key = self._get_key()

            if key == 'q' and self.mode == "category":
                break
            elif key == '\x1b':  # Escape
                if self.mode == "edit":
                    self.mode = "param_list"
                elif self.mode == "param_list":
                    self.mode = "category"
                elif self.mode == "category":
                    break
            elif key == '\x1b[A':  # Up arrow
                self._move_selection(-1)
            elif key == '\x1b[B':  # Down arrow
                self._move_selection(1)
            elif key in ('\n', '\r'):  # Enter
                self._handle_select()

    def _display(self):
        """Display current menu state."""
        os.system('clear' if os.name == 'posix' else 'cls')

        print("┌─────────────────────────────────────────────────────────────┐")
        print("│         Configuration Menu                                  │")
        print("└─────────────────────────────────────────────────────────────┘\n")

        if self.mode == "category":
            self._display_categories()
        elif self.mode == "param_list":
            self._display_parameters()
        elif self.mode == "edit":
            self._display_editor()

        print("\n" + "─" * 64)
        print("Navigation: ↑/↓ = Move | Enter = Select | Esc/q = Back/Quit")

    def _display_categories(self):
        """Display category list."""
        print("Select a category:\n")
        categories = self.registry.get_categories()

        for i, category in enumerate(categories):
            params = self.registry.get_by_category(category)
            prefix = "→ " if i == self.selected_index else "  "
            print(f"{prefix}{category.title()} ({len(params)} parameters)")

    def _display_parameters(self):
        """Display parameters in current category."""
        if not self.current_category:
            return

        print(f"Category: {self.current_category.title()}\n")
        params = self.registry.get_by_category(self.current_category)

        for i, param in enumerate(params):
            prefix = "→ " if i == self.selected_index else "  "
            value_str = param.format_value()

            # Truncate description if too long
            desc = param.description[:40] + "..." if len(param.description) > 40 else param.description

            print(f"{prefix}{param.name:25} {value_str:10} {desc}")

    def _display_editor(self):
        """Display parameter editor."""
        if not self.editing_param:
            return

        param = self.editing_param

        print(f"Editing: {param.path}\n")
        print(f"Description: {param.description}")
        print(f"Type: {param.param_type.__name__}")
        print(f"Current value: {param.current_value}\n")

        if param.choices:
            print(f"Valid choices: {', '.join(map(str, param.choices))}")
        if param.min_value is not None or param.max_value is not None:
            print(f"Range: {param.min_value} - {param.max_value}")

        print("\n" + "─" * 64)

        # Get new value from user
        # (In a real implementation, this would be more sophisticated)
        termios.tcsetattr(sys.stdin.fileno(), termios.TCSADRAIN, self._old_term_settings)

        if param.param_type == bool:
            new_value = input("\nNew value (true/false): ")
        elif param.choices:
            print("\nChoices:")
            for i, choice in enumerate(param.choices):
                print(f"  {i+1}. {choice}")
            choice_idx = input("\nSelect number: ")
            try:
                new_value = param.choices[int(choice_idx) - 1]
            except (ValueError, IndexError):
                print("Invalid choice")
                input("Press Enter to continue...")
                self.mode = "param_list"
                return
        else:
            new_value = input(f"\nNew value: ")

        # Validate and set
        success, message = self.registry.set_value(param.path, new_value)
        print(f"\n{message}")
        input("Press Enter to continue...")

        # Go back to parameter list
        self.mode = "param_list"

        # Reset terminal
        tty.setraw(sys.stdin.fileno())

    def _move_selection(self, delta: int):
        """Move selection up or down."""
        if self.mode == "category":
            max_items = len(self.registry.get_categories())
        elif self.mode == "param_list":
            max_items = len(self.registry.get_by_category(self.current_category))
        else:
            return

        self.selected_index = (self.selected_index + delta) % max_items

    def _handle_select(self):
        """Handle Enter key press."""
        if self.mode == "category":
            # Select category
            categories = self.registry.get_categories()
            self.current_category = categories[self.selected_index]
            self.mode = "param_list"
            self.selected_index = 0

        elif self.mode == "param_list":
            # Select parameter to edit
            params = self.registry.get_by_category(self.current_category)
            self.editing_param = params[self.selected_index]

            # Save terminal settings before switching to edit mode
            self._old_term_settings = termios.tcgetattr(sys.stdin.fileno())

            self.mode = "edit"

    def _get_key(self) -> str:
        """Read a single key from stdin."""
        key = sys.stdin.read(1)

        # Handle escape sequences (arrow keys, etc.)
        if key == '\x1b':
            # Read the next two characters for arrow keys
            next1 = sys.stdin.read(1)
            if next1 == '[':
                next2 = sys.stdin.read(1)
                return '\x1b[' + next2
            return '\x1b'

        return key


def auto_register_params(obj: Any, registry: ParameterRegistry,
                        path_prefix: str = "", category: str = "general"):
    """Automatically register parameters from an object's attributes.

    Args:
        obj: Object to scan for parameters
        registry: ParameterRegistry to register to
        path_prefix: Path prefix for nested objects
        category: Category to register under
    """
    # Get public attributes that aren't methods
    for attr_name in dir(obj):
        if attr_name.startswith('_'):
            continue

        try:
            attr_value = getattr(obj, attr_name)

            # Skip methods
            if callable(attr_value):
                continue

            # Skip complex objects (could recurse here if needed)
            if not isinstance(attr_value, (int, float, str, bool, list, tuple)):
                continue

            path = f"{path_prefix}.{attr_name}" if path_prefix else attr_name

            # Register with a setter callback
            def make_setter(obj, attr):
                def setter(value):
                    setattr(obj, attr, value)
                return setter

            registry.register(
                path=path,
                current_value=attr_value,
                description=f"{attr_name.replace('_', ' ').title()}",
                category=category,
                on_change=make_setter(obj, attr_name)
            )
        except Exception:
            # Skip attributes that can't be accessed
            continue
