"""
Example showing how to use the command/keybinding system in your browser.
"""

from ggene.browser.base_browser import BaseBrowser, BaseBrowserState
from dataclasses import dataclass


@dataclass
class MyBrowserState(BaseBrowserState):
    """Extended state for your browser."""
    position: int = 0
    chromosome: str = "1"
    window_size: int = 100


class MyGenomeBrowser(BaseBrowser):
    """Example genome browser using the command system."""

    def register_commands(self):
        """Register all available commands."""
        # Call parent to get base commands (quit, help)
        super().register_commands()

        # Register navigation commands
        @self.registry.register(
            "move_forward",
            "Move forward by stride",
            category="navigation"
        )
        def _move_forward(state, window):
            state.position += 100
            return state

        @self.registry.register(
            "move_backward",
            "Move backward by stride",
            category="navigation"
        )
        def _move_backward(state, window):
            state.position = max(0, state.position - 100)
            return state

        @self.registry.register(
            "move_forward_large",
            "Move forward by 10x stride",
            category="navigation"
        )
        def _move_forward_large(state, window):
            state.position += 1000
            return state

        @self.registry.register(
            "move_backward_large",
            "Move backward by 10x stride",
            category="navigation"
        )
        def _move_backward_large(state, window):
            state.position = max(0, state.position - 1000)
            return state

        # Register chromosome commands
        @self.registry.register(
            "switch_chromosome",
            "Switch to a different chromosome",
            category="chromosome",
            takes_input=True
        )
        def _switch_chromosome(state, window):
            # This would prompt for input
            new_chr = input("Enter chromosome: ")
            state.chromosome = new_chr
            return state

        @self.registry.register(
            "goto_position",
            "Jump to a specific position",
            category="navigation",
            takes_input=True
        )
        def _goto_position(state, window):
            new_pos = input("Enter position: ")
            try:
                state.position = int(new_pos)
            except ValueError:
                print("Invalid position")
            return state

        # Save keybindings command
        @self.registry.register(
            "save_keybindings",
            "Save current keybindings to file",
            category="general"
        )
        def _save_keybindings(state, window):
            self.keybindings.save_to_file("my_keybindings.yaml")
            print("Keybindings saved to my_keybindings.yaml")
            return state

    def setup_default_keybindings(self):
        """Setup default keybindings."""
        super().setup_default_keybindings()

        # Navigation - can use readable names!
        self.keybindings.bind('right_arrow', 'move_forward')
        self.keybindings.bind('left_arrow', 'move_backward')
        self.keybindings.bind('up_arrow', 'move_forward_large')
        self.keybindings.bind('down_arrow', 'move_backward_large')

        # Also bind vim keys to the same commands
        self.keybindings.bind('l', 'move_forward')
        self.keybindings.bind('h', 'move_backward')
        self.keybindings.bind('k', 'move_forward_large')
        self.keybindings.bind('j', 'move_backward_large')

        # Chromosome navigation
        self.keybindings.bind('c', 'switch_chromosome')
        self.keybindings.bind('g', 'goto_position')

        # Modifier key combinations
        self.keybindings.bind('ctrl_s', 'save_keybindings')


# Example usage:
if __name__ == "__main__":
    # Option 1: Use default keybindings
    browser = MyGenomeBrowser()

    # Option 2: Load custom keybindings from file
    browser_custom = MyGenomeBrowser(
        keybinding_file="ggene/browser/default_keybindings.yaml"
    )

    # Option 3: Programmatically set keybindings
    browser3 = MyGenomeBrowser()
    browser3.keybindings.bind('ctrl_right_arrow', 'move_forward_large')
    browser3.keybindings.bind('shift_up_arrow', 'move_forward')

    # Save your keybindings for reuse
    browser3.keybindings.save_to_file("my_custom_keys.yaml")

    # Show help (displays all keybindings organized by category)
    print(browser3.keybindings.generate_help_text())
