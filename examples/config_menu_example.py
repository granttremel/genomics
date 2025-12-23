"""
Example showing how to use the parameter configuration menu system.
"""

from ggene.browser.base_browser import BaseBrowser, BaseBrowserState
from dataclasses import dataclass


@dataclass
class MyBrowserState(BaseBrowserState):
    """Extended state with configurable parameters."""
    position: int = 0
    chromosome: str = "1"
    window_size: int = 100
    stride: int = 20
    show_features: bool = True
    show_quality: bool = False


class ConfigurableArtist:
    """Example artist with configurable parameters."""

    def __init__(self):
        # Visual parameters
        self.color_scheme = "default"
        self.show_grid = True
        self.line_height = 2

        # Rendering parameters
        self.max_features = 100
        self.feature_label_size = 20

    def render(self, state, window):
        # Rendering logic here
        return [f"Rendered with {self.color_scheme} colors"]


class MyGenomeBrowser(BaseBrowser):
    """Example browser with parameter menu."""

    def __init__(self, *args, **kwargs):
        # Initialize components that have parameters
        
        artist = ConfigurableArtist()
        state = MyBrowserState()
        
        # Call parent init (which will call register_parameters)
        super().__init__(*args, artist = artist, state = state, **kwargs)

    def register_parameters(self):
        """Register all configurable parameters."""

        # Register state parameters
        self.params.register(
            "display.window_size",
            current_value=self.state.window_size,
            param_type=int,
            description="Number of bases to display",
            category="display",
            min_value=50,
            max_value=1000,
            on_change=lambda v: setattr(self.state, 'window_size', v)
        )

        self.params.register(
            "navigation.stride",
            current_value=self.state.stride,
            param_type=int,
            description="Bases to move on arrow press",
            category="navigation",
            min_value=1,
            max_value=100,
            on_change=lambda v: setattr(self.state, 'stride', v)
        )

        self.params.register(
            "display.show_features",
            current_value=self.state.show_features,
            param_type=bool,
            description="Show genomic features",
            category="display",
            on_change=lambda v: setattr(self.state, 'show_features', v)
        )

        # Register artist parameters
        self.params.register(
            "artist.color_scheme",
            current_value=self.renderer.color_scheme,
            param_type=str,
            description="Color scheme for visualization",
            category="appearance",
            choices=["default", "dark", "light", "colorblind"],
            on_change=lambda v: setattr(self.renderer, 'color_scheme', v)
        )

        self.params.register(
            "artist.show_grid",
            current_value=self.renderer.show_grid,
            param_type=bool,
            description="Show grid lines",
            category="appearance",
            on_change=lambda v: setattr(self.renderer, 'show_grid', v)
        )

        self.params.register(
            "artist.line_height",
            current_value=self.renderer.line_height,
            param_type=int,
            description="Height of each line in characters",
            category="appearance",
            min_value=1,
            max_value=5,
            on_change=lambda v: setattr(self.renderer, 'line_height', v)
        )

        self.params.register(
            "artist.max_features",
            current_value=self.renderer.max_features,
            param_type=int,
            description="Maximum features to display",
            category="performance",
            min_value=10,
            max_value=1000,
            on_change=lambda v: setattr(self.renderer, 'max_features', v)
        )

    def reload_config(self):
        """Apply configuration changes."""
        # The on_change callbacks already updated the objects
        # But we might need to reinitialize some components

        # For example, if window size changed, update iterator
        if self.iterator:
            self.iterator.window_size = self.state.window_size

        # Re-render with new parameters
        if self.window:
            self.render()

        print("Configuration reloaded!")


# Alternative: Auto-register parameters from object attributes
class AutoConfigBrowser(BaseBrowser):
    """Browser that auto-discovers parameters."""

    def __init__(self, *args, **kwargs):
        self.renderer = ConfigurableArtist()
        super().__init__(*args, **kwargs)

    def register_parameters(self):
        """Auto-register parameters from artist."""
        from ggene.browser.config_menu import auto_register_params

        # Automatically discover and register all simple attributes
        auto_register_params(
            self.renderer,
            self.params,
            path_prefix="artist",
            category="artist_config"
        )

        # Can still manually register special ones
        self.params.register(
            "display.window_size",
            current_value=getattr(self.state, 'window_size', 100),
            description="Window size in bases",
            category="display",
            min_value=50,
            max_value=1000
        )


# Example usage demonstrating the menu workflow
def example_usage():
    """
    Demo of the configuration menu workflow:

    1. User presses Ctrl+P to open menu
    2. Menu shows categories (display, navigation, appearance, performance)
    3. User navigates with arrows, Enter to select category
    4. Menu shows parameters in that category
    5. User selects parameter to edit
    6. Menu prompts for new value with validation
    7. Value is updated via on_change callback
    8. User exits menu (Esc)
    9. Browser calls reload_config() to apply changes
    10. Browser continues with new settings
    """

    browser = MyGenomeBrowser()

    # Programmatically access and modify parameters
    window_size_param = browser.params.get("display.window_size")
    print(f"Current window size: {window_size_param.current_value}")

    # Set value programmatically (also triggers on_change)
    success, msg = browser.params.set_value("display.window_size", 200)
    print(f"Update: {msg}")

    # Get all parameters in a category
    display_params = browser.params.get_by_category("display")
    print(f"Display parameters: {[p.name for p in display_params]}")

    # Open the interactive menu (in real browser, bound to Ctrl+P)
    # browser.open_config_menu()


if __name__ == "__main__":
    example_usage()
