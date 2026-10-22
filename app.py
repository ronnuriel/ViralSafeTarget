"""Development entry point for the packaged ViralSafeTarget researcher interface."""

import sys
from pathlib import Path

SOURCE = Path(__file__).resolve().parent / "src"
if SOURCE.is_dir():  # Prefer the checkout when this development entry point is used.
    sys.path.insert(0, str(SOURCE))

def main() -> None:
    from viral_safe_target.ui_app import render_app

    render_app()


main()
