"""Launch the optional Streamlit researcher interface."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def launch_ui(*, port: int = 8501, headless: bool = False) -> int:
    """Run the packaged Streamlit app and return its process exit code."""
    try:
        import streamlit  # noqa: F401
    except ImportError as error:
        raise SystemExit(
            'The interface requires the optional UI dependency. Install it with: '
            'pip install "viral-safe-target[ui]"'
        ) from error

    app = Path(__file__).with_name("ui_app.py")
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app),
        "--server.port",
        str(port),
        "--server.headless",
        "true" if headless else "false",
    ]
    return subprocess.call(command)


def main() -> None:
    """Console-script entry point."""
    raise SystemExit(launch_ui())
