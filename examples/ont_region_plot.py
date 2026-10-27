#!/usr/bin/env python3
"""Compatibility wrapper for the platform-neutral ``region_plot.py`` example."""

from pathlib import Path

from region_plot import main

if __name__ == "__main__":
    main(Path("ont_region.png"))
