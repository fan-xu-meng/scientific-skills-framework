#!/usr/bin/env python3
"""Create TEM DM overviews and optionally Pt/C particle and FFT candidates."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys


PT_SCRIPT_ORDER = (
    "analyze_tem_particle_facet.py",
    "match_tem_fft_to_pt_pdf_card.py",
    "summarize_tem_facet_proportions.py",
    "summarize_pdf_plane_proportions.py",
)
REQUIRED_MODULES = ("ncempy", "numpy", "pandas", "scipy", "skimage", "matplotlib")


def script_command(script_path: Path) -> list[str]:
    if all(importlib.util.find_spec(module) is not None for module in REQUIRED_MODULES):
        return [sys.executable, str(script_path)]
    uv = shutil.which("uv")
    if uv is None:
        missing = [module for module in REQUIRED_MODULES if importlib.util.find_spec(module) is None]
        raise RuntimeError(
            "TEM Python dependencies are missing and uv is unavailable: " + ", ".join(missing)
        )
    return [
        uv,
        "run",
        "--quiet",
        "--python",
        "3.13",
        "--with",
        "ncempy",
        "--with",
        "scikit-image",
        "--with",
        "scipy",
        "--with",
        "pandas",
        "--with",
        "matplotlib",
        "python",
        str(script_path),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path, help="Folder containing DM3/DM4 images")
    parser.add_argument("--out-dir", required=True, type=Path, help="Separate output folder")
    parser.add_argument("--pt-candidates", action="store_true", help="Also run Pt/C-specific particle and FFT candidate checks")
    args = parser.parse_args()
    data_root = args.data_dir.expanduser().resolve()
    if not data_root.is_dir():
        raise FileNotFoundError(f"TEM input directory not found: {data_root}")
    if not any(data_root.rglob("*.dm3")) and not any(data_root.rglob("*.dm4")):
        raise FileNotFoundError(f"No DM3/DM4 files found under: {data_root}")

    environment = os.environ.copy()
    environment["TEM_DATA_DIR"] = str(data_root)
    environment["TEM_OUTPUT_DIR"] = str(args.out_dir.expanduser().resolve())
    script_root = Path(__file__).resolve().parent
    script_order = ("analyze_tem_dm.py",) + (PT_SCRIPT_ORDER if args.pt_candidates else ())
    for script_name in script_order:
        subprocess.run(
            script_command(script_root / script_name),
            cwd=data_root,
            env=environment,
            check=True,
        )


if __name__ == "__main__":
    main()
