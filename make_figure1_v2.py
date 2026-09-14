"""Compatibility entry point for Giorgio's Figure 1, now in the production workflow.

Accepts the original published workbook or its cleaned equivalent. The exact supplied
script is preserved for provenance at qc/figure1_bootstrap/giorgio_original.py.
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import data
import run_analysis
import make_figures

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", nargs="?", type=Path, default=data.WORKBOOK)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    data.WORKBOOK = args.workbook.resolve()
    run_analysis.RESULTS = make_figures.RESULTS = args.output_dir / "results"
    make_figures.FIGURES = args.output_dir / "figures"
    run_analysis.main()
