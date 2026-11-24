#!/usr/bin/env python3
"""
Master reproducible pipeline runner for the C5-MTase +1 dual-coding study.

Runs all self-contained post-census analyses directly from repository data/
(ZERO hardcoded statistics, ZERO external machine paths) and regenerates:
  1. data/selection_4d_and_stop_stats.json     (via 02_selection_4d_and_neighborhood.py)
  2. data/phylogeny/ml_bootstrap_summary.json  (via 06_analyze_ml_trees.py)
  3. data/structures/trans_vs_cis_3d_audit.json (via 07_compare_trans_vs_cis_3d.py)
  4. data/phage_p1_dmt_validation.json         (via 09_inspect_phage_p1_dmt.py)
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent

PIPELINE = [
    ("01_gtdb_census_audit.py", "GTDB R232 census verification"),
    ("02_selection_4d_and_neighborhood.py", "4D synonymous conservation, 5-frame stop purge & Test 4"),
    ("03_scan_inphared_phages.py", "INPHARED phage census & dereplicated prevalence audit"),
    ("04_verify_uniprot_dna.py", "UniProt TPCQ ENA DNA verification audit"),
    ("06_analyze_ml_trees.py", "ML phylogeny bootstrap & non-overlapping Mantel congruence"),
    ("07_compare_trans_vs_cis_3d.py", "Rigid-core 1MHT & 30-model AlphaFold 3 structural audit"),
    ("08_audit_af3_models.py", "AlphaFold 3 per-model confidence & Motif-IX lock summary"),
    ("09_inspect_phage_p1_dmt.py", "Bacteriophage P1 dmt coordinate & mutation verification"),
]


def main():
    for script, desc in PIPELINE:
        print(f"\n{'=' * 90}\n>>> Running {script}: {desc}\n{'=' * 90}")
        subprocess.run([sys.executable, str(SCRIPTS_DIR / script)], check=True)
    print("\nAll reproducible data exports and verification audits completed with exit code 0!")


if __name__ == "__main__":
    main()
