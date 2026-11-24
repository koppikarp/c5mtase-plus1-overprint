#!/usr/bin/env python3
"""
Step 3A: Publication-Grade Maximum-Likelihood Phylogeny (MAFFT L-INS-i +
VeryFastTree LG+Gamma with uniform thorough search flags across BOTH main
trees AND all 100 non-parametric ML bootstrap replicates, plus a long-branch
exclusion sensitivity check).

Runs 100% self-contained from repository data/phylogeny/ input FASTAs:
  - data/phylogeny/tree1_frame0_cat_input.faa (97 unique taxa)
  - data/phylogeny/tree2_cslsdile_domain_input.faa (67 harmonized taxa)

Outputs:
  - data/phylogeny/tree1_frame0_cat_trimmed.afa
  - data/phylogeny/tree1_frame0_cat_ml_lg.nwk
  - data/phylogeny/tree1_boot100.nwk
  - data/phylogeny/tree2_cslsdile_domain_trimmed.afa
  - data/phylogeny/tree2_cslsdile_domain_ml_lg.nwk
  - data/phylogeny/tree2_boot100.nwk
  - data/phylogeny/tree1_no_longbranch_ml_lg.nwk (+ tree1_no_longbranch_boot100.nwk)
"""
from __future__ import annotations

import os
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Tuple

import numpy as np

import lib_c5mtase as lib


VFT_THOROUGH_FLAGS = ["-lg", "-gamma", "-spr", "4", "-mlacc", "2", "-slownni"]


def _run_one_bootstrap(args: Tuple[int, List[Tuple[str, str]], np.ndarray, str]) -> str:
    b_idx, trimmed, bcols, tmpd = args
    bf = os.path.join(tmpd, f"boot_{b_idx}.afa")
    with open(bf, "w") as out:
        for h, s in trimmed:
            out.write(f">{h}\n{''.join(s[c] for c in bcols)}\n")
    res = subprocess.run(
        ["VeryFastTree", *VFT_THOROUGH_FLAGS, "-nosupport", "-threads", "1", "-quiet", bf],
        capture_output=True, text=True, check=True,
    )
    return res.stdout.strip()


def build_ml_tree_with_bootstraps(
    in_fa: Path,
    trimmed_out: Path,
    nwk_out: Path,
    boot_out: Path,
    n_boot: int = 100,
    seed: int = 42,
    exclude_taxa: set | None = None,
) -> None:
    t0 = time.time()
    records = lib.read_fasta(in_fa)
    if exclude_taxa:
        records = [r for r in records if r[0] not in exclude_taxa]

    with tempfile.TemporaryDirectory() as tmpd:
        raw_fa = os.path.join(tmpd, "raw.faa")
        aln_fa = os.path.join(tmpd, "mafft.afa")
        with open(raw_fa, "w") as out:
            for tid, _desc, seq in records:
                out.write(f">{tid}\n{seq}\n")

        # 1. MAFFT L-INS-i
        with open(aln_fa, "w") as out:
            subprocess.run(
                ["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8", "--quiet", raw_fa],
                stdout=out, check=True,
            )

        # 2. Trim columns with >=50% occupancy
        aln = lib.read_fasta(Path(aln_fa))
        L = len(aln[0][2])
        keep_cols = [c for c in range(L) if sum(1 for _, _, s in aln if s[c] != "-") / len(aln) >= 0.50]
        trimmed = [(tid, "".join(s[c] for c in keep_cols)) for tid, _d, s in aln]
        with open(trimmed_out, "w") as out:
            for tid, s in trimmed:
                out.write(f">{tid}\n{s}\n")
        print(f"  [{nwk_out.stem}] {len(trimmed)} taxa, {L} raw -> {len(keep_cols)} trimmed cols (>=50% occ)")

        # 3. Main ML tree with thorough flags
        res = subprocess.run(
            ["VeryFastTree", *VFT_THOROUGH_FLAGS, "-threads", "8", "-quiet", str(trimmed_out)],
            capture_output=True, text=True, check=True,
        )
        nwk_out.write_text(res.stdout.strip() + "\n")

        # 4. 100 ML bootstrap replicates with IDENTICAL thorough search flags!
        rng = np.random.default_rng(seed)
        Lt = len(keep_cols)
        tasks = [(b, trimmed, rng.integers(0, Lt, size=Lt), tmpd) for b in range(n_boot)]
        with ThreadPoolExecutor(max_workers=8) as ex:
            boot_nwks = list(ex.map(_run_one_bootstrap, tasks))
        boot_out.write_text("\n".join(boot_nwks) + "\n")
        print(f"  [{nwk_out.stem}] Main ML + {n_boot} thorough ML bootstraps completed in {time.time() - t0:.1f}s")


def main():
    pdir = lib.PHYLO_DIR
    print("Building Tree 1 (Frame-0 C5-MTase catalytic core, 97 unique taxa)...")
    build_ml_tree_with_bootstraps(
        pdir / "tree1_frame0_cat_input.faa",
        pdir / "tree1_frame0_cat_trimmed.afa",
        pdir / "tree1_frame0_cat_ml_lg.nwk",
        pdir / "tree1_boot100.nwk",
        n_boot=100,
        seed=42,
    )

    print("Building Tree 2 (CSLSDILE/KYFLS domain, 67 harmonized taxa)...")
    build_ml_tree_with_bootstraps(
        pdir / "tree2_cslsdile_domain_input.faa",
        pdir / "tree2_cslsdile_domain_trimmed.afa",
        pdir / "tree2_cslsdile_domain_ml_lg.nwk",
        pdir / "tree2_boot100.nwk",
        n_boot=100,
        seed=42,
    )

    print("Building Tree 1 sensitivity check (excluding extreme long branch A1_SPCT_Baci_10645, 96 taxa)...")
    build_ml_tree_with_bootstraps(
        pdir / "tree1_frame0_cat_input.faa",
        pdir / "tree1_no_longbranch_trimmed.afa",
        pdir / "tree1_no_longbranch_ml_lg.nwk",
        pdir / "tree1_no_longbranch_boot100.nwk",
        n_boot=100,
        seed=42,
        exclude_taxa={"A1_SPCT_Baci_10645"},
    )
    print("All ML trees and thorough 100-replicate bootstraps saved cleanly!")


if __name__ == "__main__":
    main()
