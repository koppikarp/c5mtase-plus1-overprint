#!/usr/bin/env python3
"""
Step 1: Pan-2,339-genome GTDB R232 C5-MTase +1 dual-coding census audit.

Default mode (100% self-contained): verifies and summarizes the bundled
data/gtdb2339_c5mtase_plus1_census.json and data/uniref100_deletion_audit.json.

Optional raw re-scan mode:
  python scripts/01_gtdb_census_audit.py --genome-dir /path/to/gtdb_fna --metadata-tsv /path/to/panel2339_summary.tsv --seeds-faa /path/to/seeds.faa
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lib_c5mtase as lib

MOTIF_IV_RE = re.compile(r"[GASCVITL][FTSAVILYCM]PC[QSN]")
MOTIF_VI_RE = re.compile(r"EN[VLIM]")
MOTIF_IX_RE = re.compile(r"[RK][LMIV][TSIVL][PRK].{1,3}[CEWYFL].[RK][LM][QEH][GSA][FYW][PRK]")
MOTIF_X_RE = re.compile(r"GN[ASGT][VMI][TAVSIC]")


def verify_bundled_census() -> None:
    gtdb = json.loads((lib.DATA_DIR / "gtdb2339_c5mtase_plus1_census.json").read_text())
    all_loci = gtdb["overprinted"] + gtdb["inframe_cterm"] + gtdb["other"]
    n_loci = len(all_loci)
    n_genomes = len({x["acc"] for x in all_loci})
    n_orders = len({(x.get("phylum"), x.get("order")) for x in all_loci})
    phyla = Counter(x["phylum"] for x in all_loci)
    arch = Counter(x.get("refined_arch", "unknown") for x in all_loci)

    print("=== GTDB R232 2,339-GENOME CENSUS VERIFICATION ===")
    print(f"  Total loci             : {n_loci}")
    print(f"  Unique genomes         : {n_genomes}")
    print(f"  Unique orders/families : {n_orders}")
    print(f"  Unique phyla           : {len(phyla)} ({sum(1 for p in phyla if p not in ('p__Nanobdellota','p__Aenigmatarchaeota'))} bacterial + {sum(1 for p in phyla if p in ('p__Nanobdellota','p__Aenigmatarchaeota'))} DPANN archaeal)")
    print(f"  Refined architectures  : {dict(arch)}")
    print("  Phylum breakdown:")
    for ph, cnt in phyla.most_common():
        sub = [x for x in all_loci if x["phylum"] == ph]
        ac = Counter(x.get("refined_arch") for x in sub)
        print(f"    {ph:24s}: {cnt:2d}  {dict(ac)}")

    uniref_path = lib.DATA_DIR / "uniref100_deletion_audit.json"
    if uniref_path.exists():
        u = json.loads(uniref_path.read_text())
        print(f"\n=== UNIREF100 / UNIPROTKB DELETION AUDIT ===")
        print(f"  Audited UniProtKB hits : {u['n_audited']} | Deleted: {u['n_deleted']} ({u['pct_deleted']}%) | Active: {u['n_active']}")
        print(f"  Pfam / InterPro hits   : {u['pfam_hits']} / {u['interpro_hits']} | PDB/CATH Foldseek E<1: {u['pdb_cath_foldseek_hits_e_lt_1']}")


def _scan_one_genome(args):
    acc, gc, use_t4, genome_dir = args
    fna = Path(genome_dir) / f"{acc}.fna.gz"
    if not fna.exists():
        return acc, []
    tbl = lib.CODON_TABLE.copy()
    if use_t4:
        tbl["TGA"] = "W"
    contigs = lib.read_fna_gz(fna)
    mtases = []
    for cid, seq in contigs.items():
        if len(seq) < 450:
            continue
        for st, s in ((1, seq), (-1, lib.revcomp(seq))):
            for fr in range(3):
                prot = lib.translate(s, tbl, offset=fr)
                pos = 0
                for seg in prot.split("*"):
                    sl = len(seg)
                    if 210 <= sl <= 850:
                        m4 = MOTIF_IV_RE.search(seg)
                        if m4 and MOTIF_VI_RE.search(seg[m4.start():m4.start() + 95]):
                            s_nt = fr + pos * 3
                            e_nt = s_nt + sl * 3
                            fnd = -1
                            for j in range(min(m4.start(), int(sl * 0.35) + 1)):
                                if s[s_nt + j * 3:s_nt + j * 3 + 3] in lib.CANON_STARTS:
                                    fnd = j
                                    break
                            if fnd != -1:
                                hprot = seg[fnd:]
                                hnt = s[s_nt + fnd * 3:e_nt]
                                pcq = m4.start() - fnd + 1
                                if pcq >= 25 and len(hprot) - pcq >= 120:
                                    p1 = lib.translate(hnt, tbl, offset=1)
                                    best_seg, best_s, best_e = "", 0, 0
                                    ppos = 0
                                    for pseg in p1.split("*"):
                                        ps, pe = ppos, ppos + len(pseg)
                                        if ps <= pcq + 95 and pe >= pcq + 70 and len(pseg) > len(best_seg):
                                            best_seg, best_s, best_e = pseg, ps, pe
                                        ppos = pe + 1
                                    if len(best_seg) >= 45:
                                        mtases.append((acc, cid, st, hprot, hnt, pcq, best_s, best_e, best_seg))
                    pos += sl + 1
    return acc, mtases


def main():
    ap = argparse.ArgumentParser(description="GTDB R232 2,339-genome C5-MTase +1 census audit")
    ap.add_argument("--genome-dir", type=str, default=None, help="Optional directory of raw GTDB .fna.gz assemblies for full re-scan")
    ap.add_argument("--metadata-tsv", type=str, default=None, help="Optional GTDB panel metadata TSV")
    ap.add_argument("--seeds-faa", type=str, default=None, help="Optional CSLSDILE/KYFLS seed FASTA for DIAMOND blastp")
    args = ap.parse_args()

    if not args.genome_dir:
        verify_bundled_census()
        return

    rows = list(csv.DictReader(open(args.metadata_tsv), delimiter="\t"))
    tasks = [(r["acc"], float(r["gc"]), float(r["f_TGA"]) < 0.05, args.genome_dir) for r in rows]
    all_mt = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for _acc, mts in ex.map(_scan_one_genome, tasks):
            all_mt.extend(mts)
    print(f"Raw GTDB scan finished: {len(tasks)} genomes -> {len(all_mt)} candidate C5-MTase +1 spans.")
    if args.seeds_faa:
        with tempfile.TemporaryDirectory() as td:
            sdb = os.path.join(td, "seeds.dmnd")
            qfa = os.path.join(td, "q.faa")
            out = os.path.join(td, "out.tsv")
            subprocess.run(["diamond", "makedb", "--in", args.seeds_faa, "-d", sdb, "--quiet"], check=True)
            with open(qfa, "w") as f:
                for i, m in enumerate(all_mt):
                    f.write(f">M_{i}\n{m[8]}\n")
            subprocess.run(["diamond", "blastp", "-d", sdb, "-q", qfa, "-o", out, "--evalue", "1e-5", "--sensitive", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "--quiet"], check=True)
            hits = sum(1 for _ in open(out))
            print(f"  DIAMOND blastp E<1e-5 hits against seeds: {hits}")


if __name__ == "__main__":
    main()
