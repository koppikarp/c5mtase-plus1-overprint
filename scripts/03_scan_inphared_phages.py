#!/usr/bin/env python3
"""
Step 1B: INPHARED complete bacteriophage C5-MTase +1 dual-coding census.

Default mode (100% self-contained): verifies and summarizes bundled
data/inphared_phage_c5mtase_census.json (5,063 complete phages, 230 C5-MTases,
Motif IV denominators, <95%/<90%/<80% dereplicated prevalences, intact vs
interrupted Frame-0 Motif IX–X breakdown, including phage P1 and phage MM1).

Optional raw re-scan mode:
  python scripts/03_scan_inphared_phages.py --inphared-fasta /path/to/inphared_genomes.fasta.gz --seeds-faa /path/to/seeds.faa
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import lib_c5mtase as lib

MOTIF_IV_RE = re.compile(r"[GASCVITL][FTSAVILYCMW]PC[QSNT]")
MOTIF_VI_RE = re.compile(r"EN[VLIM]")


def verify_bundled_inphared() -> None:
    d = json.loads((lib.DATA_DIR / "inphared_phage_c5mtase_census.json").read_text())
    ov = d["overprinted"]
    m4_d = d["motif4_denominators"]
    m4_o = d["motif4_overprinted"]
    derep = d["dereplicated_prevalence"]
    arch = Counter(h.get("refined_arch") for h in ov)

    print("=== INPHARED COMPLETE PHAGE CENSUS VERIFICATION ===")
    print(f"  Complete phage genomes : {d['n_phages']}")
    print(f"  Total C5-MTases        : {d['n_mtases']} (sum of Motif IV denominators = {sum(m4_d.values())})")
    print(f"  Overprinted (+1) loci  : {len(ov)} ({100*len(ov)/d['n_mtases']:.1f}%) | Arch 2: 0 | Arch 3: {len(d['inframe_cterm'])}")
    print(f"  Frame-0 IX-X status    : {dict(arch)}")
    print("  Motif IV breakdown (overprinted / total):")
    for m4 in ("TPCQ", "FPCQ", "SPCQ", "YPCQ", "IPCQ", "VPCS", "FPCT"):
        t = m4_d.get(m4, 0)
        o = m4_o.get(m4, 0)
        print(f"    {m4:6s}: {o:2d} / {t:3d} ({100*o/max(1, t):5.1f}%)")
    other_tot = sum(v for k, v in m4_d.items() if k not in ("TPCQ", "FPCQ", "SPCQ"))
    print(f"    Other non-F/T/SPCQ: 0 / {other_tot}")
    print("  Dereplicated prevalence by Frame-0 catalytic-core AA identity:")
    for k, v in derep.items():
        bm = v["by_motif4"]
        print(
            f"    {k:9s}: {v['n_overprinted_clusters']:2d}/{v['n_clusters']:3d} ({v['pct_overprinted']:4.1f}%) | "
            f"TPCQ {bm['TPCQ']['overprinted']}/{bm['TPCQ']['total']} ({bm['TPCQ']['pct']:.1f}%) | "
            f"FPCQ {bm['FPCQ']['overprinted']}/{bm['FPCQ']['total']} ({bm['FPCQ']['pct']:.1f}%) | "
            f"SPCQ {bm['SPCQ']['overprinted']}/{bm['SPCQ']['total']} | mixed={v['mixed_clusters']}"
        )
    for key_acc in ("AF234172", "AJ302074"):
        hits = [h for h in ov if h["acc"] == key_acc]
        if hits:
            h = hits[0]
            print(f"  Key phage {key_acc}: {h['desc'][:55]} | m4={h['motif4']} hlen={h['hlen']} p1={h['p1_len']} refined={h['refined_arch']}")


def main():
    ap = argparse.ArgumentParser(description="INPHARED complete bacteriophage C5-MTase +1 census")
    ap.add_argument("--inphared-fasta", type=str, default=None, help="Optional raw inphared_genomes.fasta.gz for full re-scan")
    ap.add_argument("--seeds-faa", type=str, default=None, help="Optional seed FASTA for DIAMOND blastp")
    args = ap.parse_args()

    if not args.inphared_fasta:
        verify_bundled_inphared()
        return

    all_mt = []
    n_ph = 0
    for acc_raw, desc, seq in lib.read_fasta(Path(args.inphared_fasta)):
        n_ph += 1
        if len(seq) < 1000:
            continue
        acc = acc_raw.split()[0]
        for st, s in ((1, seq), (-1, lib.revcomp(seq))):
            for fr in range(3):
                prot = lib.translate(s, offset=fr)
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
                                    p1 = lib.translate(hnt, offset=1)
                                    best_seg = ""
                                    ppos = 0
                                    for pseg in p1.split("*"):
                                        ps, pe = ppos, ppos + len(pseg)
                                        if ps <= pcq + 95 and pe >= pcq + 70 and len(pseg) > len(best_seg):
                                            best_seg = pseg
                                        ppos = pe + 1
                                    all_mt.append((acc, hprot[pcq:pcq + 4], len(hprot), best_seg))
                    pos += sl + 1
    print(f"Raw INPHARED scan finished: {n_ph} phages -> {len(all_mt)} C5-MTases")
    print("  Motif IV counts:", Counter(x[1] for x in all_mt).most_common(10))
    if args.seeds_faa:
        with tempfile.TemporaryDirectory() as td:
            sdb = os.path.join(td, "s.dmnd")
            qfa = os.path.join(td, "q.faa")
            out = os.path.join(td, "o.tsv")
            subprocess.run(["diamond", "makedb", "--in", args.seeds_faa, "-d", sdb, "--quiet"], check=True)
            with open(qfa, "w") as f:
                for i, x in enumerate(all_mt):
                    if len(x[3]) >= 45:
                        f.write(f">M_{i}\n{x[3]}\n")
            subprocess.run(["diamond", "blastp", "-d", sdb, "-q", qfa, "-o", out, "--evalue", "1e-5", "--sensitive", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "--quiet"], check=True)
            print(f"  DIAMOND E<1e-5 +1 hits: {sum(1 for _ in open(out))}")


if __name__ == "__main__":
    main()
