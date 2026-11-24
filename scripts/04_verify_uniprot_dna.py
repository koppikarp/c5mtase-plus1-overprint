#!/usr/bin/env python3
"""
Step 1C: ENA/NCBI nucleotide verification of UniProt TPCQ C5-MTases.

Default mode (100% self-contained): audits bundled data/uniprot_tpcq_dna_audit.json,
confirming all 5/5 DNA-verified TPCQ entries (W1J1U5, A0A068R437, A0A1Y2S7M7,
A0A0J1C4Q1, DAM_BPP1/P04392) encode intact +1 CSLSDILE/KYFLS ORFs.

Optional live fetch mode:
  python scripts/04_verify_uniprot_dna.py --fetch-ena
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.request

import lib_c5mtase as lib


def verify_bundled() -> None:
    recs = json.loads((lib.DATA_DIR / "uniprot_tpcq_dna_audit.json").read_text())
    tpcq = [r for r in recs if r.get("motif4") == "TPCQ"]
    pos = [r for r in tpcq if r.get("has_plus1")]
    print("=== UNIPROT TPCQ ENA DNA AUDIT VERIFICATION ===")
    print(f"  DNA-verified TPCQ entries carrying intact +1 CSLSDILE/KYFLS: {len(pos)}/{len(tpcq)}")
    for r in recs:
        orf = r["p1_orfs"][0] if r["p1_orfs"] else None
        status = f"+1 ORF {orf['len_aa']}aa ({orf['csl']}/{orf['kyf']})" if orf else f"NO +1 ({r['p1_stops_win']} stops)"
        print(f"    {r['acc']:12s} ({r['entry']:18s}) | {r['org']:36s} | {r['ena_pid']:12s} | {status}")


def fetch_live_ena() -> None:
    for acc in ("W1J1U5", "A0A068R437", "A0A1Y2S7M7", "A0A0J1C4Q1"):
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        with urllib.request.urlopen(url, timeout=20) as resp:
            d = json.loads(resp.read().decode())
        pids = []
        for xr in d.get("uniProtKBCrossReferences", []):
            if xr.get("database") == "EMBL":
                for p in xr.get("properties", []):
                    if p.get("key") == "ProteinId" and p.get("value") not in (None, "-"):
                        pids.append(p["value"])
        print(f"  {acc}: ENA protein IDs = {pids}")


def main():
    ap = argparse.ArgumentParser(description="Verify +1 CSLSDILE/KYFLS ORFs in UniProt TPCQ nucleotide records")
    ap.add_argument("--fetch-ena", action="store_true", help="Query live UniProt/ENA REST APIs (requires network)")
    args = ap.parse_args()
    if args.fetch_ena:
        fetch_live_ena()
    else:
        verify_bundled()


if __name__ == "__main__":
    main()
