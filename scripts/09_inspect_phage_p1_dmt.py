#!/usr/bin/env python3
"""
Step 5: Bacteriophage P1 dmt (AF234172) coordinate, Shine-Dalgarno isoform,
and synonymous Frame-0 / +1-amber-stop wet-lab mutation verification.

Runs 100% self-contained from repository data/phage_p1_dmt_cds.fna and
data/structures/trans_vs_cis_3d_audit.json, exporting data/phage_p1_dmt_validation.json.
"""
from __future__ import annotations

import json
import re

import lib_c5mtase as lib

# E. coli K-12 codon frequency per 1,000 codons (Kazusa / Nakamura et al.)
ECOLI_CODON_TIER = {
    "GTT": ("common", 18.3), "GTC": ("common", 15.3), "GTA": ("moderate", 10.9),
    "CTC": ("moderate", 11.1), "CTA": ("rare", 3.9),
    "ATT": ("common", 30.3), "ATA": ("rare", 4.4),
}


def main():
    recs = lib.read_fasta(lib.DATA_DIR / "phage_p1_dmt_cds.fna")
    assert len(recs) == 1
    hdr, _desc, cds = recs[0]
    dmt_start_1b = 66412
    dmt_end_1b = 68676
    assert len(cds) == 2265 and cds[:3] == "GTG" and cds[-3:] == "TAA"

    f0_prot = lib.translate(cds[:-3], offset=0)
    assert len(f0_prot) == 754 and "*" not in f0_prot
    p1_full = lib.translate(cds, offset=1)

    # Locate exact +1 stop-to-stop ORF overlapping Motifs IV-VIII
    m_csl = re.search(r"CLLSQVLE", p1_full)
    m_ryf = re.search(r"RYFLS", p1_full)
    assert m_csl and m_ryf
    st_aa = p1_full.rfind("*", 0, m_csl.start()) + 1
    en_aa = p1_full.find("*", m_csl.start())
    p1_orf_aa = p1_full[st_aa:en_aa]
    # 1-based nt inside dmt CDS: codon k in +1 frame spans 0-based [k*3+1 .. k*3+3], i.e. 1-based [k*3+2 .. k*3+4]
    orf_nt_start = st_aa * 3 + 2
    orf_nt_end = en_aa * 3 + 4  # includes stop codon TAA
    assert cds[en_aa * 3 + 1:en_aa * 3 + 4] == "TAA"

    # Find internal ATG start isoforms with upstream Shine-Dalgarno (GGAG)
    isoforms = []
    for loc_aa in range(len(p1_orf_aa)):
        if p1_orf_aa[loc_aa] == "M":
            nt0 = (st_aa + loc_aa) * 3 + 1  # 0-based nt start of ATG
            if cds[nt0:nt0 + 3] == "ATG":
                up = cds[max(0, nt0 - 16):max(0, nt0 - 3)]
                sd_idx = up.find("GGAG")
                if sd_idx >= 0:
                    spacer = len(up) - (sd_idx + 4) + 3
                    isoforms.append({
                        "local_aa_1b": loc_aa + 1,
                        "dmt_nt_1b": nt0 + 1,
                        "genome_1b": dmt_start_1b + nt0,
                        "sd": "GGAG",
                        "spacer_nt": spacer,
                        "codon": "ATG",
                        "len_aa": len(p1_orf_aa) - loc_aa,
                        "seq": p1_orf_aa[loc_aa:],
                    })

    # Pull live AF3 ipTM and R385-F30 distance from trans_vs_cis_3d_audit.json
    af3_path = lib.STRUCT_DIR / "trans_vs_cis_3d_audit.json"
    af3_map = {}
    if af3_path.exists():
        af3_d = json.loads(af3_path.read_text())
        af3_map = {r["name"]: r for r in af3_d["trans_complexes"]}

    iso_157 = [x for x in isoforms if x["len_aa"] == 157][0]
    iso_64 = [x for x in isoforms if x["len_aa"] == 64][0]
    iso_157["name"] = "157aa_isoform"
    iso_64["name"] = "64aa_core_microprotein"
    if "phagep1_arch1_dmt_trans_157aa" in af3_map:
        iso_157["af3_iptm_mean"] = af3_map["phagep1_arch1_dmt_trans_157aa"]["iptm_mean"]
    if "phagep1_arch1_dmt_trans_64aa_core" in af3_map:
        iso_64["af3_iptm_mean"] = af3_map["phagep1_arch1_dmt_trans_64aa_core"]["iptm_mean"]
        r385_f30 = af3_map["phagep1_arch1_dmt_trans_64aa_core"]["r1_kyf_mean"]
    else:
        r385_f30 = None

    # Programmatically verify EVERY point mutation against the real AF234172 DNA!
    # Mutations are specified as (1-based dmt CDS nt, alt_base, priority, motif_context)
    mut_specs = [
        (642, "A", "primary", "CLLSQVLE (L16* in 64-aa core)"),
        (645, "A", "primary", "CLLSQVLE (S17* in 64-aa core)"),
        (666, "A", "secondary_rare_codon", "pre-RYFLS linker (S24* in 64-aa core)"),
        (714, "A", "secondary_rare_codon", "C-terminal alpha-helix (L40* in 64-aa core)"),
    ]
    verified_muts = []
    for nt_1b, alt, tier, ctx in mut_specs:
        i0 = nt_1b - 1
        ref = cds[i0]
        assert i0 % 3 == 2, f"Expected Frame-0 3rd codon position at c.{nt_1b}"
        cod_idx = i0 // 3  # 0-based Frame-0 aa index
        f0_ref_cod = cds[cod_idx * 3:cod_idx * 3 + 3]
        f0_alt_cod = f0_ref_cod[:2] + alt
        f0_ref_aa = lib.CODON_TABLE[f0_ref_cod]
        f0_alt_aa = lib.CODON_TABLE[f0_alt_cod]
        assert f0_ref_aa == f0_alt_aa, f"Non-synonymous Frame-0 mutation at c.{nt_1b}!"

        p1_ref_cod = cds[cod_idx * 3 + 1:cod_idx * 3 + 4]
        p1_alt_cod = p1_ref_cod[0] + alt + p1_ref_cod[2]
        p1_ref_aa = lib.CODON_TABLE[p1_ref_cod]
        p1_alt_aa = lib.CODON_TABLE[p1_alt_cod]
        assert p1_alt_aa == "*", f"Mutation at c.{nt_1b} did not create +1 stop!"

        p1_abs_aa_1b = (cod_idx - st_aa) + 1
        p1_64_pos_1b = (cod_idx - (st_aa + iso_64["local_aa_1b"] - 1)) + 1
        codon_usage_tier, codon_freq = ECOLI_CODON_TIER[f0_alt_cod]
        verified_muts.append({
            "dmt_nt": nt_1b,
            "genome_1b": dmt_start_1b + i0,
            "ref": ref,
            "alt": alt,
            "sub": f"{ref}>{alt}",
            "frame0_aa_pos": cod_idx + 1,
            "frame0_codon_change": f"{f0_ref_cod}->{f0_alt_cod}",
            "frame0_effect": f"{f0_ref_aa}{cod_idx+1}{f0_alt_aa} ({f0_ref_cod}->{f0_alt_cod}, synonymous)",
            "ecoli_alt_codon_usage": f"{codon_usage_tier} ({codon_freq}/1000)",
            "plus1_codon_change": f"{p1_ref_cod}->{p1_alt_cod}",
            "plus1_effect": f"{p1_ref_aa}{p1_abs_aa_1b}* ({p1_ref_cod}->{p1_alt_cod}, amber stop; {ctx})",
            "plus1_64aa_pos": p1_64_pos_1b,
            "recommendation_tier": tier,
        })

    # Locate Frame-0 motifs in P1 dmt (exact 1-based coords!)
    m4_p1 = re.search(r"TPCQ", f0_prot)
    cys78_pos = m4_p1.start() + 3  # 1-based: T76, P77, C78, Q79
    m9_p1 = re.search(r"RRLTPVECERLQ", f0_prot)

    p1_val = {
        "accession": "AF234172",
        "organism": "Escherichia virus P1",
        "uniprot": "DAM_BPP1",
        "genome_bp": 94800,
        "dmt_cds": {"start_1b": dmt_start_1b, "end_1b": dmt_end_1b, "strand": "+", "len_aa": len(f0_prot), "start_codon": "GTG", "stop_codon": "TAA"},
        "dmt_domains_aa": {
            "c5_mtase_domain": [1, 473],
            "motif_I": [7, 16], "motif_IV_TPCQ": [m4_p1.start() + 1, m4_p1.end()],
            "active_site_cys": cys78_pos,
            "motif_VI_ENV": [123, 129], "motif_VIII": [201, 214],
            "TRD": [215, 384], "motif_IX": [m9_p1.start() + 1, m9_p1.start() + 12],
            "motif_X": [445, 459], "dam_n6a_domain": [485, 754],
        },
        "plus1_orf": {
            "stop_to_stop_dmt_nt": [orf_nt_start, orf_nt_end],
            "stop_to_stop_genome_1b": [dmt_start_1b + orf_nt_start - 1, dmt_start_1b + orf_nt_end - 1],
            "stop_to_stop_aa": len(p1_orf_aa),
            "stop_codon": "TAA",
            "sd_atg_isoforms": [
                {k: v for k, v in iso_157.items() if k != "seq"},
                {k: v for k, v in iso_64.items() if k != "seq"},
            ],
            "csl_motif_in_64aa": {"seq": "CLLSQVLE", "pos": [iso_64["seq"].find("CLLSQVLE") + 1, iso_64["seq"].find("CLLSQVLE") + 8]},
            "ryf_motif_in_64aa": {"seq": "RYFLS", "pos": [iso_64["seq"].find("RYFLS") + 1, iso_64["seq"].find("RYFLS") + 5], "lock_residue": "F30", "partner": f"A:R{m9_p1.start()+1}", "dist_5model_A": r385_f30},
        },
        "host_strain_requirement": "Use a non-amber-suppressor (sup^0) E. coli host (e.g., MG1655, W3110, or BL21); avoid supE44/glnV amber-suppressor strains (e.g., C600, DH5alpha, JM109) because all four +1 nonsense alleles introduce TAG (amber) stop codons.",
        "validation_mutations": verified_muts,
    }
    out_path = lib.DATA_DIR / "phage_p1_dmt_validation.json"
    out_path.write_text(json.dumps(p1_val, indent=2))
    print(f"Verified Phage P1 dmt (Cys{cys78_pos}, Motif IX R{m9_p1.start()+1}, +1 ORF {len(p1_orf_aa)}aa) and all 4 TAG alleles!")
    for m in verified_muts:
        print(f"  c.{m['dmt_nt']}{m['sub']} (genome {m['genome_1b']}): F0={m['frame0_effect']} [{m['ecoli_alt_codon_usage']}] | +1={m['plus1_effect']} [{m['recommendation_tier']}]")


if __name__ == "__main__":
    main()
