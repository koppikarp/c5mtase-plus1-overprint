# Pan-Prokaryotic and Viral +1 Overprinted Gene Family in C5-Cytosine DNA Methyltransferases

Code, verified datasets, maximum-likelihood phylogenies, all 30 AlphaFold 3 structural models, editable Fignimble figures, and compiled bioRxiv preprint for:

> **A Pan-Prokaryotic and Viral +1 Overprinted Gene Family Inside Full-Length C5-Cytosine DNA Methyltransferases Explains the Origin of Split-Subunit Methyltransferases**  
> Pratik Koppikar (University of Cambridge)

---

## Summary of Key Findings

1. **Widespread +1 dual-coding gene family inside full-length C5-MTases (Arch 1):** Across 2,339 representative GTDB R232 genomes (45 loci across 41 unique genomes, 41 orders/families, and 13 prokaryotic phyla spanning 11 bacterial and 2 DPANN archaeal phyla) and 5,063 complete INPHARED bacteriophages (39/230 = 17.0% raw; 15/109 = 13.8% dereplicated at <90% catalytic-core identity; 87.5% raw / 83.3% dereplicated of phage TPCQ enzymes; 0/16 SPCQ), an unannotated +1 gene carrying the conserved `CSLSDILE/KYFLS` domain sits inside full-length, Motif-IX-intact C5-MTases (35/39 continuous Motif-IX/X-intact in phages, 4/39 interrupted Frame 0; 0 Arch 2, 0 Arch 3), terminating before the target recognition domain (TRD) and Motif IX.
2. **Hidden inside model bacteriophage P1 `dmt`:** The well-studied 754-aa `dmt` methyltransferase of bacteriophage P1 (GenBank `AF234172:66412-68676`; Swiss-Prot `DAM_BPP1`) carries an unannotated 163-codon stop-free +1 ORF with two internal `GGAG` Shine-Dalgarno + `ATG` starts (157 aa and a compact 64-aa core). Prior BLAST-annotated phage MM1 `orf13`/`orf14` (`AJ302074`) is also shown to be a continuous 385-aa full-length Motif-IX-intact Arch 1 enzyme.
3. **Systematic database erasure:** 29/33 (87.9%) of audited UniProtKB matches among top MMseqs2 UniRef100 hits have been marked deleted in UniProtKB during overlapping-ORF cleanup (`data/uniref100_deletion_audit.json`), with zero Pfam/InterPro families and zero experimental PDB/CATH structural homologs.
4. **Zero-confounder statistical proof of frame-+1 purifying selection:** Whole-gene composition-corrected conservation at strictly non-stop-risk fourfold degenerate (`A/G/P/T`) sites peaks directly over `CSLSDILE/KYFLS` and collapses at the +1 stop codon (gene-level one-sided paired Wilcoxon `P = 1.22e-4` across 13 FPCQ loci [13/13 inside > outside]; strictly disjoint pairs median `P = 0.0312`; independent replication `P = 1.56e-2` across 7 TPCQ/SPCQ loci), paired with zero-filter 5-frame stop-codon purging across all 30 unique intact Arch 1 genomes (within-gene one-sided exact paired Wilcoxon max `P = 1.03e-4`; versus each gene's own 1,000-replicate flanking-codon synonymous null preserving Frame-0 sequence and codon usage: 0.0 vs. 4.7 stops, 30/30 depleted, one-sided `P = 9.31e-10`, two-sided `P = 1.86e-9`; Frame `-0` reduction confined to `+1`-shadowing `TTA`/`CTA` codons `P = 9.12e-6` while non-shadowing `TCA` is untouched `P = 0.973`) and gene-level +1 amino-acid preservation at host-mutated columns (`18.4%` vs. `10.7%`, one-sided Mann-Whitney U `P = 2.40e-6`; disjoint pairs `P = 3.75e-3`).
5. **Evolutionary direction (`Arch 1 -> Arch 2` and `Arch 1 -> Arch 3`):** DNMT-rooted LG+Gamma ML phylogenetics with uniform thorough search settings across main trees and 100 bootstrap replicates (97 unique taxa Tree 1; 67 harmonized taxa Tree 2) shows that single-chain C-terminal fusions (`A0A6N9NLZ4`, Arch 3, ML 97%/SH 99.5% in Tree 1 and ML 68%/SH 94.6% in Tree 2, nested inside a 29-taxon TPCQ superclade at ML 95%/SH 89.8%) evolved once from an Arch 1 TPCQ ancestor, while split-Motif-IX enzymes (`M.EcoHK31I`/`M.EaeI`, canonical Arch 2, ML 62%/SH 81.5% nested in Arch 1 FPCQ at inner ML 77%/SH 91.9%) and the independent *Bacillota* split-TRD/IX/X lineage (`GCA_035934735.1`, nested inside the *Bacillota* Arch 1 SPCQ clade at ML 71%/SH 89.4%) represent derived split-frame transitions. Patristic congruence holds across all 66 shared taxa (`r = 0.56` [`0.67` excl. one long-branch outlier], Mantel `P = 0.0002`) and strictly within the zero-DNA-overlap Arch 3 subset (`n = 7`, `r = 0.75`, `rho = 0.72`, Mantel `P = 0.0004`).
6. **AlphaFold 3 *trans* structural convergence:** Across all 30 diffusion models (6 complexes × 5 models, all bundled in `data/structures/af3_models/`) analyzed with strictly rigid Rossmann-core superposition (`1MHT` to `AF-A0A6N9NLZ4-F1` core RMSD `0.84 A`, excluding the mobile Motif IV base-flipping loop), the unlinked +1 polypeptide docks onto the exact Motif VIII/TRD/Motif IX collar occupied in *cis* by Arch 3 fusions (`ipTM 0.63-0.83`; Phage P1 64-aa core `0.830 ± 0.006`, COM `1.20 ± 0.28 A`, Dock RMSD `2.06 ± 0.26 A`, matching positive control `M.EaeI` `0.840`; negative control `+1`-less `M.HhaI` `0.160`), preserving a conserved `3.69-4.04 A` Motif IX `R1` to `KYFLS F/Y` contact across all three architectures without occluding the catalytic Cys (`~22-28 A`) or SAH (`~19 A`).

---

## Repository Layout

```text
c5mtase-plus1-overprint/
├── preprint/
│   ├── c5mtase_plus1_preprint.tex         # Preprint LaTeX manuscript
│   ├── c5mtase_plus1_preprint.pdf         # Compiled vector PDF
│   └── figures/                           # Native Fignimble (pachterlab/Fignimble) project
│       ├── style.py                       # Shared figure styling and global matte color palette
│       ├── graphical_abstract.py          # Graphical Abstract
│       ├── fig1.py ... fig5.py            # Main Figures 1-5 (183 mm double-column layout)
│       ├── panels/                        # Modular per-panel Matplotlib scripts (data-driven)
│       ├── rendered/                      # Exported vector PDFs
│       └── build_all.sh                   # One-command preflight + export script
├── data/
│   ├── gtdb2339_c5mtase_plus1_census.json # 45 verified GTDB R232 loci across 13 phyla
│   ├── gtdb2339_ctrl_c5mtases.json        # 26 unbiased homology-negative low-GC control C5-MTases
│   ├── inphared_phage_c5mtase_census.json # 39/230 INPHARED phage hits, Motif IV denominators & derep stats
│   ├── uniref100_deletion_audit.json      # 33 audited UniProtKB accessions (29/33 deleted, 0 Pfam/InterPro)
│   ├── uniprot_tpcq_dna_audit.json        # Nucleotide verification of all 5/5 UniProt TPCQ entries
│   ├── selection_4d_and_stop_stats.json   # Unified 4D, genome-bootstrap, 5-frame stops & Test 4 stats
│   ├── phage_p1_dmt_cds.fna               # Exact AF234172 dmt coding DNA (2,265 bp)
│   ├── phage_p1_dmt_validation.json       # Verified P1 coordinates, SD+ATG isoforms & TAG switches
│   ├── phylogeny/                         # Harmonized FASTAs, trimmed alignments, ML & 100-boot Newick trees
│   └── structures/                        # AF3 batch JSON, 30-model audit, all 30 AF3 models, 1MHT, AF-A0A6N9NLZ4
└── scripts/                               # 100% self-contained, repo-relative analysis scripts (00..09 + lib)
```

---

## Reproducing All Analyses, Figures, and the Preprint

### 1. Re-run the full self-contained analysis pipeline (zero external files required)

```bash
python scripts/00_export_all_figure_data.py
```

*(Requires Python 3.10+ with `numpy` and `scipy`; optionally run `python scripts/05_build_ml_phylogeny.py` with `mafft` and `VeryFastTree` on `PATH` to regenerate all ML trees and 100-replicate bootstraps from `data/phylogeny/*_input.faa`.)*

### 2. Re-render or visually edit figures with Fignimble

```bash
uv tool install "git+https://github.com/pachterlab/Fignimble.git"
bash preprint/figures/build_all.sh
```

### 3. Compile the preprint PDF

```bash
cd preprint
tectonic c5mtase_plus1_preprint.tex
```

---

## License and Data Provenance

- Code, figures, and manuscript text: **MIT License** (`LICENSE`).
- 3D *trans* complexes were generated with **AlphaFold Server (AlphaFold 3)** under the Google DeepMind AlphaFold 3 Non-Commercial Terms of Use; *cis* monomer models were retrieved from AlphaFold DB (`CC-BY-4.0`); experimental crystal structure `1MHT` was retrieved from RCSB PDB (`CC0`).
