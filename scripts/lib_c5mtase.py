#!/usr/bin/env python3
"""
Shared, self-contained utilities for the C5-MTase +1 dual-coding repository.
All paths are resolved strictly relative to the repository root (zero external
machine-specific paths or outside module imports).
"""
from __future__ import annotations

import gzip
import math
import random
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# ------------------------------------------------------------------------------
# Repository paths (portable, relative to this file)
# ------------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
PHYLO_DIR = DATA_DIR / "phylogeny"
STRUCT_DIR = DATA_DIR / "structures"
PREPRINT_DIR = REPO_ROOT / "preprint"

# ------------------------------------------------------------------------------
# Genetic code & sequence helpers
# ------------------------------------------------------------------------------
COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S", "CCT": "P", "CCC": "P",
    "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T",
    "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N",
    "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D",
    "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CANON_STARTS = {"ATG", "GTG", "TTG"}
FOURFOLD_NOSTOP = {"A", "G", "P", "T"}  # Excludes V (GTN can generate +1 TGA/TAG/TAA)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def revcomp(seq: str) -> str:
    return seq.translate(COMP_TABLE)[::-1]


def translate(nt: str, table: Optional[Dict[str, str]] = None, offset: int = 0) -> str:
    tbl = table or CODON_TABLE
    n = (len(nt) - offset) // 3
    return "".join(tbl.get(nt[k:k + 3], "X") for k in range(offset, offset + n * 3, 3))


def read_fasta(path: Path) -> List[Tuple[str, str, str]]:
    """Read plain or .gz FASTA into list of (id, desc, seq)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    out: List[Tuple[str, str, str]] = []
    hdr: Optional[str] = None
    buf: List[str] = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if hdr is not None:
                    parts = hdr.split(None, 1)
                    out.append((parts[0], parts[1] if len(parts) > 1 else "", "".join(buf).upper()))
                hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if hdr is not None:
            parts = hdr.split(None, 1)
            out.append((parts[0], parts[1] if len(parts) > 1 else "", "".join(buf).upper()))
    return out


def read_fna_gz(path: Path) -> Dict[str, str]:
    return {cid: seq for cid, _desc, seq in read_fasta(path)}


def nw_align(s1: str, s2: str, match: int = 2, mis: int = -1, gap: int = -4) -> Tuple[str, str]:
    """Global Needleman-Wunsch alignment for protein sequences."""
    n, m = len(s1), len(s2)
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    tb = np.zeros((n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        dp[i, 0] = i * gap
        tb[i, 0] = 1
    for j in range(1, m + 1):
        dp[0, j] = j * gap
        tb[0, j] = 2
    for i in range(1, n + 1):
        c1 = s1[i - 1]
        dp_im1 = dp[i - 1]
        dp_i = dp[i]
        tb_i = tb[i]
        for j in range(1, m + 1):
            d = dp_im1[j - 1] + (match if c1 == s2[j - 1] else mis)
            u = dp_im1[j] + gap
            l = dp_i[j - 1] + gap
            if d >= u and d >= l:
                dp_i[j] = d
                tb_i[j] = 0
            elif u >= l:
                dp_i[j] = u
                tb_i[j] = 1
            else:
                dp_i[j] = l
                tb_i[j] = 2
    a1: List[str] = []
    a2: List[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        t = tb[i, j]
        if t == 0:
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif t == 1:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        else:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def gc3_4d(hprot: str, hnt: str, lo: int = 0, hi: Optional[int] = None) -> float:
    """GC3 at non-stop-risk 4-fold degenerate sites (A, G, P, T)."""
    hi_i = len(hprot) if hi is None else min(len(hprot), hi)
    lo_i = max(0, lo)
    b3 = Counter(
        hnt[k * 3 + 2]
        for k in range(lo_i, hi_i)
        if hprot[k] in FOURFOLD_NOSTOP and k * 3 + 2 < len(hnt)
    )
    tot = sum(b3[b] for b in "ACGT")
    return (b3["G"] + b3["C"]) / max(1, tot)


def whole_gene_4d_exp(hprot1: str, hnt1: str, hprot2: str, hnt2: str) -> float:
    """Expected 3rd-position identity at A/G/P/T 4D sites from whole-gene base frequencies."""
    b1 = Counter(hnt1[k * 3 + 2] for k in range(len(hprot1)) if hprot1[k] in FOURFOLD_NOSTOP and k * 3 + 2 < len(hnt1))
    b2 = Counter(hnt2[k * 3 + 2] for k in range(len(hprot2)) if hprot2[k] in FOURFOLD_NOSTOP and k * 3 + 2 < len(hnt2))
    t1 = max(1, sum(b1[b] for b in "ACGT"))
    t2 = max(1, sum(b2[b] for b in "ACGT"))
    return sum((b1[b] / t1) * (b2[b] / t2) for b in "ACGT")


def count_5frames_stops(hnt: str, pcq: int, rel_lo: int = -25, rel_hi: int = 110, table: Optional[Dict[str, str]] = None) -> Dict[str, int]:
    """Count stop codons in all 5 alternative reading frames across [pcq+rel_lo .. pcq+rel_hi]."""
    tbl = table or CODON_TABLE
    lo_nt = max(0, (pcq + rel_lo) * 3)
    hi_nt = min(len(hnt), (pcq + rel_hi) * 3)
    sub = hnt[lo_nt:hi_nt]
    rc = revcomp(sub)
    out: Dict[str, int] = {}
    for ft, s, off in (("+1", sub, 1), ("+2", sub, 2), ("-0", rc, 0), ("-1", rc, 1), ("-2", rc, 2)):
        nc = (len(s) - off) // 3
        out[ft] = sum(1 for k in range(off, off + nc * 3, 3) if tbl.get(s[k:k + 3], "X") == "*")
    return out


# ------------------------------------------------------------------------------
# Newick / phylogenetic tree helpers
# ------------------------------------------------------------------------------
class Node:
    def __init__(self, name: str = "", length: float = 0.0):
        self.name = name
        self.length = length
        self.children: List[Node] = []
        self.parent: Optional[Node] = None
        self.support: Optional[float] = None

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def leaves(self) -> List[Node]:
        if self.is_leaf():
            return [self]
        out: List[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out

    def leaf_names(self) -> set:
        return {l.name for l in self.leaves()}


def parse_newick(s: str) -> Node:
    s = s.strip().rstrip(";")
    idx = 0
    n = len(s)

    def parse_node() -> Node:
        nonlocal idx
        node = Node()
        if idx < n and s[idx] == "(":
            idx += 1
            while True:
                ch = parse_node()
                ch.parent = node
                node.children.append(ch)
                if idx < n and s[idx] == ",":
                    idx += 1
                elif idx < n and s[idx] == ")":
                    idx += 1
                    break
                else:
                    break
        start = idx
        while idx < n and s[idx] not in ",):":
            idx += 1
        lbl = s[start:idx].strip()
        if lbl:
            if node.is_leaf():
                node.name = lbl
            else:
                try:
                    node.support = float(lbl)
                except ValueError:
                    node.name = lbl
        if idx < n and s[idx] == ":":
            idx += 1
            bstart = idx
            while idx < n and s[idx] not in ",)":
                idx += 1
            try:
                node.length = float(s[bstart:idx].strip())
            except ValueError:
                node.length = 0.0
        return node

    return parse_node()


def reroot_on_clade(root: Node, outgroup_names: Iterable[str]) -> Node:
    """Reroot unrooted/rooted Newick tree on the edge separating outgroup_names from the rest."""
    og = set(outgroup_names)
    adj: Dict[Node, List[Tuple[Node, float, Optional[float]]]] = defaultdict(list)

    def collect(u: Node):
        for c in u.children:
            adj[u].append((c, c.length, c.support))
            adj[c].append((u, c.length, c.support))
            collect(c)

    collect(root)
    all_leaves = {nd.name for nd in adj if len(adj[nd]) == 1 and nd.name}
    target_edge = None
    for u in list(adj.keys()):
        for v, bl, sup in adj[u]:
            if id(u) >= id(v):
                continue
            seen = {u}
            stack = [v]
            lvs = set()
            while stack:
                cur = stack.pop()
                if cur in seen:
                    continue
                seen.add(cur)
                if len(adj[cur]) == 1 and cur.name:
                    lvs.add(cur.name)
                for nb, _, _ in adj[cur]:
                    if nb not in seen:
                        stack.append(nb)
            if lvs == og or (all_leaves - lvs) == og:
                target_edge = (u, v, bl, sup)
                break
        if target_edge:
            break
    if not target_edge:
        return root
    u, v, bl, _ = target_edge
    adj[u] = [x for x in adj[u] if x[0] is not v]
    adj[v] = [x for x in adj[v] if x[0] is not u]
    new_root = Node("ROOT", 0.0)
    adj[new_root].append((u, bl / 2.0, None))
    adj[u].append((new_root, bl / 2.0, None))
    adj[new_root].append((v, bl / 2.0, None))
    adj[v].append((new_root, bl / 2.0, None))

    def build(cur: Node, par: Optional[Node], bl_in: float, sup_in: Optional[float]) -> Node:
        nd = Node(cur.name if len(adj[cur]) == 1 else "", bl_in)
        nd.support = sup_in
        for nb, w, sup in adj[cur]:
            if nb is par:
                continue
            ch = build(nb, cur, w, sup)
            ch.parent = nd
            nd.children.append(ch)
        return nd

    return build(new_root, None, 0.0, None)


def midpoint_reroot(root: Node) -> Node:
    """Reroot tree at the midpoint of the longest leaf-to-leaf path."""
    adj: Dict[Node, List[Tuple[Node, float, Optional[float]]]] = defaultdict(list)

    def collect(u: Node):
        for c in u.children:
            adj[u].append((c, max(0.0, c.length), c.support))
            adj[c].append((u, max(0.0, c.length), c.support))
            collect(c)

    collect(root)
    leaves = [nd for nd in adj if len(adj[nd]) == 1]

    def farthest(start: Node) -> Tuple[Node, float, Dict[Node, Optional[Node]], Dict[Node, float]]:
        dist = {start: 0.0}
        par: Dict[Node, Optional[Node]] = {start: None}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, w, _ in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + w
                    par[v] = u
                    stack.append(v)
        best = max(leaves, key=lambda x: dist[x])
        return best, dist[best], par, dist

    a, _, _, _ = farthest(leaves[0])
    b, diameter, par, dist = farthest(a)
    half = diameter / 2.0
    cur = b
    while par[cur] is not None and dist[par[cur]] > half:
        cur = par[cur]
    p = par[cur]
    if p is None:
        return root
    edge_len = dist[cur] - dist[p]
    rem_cur = dist[cur] - half
    rem_p = half - dist[p]
    adj[p] = [x for x in adj[p] if x[0] is not cur]
    adj[cur] = [x for x in adj[cur] if x[0] is not p]
    nr = Node("MIDROOT", 0.0)
    adj[nr].append((cur, rem_cur, None))
    adj[cur].append((nr, rem_cur, None))
    adj[nr].append((p, rem_p, None))
    adj[p].append((nr, rem_p, None))

    def build(u: Node, parent: Optional[Node], bl_in: float, sup_in: Optional[float]) -> Node:
        nd = Node(u.name if len(adj[u]) == 1 else "", bl_in)
        nd.support = sup_in
        for nb, w, sup in adj[u]:
            if nb is parent:
                continue
            ch = build(nb, u, w, sup)
            ch.parent = nd
            nd.children.append(ch)
        return nd

    return build(nr, None, 0.0, None)


def patristic_matrix(root: Node) -> Tuple[List[str], np.ndarray]:
    adj: Dict[Node, List[Tuple[Node, float]]] = defaultdict(list)

    def collect(u: Node):
        for c in u.children:
            adj[u].append((c, max(0.0, c.length)))
            adj[c].append((u, max(0.0, c.length)))
            collect(c)

    collect(root)
    leaves = sorted([nd for nd in adj if len(adj[nd]) == 1 and nd.name], key=lambda x: x.name)
    n = len(leaves)
    idx = {nd: i for i, nd in enumerate(leaves)}
    mat = np.zeros((n, n), dtype=float)
    for i, start in enumerate(leaves):
        dist = {start: 0.0}
        stack = [start]
        while stack:
            u = stack.pop()
            du = dist[u]
            for v, w in adj[u]:
                if v not in dist:
                    dist[v] = du + w
                    stack.append(v)
        for nd, j in idx.items():
            mat[i, j] = dist[nd]
    return [l.name for l in leaves], mat


def pearson_r(x: Sequence[float], y: Sequence[float]) -> float:
    ax = np.asarray(x, dtype=float)
    ay = np.asarray(y, dtype=float)
    dx = ax - ax.mean()
    dy = ay - ay.mean()
    den = float(np.sqrt(np.sum(dx * dx) * np.sum(dy * dy)))
    return float(np.sum(dx * dy) / den) if den > 0 else 0.0


def rank_vec(vals: Sequence[float]) -> List[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg_r = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg_r
        i = j + 1
    return ranks


def mantel_test(m1: np.ndarray, m2: np.ndarray, n_perm: int = 5000, seed: int = 42) -> Dict[str, float]:
    """Exact permutation Mantel test (one-sided greater) returning Pearson r, Spearman rho, and P."""
    n = m1.shape[0]
    iu = np.triu_indices(n, k=1)
    v1 = m1[iu]
    v2 = m2[iu]
    r_obs = pearson_r(v1, v2)
    rho_obs = pearson_r(rank_vec(v1), rank_vec(v2))
    rng = np.random.default_rng(seed)
    ge_r = 0
    ge_rho = 0
    rv1 = np.asarray(rank_vec(v1), dtype=float)
    for _ in range(n_perm):
        perm = rng.permutation(n)
        vp = m2[np.ix_(perm, perm)][iu]
        if pearson_r(v1, vp) >= r_obs:
            ge_r += 1
        if pearson_r(rv1, rank_vec(vp)) >= rho_obs:
            ge_rho += 1
    return {
        "n_taxa": int(n),
        "n_pairs": int(len(v1)),
        "pearson_r": round(r_obs, 4),
        "spearman_rho": round(rho_obs, 4),
        "mantel_p_pearson": round((ge_r + 1) / (n_perm + 1), 6),
        "mantel_p_spearman": round((ge_rho + 1) / (n_perm + 1), 6),
        "n_perm": int(n_perm),
    }


def norm_id(name: str) -> str:
    if name in ("A2_UP_M_EaeI", "A2_UP_O52510_M_EaeI_beta"):
        return "A2_EaeI"
    return name


def classify_taxon(name: str) -> str:
    if name.startswith("OUT_Arch"):
        return "NOp1_ARCH_OUT"
    if name.startswith("OUT_"):
        return "OUT_EUK"
    if name.startswith("A3_"):
        return "ARCH3"
    if name.startswith("A2_"):
        return "ARCH2"
    if name.startswith("A1_TPCQ_"):
        return "ARCH1_TPCQ"
    if name.startswith("A1_FPCQ_"):
        return "ARCH1_FPCQ"
    if name.startswith("A1_SPCQ_") or name.startswith("A1_SPCT_"):
        return "ARCH1_SPCQ"
    if name.startswith("NOp1_ARCH_"):
        return "NOp1_ARCH_OUT"
    if name.startswith("NOp1_TPCQ"):
        return "NOp1_TPCQ"
    if name.startswith("NOp1_SPCQ") or name.startswith("NOp1_SPCT"):
        return "NOp1_SPCQ"
    if name.startswith("NOp1_"):
        return "NOp1_FPCQ"
    return "OTHER"


def get_unrooted_splits(root: Node) -> set:
    """Return canonical set of unrooted non-trivial bipartitions (frozenset of smaller/lexicographical side)."""
    all_lvs = frozenset(root.leaf_names())
    splits = set()

    def walk(nd: Node):
        lvs = frozenset(nd.leaf_names())
        comp = all_lvs - lvs
        if 1 < len(lvs) < len(all_lvs) - 1:
            key = lvs if (len(lvs) < len(comp) or (len(lvs) == len(comp) and min(lvs) < min(comp))) else comp
            splits.add(key)
        for c in nd.children:
            walk(c)

    walk(root)
    return splits


def canonical_split(clade: Iterable[str], all_leaves: Iterable[str]) -> frozenset:
    s = frozenset(clade)
    all_s = frozenset(all_leaves)
    comp = all_s - s
    return s if (len(s) < len(comp) or (len(s) == len(comp) and min(s) < min(comp))) else comp


# ------------------------------------------------------------------------------
# 3D structure (mmCIF) & rigid-core C5-MTase superposition helpers
# ------------------------------------------------------------------------------
# Rigid-core anchors EXCLUDE the mobile base-flipping catalytic loop M4[+2..+8]
# (residues 81-88 in M.HhaI 1MHT, which displaces ~27 A upon DNA binding!).
RIGID_MT_ANCHORS = [
    ("M1", -4, 12),   # Motif I SAM-binding Rossmann strand + GxG loop
    ("M4", -4, 2),    # Motif IV rigid pre-catalytic beta-strand [-4..+1] ONLY!
    ("M6", -6, 10),   # Motif VI ENV catalytic core strand
    ("M8", -4, 8),    # Motif VIII QxRxR core helix
    ("M10", -4, 10),  # Motif X GNxxV C-terminal Rossmann helix
]


def parse_cif_atoms(path: Path) -> Dict[str, List[dict]]:
    """Parse mmCIF ATOM/HETATM records grouped by auth_asym_id."""
    chains: Dict[str, List[dict]] = defaultdict(list)
    cols: List[str] = []
    in_atom = False
    with open(path) as f:
        for line in f:
            ls = line.strip()
            if ls.startswith("_atom_site."):
                cols.append(ls.split(".")[1].strip())
                in_atom = True
                continue
            if in_atom and (ls.startswith("ATOM ") or ls.startswith("HETATM ")):
                parts = ls.split()
                if len(parts) < len(cols):
                    continue
                rec = dict(zip(cols, parts))
                ch = rec.get("auth_asym_id") or rec.get("label_asym_id")
                seq_s = rec.get("auth_seq_id") or rec.get("label_seq_id")
                try:
                    resseq = int(seq_s)
                except (TypeError, ValueError):
                    resseq = -1
                chains[ch].append({
                    "group": rec.get("group_PDB", "ATOM"),
                    "atom": rec.get("label_atom_id", "").strip('"'),
                    "comp": rec.get("label_comp_id", ""),
                    "resseq": resseq,
                    "xyz": np.array([
                        float(rec["Cartn_x"]),
                        float(rec["Cartn_y"]),
                        float(rec["Cartn_z"]),
                    ], dtype=float),
                    "bfac": float(rec.get("B_iso_or_equiv", 0.0)),
                })
            elif in_atom and (ls.startswith("#") or ls.startswith("loop_")):
                if cols:
                    in_atom = False
    return dict(chains)


def chain_ca_sequence(atoms: List[dict]) -> Tuple[List[int], str, np.ndarray, np.ndarray]:
    """Extract sorted unique CA residues: (resseq_list, seq_1letter, ca_coords[N,3], bfac[N])."""
    seen = set()
    ids, seq, coords, bfs = [], [], [], []
    for a in atoms:
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["comp"] in AA3TO1:
            if a["resseq"] in seen:
                continue
            seen.add(a["resseq"])
            ids.append(a["resseq"])
            seq.append(AA3TO1[a["comp"]])
            coords.append(a["xyz"])
            bfs.append(a["bfac"])
    return ids, "".join(seq), np.asarray(coords, dtype=float), np.asarray(bfs, dtype=float)


def kabsch_align(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Find R, t minimizing || (P @ R.T + t) - Q ||_F, returning (R, t, rmsd)."""
    p_c = P.mean(axis=0)
    q_c = Q.mean(axis=0)
    P0 = P - p_c
    Q0 = Q - q_c
    H = P0.T @ Q0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, -1.0 if d < 0 else 1.0])
    R = Vt.T @ D @ U.T
    t = q_c - p_c @ R.T
    P_al = P @ R.T + t
    rmsd = float(np.sqrt(np.mean(np.sum((P_al - Q) ** 2, axis=1))))
    return R, t, rmsd


def find_mtase_landmarks(seq: str) -> Dict[str, int]:
    """Locate canonical C5-MTase core Rossmann motifs I, IV (P of PCQ), VI, VIII, IX, X."""
    lm: Dict[str, int] = {}
    m4 = re.search(r"[FTSAVILYCMW]PC[QSNT]", seq)
    if m4:
        lm["M4"] = m4.start() + 1  # index of P
    m1 = re.search(r"[LIVMFY][DF][LIVMFY][FWY][ASCG]G[CG][GSA]", seq)
    if not m1:
        m1 = re.search(r"[LIVMFY].[LIVMFY][FWY][ASCG]G[CG]", seq[:120])
    if m1:
        lm["M1"] = m1.start()
    m6 = re.search(r"EN[VLIM]", seq[lm.get("M4", 0):]) if "M4" in lm else re.search(r"EN[VLIM]", seq)
    if m6:
        lm["M6"] = m6.start() + lm.get("M4", 0)
    m8 = re.search(r"Q.[RK].R[LIVMFY]", seq[lm.get("M6", 0):]) if "M6" in lm else None
    if m8:
        lm["M8"] = m8.start() + lm.get("M6", 0)
    m9 = re.search(r"R[RK][LMIFV][TIVSLH]P.{2,4}R[LVIM][QMV]", seq) or re.search(r"RKLHPRECAR", seq)
    if m9:
        lm["M9"] = m9.start()
    m10 = re.search(r"GN[ASGT][VMI][TAVSIC]", seq[lm.get("M8", 0):]) if "M8" in lm else re.search(r"GN[ASGT][VMI]", seq[-150:])
    if m10:
        lm["M10"] = m10.start() + (lm.get("M8", 0) if "M8" in lm else max(0, len(seq) - 150))
    return lm


def superpose_mtase_rigid_core(
    seq_mob: str, ca_mob: np.ndarray, seq_ref: str, ca_ref: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """Superpose mobile C5-MTase onto reference using strictly rigid Rossmann core motifs."""
    lm_m = find_mtase_landmarks(seq_mob)
    lm_r = find_mtase_landmarks(seq_ref)
    pts_m, pts_r = [], []
    for key, lo, hi in RIGID_MT_ANCHORS:
        if key in lm_m and key in lm_r:
            im, ir = lm_m[key], lm_r[key]
            for off in range(lo, hi):
                pm, pr = im + off, ir + off
                if 0 <= pm < len(ca_mob) and 0 <= pr < len(ca_ref):
                    pts_m.append(ca_mob[pm])
                    pts_r.append(ca_ref[pr])
    R, t, rmsd = kabsch_align(np.asarray(pts_m, dtype=float), np.asarray(pts_r, dtype=float))
    return R, t, rmsd, len(pts_m)
