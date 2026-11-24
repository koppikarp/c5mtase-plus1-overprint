"""Figure 4, Panel c: Real atomic-coordinate 3D superposition of Phage P1 trans AF3 complex onto cis A0A6N9NLZ4 & 1MHT DNA."""
import json
import re
import numpy as np
from matplotlib.collections import LineCollection
from style import DATA_DIR, PAL

AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}


def _parse_cif(path):
    lines = path.read_text().splitlines()
    atoms = []
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            i += 1
            cols = []
            while i < len(lines) and lines[i].strip().startswith("_"):
                cols.append(lines[i].strip()); i += 1
            if cols and cols[0].startswith("_atom_site."):
                ci = {c: j for j, c in enumerate(cols)}
                while i < len(lines):
                    s = lines[i].strip()
                    if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                        break
                    p = s.split()
                    if len(p) >= len(cols):
                        rec = p[ci["_atom_site.group_PDB"]]
                        if rec in ("ATOM", "HETATM"):
                            an = p[ci["_atom_site.label_atom_id"]].strip('"')
                            ch = p[ci.get("_atom_site.auth_asym_id", ci["_atom_site.label_asym_id"])]
                            comp = p[ci["_atom_site.label_comp_id"]]
                            ss = p[ci.get("_atom_site.auth_seq_id", ci["_atom_site.label_seq_id"])]
                            if ss in (".", "?"): ss = p[ci["_atom_site.label_seq_id"]]
                            try:
                                atoms.append((rec, ch, int(ss), comp, an,
                                              np.array([float(p[ci["_atom_site.Cartn_x"]]),
                                                        float(p[ci["_atom_site.Cartn_y"]]),
                                                        float(p[ci["_atom_site.Cartn_z"]])])))
                            except ValueError:
                                pass
                    i += 1
                break
        else:
            i += 1
    return atoms


def _chain_ca(atoms, chain):
    d = {}
    for rec, ch, resi, comp, an, xyz in atoms:
        if rec == "ATOM" and ch == chain and an == "CA" and comp in AA3 and resi not in d:
            d[resi] = (AA3[comp], xyz)
    ids = sorted(d.keys())
    return ids, "".join(d[r][0] for r in ids), np.array([d[r][1] for r in ids])


def _kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cQ - cP @ R.T


def _landmarks(seq):
    lm = {}
    m4 = re.search(r"[GASCVITLFY][FTSAVILYCMWP]PC[QSNTG]", seq)
    if m4: lm["M4"] = m4.start() + 2
    m6 = re.search(r"[IVLMFWYC]ENV", seq)
    if m6: lm["M6"] = m6.start() + 1
    m10 = re.search(r"GN[SA][VWILMT][ATVIL]", seq)
    if m10: lm["M10"] = m10.start()
    return lm


def _smooth_catmull(pts, subdiv=3):
    if len(pts) < 4:
        return pts
    out = []
    p = np.vstack([pts[0], pts, pts[-1]])
    for i in range(1, len(p) - 2):
        p0, p1, p2, p3 = p[i - 1], p[i], p[i + 1], p[i + 2]
        for s in range(subdiv):
            t = s / subdiv
            t2, t3 = t * t, t * t * t
            pt = 0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t3)
            out.append(pt)
    out.append(pts[-1])
    return np.array(out)


def draw(fig):
    ax = fig.add_subplot(111)

    # Load real coordinates: cis A0A6N9NLZ4, trans Phage P1 64aa AF3, and 1MHT DNA
    a0_at = _parse_cif(DATA_DIR / "structures/AF-A0A6N9NLZ4-F1-model_v6.cif")
    p1_at = _parse_cif(DATA_DIR / "structures/af3_models/fold_phagep1_arch1_dmt_trans_64aa_core_model_0.cif")
    mht_at = _parse_cif(DATA_DIR / "structures/1MHT.cif")
    audit = json.loads((DATA_DIR / "structures/trans_vs_cis_3d_audit.json").read_text())
    p1_rec = [r for r in audit["trans_complexes"] if r["name"] == "phagep1_arch1_dmt_trans_64aa_core"][0]

    _, a0_seq, a0_ca = _chain_ca(a0_at, "A")
    _, p1A_seq, p1A_ca = _chain_ca(p1_at, "A")
    _, p1B_seq, p1B_ca = _chain_ca(p1_at, "B")
    _, mht_seq, mht_ca = _chain_ca(mht_at, "A")

    lm0, lm1, lmH = _landmarks(a0_seq), _landmarks(p1A_seq), _landmarks(mht_seq)

    def align_to_a0(seq, ca, lm):
        # Strictly rigid Rossmann core anchors: M4[-4..+1] EXCLUDES the mobile catalytic loop M4[+2..+8]!
        P, Q = [], []
        for k, a, b in [("M4", -4, 2), ("M6", -3, 7), ("M10", -1, 7)]:
            i1, i2 = lm[k], lm0[k]
            for off in range(a, b):
                P.append(ca[i1 + off]); Q.append(a0_ca[i2 + off])
        return _kabsch(np.array(P), np.array(Q))

    R1, t1 = align_to_a0(p1A_seq, p1A_ca, lm1)
    p1A_aln = p1A_ca @ R1.T + t1
    p1B_aln = p1B_ca @ R1.T + t1

    RH, tH = align_to_a0(mht_seq, mht_ca, lmH)
    dna_by_ch = {}
    for rec, ch, resi, comp, an, xyz in mht_at:
        if comp in ("DC", "DG", "DA", "DT", "5FC") and an == "P":
            dna_by_ch.setdefault(ch, []).append(xyz)
    dna_chains = [np.array(v) @ RH.T + tH for _, v in sorted(dna_by_ch.items()) if len(v) >= 4]
    dna_C, dna_D = dna_chains[0], dna_chains[1]

    # Key landmark coordinates in a0 frame
    cys_pt = p1A_aln[77]     # Cys78 (0-based 77)
    r385_pt = p1A_aln[384]   # Motif IX R385 (0-based 384)
    f30_pt = p1B_aln[29]     # RYFLS F30 (0-based 29)
    dna_mid = 0.5 * (dna_C.mean(0) + dna_D.mean(0))

    # Build a clean camera basis: x-axis from Cys78 -> F30, y-axis toward DNA duplex!
    vx = f30_pt - cys_pt
    vx /= np.linalg.norm(vx)
    vy = (dna_mid - 0.5 * (cys_pt + f30_pt))
    vy = vy - np.dot(vy, vx) * vx
    vy /= np.linalg.norm(vy)
    vz = np.cross(vx, vy)
    # Slight tilt for clean 3D depth separation
    ang = np.deg2rad(16)
    vy2 = np.cos(ang) * vy + np.sin(ang) * vz
    vz2 = -np.sin(ang) * vy + np.cos(ang) * vz
    basis = np.vstack([vx, vy2, vz2])  # 3x3

    origin = p1A_aln.mean(0)
    def proj(pts):
        return (pts - origin) @ basis.T

    P_A = proj(p1A_aln)
    P_B = proj(p1B_aln)
    P_cis = proj(a0_ca[387:446])  # cis CSLSDILE/KYFLS domain
    P_dC = proj(dna_C)
    P_dD = proj(dna_D)

    # Collect depth-sorted ribbon segments: (z_mid, x0, y0, x1, y1, color, lw, alpha)
    segs = []
    def add_ribbon(pts3, col, lw=2.0, alpha=0.92, subdiv=2):
        sm = _smooth_catmull(pts3, subdiv=subdiv)
        for k in range(len(sm) - 1):
            zm = 0.5 * (sm[k, 2] + sm[k + 1, 2])
            segs.append((zm, sm[k, 0], sm[k, 1], sm[k + 1, 0], sm[k + 1, 1], col, lw, alpha))

    # Chain A catalytic core (subtle matte slate) + Motif VIII (brighter slate) + Motif IX (sage green!)
    add_ribbon(P_A[:200], "#9BB1CE", lw=1.55, alpha=0.68, subdiv=2)
    add_ribbon(P_A[200:215], PAL["frame0"], lw=2.4, alpha=0.95, subdiv=2)
    add_ribbon(P_A[215:383], "#B5C5DA", lw=1.45, alpha=0.58, subdiv=2)
    add_ribbon(P_A[383:398], PAL["motif9"], lw=2.8, alpha=0.98, subdiv=3)
    add_ribbon(P_A[398:], "#9BB1CE", lw=1.55, alpha=0.68, subdiv=2)

    # DNA phosphate strands + base-pair rungs (matte charcoal-sand)
    add_ribbon(P_dC, PAL["dna"], lw=2.2, alpha=0.72, subdiv=2)
    add_ribbon(P_dD, PAL["dna"], lw=2.2, alpha=0.72, subdiv=2)
    n_rungs = min(len(P_dC), len(P_dD))
    for k in range(n_rungs):
        p1, p2 = P_dC[k], P_dD[n_rungs - 1 - k]
        segs.append((0.5 * (p1[2] + p2[2]), p1[0], p1[1], p2[0], p2[1], "#9E958C", 0.9, 0.55))

    # Cis A0A6N9NLZ4 CSLSDILE/KYFLS overlay (petrol teal outline!)
    add_ribbon(P_cis, PAL["arch3"], lw=1.9, alpha=0.78, subdiv=2)

    # Trans Phage P1 64aa +1 microprotein (bold matte terracotta!)
    add_ribbon(P_B, PAL["plus1"], lw=2.9, alpha=0.98, subdiv=3)

    # Sort back-to-front by z-depth and draw with white PyMOL-style casing
    segs.sort(key=lambda s: s[0])
    for zm, x0, y0, x1, y1, col, lw, alp in segs:
        ax.plot([x0, x1], [y0, y1], color="white", lw=lw + 1.1, alpha=alp * 0.75, solid_capstyle="round", zorder=2)
        ax.plot([x0, x1], [y0, y1], color=col, lw=lw, alpha=alp, solid_capstyle="round", zorder=3)

    # Highlight key atomic lock spheres & distances
    p_cys = proj(cys_pt[None])[0]
    p_r385 = proj(r385_pt[None])[0]
    p_f30 = proj(f30_pt[None])[0]
    p_dna_lbl = P_dC[np.argmax(P_dC[:, 1])]

    ax.scatter([p_cys[0]], [p_cys[1]], s=44, fc=PAL["motif4"], ec=PAL["ink"], lw=0.7, zorder=8)
    ax.scatter([p_r385[0]], [p_r385[1]], s=52, fc=PAL["motif9"], ec=PAL["ink"], lw=0.8, zorder=9)
    ax.scatter([p_f30[0]], [p_f30[1]], s=52, fc=PAL["plus1"], ec=PAL["ink"], lw=0.8, zorder=9)
    ax.plot([p_r385[0], p_f30[0]], [p_r385[1], p_f30[1]], ls="--", lw=1.1, color=PAL["ink"], zorder=10)

    # Clean sparse callouts with leader arrows
    all_xy = np.vstack([P_A[:, :2], P_B[:, :2], P_cis[:, :2], P_dC[:, :2], P_dD[:, :2]])
    xmin, ymin = all_xy.min(0) - 6.0
    xmax, ymax = all_xy.max(0) + 6.5

    ax.annotate(f"DNA duplex\n(1MHT, {p1_rec['dna_min_dist']:.1f} Å Cα)", xy=(p_dna_lbl[0], p_dna_lbl[1]),
                xytext=(xmin + 2.0, ymax - 2.0), fontsize=5.5, color=PAL["dna"], fontweight="bold",
                ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=PAL["dna"], lw=0.6))
    ax.annotate(f"Cys78\n(Motif IV, {p1_rec['cys_min_dist_ca_mean']:.0f} Å)", xy=(p_cys[0], p_cys[1]),
                xytext=(xmin + 1.5, ymin + 3.0), fontsize=5.5, color=PAL["motif4"], fontweight="bold",
                ha="left", va="bottom",
                arrowprops=dict(arrowstyle="-", color=PAL["motif4"], lw=0.6))
    ax.annotate(f"Motif IX R385\n↔ F30 ({p1_rec['r1_kyf_mean']:.2f} Å)", xy=(0.5 * (p_r385[0] + p_f30[0]), 0.5 * (p_r385[1] + p_f30[1])),
                xytext=(xmax - 1.5, ymax - 2.0), fontsize=5.6, color=PAL["ink"], fontweight="bold",
                ha="right", va="top",
                arrowprops=dict(arrowstyle="-|>", color=PAL["ink"], lw=0.7, mutation_scale=7))
    ax.annotate(f"Trans +1 (P1)\n+ cis Arch 3\n(COM {p1_rec['com_dist_vs_cis']:.2f} Å)", xy=(P_B[45, 0], P_B[45, 1]),
                xytext=(xmax - 1.5, ymin + 2.5), fontsize=5.5, color=PAL["plus1"], fontweight="bold",
                ha="right", va="bottom",
                arrowprops=dict(arrowstyle="-", color=PAL["plus1"], lw=0.6))

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("3D superposition: Phage P1 trans vs. cis Arch 3", fontsize=6.6, fontweight="bold", pad=3)
