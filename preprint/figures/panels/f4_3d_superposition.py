"""Figure 4, Panel c: py2Dmol-style flat-shaded Richardson cartoon of Phage P1 trans AF3 complex onto cis A0A6N9NLZ4 & 1MHT DNA."""
import json
import re
import numpy as np
from style import DATA_DIR, PAL
from panels.struct3d import (
    _catmull,
    _depth_shade,
    build_cartoon_primitives,
    draw_outlined_atom,
    kabsch,
    parse_cif,
    render_cartoon,
)


def _landmarks(chain_d):
    ids = sorted(r for r in chain_d if chain_d[r]["ca"] is not None)
    seq = "".join(chain_d[r]["resn"] for r in ids)
    lm = {}
    m4 = re.search(r"[GASCVITLFY][FTSAVILYCMWP]PC[QSNTG]", seq)
    if m4:
        lm["M4"] = ids[m4.start() + 2]
    m6 = re.search(r"[IVLMFWYC]ENV", seq)
    if m6:
        lm["M6"] = ids[m6.start() + 1]
    m10 = re.search(r"GN[SA][VWILMT][ATVIL]", seq)
    if m10:
        lm["M10"] = ids[m10.start()]
    return lm


def _build_dna_primitives(P_c, P_d, zmin, zmax, backbone_col="#7C736B", rung_col="#A89F96"):
    """Build depth-sorted py2Dmol-style flat-shaded ribbons for DNA phosphate strands + base-pair rungs."""
    prims = []
    zspan = max(1e-3, zmax - zmin)

    def add_tube(pts3, width, col, outline="#1C1C1C", lw=0.42, alpha=0.84, subdiv=3):
        curve = _catmull(pts3, subdiv=subdiv)
        m = len(curve)
        if m < 2:
            return
        xy = curve[:, :2]
        tang = np.zeros_like(xy)
        tang[1:-1] = xy[2:] - xy[:-2]
        tang[0] = xy[1] - xy[0]
        tang[-1] = xy[-1] - xy[-2]
        tang = tang / np.maximum(np.linalg.norm(tang, axis=1, keepdims=True), 1e-6)
        perp = np.column_stack([-tang[:, 1], tang[:, 0]])
        left = xy + perp * width
        right = xy - perp * width
        for k in range(m - 1):
            zm = 0.5 * (curve[k, 2] + curve[k + 1, 2])
            zn = (zm - zmin) / zspan
            fc = _depth_shade(col, zn, face_dot=0.7, fog=0.25)
            oc = _depth_shade(outline, zn, face_dot=1.0, fog=0.18) if outline else None
            prims.append({
                "z": zm,
                "poly": np.array([left[k], left[k + 1], right[k + 1], right[k]]),
                "edge_L": (left[k], left[k + 1]),
                "edge_R": (right[k], right[k + 1]),
                "cap0": (left[0], right[0]) if k == 0 else None,
                "cap1": (left[-1], right[-1]) if k == m - 2 else None,
                "color": fc,
                "outline": oc,
                "lw": lw,
                "alpha": alpha,
            })

    add_tube(P_c, 0.72, backbone_col, outline="#2A2623", lw=0.45, alpha=0.84, subdiv=3)
    add_tube(P_d, 0.72, backbone_col, outline="#2A2623", lw=0.45, alpha=0.84, subdiv=3)

    n_rungs = min(len(P_c), len(P_d))
    for k in range(n_rungs):
        p1, p2 = P_c[k], P_d[n_rungs - 1 - k]
        add_tube(np.vstack([p1, p2]), 0.28, rung_col, outline="#5A524C", lw=0.25, alpha=0.68, subdiv=1)

    return prims


def draw(fig):
    ax = fig.add_subplot(111)

    # Load real coordinates: cis A0A6N9NLZ4, trans Phage P1 64aa AF3, and 1MHT DNA
    a0_ch = parse_cif(DATA_DIR / "structures/AF-A0A6N9NLZ4-F1-model_v6.cif")["A"]
    p1_all = parse_cif(DATA_DIR / "structures/af3_models/fold_phagep1_arch1_dmt_trans_64aa_core_model_0.cif")
    p1A, p1B = p1_all["A"], p1_all["B"]
    mht_all = parse_cif(DATA_DIR / "structures/1MHT.cif")
    mhtA = mht_all["A"]

    audit = json.loads((DATA_DIR / "structures/trans_vs_cis_3d_audit.json").read_text())
    p1_rec = [r for r in audit["trans_complexes"] if r["name"] == "phagep1_arch1_dmt_trans_64aa_core"][0]

    lm0 = _landmarks(a0_ch)
    lm1 = _landmarks(p1A)
    lmH = _landmarks(mhtA)

    def align_to_a0(chain_d, lm):
        # Strictly rigid Rossmann core anchors: M4[-4..+1] EXCLUDES mobile catalytic loop M4[+2..+8]
        P, Q = [], []
        for k, a, b in [("M4", -4, 2), ("M6", -3, 7), ("M10", -1, 7)]:
            i1, i2 = lm[k], lm0[k]
            for off in range(a, b):
                P.append(chain_d[i1 + off]["ca"])
                Q.append(a0_ch[i2 + off]["ca"])
        return kabsch(np.array(P), np.array(Q))

    R1, t1 = align_to_a0(p1A, lm1)
    RH, tH = align_to_a0(mhtA, lmH)
    xf1 = lambda pt: pt @ R1.T + t1
    xfH = lambda pt: pt @ RH.T + tH

    # Extract 1MHT DNA phosphate strands (chains B and C)
    dna_chains = []
    for ch_id in sorted(mht_all.keys()):
        if ch_id == "A":
            continue
        pts = [xfH(mht_all[ch_id][r]["atoms"]["P"]) for r in sorted(mht_all[ch_id]) if "P" in mht_all[ch_id][r]["atoms"]]
        if len(pts) >= 4:
            dna_chains.append(np.array(pts))
    dna_C, dna_D = dna_chains[0], dna_chains[1]

    # Key landmark coordinates in A0A6N9NLZ4 reference frame
    cys_pt = xf1(p1A[78]["ca"])    # Motif IV Cys78
    r385_pt = xf1(p1A[385]["ca"])  # Motif IX R385
    f30_pt = xf1(p1B[30]["ca"])    # RYFLS F30
    dna_mid = 0.5 * (dna_C.mean(0) + dna_D.mean(0))

    # Clean camera basis: x-axis Cys78 -> F30, y-axis toward DNA duplex, slight 3D tilt
    vx = f30_pt - cys_pt
    vx /= np.linalg.norm(vx)
    vy = dna_mid - 0.5 * (cys_pt + f30_pt)
    vy = vy - np.dot(vy, vx) * vx
    vy /= np.linalg.norm(vy)
    vz = np.cross(vx, vy)
    ang = np.deg2rad(16)
    vy2 = np.cos(ang) * vy + np.sin(ang) * vz
    vz2 = -np.sin(ang) * vy + np.cos(ang) * vz
    basis = np.vstack([vx, vy2, vz2])

    p1A_res = sorted(r for r in p1A if p1A[r]["ca"] is not None)
    p1B_res = sorted(r for r in p1B if p1B[r]["ca"] is not None)
    cis_res = [r for r in range(388, 447) if r in a0_ch and a0_ch[r]["ca"] is not None]

    origin = np.mean([xf1(p1A[r]["ca"]) for r in p1A_res], axis=0)
    proj_p1 = lambda pt: (xf1(pt) - origin) @ basis.T
    proj_a0 = lambda pt: (pt - origin) @ basis.T
    proj_arr = lambda arr: (arr - origin) @ basis.T

    P_dC = proj_arr(dna_C)
    P_dD = proj_arr(dna_D)

    # Color function for full Phage P1 C5-MTase (Chain A: soft slate background, bold Motif VIII & IX)
    def p1A_col(r):
        if 384 <= r <= 398:
            return PAL["motif9"]
        if 201 <= r <= 215:
            return PAL["frame0"]
        return "#C7D6E7"

    # Build py2Dmol flat-shaded Richardson cartoon primitives
    prims = []
    # Background C5-MTase catalytic core & TRD (soft outlined ribbons so secondary structure reads cleanly)
    core_bg = [r for r in p1A_res if not (201 <= r <= 215 or 384 <= r <= 398)]
    # Render Chain A in contiguous chunks so Motif VIII and Motif IX get crisp dark outlines while core stays subtle
    for seg_rs, is_hi in [
        ([r for r in p1A_res if r <= 201], False),
        ([r for r in p1A_res if 201 <= r <= 215], True),
        ([r for r in p1A_res if 215 <= r <= 384], False),
        ([r for r in p1A_res if 384 <= r <= 398], True),
        ([r for r in p1A_res if r >= 398], False),
    ]:
        if len(seg_rs) < 2:
            continue
        prims.extend(build_cartoon_primitives(
            p1A, seg_rs, proj_p1, p1A_col,
            helix_w=1.22 if is_hi else 1.04,
            strand_w=1.14 if is_hi else 0.96,
            loop_w=0.44 if is_hi else 0.34,
            alpha=0.98 if is_hi else 0.62,
            outline="#1C1C1C" if is_hi else "#788DA8",
            outline_lw=0.48 if is_hi else 0.28,
            subdiv=4 if is_hi else 3,
            fog=0.18 if is_hi else 0.42,
            style="ribbon",
        ))

    # Cis A0A6N9NLZ4 CSLSDILE/KYFLS domain overlay (petrol teal flat-shaded cartoon)
    prims.extend(build_cartoon_primitives(
        a0_ch, cis_res, proj_a0, lambda r: PAL["arch3"],
        helix_w=1.02, strand_w=0.98, loop_w=0.38,
        alpha=0.76, outline="#134E5C", outline_lw=0.36,
        subdiv=4, fog=0.18, style="ribbon",
    ))

    # Trans Phage P1 64-aa +1 microprotein (bold matte terracotta hero cartoon)
    prims.extend(build_cartoon_primitives(
        p1B, p1B_res, proj_p1, lambda r: PAL["plus1"],
        helix_w=1.26, strand_w=1.18, loop_w=0.46,
        alpha=0.98, outline="#1C1C1C", outline_lw=0.52,
        subdiv=4, fog=0.14, style="ribbon",
    ))

    # DNA duplex phosphate strands + rungs
    all_z = [p["z"] for p in prims]
    zmin_g, zmax_g = min(all_z), max(all_z)
    prims.extend(_build_dna_primitives(P_dC, P_dD, zmin_g, zmax_g))

    render_cartoon(ax, prims)

    # Highlight key atomic lock spheres & dashed contact
    p_cys = proj_p1(p1A[78]["ca"])
    p_r385 = proj_p1(p1A[385]["ca"])
    p_f30 = proj_p1(p1B[30]["ca"])
    p_dna_lbl = P_dC[np.argmax(P_dC[:, 1])]
    p_p1b_lbl = proj_p1(p1B[46]["ca"])

    ax.plot([p_r385[0], p_f30[0]], [p_r385[1], p_f30[1]], ls="--", lw=1.25, color=PAL["ink"], zorder=9.5)
    draw_outlined_atom(ax, p_cys, 1.35, PAL["motif4"], zord=10.0)
    draw_outlined_atom(ax, p_r385, 1.45, PAL["motif9"], zord=10.2)
    draw_outlined_atom(ax, p_f30, 1.45, PAL["plus1"], zord=10.4)

    # Compute exact unclipped bounds across ALL full chains + DNA (zero top/bottom/side cutoffs!)
    all_pts = np.vstack([
        np.array([proj_p1(p1A[r]["ca"])[:2] for r in p1A_res]),
        np.array([proj_p1(p1B[r]["ca"])[:2] for r in p1B_res]),
        np.array([proj_a0(a0_ch[r]["ca"])[:2] for r in cis_res]),
        P_dC[:, :2],
        P_dD[:, :2],
    ])
    x0, y0 = all_pts.min(0)
    x1, y1 = all_pts.max(0)
    xmin, xmax = x0 - 11.5, x1 + 13.5
    ymin, ymax = y0 - 9.0, y1 + 8.5

    bbox_kw = dict(fc="white", ec="none", pad=0.45, alpha=0.88)

    # Clean sparse callouts anchored safely in the margins
    ax.annotate(
        f"DNA duplex\n(1MHT, {p1_rec['dna_min_dist']:.1f} Å Cα)",
        xy=(p_dna_lbl[0], p_dna_lbl[1]),
        xytext=(xmin + 1.5, ymax - 1.5),
        fontsize=5.2, color=PAL["dna"], fontweight="bold", ha="left", va="top",
        bbox=bbox_kw,
        arrowprops=dict(arrowstyle="-", color=PAL["dna"], lw=0.6, shrinkB=2),
        zorder=11,
    )
    ax.annotate(
        f"Cys78\n(Motif IV, {p1_rec['cys_min_dist_ca_mean']:.0f} Å)",
        xy=(p_cys[0], p_cys[1]),
        xytext=(xmin + 1.5, ymin + 5.0),
        fontsize=5.2, color=PAL["motif4"], fontweight="bold", ha="left", va="bottom",
        bbox=bbox_kw,
        arrowprops=dict(arrowstyle="-", color=PAL["motif4"], lw=0.6, shrinkB=3),
        zorder=11,
    )
    ax.annotate(
        f"Motif IX R385\n↔ F30 ({p1_rec['r1_kyf_mean']:.2f} Å)",
        xy=(0.5 * (p_r385[0] + p_f30[0]), 0.5 * (p_r385[1] + p_f30[1])),
        xytext=(xmax - 1.5, ymax - 1.5),
        fontsize=5.3, color=PAL["ink"], fontweight="bold", ha="right", va="top",
        bbox=bbox_kw,
        arrowprops=dict(arrowstyle="-|>", color=PAL["ink"], lw=0.7, mutation_scale=7, shrinkB=3),
        zorder=11,
    )
    ax.annotate(
        f"Trans +1 (P1)\n+ cis Arch 3\n(COM {p1_rec['com_dist_vs_cis']:.2f} Å)",
        xy=(p_p1b_lbl[0], p_p1b_lbl[1]),
        xytext=(xmax - 1.5, ymin + 5.0),
        fontsize=5.2, color=PAL["plus1"], fontweight="bold", ha="right", va="bottom",
        bbox=bbox_kw,
        arrowprops=dict(arrowstyle="-", color=PAL["plus1"], lw=0.6, shrinkB=2),
        zorder=11,
    )

    # Compact legend matching TMEM168 Fig. 1b aesthetic
    ax.plot([], [], color="#A6BDD8", lw=2.6, label="P1 C5-MTase (1–473)")
    ax.plot([], [], color=PAL["motif9"], lw=2.6, label="Motif IX")
    ax.plot([], [], color=PAL["plus1"], lw=2.6, label="Trans +1 (64 aa)")
    ax.plot([], [], color=PAL["arch3"], lw=2.4, label="Cis Arch 3")
    ax.legend(loc="lower center", bbox_to_anchor=(0.48, -0.04), ncol=2, frameon=False,
              fontsize=4.7, handlelength=1.0, columnspacing=0.9, labelspacing=0.18, borderpad=0.1)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin - 5.5, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("3D superposition: Phage P1 trans vs. cis Arch 3", fontsize=6.6, fontweight="bold", pad=3)
