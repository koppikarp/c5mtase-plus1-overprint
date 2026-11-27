"""py2Dmol / Richardson-inspired flat-shaded, black-outlined pure-vector 3D cartoon renderer.

Computes real secondary structure (alpha-helices, beta-strands with arrowheads, smooth loops)
from atomic coordinates, renders depth-sorted outlined polygons/ribbons with subtle
front-to-back depth cueing, and outputs 100% vector Matplotlib/Fignimble geometry.
"""
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, Polygon

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_pdb(path, chain_filter=None):
    res = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM  "):
                continue
            ch = line[21].strip() or "A"
            if chain_filter and ch != chain_filter:
                continue
            an = line[12:16].strip()
            comp = line[17:20].strip()
            resi = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            rd = res.setdefault(resi, {"resn": AA3.get(comp, "X"), "ca": None, "atoms": {}})
            rd["atoms"][an] = xyz
            if an == "CA":
                rd["ca"] = xyz
    return res


def parse_cif(path):
    chains, cols, in_atom = {}, [], False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1]); in_atom = True; continue
            if in_atom and (s.startswith("ATOM ") or s.startswith("HETATM ")):
                p = s.split()
                if len(p) < len(cols):
                    continue
                ch = p[cols.index("auth_asym_id")]
                seq = p[cols.index("label_seq_id")]
                if not seq.isdigit():
                    continue
                r = int(seq)
                an = p[cols.index("label_atom_id")].strip('"')
                comp = p[cols.index("label_comp_id")]
                xyz = np.array([float(p[cols.index("Cartn_x")]), float(p[cols.index("Cartn_y")]), float(p[cols.index("Cartn_z")])])
                rd = chains.setdefault(ch, {}).setdefault(r, {"resn": AA3.get(comp, "X"), "ca": None, "atoms": {}})
                rd["atoms"][an] = xyz
                if an == "CA":
                    rd["ca"] = xyz
            elif in_atom and s.startswith("#"):
                in_atom = False
    return chains


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cQ - cP @ R.T


def assign_ss(chain_d, res_list):
    """Fast backbone-geometry secondary structure assignment (H=helix, E=strand, C=coil)."""
    n = len(res_list)
    cas = np.array([chain_d[r]["ca"] for r in res_list])
    ss = ["C"] * n
    # Alpha-helix: d(i,i+3) in [4.5, 6.0] and d(i,i+4) in [5.4, 7.0]
    for i in range(n - 4):
        d2 = np.linalg.norm(cas[i + 2] - cas[i])
        d3 = np.linalg.norm(cas[i + 3] - cas[i])
        d4 = np.linalg.norm(cas[i + 4] - cas[i])
        if 4.6 <= d2 <= 6.2 and 4.5 <= d3 <= 6.0 and 5.3 <= d4 <= 7.0:
            for k in range(i, i + 4):
                ss[k] = "H"
    # Beta-strand: extended CA(i,i+2) >= 6.1 and CA(i,i+3) >= 8.9
    for i in range(n - 3):
        if ss[i] == "H":
            continue
        d2 = np.linalg.norm(cas[i + 2] - cas[i])
        d3 = np.linalg.norm(cas[i + 3] - cas[i])
        if d2 >= 6.15 and d3 >= 8.95:
            for k in range(i, i + 3):
                if ss[k] == "C":
                    ss[k] = "E"
    # Clean isolated 1-2 residue blips
    for i in range(1, n - 1):
        if ss[i] != ss[i - 1] and ss[i] != ss[i + 1]:
            ss[i] = ss[i - 1]
    return ss


def _catmull(pts, subdiv=4):
    pts = np.asarray(pts, dtype=float)
    if len(pts) < 3:
        return pts
    out = []
    p = np.vstack([pts[0], pts, pts[-1]])
    for i in range(1, len(p) - 2):
        p0, p1, p2, p3 = p[i - 1], p[i], p[i + 1], p[i + 2]
        for s in range(subdiv):
            t = s / subdiv
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(pts[-1])
    return np.array(out)


def _blend_rgb(hex_col, target_rgb, frac):
    c = np.array(to_rgb(hex_col))
    t = np.array(target_rgb)
    return tuple(np.clip((1.0 - frac) * c + frac * t, 0.0, 1.0))


def _depth_shade(hex_col, z_norm, face_dot=0.6, fog=0.36):
    """py2Dmol / Richardson flat shading: subtle depth fog + face brightness."""
    # z_norm in [0, 1] where 1 = front, 0 = back
    c = np.array(to_rgb(hex_col))
    # Subtle facet brightening/darkening
    shade = 0.86 + 0.24 * np.clip(abs(face_dot), 0.0, 1.0)
    c = np.clip(c * shade, 0.0, 1.0)
    # Depth fog toward crisp white for distant segments
    fog_amt = fog * (1.0 - np.clip(z_norm, 0.0, 1.0))
    c = (1.0 - fog_amt) * c + fog_amt * np.array([1.0, 1.0, 1.0])
    return tuple(np.clip(c, 0.0, 1.0))


def build_cartoon_primitives(chain_d, res_list, proj_fn, color_fn,
                             helix_w=1.15, strand_w=1.05, loop_w=0.42,
                             alpha=1.0, outline="#1C1C1C", outline_lw=0.45,
                             subdiv=4, ss_override=None, fog=0.34, style="ribbon"):
    """
    Build depth-sorted py2Dmol-style flat-shaded cartoon primitives:
    returns list of dicts: {z, poly (Nx2), color, edge, lw, alpha}
    style='ribbon': Richardson coiled ribbons for helices, arrowheaded bands for strands, tubes for loops
    style='cylinder': Smooth rounded cylindrical pipes for TM helices (extra clean for 12-TM bundles!)
    """
    n = len(res_list)
    if n < 2:
        return []
    ss = ss_override if ss_override is not None else assign_ss(chain_d, res_list)
    pts3 = np.array([proj_fn(chain_d[r]["ca"]) for r in res_list])
    cols = [color_fn(r) for r in res_list]

    # Segment into contiguous runs of same SS
    runs = []
    start = 0
    for i in range(1, n):
        if ss[i] != ss[start]:
            runs.append((ss[start], start, i))
            start = i
    runs.append((ss[start], start, n))

    # Global z range for depth cueing
    zmin, zmax = pts3[:, 2].min(), pts3[:, 2].max()
    zspan = max(1e-3, zmax - zmin)

    prims = []

    def add_quad_ribbon(curve, widths, col_list, is_arrow_end=False):
        m = len(curve)
        if m < 2:
            return
        # Compute 2D perpendicular normal smoothly along screen projection
        xy = curve[:, :2]
        tang = np.zeros_like(xy)
        tang[1:-1] = xy[2:] - xy[:-2]
        tang[0] = xy[1] - xy[0]
        tang[-1] = xy[-1] - xy[-2]
        tlen = np.linalg.norm(tang, axis=1, keepdims=True)
        tang = tang / np.maximum(tlen, 1e-6)
        perp = np.column_stack([-tang[:, 1], tang[:, 0]])

        # 3D tangent for face shading
        t3 = np.zeros_like(curve)
        t3[1:-1] = curve[2:] - curve[:-2]
        t3[0] = curve[1] - curve[0]
        t3[-1] = curve[-1] - curve[-2]
        t3 = t3 / np.maximum(np.linalg.norm(t3, axis=1, keepdims=True), 1e-6)

        w_arr = np.asarray(widths, dtype=float)
        if is_arrow_end and m >= 5:
            # Taper last 4 segments into a crisp Richardson beta-strand arrowhead!
            ak = min(5, m - 1)
            w_arr[-ak] *= 1.65
            for s in range(1, ak):
                w_arr[-ak + s] = w_arr[-ak] * (1.0 - s / ak)

        left = xy + perp * w_arr[:, None]
        right = xy - perp * w_arr[:, None]

        # Build overlapping micro-quads so depth-sorting handles helical crossings cleanly!
        for k in range(m - 1):
            zm = 0.5 * (curve[k, 2] + curve[k + 1, 2])
            zn = (zm - zmin) / zspan
            # Slight overlap along tangent eliminates hairline gaps between quads
            poly = np.array([left[k], left[k + 1], right[k + 1], right[k]])
            face_z = 1.0 - abs(t3[k, 2])
            fc = _depth_shade(col_list[k], zn, face_dot=face_z, fog=fog)
            # Draw side edges ONLY (left edge & right edge) via custom edge lines, plus seamless fill!
            prims.append({
                "z": zm,
                "poly": poly,
                "edge_L": (left[k], left[k + 1]),
                "edge_R": (right[k], right[k + 1]),
                "cap0": (left[0], right[0]) if k == 0 else None,
                "cap1": (left[-1], right[-1]) if k == m - 2 else None,
                "color": fc,
                "outline": _depth_shade(outline, zn, face_dot=1.0, fog=fog * 0.55) if outline else None,
                "lw": outline_lw,
                "alpha": alpha,
            })

    for stype, a, b in runs:
        # Extend 1 residue on each side for seamless joins
        ia = max(0, a - 1)
        ib = min(n, b + 1)
        sub_pts = pts3[ia:ib]
        sub_cols = cols[ia:ib]
        if len(sub_pts) < 2:
            continue

        if stype == "H" and style == "cylinder" and (b - a) >= 5:
            # Straight/slightly curved helical axis cylinder (ultra-clean for 12-TM bundles!)
            ax_pts = []
            ax_cols = []
            for k in range(a, b):
                lo, hi = max(a, k - 2), min(b, k + 3)
                ax_pts.append(pts3[lo:hi].mean(0))
                ax_cols.append(cols[k])
            ax_pts = np.array(ax_pts)
            curve = _catmull(ax_pts, subdiv=3)
            c_idx = np.clip(np.linspace(0, len(ax_cols) - 1, len(curve)).round().astype(int), 0, len(ax_cols) - 1)
            c_list = [ax_cols[j] for j in c_idx]
            add_quad_ribbon(curve, np.full(len(curve), helix_w * 1.35), c_list, is_arrow_end=False)
        elif stype == "H":
            # Smooth Richardson coiled helix ribbon
            curve = _catmull(sub_pts, subdiv=subdiv)
            c_idx = np.clip(np.linspace(0, len(sub_cols) - 1, len(curve)).round().astype(int), 0, len(sub_cols) - 1)
            c_list = [sub_cols[j] for j in c_idx]
            add_quad_ribbon(curve, np.full(len(curve), helix_w), c_list, is_arrow_end=False)
        elif stype == "E" and (b - a) >= 3:
            # Smooth flattened beta-strand with Richardson arrowhead!
            sm_raw = sub_pts.copy()
            # Gently smooth beta zig-zag into a crisp flat strand
            if len(sm_raw) >= 3:
                sm_raw[1:-1] = 0.25 * sub_pts[:-2] + 0.50 * sub_pts[1:-1] + 0.25 * sub_pts[2:]
            curve = _catmull(sm_raw, subdiv=subdiv)
            c_idx = np.clip(np.linspace(0, len(sub_cols) - 1, len(curve)).round().astype(int), 0, len(sub_cols) - 1)
            c_list = [sub_cols[j] for j in c_idx]
            add_quad_ribbon(curve, np.full(len(curve), strand_w), c_list, is_arrow_end=True)
        else:
            # Smooth rounded loop tube
            curve = _catmull(sub_pts, subdiv=max(2, subdiv - 1))
            c_idx = np.clip(np.linspace(0, len(sub_cols) - 1, len(curve)).round().astype(int), 0, len(sub_cols) - 1)
            c_list = [sub_cols[j] for j in c_idx]
            add_quad_ribbon(curve, np.full(len(curve), loop_w), c_list, is_arrow_end=False)

    return prims


def render_cartoon(ax, prims):
    """Draw depth-sorted py2Dmol flat-shaded polygons with crisp continuous side outlines."""
    prims = sorted(prims, key=lambda p: p["z"])
    for idx, p in enumerate(prims):
        zord = 2.0 + 3.0 * (idx / max(1, len(prims)))
        # Fill quad (tiny facecolor-matched edgecolor seals sub-pixel seams between adjacent quads)
        ax.add_patch(Polygon(p["poly"], closed=True, facecolor=p["color"],
                             edgecolor=p["color"], linewidth=0.22, alpha=p["alpha"], zorder=zord))
        if p["outline"]:
            el, er = p["edge_L"], p["edge_R"]
            ax.plot([el[0][0], el[1][0]], [el[0][1], el[1][1]],
                    color=p["outline"], lw=p["lw"], alpha=p["alpha"], solid_capstyle="round", zorder=zord + 0.001)
            ax.plot([er[0][0], er[1][0]], [er[0][1], er[1][1]],
                    color=p["outline"], lw=p["lw"], alpha=p["alpha"], solid_capstyle="round", zorder=zord + 0.001)
            if p["cap0"] is not None:
                c0 = p["cap0"]
                ax.plot([c0[0][0], c0[1][0]], [c0[0][1], c0[1][1]],
                        color=p["outline"], lw=p["lw"], alpha=p["alpha"], solid_capstyle="round", zorder=zord + 0.001)
            if p["cap1"] is not None:
                c1 = p["cap1"]
                ax.plot([c1[0][0], c1[1][0]], [c1[0][1], c1[1][1]],
                        color=p["outline"], lw=p["lw"], alpha=p["alpha"], solid_capstyle="round", zorder=zord + 0.001)


def draw_outlined_stick(ax, p1, p2, col, lw=2.3, outline="#1C1C1C", zord=8.0, ls="-", alpha=0.96):
    """py2Dmol-style black-outlined covalent bond stick."""
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=outline, lw=lw + 1.1, ls=ls, alpha=alpha, solid_capstyle="round", zorder=zord)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=col, lw=lw, ls=ls, alpha=alpha, solid_capstyle="round", zorder=zord + 0.05)


def draw_outlined_atom(ax, pt, radius, col, outline="#1C1C1C", zord=9.0):
    """py2Dmol-style flat-shaded sphere with crisp dark outline and subtle highlight."""
    ax.add_patch(Circle((pt[0], pt[1]), radius, facecolor=col, edgecolor=outline, linewidth=0.65, zorder=zord))
    ax.add_patch(Circle((pt[0] - radius * 0.28, pt[1] + radius * 0.28), radius * 0.24,
                        facecolor="white", edgecolor="none", alpha=0.55, zorder=zord + 0.05))
