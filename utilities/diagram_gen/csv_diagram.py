#!/usr/bin/env python3
"""
csv_diagram.py — parse a SeQuant CSV-CCSD(-DF) term string (format of pao-ccsd-df.md)
and render the corresponding CSV-CC Goldstone diagram with matplotlib.

Conventions implemented (see csv-cc-diagram-rules.md in this directory for citations):
  * hole line  = downward arrow (blue), occupied labels i_n
  * particle line = upward arrow (red), PNO labels a_n^{dom}
  * two-electron interaction = ONE vertex = two half-vertices (dots) joined by a
    dashed line; under DF the dashed line carries the auxiliary label K
  * one-body (Fock) interaction = a dot joined by a dashed line to a "×" marker
    (A2: "the point of action ... marked by the interaction line —×")
  * bra index = outgoing (creation); ket index = incoming (annihilation);
    every half-vertex has exactly 1 incoming + 1 outgoing line
  * fermion flow of a line: from the terminal node where its index sits in a BRA
    slot to the terminal node where it sits in a KET slot
  * amplitude vertex = solid horizontal bar at the bottom, label \\hat{t}_{occ}^{csv-virts}
  * CSV overlay: PNO domain as superscript a^{ij}; PAO->PNO transform on a g/f leg
    annotated "C"; a domain-changing bridge C·s·C annotated "CsC"; DF split of the
    interaction annotated "K" on the dashed line
  * external lines = indices listed in the projector Ŝ{occ;virt}; drawn to the top

CSV OVERLAY KNOB (rule B0, skeleton invariance)
  draw(..., show_csv=False) removes every CSV/DF annotation — domain superscripts,
  C / CsC / s boxes, the K label and the g_L/g_R half-vertex tags — leaving the
  plain canonical CC Goldstone diagram.  show_aux=<bool> overrides the K/g_L/g_R
  part alone (it defaults to show_csv).

INPUT
  Any SeQuant term string is accepted — inline, from a file, or on stdin.  The
  ':N-C-S' suffixes, markdown list numbering and backticks are stripped, so a
  line copied straight out of pao-ccsd*.md works as-is:

    python3 csv_diagram.py '<term>' out.png                # inline
    python3 csv_diagram.py --file term.txt out.png         # from a file
    sed -n '281p' ../pao-ccsdt-df.md | python3 csv_diagram.py - out.png
    python3 csv_diagram.py --file many.txt outdir/         # N terms -> outdir/
    python3 render_all.py --src ../pao-ccsdtq-df.md --only R4   # whole section

LAYOUT
  * geometry is computed analytically (no renderer round-trip): text extents come
    from a scratch canvas in DATA units, so the figure is built at exactly 1 data
    unit = UNIT inches and every label is placed by a collision-avoiding greedy
    search against vertices, bars, other labels and the fermion lines themselves.
  * amplitude bars, half-vertices and external stubs are all centred on a common
    axis, so mirror-image terms render as mirror images.
  * two fermion lines joining the SAME pair of endpoints (e.g. the particle+hole
    pair of a T1 hanging off a single interaction dot) are bowed apart into a
    slim pointy ellipse (the canonical Goldstone bubble); lines that merely run
    together (a shared apex, an external stub grazing a half-vertex) are fanned
    apart by the same mechanism, so no two lines are ever drawn on top of one
    another.
  * C / CsC / s ride ON their fermion line and slide along it to dodge obstacles.

VERIFIED ON
  pao-ccsd-df.md (86), pao-ccsdt-df.md (201) and pao-ccsdtq-df.md (373 terms,
  T1..T4 amplitude bars) — see audit_diagrams.py: no label/label, label/vertex or
  line/line overlaps in any of them, in either variant.

LIMITATIONS (not yet handled):
  * line-crossing minimisation is only heuristic (amplitude/vertex ordering by
    mean partner position); very crowded terms may still cross lines, which is
    legal but not always pretty.  Lines never coincide and labels never overlap,
    but a label box may sit on top of one foreign line (1/86 CCSD, 8/201 CCSDT
    and 32/373 CCSDTQ terms); it keeps an opaque background and stays readable.
  * assumes at most ONE interaction (<=2 g half-vertices sharing one Κ, or one f).
    True for every CCSD/CCSDT/CCSDTQ residual term in the files above.
  * does not compute the diagram sign/weight (h, loops, equivalent lines);
    it draws and annotates only. The printed coefficient is taken verbatim.
  * an s{μ̃;μ̃} that is NOT flanked by two C's (never occurs in the files) would
    be annotated as a bare "s".
"""

import argparse
import os
import re
import sys
import unicodedata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HOLE = "#1f4fd8"   # blue, arrow down
PART = "#c0281f"   # red,  arrow up
INT  = "#555555"   # gray dashed interaction

Y_T, Y_HV, Y_EXT = 1.0, 3.9, 6.5

UNIT = 1.0           # inches per data unit  => 72 pt per data unit
FS_IDX, FS_ANN, FS_T, FS_V, FS_TITLE, FS_LEG = 12, 11, 11, 12, 13, 10

APEX_STEP = 0.85     # spacing of the excitation apexes on one amplitude bar
BAR_HALF  = 0.30     # bar overhang beyond the outermost apex
T_GAP     = 0.75     # clear space between neighbouring amplitude bars/labels
HV_GAP    = 3.4      # minimum separation of the two interaction half-vertices
EXT_BIAS  = 0.14     # tie-break only: holes lean left of their anchor, particles right
EXT_GAP   = 0.55     # clear space between neighbouring external labels
STUB      = 1.15     # length of the Fock  ---x  stub

LAST_AUDIT = None    # geometry of the most recent draw() — used by audit_diagrams.py

# ----------------------------------------------------------------------------
# parsing
# ----------------------------------------------------------------------------

TENSOR_RE = re.compile(r"([^\s*(){}]+)\{([^{}]*)\}(?::[A-ZN\-]+)?")
COEFF_RE  = re.compile(r"^\s*([+-]?\d+(?:/\d+)?)\s+")
MD_ITEM   = re.compile(r"^\s*\d+\.\s+(.*)$")     # "12. `<term>`" markdown list line


def clean_term(s):
    """Accept a raw SeQuant term as pasted from anywhere: a bare expression, a
    markdown list item (``12. `...` ``), backtick-quoted, or with stray
    surrounding whitespace / a leading '+'."""
    s = s.strip()
    m = MD_ITEM.match(s)
    if m and "{" in m.group(1):
        s = m.group(1).strip()
    s = s.strip("`").strip()
    if s.startswith("+") and not s[1:2].isdigit():
        s = s[1:].strip()
    return s


def _norm(s):
    return unicodedata.normalize("NFC", s)


def kind_of(idx):
    """occ / pno / pao / aux from the index name."""
    if idx.startswith("i"):
        return "occ"
    if idx.startswith("a"):
        return "pno"
    if "μ" in idx:
        return "pao"
    if "Κ" in idx or "K" in idx:
        return "aux"
    raise ValueError(f"unknown index kind: {idx!r}")


def split_outside(s, sep):
    """split on sep, but not inside <...> (PNO domain brackets)."""
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(cur)); cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


def parse_index(tok):
    """'a_3<i_1,i_2>' -> ('a_3', ['i_1','i_2']);  'i_1' -> ('i_1', [])"""
    tok = tok.strip()
    m = re.match(r"([^<>]+)(?:<([^<>]*)>)?$", tok)
    name = m.group(1).strip()
    dom = [d.strip() for d in m.group(2).split(",")] if m.group(2) else []
    return name, dom


class Tensor:
    def __init__(self, label, slots):
        self.label = label                      # 'g','f','t','s','C','S^'
        self.bra, self.ket = slots[0], slots[1]  # lists of (name, dom)
        self.aux = slots[2] if len(slots) > 2 else []

    def __repr__(self):
        return f"{self.label}{{{self.bra};{self.ket};{self.aux}}}"


def parse_term(term):
    """Return (coeff:str, tensors:[Tensor]). Parens are irrelevant (associativity)."""
    term = _norm(clean_term(term))
    if "{" not in term:
        raise ValueError(f"no tensor found in input: {term[:80]!r}")
    m = COEFF_RE.match(term)
    coeff = m.group(1) if m else "1"
    tensors = []
    for lab, body in TENSOR_RE.findall(term):
        lab = _norm(lab)
        if lab in ("Ŝ", "Ŝ", "S^"):
            lab = "S^"
        slots = [[parse_index(t) for t in split_outside(part, ",") if t.strip()]
                 for part in body.split(";")]
        while len(slots) < 2:
            slots.append([])
        tensors.append(Tensor(lab, slots))
    return coeff, tensors


# ----------------------------------------------------------------------------
# graph construction: terminals + C/s connector chains -> fermion lines
# ----------------------------------------------------------------------------

TERMINALS = {"g", "f", "t", "S^"}


class Node:
    """A drawable vertex: half-vertex ('hv'), amplitude ('t'), or external ('ext')."""
    def __init__(self, kind, tag, tensor=None):
        self.kind, self.tag, self.tensor = kind, tag, tensor
        self.x = self.y = 0.0

    def __repr__(self):
        return f"<{self.kind}:{self.tag}>"


class FermionLine:
    def __init__(self, src, dst, src_idx, dst_idx, nC, nS):
        self.src, self.dst = src, dst           # Node (bra end -> ket end)
        self.src_idx, self.dst_idx = src_idx, dst_idx  # (name, dom) at each end
        self.nC, self.nS = nC, nS               # of C / s tensors on the chain
        # hole vs particle from terminal index kinds
        kinds = {kind_of(src_idx[0]), kind_of(dst_idx[0])}
        self.is_hole = "occ" in kinds

    @property
    def annotation(self):
        if self.nS:
            return "CsC" if self.nC >= 2 else "sC" if self.nC else "s"
        return "C" if self.nC == 1 else ("CC" if self.nC else None)


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def build_graph(tensors):
    """Returns (nodes, lines, hvs, hv_pair_or_None, externals dict)."""
    uf = UF()
    conn_count = {}          # root -> [nC, nS] (filled after unions)
    connectors = []          # (label, p, q)
    terminal_occ = {}        # index name -> list of (node, 'bra'|'ket', (name,dom))
    nodes, hvs, tnodes = [], [], []
    S_hat = None

    for T in tensors:
        if T.label in ("C", "s"):
            p, q = T.bra[0][0], T.ket[0][0]
            uf.union(p, q)
            connectors.append((T.label, p, q))
        elif T.label in ("g", "f"):
            n = Node("hv", f"{T.label}{len(hvs)}", T)
            hvs.append(n); nodes.append(n)
        elif T.label == "t":
            n = Node("t", f"t{len(tnodes)}", T)
            tnodes.append(n); nodes.append(n)
        elif T.label == "S^":
            S_hat = T
        else:
            raise ValueError(f"unhandled tensor {T}")

    def add_term(node, side, idx):
        terminal_occ.setdefault(uf.find(idx[0]), []).append((node, side, idx))

    for n in nodes:
        for idx in n.tensor.bra:
            add_term(n, "bra", idx)
        for idx in n.tensor.ket:
            add_term(n, "ket", idx)

    # externals: one ext node per Ŝ index. Ŝ bra=occ (creation at top edge),
    # Ŝ ket=virt. Fermion flow: occ Ŝ->partner (down); virt partner->Ŝ (up).
    ext_nodes = {}
    if S_hat is not None:
        for idx in S_hat.bra:                    # occupied externals
            en = Node("ext", idx[0]); en.ext_idx = idx
            ext_nodes[idx[0]] = en
            add_term(en, "bra", idx)
        for idx in S_hat.ket:                    # virtual externals
            en = Node("ext", idx[0]); en.ext_idx = idx
            ext_nodes[idx[0]] = en
            add_term(en, "ket", idx)

    # count connectors per chain
    for lab, p, q in connectors:
        r = uf.find(p)
        c = conn_count.setdefault(r, [0, 0])
        c[0 if lab == "C" else 1] += 1

    lines = []
    for root, occs in terminal_occ.items():
        if len(occs) != 2:
            raise ValueError(f"chain {root} has {len(occs)} terminals (want 2): {occs}")
        (n1, s1, i1), (n2, s2, i2) = occs
        if s1 == s2:
            raise ValueError(f"chain {root}: both terminals are {s1} slots")
        if s1 == "ket":
            (n1, s1, i1), (n2, s2, i2) = (n2, s2, i2), (n1, s1, i1)
        nC, nS = conn_count.get(root, [0, 0])
        lines.append(FermionLine(n1, n2, i1, i2, nC, nS))

    # interaction pairing: g's sharing an aux index (DF), else lone f/g
    hv_pair = None
    if len(hvs) == 2:
        hv_pair = (hvs[0], hvs[1])
    elif len(hvs) > 2:
        raise ValueError(">2 interaction half-vertices not supported")
    return nodes + list(ext_nodes.values()), lines, hvs, hv_pair, ext_nodes


# ----------------------------------------------------------------------------
# labels
# ----------------------------------------------------------------------------

class Opts:
    """Rendering knobs.  show_csv=False => canonical CC Goldstone skeleton (rule B0).

    The three overlay layers can also be toggled independently:
      csv  C / CsC / s boxes on the fermion lines
      dom  PNO domain superscripts a^{ij}
      aux  the DF label K on the interaction line (+ the g_L/g_R half-vertex tags)
    """
    def __init__(self, show_csv=True, show_aux=None, show_dom=None):
        self.csv = bool(show_csv)
        self.aux = self.csv if show_aux is None else bool(show_aux)
        self.dom = self.csv if show_dom is None else bool(show_dom)


def _sub(name):
    """'a_3' -> 'a_3' latex; 'μ̃_1' -> '\\tilde{\\mu}_1'"""
    if "μ" in name:
        n = name.split("_")[-1]
        return rf"\tilde{{\mu}}_{{{n}}}"
    base, _, n = name.partition("_")
    return rf"{base}_{{{n}}}" if n else base


def idx_latex(idx, opts):
    name, dom = idx
    s = _sub(name)
    if dom and opts.dom:
        s += "^{\\," + "\\,".join(_sub(d) for d in dom) + "}"
    return f"${s}$"


def t_label(tensor, opts):
    occ = "\\,".join(_sub(n) for n, _ in tensor.ket)
    vir = "\\,".join(_sub(n) + ("^{" + "\\,".join(_sub(d) for d in dom) + "}"
                                if (dom and opts.dom) else "")
                     for n, dom in tensor.bra)
    return rf"$\hat{{t}}_{{\,{occ}}}^{{\,{vir}}}$"


# ----------------------------------------------------------------------------
# text metrics (DATA units) — measured once on a scratch canvas
# ----------------------------------------------------------------------------

_SCRATCH = None
_SIZE_CACHE = {}


def text_size(s, fontsize):
    """(width, height) of a rendered (math)text string, in data units."""
    key = (s, fontsize)
    hit = _SIZE_CACHE.get(key)
    if hit is not None:
        return hit
    global _SCRATCH
    if _SCRATCH is None:
        _SCRATCH = plt.figure(figsize=(2, 2), dpi=100)
        _SCRATCH.canvas.draw()
    art = _SCRATCH.text(0.5, 0.5, s, fontsize=fontsize)
    bb = art.get_window_extent(_SCRATCH.canvas.get_renderer())
    art.remove()
    res = (bb.width / _SCRATCH.dpi / UNIT, bb.height / _SCRATCH.dpi / UNIT)
    _SIZE_CACHE[key] = res
    return res


def box_of(x, y, w, h, ha="center", va="center", pad=0.06):
    if ha == "center":
        x0, x1 = x - w / 2, x + w / 2
    elif ha == "left":
        x0, x1 = x, x + w
    else:
        x0, x1 = x - w, x
    if va == "center":
        y0, y1 = y - h / 2, y + h / 2
    elif va == "bottom":
        y0, y1 = y, y + h
    else:
        y0, y1 = y - h, y
    return (x0 - pad, y0 - pad, x1 + pad, y1 + pad)


# ----------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------

def _perp(u):
    return np.array([u[1], -u[0]])


def bow_path(p0, p1, bow, n=41):
    """Quadratic-Bezier arc from p0 to p1 bulging by `bow` (perp, +=right of p0->p1).
    Two opposite bows between the same endpoints make a pointy ellipse (lens)."""
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    d = p1 - p0
    L = float(np.hypot(*d))
    if L < 1e-9 or abs(bow) < 1e-9:
        return np.array([p0, p1])
    ctrl = (p0 + p1) / 2 + _perp(d / L) * 2.0 * bow
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * ctrl + t ** 2 * p1


def path_at(pts, f):
    """point + unit tangent at arc-length fraction f of a polyline."""
    seg = np.diff(pts, axis=0)
    ln = np.hypot(seg[:, 0], seg[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(ln)])
    tot = cum[-1]
    if tot < 1e-9:
        return pts[0].copy(), np.array([0.0, 1.0])
    s = float(np.clip(f, 0.0, 1.0)) * tot
    k = int(np.clip(np.searchsorted(cum, s) - 1, 0, len(seg) - 1))
    u = (s - cum[k]) / max(ln[k], 1e-12)
    return pts[k] + seg[k] * u, seg[k] / max(ln[k], 1e-12)


def _overlap(a, b):
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if (dx > 0 and dy > 0) else 0.0


def _seg_hits_box(p, q, box):
    """Liang-Barsky segment/AABB test."""
    x0, y0, x1, y1 = box
    tmin, tmax = 0.0, 1.0
    dx, dy = q[0] - p[0], q[1] - p[1]
    for pp, qq in ((-dx, p[0] - x0), (dx, x1 - p[0]),
                   (-dy, p[1] - y0), (dy, y1 - p[1])):
        if abs(pp) < 1e-12:
            if qq < 0:
                return False
        else:
            r = qq / pp
            if pp < 0:
                if r > tmax:
                    return False
                tmin = max(tmin, r)
            else:
                if r < tmin:
                    return False
                tmax = min(tmax, r)
    return tmin <= tmax


class Placer:
    """Greedy collision-avoiding label placement against boxes and line segments."""

    def __init__(self):
        self.boxes = []                  # occupied rectangles
        self.segs = []                   # (p, q, owner)
        self.tags = []                   # parallel to boxes: (tag, is_label)

    def block(self, box, tag="obstacle", is_label=False, owner=None):
        self.boxes.append(box)
        self.tags.append((tag, is_label, id(owner) if owner is not None else None))
        return box

    def block_text(self, x, y, s, fs, ha="center", va="center", pad=0.06):
        w, h = text_size(s, fs)
        return self.block(box_of(x, y, w, h, ha, va, pad))

    def add_path(self, pts, owner, k=12):
        pts = np.asarray(pts, float)
        idx = np.unique(np.linspace(0, len(pts) - 1, min(k + 1, len(pts))).astype(int))
        for a, b in zip(idx[:-1], idx[1:]):
            self.segs.append((pts[a], pts[b], owner))

    def cost(self, box, owner=None):
        """Any overlap with an existing box is disqualifying-expensive; sitting on
        a foreign fermion line is merely ugly."""
        c = 0.0
        for b in self.boxes:
            ov = _overlap(box, b)
            if ov > 1e-9:
                c += 3.0 + 25.0 * ov
        for p, q, o in self.segs:
            if _seg_hits_box(p, q, box):
                c += 0.30 if o is owner else 0.85
        return c

    def place(self, cands, w, h, owner=None, pad=0.06, tag="label"):
        """cands: [((x, y), preference_penalty)] — lower penalty = more wanted."""
        best = None
        for x, y, pen in cands:
            box = box_of(x, y, w, h, "center", "center", pad)
            c = self.cost(box, owner) + pen
            if best is None or c < best[0]:
                best = (c, (x, y), box)
        self.block(best[2], tag, True, owner)
        return best[1]


# ----------------------------------------------------------------------------
# layout
# ----------------------------------------------------------------------------

def _bar_half(tensor):
    m = max(len(tensor.bra), len(tensor.ket), 1)
    return APEX_STEP * (m - 1) / 2 + BAR_HALF


def assign_apexes(tnodes, lines):
    """Map (t-node, line) -> apex point on the bar.  The legs of one excitation
    (particle+hole pair) MEET at a shared apex: a T_m vertex is m V-shapes."""
    slots = {}
    for n in tnodes:
        bra_names = [nm for nm, _ in n.tensor.bra]   # virtuals -> particle legs
        ket_names = [nm for nm, _ in n.tensor.ket]   # occupieds -> hole legs
        m = max(len(bra_names), len(ket_names), 1)
        bra_line = [None] * len(bra_names)
        ket_line = [None] * len(ket_names)
        for L in lines:
            if L.src is n:
                idx = L.src_idx
            elif L.dst is n:
                idx = L.dst_idx
            else:
                continue
            nm = idx[0]
            if nm in bra_names:
                bra_line[bra_names.index(nm)] = L
            elif nm in ket_names:
                ket_line[ket_names.index(nm)] = L

        def other_x(L):
            if L is None:
                return n.x
            o = L.dst if L.src is n else L.src
            return o.x

        pairs = [(bra_line[k] if k < len(bra_line) else None,
                  ket_line[k] if k < len(ket_line) else None) for k in range(m)]
        pairs.sort(key=lambda pr: float(np.mean([other_x(L) for L in pr
                                                 if L is not None] or [n.x])))
        if m == 1:
            apex_xs = [n.x]
        else:
            apex_xs = [n.x - APEX_STEP * (m - 1) / 2 + APEX_STEP * k for k in range(m)]
        for k, (bl, kl) in enumerate(pairs):
            for L in (bl, kl):
                if L is not None:
                    slots[(n, id(L))] = (apex_xs[k], n.y)
        for L in lines:                              # safety: stray leg -> centre
            if (L.src is n or L.dst is n) and (n, id(L)) not in slots:
                slots[(n, id(L))] = (n.x, n.y)
        n.apex_xs = apex_xs
    return slots


def _spread(pos, half, gap):
    """Push overlapping items apart symmetrically (centroid preserving)."""
    pos = list(pos)
    for _ in range(400):
        moved = False
        for k in range(len(pos) - 1):
            need = half[k] + half[k + 1] + gap
            d = pos[k + 1] - pos[k]
            if d < need - 1e-9:
                push = (need - d) / 2.0
                pos[k] -= push
                pos[k + 1] += push
                moved = True
        if not moved:
            break
    return pos


def _symmetrize(hvs):
    if len(hvs) != 2:
        return
    a, b = hvs
    if a.x > b.x:
        a.x, b.x = b.x, a.x
    mid = (a.x + b.x) / 2.0
    halfgap = max(HV_GAP, b.x - a.x) / 2.0
    a.x, b.x = mid - halfgap, mid + halfgap


def layout(nodes, lines, hvs, opts):
    """Place every node; returns the apex map."""
    tnodes = [n for n in nodes if n.kind == "t"]
    exts = [n for n in nodes if n.kind == "ext"]

    # which half-vertices does each amplitude touch? (left / straddling / right)
    touch = {n: set() for n in tnodes}
    for L in lines:
        for a, b in ((L.src, L.dst), (L.dst, L.src)):
            if a.kind == "t" and b.kind == "hv":
                touch[a].add(b.tag)

    def group(n):
        s = touch[n]
        if len(hvs) == 2:
            if s == {hvs[0].tag}:
                return 0
            if s == {hvs[0].tag, hvs[1].tag}:
                return 1
            if s == {hvs[1].tag}:
                return 2
        return 1
    tnodes.sort(key=group)

    # --- bottom row: spacing driven by bar width AND amplitude-label width
    half = [max(text_size(t_label(n.tensor, opts), FS_T)[0] / 2, _bar_half(n.tensor))
            for n in tnodes]
    x = 0.0
    for k, n in enumerate(tnodes):
        if k:
            x += half[k - 1] + half[k] + T_GAP
        n.x, n.y = x, Y_T
    if tnodes:                                        # centre the row on 0
        c = (tnodes[0].x + tnodes[-1].x) / 2
        for n in tnodes:
            n.x -= c

    def hv_partner_xs(hv, slots=None):
        out = []
        for L in lines:
            for a, b in ((L.src, L.dst), (L.dst, L.src)):
                if a is hv and b.kind == "t":
                    out.append(slots[(b, id(L))][0] if slots else b.x)
        return out

    partner_of = {}
    line_of = {}
    for en in exts:
        for L in lines:
            if L.src is en or L.dst is en:
                partner_of[en] = L.dst if L.src is en else L.src
                line_of[en] = L

    def place_hvs(slots=None):
        for j, hv in enumerate(hvs):
            px = hv_partner_xs(hv, slots)
            hv.x = float(np.mean(px)) if px else (HV_GAP * (j - 0.5))
            hv.y = Y_HV
        _symmetrize(hvs)

    def ext_pref(en, slots=None):
        p = partner_of.get(en)
        if p is None:
            return 0.0
        if p.kind == "t" and slots is not None:
            px = slots[(p, id(line_of[en]))][0]
        else:
            px = p.x
        return px + (-EXT_BIAS if kind_of(en.tag) == "occ" else EXT_BIAS)

    # provisional -> apexes -> refined half-vertices -> apexes -> externals
    place_hvs()
    for en in exts:
        en.x, en.y = ext_pref(en), Y_EXT
    slots = assign_apexes(tnodes, lines)
    place_hvs(slots)
    slots = assign_apexes(tnodes, lines)

    if exts:
        exts.sort(key=lambda e: (ext_pref(e, slots), e.tag))
        halfw = [text_size(idx_latex(e.ext_idx, opts), FS_IDX)[0] / 2 for e in exts]
        prefs = [ext_pref(e, slots) for e in exts]
        for _ in range(3):
            pos = _spread(prefs, halfw, EXT_GAP)
            # keep external stubs from spearing an interaction dot
            bumped = False
            for k, en in enumerate(exts):
                p = partner_of.get(en)
                if p is None or p.kind == "hv":
                    continue
                p0 = slots[(p, id(line_of[en]))] if p.kind == "t" else (p.x, p.y)
                if abs(Y_EXT - p0[1]) < 1e-9:
                    continue
                f = (Y_HV - p0[1]) / (Y_EXT - p0[1])
                xh = p0[0] + (pos[k] - p0[0]) * f
                for hv in hvs:
                    if abs(xh - hv.x) < 0.45:
                        d = 0.45 - abs(xh - hv.x)
                        s = 1.0 if pos[k] >= hv.x else -1.0
                        prefs[k] = pos[k] + s * (d / max(f, 0.2) + 0.15)
                        bumped = True
            if not bumped:
                break
        pos = _spread(prefs, halfw, EXT_GAP)
        for en, x in zip(exts, pos):
            en.x, en.y = x, Y_EXT

    # --- global centring: mirror-image terms come out as mirror images
    allx = [n.x for n in tnodes + exts + hvs]
    for n in tnodes:
        allx += list(getattr(n, "apex_xs", []))
    if allx:
        c = (min(allx) + max(allx)) / 2.0
        for n in tnodes + exts + hvs:
            n.x -= c
        slots = assign_apexes(tnodes, lines)
    return slots


# ----------------------------------------------------------------------------
# drawing
# ----------------------------------------------------------------------------

def _sides(bow, p, q):
    """Perpendicular sides to try first: the bulge side of a bowed line, else the
    side pointing away from the diagram axis (keeps mirror terms mirror-imaged)."""
    if abs(bow) > 1e-9:
        return (1, -1) if bow > 0 else (-1, 1)
    want = 1.0 if p[0] >= 0 else -1.0
    first = 1 if q[0] * want >= 0 else -1
    return (first, -first)


def _line_label_pairs(L, opts):
    """(text_near_src, text_near_dst) — a bridge shows a different index at each
    end, an ordinary line shows one label.  PAO ends borrow the other end's name."""
    ks, kd = kind_of(L.src_idx[0]), kind_of(L.dst_idx[0])
    a = L.src_idx if ks != "pao" else (L.dst_idx if kd != "pao" else L.src_idx)
    b = L.dst_idx if kd != "pao" else (L.src_idx if ks != "pao" else L.dst_idx)
    ta, tb = idx_latex(a, opts), idx_latex(b, opts)
    return (ta, tb) if ta != tb else (ta, None)


def draw(term, title=None, out="diagram.png", verbose=True,
         show_csv=True, show_aux=None, show_dom=None, legend=True, dpi=170):
    opts = Opts(show_csv, show_aux, show_dom)
    coeff, tensors = parse_term(term)
    nodes, lines, hvs, hv_pair, ext_nodes = build_graph(tensors)
    slots = layout(nodes, lines, hvs, opts)
    tnodes = [n for n in nodes if n.kind == "t"]

    if verbose:
        print(f"coeff={coeff}; {len(hvs)} half-vertices; {len(tnodes)} amplitudes; "
              f"{len(ext_nodes)} external lines")
        for L in lines:
            print(f"  {'hole' if L.is_hole else 'part'} "
                  f"{L.src}({L.src_idx[0]}) -> {L.dst}({L.dst_idx[0]})"
                  f"  ann={L.annotation}")

    def endpoint(node, L):
        return np.array(slots[(node, id(L))] if node.kind == "t" else (node.x, node.y),
                        float)

    # ---- fermion-line geometry (bowing lines that share both endpoints) -----
    geo = {}
    groups = {}
    for L in lines:
        p, q = endpoint(L.src, L), endpoint(L.dst, L)
        key = tuple(sorted([tuple(np.round(p, 4)), tuple(np.round(q, 4))]))
        groups.setdefault(key, []).append(L)
    def _bow_line(L, bw):
        """(re)build L's polyline as an arc bulging `bw` off the straight chord."""
        p, q = endpoint(L.src, L), endpoint(L.dst, L)
        lo, hi = (p, q) if (p[1], p[0]) <= (q[1], q[0]) else (q, p)
        pts = bow_path(lo, hi, bw)
        if np.hypot(*(pts[0] - p)) > np.hypot(*(pts[-1] - p)):
            pts = pts[::-1]
        geo[id(L)] = dict(pts=pts, bow=bw, color=HOLE if L.is_hole else PART)

    # (a) lines sharing BOTH endpoints -> pointy ellipse (Goldstone bubble)
    lens = set()
    for key, grp in groups.items():
        n = len(grp)
        if n == 1:
            _bow_line(grp[0], 0.0)
            continue
        grp.sort(key=lambda L: (not L.is_hole, L.src_idx[0]))
        chord = np.hypot(*(np.array(key[1]) - np.array(key[0])))
        span = float(np.clip(0.10 * chord, 0.14, 0.30))   # slim lens
        for k, L in enumerate(grp):
            _bow_line(L, -span + 2 * span * k / (n - 1))
            lens.add(id(L))

    # (b) lines that merely RUN TOGETHER (sharing one endpoint, or an external
    #     stub grazing a half-vertex) are bowed apart too — just enough to read
    #     as two lines.  Iterated, because bowing can uncover a new near-miss.
    def _samples(pts, n=24):
        d = np.cumsum(np.r_[0.0, np.hypot(*np.diff(pts, axis=0).T)])
        if d[-1] < 1e-9:
            return np.repeat(pts[:1], n, axis=0)
        u = np.linspace(0, d[-1], n)
        return np.column_stack([np.interp(u, d, pts[:, 0]),
                                np.interp(u, d, pts[:, 1])])

    def _runs_together(A, B, smp, hair=0.11, frac=0.30):
        a, b = smp[id(A)], smp[id(B)]
        d = np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1])
        return max((d.min(axis=1) < hair).mean(), (d.min(axis=0) < hair).mean()) > frac

    def _key(L):
        return (not L.is_hole, L.src_idx[0], L.dst_idx[0])

    # A bundle is fanned as a whole; the bundle graph only ever GROWS (so the
    # fan never oscillates), and the fan widens only once membership is stable.
    uf2, sep = UF(), 0.48
    for _ in range(6):
        smp = {id(L): _samples(geo[id(L)]["pts"]) for L in lines}
        hits = [(A, B) for i, A in enumerate(lines) for B in lines[i + 1:]
                if _runs_together(A, B, smp)]
        if not hits:
            break
        grew = False
        for A, B in hits:
            if uf2.find(id(A)) != uf2.find(id(B)):
                uf2.union(id(A), id(B))
                grew = True
        if not grew:
            sep = min(sep * 1.3, 0.95)
        comps = {}
        for L in lines:
            comps.setdefault(uf2.find(id(L)), []).append(L)
        for grp in comps.values():
            if len(grp) < 2:
                continue
            free = [L for L in sorted(grp, key=_key) if id(L) not in lens]
            if len(free) == 1:                 # push it off whatever it runs with
                other = [L for L in grp if L is not free[0]]
                b0 = np.sign(geo[id(other[0])]["bow"] or 1.0) if other else 1.0
                _bow_line(free[0], -b0 * sep)
            elif free:
                bows = [(k - (len(free) - 1) / 2) * sep for k in range(len(free))]
                cap = max(abs(b) for b in bows)          # keep arcs from wandering
                if cap > 0.55:
                    bows = [b * 0.55 / cap for b in bows]
                for L, bw in zip(free, bows):
                    _bow_line(L, bw)

    # ---- obstacles ---------------------------------------------------------
    P = Placer()
    texts = []      # (x, y, s, color, fs, ha, va, style)

    for n in tnodes:
        x0, x1 = n.apex_xs[0] - BAR_HALF, n.apex_xs[-1] + BAR_HALF
        P.block((x0 - 0.10, n.y - 0.14, x1 + 0.10, n.y + 0.14), "bar")
        s = t_label(n.tensor, opts)
        w, h = text_size(s, FS_T)
        P.block(box_of(n.x, n.y - 0.34, w, h, "center", "top"), f"t:{s}", True)
        texts.append((n.x, n.y - 0.34, s, "black", FS_T, "center", "top", "t"))

    for hv in hvs:
        P.block((hv.x - 0.20, hv.y - 0.20, hv.x + 0.20, hv.y + 0.20), "dot")

    for en in [n for n in nodes if n.kind == "ext"]:
        s = idx_latex(en.ext_idx, opts)
        w, h = text_size(s, FS_IDX)
        col = HOLE if kind_of(en.tag) == "occ" else PART
        P.block(box_of(en.x, en.y + 0.22, w, h, "center", "bottom"), f"ext:{s}", True)
        texts.append((en.x, en.y + 0.22, s, col, FS_IDX, "center", "bottom", "ext"))

    for L in lines:
        P.add_path(geo[id(L)]["pts"], owner=L)

    # ---- interaction -------------------------------------------------------
    dash_seg = None
    xmark = None
    if len(hvs) == 2:
        a, b = hvs
        dash_seg = ((a.x, a.y), (b.x, b.y))
        P.segs.append((np.array(dash_seg[0]), np.array(dash_seg[1]), "dash"))
        if opts.aux:
            aux = a.tensor.aux[0][0] if a.tensor.aux else None
            if aux:
                w, h = text_size("$K$", FS_V)
                # K belongs ON the dashed line; step sideways before stepping off it
                cands = [((a.x + b.x) / 2 + d, a.y + dy, 0.30 * abs(d) + 1.30 * abs(dy))
                         for d in (0, -0.55, 0.55, -1.1, 1.1, -1.65, 1.65)
                         for dy in (0, 0.42, -0.42)]
                x, y = P.place(cands, w, h, owner="dash", tag="K")
                texts.append((x, y, "$K$", INT, FS_V, "center", "center", "anno_int"))
            for hv, lbl, side in ((a, "$g_L$", -1.0), (b, "$g_R$", 1.0)):
                w, h = text_size(lbl, FS_V)
                r = 0.30 + w / 2
                cands = []
                for k, (dx, dy) in enumerate(((side * r, 0.44), (side * r, -0.44),
                                              (-side * r, 0.44), (-side * r, -0.44))):
                    for grow in (0.0, 0.45, 0.95):
                        cands.append((hv.x + dx + np.sign(dx) * grow, hv.y + dy,
                                      0.18 * k + 0.55 * grow))
                x, y = P.place(cands, w, h, owner="dash", tag=lbl)
                texts.append((x, y, lbl, INT, FS_V, "center", "center", "hv"))
    elif len(hvs) == 1:
        hv = hvs[0]
        best, bestc = (1.0, STUB), None
        for s in (1.0, -1.0):
            for ln in (STUB, STUB + 0.45, STUB + 0.9):
                tip = hv.x + s * ln
                box = (min(hv.x, tip) - 0.1, hv.y - 0.28,
                       max(hv.x, tip) + 0.24, hv.y + 0.28)
                # the x marker itself must sit clear of everything
                c = P.cost(box) + 3.0 * P.cost((tip - 0.26, hv.y - 0.26,
                                                tip + 0.26, hv.y + 0.26))
                if bestc is None or c < bestc - 1e-9:
                    best, bestc = (s, ln), c
        s, ln = best
        dash_seg = ((hv.x, hv.y), (hv.x + s * ln, hv.y))
        xmark = (hv.x + s * ln, hv.y)
        P.block((min(hv.x, xmark[0]) - 0.05, hv.y - 0.22,
                 max(hv.x, xmark[0]) + 0.22, hv.y + 0.22))

    # ---- CSV annotations (C / CsC / s) -------------------------------------
    for L in lines:
        ann = L.annotation if opts.csv else None
        if not ann:
            continue
        g = geo[id(L)]
        pts = g["pts"]
        w, h = text_size(f"${ann}$", FS_ANN)
        # "C" belongs next to the g/f end; bridges sit mid-line
        # C / CsC / s ride ON the line (a bead threaded on the fermion line);
        # they slide ALONG it to dodge obstacles, and only step aside as a
        # last resort.  "C" stays near its g/f end (rule B4).
        if ann == "C" and (L.src.kind == "hv" or L.dst.kind == "hv"):
            near_src = L.src.kind == "hv"
            fr = [0.18, 0.26, 0.34, 0.42, 0.12, 0.50, 0.58, 0.66]
            if not near_src:
                fr = [1.0 - f for f in fr]
        else:
            fr = [0.50, 0.42, 0.58, 0.34, 0.66, 0.26, 0.74, 0.18, 0.82]
        cands = [(*path_at(pts, f)[0], 0.12 * j) for j, f in enumerate(fr)]
        for j, f in enumerate(fr[:4]):                   # last-resort step-aside
            p, t = path_at(pts, f)
            q = _perp(t)
            off = abs(q[0]) * w / 2 + abs(q[1]) * h / 2 + 0.16
            for k, s in enumerate(_sides(g["bow"], p, q)):
                cands.append((p[0] + s * q[0] * off, p[1] + s * q[1] * off,
                              2.2 + 0.12 * j + 0.14 * k))
        x, y = P.place(cands, w, h, owner=L, tag=f"ann:{ann}")
        texts.append((x, y, f"${ann}$", g["color"], FS_ANN, "center", "center", "anno"))

    # ---- index labels ------------------------------------------------------
    for L in lines:
        g = geo[id(L)]
        pts = g["pts"]
        ta, tb = _line_label_pairs(L, opts)
        jobs = []
        src_ext = L.src.kind == "ext"
        dst_ext = L.dst.kind == "ext"
        if tb is None:                                  # one name for the line
            if not (src_ext or dst_ext):                # externals label at the top
                jobs.append((ta, [0.38, 0.52, 0.26, 0.64, 0.14, 0.76, 0.88]))
        else:                                           # bridge: name at each end
            if not src_ext:
                jobs.append((ta, [0.24, 0.34, 0.16, 0.44, 0.08, 0.52]))
            if not dst_ext:
                jobs.append((tb, [0.76, 0.66, 0.84, 0.56, 0.92, 0.48]))
        for s, fr in jobs:
            w, h = text_size(s, FS_IDX)
            cands = []
            for extra in (0.17, 0.52, 0.95, 1.55):      # progressively roomier
                for j, f in enumerate(fr):
                    p, t = path_at(pts, f)
                    q = _perp(t)
                    off = abs(q[0]) * w / 2 + abs(q[1]) * h / 2 + extra
                    for k, sg in enumerate(_sides(g["bow"], p, q)):
                        cands.append((p[0] + sg * q[0] * off, p[1] + sg * q[1] * off,
                                      1.40 * (extra - 0.17) + 0.09 * j + 0.13 * k))
            x, y = P.place(cands, w, h, owner=L, tag=f"idx:{s}")
            texts.append((x, y, s, g["color"], FS_IDX, "center", "center", "idx"))

    # ---- audit record (see audit_diagrams.py) ------------------------------
    global LAST_AUDIT
    LAST_AUDIT = dict(boxes=list(P.boxes), tags=list(P.tags),
                      segs=[(tuple(p), tuple(q), o) for p, q, o in P.segs])

    # ---- extents -----------------------------------------------------------
    xs, ys = [], []
    for L in lines:
        pts = geo[id(L)]["pts"]
        xs += [pts[:, 0].min(), pts[:, 0].max()]
        ys += [pts[:, 1].min(), pts[:, 1].max()]
    for b in P.boxes:
        xs += [b[0], b[2]]; ys += [b[1], b[3]]
    if dash_seg:
        xs += [dash_seg[0][0], dash_seg[1][0]]; ys += [dash_seg[0][1], dash_seg[1][1]]
    xs = xs or [-1.0, 1.0]; ys = ys or [0.0, 1.0]
    m = 0.30
    xr = max(abs(min(xs) - m), abs(max(xs) + m))
    xmin, xmax = -xr, xr                    # symmetric frame about the diagram axis
    ymin, ymax = min(ys) - m, max(ys) + m
    if ext_nodes:
        ymax = max(ymax, Y_EXT + 0.60)

    # ---- figure ------------------------------------------------------------
    LM = RM = 0.14
    TM = 0.52 if title else 0.14
    BM = 0.62 if legend else 0.14
    figw = (xmax - xmin) * UNIT + LM + RM
    figh = (ymax - ymin) * UNIT + TM + BM

    need = 0.6
    if title:
        need = max(need, text_size(title, FS_TITLE)[0] * UNIT + 0.4)
    if legend:
        need = max(need, 6.1)
    if figw < need:                          # widen symmetrically
        extra = (need - figw) / UNIT
        xmin -= extra / 2; xmax += extra / 2
        figw = (xmax - xmin) * UNIT + LM + RM

    fig = plt.figure(figsize=(figw, figh), dpi=dpi)
    ax = fig.add_axes([LM / figw, BM / figh,
                       (xmax - xmin) * UNIT / figw, (ymax - ymin) * UNIT / figh])
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.axis("off")

    if ext_nodes:                            # projection boundary (no caption)
        ax.axhline(Y_EXT, color="black", lw=0.6, ls=(0, (1, 4)), alpha=0.35, zorder=1)

    # amplitude bars
    for n in tnodes:
        ax.plot([n.apex_xs[0] - BAR_HALF, n.apex_xs[-1] + BAR_HALF], [n.y, n.y],
                color="black", lw=4, solid_capstyle="butt", zorder=5)

    # interaction
    if dash_seg:
        ax.plot([dash_seg[0][0], dash_seg[1][0]], [dash_seg[0][1], dash_seg[1][1]],
                color=INT, lw=1.6, ls=(0, (6, 4)), zorder=2)
    for hv in hvs:
        ax.plot([hv.x], [hv.y], "o", color=INT, ms=11, zorder=6)
    if xmark is not None:                    # one-body vertex:  ---x
        ax.plot([xmark[0]], [xmark[1]], marker="x", color=INT, ms=13, mew=2.6,
                zorder=6, ls="none")

    # fermion lines + arrows
    for L in lines:
        g = geo[id(L)]
        pts = g["pts"]
        ax.plot(pts[:, 0], pts[:, 1], color=g["color"], lw=1.8, zorder=3,
                solid_capstyle="round")
        p, t = path_at(pts, 0.5)
        ax.annotate("", xy=p + t * 0.28, xytext=p - t * 0.28,
                    arrowprops=dict(arrowstyle="-|>", color=g["color"], lw=0,
                                    mutation_scale=22), zorder=6)

    # labels
    for x, y, s, col, fs, ha, va, style in texts:
        if style in ("anno", "anno_int"):
            ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col, zorder=9,
                    bbox=dict(boxstyle="round,pad=0.16", fc="#fff7d6", ec=col,
                              alpha=0.96))
        elif style == "t":
            ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col,
                    fontweight="bold", zorder=9)
        elif style == "hv":
            ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col,
                    fontweight="bold", zorder=9)
        else:
            ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col, zorder=8,
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none",
                              alpha=0.80))

    if legend:
        if len(hvs) == 2:
            ilab = "two-body vertex" + (" ($K$)" if opts.aux else "")
            ihandle = Line2D([0], [0], color=INT, lw=1.6, ls="--")
        else:
            ilab = "one-body vertex ($-\\!-\\!\\times$)"
            ihandle = Line2D([0], [0], color=INT, lw=1.6, ls="--")
        leg = [Line2D([0], [0], color=PART, lw=2, label="particle line (up)"),
               Line2D([0], [0], color=HOLE, lw=2, label="hole line (down)"),
               (ihandle, ilab)]
        fig.legend(handles=[leg[0], leg[1], leg[2][0]],
                   labels=[leg[0].get_label(), leg[1].get_label(), leg[2][1]],
                   loc="lower center", ncol=3, frameon=False, fontsize=FS_LEG,
                   bbox_to_anchor=(0.5, 0.004))
    if title:
        ax.set_title(f"{title}   (coefficient {coeff})", fontsize=FS_TITLE, pad=8)

    out = os.path.expanduser(out)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    if verbose:
        print("saved:", out)
    return out


# ----------------------------------------------------------------------------
# demo: R2 terms from pao-ccsd-df.md
# ----------------------------------------------------------------------------

TERM43 = ("Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:N-C-S * g{i_3;μ̃_1;Κ_1}:N-C-S * "
          "C{μ̃_1;a_3<i_1>}:N-C-S * g{i_4;μ̃_2;Κ_1}:N-C-S * C{μ̃_2;a_4<i_2>}:N-C-S * "
          "s{μ̃_3;μ̃_4}:N-C-S * C{a_1<i_1,i_2>;μ̃_3}:N-C-S * C{μ̃_4;a_5<i_3>}:N-C-S * "
          "s{μ̃_5;μ̃_6}:N-C-S * C{a_2<i_1,i_2>;μ̃_5}:N-C-S * C{μ̃_6;a_6<i_4>}:N-C-S * "
          "t{a_3<i_1>;i_1}:N-C-S * t{a_4<i_2>;i_2}:N-C-S * t{a_5<i_3>;i_3}:N-C-S * "
          "t{a_6<i_4>;i_4}:N-C-S")

# R2 term 30: (oo|oo) interaction — both half-vertices carry two HOLE legs;
# both external particles reach the amplitudes only through CsC bridges.
TERM30 = ("((Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:N-C-S * ((s{μ̃_1;μ̃_2}:N-C-S * "
          "(C{μ̃_2;a_3<i_3>}:N-C-S * t{a_3<i_3>;i_3}:N-C-S)) * "
          "C{a_1<i_1,i_2>;μ̃_1}:N-C-S)) * (g{i_3;i_1;Κ_1}:N-C-S * "
          "g{i_4;i_2;Κ_1}:N-C-S)) * ((s{μ̃_3;μ̃_4}:N-C-S * (C{μ̃_4;a_4<i_4>}:N-C-S * "
          "t{a_4<i_4>;i_4}:N-C-S)) * C{a_2<i_1,i_2>;μ̃_3}:N-C-S)")

# R2 term 49: (vv|vv) interaction + one T2 — each half-vertex has particle-in
# (from the T2, via C) and particle-out (external, via C); the T2 also carries
# both external hole lines directly.
TERM49 = ("(Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:N-C-S * (((g{μ̃_1;μ̃_2;Κ_1}:N-C-S * "
          "C{a_1<i_1,i_2>;μ̃_1}:N-C-S) * C{μ̃_2;a_3<i_1,i_2>}:N-C-S) * "
          "((g{μ̃_3;μ̃_4;Κ_1}:N-C-S * C{a_2<i_1,i_2>;μ̃_3}:N-C-S) * "
          "C{μ̃_4;a_4<i_1,i_2>}:N-C-S))) * t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}:N-C-S")


def _cli():
    ap = argparse.ArgumentParser(
        description="render CSV-CC Goldstone diagram(s) from SeQuant term string(s)",
        epilog="the term may be given inline, with --file, or on stdin; markdown "
               "list numbering and backticks are stripped automatically")
    ap.add_argument("term", nargs="?", help="SeQuant term string ('-' = read stdin)")
    ap.add_argument("out", nargs="?", default="diagram.png",
                    help="output PNG (or output DIRECTORY when --file holds many terms)")
    ap.add_argument("--file", metavar="PATH",
                    help="read the term(s) from a file ('-' = stdin); every line "
                         "holding a tensor expression is rendered")
    ap.add_argument("--title", default=None)
    ap.add_argument("--no-csv", action="store_true",
                    help="canonical CC skeleton: drop C/CsC/s, domains and K")
    ap.add_argument("--no-dom", action="store_true", help="drop domain superscripts only")
    ap.add_argument("--no-aux", action="store_true", help="drop the K / g_L,g_R tags only")
    ap.add_argument("--no-legend", action="store_true")
    args = ap.parse_args()
    kw = dict(show_csv=not args.no_csv,
              show_dom=False if args.no_dom else None,
              show_aux=False if args.no_aux else None,
              legend=not args.no_legend)

    src = args.file or ("-" if args.term == "-" else None)
    if src is not None:                                  # file / stdin: 1..N terms
        # with --file the single positional is the OUTPUT ("--file in.txt out.png")
        out = args.out if (args.out != "diagram.png" or not args.file) else \
            (args.term or args.out)
        text = sys.stdin.read() if src == "-" else open(src, encoding="utf-8").read()
        terms = [ln for ln in text.splitlines() if "{" in ln and "}" in ln]
        if not terms:
            ap.error(f"no tensor expression found in {src}")
        if len(terms) == 1:
            draw(terms[0], title=args.title, out=out, **kw)
        else:
            outdir = out if (out.endswith("/") or "." not in os.path.basename(out)) \
                else (os.path.dirname(out) or ".")
            os.makedirs(outdir, exist_ok=True)
            print(f"{len(terms)} terms -> {outdir}/")
            for k, t in enumerate(terms, 1):
                draw(t, title=(args.title or "term") + f" {k}", verbose=False,
                     out=os.path.join(outdir, f"term{k:03d}.png"), **kw)
    elif args.term:
        draw(args.term, title=args.title, out=args.out, **kw)
    else:
        print("no term given — drawing the three built-in demo terms here:")
        for nm, tm, ti in (("term43", TERM43, r"CSV-CCSD $R_2$ term 43"),
                           ("term30", TERM30, r"CSV-CCSD $R_2$ term 30 (oo|oo)"),
                           ("term49", TERM49, r"CSV-CCSD $R_2$ term 49 (vv|vv)")):
            draw(tm, title=ti, out=f"csv_ccsd_{nm}.png", verbose=False, **kw)


if __name__ == "__main__":
    _cli()
