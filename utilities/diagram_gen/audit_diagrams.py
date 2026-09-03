#!/usr/bin/env python3
"""
audit_diagrams.py — machine check that every rendered CSV-CC diagram is clean.

For each term of pao-ccsd-df.md it re-runs the layout and reports
  1. label/label overlaps        (two text boxes sharing area)
  2. label/vertex-or-bar overlaps
  3. labels sitting on top of a fermion line that is not their own
  4. fermion lines drawn on top of each other (identical endpoints not bowed
     apart, or two segments nearly collinear over a long stretch)

Exit status 0 iff nothing is reported.  Run for both variants:
    python3 audit_diagrams.py            # CSV overlay
    python3 audit_diagrams.py --no-csv   # canonical skeleton
"""

import argparse
import os
import sys
import tempfile

import numpy as np

import csv_diagram as cd
from render_all import SRC, extract_terms, filter_sections, section_order


_TMP = os.path.join(tempfile.gettempdir(), "_csvdiag_audit.png")


def overlaps(a, b):
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if (dx > 0 and dy > 0) else 0.0


def audit_term(expr, show_csv):
    cd.draw(expr, title=None, out=_TMP, verbose=False,
            show_csv=show_csv, legend=False)
    A = cd.LAST_AUDIT
    boxes, tags, segs = A["boxes"], A["tags"], A["segs"]
    issues = []

    lab = [(b, t, ow) for b, (t, is_lab, ow) in zip(boxes, tags) if is_lab]
    obs = [(b, t) for b, (t, is_lab, ow) in zip(boxes, tags) if not is_lab]

    for i in range(len(lab)):
        for j in range(i + 1, len(lab)):
            ov = overlaps(lab[i][0], lab[j][0])
            if ov > 1e-4:
                issues.append(f"label/label {lab[i][1]} ~ {lab[j][1]}  ({ov:.3f})")
    for b, t, _ in lab:
        for ob, ot in obs:
            ov = overlaps(b, ob)
            if ov > 1e-4:
                issues.append(f"label/{ot} {t}  ({ov:.3f})")

    # a label box sitting on a fermion line OTHER than the one it belongs to
    for b, t, own in lab:
        foreign = set()
        for p, q, o in segs:
            if o in ("dash", None) or id(o) == own:
                continue
            if cd._seg_hits_box(np.array(p), np.array(q), b):
                foreign.add(id(o))
        if foreign:
            issues.append(f"label sits on {len(foreign)} foreign line(s): {t}")

    # duplicated / overlapping fermion strokes: sample each line's polyline and
    # measure how much of it runs within HAIR of another line
    paths = {}
    for p, q, o in segs:
        if o in ("dash", None):
            continue
        paths.setdefault(id(o), []).append((np.array(p), np.array(q)))
    sampled = {}
    for k, ss in paths.items():
        pts = [ss[0][0]] + [b for _, b in ss]
        pts = np.array(pts)
        d = np.cumsum(np.r_[0, np.hypot(*np.diff(pts, axis=0).T)])
        if d[-1] < 1e-9:
            continue
        u = np.linspace(0, d[-1], 60)
        sampled[k] = np.column_stack([np.interp(u, d, pts[:, 0]),
                                      np.interp(u, d, pts[:, 1])])
    HAIR = 0.07
    keys = list(sampled)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            A, B = sampled[keys[i]], sampled[keys[j]]
            dist = np.hypot(A[:, None, 0] - B[None, :, 0], A[:, None, 1] - B[None, :, 1])
            frac = float((dist.min(axis=1) < HAIR).mean())
            if frac > 0.35:
                issues.append(f"two fermion lines coincide over {frac:.0%} of a line")
    return issues


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=SRC, help="SeQuant equation markdown file")
    ap.add_argument("--only", default=None, help="e.g. R3 or R4:1,2")
    ap.add_argument("--no-csv", action="store_true")
    args = ap.parse_args()
    show_csv = not args.no_csv

    src = os.path.abspath(os.path.expanduser(args.src))
    if not os.path.exists(src):
        ap.error(f"equation file not found: {src} (pass one with --src <file.md>)")
    sections = extract_terms(src)
    if args.only:
        sections = filter_sections(sections, args.only)
    bad = n = 0
    for tag in section_order(sections):
        for num, expr in sections.get(tag, []):
            n += 1
            try:
                issues = audit_term(expr, show_csv)
            except Exception as exc:  # noqa: BLE001
                print(f"{tag} term {num}: RENDER FAILED {type(exc).__name__}: {exc}")
                bad += 1
                continue
            if issues:
                bad += 1
                print(f"{tag} term {num}:")
                for it in issues:
                    print(f"    {it}")
    variant = "CSV overlay" if show_csv else "canonical"
    print(f"\n[{variant}] {os.path.basename(args.src)}: "
          f"{bad} of {n} term(s) with issues")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
