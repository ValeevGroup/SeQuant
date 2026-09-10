#!/usr/bin/env python3
"""
render_all.py — batch-render every term of a SeQuant equation file (sections
E / R1 / R2 / R3 / R4 of pao-ccsd-df.md, pao-ccsdt-df.md, pao-ccsdtq-df.md, ...)
to a CSV-CC Goldstone diagram PNG using the csv_diagram renderer.

Two variants (rule B0, skeleton invariance):
  --variant csv        (default)  full CSV/DF overlay: domains, C/CsC boxes, K
  --variant canonical             plain CC Goldstone skeleton (no C/s, no K)
  --variant both                  render both sets

  python3 render_all.py                                   # CCSD, CSV overlay
  python3 render_all.py --variant both
  python3 render_all.py --src ../pao-ccsdt-df.md --only R3
  python3 render_all.py --src ../pao-ccsdtq-df.md --only R4:1,2,3 --outdir /tmp/q

Outputs go to diagrams/<outdir>/<section>_term<NN>.png (zero-padded); the default
outdir is derived from the source file name and the variant.  A manifest INDEX.md
is written summarizing successes/failures per section.
"""

import argparse
import os
import re

from csv_diagram import draw  # parse_term is used internally by draw

HERE = os.path.dirname(os.path.abspath(__file__))

# Default equation file: the CCSD one, wherever this script happens to live.
SRC_CANDIDATES = [
    os.path.join(HERE, "pao-ccsd-df.md"),                    # alongside the script
    os.path.join(os.path.dirname(HERE), "pao-ccsd-df.md"),   # one directory up
]
SRC = next((p for p in SRC_CANDIDATES if os.path.exists(p)), SRC_CANDIDATES[0])

# "## E — ...", "## R1 — ...", ... "## R4 — ..." (CCSD / CCSDT / CCSDTQ files)
SECTION_RE = re.compile(r"^##\s+(E|R\d+)\b")
TERM_RE = re.compile(r"^\s*(\d+)\.\s+`(.+?)`\s*$")


def extract_terms(path):
    """Return {section_tag: [(num, expr), ...]} in file order."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()

    sections = {}
    current = None
    for line in lines:
        stripped = line.rstrip("\n")
        if stripped.startswith("## "):
            m = SECTION_RE.match(stripped)
            current = m.group(1) if m else None      # Notation etc. -> no section
            if current:
                sections.setdefault(current, [])
            continue
        if current is None:
            continue
        m = TERM_RE.match(stripped)
        if m:
            sections[current].append((int(m.group(1)), m.group(2)))
    return sections


def section_order(sections):
    """E first, then R1, R2, R3, ... in numeric order."""
    return ([t for t in ("E",) if t in sections] +
            sorted((t for t in sections if t != "E"), key=lambda t: int(t[1:])))


def render_set(sections, outdir, show_csv, kind="CSV-CCSD"):
    os.makedirs(outdir, exist_ok=True)
    results = {}  # tag -> list of (num, filename_or_None, error_or_None)
    for tag in section_order(sections):
        terms = sections.get(tag, [])
        results[tag] = []
        for num, expr in terms:
            fname = f"{tag}_term{num:02d}.png"
            out = os.path.join(outdir, fname)
            title = f"{kind} {tag} term {num}"
            try:
                draw(expr, title=title, out=out, verbose=False, show_csv=show_csv)
                results[tag].append((num, fname, None))
            except Exception as exc:  # noqa: BLE001 - want to keep going
                msg = f"{type(exc).__name__}: {exc}"
                results[tag].append((num, None, msg))
                print(f"FAILED {tag} term {num}: {msg}")
    write_index(results, outdir, show_csv, kind)
    return results


def write_index(results, outdir, show_csv, kind="CSV-CCSD"):
    head = f"# {kind} diagram batch — INDEX"
    if not show_csv:
        head += "  (canonical skeleton)"
    lines = [head, ""]
    grand_ok = grand_fail = 0
    for tag in section_order(results):
        rows = results.get(tag, [])
        ok = [r for r in rows if r[2] is None]
        fail = [r for r in rows if r[2] is not None]
        grand_ok += len(ok)
        grand_fail += len(fail)
        lines.append(f"## {tag}  ({len(ok)} rendered / {len(fail)} failed of {len(rows)})")
        lines.append("")
        for num, fname, err in rows:
            if err is None:
                lines.append(f"- term {num}: `{fname}`")
            else:
                lines.append(f"- term {num}: FAILED: {err}")
        lines.append("")
    lines.insert(1, "")
    lines.insert(2, f"**Totals:** {grand_ok} rendered / {grand_fail} failed "
                    f"(of {grand_ok + grand_fail} terms).")
    with open(os.path.join(outdir, "INDEX.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{outdir}: {grand_ok} rendered / {grand_fail} failed")


def filter_sections(sections, only):
    """only = 'R2' or 'R2:16' or 'R2:16,49'"""
    tag, _, nums = only.partition(":")
    keep = {int(n) for n in nums.split(",") if n} if nums else None
    return {t: [(n, e) for n, e in v if keep is None or n in keep]
            for t, v in sections.items() if t == tag}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--src", default=SRC, help="SeQuant equation markdown file")
    ap.add_argument("--variant", choices=("csv", "canonical", "both"), default="csv")
    ap.add_argument("--outdir", default=None,
                    help="output directory (single-variant runs only)")
    ap.add_argument("--only", default=None,
                    help="restrict to e.g. R2 or R2:16 or R2:16,49")
    ap.add_argument("--kind", default=None,
                    help="title prefix (default: derived from --src)")
    args = ap.parse_args()

    src = os.path.abspath(os.path.expanduser(args.src))
    if not os.path.exists(src):
        ap.error("equation file not found: %s\nlooked for a default in:\n  %s\n"
                 "pass one explicitly with --src <file.md>"
                 % (src, "\n  ".join(SRC_CANDIDATES)))
    sections = extract_terms(src)
    if args.only:
        sections = filter_sections(sections, args.only)

    stem = os.path.splitext(os.path.basename(src))[0]        # e.g. pao-ccsdt-df
    theory = stem.replace("pao-", "").replace("-df", "").upper()  # CCSDT
    kind = args.kind or f"CSV-{theory}"
    base = "ccsd_all" if stem == "pao-ccsd-df" else f"{theory.lower()}_all"

    out = {}
    if args.variant in ("csv", "both"):
        d = (args.outdir if (args.outdir and args.variant == "csv")
             else os.path.join(HERE, base))
        out["csv"] = render_set(sections, d, True, kind)
    if args.variant in ("canonical", "both"):
        d = (args.outdir if (args.outdir and args.variant == "canonical")
             else os.path.join(HERE, base + "_canonical"))
        out["canonical"] = render_set(sections, d, False, theory)
    return out


if __name__ == "__main__":
    main()
