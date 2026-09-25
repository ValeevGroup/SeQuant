# diagram_gen — Goldstone diagrams from SeQuant term strings

Turn a SeQuant coupled-cluster term into a publication-style Goldstone diagram (PNG).

```
python3 csv_diagram.py '<term string>' out.png
```

Everything is driven by the **term string** that SeQuant prints — you paste it in, you
get a diagram out. The renderer understands the CSV/PNO–DF flavour of those strings
(PAO→PNO transforms `C`, PAO overlaps `s`, density-fitting auxiliary index `Κ`) and can
either draw them as annotations or drop them to leave a plain textbook CC diagram.

---

## 1. Requirements

Python 3 with `matplotlib` and `numpy`. Nothing else; no SeQuant build, no LaTeX
installation (matplotlib's built-in mathtext renders the labels).

```bash
python3 -c "import matplotlib, numpy; print('ok')"
```

---

## 2. Quick start

```bash
cd utilities/diagram_gen

# a) one term, given inline
python3 csv_diagram.py '2 (f{i_1;μ̃_1} * C{μ̃_1;a_1<i_1>}) * t{a_1<i_1>;i_1}' bubble.png

# b) one term, from a file
python3 csv_diagram.py --file term.txt out.png

# c) many terms: one file, one PNG per line -> term001.png, term002.png, ...
python3 csv_diagram.py --file sample_terms.txt outdir/

# d) the canonical CC skeleton instead of the CSV/DF overlay
python3 csv_diagram.py --file term.txt out.png --no-csv
```

`sample_terms.txt` in this folder holds four real terms to try; `examples/` shows what
they render to.

| example | what it shows |
|---|---|
| `examples/ccsd_E3_fock_bubble.png` | one-body `—×` vertex, T1, particle+hole drawn as a lens |
| `examples/ccsd_R2t49_vvvv_ladder.png` | (vv\|vv) ladder + T2, four `C` boxes, `K` on the interaction |
| `examples/ccsd_R2t30_oooo_csc_bridge.png` | (oo\|oo) + two T1, `CsC` domain-changing bridges |
| `examples/ccsdt_R1t1_with_T3.png` | a T3 amplitude bar (three excitation apexes) |
| `examples/ccsd_R2t49_canonical.png` | same term as the ladder above, with `--no-csv` |

---

## 3. The input: what a term string looks like

A term is a product of tensors, written `label{bra;ket}` (or `label{bra;ket;aux}`),
joined by `*`. Parentheses are irrelevant — the renderer ignores the association.

```
1/2 Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>} * g{μ̃_1;μ̃_2;Κ_1} * C{a_1<i_1,i_2>;μ̃_1} * … * t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}
     ↑ coefficient    ↑ external indices (projector)      ↑ interaction   ↑ PAO→PNO transform      ↑ amplitude
```

### Tensors understood

| label | slots | meaning | drawn as |
|---|---|---|---|
| `g{p;q;Κ}` | bra, ket, auxiliary | one half of a two-body interaction (a DF three-centre factor) | a dot; the two `g`'s sharing a `Κ` are joined by the dashed interaction line |
| `f{p;q}` | bra, ket | one-body (Fock) operator | a dot joined by a dashed line to a `×` |
| `t{virt…;occ…}` | bra=virtuals, ket=occupieds | cluster amplitude T₁…T₄ | solid horizontal bar at the bottom, one V-shaped particle/hole pair per excitation |
| `C{p;q}` | bra, ket | PAO→PNO transform | *not* a vertex — it splices two indices into one fermion line and puts a `C` box on it |
| `s{μ̃;μ̃}` | bra, ket | PAO–PAO overlap metric | same, part of a `CsC` bridge |

### Index naming (this is how line types are decided)

| starts with | kind | drawn as |
|---|---|---|
| `i…` | occupied | **hole** line — blue, arrow pointing down |
| `a…` | PNO virtual | **particle** line — red, arrow pointing up |
| `μ̃…` | PAO virtual | internal to a `C`/`s` chain; never a line end in the picture |
| `Κ…` | DF auxiliary | the `K` label on the dashed interaction line |

`a_3<i_1,i_2>` means "virtual `a_3` restricted to pair domain (i₁i₂)"; the `<…>` part is
printed as the superscript `a₃^{i₁i₂}`.

### What gets cleaned up automatically

all of these are stripped for you:

* `:N-C-S` (and similar) suffixes after a tensor
* markdown list numbering and backticks: ``12. `<term>` ``
* surrounding whitespace, a leading `+`
* any parenthesisation

So this works as a single argument, unchanged:

```
1. `-6 Ŝ{i_1,i_2,i_3;a_1<i_1,i_2,i_3>,…}:N-C-S * g{i_4;μ̃_1;Κ_1}:N-C-S * …`
```

### Rules the term must satisfy

1. **Every index appears exactly twice** — once in a bra slot, once in a ket slot.
   A truncated paste is the usual cause of failure here.
2. **At most one interaction**: either two `g`'s sharing one `Κ`, or one `f`.
   True of every CC residual term SeQuant emits.
3. Amplitudes may be T₁…T₄ (verified up to T₄; the bar just grows an apex per excitation).

---

## 4. `csv_diagram.py` — render one term (or a small batch)

```
python3 csv_diagram.py [TERM] [OUT.png] [options]
python3 csv_diagram.py --file PATH [OUT.png|OUTDIR/] [options]
```

| option | effect |
|---|---|
| `TERM` | the term string; use `-` to read it from stdin |
| `OUT.png` | output path (default `diagram.png`). With `--file` and several terms, this is a **directory** and files are named `term001.png`, … |
| `--file PATH` | read term(s) from a file (`-` = stdin). Every line containing `{…}` is rendered; other lines (comments, blanks) are skipped |
| `--title TEXT` | figure title. The coefficient is appended automatically |
| `--no-csv` | **canonical CC skeleton**: drop the `C`/`CsC`/`s` boxes, the domain superscripts and the `K` / `g_L,g_R` tags |
| `--no-dom` | drop the domain superscripts only (`a₃^{i₁i₂}` → `a₃`) |
| `--no-aux` | drop the `K` label and the `g_L`/`g_R` tags only |
| `--no-legend` | omit the legend strip at the bottom |

With no arguments at all it renders three built-in demo terms into the current directory.

### The overlay knob

`--no-csv` implements rule **B0** of `csv-cc-diagram-rules.md` (skeleton invariance):
deleting every domain superscript, `C`/`CsC` box and the `K` label must leave a textbook
Goldstone diagram. Compare `examples/ccsd_R2t49_vvvv_ladder.png` with
`examples/ccsd_R2t49_canonical.png` — same skeleton, no annotations. The three layers
(`--no-csv`, `--no-dom`, `--no-aux`) can be mixed freely.

---

## 5. `render_all.py` — batch a whole equation file

Renders every numbered term of an equation markdown file (sections `E`, `R1`, `R2`,
`R3`, `R4`) and writes an `INDEX.md` manifest next to the PNGs.

```bash
# whole CCSD file, CSV overlay + canonical skeleton
python3 render_all.py --src pao-ccsd-df.md --variant both

# just the triples section of the CCSDT file
python3 render_all.py --src pao-ccsdt-df.md --only R3

# three specific R4 terms into a scratch directory
python3 render_all.py --src pao-ccsdtq-df.md --only R4:1,2,3 --outdir /tmp/q
```

| option | effect |
|---|---|
| `--src FILE` | the equation markdown file. Defaults to `pao-ccsd-df.md` next to the script or one directory up |
| `--variant csv\|canonical\|both` | overlay, plain skeleton, or both (default `csv`) |
| `--only TAG[:N,N,…]` | restrict to a section (`R3`) or specific terms (`R2:16,49`) |
| `--outdir DIR` | output directory (single-variant runs). Default: `ccsd_all`, `ccsdt_all`, `ccsdtq_all`, … derived from the file name, with `_canonical` appended for the skeleton variant |
| `--kind TEXT` | title prefix (default derived from the file name, e.g. `CSV-CCSDT`) |

It expects terms formatted as a markdown numbered list under a `## E`/`## R1`/… heading:

```markdown
## R2 — doubles residual

1. `<term string>`
2. `<term string>`
```

Failures never abort the run — the offending term is listed in `INDEX.md` with its error.

---

## 6. `audit_diagrams.py` — check the output is clean

Re-runs the layout for every term of a file and reports anything unreadable: overlapping
labels, labels sitting on a vertex or a bar, labels on a foreign line, and fermion lines
drawn on top of each other.

```bash
python3 audit_diagrams.py --src pao-ccsdtq-df.md
python3 audit_diagrams.py --src pao-ccsd-df.md --no-csv
```

Over the CCSD/CCSDT/CCSDTQ equation files the renderer was developed against (660 terms,
both variants) the audit reports **no label/label, label/vertex or line/line overlaps**.
The only remaining reports are labels sitting on top of one *foreign* line — 1 of 86
CCSD terms, 8 of 201 CCSDT, 32 of 373 CCSDTQ. Those labels keep an opaque background,
so they stay readable.

---

## 7. Using it from Python

```python
import csv_diagram as cd

cd.draw(term_string, title="R2 term 49", out="out.png",
        show_csv=True,     # C / CsC / s boxes
        show_dom=None,     # domain superscripts   (None = follow show_csv)
        show_aux=None,     # K and g_L/g_R tags    (None = follow show_csv)
        legend=True, dpi=170, verbose=False)      # -> returns the output path
```

The parsing/graph layer is reusable on its own:

```python
coeff, tensors            = cd.parse_term(term_string)
nodes, lines, hvs, pair, ext = cd.build_graph(tensors)
for L in lines:
    print("hole" if L.is_hole else "particle",
          L.src, L.src_idx[0], "->", L.dst, L.dst_idx[0], L.annotation)
```

`build_graph` collapses every `C`/`s` chain, so `lines` is the physical fermion-line
graph: each `FermionLine` runs from the vertex where its index sits in a **bra** slot to
the one where it sits in a **ket** slot, and `.annotation` is `None`, `"C"` or `"CsC"`
depending on what the chain contained.

---

## 8. How to read the diagram

| element | meaning |
|---|---|
| red line, arrow up | particle line (virtual index) |
| blue line, arrow down | hole line (occupied index) |
| solid black bar | a cluster amplitude T_m; each excitation is one V of a particle + a hole leg meeting at an apex on the bar. The label under the bar is `t̂` with its occupied indices as subscripts and virtuals as superscripts |
| two grey dots joined by a dashed line | **one** two-body interaction vertex (its two half-vertices); `g_L`/`g_R` tag the left/right half |
| `K` on the dashed line | the density-fitting auxiliary index shared by the two `g` factors |
| dot + dashed line + `×` | the one-body (Fock) vertex |
| `C` box on a line | a PAO→PNO transform on that interaction leg |
| `CsC` box on a line | a domain-changing bridge C·s·C joining two PNOs of different domains |
| two lines bowed into a slim pointy ellipse | a particle and a hole line connecting the *same* two points — the canonical Goldstone bubble (they would otherwise be drawn on top of each other) |
| faint dotted line across the top | the projection boundary; lines crossing it are the external indices of `Ŝ` |

Conventions and their textbook citations are in `csv-cc-diagram-rules.md`. Sign/weight
rules are **not** applied — the coefficient printed in the title is whatever the term
carried.

---

## 9. When it fails

Errors are raised as `ValueError` with a specific message:

| message | cause |
|---|---|
| `no tensor found in input: …` | the string had no `{…}` — wrong argument, or the shell ate it (quote with single quotes) |
| `chain <idx> has N terminals (want 2)` | an index appears once, or three or more times — usually a truncated paste |
| `chain <idx>: both terminals are bra slots` | the same index used twice on the same side |
| `unknown index kind: 'x_1'` | an index name that starts with something other than `i`, `a`, `μ̃`, `Κ` |
| `>2 interaction half-vertices not supported` | more than one interaction in the term |

Shell tip: always wrap the term in **single** quotes — it contains `*`, `{`, `}` and
backticks, all of which the shell would otherwise interpret. If quoting gets awkward,
put the term in a file and use `--file`.

---

## 10. Files

| file | role |
|---|---|
| `csv_diagram.py` | the renderer: parser, fermion-line graph, layout, drawing, CLI. Self-contained |
| `render_all.py` | batch driver over an equation markdown file; writes `INDEX.md` |
| `audit_diagrams.py` | machine check for overlapping labels/lines (imports the other two) |
| `csv-cc-diagram-rules.md` | the diagram conventions with citations (Shavitt & Bartlett, Crawford & Schaefer) |
| `sample_terms.txt` | four real terms to try |
| `examples/` | what those terms render to |

`render_all.py` and `audit_diagrams.py` operate on equation markdown files
(`pao-ccsd-df.md`, `pao-ccsdt-df.md`, `pao-ccsdtq-df.md`) that are not part of this
repository; point them at your own with `--src`.
