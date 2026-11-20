# CSV-CC Goldstone diagram rules

Rule set for drawing and interpreting CSV-CCSD (cluster-specific-virtual, locally
density-fitted) Goldstone diagrams. Part A is the canonical textbook layer; Part B is
exactly what the CSV/DF overlay adds. Strip every Part-B annotation and a valid
CSV-CC diagram must reduce to the textbook diagram of Part A.

Citations:
- **C&S** = Crawford & Schaefer 2000, *An Introduction to Coupled Cluster Theory*,
  chapter "An Introduction to Coupled Cluster Diagrams", printed pp. 77–95
  (PDF page = printed − 32; formalism credited to Kucharski & Bartlett, p. 77).
- **S&B** = Shavitt & Bartlett 2009, *Many-Body Methods in Chemistry and Physics*,
  Ch. 4 (printed page = PDF − 16), plus the CC rules of Ch. 9–10 (Fig. 10.1, p. 296).
- **KST** = K. S. Thanthiriwatte, *Diagrammatic techniques in Coupled-Cluster
  Theory*, slides, Aug 2010 (18 slides).

---

## Part A — canonical skeleton rules

### A1. Line direction and labels
Upward-directed lines are **particle** lines (virtuals a, b, c, …); downward-directed
lines are **hole** lines (occupieds i, j, k, …). "Downward-directed lines represent
hole states (orbitals occupied in the reference) and upward directed lines represent
particle states" (C&S p. 78, Fig. 1; S&B p. 91: "pointing upward for particles and
downward for holes"; KST slide 4). The horizontal position of a line "has no
significance" (S&B p. 91). The reference |0⟩ is drawn as nothing (C&S Fig. 1c;
KST slide 4).

### A2. Vertex inventory
- **One-body operator (Fock)**: a dot joined by a horizontal *dashed* line to a
  marker "×" (C&S Fig. 2, p. 79; S&B p. 93 "the point of action … marked by the
  interaction line —×"; KST slide 5). Four fragments f_ab, f_ij, f_ia, f_ai with
  excitation levels 0, 0, −1, +1.
- **Two-body operator**: ONE vertex = **two half-vertices at the same level joined by
  a horizontal dashed interaction line** ("The two half-vertices and the interaction
  line constitute a single vertex", S&B p. 111, Eq. 4.28; C&S Fig. 3, p. 80: the nine
  V_N fragments; KST slide 6). Electron assignment: "electron 1 ↔ left half-vertex;
  electron 2 ↔ right half-vertex" (S&B p. 111, boxed).
- **Cluster operators T_m**: a **solid horizontal bar** with m hole–particle line
  pairs rising from it; no lines below the bar ("the cluster operators contain only
  q-creation strings … they contain no lines below the horizontal bar", C&S p. 80,
  Fig. 4; S&B pp. 273–274, Eqs. 9.84–9.87 "solid horizontal lines to represent the
  T̂_m operators"; KST slides 4, 8).

### A3. Half-vertex degree rule
"**Each individual half-vertex will have one incoming and one outgoing line**, each of
which may be a particle line or a hole line" (S&B p. 111, verbatim). This holds for
every interaction dot (one-body or half of two-body).

### A4. Bra = outgoing = creation; ket = incoming = annihilation
"incoming line ↔ annihilation operator ↔ ket state (IAK); outgoing line ↔ creation
operator ↔ bra state (OCB)" (S&B p. 98, mnemonics; also p. 95 rule 1 and the boxed
table on p. 111). Creation operators lie above the interaction line, annihilation
below (C&S p. 78; KST slide 4). Consequently: a particle line above a vertex with
arrow leaving it = created particle; a hole line below a vertex with arrow leaving
it = created hole; etc.

### A5. Integral assignment
- One-body: **f = ⟨out | f | in⟩** — "(out|f|in), where out indicates the index of the
  outgoing directed line and in the incoming" (C&S p. 84; KST slide 14 rule 2;
  S&B p. 93: bra = leaving line, ket = entering line).
- Two-body: **⟨left-out right-out ‖ left-in right-in⟩**, with operator product
  {(left-out)† (right-out)† (right-in)(left-in)} (S&B p. 112, boxed, verbatim;
  C&S p. 85; KST slide 6, verbatim). Goldstone diagrams proper use the
  non-antisymmetrized ⟨pq|rs⟩ (S&B p. 112); antisymmetrized-Goldstone (Brandow/ASG)
  and Hugenholtz variants use ⟨pq‖rs⟩ (S&B pp. 118–123). *The CSV-CCSD spin-traced
  equations are in the non-antisymmetrized (Goldstone-proper) regime.*

### A6. Amplitude assignment
"With every T_m vertex associate an amplitude t^{ab…}_{ij…}" where the (hole,
particle) pairs are read off the bar **left to right**, each vertical pair supplying
one (subscript, superscript) column (KST slide 14 rule 4; C&S pp. 84–85 "hole and
particle indices in their left-to-right order in the diagram"; S&B p. 282 and
Fig. 10.1 rule 4). External pairing matches the bra determinant: "surviving hole and
particle lines at the top are paired vertically in |Φ…⟩ if they are part of the same
path" (S&B p. 104).

### A7. Vertical placement (fictitious time)
Operators act bottom-to-top: rightmost operator (a T) at the bottom, Hamiltonian
fragment above it, resulting/bra state at the top (C&S p. 81, p. 83 "the rightmost
operator's interaction line must lie at the bottom"; S&B p. 91; KST slide 9).
T1 vertices commute, so the relative order/side-by-side placement of T bars is
immaterial (C&S p. 86; S&B p. 274, Eq. 9.88). External lines run to the top edge;
energy diagrams "can contain no external lines" (C&S p. 83); the R1/R2 equations have
exactly 2/4 external lines (C&S p. 88).

### A8. Internal summation
"Summations are included over all 'internal' indices — lines that begin and end at
operator interaction lines and do not extend to infinity above or below" (C&S p. 84;
KST slide 14 rule 5: "Sum over all internal line labels, i.e. lines terminating
below the last Ĥ_N"; S&B p. 100, Eq. 4.19).

### A9. Sign rule: (−1)^{h+l} with quasi-loops
Sign = (−1)^{h+l}, h = number of hole lines, l = number of loops. "A loop is a route
along a series of directed lines that either returns to its beginning or begins at
one external line and ends at another" (C&S p. 84). Open diagrams: "for open diagrams
fictitious external loop should be included" (KST slide 14 rule 6); S&B formalizes
this as **quasiloops** — "open lines … paired at the top according to the pairing
pattern in the corresponding Slater determinant" (S&B p. 135, Eq. 5.15; Fig. 10.1
rule 8, p. 296). S&B writes the equivalent form (−1)^{h−l} (p. 109, Eq. 4.26).

### A10. Equivalent lines → factor 1/2
"A pair of lines is equivalent if they connect the same pair of vertices in the same
direction. Each pair … contributes a factor 1/2" (S&B p. 120; C&S p. 85: lines
"beginning at the same operator interaction line and ending at the same interaction
line"; KST slide 15 rule 7: weight (1/2)^m). Groups of n ≥ 3 equivalent lines
(MBPT/UCC only) give 1/n! (C&S p. 85, footnote q).

### A11. Equivalent T vertices → factor 1/n!
"If there are n equivalent vertices in the diagram, they contribute a prefactor of
1/n!" (C&S p. 87). Two T vertices are equivalent when they have the same m and
"are connected to the interaction vertex with the same number … of particle lines
and the same number of hole lines" (S&B p. 287; KST slide 15 rule 9).

### A12. Permutation operators for external lines
"Each pair of unique, external hole or particle lines introduces a permutation
function P(pq)", P(pq) f(p,q) = f(p,q) − f(q,p); pairs originating from the *same*
interaction line are not unique (already antisymmetric) (C&S p. 91, Eq. 154 p. 76;
KST slide 15 rule 8; S&B Fig. 10.1 rule 10). *In the spin-traced biorthogonal
CSV equations this role is played by the projector Ŝ at the head of each term.*

### A13. Connected-cluster requirement
Every Hamiltonian fragment "must share at least one index with every cluster operator
on its right" (C&S p. 92); diagrams in which any T bar fails to contract with the
interaction are discarded. All T-vertex lines that do not reach the interaction must
be external.

### A14. Kucharski–Bartlett sign sequences (enumeration bookkeeping)
Unique connectivities are enumerated by assigning "+" to particle and "−" to hole
lines below the Hamiltonian interaction line / above the cluster bars, and combining
the signs in all unique ways, e.g. the five unique sequences of (V_N T1 T1 T2)_c
in the T2 equation (C&S pp. 93–95; same idea as KST slides 17–18 "signatures":
T1 = +−, T2 = ++−−, V fragments = ++−−, etc.).

---

## Part B — what CSV + DF adds (overlay annotations only)

**B0. Skeleton invariance.** The CSV-CC diagram is a canonical Goldstone diagram
(Part A) plus annotations. Deleting all domain superscripts, C / CsC boxes and the
K label must leave a textbook diagram (e.g. R2 term 43 → the (V_N T1^4)_C doubles
diagram, C&S pp. 93–95).

**B1. Domain superscripts.** A PNO virtual is written a^{ij} = virtual a restricted
to pair domain ⟨ij⟩ (SeQuant `a_1<i_1,i_2>`; pao-ccsd-df.md Notation). Particle-line
labels always carry their domain. Occupied indices are ordinary MOs (no domain).

**B2. DF split of the interaction (K).** Under density fitting
g_{pq,rs} = Σ_K B^K_{pq} B^K_{rs}: each half-vertex is one three-center factor
`g{bra;ket;Κ}` and **the dashed interaction line between the two half-vertices
carries the auxiliary label K** (write just "K"). The two `g` factors of a term are
identified by their shared Κ index. A term with `f{...}` instead has the one-body
×-vertex and no K.

**B3. Slot → line mapping.** Every SeQuant tensor prints as `label{bra;ket[;aux]}`.
Apply A4 directly: **bra slot = outgoing line, ket slot = incoming line** at that
vertex. For amplitudes `t{virtuals;occupieds}`: virtuals are bra (created particles,
arrows leaving the bar upward), occupieds are ket (created holes: downward arrows
entering the bar from above). For the projector `Ŝ{occ;virt}` the occupieds are
bra-side (external hole lines flow from the top edge down into the diagram) and the
virtuals ket-side (external particle lines flow up out of the diagram).

**B4. C annotation (PAO→PNO transform on an interaction leg).** g and f store their
*virtual* indices in the PAO basis (μ̃); a factor `C{μ̃;a<dom>}` (or transposed)
converts that leg to a PNO. Draw ONE fermion line for the chain g—C—a and put a
boxed **C** on the arrow near the g end: g —C→ a^{ij}. Occupied legs of g/f are MO
indices — never annotated. A leg with C is always a particle line.

**B5. CsC annotation (domain-changing bridge).** A chain
`C{a<dom1>;μ̃} · s{μ̃;μ̃'} · C{μ̃';a'<dom2>}` (s = the PAO–PAO overlap metric) joins
two PNO indices of *different* domains into one particle line. Draw one line and put
a boxed **CsC** at mid-height. Typical use: connecting an external a_1^{i1 i2} to an
amplitude virtual a_5^{i3}.

**B6. Fermion lines = connector chains.** Formally: C and s are not vertices. A
fermion line is a maximal chain of indices linked through C and s factors whose two
endpoints are slots of *terminal* tensors (g, f, t, Ŝ) — always exactly one bra slot
and one ket slot (this is the CSV image of A3/A4). Chain content ⇒ annotation:
one C → "C"; C·s·C → "CsC"; none → bare line.

**B7. Amplitude vertex label.** `\hat{t}_{occ…}^{csv-virt…}` under the bar, e.g.
t{a_3<i_1>;i_1} → t̂_{i1}^{a3^{i1}}; t{a_1<i1,i2>,a_2<i1,i2>;i_1,i_2} →
t̂_{i1 i2}^{a1^{i1 i2} a2^{i1 i2}}. Amplitudes are kept in CSV/PNO form — never
expanded into PAO (pao-ccsd-df.md Notation) — so C/CsC boxes never sit *between* a
T bar and its own index label; they sit on the line above it.

**B8. External lines and Ŝ.** The external indices are exactly those listed in
`Ŝ{i…;a<dom>…}`; Ŝ is the biorthogonal particle-symmetrizer/projector marking
R = ⟨P|H̄|0⟩ — not a contracted tensor (pao-ccsd-df.md Notation). It replaces the
P(pq) machinery of A12 in the spin-traced biorthogonal formulation.

**B9. Signs/weights.** The printed rational coefficient of each SeQuant term already
absorbs the (−1)^{h+l} sign (A9), the equivalent-line/vertex weights (A10/A11) and
the closed-shell spin-tracing factors; the diagram is drawn for structure, and the
A9–A11 count can be used as a consistency check (e.g. term 43: h = 4 hole lines,
l = 2 quasiloops pairing (i1,a1),(i2,a2) ⇒ (−1)^6 = +1, matching its coefficient +1).

---

## Quick checklist for auditing a drawn CSV-CC diagram

1. Every interaction half-vertex: exactly 1 in + 1 out (A3).
2. Every T1 bar: one particle up-arrow leaving + one hole down-arrow entering;
   T_m bar: m of each (A2/A6).
3. Arrow directions: particles up, holes down, everywhere along a path (A1).
4. bra/ket of every SeQuant slot honored as out/in (A4/B3).
5. External lines = Ŝ indices only, run to the top edge, paired (i_k, a_k)
   per the bra determinant (A7/B8).
6. Connectedness: every T bar shares ≥1 line with the interaction (A13).
7. K on the dashed line iff the term has two `g` factors sharing a Κ (B2).
8. C on every g/f particle leg; CsC on every C·s·C bridge; domains on every
   particle label; no annotation on occupied lines (B1, B4, B5).
9. Strip annotations → a textbook Goldstone diagram (B0).
