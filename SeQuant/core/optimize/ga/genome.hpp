#ifndef SEQUANT_CORE_OPTIMIZE_GA_GENOME_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_GENOME_HPP

// Leaf-insertion codec for binary trees over n leaves, represented as laminar
// families of leaf-set bitmasks. A tree on leaves {0..n-1} <-> a laminar
// family of 2n-1 masks (all singletons, all internal clusters, the full set).
// Gene k (k = 2..n) has range 2k-3; the code is a bijection onto the
// (2n-3)!! binary trees. Ported 1:1 from the Python prototype
// (proto_csv_core.py: canon_order/ins/decode_tree/encode_tree/children/
// nni_moves) -- the gene values index into the SAME canonical ordering, so a
// genome decodes to the same tree in both implementations.

#include <SeQuant/core/container.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <utility>
#include <vector>

namespace sequant::opt::ga {

/// Cluster of leaves as a bitmask. The same codec runs the per-term L1
/// trees (leaves = one term's tensors, CCSD needs n <= 14) AND the
/// per-target L2 trees (leaves = a target's summands: CSV-CCSD R2 has 40),
/// so the width must cover max(factors per term, summands per target);
/// build_key_table rejects inputs that don't fit.
using NodeMask = std::uint64_t;

/// Laminar family: sorted (see canon_less) vector of distinct masks.
using Laminar = container::svector<NodeMask>;

/// One integer gene per leaf insertion; gene k in [0, 2k-3).
using TreeCode = container::svector<int>;

/// Python's canon_order key: (|S|, tuple(sorted(S))) -- popcount first, then
/// lexicographic on the ascending element list (= successive lowest set bits).
///
/// At equal popcount the lexicographic comparison is decided entirely by the
/// LOWEST differing bit l: below l the two element lists agree, and whichever
/// mask owns l contributes l at that position while the other contributes
/// something strictly greater. So a < b iff l belongs to a. Checked
/// exhaustively against the bit-stripping loop it replaces over all pairs of
/// 10-bit masks in the "ga genome codec" test.
inline bool canon_less(NodeMask a, NodeMask b) {
  const int pa = std::popcount(a), pb = std::popcount(b);
  if (pa != pb) return pa < pb;
  const NodeMask d = a ^ b;
  return d && (a & (d & (~d + 1)));  // lowest differing bit belongs to a
}

/// `std::sort` is unstable, which is harmless here and is what lets the merge
/// rewrites below (T-A2, T-A2b) be exact: `canon_less` is a strict TOTAL order
/// on masks, not merely a weak one. For a != b of equal popcount the XOR is
/// nonzero and its lowest set bit belongs to exactly one of them, so exactly
/// one of canon_less(a,b) / canon_less(b,a) holds; unequal popcounts are
/// ordered by popcount. Hence two masks are canon_less-equivalent iff they are
/// the same value, the sorted image of a multiset of masks is unique, and any
/// procedure that emits the elements in nondecreasing order emits exactly what
/// `canon_sort` would -- there is no permutation freedom among "ties".
inline void canon_sort(Laminar& fam) {
  std::sort(fam.begin(), fam.end(), canon_less);
}

/// Python's `tuple(sorted(S))` order: lexicographic on the ascending element
/// list, WITHOUT the popcount tie applied first (a prefix sorts before its
/// extensions). Used for the child ordering inside decoded trees and the
/// producer-fibre ordering; distinct from canon_less.
inline bool lex_less(NodeMask a, NodeMask b) {
  while (a && b) {
    const int la = std::countr_zero(a), lb = std::countr_zero(b);
    if (la != lb) return la < lb;
    a &= a - 1;
    b &= b - 1;
  }
  return !a && b;
}

/// Insert leaf `l` (a single-bit mask, absent from every member of `fam`) as a
/// sibling of the subtree rooted at family member T: every proper superset of T
/// widens by l, T and its copies stay, and {l} joins. Keeps the family
/// canon-sorted.
///
/// The output is the union of three runs, each already canon-sorted, so it is
/// produced by a 3-way merge rather than a full sort (T-A2 -- the old
/// canon_sort was Theta(n^2 log n) comparisons per decoded tree, and the L2
/// block does 54 insertions each sorting up to 109 masks):
///
///  (a) the kept members -- a subsequence of `fam`, hence still sorted;
///  (b) the widened supersets `C | l` for `C` superset of `T`, emitted in
///      `fam` order;
///  (c) the singleton `{l}`.
///
/// Run (b) is sorted because `C -> C | l` is order-preserving for `canon_less`
/// whenever `l` is a bit no `C` already has. Proof: popcount(C|l) =
/// popcount(C)+1 for every member of the run, so the popcount comparison is
/// unchanged; and at equal popcount `canon_less` looks only at the lowest bit
/// of the XOR, with (C|l)^(C'|l) = C^C' since bit l cancels -- so both the
/// deciding bit and its owner are the same as for C vs C'. (For this codec `l`
/// is in fact the new highest bit, which also says the widened masks merely
/// append the same top element to each ascending element list.)
///
/// The three runs are disjoint: (a) has no member containing `l`, (b)'s all
/// contain `l` and have popcount >= 2, and (c) has popcount 1 -- so no ties
/// arise and the merge reproduces the sort exactly.
inline void ins(Laminar& fam, NodeMask T, NodeMask l) {
  Laminar out;
  out.reserve(fam.size() + 2);
  const auto end = fam.end();
  auto a = fam.begin(), b = fam.begin();  // cursors over runs (a) and (b)
  const auto skip_a = [&] {               // (a): skip proper supersets of T
    while (a != end && (*a & T) == T && *a != T) ++a;
  };
  const auto skip_b = [&] {  // (b): skip non-supersets of T
    while (b != end && (*b & T) != T) ++b;
  };
  skip_a();
  skip_b();
  bool have_l = true;  // run (c), one element
  while (a != end || b != end || have_l) {
    NodeMask best = 0;
    int src = -1;
    if (a != end) {
      best = *a;
      src = 0;
    }
    if (b != end) {
      const NodeMask w = *b | l;
      if (src < 0 || canon_less(w, best)) {
        best = w;
        src = 1;
      }
    }
    if (have_l && (src < 0 || canon_less(l, best))) {
      best = l;
      src = 2;
    }
    out.push_back(best);
    if (src == 0) {
      ++a;
      skip_a();
    } else if (src == 1) {
      ++b;
      skip_b();
    } else {
      have_l = false;
    }
  }
  fam = std::move(out);
}

/// Decode `code` into the laminar family of a binary tree on `n` leaves.
inline Laminar decode_tree(TreeCode const& code, int n) {
  Laminar fam{NodeMask{1}};
  for (int k = 2; k <= n; ++k) ins(fam, fam[code[k - 2]], NodeMask{1} << (k - 1));
  return fam;
}

/// Inverse of decode_tree; decode_tree(encode_tree(fam, n), n) == fam.
///
/// `fam` must be the (canon-sorted) laminar family of a binary tree on leaves
/// {0..n-1}: all n singletons, all internal clusters, the full set, 2n-1
/// distinct masks. Every caller supplies one -- `decode_tree` / `ins` output,
/// an `nni_moves` neighbour, or `sequence_to_family`'s postfix walk.
///
/// The peeling step is the exact mirror of `ins`, so it is a 2-way merge
/// rather than a sort + unique (T-A2b -- one canon_sort of up to 2n-1 masks
/// per leaf, i.e. Theta(n^2 log n) comparisons per encode, and `hill_climb`
/// encodes every NNI neighbour it tries). `next` is the union of two runs:
///
///  (a) the members that do not contain `sl`, passed through unchanged -- a
///      subsequence of `fam`, hence still sorted;
///  (b) `C & ~sl` for the members that DO contain `sl` other than `sl` itself
///      and P, emitted in `fam` order.
///
/// Run (b) is sorted by the mirror of T-A2's lemma: the members containing
/// `sl` form a chain ({sl} < P < Q1 < ... under inclusion), so (b) is exactly
/// the STRICT supersets of P; clearing `sl` drops every popcount by 1, leaving
/// the popcount comparison unchanged, and at equal popcount `canon_less` reads
/// only the lowest bit of the XOR, with (Q&~sl)^(Q'&~sl) = Q^Q' since bit `sl`
/// cancels -- same deciding bit, same owner.
///
/// The two runs are disjoint, so no ties arise, `std::unique` had nothing to
/// remove, and the merge reproduces sort + unique exactly. Proof: T = P & ~sl
/// is P's sibling child of `sl`, so T is a member; every element of (b) is
/// Q & ~sl with Q a strict superset of P, and since Q and P both contain `sl`
/// this strictly contains P & ~sl = T. A collision with run (a) would then be
/// a member D strictly containing T with `sl` not in D. But `fam` is laminar,
/// so the members containing T are T's ancestors -- T, then P, then P's
/// ancestors -- and every proper one contains `sl`. Contradiction. (Within a
/// run there are no duplicates either: `fam`'s masks are distinct and
/// C -> C & ~sl is injective on masks that all contain `sl`.)
///
/// `canon_less` is a total order (see `canon_sort`), so the sorted sequence of
/// this multiset is unique and the merge output is element-for-element what
/// the sort produced.
inline TreeCode encode_tree(Laminar const& fam_in, int n) {
  TreeCode code(n - 1);
  // Ping-pong `std::vector` working buffers: `std::swap` at the bottom is a
  // pointer exchange, so the n-1 peels allocate twice total instead of once
  // each (a `Laminar` is a small_vector, whose swap/move cannot pointer-swap
  // around the inline storage). Element values and order are untouched.
  std::vector<NodeMask> fam(fam_in.begin(), fam_in.end()), next;
  next.reserve(fam.size());
  for (int k = n; k >= 2; --k) {
    const NodeMask sl = NodeMask{1} << (k - 1);
    // Parent of leaf k-1: the members containing bit `sl` form a CHAIN under
    // inclusion (laminarity), and along a chain canon order is strictly
    // ascending popcount -- so the FIRST member met that contains `sl`,
    // besides the singleton itself, IS the popcount-minimal one the old
    // whole-family scan selected. (Uniqueness: nested distinct masks cannot
    // tie on popcount, so "first" and "strict minimum" pick the same member.)
    NodeMask P = 0;
    for (NodeMask C : fam)
      if ((C & sl) && C != sl) {
        P = C;
        break;
      }
    const NodeMask T = P & ~sl;
    next.clear();
    const auto end = fam.end();
    auto a = fam.begin(), b = fam.begin();  // cursors over runs (a) and (b)
    const auto skip_a = [&] {               // (a): members without `sl`
      while (a != end && (*a & sl)) ++a;
    };
    const auto skip_b = [&] {  // (b): the strict ancestors of P
      while (b != end && (!(*b & sl) || *b == sl || *b == P)) ++b;
    };
    skip_a();
    skip_b();
    while (a != end || b != end) {
      if (b == end || (a != end && canon_less(*a, *b & ~sl))) {
        next.push_back(*a);
        ++a;
        skip_a();
      } else {
        next.push_back(*b & ~sl);
        ++b;
        skip_b();
      }
    }
    const auto it = std::lower_bound(next.begin(), next.end(), T, canon_less);
    code[k - 2] = static_cast<int>(it - next.begin());
    std::swap(fam, next);
  }
  return code;
}

/// The two children of internal cluster S: its maximal proper subsets in fam.
///
/// RETURN ORDER IS PART OF THE CONTRACT: the pair comes back canon-LARGER
/// child first. The descending scan's first hit is the canon-maximum over all
/// of S's strict subsets in `fam`, and every such subset is a descendant,
/// canon-less than the child above it (a strict subset has strictly smaller
/// popcount) -- so the maximum is the canon-larger child; the cover mask then
/// skips that child's whole subtree, and the same argument hands back the
/// other child next. `nni_each_slot`'s slot order is built on this.
inline std::pair<NodeMask, NodeMask> tree_children(Laminar const& fam,
                                                   NodeMask S) {
  NodeMask c1 = 0, c2 = 0;
  // fam is canon-sorted (ascending popcount), so scanning downward meets each
  // child before any of that child's descendants; skip masks covered so far.
  for (auto it = fam.rbegin(); it != fam.rend(); ++it) {
    const NodeMask C = *it;
    if ((C & S) != C || C == S || (C & ~(c1 | c2)) != C) continue;
    (c1 ? c2 : c1) = C;
    if ((c1 | c2) == S) break;
  }
  return {c1, c2};
}

/// Per internal cluster of a decoded tree: (S, c1, c2) with c1 < c2
/// lexicographically, in `fam` order. A tree on n leaves has exactly n-1.
/// (Lives here, next to the codec, because the NNI walk below reads it; the
/// decode memo and `ForestState` in fitness.hpp consume the same definition.)
using ChildTable = container::svector<std::array<NodeMask, 3>>;

/// The children table of an already-decoded family, in `fam` order. The
/// canonical definition -- `decode_forest` and `TreeMemo` both go through here
/// so a memo hit is byte-identical to a fresh decode.
inline void build_children(Laminar const& fam, ChildTable& ch) {
  ch.clear();
  for (NodeMask S : fam) {
    if (std::popcount(S) < 2) continue;
    auto [c1, c2] = tree_children(fam, S);
    if (lex_less(c2, c1)) std::swap(c1, c2);
    ch.push_back({S, c1, c2});
  }
}

/// `tree_children(fam, S)` read out of the ChildTable instead of scanned for:
/// same pair in the same (canon-larger first) order, O(log n).
///
/// The table stores each pair lex-ordered, but the raw order is recoverable
/// without touching `fam`: `tree_children` returns the canon-larger child
/// first (see its contract above), and `canon_less` is a strict total order,
/// so "canon-larger of the two" names the same mask however the pair is
/// stored. `S` must be an internal member of the family the table was built
/// from; the table is in fam order, whose restriction to internals is still
/// canon-sorted, hence the binary search.
inline std::pair<NodeMask, NodeMask> tree_children(ChildTable const& ch,
                                                   NodeMask S) {
  auto const it = std::lower_bound(
      ch.begin(), ch.end(), S,
      [](std::array<NodeMask, 3> const& e, NodeMask s) {
        return canon_less(e[0], s);
      });
  const NodeMask c1 = (*it)[1], c2 = (*it)[2];
  return canon_less(c1, c2) ? std::pair{c2, c1} : std::pair{c1, c2};
}

/// The NNI neighbour of `fam` in which family member `A` is replaced by the
/// cluster `W`.
///
/// Each neighbour differs from `fam` in exactly ONE mask -- A leaves, W =
/// keep | A3 joins -- so it is built by splicing rather than by re-sorting
/// (T-A2b; `hill_climb` rebuilds the whole move list for every block of every
/// sweep, and the 55-leaf L2 block sorts 109 masks per move). `fam` is
/// canon-sorted on entry (the codec's standing invariant: `decode_tree` /
/// `ins` emit sorted families, as does this function), so:
///
///   * `fam` with A removed is a subsequence of a sorted sequence, hence
///     sorted;
///   * inserting W at `lower_bound(W)` of that subsequence therefore yields
///     the sorted sequence of the same multiset;
///   * and because `canon_less` is a total order (see `canon_sort`), that
///     sorted sequence is UNIQUE -- so it is element-for-element what the old
///     `canon_sort` produced, whatever permutation the unstable sort chose.
///
/// The index arithmetic below just fuses the erase and the insert into one
/// shift. `q` is `lower_bound(W)` taken over `fam` WITH A still in it: if A
/// precedes W (q > p) then dropping A slides W's target down to q-1 and the
/// span (p, q) moves one slot left; otherwise indices below p are untouched
/// and the span [q, p) moves one slot right.
inline Laminar nni_neighbour(Laminar const& fam, NodeMask A, NodeMask W) {
  Laminar next = fam;
  const auto beg = next.begin();
  const auto p = std::find(beg, next.end(), A) - beg;
  const auto q = std::lower_bound(beg, next.end(), W, canon_less) - beg;
  if (q > p) {
    std::move(beg + p + 1, beg + q, beg + p);
    next[q - 1] = W;
  } else {
    std::move_backward(beg + q, beg + p, beg + p + 1);
    next[q] = W;
  }
  return next;
}

/// THE nearest-neighbour-interchange enumeration, in one place (T-A2c).
///
/// A *slot* is a (grandparent B, child A) pair with both internal: for every
/// member B of `fam` with popcount >= 2, in `fam` order, its children
/// (b1, b2) = `tree_children(fam, B)` are visited as A = b1 (sibling A3 = b2)
/// then A = b2 (sibling A3 = b1), keeping only those with popcount(A) >= 2.
/// `fn(A, A3)` is called once per slot; returning false stops the walk.
///
/// Slot s owns the two consecutive move indices 2s and 2s+1, because the move
/// a slot produces is `A` replaced by `keep | A3` for `keep` running over
/// `{A2, A1}` = `tree_children(fam, A)` swapped -- A2 first. So, writing S for
/// the number of slots:
///
///     nni_moves(fam).size() == 2 * S
///     nni_moves(fam)[k]     == slot k/2 with keep = (k odd ? A1 : A2)
///
/// The walk is driven by the family's ChildTable (`ch` must be
/// `build_children(fam)` -- `nni_kick` and `hill_climb` get it from the
/// decode memo for free), which replaces one O(|fam|) `tree_children` scan
/// per internal member with a table read. Equivalence to the scan-driven walk
/// it replaced, piece by piece:
///
///   * `ch`'s entries are exactly the internal members, in `fam` order --
///     `build_children` filters the same popcount >= 2 test out of the same
///     iteration, so the B sequence is unchanged;
///   * per B, the scan produced (b1, b2) canon-larger child first (the
///     `tree_children` contract); the table stores the pair lex-ordered, and
///     "canon-larger of the two" recovers the same first element because
///     `canon_less` is a strict total order on distinct masks.
///
/// This is why `nni_move_count` (closed form) and `nni_move_at` may be used
/// as a count-draw-construct pair in place of `nni_moves` + an index into it,
/// and why the fam-only overloads below may delegate here through a freshly
/// built table: every entry point runs the SAME walk, and all of them are
/// pinned together against the pre-T-A2b reference enumeration exhaustively
/// for n <= 7 in the "ga genome codec" test.
template <typename Fn>
inline void nni_each_slot(ChildTable const& ch, Fn&& fn) {
  for (auto const& e : ch) {
    const NodeMask b1 = canon_less(e[1], e[2]) ? e[2] : e[1];
    const NodeMask b2 = b1 == e[1] ? e[2] : e[1];
    for (int i = 0; i < 2; ++i) {
      const NodeMask A = i ? b2 : b1, A3 = i ? b1 : b2;
      if (std::popcount(A) < 2) continue;
      if (!fn(A, A3)) return;
    }
  }
}

/// Number of NNI neighbours of `fam`, without building any of them.
/// Equals `nni_moves(fam).size()` -- see `nni_each_slot`. This is the value
/// `nni_kick` draws its move index below, so it has to be exact, and it is
/// pinned against the pre-T-A2b reference body exhaustively for n <= 7 in the
/// "ga genome codec" test.
///
/// Closed form, no walk: a slot is an (internal B, internal child A) pair,
/// and every internal node except the root is the child of exactly one
/// internal node (all parents are internal), so on n >= 2 leaves there are
/// exactly (n-1) - 1 slots -- the internal non-root nodes -- each owning two
/// moves. With |fam| = 2n - 1 that is 2(n-2) = |fam| - 3. The walk this
/// replaces was O(n^2) (one O(|fam|) `tree_children` scan per internal node)
/// and ran once per NNI kick inside `ga_once`'s SERIAL breeding phase.
///
/// The leaf-count form below is the same constant BEFORE any decode: it is
/// what lets `ga_once` draw a kick's move index without materializing the
/// kicked family at all (the count is the one family-dependent value the
/// rng stream used to need).
inline std::size_t nni_move_count(int n_leaves) {
  return n_leaves < 3 ? 0 : 2 * static_cast<std::size_t>(n_leaves - 2);
}
inline std::size_t nni_move_count(Laminar const& fam) {
  // |fam| = 2n - 1
  return nni_move_count(static_cast<int>((fam.size() + 1) / 2));
}

/// The k-th NNI neighbour of `fam`, k < `nni_move_count(fam)`, built without
/// materializing any of the others -- `nni_moves(fam)[k]` exactly (T-A2c).
/// `ch` must be `build_children(fam)`; A's children are read from it in
/// `tree_children` order (canon-larger first), see the ChildTable overload of
/// `tree_children`.
///
/// `nni_kick` draws one move out of the ~104 the 55-leaf L2 block offers and
/// throws the rest away; building all of them was 61 % of the serial breeding
/// phase, which is what Amdahl-caps T-C3's parallel evaluation.
inline Laminar nni_move_at(Laminar const& fam, ChildTable const& ch,
                           std::size_t k) {
  Laminar out;
  std::size_t s = 0;
  nni_each_slot(ch, [&](NodeMask A, NodeMask A3) {
    if (s++ != k / 2) return true;
    const auto [A1, A2] = tree_children(ch, A);
    out = nni_neighbour(fam, A, ((k & 1) ? A1 : A2) | A3);
    return false;
  });
  return out;
}

/// The fam-only form: builds the table and takes the same walk.
inline Laminar nni_move_at(Laminar const& fam, std::size_t k) {
  ChildTable ch;
  build_children(fam, ch);
  return nni_move_at(fam, ch, k);
}

/// All NNI neighbours of `fam`, in the order `nni_move_at` indexes.
/// (The list's ORDER is part of the result: `hill_climb` breaks ties
/// first-wins on a strict `<`.) `ch` must be `build_children(fam)`.
inline container::svector<Laminar> nni_moves(Laminar const& fam,
                                             ChildTable const& ch) {
  container::svector<Laminar> out;
  nni_each_slot(ch, [&](NodeMask A, NodeMask A3) {
    const auto [A1, A2] = tree_children(ch, A);
    for (NodeMask keep : {A2, A1})  // swap the other child out
      out.push_back(nni_neighbour(fam, A, keep | A3));
    return true;
  });
  return out;
}

/// The fam-only form: builds the table and takes the same walk.
inline container::svector<Laminar> nni_moves(Laminar const& fam) {
  ChildTable ch;
  build_children(fam, ch);
  return nni_moves(fam, ch);
}

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_GENOME_HPP
