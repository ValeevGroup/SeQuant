#ifndef SEQUANT_CORE_OPTIMIZE_GA_GENOME_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_GENOME_HPP

// Leaf-insertion codec for binary trees over n leaves, represented as laminar
// families of leaf-set bitmasks: a tree on {0..n-1} <-> 2n-1 masks (all
// singletons, all internal clusters, the full set). Gene k (k = 2..n) has
// range 2k-3; the code is a bijection onto the (2n-3)!! binary trees. Ported
// 1:1 from the Python prototype (proto_csv_core.py) -- the gene values index
// the SAME canonical ordering, so a genome decodes to the same tree in both.

#include <SeQuant/core/container.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <utility>
#include <vector>

namespace sequant::opt::ga {

/// Cluster of leaves as a bitmask. One codec runs both the per-term L1 trees
/// (leaves = a term's tensors) and the per-target L2 trees (leaves = a target's
/// summands), so the width must cover both; build_key_table rejects the rest.
using NodeMask = std::uint64_t;

/// Laminar family: sorted (see canon_less) vector of distinct masks.
using Laminar = container::svector<NodeMask>;

/// One integer gene per leaf insertion; gene k in [0, 2k-3).
using TreeCode = container::svector<int>;

/// Python's canon_order key: (|S|, tuple(sorted(S))) -- popcount first, then
/// lexicographic on the ascending element list, which at equal popcount is
/// decided by the LOWEST differing bit. Pinned by the "ga genome codec" test.
inline bool canon_less(NodeMask a, NodeMask b) {
  const int pa = std::popcount(a), pb = std::popcount(b);
  if (pa != pb) return pa < pb;
  const NodeMask d = a ^ b;
  return d && (a & (d & (~d + 1)));  // lowest differing bit belongs to a
}

/// **`canon_less` is a strict TOTAL order**: distinct masks never tie, so the
/// sorted image of a set of masks is UNIQUE and every merge/splice rewrite
/// below reproduces `canon_sort` exactly despite `std::sort` being unstable.
inline void canon_sort(Laminar& fam) {
  std::sort(fam.begin(), fam.end(), canon_less);
}

/// Python's `tuple(sorted(S))` order: lexicographic on the ascending element
/// list, WITHOUT the popcount tie first. Used for the child ordering inside
/// decoded trees and the producer-fibre ordering; distinct from canon_less.
inline bool lex_less(NodeMask a, NodeMask b) {
  while (a && b) {
    const int la = std::countr_zero(a), lb = std::countr_zero(b);
    if (la != lb) return la < lb;
    a &= a - 1;
    b &= b - 1;
  }
  return !a && b;
}

/// Insert leaf `l` (a single-bit mask absent from `fam`) as a sibling of the
/// subtree rooted at family member T: proper supersets of T widen by l, T and
/// its copies stay, {l} joins. KEEPS THE FAMILY CANON-SORTED: the 3-way merge
/// of the kept members, the widened supersets and `{l}` runs over three
/// disjoint already-sorted runs, so its output is exactly what sorting the
/// union gives. Pinned exhaustively for n <= 7 by the "ga genome codec" test.
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

/// Inverse of decode_tree; decode_tree(encode_tree(fam, n), n) == fam. `fam`
/// must be the canon-sorted laminar family of a binary tree on {0..n-1}.
/// Each peel is the exact mirror of `ins`: a 2-way merge of the members
/// without `sl` and the strict ancestors of P with `sl` cleared, two disjoint
/// canon-sorted runs, so it reproduces sort + unique element for element.
/// Pinned exhaustively for n <= 7 by the "ga genome codec" test.
inline TreeCode encode_tree(Laminar const& fam_in, int n) {
  TreeCode code(n - 1);
  // Ping-pong `std::vector` buffers so the bottom `std::swap` is a pointer
  // exchange (a small_vector's is not). Element values and order untouched.
  std::vector<NodeMask> fam(fam_in.begin(), fam_in.end()), next;
  next.reserve(fam.size());
  for (int k = n; k >= 2; --k) {
    const NodeMask sl = NodeMask{1} << (k - 1);
    // Parent of leaf k-1: the members containing `sl` form a chain under
    // inclusion (laminarity) and canon order ascends popcount along a chain,
    // so the first such member past the singleton is the popcount-minimal one.
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
/// RETURN ORDER IS PART OF THE CONTRACT -- canon-LARGER child first (the
/// descending scan's first hit is the canon-maximum over S's strict subsets,
/// and the cover mask then skips its subtree). `nni_each_slot`'s slot order,
/// hence which NNI move a drawn index names, is built on it.
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
using ChildTable = container::svector<std::array<NodeMask, 3>>;

/// The children table of an already-decoded family, in `fam` order. The one
/// definition: `decode_forest` and `TreeMemo` both go through here, so a memo
/// hit is byte-identical to a fresh decode.
inline void build_children(Laminar const& fam, ChildTable& ch) {
  ch.clear();
  for (NodeMask S : fam) {
    if (std::popcount(S) < 2) continue;
    auto [c1, c2] = tree_children(fam, S);
    if (lex_less(c2, c1)) std::swap(c1, c2);
    ch.push_back({S, c1, c2});
  }
}

/// `tree_children(fam, S)` read out of the ChildTable: same pair in the same
/// (canon-larger first) order. The table stores each pair lex-ordered, but
/// "canon-larger of the two" recovers the contractual first element because
/// `canon_less` is a strict total order. `S` must be an internal member.
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
/// cluster `W`. Exactly one mask changes, so this is a splice, not a re-sort:
/// `fam` is canon-sorted on entry (the codec's standing invariant), and
/// inserting W at its `lower_bound` gives the unique sorted sequence of the
/// new set. `q` is that `lower_bound` taken with A still in place.
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

/// THE nearest-neighbour-interchange enumeration, in one place. A *slot* is a
/// (grandparent B, child A) pair with both internal: for every member B of
/// `fam` with popcount >= 2, in `fam` order, its children (b1, b2) =
/// `tree_children(fam, B)` are visited as A = b1 (sibling A3 = b2) then A = b2,
/// keeping only those with popcount(A) >= 2. `fn(A, A3)` runs once per slot;
/// returning false stops the walk. `ch` must be `build_children(fam)`.
///
/// **This walk defines the move INDEX**: slot s owns 2s and 2s+1, keep running
/// over `tree_children(fam, A)` swapped (A2 first), which is what makes
/// `nni_move_count` + `nni_move_at` a valid count-draw-construct pair for
/// `nni_moves` + an index. Every entry point below runs THIS walk, and all are
/// pinned together for n <= 7 in the "ga genome codec" test.
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

/// Number of NNI neighbours of `fam` == `nni_moves(fam).size()`, in closed
/// form: every internal non-root node is one slot owning two moves, i.e.
/// 2(n-2). The search draws a kick's move index against this, so it must be
/// EXACT; pinned for n <= 7 in the "ga genome codec" test. Depending only on
/// the leaf count is what lets a kick be drawn before its block is decoded.
inline std::size_t nni_move_count(int n_leaves) {
  return n_leaves < 3 ? 0 : 2 * static_cast<std::size_t>(n_leaves - 2);
}
inline std::size_t nni_move_count(Laminar const& fam) {
  // |fam| = 2n - 1
  return nni_move_count(static_cast<int>((fam.size() + 1) / 2));
}

/// The k-th NNI neighbour of `fam`, k < `nni_move_count(fam)`, built without
/// materializing any of the others -- `nni_moves(fam)[k]` exactly. `ch` must
/// be `build_children(fam)`; A's children are read from it in `tree_children`
/// order (canon-larger first).
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
