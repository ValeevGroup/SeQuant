#include <SeQuant/domain/mbpt/spinor.hpp>
#include <regex>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <bit>
#include <cstdint>
#include <numeric>
#include <optional>
#include <string_view>
#include <utility>

namespace sequant::mbpt {

namespace {

/// Expands the bra antisymmetry of every tensor named @p label into the
/// NonSymm form, leaving all other tensors untouched. Intended to be run on a
/// Kramers-FREE expression (all indices spin-`any`): `expand_antisymm`'s
/// `ms_conserving_columns` guard then keeps every permutation (no Ms filter),
/// which is what we want — Kramers is not an Ms-conserved label.
ExprPtr expand_label_antisymm(const ExprPtr& expr, std::wstring_view label) {
  if (expr->is<Tensor>()) {
    const auto& t = expr->as<Tensor>();
    return t.label() == label ? expand_antisymm(t) : expr;
  } else if (expr->is<Product>()) {
    const auto& p = expr->as<Product>();
    auto result = std::make_shared<Product>();
    result->scale(p.scalar());
    for (const auto& factor : p) {
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == label) {
        result->append(1, expand_antisymm(factor->as<Tensor>()),
                       Product::Flatten::No);
      } else {
        result->append(1, factor, Product::Flatten::No);
      }
    }
    ExprPtr r = result;
    expand(r);  // distribute the antisymmetrizer-expansion Sum factor
    rapid_simplify(r);
    return r;
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (const auto& summand : *expr) {
      result->append(expand_label_antisymm(summand, label));
    }
    return result;
  }
  return expr;
}

/// Returns the first antisymmetrizer (Â) tensor found anywhere in @p e, or
/// nullopt. In a CC residual term Â is a leading factor whose bra/ket are the
/// external (virtual/occupied) index groups.
std::optional<Tensor> find_antisymmetrizer(const ExprPtr& e) {
  if (e->is<Tensor>()) {
    const auto& t = e->as<Tensor>();
    if (t.label() == reserved::antisymm_label()) return t;
    return std::nullopt;
  } else if (e->is<Product>()) {
    for (const auto& f : e->as<Product>()) {
      if (auto r = find_antisymmetrizer(f)) return r;
    }
    return std::nullopt;
  } else if (e->is<Sum>()) {
    for (const auto& s : *e) {
      if (auto r = find_antisymmetrizer(s)) return r;
    }
    return std::nullopt;
  }
  return std::nullopt;
}

/// Applies a bit permutation @p perm (bit p -> bit perm[p]) to config @p m.
std::uint64_t apply_bit_perm(const container::svector<std::size_t>& perm,
                             std::uint64_t m, std::size_t n) {
  std::uint64_t out = 0;
  for (std::size_t p = 0; p < n; ++p)
    if ((m >> p) & 1u) out |= (std::uint64_t{1} << perm[p]);
  return out;
}

}  // namespace

container::svector<container::svector<std::size_t>> kramers_external_generators(
    std::size_t rank) {
  const std::size_t n = 2 * rank;
  container::svector<container::svector<std::size_t>> gens;
  // S_rank adjacent transpositions of each external group: virtual [0,rank),
  // occupied [rank,n).
  auto add = [&](std::size_t lo, std::size_t hi) {
    for (std::size_t p = lo; p + 1 < hi; ++p) {
      container::svector<std::size_t> perm(n);
      std::iota(perm.begin(), perm.end(), std::size_t{0});
      std::swap(perm[p], perm[p + 1]);
      gens.push_back(std::move(perm));
    }
  };
  add(0, rank);
  add(rank, n);
  return gens;
}

container::svector<container::svector<std::uint64_t>> kramers_config_orbits(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& bit_perms,
    bool use_T) {
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;

  container::svector<container::svector<std::uint64_t>> orbits;
  container::set<std::uint64_t> seen;
  for (std::uint64_t m = 0; m <= full_mask; ++m) {
    if (seen.find(m) != seen.end()) continue;

    // BFS the orbit of m under the generators.
    container::set<std::uint64_t> orbit{m};
    container::svector<std::uint64_t> frontier{m};
    while (!frontier.empty()) {
      container::svector<std::uint64_t> next;
      for (const std::uint64_t x : frontier) {
        container::svector<std::uint64_t> images;
        for (const auto& perm : bit_perms)
          images.push_back(apply_bit_perm(perm, x, n));
        if (use_T) images.push_back((~x) & full_mask);
        for (const std::uint64_t y : images)
          if (orbit.insert(y).second) next.push_back(y);
      }
      frontier = std::move(next);
    }

    // container::set iterates ascending, so front() is the orbit minimum.
    container::svector<std::uint64_t> members(orbit.begin(), orbit.end());
    for (const std::uint64_t x : members) seen.insert(x);
    orbits.push_back(std::move(members));
  }
  return orbits;
}

container::svector<KramersBlock> kramers_external_blocks(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& antisym_perms,
    bool use_T,
    const container::svector<container::svector<std::size_t>>& symm_perms) {
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;

  using Perm = container::svector<std::size_t>;
  Perm identity(n);
  std::iota(identity.begin(), identity.end(), std::size_t{0});

  // transform from canonical to a config: block(cfg) = sign*[conj]*perm(canon),
  // where perm acts as block(cfg)[v] = block(canon)[v_perm] (axis q of the
  // output reads axis perm[q] of the input — the same convention kr_recon /
  // kramers_leaf reconstruct with: `blk(std) = rep(perm)`).
  struct Xform {
    int sign;
    bool conj;
    Perm perm;
  };
  // compose: apply generator g (cfg x -> y) on top of transform tx (canon ->
  // x). A generator acts as block(y) = g_sign*[g_conj]*g_perm(block(x)).
  // Substituting block(x) = tx_sign*[tx_conj]*tx_perm(block(canon)) and using
  // (p2 . p1)[k] = p2[p1[k]] for the perm convention above gives the COMPOSED
  // perm g_perm . tx_perm, i.e. ty.perm[k] = g_perm[tx.perm[k]] (NOT
  // tx.perm[g_perm[k]]). The two orders coincide for rank <= 2 (the external
  // groups are S_1/S_2, abelian and all-involution, so CCD/CCSD never exposed
  // the bug) but differ for the non-abelian S_3 (and higher) external groups of
  // triples and beyond — 3-cycles reconstruct to the wrong tensor with the
  // reversed order. (Validated offline at ranks 1-4 against a fully
  // antisymmetric + TRS reference: reversed order FAILs at rank 3, this order
  // PASSes to machine precision at every rank.)
  auto compose = [n](const Xform& tx, int g_sign, bool g_conj,
                     const Perm& g_perm) {
    Xform ty;
    ty.sign = tx.sign * g_sign;
    ty.conj = tx.conj ^ g_conj;
    ty.perm.resize(n);
    for (std::size_t k = 0; k < n; ++k) ty.perm[k] = g_perm[tx.perm[k]];
    return ty;
  };

  container::svector<KramersBlock> blocks;
  container::set<std::uint64_t> seen;
  for (std::uint64_t m = 0; m <= full_mask; ++m) {
    if (seen.find(m) != seen.end()) continue;

    // BFS the orbit, tracking each member's transform from the canonical m.
    container::map<std::uint64_t, Xform> xform;
    xform.emplace(m, Xform{+1, false, identity});
    container::svector<std::uint64_t> frontier{m};
    while (!frontier.empty()) {
      container::svector<std::uint64_t> next;
      for (const std::uint64_t x : frontier) {
        const Xform tx =
            xform.at(x);  // copy: emplace below may invalidate refs
        // antisymmetric generators (residual external antisymmetry): sign -1,
        // no conj, permutation = perm
        for (const auto& perm : antisym_perms) {
          const std::uint64_t y = apply_bit_perm(perm, x, n);
          if (xform.find(y) == xform.end()) {
            xform.emplace(y, compose(tx, -1, false, perm));
            next.push_back(y);
          }
        }
        // symmetric generators (e.g. particle interchange sigma on a raw,
        // non-antisymmetrized integral): sign +1, no conj, permutation = perm
        for (const auto& perm : symm_perms) {
          const std::uint64_t y = apply_bit_perm(perm, x, n);
          if (xform.find(y) == xform.end()) {
            xform.emplace(y, compose(tx, +1, false, perm));
            next.push_back(y);
          }
        }
        // global time reversal: complement, conj, sign (-1)^(#down of x)
        if (use_T) {
          const std::uint64_t y = (~x) & full_mask;
          if (xform.find(y) == xform.end()) {
            const int t_sign = (std::popcount(x) & 1) ? -1 : +1;
            xform.emplace(y, compose(tx, t_sign, true, identity));
            next.push_back(y);
          }
        }
      }
      frontier = std::move(next);
    }

    KramersBlock block;
    block.canonical = m;  // orbit-min (m increases, first unseen is the min)
    for (const auto& [cfg, t] : xform) {  // map iterates ascending by cfg
      block.members.push_back({cfg, t.sign, t.conj, t.perm});
      seen.insert(cfg);
    }
    blocks.push_back(std::move(block));
  }
  return blocks;
}

container::svector<container::svector<std::size_t>> kramers_external_groups(
    std::size_t rank) {
  container::svector<container::svector<std::size_t>> groups(2);
  for (std::size_t s = 0; s < rank; ++s) groups[0].push_back(s);
  for (std::size_t s = rank; s < 2 * rank; ++s) groups[1].push_back(s);
  return groups;
}

std::optional<KramersBlockMember> kramers_transform(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& groups,
    bool use_T, std::uint64_t canon, std::uint64_t cfg) {
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;
  container::svector<bool> in_group(n, false);
  for (auto const& g : groups)
    for (auto s : g) in_group[s] = true;

  // Order-preserving within-group matching with dst_bit[q] = src_bit[perm[q]]
  // (the KramersBlockMember convention: output slot q reads input slot
  // perm[q]); sign = product of the per-group permutation parities. nullopt
  // when a group's bit multiset differs or an ungrouped slot mismatches.
  auto match = [&](std::uint64_t src, std::uint64_t dst)
      -> std::optional<std::pair<int, container::svector<std::size_t>>> {
    for (std::size_t q = 0; q < n; ++q)
      if (!in_group[q] && (((src >> q) ^ (dst >> q)) & 1u)) return std::nullopt;
    container::svector<std::size_t> perm(n);
    std::iota(perm.begin(), perm.end(), std::size_t{0});
    int sign = +1;
    for (auto const& g : groups) {
      container::svector<std::size_t> src1, src0, dst1, dst0;
      for (auto s : g) (((src >> s) & 1u) ? src1 : src0).push_back(s);
      for (auto s : g) (((dst >> s) & 1u) ? dst1 : dst0).push_back(s);
      if (src1.size() != dst1.size()) return std::nullopt;
      // forward slot map (the KramersBlockMember convention validated by the
      // [orbit-transform] value test): input slot p feeds output slot perm[p]
      for (std::size_t k = 0; k < dst1.size(); ++k) perm[src1[k]] = dst1[k];
      for (std::size_t k = 0; k < dst0.size(); ++k) perm[src0[k]] = dst0[k];
      // parity of perm restricted to g
      container::svector<std::size_t> pos(n, 0);
      for (std::size_t k = 0; k < g.size(); ++k) pos[g[k]] = k;
      container::svector<std::size_t> seq;
      for (auto s : g) seq.push_back(pos[perm[s]]);
      for (std::size_t a = 0; a < seq.size(); ++a)
        for (std::size_t b = a + 1; b < seq.size(); ++b)
          if (seq[a] > seq[b]) sign = -sign;
    }
    return std::make_pair(sign, std::move(perm));
  };

  if (auto direct = match(canon, cfg))
    return KramersBlockMember{cfg, direct->first, false,
                              std::move(direct->second)};
  if (use_T) {
    // cfg = T(pre): block(cfg) = (-1)^(#down of pre) * conj(block(pre)); the
    // canon -> pre permutation carries through unchanged (T is slot-identity).
    const std::uint64_t pre = (~cfg) & full_mask;
    if (auto m = match(canon, pre)) {
      const int t_sign = (std::popcount(pre) & 1) ? -1 : +1;
      return KramersBlockMember{cfg, m->first * t_sign, true,
                                std::move(m->second)};
    }
  }
  return std::nullopt;
}

bool has_antisymmetrizer(const ExprPtr& expr) {
  return find_antisymmetrizer(expr).has_value();
}

namespace {

/// true if @p t is a Fock tensor whose indices span BOTH Kramers labels
bool is_mixed_kramers_fock(const Tensor& t) {
  std::wstring label{t.label()};
  if (has_conj_suffix(label)) label.pop_back();
  if (label != L"f") return false;
  bool up = false, dn = false;
  for (const auto& ix : t.const_braket()) {
    const auto s = to_spin(ix.space().qns());
    if (s == Spin::alpha)
      up = true;
    else if (s == Spin::beta)
      dn = true;
  }
  return up && dn;
}

}  // namespace

ExprPtr drop_mixed_kramers_fock_terms(const ExprPtr& expr) {
  if (!expr) return expr;
  if (expr->is<Tensor>())
    return is_mixed_kramers_fock(expr->as<Tensor>()) ? ex<Constant>(0) : expr;
  if (expr->is<Product>()) {
    const auto& p = expr->as<Product>();
    auto result = std::make_shared<Product>();
    result->scale(p.scalar());
    for (const auto& factor : p) {
      auto f = drop_mixed_kramers_fock_terms(factor);
      // one vanishing factor kills the whole product
      if (f->is<Constant>() && f->as<Constant>().is_zero())
        return ex<Constant>(0);
      result->append(1, f);
    }
    return result;
  }
  if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (const auto& summand : expr->as<Sum>().summands()) {
      auto s = drop_mixed_kramers_fock_terms(summand);
      if (s->is<Constant>() && s->as<Constant>().is_zero()) continue;
      result->append(s);
    }
    if (result->summands().empty()) return ex<Constant>(0);
    return result;
  }
  // Constant / Variable / opaque wrappers (e.g. RealPart) pass through
  return expr;
}

container::svector<ExprPtr> closed_shell_kramers_CC_trace(
    const ExprPtr& expr, bool expand_g, bool use_T,
    bool drop_mixed_kramers_fock) {
  // Stage 1: factor out the antisymmetrizer Â (kept, not expanded). Its bra/ket
  // are the external virtual/occupied index groups.
  auto A = find_antisymmetrizer(expr);
  SEQUANT_ASSERT(A.has_value() &&
                 "closed_shell_kramers_CC_trace: residual must carry Â");

  // External indices in bit order: bra group (virtuals) then ket group
  // (occupieds). Bit k = Kramers label of ext[k] (0 = up, 1 = down).
  container::svector<Index> ext;
  for (const auto& idx : A->bra()) ext.push_back(idx);
  for (const auto& idx : A->ket()) ext.push_back(idx);
  const std::size_t n_bra = A->bra_rank();
  const std::size_t n = ext.size();

  // Stage 2 (external fold): generators = the S_k adjacent transpositions of
  // each external group (rank-general: the residual's antisymmetry in each
  // external pair/tuple) + global time reversal T. Doubles -> 5 blocks.
  SEQUANT_ASSERT(n == 2 * n_bra &&
                 "closed_shell_kramers_CC_trace: expected equal-rank external "
                 "virtual/occupied groups");
  const auto gens = kramers_external_generators(n_bra);
  const auto orbits = kramers_config_orbits(n, gens, use_T);

  // A-expand: expand the antisymmetrizer Â into its explicit signed external
  // index permutations (the round-1 validated approach). This makes the FULL
  // external antisymmetry — INCLUDING cross-Kramers pairs — explicit in every
  // term. The earlier strip-Â + within-block bitwise antisymmetrization could
  // not supply this for the t-dependent terms: bitwise only antisymmetrizes
  // SAME-Kramers external pairs, but the cross-Kramers external antisymmetry of
  // a t-term needs its signed partner permutation, which only A-expand provides
  // (a post-hoc block-global sign merely flips, it does not antisymmetrize).
  // With ḡ still antisymmetric here, the driver Â[ḡ] reduces to ḡ: the bra!ket!
  // signed permutations of an antisymmetric tensor sum to bra!ket!·ḡ,
  // cancelling expand_A_op's 1/(bra!ket!) normalization.
  const bool kram_stage = std::getenv("MPQC_KRAM_PROTO_CHECK") != nullptr;
  auto stage = [&](const char* s) {
    if (kram_stage) std::wcout << L"[kram-stage] " << s << L"\n" << std::flush;
  };
  ExprPtr inner = expand_A_op(expr);
  stage("A-expand done");
  expand(inner);
  stage("expand done");
  rapid_simplify(inner);
  stage("rapid_simplify(inner) done -> entering block loop");

  // Optional g-expansion (AFTER A-expand, while still Kramers-FREE so
  // expand_antisymm's Ms-conserving guard keeps every permutation): replace
  // each ḡ with its raw NonSymm expansion. The factory [as] block omits the
  // cross- Kramers swap for any mixed-Kramers pair (INTERNAL pairs included),
  // so the [as] form is wrong there; raw-g leaves are each a correct direct
  // ⟨..|..⟩ fetch and the antisymmetrization becomes explicit config terms.
  if (expand_g) {
    inner = expand_label_antisymm(inner, L"g");
  }

  // Emit one block per external representative.
  container::svector<ExprPtr> blocks;
  blocks.reserve(orbits.size());

  // Rebuild each proto-indexed replacement value's proto-bundle by mapping the
  // ORIGINAL proto-indices through the same config map. make_spin{alpha,beta}
  // labels an index AND (via make_index_with_spincase's recursion) its proto-
  // indices with the *outer* index's single Kramers spin — which is wrong for a
  // CSV/PNS virtual whose proto-bundle is the occ pair (each occ carries its
  // own config spin, e.g. i↓ vs i↑). Left unfixed, canonicalize reconciles the
  // inconsistent proto vs occ labels by collapsing the two proto-indices to the
  // same dummy (a↑^{i↑ i↑}), which trips TensorNetworkV3::Edge::add_vertex. A
  // proto-index absent from the map (e.g. an already-labeled external occ
  // during the internal fold) is kept as-is.
  auto fix_proto_bundles = [](container::map<Index, Index>& repl) {
    for (auto& kv : repl) {
      const Index& orig = kv.first;
      if (!orig.has_proto_indices()) continue;
      auto proto = orig.proto_indices();
      for (auto& p : proto) {
        auto it = repl.find(p);
        if (it != repl.end() && !(it->second == p)) p = it->second;
      }
      // ALWAYS rebuild kv.second's proto-bundle from orig's (config-mapped)
      // proto, never keep make_spin{alpha,beta}'s: that helper stamps the WHOLE
      // proto-bundle with the outer virtual's single Kramers spin, which is
      // wrong for a CSV/PNS virtual whose proto is the occ pair. The correct
      // proto is orig's, with each proto-index carrying its OWN config spin
      // (mapped via repl if that occ is in this config map, else kept at its
      // already-labeled external value). Rebuilding only when some proto-index
      // was in the map (the old `if (changed)`) missed the case where EVERY
      // proto-index is external -- e.g. a DOWN internal virtual of an UP
      // external pair in the vv-Fock coupling f_a^c t̄: make_spinbeta stamped
      // i↓,i↓ but the pair (and thus the correct proto) is i↑,i↑. Left unfixed,
      // that virtual carries the down-pair proto tuple, so the amplitude spans
      // two distinct proto tuples and its ToT gains a spurious second pair of
      // outer occ axes that no longer matches the residual head.
      kv.second =
          Index(kv.second.space(), kv.second.ordinal(), std::move(proto));
    }
  };
  std::size_t block_idx = 0;
  for (const auto& orbit : orbits) {
    const std::uint64_t cfg = orbit.front();  // canonical (orbit-min)
    stage((std::string("block ") + std::to_string(block_idx++) +
           " start (cfg=" + std::to_string(cfg) + ")")
              .c_str());

    // Stage 2a: label the external indices for this block.
    container::map<Index, Index> ext_repl;
    for (std::size_t k = 0; k < n; ++k) {
      const bool down = (cfg >> k) & 1u;
      ext_repl.emplace(ext[k],
                       down ? make_spinbeta(ext[k]) : make_spinalpha(ext[k]));
    }
    fix_proto_bundles(ext_repl);
    ExprPtr block = append_spin(inner, ext_repl);

    // Stage 2b (internal fold): the still-unlabeled (spin-any) indices of a
    // term are its contracted/internal ones. Enumerate all their Kramers
    // configurations and sum (no Ms filter — relativistic); canonicalize then
    // performs the particle-interchange (sigma) merge. The internal-T-reach
    // (extra fold in the self-complementary external block) is deferred.
    //
    // The fold MUST be applied per term (summand): different terms of one
    // residual carry different internal-index sets, and a term invariant under
    // a given internal index (e.g. the driver, with no internal index) would be
    // over-counted 2^k if the whole block were folded over the union of all
    // terms' internal indices.
    auto fold_term = [&](const ExprPtr& term) -> ExprPtr {
      const auto used =
          get_used_indices<container::set<Index, Index::LabelCompare>>(term);
      container::svector<Index> internal;
      for (const auto& idx : used)
        if (to_spin(idx.space().qns()) == Spin::any) internal.push_back(idx);
      if (internal.empty()) return term;
      const std::size_t ni = internal.size();
      SEQUANT_ASSERT(ni <= 62);
      const std::uint64_t nconfigs = pow2(ni);
      auto sum = std::make_shared<Sum>();
      for (std::uint64_t ic = 0; ic < nconfigs; ++ic) {
        container::map<Index, Index> int_repl;
        for (std::size_t k = 0; k < ni; ++k) {
          const bool down = (ic >> k) & 1u;
          int_repl.emplace(internal[k], down ? make_spinbeta(internal[k])
                                             : make_spinalpha(internal[k]));
        }
        fix_proto_bundles(int_repl);
        sum->append(append_spin(term, int_repl));
      }
      return sum;
    };

    if (block->is<Sum>()) {
      auto folded = std::make_shared<Sum>();
      for (const auto& summand : *block) folded->append(fold_term(summand));
      block = folded;
    } else {
      block = fold_term(block);
    }

    // Diagnostic (MPQC_KRAM_PROTO_CHECK): detect any index whose proto-index
    // bundle contains a repeated index (e.g. a↑_1^{i↑_2 i↑_2}), which trips
    // TensorNetworkV3 during canonicalize. Reports whether the repeat is
    // present BEFORE canonicalize (from append_spin/fold) or introduced BY
    // canonicalize.
    auto proto_repeat_report = [](const ExprPtr& e, const char* when) {
      if (!std::getenv("MPQC_KRAM_PROTO_CHECK")) return;
      // String-based scan of the latex (robust vs LabelCompare dedup): look for
      // any tensor index whose proto-bundle repeats an index, i.e. `^{{X}{X}}`.
      const std::wstring tex = to_latex(e);
      static const std::wregex re(LR"(\^\{\{([^}]+)\}\{([^}]+)\}\})");
      for (auto it = std::wsregex_iterator(tex.begin(), tex.end(), re);
           it != std::wsregex_iterator(); ++it)
        if ((*it)[1].str() == (*it)[2].str())
          std::wcout << L"[proto-check " << when << L"] repeated proto bundle {"
                     << (*it)[1].str() << L"}{" << (*it)[2].str() << L"}\n";
    };
    proto_repeat_report(block, "pre-canon");
    if (std::getenv("MPQC_KRAM_BLOCK_DUMP"))
      std::wcout << L"[kram-block] block " << (block_idx - 1)
                 << L" pre-canon:\n"
                 << to_latex(block) << L"\n"
                 << std::flush;
    stage("  pre-canon check done, calling canonicalize");
    canonicalize(block);  // sigma merge + dummy canonicalization
    stage("  canonicalize done");
    proto_repeat_report(block, "post-canon");
    rapid_simplify(block);
    if (drop_mixed_kramers_fock) {
      auto filtered = drop_mixed_kramers_fock_terms(block);
      // Never annihilate a whole residual block: downstream consumers derive
      // the block's external (bra/ket) layout from its expression, and a bare
      // Constant{0} carries none. Such a block evaluates to zero anyway — at
      // zero cost, since its vanishing Fock leaf is a structurally empty
      // tensor — so keep it verbatim.
      if (!(filtered->is<Constant>() && filtered->as<Constant>().is_zero()))
        block = std::move(filtered);
    }
    blocks.push_back(block);
  }
  return blocks;
}

ExprPtr closed_shell_kramers_trace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool fold_T, bool expand_g, bool drop_mixed_kramers_fock) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Step 0/6: optionally expand the integral `g`'s antisymmetry (Kramers-free,
  // so no Ms filtering), keeping the amplitude `t` antisymmetric (the level-1
  // form). When expand_g is false, `g` stays antisymmetric (ḡ) so the evaluator
  // fetches the factory [as] block directly — the cross-Kramers antisymmetry is
  // then handled inside the integral, not by a raw config sum. (Expanding g is
  // a later optimization stage.)
  ExprPtr traced_input =
      expand_g ? expand_label_antisymm(expr, L"g") : expr->clone();
  canonicalize(traced_input);
  rapid_simplify(traced_input);

  // A SUM must be traced per summand: the configuration enumeration below
  // runs over the indices of the traced expression, so a summand that does
  // not use every index of the sum (a mixed-rank sum, e.g. the CCSD energy's
  // 2-slot f.t1 next to the 4-slot g.t2) would be replicated once per
  // configuration of the indices it does NOT carry -- 2^(n_missing) phantom
  // duplicates (measured: the CCSD energy's f.t1 config terms appeared 4x).
  // Tracing each summand over its own index set has no phantom bits; for
  // uniform-rank sums (MP2/CCD energies) the result is term-for-term
  // equivalent to the whole-sum enumeration.
  if (traced_input->is<Sum>()) {
    auto out = std::make_shared<Sum>();
    for (const auto& summand : traced_input->as<Sum>().summands())
      out->append(closed_shell_kramers_trace(summand, ext_index_groups, fold_T,
                                             /*expand_g=*/false));
    ExprPtr result{out};
    flatten(result);
    return result;
  }

  // Steps 1-2: classify indices and build groups (each internal index its own
  // group; external groups appended verbatim). No Ms / rank assumptions.
  const auto all_indices =
      get_used_indices<container::set<Index, Index::LabelCompare>>(
          traced_input);

  container::set<Index> ext_indices;
  for (const auto& group : ext_index_groups) {
    for (const auto& idx : group) {
      Index x = idx;
      x.reset_tag();
      ext_indices.insert(std::move(x));
    }
  }

  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> index_groups;
  for (const auto& idx : all_indices) {
    if (ext_indices.find(idx) == ext_indices.end())
      index_groups.emplace_back(IndexGroup(1, idx));
  }
  for (const auto& group : ext_index_groups) {
    index_groups.emplace_back(IndexGroup(group.begin(), group.end()));
  }

  // Step 3: enumerate the 2^n Kramers configurations and fold conjugate
  // (global time-reversal) partners. The config integer assigns bit k to
  // index_groups[k]; bit 0 = Kramers-up (Spin::alpha), bit 1 = down
  // (Spin::beta). The global-T partner is the full one's-complement; for a
  // closed contraction block(comp) = conj(block(cfg)), so the pair sums to
  // 2 Re(block). We keep the lexicographically smaller of each T-pair.
  const std::size_t n = index_groups.size();
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t nconfigs = pow2(n);
  const std::uint64_t full_mask = nconfigs - 1;

  // Accumulate canonical representatives, merging configurations that
  // `canonicalize` made identical (the particle-interchange / sigma fold) and
  // summing their T-fold multiplicities. RealPart is opaque to
  // `simplify`/`canonicalize`, so we merge explicitly here (exact Expr equality
  // of the canonicalized blocks) rather than relying on the Sum machinery.
  container::svector<std::pair<ExprPtr, std::int64_t>> reps;
  for (std::uint64_t cfg = 0; cfg < nconfigs; ++cfg) {
    if (fold_T) {
      const std::uint64_t comp = cfg ^ full_mask;
      if (cfg > comp) continue;  // keep canonical T-rep (cfg < comp for n >= 1)
    }

    container::map<Index, Index> replacements;
    for (std::size_t k = 0; k < n; ++k) {
      const bool down = (cfg >> k) & 1u;
      for (const auto& idx : index_groups[k]) {
        replacements.emplace(idx,
                             down ? make_spinbeta(idx) : make_spinalpha(idx));
      }
    }

    ExprPtr block = append_spin(traced_input, replacements);
    canonicalize(block);  // particle-interchange (sigma) merge within the block
    rapid_simplify(block);
    if (drop_mixed_kramers_fock) {
      block = drop_mixed_kramers_fock_terms(block);
      // a configuration left with nothing contributes no representative
      if (block->is<Constant>() && block->as<Constant>().is_zero()) continue;
    }

    bool merged = false;
    for (auto& [rep, mult] : reps) {
      if (*rep == *block) {
        ++mult;
        merged = true;
        break;
      }
    }
    if (!merged) reps.emplace_back(block, std::int64_t{1});
  }

  // Emit one representative per orbit. With the T-fold, coefficient = 2 x
  // multiplicity and each block is RealPart-wrapped (A + A* = 2 Re A). Without
  // it, every configuration is emitted verbatim (coefficient = multiplicity, no
  // RealPart) — the complex sum whose real part the caller takes.
  auto result = std::make_shared<Sum>();
  for (const auto& [block, mult] : reps) {
    if (fold_T)
      result->append(ex<Constant>(2 * mult) * real_part(block));
    else
      result->append(ex<Constant>(mult) * block);
  }
  return result;
}

namespace {

using sequant::mbpt::Spin;

// Kramers flavor of an index; Spin::any = unflavored.
Spin kr_flavor(const Index& idx) {
  return sequant::mbpt::to_spin(idx.space().qns());
}

// Flip the flavor of a PLAIN (proto-free) index; identity on unflavored.
Index kr_flip_plain(const Index& idx) {
  SEQUANT_ASSERT(!idx.has_proto_indices());
  const Spin s = kr_flavor(idx);
  if (s == Spin::any) return idx;
  return s == Spin::alpha ? make_spinbeta(make_spinfree(idx))
                          : make_spinalpha(make_spinfree(idx));
}

// Rebase one flat term (Product of Tensors/Constants/Variables, or a lone
// Tensor). Returns the rebased term; nullptr if the term has a shape this
// pass does not handle (caller keeps the original).
ExprPtr kramers_rebase_term(const ExprPtr& term,
                            const container::set<Index>& externals) {
  // --- collect the tensor factors
  container::svector<Tensor> leaves;
  ExprPtr scalar_tail;  // non-tensor factors (Constants/Variables), in order
  container::svector<std::size_t> leaf_pos;  // factor position of each leaf
  container::svector<ExprPtr> factors;
  if (term->is<Tensor>()) {
    factors.push_back(term);
  } else if (term->is<Product>()) {
    for (auto const& f : term->as<Product>().factors()) factors.push_back(f);
  } else {
    return nullptr;
  }
  for (std::size_t i = 0; i < factors.size(); ++i) {
    if (factors[i]->is<Tensor>()) {
      leaves.push_back(factors[i]->as<Tensor>());
      leaf_pos.push_back(i);
    } else if (factors[i]->is<Constant>() || factors[i]->is<Variable>()) {
      continue;
    } else {
      return nullptr;  // nested Sum etc.: leave the term untouched
    }
  }
  if (leaves.empty()) return nullptr;

  // --- flavored slot indices per leaf, and index -> leaves incidence
  auto slot_indices = [](const Tensor& t) {
    container::svector<Index> out;
    for (auto const& idx : t.bra()) out.push_back(idx);
    for (auto const& idx : t.ket()) out.push_back(idx);
    for (auto const& idx : t.aux()) out.push_back(idx);
    return out;
  };
  container::map<Index, container::svector<std::size_t>> incidence;
  for (std::size_t l = 0; l < leaves.size(); ++l)
    for (auto const& idx : slot_indices(leaves[l]))
      if (kr_flavor(idx) != Spin::any) incidence[idx].push_back(l);

  // --- union-find over leaves
  container::svector<std::size_t> parent(leaves.size());
  for (std::size_t l = 0; l < leaves.size(); ++l) parent[l] = l;
  auto find = [&parent](std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  for (auto const& [idx, ls] : incidence)
    for (std::size_t k = 1; k < ls.size(); ++k) {
      const auto a = find(ls[0]), b = find(ls[k]);
      if (a != b) parent[a] = b;
    }

  // --- freeze components touching externals or dangling flavored indices
  container::set<std::size_t> frozen;
  for (auto const& [idx, ls] : incidence) {
    const bool ext = externals.find(idx) != externals.end();
    const bool dangling = ls.size() == 1;
    if (ext || dangling)
      for (auto l : ls) frozen.insert(find(l));
  }

  // --- per component: decide the canonical orientation
  // fingerprint entry: (label, flavor bits over the leaf's flavored slots)
  auto leaf_print = [&](const Tensor& t, bool flipped) {
    std::pair<std::wstring, container::svector<int>> p{std::wstring{t.label()},
                                                       {}};
    for (auto const& idx : slot_indices(t)) {
      const Spin s = kr_flavor(idx);
      if (s == Spin::any) continue;
      const int down = (s == Spin::beta) ? 1 : 0;
      p.second.push_back(flipped ? 1 - down : down);
    }
    return p;
  };
  container::set<std::size_t> flip_roots;
  container::map<std::size_t, container::svector<std::size_t>> comp_leaves;
  for (std::size_t l = 0; l < leaves.size(); ++l)
    comp_leaves[find(l)].push_back(l);
  for (auto const& [root, ls] : comp_leaves) {
    if (frozen.find(root) != frozen.end()) continue;
    std::vector<std::pair<std::wstring, container::svector<int>>> fp, fp_flip;
    for (auto l : ls) {
      fp.push_back(leaf_print(leaves[l], false));
      fp_flip.push_back(leaf_print(leaves[l], true));
    }
    std::sort(fp.begin(), fp.end());
    std::sort(fp_flip.begin(), fp_flip.end());
    if (fp_flip < fp) flip_roots.insert(root);
  }
  if (flip_roots.empty()) return nullptr;  // nothing to do

  // --- build the replacement map: every flavored slot index of a flipped
  // component maps to its flavor-flipped spelling. Two passes: plain indices
  // first, then composites (base flavor flipped, protos mapped through the
  // plain entries so decoration follows its referent).
  container::map<Index, Index> repl;
  for (auto const& [idx, ls] : incidence) {
    if (flip_roots.find(find(ls.front())) == flip_roots.end()) continue;
    if (!idx.has_proto_indices()) repl.emplace(idx, kr_flip_plain(idx));
  }
  for (auto const& [idx, ls] : incidence) {
    if (flip_roots.find(find(ls.front())) == flip_roots.end()) continue;
    if (!idx.has_proto_indices()) continue;
    Index base(idx.space(), idx.ordinal());
    Index base_flipped = kr_flip_plain(base);
    auto protos = idx.proto_indices();
    for (auto& p : protos) {
      auto it = repl.find(p);
      if (it != repl.end()) p = it->second;
    }
    repl.emplace(idx,
                 Index(base_flipped.space(), base_flipped.ordinal(), protos));
  }

  // --- apply: transform every factor's indices (decoration everywhere
  // follows flipped referents); conjugate-mark the flipped components' leaves
  auto result = std::make_shared<Product>();
  if (term->is<Product>()) result->scale(term->as<Product>().scalar());
  std::size_t li = 0;
  for (std::size_t i = 0; i < factors.size(); ++i) {
    if (li < leaf_pos.size() && leaf_pos[li] == i) {
      Tensor t{leaves[li]};
      t.transform_indices(repl);
      if (flip_roots.find(find(li)) != flip_roots.end()) t.conjugate();
      result->append(1, ex<Tensor>(std::move(t)), Product::Flatten::No);
      ++li;
    } else {
      result->append(1, factors[i], Product::Flatten::No);
    }
  }
  return result;
}

}  // namespace

ExprPtr kramers_internal_rebase(const ExprPtr& expr,
                                const container::set<Index>& externals) {
  if (!expr) return expr;
  if (expr->is<Sum>()) {
    auto out = std::make_shared<Sum>();
    for (auto const& s : expr->as<Sum>().summands()) {
      auto r = kramers_rebase_term(s, externals);
      out->append(r ? r : s);
    }
    return out;
  }
  auto r = kramers_rebase_term(expr, externals);
  return r ? r : expr;
}

}  // namespace sequant::mbpt
