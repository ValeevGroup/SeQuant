//
// Spinor decomposition passes for relativistic 2-component theories.
// See spinor.hpp for design overview.
//

#include <SeQuant/domain/mbpt/spinor.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/core/tensor_network/v1.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace sequant::mbpt {

namespace {

/// Construct a new Index in the same IndexSpace type as @p idx, with
/// the @p k Kramers QN bit set (replacing any existing Kramers state)
/// and the label decorated with the matching ⇑/⇓ glyph.
///
/// Mirrors `detail::make_index_with_spincase` in `spin.cpp` so that
/// SeQuant's IndexSpace registry is reused: if a space with the
/// resulting label and qns already exists, we adopt it.
Index make_index_with_kramers(const Index& idx, mbpt::Kramers k) {
  // sanity check: at most one Kramers annotation already
  SEQUANT_ASSERT(!(idx.label().find(L'⇑') != std::wstring::npos &&
                   idx.label().find(L'⇓') != std::wstring::npos));

  auto qns = mbpt::kramers_annotation_remove(idx.space().qns()).unIon(k);
  IndexSpace space;
  const auto label =
      mbpt::kramers_annotation_replace(idx.space().base_key(), k);
  if (auto isr = get_default_context().index_space_registry()) {
    auto* space_ptr = isr->retrieve_ptr(label);
    if (space_ptr && space_ptr->type() == idx.space().type() &&
        space_ptr->qns() == qns) {
      space = *space_ptr;
    }
  }
  if (!space) {
    space = IndexSpace{label, idx.space().type(), qns,
                       idx.space().approximate_size()};
  }
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_kramers(pidx, k);
  return Index{space, idx.ordinal(), protoindices};
}

/// Standard SeQuant post-trace pipeline: reset index tags then
/// `simplify` (= expand + rapid_simplify + canonicalize +
/// rapid_simplify, see `expr_algorithms.cpp`). Mirrors the spin tracer
/// idiom at `spin.cpp:953-956`.
ExprPtr& simplify_with_reset(ExprPtr& e) {
  detail::reset_idx_tags(e);
  return simplify(e);
}

/// Read the specialized Kramers state (`up`/`down`/`any`) of an index.
Kramers index_kramers_state(const Index& idx) {
  const auto k_bits = idx.space().qns().to_int32() & mask_v<Kramers>;
  if (k_bits == static_cast<bitset_t>(Kramers::up)) return Kramers::up;
  if (k_bits == static_cast<bitset_t>(Kramers::down)) return Kramers::down;
  return Kramers::any;
}

/// Return a clone of @p t with `set_label(label)` applied. Used by the
/// label-rewriting passes (conjugation toggle, complex split,
/// antisymm-recombine) to avoid the noisier full-Tensor reconstruction.
ExprPtr tensor_with_label(const Tensor& t, std::wstring label) {
  auto out = std::static_pointer_cast<Tensor>(t.clone());
  out->set_label(std::move(label));
  return out;
}

/// Walk @p expr and collect the unique Kramers-eligible (`Kramers::any`)
/// indices in encounter order. Indices that already carry a specialized
/// Kramers state are skipped.
container::svector<Index> collect_kramers_indices(const ExprPtr& expr) {
  container::svector<Index> seen;
  container::set<Index> seen_set;
  expr->visit(
      [&](const ExprPtr& e) {
        if (!e->is<Tensor>()) return;
        for (const auto& idx : e->as<Tensor>().const_braket()) {
          if (index_kramers_state(idx) != Kramers::any) continue;
          if (seen_set.insert(idx).second) seen.push_back(idx);
        }
      },
      /* recursive = */ true);
  return seen;
}

/// Build the Index → Kramers-specialized-Index replacement for one
/// configuration (bit `k` of @p cfg → Kramers::down, else Kramers::up).
container::map<Index, Index> kramers_replacement_map(
    const container::svector<Index>& indices, std::uint64_t cfg) {
  container::map<Index, Index> repl;
  for (std::size_t k = 0; k < indices.size(); ++k) {
    const bool is_dn = (cfg >> k) & 1u;
    repl[indices[k]] =
        is_dn ? make_kramers_dn(indices[k]) : make_kramers_up(indices[k]);
  }
  return repl;
}

/// Substitute Kramers labels onto @p expr per @p repl, expand
/// antisymmetrizers, then run the canonical SeQuant simplify pass.
ExprPtr specialize_and_simplify(const ExprPtr& expr,
                                const container::map<Index, Index>& repl) {
  auto term = append_kramers(expr, repl);
  term = expand_antisymm(term, /* skip_spinsymm */ false);
  return simplify_with_reset(term);
}

}  // namespace

Index make_kramers_up(const Index& idx) {
  return make_index_with_kramers(idx, mbpt::Kramers::up);
}

Index make_kramers_dn(const Index& idx) {
  return make_index_with_kramers(idx, mbpt::Kramers::down);
}

Index make_kramers_free(const Index& idx) {
  return make_index_with_kramers(idx, mbpt::Kramers::any);
}

std::wstring toggle_conj_suffix(std::wstring_view label) {
  if (has_conj_suffix(label))
    return std::wstring{label.substr(0, label.size() - 1)};
  return std::wstring{label} + L'*';
}

ExprPtr append_kramers(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements) {
  // The transformation is identical to spin's `append_spin`: in both
  // cases the only operation is `Tensor::transform_indices`, which is
  // QN-agnostic. Delegate to avoid duplicating the recursion logic.
  return append_spin(expr, index_replacements);
}

// =====================================================================
// kramers_trace: 2^n exhaustive enumeration.
// =====================================================================

ExprPtr kramers_trace(const ExprPtr& expr) {
  // Enumerate 2^N Kramers configurations of @p expr's spinor
  // (Kramers::any) indices, then sum. TRS pairs (X, X̄) related by full
  // bar-flip are complex conjugates of each other under the (-1)^k
  // identity; per-tensor signs cancel iff every barred index appears in
  // an even number of tensors — true for any all-internal-index
  // contraction. CC residuals carry external bars and need explicit
  // (-1)^(external bars) bookkeeping; not yet implemented.

  detail::reset_idx_tags(expr);
  const auto indices = collect_kramers_indices(expr);
  const std::size_t n = indices.size();
  if (n == 0) return expr;
  if (n > 30) {
    throw std::runtime_error(
        "kramers_trace: too many Kramers indices (would emit > 2^30 terms)");
  }

  auto result = std::make_shared<Sum>();
  const std::uint64_t n_configs = std::uint64_t{1} << n;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    result->append(
        specialize_and_simplify(expr, kramers_replacement_map(indices, cfg)));
  }
  ExprPtr r = result;
  return simplify_with_reset(r);
}

container::svector<ResultExpr> kramers_trace(const ResultExpr& expr,
                                             bool fold_klein) {
  // Collect external indices from LHS (bra + ket; aux ignored — the
  // pair-density use case doesn't have aux slots).
  container::svector<Index> ext_indices;
  for (const auto& idx : expr.bra()) ext_indices.push_back(idx);
  for (const auto& idx : expr.ket()) ext_indices.push_back(idx);

  const std::size_t n_ext = ext_indices.size();
  if (n_ext > 30) {
    throw std::runtime_error(
        "kramers_trace(ResultExpr): too many external indices");
  }
  if (n_ext == 0) {
    // Scalar case: just run the regular tracer + pipeline and return
    // a single ResultExpr (variable-shape; no LHS slots).
    auto traced = kramers_trace(expr.expression());
    traced = fold_conj_pairs(traced);
    traced = antisymm_recombine(traced);
    container::svector<ResultExpr> single;
    single.emplace_back(
        bra(container::svector<Index>{}), ket(container::svector<Index>{}),
        aux(container::svector<Index>{}), expr.symmetry(),
        expr.braket_symmetry(), expr.column_symmetry(),
        expr.has_label() ? std::optional<std::wstring>(expr.label())
                         : std::nullopt,
        traced);
    return single;
  }

  // Process only canonical representatives under the Kramers
  // symmetry group of the LHS. For a rank-(p,q) result the symmetry
  // group is the Klein 4-group generated by:
  //   α: flip all bra-side Kramers states  (partial flip on virtuals)
  //   β: flip all ket-side Kramers states  (partial flip on occupieds)
  //   αβ: flip all Kramers states           (full TRS)
  // For Kramers-restricted closed-shell systems these are all
  // symmetries of the trace (verifiable: traced expressions in one
  // orbit are equal up to Kramers-conjugation relabeling). The orbit
  // size is up to 4; emitting only the lex-smallest cfg in each orbit
  // reduces 2^N raw blocks by a factor of 4 (typically).
  const std::uint64_t n_configs = std::uint64_t{1} << n_ext;
  const std::uint64_t mask = n_configs - 1;
  const std::size_t n_bra = expr.bra().size();
  const std::uint64_t bra_mask = ((std::uint64_t{1} << n_bra) - 1);
  const std::uint64_t ket_mask = mask ^ bra_mask;

  container::svector<ResultExpr> out;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    const std::uint64_t cfg_alpha = cfg ^ bra_mask;
    const std::uint64_t cfg_beta = cfg ^ ket_mask;
    const std::uint64_t cfg_full = cfg ^ mask;
    const std::uint64_t canonical =
        fold_klein ? std::min({cfg, cfg_alpha, cfg_beta, cfg_full})
                   : std::min(cfg, cfg_full);
    if (cfg != canonical) continue;  // skip non-canonical

    // Build Kramers-specialized externals for this config.
    container::map<Index, Index> ext_repl;
    for (std::size_t k = 0; k < n_ext; ++k) {
      const bool is_dn = (cfg >> k) & 1u;
      ext_repl[ext_indices[k]] = is_dn ? make_kramers_dn(ext_indices[k])
                                       : make_kramers_up(ext_indices[k]);
    }

    // Substitute externals on the RHS, then run the FULL kramers
    // pipeline on the substituted expression:
    //   1) kramers_trace — enumerates remaining Kramers::any indices
    //      (the internal dummies) + expand_antisymm + canonicalize
    //      per config + final cross-product canonicalize. After this
    //      the RHS is a Sum of NonSymm products.
    //   2) antisymm_recombine — folds direct/exchange NonSymm pairs
    //      back into Kramers-restricted Antisymm tensors, restoring
    //      the mixed antisym/non-antisym MP2-style structure.
    //
    // Install a per-block `named_indices` context override carrying
    // the kramers-substituted externals — every external appears
    // exactly twice in the RHS (once via the LHS-side slot, once via
    // the dummy contraction), so without the override SeQuant's
    // canonicalize would treat them as dummies and rename them. The
    // override keeps the externals fixed and consistent with the
    // kramers-specialized LHS we build below.
    //
    // We do NOT call fold_conj_pairs here: it folds TRS pairs WITHIN
    // a Sum, but within a single external-fixed block the externals
    // are already specialized — flip_kramers would change them too
    // and the partner would land in a different block. External TRS
    // folding is handled at the block level above (cfg <= cfg_flipped
    // canonicalization).
    auto rhs_subs = append_kramers(expr.expression(), ext_repl);
    container::set<Index> specialized_externals;
    for (const auto& kv : ext_repl) specialized_externals.insert(kv.second);
    auto saved_ctx = get_default_context();
    auto block_ctx = saved_ctx;
    CanonicalizeOptions copts;
    copts.named_indices = specialized_externals;
    block_ctx.set(copts);
    set_default_context(block_ctx);
    auto traced = kramers_trace(rhs_subs);
    traced = antisymm_recombine(traced);
    set_default_context(saved_ctx);

    // Drop blocks that vanished.
    if (traced->is<Constant>() && traced->as<Constant>().is_zero()) continue;
    if (traced->is<Sum>() && traced->as<Sum>().summands().empty()) continue;

    // Build the kramers-specialized LHS slots in input order.
    container::svector<Index> spec_bra, spec_ket;
    for (const auto& idx : expr.bra()) spec_bra.push_back(ext_repl.at(idx));
    for (const auto& idx : expr.ket()) spec_ket.push_back(ext_repl.at(idx));

    out.emplace_back(bra(std::move(spec_bra)), ket(std::move(spec_ket)),
                     aux(container::svector<Index>{}), expr.symmetry(),
                     expr.braket_symmetry(), expr.column_symmetry(),
                     expr.has_label()
                         ? std::optional<std::wstring>(expr.label())
                         : std::nullopt,
                     traced);
  }
  return out;
}

// =====================================================================
// flip_kramers + fold_conj_pairs.
// =====================================================================

ExprPtr flip_kramers(const ExprPtr& expr) {
  // Walk every tensor, build the Index → flipped-Index map for any
  // index already carrying a specialized Kramers state, then apply via
  // append_kramers.
  container::map<Index, Index> repl;
  expr->visit(
      [&](const ExprPtr& e) {
        if (!e->is<Tensor>()) return;
        for (const auto& idx : e->as<Tensor>().const_braket()) {
          if (repl.contains(idx)) continue;
          switch (index_kramers_state(idx)) {
            case Kramers::up:
              repl[idx] = make_kramers_dn(idx);
              break;
            case Kramers::down:
              repl[idx] = make_kramers_up(idx);
              break;
            default:
              break;  // Kramers::any or unset — leave alone
          }
        }
      },
      /* recursive = */ true);
  return repl.empty() ? expr : append_kramers(expr, repl);
}

namespace {

/// Extract the leading scalar of a Product, or 1 for a non-Product.
Constant::scalar_type product_scalar(const ExprPtr& e) {
  return e->is<Product>() ? e->as<Product>().scalar()
                          : Constant::scalar_type{1};
}

/// Build a copy of @p e with the leading Product scalar dropped.
/// Returns a clone if @p e is not a Product.
ExprPtr product_factors_only(const ExprPtr& e) {
  if (!e->is<Product>()) return e->clone();
  auto p = std::make_shared<Product>();
  for (auto&& f : e->as<Product>()) p->append(1, f->clone());
  return p;
}

}  // namespace

ExprPtr fold_conj_pairs(const ExprPtr& expr) {
  detail::reset_idx_tags(expr);
  if (!expr->is<Sum>()) return expr;
  auto const& sum = expr->as<Sum>();
  const std::size_t n = sum.size();

  // O(n) canonical-form hash bucketing: for each summand A, the
  // canonical key is the lex-min of (normalize(A), normalize(flip(A))).
  // Both A and its full-bar partner map to the same bucket. Per-bucket
  // work is O(1) on average since TRS orbits have size ≤ 2.
  struct BucketEntry {
    std::size_t idx;
    ExprPtr factors_normalized;
    Constant::scalar_type scalar;
    bool self_conj;
  };
  container::map<Expr::hash_type, container::svector<BucketEntry>> buckets;

  std::vector<ExprPtr> originals;
  originals.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    originals.push_back(sum.summand(i)->clone());
    auto factors = product_factors_only(originals[i]);
    simplify_with_reset(factors);
    auto factors_flipped = flip_kramers(factors);
    simplify_with_reset(factors_flipped);
    const bool self_conj = (*factors == *factors_flipped);
    ExprPtr canon = (*factors < *factors_flipped) ? factors : factors_flipped;
    buckets[canon->hash_value()].push_back(
        {i, factors, product_scalar(originals[i]), self_conj});
  }

  auto out = std::make_shared<Sum>();
  std::vector<bool> used(n, false);

  // Emit the folded form of (A, B) per their scalar relationship.
  auto emit_pair = [&](const BucketEntry& a, const BucketEntry& b) {
    using cscalar = Constant::scalar_type;
    if (a.scalar == b.scalar) {
      out->append(ex<Constant>(rational{2}) * real_part(originals[a.idx]));
    } else if (a.scalar == -b.scalar) {
      out->append(ex<Constant>(cscalar{Complex<rational>{0, 2}}) *
                  imaginary_part(originals[a.idx]));
    } else {
      out->append(originals[a.idx]);
      out->append(originals[b.idx]);
    }
    used[a.idx] = used[b.idx] = true;
  };

  for (auto& [_, entries] : buckets) {
    std::vector<bool> taken(entries.size(), false);
    for (std::size_t i = 0; i < entries.size(); ++i) {
      if (taken[i] || used[entries[i].idx]) continue;
      const auto& a = entries[i];
      if (a.self_conj) {
        // Self-conjugate: pair with another self-conj if present, else
        // emit the lone Re(A).
        bool merged = false;
        for (std::size_t j = i + 1; j < entries.size(); ++j) {
          if (taken[j] || used[entries[j].idx]) continue;
          const auto& b = entries[j];
          if (b.self_conj && *b.factors_normalized == *a.factors_normalized) {
            emit_pair(a, b);
            taken[i] = taken[j] = true;
            merged = true;
            break;
          }
        }
        if (!merged) {
          out->append(real_part(originals[a.idx]));
          used[a.idx] = true;
          taken[i] = true;
        }
        continue;
      }
      // Non-self-conj: confirm partner via structural equality with
      // flip(A) (defending against hash collisions).
      bool paired = false;
      for (std::size_t j = i + 1; j < entries.size(); ++j) {
        if (taken[j] || used[entries[j].idx]) continue;
        const auto& b = entries[j];
        auto a_flip = flip_kramers(a.factors_normalized);
        simplify_with_reset(a_flip);
        if (*b.factors_normalized == *a_flip) {
          emit_pair(a, b);
          taken[i] = taken[j] = true;
          paired = true;
          break;
        }
      }
      if (!paired) {
        out->append(originals[a.idx]);
        used[a.idx] = true;
        taken[i] = true;
      }
    }
  }

  return out;
}

// =====================================================================
// Cycle decomposition (used by kramers_trace_cycles).
// =====================================================================

namespace {

/// One step of a contraction-cycle walk.
struct CycleNode {
  std::size_t tensor_idx;
  std::size_t braket_pos;
  Index idx;
};

/// One contraction cycle: alternating intra-tensor (bra[k] ↔ ket[k] of
/// the same tensor) and inter-tensor (same Index in two tensors) edges.
struct ContractionCycle {
  container::svector<CycleNode> nodes;
};

/// Decompose a Product's index-contraction graph into cycles.
///
/// @pre every Index in @p product appears in exactly two tensor slot
///      positions; each tensor's bra[k]/ket[k] positions pair as
///      intra-tensor edges
/// @todo CC support: the k ↔ k+rank/2 pairing assumes a rectangular
///       tensor (bra_rank == ket_rank); some CC intermediates aren't.
container::svector<ContractionCycle> decompose_cycles(const Product& product) {
  struct Slot {
    std::size_t tensor_idx;
    std::size_t braket_pos;
    Index idx;
  };
  container::svector<Slot> slots;
  std::size_t tensor_counter = 0;
  for (auto&& f : product) {
    if (!f->is<Tensor>()) {
      ++tensor_counter;
      continue;
    }
    auto const& t = f->as<Tensor>();
    std::size_t pos = 0;
    for (auto const& idx : t.const_braket())
      slots.push_back(Slot{tensor_counter, pos++, idx});
    ++tensor_counter;
  }
  const std::size_t N = slots.size();
  container::svector<bool> visited(N, false);
  container::svector<std::size_t> intra_partner(N, N);
  container::svector<std::size_t> inter_partner(N, N);

  // Intra-tensor pairing per tensor (k ↔ k+rank/2 within each
  // contiguous slot run, assuming rectangular).
  for (std::size_t i = 0; i < N;) {
    std::size_t j = i;
    while (j < N && slots[j].tensor_idx == slots[i].tensor_idx) ++j;
    const std::size_t tensor_rank = j - i;
    SEQUANT_ASSERT(tensor_rank % 2 == 0);
    const std::size_t half = tensor_rank / 2;
    for (std::size_t k = 0; k < half; ++k) {
      intra_partner[i + k] = i + half + k;
      intra_partner[i + half + k] = i + k;
    }
    i = j;
  }

  // Inter-tensor pairing: same Index in two distinct slots.
  for (std::size_t i = 0; i < N; ++i) {
    if (inter_partner[i] != N) continue;
    for (std::size_t j = i + 1; j < N; ++j) {
      if (inter_partner[j] == N && slots[j].idx == slots[i].idx) {
        inter_partner[i] = j;
        inter_partner[j] = i;
        break;
      }
    }
  }

  // Walk cycles: alternate intra → inter → intra → ... edges.
  container::svector<ContractionCycle> cycles;
  for (std::size_t start = 0; start < N; ++start) {
    if (visited[start]) continue;
    ContractionCycle cyc;
    std::size_t cur = start;
    bool take_intra = true;
    while (!visited[cur]) {
      visited[cur] = true;
      cyc.nodes.push_back(
          {slots[cur].tensor_idx, slots[cur].braket_pos, slots[cur].idx});
      const std::size_t next =
          take_intra ? intra_partner[cur] : inter_partner[cur];
      if (next == N) break;
      cur = next;
      take_intra = !take_intra;
    }
    cycles.push_back(std::move(cyc));
  }
  return cycles;
}

/// Distinct Kramers-eligible (`Kramers::any`) indices visited by a
/// cycle, in first-encounter order.
container::svector<Index> cycle_distinct_indices(const ContractionCycle& c) {
  container::svector<Index> out;
  container::set<Index> seen;
  for (auto const& n : c.nodes) {
    if (index_kramers_state(n.idx) != Kramers::any) continue;
    if (seen.insert(n.idx).second) out.push_back(n.idx);
  }
  return out;
}

/// Rotation-invariant shape signature of a cycle: the lex-min rotation
/// of the IndexSpace-type sequence along the walk. Two cycles with the
/// same shape are interchangeable under the expression's tensor
/// (anti)symmetries.
std::wstring cycle_shape(const ContractionCycle& c) {
  std::wstring raw;
  for (auto const& n : c.nodes)
    raw += static_cast<wchar_t>(L'a' + (n.idx.space().type().to_int32() & 0xF));
  if (raw.empty()) return raw;
  std::wstring best = raw, rot = raw;
  for (std::size_t k = 1; k < raw.size(); ++k) {
    std::rotate(rot.begin(), rot.begin() + 1, rot.end());
    if (rot < best) best = rot;
  }
  return best;
}

/// Multinomial coefficient (Σ counts)! / ∏ count_i!.
double multiset_multiplicity(const container::svector<std::size_t>& counts) {
  auto factorial = [](std::size_t k) {
    double f = 1;
    for (std::size_t i = 2; i <= k; ++i) f *= static_cast<double>(i);
    return f;
  };
  std::size_t total = 0;
  for (auto c : counts) total += c;
  double num = factorial(total);
  for (auto c : counts) num /= factorial(c);
  return num;
}

}  // namespace

ExprPtr kramers_trace_cycles(const ExprPtr& expr) {
  detail::reset_idx_tags(expr);
  if (expr->is<Sum>()) {
    auto sum_result = std::make_shared<Sum>();
    for (auto&& s : *expr) sum_result->append(kramers_trace_cycles(s));
    ExprPtr r = sum_result;
    return simplify_with_reset(r);
  }
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  Product prod;
  if (expr->is<Product>()) {
    prod = expr->as<Product>();
  } else if (expr->is<Tensor>()) {
    prod.append(1, expr, Product::Flatten::No);
  } else {
    throw std::runtime_error("kramers_trace_cycles: unsupported Expr type");
  }

  // 1. Decompose into cycles + collect per-cycle Kramers-eligible indices.
  const auto cycles = decompose_cycles(prod);
  if (cycles.empty()) return expr;
  const std::size_t n_cycles = cycles.size();
  std::vector<container::svector<Index>> cyc_indices(n_cycles);
  for (std::size_t c = 0; c < n_cycles; ++c)
    cyc_indices[c] = cycle_distinct_indices(cycles[c]);

  // 2. Group cycles into shape-equivalence classes.
  container::map<std::wstring, container::svector<std::size_t>> shape_classes;
  for (std::size_t c = 0; c < n_cycles; ++c)
    shape_classes[cycle_shape(cycles[c])].push_back(c);

  // 3. Per class: enumerate Kramers-labeling MULTISETS (sorted m-tuples
  //    of per-cycle 2^d configs) with their multinomial multiplicity.
  struct ClassConfig {
    container::svector<std::size_t> class_ids;
    container::svector<std::pair<container::svector<std::uint32_t>, double>>
        multisets;
  };
  container::svector<ClassConfig> class_configs;
  for (auto& [shape, cyc_ids] : shape_classes) {
    ClassConfig cc;
    cc.class_ids = cyc_ids;
    const std::size_t m = cyc_ids.size();
    const std::size_t L = std::size_t{1} << cyc_indices[cyc_ids.front()].size();
    container::svector<std::size_t> pick(m, 0);
    std::function<void(std::size_t, std::size_t)> rec = [&](std::size_t pos,
                                                            std::size_t start) {
      if (pos == m) {
        container::svector<std::uint32_t> assignment(m);
        for (std::size_t i = 0; i < m; ++i)
          assignment[i] = static_cast<std::uint32_t>(pick[i]);
        container::svector<std::size_t> counts;
        std::size_t run = 1;
        for (std::size_t i = 1; i <= m; ++i) {
          if (i < m && pick[i] == pick[i - 1]) {
            ++run;
          } else {
            counts.push_back(run);
            run = 1;
          }
        }
        cc.multisets.emplace_back(assignment, multiset_multiplicity(counts));
        return;
      }
      for (std::size_t v = start; v < L; ++v) {
        pick[pos] = v;
        rec(pos + 1, v);
      }
    };
    rec(0, 0);
    class_configs.push_back(std::move(cc));
  }

  // 4. Cartesian product across classes: each global config picks one
  //    multiset per class; build the per-index Kramers replacement,
  //    specialize, and scale by the product of multiplicities.
  auto result = std::make_shared<Sum>();
  const std::size_t n_classes = class_configs.size();
  container::svector<std::size_t> sel(n_classes, 0);
  std::function<void(std::size_t)> build = [&](std::size_t ci) {
    if (ci == n_classes) {
      container::map<Index, Index> repl;
      double mult = 1.0;
      for (std::size_t c = 0; c < n_classes; ++c) {
        const auto& cc = class_configs[c];
        const auto& [assignment, m] = cc.multisets[sel[c]];
        mult *= m;
        for (std::size_t k = 0; k < cc.class_ids.size(); ++k) {
          const std::size_t cyc = cc.class_ids[k];
          const std::uint32_t labeling = assignment[k];
          for (std::size_t d = 0; d < cyc_indices[cyc].size(); ++d) {
            const bool is_dn = (labeling >> d) & 1u;
            const Index& src = cyc_indices[cyc][d];
            repl[src] = is_dn ? make_kramers_dn(src) : make_kramers_up(src);
          }
        }
      }
      auto term = specialize_and_simplify(prod.clone(), repl);
      if (mult != 1.0) {
        term = ex<Constant>(rational{static_cast<std::intmax_t>(mult)}) * term;
        non_canon_simplify(term);
      }
      result->append(term);
      return;
    }
    for (std::size_t s = 0; s < class_configs[ci].multisets.size(); ++s) {
      sel[ci] = s;
      build(ci + 1);
    }
  };
  build(0);

  ExprPtr r = result;
  return simplify_with_reset(r);
}

// =====================================================================
// Antisymmetrizer recombination (inverse of expand_antisymm).
// =====================================================================

namespace {

/// Parsed Sum summand: total scalar, wrapper kind (0 raw / 1 Re / 2
/// Im), and the underlying tensor factors.
struct RecombEntry {
  Constant::scalar_type scalar{1};
  int wrapper = 0;
  container::svector<ExprPtr> tensors;
  bool alive = true;
};

/// Parse a Sum summand into a RecombEntry. Recognized shapes:
///   Tensor | Product(s; Tensor...) | Re/Im(inner)
///                                  | Product(s; Re/Im(inner))
/// where inner = Tensor | Product(s; Tensor...).
bool parse_recomb_entry(const ExprPtr& summand, RecombEntry& e) {
  e = RecombEntry{};
  ExprPtr cur = summand;
  if (cur->is<Product>() && cur->as<Product>().size() == 1) {
    auto const& p = cur->as<Product>();
    auto const& f0 = p.factor(0);
    if (f0->is<RealPart>() || f0->is<ImagPart>()) {
      e.scalar = e.scalar * p.scalar();
      cur = f0;
    }
  }
  if (cur->is<RealPart>()) {
    e.wrapper = 1;
    cur = cur->as<RealPart>().inner();
  } else if (cur->is<ImagPart>()) {
    e.wrapper = 2;
    cur = cur->as<ImagPart>().inner();
  }
  if (cur->is<Tensor>()) {
    e.tensors.push_back(cur);
    return true;
  }
  if (cur->is<Product>()) {
    auto const& p = cur->as<Product>();
    e.scalar = e.scalar * p.scalar();
    for (auto&& f : p) {
      if (!f->is<Tensor>()) return false;
      e.tensors.push_back(f);
    }
    return true;
  }
  return false;
}

/// Detect whether @p t2 is @p t1 with a single column transposition.
/// Returns 1 if t2 = bra-swap(t1), 2 if t2 = ket-swap(t1), 0 otherwise.
///
/// @todo CC support: rank-(2,2) only. CCSDT triples residuals need the
///       single-column-swap detection generalized to rank-(n,n).
int detect_single_column_swap(const Tensor& t1, const Tensor& t2) {
  if (t1.label() != t2.label()) return 0;
  if (t1.bra_rank() != 2 || t1.ket_rank() != 2) return 0;
  if (t2.bra_rank() != 2 || t2.ket_rank() != 2) return 0;
  if (t1.symmetry() != t2.symmetry()) return 0;
  const auto& b1 = t1.bra();
  const auto& k1 = t1.ket();
  const auto& b2 = t2.bra();
  const auto& k2 = t2.ket();
  const bool bra_same = (b1.at(0) == b2.at(0) && b1.at(1) == b2.at(1));
  const bool bra_swap = (b1.at(0) == b2.at(1) && b1.at(1) == b2.at(0));
  const bool ket_same = (k1.at(0) == k2.at(0) && k1.at(1) == k2.at(1));
  const bool ket_swap = (k1.at(0) == k2.at(1) && k1.at(1) == k2.at(0));
  if (bra_swap && ket_same) return 1;
  if (ket_swap && bra_same) return 2;
  return 0;
}

/// Re-tag a Tensor with `Symmetry::Antisymm` (Tensor::symmetry is
/// constructor-set; we rebuild via the existing Tensor ctor).
ExprPtr make_antisymm(const Tensor& t) {
  return ex<Tensor>(
      t.label(), bra(container::svector<Index>{t.bra().begin(), t.bra().end()}),
      ket(container::svector<Index>{t.ket().begin(), t.ket().end()}),
      aux(container::svector<Index>{t.aux().begin(), t.aux().end()}),
      Symmetry::Antisymm, t.braket_symmetry(), t.column_symmetry());
}

/// Try to recombine @p A and @p B into one antisymmetrized entry.
/// Succeeds when, for some tensor position p: all other positions
/// match structurally, t2[p] is a single column-swap of t1[p], and
/// `B.scalar == -A.scalar` (single transposition → odd → sign -1).
std::optional<RecombEntry> try_recombine(const RecombEntry& A,
                                         const RecombEntry& B) {
  if (A.wrapper != B.wrapper) return std::nullopt;
  if (A.tensors.size() != B.tensors.size()) return std::nullopt;
  if (!(B.scalar == -A.scalar)) return std::nullopt;
  const std::size_t n = A.tensors.size();
  for (std::size_t p = 0; p < n; ++p) {
    bool others_match = true;
    for (std::size_t q = 0; q < n; ++q) {
      if (q == p) continue;
      if (!(*A.tensors[q] == *B.tensors[q])) {
        others_match = false;
        break;
      }
    }
    if (!others_match) continue;
    if (!A.tensors[p]->is<Tensor>() || !B.tensors[p]->is<Tensor>()) continue;
    if (detect_single_column_swap(A.tensors[p]->as<Tensor>(),
                                  B.tensors[p]->as<Tensor>()) == 0)
      continue;
    RecombEntry R = A;
    R.tensors[p] = make_antisymm(A.tensors[p]->as<Tensor>());
    return R;
  }
  return std::nullopt;
}

}  // namespace

ExprPtr antisymm_recombine(const ExprPtr& expr) {
  detail::reset_idx_tags(expr);
  if (!expr->is<Sum>()) return expr;
  auto const& sum = expr->as<Sum>();

  std::vector<RecombEntry> entries;
  entries.reserve(sum.size());
  for (auto&& s : sum) {
    RecombEntry e;
    if (!parse_recomb_entry(s, e)) return expr;  // unrecognized — bail safely
    entries.push_back(std::move(e));
  }

  // Iterate to a fixed point: each pass merges one recombinable pair.
  // A second pass over once-recombined terms catches doubly-antisymm
  // blocks (the merged antisymm tensor becomes the "fixed" tensor).
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t i = 0; i < entries.size() && !changed; ++i) {
      if (!entries[i].alive) continue;
      for (std::size_t j = i + 1; j < entries.size() && !changed; ++j) {
        if (!entries[j].alive) continue;
        auto rec = try_recombine(entries[i], entries[j]);
        if (!rec) rec = try_recombine(entries[j], entries[i]);
        if (rec) {
          entries[i] = *rec;
          entries[j].alive = false;
          changed = true;
        }
      }
    }
  }

  auto out = std::make_shared<Sum>();
  for (auto const& e : entries) {
    if (!e.alive) continue;
    auto prod = std::make_shared<Product>();
    for (auto const& t : e.tensors) prod->append(1, t, Product::Flatten::No);
    ExprPtr inner = prod;
    if (e.wrapper == 1)
      inner = ex<RealPart>(inner);
    else if (e.wrapper == 2)
      inner = ex<ImagPart>(inner);
    out->append(ex<Constant>(e.scalar) * inner);
  }
  return out;
}

// =====================================================================
// Burnside-orbit Kramers tracer.
// =====================================================================

namespace {

/// True iff @p p and @p q both sit in the same antisymmetric bra/ket
/// group of @p t (i.e. t is Antisymm and {p,q} ⊆ bra or {p,q} ⊆ ket).
bool same_antisymm_group(const Tensor& t, const Index& p, const Index& q) {
  if (t.symmetry() != Symmetry::Antisymm) return false;
  auto in = [](const auto& range, const Index& x) {
    return std::find(range.begin(), range.end(), x) != range.end();
  };
  const bool p_bra = in(t.bra(), p), q_bra = in(t.bra(), q);
  const bool p_ket = in(t.ket(), p), q_ket = in(t.ket(), q);
  return (p_bra && q_bra) || (p_ket && q_ket);
}

bool tensor_has_index(const Tensor& t, const Index& p) {
  for (auto&& idx : t.const_braket())
    if (idx == p) return true;
  return false;
}

container::svector<const Tensor*> product_tensors(const Product& p) {
  container::svector<const Tensor*> out;
  for (auto&& f : p)
    if (f->is<Tensor>()) out.push_back(&f->as<Tensor>());
  return out;
}

}  // namespace

ExprPtr kramers_trace_burnside(const ExprPtr& expr) {
  detail::reset_idx_tags(expr);
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& t : *expr) s->append(kramers_trace_burnside(t));
    ExprPtr r = s;
    return simplify_with_reset(r);
  }
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  Product prod;
  if (expr->is<Product>()) {
    prod = expr->as<Product>();
  } else if (expr->is<Tensor>()) {
    prod.append(1, expr, Product::Flatten::No);
  } else {
    throw std::runtime_error("kramers_trace_burnside: unsupported Expr type");
  }

  const auto indices = collect_kramers_indices(expr);
  const std::size_t n = indices.size();
  if (n == 0) return expr;
  if (n > 30) {
    throw std::runtime_error(
        "kramers_trace_burnside: too many Kramers indices (> 2^30)");
  }
  const auto tensors = product_tensors(prod);

  // 1. Find transposition generators of the index-permutation symmetry
  //    group: (a,b) qualifies iff in every tensor that contains both
  //    they share an antisymmetric bra/ket group (sign +1 net).
  container::svector<std::pair<std::size_t, std::size_t>> generators;
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = a + 1; b < n; ++b) {
      std::size_t n_common = 0;
      bool ok = true;
      for (const Tensor* t : tensors) {
        const bool ha = tensor_has_index(*t, indices[a]);
        const bool hb = tensor_has_index(*t, indices[b]);
        if (ha && hb) {
          ++n_common;
          if (!same_antisymm_group(*t, indices[a], indices[b])) {
            ok = false;
            break;
          }
        } else if (ha != hb) {
          ok = false;
          break;
        }
      }
      if (ok && n_common > 0 && n_common % 2 == 0)
        generators.emplace_back(a, b);
    }
  }

  // 2. BFS-close the transposition generators into the full group.
  container::svector<container::svector<std::size_t>> group;
  {
    container::svector<std::size_t> identity(n);
    for (std::size_t k = 0; k < n; ++k) identity[k] = k;
    container::set<container::svector<std::size_t>> seen;
    container::svector<container::svector<std::size_t>> frontier{identity};
    seen.insert(identity);
    while (!frontier.empty()) {
      container::svector<container::svector<std::size_t>> next;
      for (auto const& g : frontier) {
        group.push_back(g);
        for (auto const& [a, b] : generators) {
          auto ng = g;
          std::swap(ng[a], ng[b]);
          if (seen.insert(ng).second) next.push_back(ng);
        }
      }
      frontier = std::move(next);
    }
  }

  // 3. Burnside orbit enumeration over 2^n bit-string configs.
  auto apply_perm = [](std::uint64_t cfg,
                       const container::svector<std::size_t>& perm) {
    std::uint64_t out = 0;
    for (std::size_t k = 0; k < perm.size(); ++k)
      if ((cfg >> k) & 1u) out |= (std::uint64_t{1} << perm[k]);
    return out;
  };
  const std::uint64_t n_configs = std::uint64_t{1} << n;
  container::svector<bool> visited(n_configs, false);
  container::svector<std::pair<std::uint64_t, std::size_t>> reps;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    if (visited[cfg]) continue;
    container::set<std::uint64_t> orbit;
    for (auto const& g : group) orbit.insert(apply_perm(cfg, g));
    for (auto o : orbit) visited[o] = true;
    reps.emplace_back(*orbit.begin(), orbit.size());
  }

  // 4. Per orbit rep: substitute, simplify, scale by orbit size.
  auto result = std::make_shared<Sum>();
  for (auto const& [cfg, orbit_size] : reps) {
    auto term =
        specialize_and_simplify(expr, kramers_replacement_map(indices, cfg));
    if (orbit_size != 1) {
      term =
          ex<Constant>(rational{static_cast<std::intmax_t>(orbit_size)}) * term;
      non_canon_simplify(term);
    }
    result->append(term);
  }
  ExprPtr r = result;
  return simplify_with_reset(r);
}

// =====================================================================
// Complex (Re/Im) split. See spinor.hpp for the design overview.
// =====================================================================

namespace {

/// Returns (re, im) such that @p expr == re + i*im, with both `re` and
/// `im` expressed over real-valued tensors (labels carrying `~r`/`~i`).
std::pair<ExprPtr, ExprPtr> split_re_im(const ExprPtr& expr) {
  using S = Constant::scalar_type;
  if (expr->is<Constant>()) {
    const auto v = expr->as<Constant>().value();
    return {ex<Constant>(S{v.real()}), ex<Constant>(S{v.imag()})};
  }
  if (expr->is<RealPart>()) {
    auto inner = split_re_im(expr->as<RealPart>().inner());
    return {inner.first, ex<Constant>(S{0})};
  }
  if (expr->is<ImagPart>()) {
    auto inner = split_re_im(expr->as<ImagPart>().inner());
    return {inner.second, ex<Constant>(S{0})};
  }
  if (expr->is<Tensor>()) {
    auto const& t = expr->as<Tensor>();
    const bool conj = has_conj_suffix(t.label());
    const std::wstring core =
        conj ? std::wstring{t.label().substr(0, t.label().size() - 1)}
             : std::wstring{t.label()};
    ExprPtr re = tensor_with_label(t, complex_label_add(core, /*imag*/ false));
    ExprPtr im = tensor_with_label(t, complex_label_add(core, /*imag*/ true));
    if (conj) im = ex<Constant>(S{-1}) * im;  // g* = g~r − i·g~i
    return {re, im};
  }
  if (expr->is<Sum>()) {
    auto re = std::make_shared<Sum>();
    auto im = std::make_shared<Sum>();
    for (auto&& s : *expr) {
      auto part = split_re_im(s);
      re->append(part.first);
      im->append(part.second);
    }
    return {re, im};
  }
  if (expr->is<Product>()) {
    auto const& p = expr->as<Product>();
    auto acc = split_re_im(ex<Constant>(p.scalar()));
    for (auto&& f : p) {
      auto fac = split_re_im(f);
      // (ar + i·ai)(fr + i·fi) = (ar·fr − ai·fi) + i(ar·fi + ai·fr)
      ExprPtr nr = acc.first * fac.first +
                   ex<Constant>(S{-1}) * (acc.second * fac.second);
      ExprPtr ni = acc.first * fac.second + acc.second * fac.first;
      acc = {nr, ni};
    }
    return acc;
  }
  if (expr->is<Variable>())
    throw std::runtime_error("complex_split: Variable not supported");
  throw std::runtime_error("complex_split: unsupported Expr type");
}

/// Distribute and constant-fold a real-tensor scalar expression.
///
/// Deliberately does NOT canonicalize (would permute Antisymm `~r`/`~i`
/// indices and emit signs the evaluator cannot see). `non_canon_simplify`
/// only distributes and constant-folds; it never permutes tensor indices.
ExprPtr clean_real(ExprPtr e) { return non_canon_simplify(e); }

}  // namespace

ExprPtr complex_split(const ExprPtr& expr) {
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& summand : *expr) s->append(complex_split(summand));
    return s;
  }
  if (expr->is<RealPart>()) {
    return ex<RealPart>(
        clean_real(split_re_im(expr->as<RealPart>().inner()).first));
  }
  if (expr->is<ImagPart>()) {
    // Im(inner) is itself real-valued; its real-tensor form is
    // `.second` of the split — wrap in RealPart for the evaluator.
    return ex<RealPart>(
        clean_real(split_re_im(expr->as<ImagPart>().inner()).second));
  }
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;
  if (expr->is<Product>()) {
    auto const& p = expr->as<Product>();
    if (p.size() == 1 &&
        (p.factor(0)->is<RealPart>() || p.factor(0)->is<ImagPart>())) {
      auto prod = std::make_shared<Product>();
      prod->scale(p.scalar());
      prod->append(1, complex_split(p.factor(0)), Product::Flatten::No);
      return prod;
    }
    throw std::runtime_error(
        "complex_split: expected `Constant * (Re|Im)(...)` summand — run "
        "fold_conj_pairs first so every term is real-wrapped");
  }
  throw std::runtime_error("complex_split: unsupported top-level Expr type");
}

// =====================================================================
// V1-canonicalization variant of the open-shell tracer.
//
// Mirrors `Product::canonicalize_impl` (expr.cpp:174-254) but forces
// the all-tensor sub-block through `TensorNetworkV1` rather than the
// hard-coded `using TN = TensorNetwork (= V3)`. This is the minimum
// surface area needed for V1 canonicalization without modifying
// SeQuant core; the rest of the pipeline (rapid_simplify, expand,
// flattening, hashing) is V-agnostic.
// =====================================================================

namespace {

/// Canonicalize one Product through TensorNetworkV1. Mirrors
/// `Product::canonicalize_impl`'s all-tensor branch (expr.cpp:216-254)
/// with the TN typedef forced to V1. Falls back to the input unchanged
/// for Products that contain non-tensor non-scalar factors (e.g.
/// embedded RealPart/ImagPart wrappers — the open-shell pipeline never
/// produces these mid-pass).
ExprPtr canonicalize_v1_product(Product& product, CanonicalizeOptions opts) {
  using NamedIndexSet = tensor_network::NamedIndexSet;

  // recursively canonicalize non-tensor subfactors first
  for (auto& factor : product.factors()) {
    if (factor->is<AbstractTensor>()) continue;
    auto bp = factor->canonicalize(opts);
    if (bp) {
      SEQUANT_ASSERT(bp->is<Constant>());
      product.scale(std::static_pointer_cast<Constant>(bp)->value());
    }
  }

  // bail if any non-tensor non-scalar factors remain
  for (const auto& factor : product.factors()) {
    if (!factor->is<AbstractTensor>() && !factor->is_scalar()) return {};
  }

  // pull tensor factors out into their own range for TN construction.
  // We ignore scalar factors (Constants embedded as factors): they are
  // multiplied into the scalar at the end of this pass.
  container::svector<ExprPtr> tensor_factors;
  Constant::scalar_type folded_scalar{1};
  for (const auto& factor : product.factors()) {
    if (factor->is<AbstractTensor>()) {
      tensor_factors.push_back(factor);
    } else if (factor->is<Constant>()) {
      folded_scalar = folded_scalar * factor->as<Constant>().value();
    }
  }

  if (tensor_factors.empty()) return {};

  TensorNetworkV1 tn(tensor_factors);

  std::shared_ptr<NamedIndexSet> named_indices;
  if (opts.named_indices) {
    named_indices = std::make_shared<NamedIndexSet>(opts.named_indices->begin(),
                                                    opts.named_indices->end());
  }

  ExprPtr canon_factor =
      tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                      /*fast=*/opts.method == CanonicalizationMethod::Rapid,
                      named_indices.get());

  // Rebuild factors_ in canonical order.
  auto& factors = product.factors();
  factors.clear();
  for (const auto& tptr : tn.tensors()) {
    auto exprptr = std::dynamic_pointer_cast<Expr>(tptr);
    SEQUANT_ASSERT(exprptr);
    factors.push_back(exprptr);
  }

  product.scale(folded_scalar);
  if (canon_factor) product.scale(canon_factor->as<Constant>().value());
  return {};
}

/// Recurse Sum → for each summand → if Product, canonicalize via V1;
/// then rebuild the Sum so equivalent summands fold by hash. Mirrors
/// the structure of `Sum::canonicalize_impl` (expr.cpp:374-432) but
/// routes Product canonicalization through V1 rather than V3.
ExprPtr canonicalize_v1_recurse(const ExprPtr& expr, CanonicalizeOptions opts) {
  if (expr->is<Sum>()) {
    auto out = std::make_shared<Sum>();
    for (const auto& s : expr->as<Sum>().summands()) {
      out->append(canonicalize_v1_recurse(s, opts));
    }
    return out;
  }
  if (expr->is<Product>()) {
    auto product = std::static_pointer_cast<Product>(expr->clone());
    canonicalize_v1_product(*product, opts);
    return product;
  }
  // Tensor / Constant / Variable / wrapper nodes: nothing to do.
  return expr;
}

}  // namespace

ExprPtr& canonicalize_v1(ExprPtr& expr, CanonicalizeOptions opts) {
  expr = canonicalize_v1_recurse(expr, opts);
  return expr;
}

ExprPtr& simplify_v1(ExprPtr& expr) {
  expand(expr);
  rapid_simplify(expr);
  canonicalize_v1(expr);
  rapid_simplify(expr);
  return expr;
}

// =====================================================================
// open_shell_spintrace_v1: V1-driven clone of
// `mbpt::open_shell_spintrace_impl` (spin.cpp:1367).
//
// Both V1 and V3 enforce bra↔ket-only edge connectivity in their
// graph builders, so an expression with a dummy summed across two
// bras (or two kets) — like the PNO pair density
// `t̄^{ac}_{ij} t̄^{bc}_{ij}` (`c` in two bras) — trips the assertion
// regardless of which TN version we use. The closed-shell tracer
// `mbpt::spintrace` (spin.cpp:1626) avoids the issue by NEVER
// canonicalizing the spin-traced products: it only does `expand` +
// `rapid_simplify`, both of which are TN-agnostic. We mirror that
// strategy here — `simplify_v1` and `canonicalize_v1` still exist as
// useful primitives for callers whose expressions DO have CC-style
// wiring, but in this open-shell pipeline we deliberately skip those
// steps and stay TN-free.
// =====================================================================

namespace {

/// True if a fully-flattened Product is spin-symmetric, i.e. has the
/// concatenated bra and ket lists carrying the same per-slot QN.
/// Mirrors `spin_symm_product` lambda in `open_shell_spintrace_impl`.
bool spin_symm_product_v1(const Product& product) {
  container::svector<Index> cBra, cKet;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto const& t = term->as<Tensor>();
      cBra.insert(cBra.end(), t.bra().begin(), t.bra().end());
      cKet.insert(cKet.end(), t.ket().begin(), t.ket().end());
    } else if (term->is<Product>() || term->is<Sum>()) {
      throw Exception(
          "open_shell_spintrace_v1: nested Product/Sum unsupported in "
          "spin_symm_product");
    }
  }
  if (cKet.size() != cBra.size()) return false;
  auto i_ket = cKet.begin();
  for (auto& b : cBra) {
    if (b.space().qns() != i_ket->space().qns()) return false;
    ++i_ket;
  }
  return true;
}

}  // namespace

std::vector<ExprPtr> open_shell_spintrace_v1(
    const ExprPtr& expr,
    container::svector<container::svector<Index>> ext_index_groups,
    std::optional<int> target_spin_case) {
  if (expr->is<Constant>() || expr->is<Variable>())
    return std::vector<ExprPtr>{expr};

  // Internal vs external indices.
  container::set<Index, Index::LabelCompare> grand_idxlist =
      get_used_indices<container::set<Index, Index::LabelCompare>>(expr);
  container::set<Index> ext_idxlist;
  for (const auto& grp : ext_index_groups) {
    for (const auto& idx : grp) {
      Index ix = idx;
      ix.reset_tag();
      ext_idxlist.insert(std::move(ix));
    }
  }
  container::set<Index> int_idxlist;
  for (auto&& gidx : grand_idxlist) {
    if (ext_idxlist.find(gidx) == ext_idxlist.end()) int_idxlist.insert(gidx);
  }

  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> int_index_groups;
  for (auto&& i : int_idxlist) int_index_groups.emplace_back(IndexGroup(1, i));

  SEQUANT_ASSERT(grand_idxlist.size() ==
                 int_idxlist.size() + ext_idxlist.size());

  auto make_spinspecific = [](const Index& idx, long int spin_bit) {
    return spin_bit == 0 ? make_spinalpha(idx) : make_spinbeta(idx);
  };

  // Internal-index replacements: 2^(n_internal_groups) configurations.
  auto spin_cases = [&](const container::svector<IndexGroup>& idx_group) {
    const auto ncases = std::uint64_t{1} << idx_group.size();
    container::svector<container::map<Index, Index>> all(ncases);
    for (std::uint64_t i = 0; i < ncases; ++i) {
      container::map<Index, Index> rep;
      for (std::size_t g = 0; g < idx_group.size(); ++g) {
        const auto bit = (i >> g) & 1u;
        for (auto& idx : idx_group[g])
          rep.emplace(idx, make_spinspecific(idx, bit));
      }
      all[i] = rep;
    }
    return all;
  };

  // External-index replacements: pad the alpha-prefix with k betas for
  // k = 0..n_groups (n_groups+1 distinct M_S sectors).
  auto ext_spin_cases = [&](const auto& idx_groups) {
    container::svector<container::map<Index, Index>> all;
    for (std::size_t i = 0; i <= idx_groups.size(); ++i) {
      container::svector<int> spins(idx_groups.size(), 0);
      std::fill(spins.end() - i, spins.end(), 1);
      container::map<Index, Index> rep;
      for (std::size_t j = 0; j < idx_groups.size(); ++j) {
        for (const auto& idx : idx_groups[j])
          rep.emplace(idx, make_spinspecific(idx, spins[j]));
      }
      all.push_back(rep);
    }
    return all;
  };

  auto i_rep = spin_cases(int_index_groups);
  auto e_rep = ext_spin_cases(ext_index_groups);
  if (target_spin_case) {
    auto pick = e_rep.at(*target_spin_case);
    e_rep.clear();
    e_rep.push_back(pick);
  }

  // Expand antisymmetrizers + flatten — but NO canonicalize. The
  // expression may contain bra↔bra (or ket↔ket) summed dummies that
  // both TN versions reject; staying TN-free keeps the pipeline
  // tolerant. Mirrors `spintrace_impl` (spin.cpp:1626).
  auto expanded = expand_A_op(expr);
  detail::reset_idx_tags(expanded);
  expand(expanded);
  rapid_simplify(expanded);

  std::vector<ExprPtr> result;
  for (auto& e : e_rep) {
    auto spin_expr = append_spin(expanded, e);
    detail::reset_idx_tags(spin_expr);
    Sum e_result{};
    for (auto& i : i_rep) {
      ExprPtr spin_expr_i = append_spin(spin_expr, i);
      spin_expr_i = expand_antisymm(spin_expr_i, /*skip_spinsymm=*/true);
      expand(spin_expr_i);
      rapid_simplify(spin_expr_i);
      detail::reset_idx_tags(spin_expr_i);
      Sum i_result{};
      if (spin_expr_i->is<Tensor>() || spin_expr_i->is<Constant>() ||
          spin_expr_i->is<Variable>()) {
        e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Product>()) {
        if (spin_symm_product_v1(spin_expr_i->as<Product>()))
          e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Sum>()) {
        for (auto& pr : *spin_expr_i) {
          if (pr->is<Product>()) {
            if (spin_symm_product_v1(pr->as<Product>())) i_result.append(pr);
          } else if (pr->is<Tensor>()) {
            if (ms_conserving_columns(pr->as<Tensor>())) i_result.append(pr);
          } else if (pr->is<Constant>() || pr->is<Variable>()) {
            i_result.append(pr);
          } else {
            throw Exception(
                "open_shell_spintrace_v1: unsupported summand type");
          }
        }
        e_result.append(std::make_shared<Sum>(i_result));
      }
    }
    result.push_back(std::make_shared<Sum>(e_result));
  }

  if (target_spin_case) {
    SEQUANT_ASSERT(result.size() == 1 &&
                   "Spin-specific case must return one expression");
  }

  // Final per-spin-block clean-up: rapid_simplify only (no canonicalize).
  for (auto& expression : result) {
    detail::reset_idx_tags(expression);
    rapid_simplify(expression);
  }
  return result;
}

}  // namespace sequant::mbpt
