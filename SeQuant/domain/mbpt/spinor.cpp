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
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace sequant::mbpt {

namespace detail {

/// Construct a new Index in the same IndexSpace type as @p idx, with
/// the @p k Kramers QN bit set (replacing any existing Kramers state)
/// and the label decorated with the matching ⇑/⇓ glyph. Spin bits
/// (alpha/beta) are NOT touched here — Spin and Kramers are mutually
/// exclusive on a single index, but the per-index check is left to
/// the caller (a typical workflow only ever sets one of them).
///
/// Mirrors `detail::make_index_with_spincase` in `spin.cpp` so that
/// SeQuant's IndexSpace registry is reused: if a space with the
/// resulting label and qns already exists, we adopt it.
Index make_index_with_kramers(const Index& idx, mbpt::Kramers k) {
  // sanity check: at most one Kramers annotation already
  SEQUANT_ASSERT(!(idx.label().find(L'⇑') != std::wstring::npos &&
                   idx.label().find(L'⇓') != std::wstring::npos));

  // Strip any existing Kramers bits, then OR in the requested state.
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
                       // assume size does not depend on Kramers state
                       idx.space().approximate_size()};
  }
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_kramers(pidx, k);
  return Index{space, idx.ordinal(), protoindices};
}

}  // namespace detail

Index make_kramers_up(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::up);
}

Index make_kramers_dn(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::down);
}

Index make_kramers_free(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::any);
}

// =====================================================================
// Conjugation infrastructure (groundwork for Kramers tracer + quaternion
// decomposer).
// =====================================================================

std::wstring toggle_conj_suffix(std::wstring_view label) {
  if (has_conj_suffix(label))
    return std::wstring{label.substr(0, label.size() - 1)};
  return std::wstring{label} + L'*';
}

namespace {

/// Apply conjugation to a single Tensor: rebuild it with toggled label
/// and identical bra/ket/aux/symmetry/etc.
ExprPtr conjugate_tensor(const Tensor& t) {
  auto new_label = toggle_conj_suffix(t.label());
  // Tensor copy ctor + label edit isn't directly exposed; reconstruct.
  return ex<Tensor>(
      new_label, bra(container::svector<Index>{t.bra().begin(), t.bra().end()}),
      ket(container::svector<Index>{t.ket().begin(), t.ket().end()}),
      aux(container::svector<Index>{t.aux().begin(), t.aux().end()}),
      t.symmetry(), t.braket_symmetry(), t.column_symmetry());
}

}  // namespace

ExprPtr conjugate(const ExprPtr& expr) {
  if (expr->is<Tensor>()) {
    return conjugate_tensor(expr->as<Tensor>());
  }
  if (expr->is<Constant>()) {
    auto v = expr->as<Constant>().value();
    return ex<Constant>(sequant::conj(v));
  }
  if (expr->is<Variable>()) {
    throw std::runtime_error(
        "conjugate: Variable conjugation not yet supported");
  }
  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    auto p = std::make_shared<Product>();
    p->scale(sequant::conj(prod.scalar()));
    for (auto&& f : prod) p->append(1, conjugate(f));
    return p;
  }
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& summand : *expr) s->append(conjugate(summand));
    return s;
  }
  throw std::runtime_error("conjugate: unsupported Expr type");
}

ExprPtr append_kramers(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements) {
  auto add_to_tensor = [&](const Tensor& tensor) {
    auto t = std::make_shared<Tensor>(tensor);
    t->transform_indices(index_replacements);
    return t;
  };
  auto add_to_product = [&](const Product& product) {
    auto p = std::make_shared<Product>();
    p->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        p->append(1, add_to_tensor(term->as<Tensor>()));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        p->append(1, term);
      } else {
        throw std::runtime_error(
            "append_kramers: unsupported Expr type in product: " +
            term->type_name());
      }
    }
    return p;
  };

  if (expr->is<Tensor>()) return add_to_tensor(expr->as<Tensor>());
  if (expr->is<Product>()) return add_to_product(expr->as<Product>());
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      s->append(append_kramers(summand, index_replacements));
    }
    return s;
  }
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;
  throw std::runtime_error("append_kramers: unsupported Expr type");
}

namespace {

/// Walk @p expr and collect the unique Kramers-eligible (i.e.,
/// `Kramers::any`) indices in encounter order. Indices that already
/// have a specific Kramers state (up/down) are skipped.
container::svector<Index> collect_kramers_indices(const ExprPtr& expr) {
  container::svector<Index> seen;
  container::set<Index> seen_set;
  auto consider = [&](const Index& idx) {
    auto qns = idx.space().qns().to_int32();
    if ((qns & mask_v<Kramers>) != static_cast<bitset_t>(Kramers::any)) {
      // Either no Kramers state set or already specialized — skip.
      return;
    }
    if (seen_set.insert(idx).second) seen.push_back(idx);
  };
  expr->visit(
      [&](const ExprPtr& e) {
        if (e->is<Tensor>()) {
          for (const auto& idx : e->as<Tensor>().const_braket()) {
            consider(idx);
          }
        }
      },
      /* recursive = */ true);
  return seen;
}

}  // namespace

ExprPtr kramers_trace(const ExprPtr& expr) {
  // Enumerate the 2^N Kramers configurations of the expression's spinor
  // (Kramers::any) indices, then apply TRS pairing: configs (X, X̄) related
  // by full bar-flip are complex conjugates of each other (per the
  // (-1)^k full-bar identity, with the per-tensor signs cancelling for
  // products in which every barred index appears in an even number of
  // tensors — true for any all-internal-index contraction).
  //
  // For each TRS pair we materialize only the "lex-smaller" configuration
  // X explicitly; the X̄ partner is emitted as `conjugate(T_X)`. Tensors
  // in the conj branch carry a `*` label suffix (per the convention in
  // the conjugate/real_part/imaginary_part scaffolding above) which the
  // evaluator dispatches into a `.conj()` call on the underlying numeric
  // tensor. Term count is unchanged, but the symbolic structure exposes
  // the conjugate-pair relationship for downstream simplification (e.g.,
  // `real_part` collapse in real-scalar callers).
  //
  // Standard SeQuant simplifications (dummy renaming, antisym fold, etc.)
  // are NOT applied here — caller invokes `sequant::canonicalize` if it
  // wants those.

  const auto indices = collect_kramers_indices(expr);
  const std::size_t n = indices.size();
  if (n == 0) return expr;
  if (n > 30) {
    throw std::runtime_error(
        "kramers_trace: too many Kramers indices (would emit > 2^30 terms)");
  }

  const std::uint64_t n_configs = 1ull << n;

  // Enumerate all 2^n Kramers configurations as separate Products. We
  // mirror open_shell_spintrace's pipeline: substitute → expand_antisymm →
  // expand → simplify → canonicalize, working entirely in expanded
  // (NonSymm) tensor form so the standard SeQuant simplifier sees full
  // index permutability for fold detection.
  auto result = std::make_shared<Sum>();
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    container::map<Index, Index> repl;
    for (std::size_t i = 0; i < n; ++i) {
      const bool is_down = (cfg >> i) & 1u;
      Index new_idx =
          is_down ? make_kramers_dn(indices[i]) : make_kramers_up(indices[i]);
      repl[indices[i]] = new_idx;
    }
    auto term = append_kramers(expr, repl);
    term = expand_antisymm(term, /*skip_spinsymm*/ false);
    expand(term);
    rapid_simplify(term);
    canonicalize(term);
    rapid_simplify(term);
    result->append(term);
  }
  // Final cross-product canonicalize: now that every term is in NonSymm
  // expanded form, dummy renaming can fold equivalents across Kramers
  // configurations.
  ExprPtr r = result;
  expand(r);
  rapid_simplify(r);
  canonicalize(r);
  rapid_simplify(r);
  return r;
}

// =====================================================================
// Conjugate-pair fold pass: detects TRS-related Products in a Sum and
// rewrites (A + B) → 2·Re(A) where B = flip_kramers(A).
// =====================================================================

ExprPtr flip_kramers(const ExprPtr& expr) {
  // Walk all tensors and collect every Kramers-tagged index. For each
  // distinct (label, ordinal, space-type) seen, build a partner index
  // with the opposite Kramers state. Then apply the swap via the existing
  // append_kramers replacement machinery.
  container::map<Index, Index> repl;
  expr->visit(
      [&](const ExprPtr& e) {
        if (!e->is<Tensor>()) return;
        for (const auto& idx : e->as<Tensor>().const_braket()) {
          if (repl.contains(idx)) continue;
          const auto qns_int = idx.space().qns().to_int32();
          const auto k_bits = qns_int & mask_v<Kramers>;
          if (k_bits == static_cast<bitset_t>(Kramers::up)) {
            repl[idx] = make_kramers_dn(idx);
          } else if (k_bits == static_cast<bitset_t>(Kramers::down)) {
            repl[idx] = make_kramers_up(idx);
          }
          // Kramers::any or no Kramers state — leave alone.
        }
      },
      /* recursive = */ true);
  if (repl.empty()) return expr;
  return append_kramers(expr, repl);
}

namespace {

/// Normalize a Product-like expression for structural comparison: clone,
/// apply canonicalize + rapid_simplify so dummy renaming and ordering
/// match other already-canonicalized summands.
ExprPtr normalize(const ExprPtr& e) {
  ExprPtr c = e->clone();
  expand(c);
  rapid_simplify(c);
  canonicalize(c);
  rapid_simplify(c);
  return c;
}

/// Compare two scalars as complex rationals.
bool scalar_eq(const Constant::scalar_type& a, const Constant::scalar_type& b) {
  return a == b;
}

/// Extract the leading scalar of a Product, or 1 if the expr is a bare
/// Tensor.
Constant::scalar_type product_scalar(const ExprPtr& e) {
  if (e->is<Product>()) return e->as<Product>().scalar();
  return Constant::scalar_type{1};
}

/// Build a "factor-only" copy of a Product (drop the leading scalar).
/// Returns @p e unchanged if it isn't a Product.
ExprPtr product_factors_only(const ExprPtr& e) {
  if (!e->is<Product>()) return e->clone();
  auto p = std::make_shared<Product>();
  p->scale(Constant::scalar_type{1});
  for (auto&& f : e->as<Product>()) p->append(1, f->clone());
  return p;
}

}  // namespace

ExprPtr fold_conj_pairs(const ExprPtr& expr) {
  if (!expr->is<Sum>()) return expr;
  auto const& sum = expr->as<Sum>();
  const std::size_t n = sum.size();

  // O(n) canonical-form hash bucketing. For each summand A, compute a
  // TRS-invariant key by taking the lex-smaller of normalize(A) and
  // normalize(flip_kramers(A)). A and its full-bar partner B = flip(A)
  // both map to the same key, so they land in the same bucket.
  //
  // Within a bucket we re-verify structural equality (hash collisions are
  // possible; equality is O(tensor-count) per pair). Because TRS pairs
  // form orbits of size ≤ 2, buckets are tiny and per-bucket work is O(1)
  // on average. Total cost: O(n) normalizations + O(n) flips + O(n) hash
  // inserts, vs. O(n²) for naive pairwise search.
  //
  // Scalar handling: we also bucket by scalar relationship. Two summands
  // with matching scalars (a, a) collapse to `2·Re(A)`; opposite scalars
  // (a, -a) collapse to `2i·Im(A)`. Self-conjugate singletons emit
  // `Re(A)` (real-valued).
  struct BucketEntry {
    std::size_t idx;
    ExprPtr factors_normalized;    // factor part, post-normalize
    Constant::scalar_type scalar;  // leading scalar prefactor
    bool self_conj;                // factors == flip(factors)
  };
  container::map<Expr::hash_type, container::svector<BucketEntry>> buckets;

  std::vector<ExprPtr> originals;
  originals.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    originals.push_back(sum.summand(i)->clone());
    auto factors = normalize(product_factors_only(originals[i]));
    auto factors_flipped = normalize(flip_kramers(factors));
    const bool self_conj = (*factors == *factors_flipped);
    // Canonical factor key: lex-min of (factors, factors_flipped). Both A
    // and flip(A) produce the same key, so they collide in the same bucket.
    ExprPtr canon = (*factors < *factors_flipped) ? factors : factors_flipped;
    const auto h = canon->hash_value();
    buckets[h].push_back(
        BucketEntry{i, factors, product_scalar(originals[i]), self_conj});
  }

  auto out = std::make_shared<Sum>();
  std::vector<bool> used(n, false);

  // Helper to emit a folded summand based on scalar relationship.
  auto emit_pair = [&](const BucketEntry& a, const BucketEntry& b) {
    if (scalar_eq(a.scalar, b.scalar)) {
      // (A + B) = 2 Re(A)
      out->append(ex<Constant>(rational{2}) * real_part(originals[a.idx]));
    } else if (scalar_eq(a.scalar, -b.scalar)) {
      // (A - B) = 2i Im(A)  (B carries -A's scalar)
      using cscalar = Constant::scalar_type;
      out->append(ex<Constant>(cscalar{Complex<rational>{0, 2}}) *
                  imaginary_part(originals[a.idx]));
    } else {
      // Same factor-canonical key but unrelated scalars — emit both raw.
      out->append(originals[a.idx]);
      out->append(originals[b.idx]);
    }
    used[a.idx] = used[b.idx] = true;
  };

  for (auto& [_, entries] : buckets) {
    // Verify hash-bucket entries actually share the same canonical factor
    // form (defend against hash collisions): partition by structural
    // equality on factors_normalized OR its flip.
    std::vector<bool> taken(entries.size(), false);
    for (std::size_t i = 0; i < entries.size(); ++i) {
      if (taken[i] || used[entries[i].idx]) continue;
      const auto& a = entries[i];
      // Self-conjugate singleton in this bucket?
      if (a.self_conj) {
        // Look for another self-conj with same scalar — if none, emit Re(A).
        bool merged = false;
        for (std::size_t j = i + 1; j < entries.size(); ++j) {
          if (taken[j] || used[entries[j].idx]) continue;
          const auto& b = entries[j];
          if (b.self_conj && *b.factors_normalized == *a.factors_normalized) {
            // Two self-conj copies → 2·Re(A) (or scalar-aware emit).
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
      // Look for the partner: another entry whose factors equal flip(A).
      // We don't have flip(A) materialized here — use the structural test
      // (a.factors and b.factors land in the same bucket only when one is
      // the flip of the other).
      bool paired = false;
      for (std::size_t j = i + 1; j < entries.size(); ++j) {
        if (taken[j] || used[entries[j].idx]) continue;
        const auto& b = entries[j];
        // Confirm they're actually a TRS pair: b.factors should equal
        // flip(a.factors).
        auto a_flip = normalize(flip_kramers(a.factors_normalized));
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
// Cycle decomposition. See spinor.hpp for the design summary.
// =====================================================================

container::svector<ContractionCycle> kramers_cycles(const Product& product) {
  // Collect (tensor_idx, braket_pos, Index) for every position in every
  // Tensor factor. Indices that aren't Tensor are skipped (Constants,
  // Variables don't contribute to cycles).
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
    for (auto const& idx : t.const_braket()) {
      slots.push_back(Slot{tensor_counter, pos, idx});
      ++pos;
    }
    ++tensor_counter;
  }
  const std::size_t N = slots.size();
  container::svector<bool> visited(N, false);

  // Pre-build two adjacency maps:
  // - intra-tensor: bra[k] ↔ ket[k] of the same tensor (positions 0↔rank,
  //   1↔rank+1, ... where rank = bra_rank).
  // - inter-tensor: same Index occurrence in two distinct slots.
  container::svector<std::size_t> intra_partner(N, N);  // N = "no partner"
  container::svector<std::size_t> inter_partner(N, N);

  // intra-tensor pairing per tensor: walk positions and pair k ↔ k+rank/2
  // within each tensor's slots. We rebuild rank by counting consecutive
  // slots with the same tensor_idx.
  for (std::size_t i = 0; i < N;) {
    std::size_t j = i;
    while (j < N && slots[j].tensor_idx == slots[i].tensor_idx) ++j;
    const std::size_t tensor_rank = j - i;  // total bra+ket count
    SEQUANT_ASSERT(tensor_rank % 2 == 0);
    const std::size_t half = tensor_rank / 2;
    for (std::size_t k = 0; k < half; ++k) {
      intra_partner[i + k] = i + half + k;
      intra_partner[i + half + k] = i + k;
    }
    i = j;
  }

  // inter-tensor pairing: same Index in two distinct positions across the
  // whole Product. Use Index-equality (which compares full label).
  for (std::size_t i = 0; i < N; ++i) {
    if (inter_partner[i] != N) continue;
    for (std::size_t j = i + 1; j < N; ++j) {
      if (inter_partner[j] != N) continue;
      if (slots[j].idx == slots[i].idx) {
        inter_partner[i] = j;
        inter_partner[j] = i;
        break;
      }
    }
  }

  // Walk cycles. From each unvisited slot, alternate intra → inter →
  // intra → ... edges until returning to the start.
  container::svector<ContractionCycle> cycles;
  for (std::size_t start = 0; start < N; ++start) {
    if (visited[start]) continue;
    ContractionCycle cyc;
    std::size_t cur = start;
    bool take_intra = true;  // first hop is intra-tensor
    while (!visited[cur]) {
      visited[cur] = true;
      cyc.nodes.push_back(ContractionCycle::Node{
          slots[cur].tensor_idx, slots[cur].braket_pos, slots[cur].idx});
      const std::size_t next =
          take_intra ? intra_partner[cur] : inter_partner[cur];
      if (next == N) break;  // dangling (shouldn't happen for closed expr)
      cur = next;
      take_intra = !take_intra;
    }
    cycles.push_back(std::move(cyc));
  }

  return cycles;
}

namespace {

/// Format a cycle's Kramers labelling as a wstring like "UAUA" or "UDDU"
/// where each char is the visited index's Kramers state.
std::wstring cycle_kramers_label(const ContractionCycle& c) {
  std::wstring s;
  s.reserve(c.nodes.size());
  for (auto const& n : c.nodes) {
    const auto qns_int = n.idx.space().qns().to_int32();
    const auto k_bits = qns_int & mask_v<Kramers>;
    if (k_bits == static_cast<bitset_t>(Kramers::up)) {
      s += L'U';
    } else if (k_bits == static_cast<bitset_t>(Kramers::down)) {
      s += L'B';
    } else {
      s += L'?';
    }
  }
  return s;
}

/// Number of canonical Kramers patterns of a length-L cycle under
/// cyclic-rotation symmetry. By Burnside: (1/L) Σ_{d | L} φ(L/d) · 2^d.
std::size_t canonical_pattern_count_cyclic(std::size_t L) {
  if (L == 0) return 1;
  // φ(n)
  auto phi = [](std::size_t n) {
    std::size_t result = n;
    for (std::size_t p = 2; p * p <= n; ++p) {
      if (n % p == 0) {
        while (n % p == 0) n /= p;
        result -= result / p;
      }
    }
    if (n > 1) result -= result / n;
    return result;
  };
  std::size_t total = 0;
  for (std::size_t d = 1; d <= L; ++d) {
    if (L % d == 0) {
      total += phi(L / d) * (std::size_t{1} << d);
    }
  }
  return total / L;
}

}  // namespace

void kramers_cycle_dump(const ExprPtr& expr, std::wostream& os) {
  auto dump_product = [&](const Product& p, std::size_t pidx) {
    os << L"  Product[" << pidx << L"]:";
    for (auto&& f : p) {
      if (f->is<Tensor>()) {
        auto const& t = f->as<Tensor>();
        os << L" " << t.label();
      }
    }
    os << L"\n";
    auto cycles = kramers_cycles(p);
    os << L"    " << cycles.size() << L" cycle(s):\n";
    for (std::size_t ci = 0; ci < cycles.size(); ++ci) {
      auto const& c = cycles[ci];
      // Distinct indices visited.
      container::set<Index> distinct;
      for (auto const& n : c.nodes) distinct.insert(n.idx);
      const auto canon = canonical_pattern_count_cyclic(c.nodes.size());
      os << L"      cycle[" << ci << L"] len=" << c.nodes.size()
         << L" distinctIdx=" << distinct.size() << L" labels=`"
         << cycle_kramers_label(c) << L"`" << L"  cyclic-canonical-#patterns="
         << canon << L"  walk=";
      bool first = true;
      for (auto const& n : c.nodes) {
        if (!first) os << L"→";
        os << L"t" << n.tensor_idx << L"/" << n.braket_pos << L":"
           << n.idx.label();
        first = false;
      }
      os << L"\n";
    }
  };

  if (expr->is<Product>()) {
    dump_product(expr->as<Product>(), 0);
  } else if (expr->is<Sum>()) {
    auto const& s = expr->as<Sum>();
    for (std::size_t i = 0; i < s.size(); ++i) {
      auto const& term = s.summand(i);
      if (term->is<Product>()) dump_product(term->as<Product>(), i);
    }
  }
}

std::wstring cycle_canonical_label(const ContractionCycle& c) {
  // Build raw label, then take lex-min over all cyclic rotations.
  const auto raw = cycle_kramers_label(c);
  if (raw.empty()) return raw;
  std::wstring best = raw;
  std::wstring rot = raw;
  for (std::size_t k = 1; k < raw.size(); ++k) {
    // rotate by 1: move first char to end
    std::rotate(rot.begin(), rot.begin() + 1, rot.end());
    if (rot < best) best = rot;
  }
  return best;
}

std::wstring cycle_canonical_signature(const Product& product) {
  auto cycles = kramers_cycles(product);
  container::svector<std::wstring> labels;
  labels.reserve(cycles.size());
  for (auto const& c : cycles) labels.push_back(cycle_canonical_label(c));
  std::sort(labels.begin(), labels.end());
  std::wstring sig;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    if (i > 0) sig += L'|';
    sig += labels[i];
  }
  return sig;
}

ExprPtr cycle_canonical_fold(const ExprPtr& expr) {
  if (!expr->is<Sum>()) return expr;
  auto const& sum = expr->as<Sum>();
  // Bucket Sum summands by cycle_canonical_signature. Within each bucket
  // the contractions are numerically equal up to scalar prefactor; sum
  // the prefactors and emit a single representative per bucket.
  container::map<std::wstring, container::svector<std::size_t>> buckets;
  std::vector<ExprPtr> originals;
  originals.reserve(sum.size());
  for (std::size_t i = 0; i < sum.size(); ++i) {
    originals.push_back(sum.summand(i)->clone());
    if (!originals[i]->is<Product>()) {
      // Non-Product summands (Tensor, Constant, Variable) — emit a
      // distinct bucket per index so they pass through unchanged.
      buckets[L"__pass" + std::to_wstring(i)].push_back(i);
      continue;
    }
    const auto sig = cycle_canonical_signature(originals[i]->as<Product>());
    buckets[sig].push_back(i);
  }

  auto out = std::make_shared<Sum>();
  for (auto& [sig, idxs] : buckets) {
    if (idxs.size() == 1) {
      out->append(originals[idxs[0]]);
      continue;
    }
    // Sum scalars; reuse first member's tensor structure as the
    // representative. NOTE: this assumes the bucketed Products are
    // numerically equal as contractions, which holds when the cycles'
    // canonical Kramers labels match — that's the bucket-key meaning.
    Constant::scalar_type total{0};
    for (auto i : idxs) {
      total += product_scalar(originals[i]);
    }
    if (total == Constant::scalar_type{0}) continue;
    auto rep = product_factors_only(originals[idxs.front()]);
    out->append(ex<Constant>(total) * rep);
  }
  return out;
}

// =====================================================================
// Standalone cycle-driven Kramers tracer. See spinor.hpp.
// =====================================================================

namespace {

/// Distinct Kramers-eligible (`Kramers::any`) indices visited by a cycle,
/// in first-encounter order.
container::svector<Index> cycle_distinct_indices(const ContractionCycle& c) {
  container::svector<Index> out;
  container::set<Index> seen;
  for (auto const& n : c.nodes) {
    const auto qns = n.idx.space().qns().to_int32();
    if ((qns & mask_v<Kramers>) != static_cast<bitset_t>(Kramers::any))
      continue;  // already specialized or not Kramers — skip
    if (seen.insert(n.idx).second) out.push_back(n.idx);
  }
  return out;
}

/// Structural shape signature of a cycle: a string encoding the walk's
/// IndexSpace-type sequence, made rotation-invariant by taking the
/// lex-min rotation. Two cycles with the same shape are interchangeable
/// under the expression's tensor (anti)symmetries — the caller is
/// responsible for ensuring those symmetries actually hold.
std::wstring cycle_shape(const ContractionCycle& c) {
  std::wstring raw;
  for (auto const& n : c.nodes) {
    // encode IndexSpace base type as a hex digit (stable per registry)
    raw += static_cast<wchar_t>(L'a' + (n.idx.space().type().to_int32() & 0xF));
  }
  if (raw.empty()) return raw;
  std::wstring best = raw, rot = raw;
  for (std::size_t k = 1; k < raw.size(); ++k) {
    std::rotate(rot.begin(), rot.begin() + 1, rot.end());
    if (rot < best) best = rot;
  }
  return best;
}

/// All 2^n Kramers labelings of a cycle's @p n distinct indices, each a
/// bitset (bit i = 1 ⇔ index i is Kramers-down).
container::svector<std::uint32_t> cycle_labelings(std::size_t n_distinct) {
  container::svector<std::uint32_t> out;
  const std::uint32_t n_cfg = 1u << n_distinct;
  out.reserve(n_cfg);
  for (std::uint32_t c = 0; c < n_cfg; ++c) out.push_back(c);
  return out;
}

/// Multinomial coefficient m! / ∏ count_i! for a multiset assignment.
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
  // Sum: recurse over summands.
  if (expr->is<Sum>()) {
    auto sum_result = std::make_shared<Sum>();
    for (auto&& s : *expr) sum_result->append(kramers_trace_cycles(s));
    ExprPtr r = sum_result;
    expand(r);
    rapid_simplify(r);
    canonicalize(r);
    rapid_simplify(r);
    return r;
  }
  // Constant / Variable: no Kramers indices, pass through.
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Normalize to a Product.
  Product prod;
  if (expr->is<Product>()) {
    prod = expr->as<Product>();
  } else if (expr->is<Tensor>()) {
    prod.append(1, expr, Product::Flatten::No);
  } else {
    throw std::runtime_error("kramers_trace_cycles: unsupported Expr type");
  }

  // 1. Decompose into cycles.
  const auto cycles = kramers_cycles(prod);
  const std::size_t n_cycles = cycles.size();
  if (n_cycles == 0) return expr;

  // Per-cycle distinct Kramers indices + labelings.
  std::vector<container::svector<Index>> cyc_indices(n_cycles);
  std::vector<container::svector<std::uint32_t>> cyc_labelings(n_cycles);
  for (std::size_t c = 0; c < n_cycles; ++c) {
    cyc_indices[c] = cycle_distinct_indices(cycles[c]);
    cyc_labelings[c] = cycle_labelings(cyc_indices[c].size());
  }

  // 2. Group cycles into structural-equivalence classes by shape.
  container::map<std::wstring, container::svector<std::size_t>> shape_classes;
  for (std::size_t c = 0; c < n_cycles; ++c)
    shape_classes[cycle_shape(cycles[c])].push_back(c);

  // 3. Per class, enumerate Kramers-labeling MULTISETS (sorted tuples of
  //    per-cycle labelings) with their multiplicity. A class of m cycles
  //    each with L labelings yields C(L+m-1, m) multisets.
  //
  //    Represent a class assignment as a sorted vector<uint32_t> of
  //    length m (one labeling index per cycle in the class). Multiplicity
  //    = multinomial(counts of each distinct labeling).
  struct ClassConfig {
    container::svector<std::size_t> class_ids;  // cycle indices in this class
    // each entry: (sorted labeling-assignment, multiplicity)
    container::svector<std::pair<container::svector<std::uint32_t>, double>>
        multisets;
  };
  container::svector<ClassConfig> class_configs;
  for (auto& [shape, cyc_ids] : shape_classes) {
    ClassConfig cc;
    cc.class_ids = cyc_ids;
    const std::size_t m = cyc_ids.size();
    // All cycles in a class share labeling-set size (same shape ⇒ same
    // #distinct indices). Use the first cycle's labeling set.
    const auto& labels = cyc_labelings[cyc_ids.front()];
    const std::size_t L = labels.size();
    // Enumerate sorted m-tuples (multisets) from L labelings.
    container::svector<std::size_t> pick(m, 0);
    std::function<void(std::size_t, std::size_t)> rec = [&](std::size_t pos,
                                                            std::size_t start) {
      if (pos == m) {
        // pick[] holds indices into `labels`, non-decreasing.
        container::svector<std::uint32_t> assignment(m);
        for (std::size_t i = 0; i < m; ++i) assignment[i] = labels[pick[i]];
        // multiplicity = multinomial over runs of equal pick values
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
  //    multiset per class. Build the per-index Kramers replacement,
  //    substitute, expand_antisymm, canonicalize, scale by the product
  //    of per-class multiplicities.
  auto result = std::make_shared<Sum>();
  const std::size_t n_classes = class_configs.size();
  container::svector<std::size_t> sel(n_classes, 0);
  std::function<void(std::size_t)> build = [&](std::size_t ci) {
    if (ci == n_classes) {
      // Assemble the full per-index Kramers replacement.
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
      auto term = append_kramers(prod.clone(), repl);
      term = expand_antisymm(term, /*skip_spinsymm*/ false);
      expand(term);
      rapid_simplify(term);
      canonicalize(term);
      rapid_simplify(term);
      if (mult != 1.0) {
        term = ex<Constant>(rational{static_cast<std::intmax_t>(mult)}) * term;
        expand(term);
        rapid_simplify(term);
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

  // Final cross-term canonicalize so equivalent multiset terms fold.
  ExprPtr r = result;
  expand(r);
  rapid_simplify(r);
  canonicalize(r);
  rapid_simplify(r);
  return r;
}

// =====================================================================
// Antisymmetrizer recombination (Phase 2.11). Inverse of expand_antisymm.
// See spinor.hpp for the algorithm.
// =====================================================================

namespace {

/// A parsed summand: total scalar, wrapper kind (0 raw / 1 Re / 2 Im), and
/// the tensor factors. The inner Product's scalar is folded into `scalar`.
struct RecombEntry {
  Constant::scalar_type scalar{1};
  int wrapper = 0;
  container::svector<ExprPtr> tensors;  // each is a Tensor ExprPtr
  bool alive = true;
};

/// Decompose a Sum summand into a RecombEntry. Recognized shapes:
///   Tensor | Product(s; Tensor...) | Re/Im(inner) | Product(s; Re/Im(inner))
/// where inner = Tensor | Product(s; Tensor...).
bool parse_recomb_entry(const ExprPtr& summand, RecombEntry& e) {
  e = RecombEntry{};
  ExprPtr cur = summand;
  // optional outer `Constant * (Re|Im wrapper)`
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
/// Rank-4 (2+2) only for now. Both indices compared by full Index identity
/// (recombinable terms share dummy names post-canonicalize).
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

/// Re-tag a Tensor with Symmetry::Antisymm (keeping label/indices/other
/// symmetries). Used to materialize a recombined antisymmetrized tensor.
ExprPtr make_antisymm(const Tensor& t) {
  return ex<Tensor>(
      t.label(), bra(container::svector<Index>{t.bra().begin(), t.bra().end()}),
      ket(container::svector<Index>{t.ket().begin(), t.ket().end()}),
      aux(container::svector<Index>{t.aux().begin(), t.aux().end()}),
      Symmetry::Antisymm, t.braket_symmetry(), t.column_symmetry());
}

/// Try to recombine two entries into one antisymmetrized entry.
/// Succeeds when, for some tensor position p: all other positions match
/// structurally, t2[p] is a single column-swap of t1[p], and the scalars
/// satisfy B.scalar == -A.scalar (single transposition → odd → sign -1).
/// Result: A.scalar * (... A.tensors[p] tagged Antisymm ...).
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
    const int swap = detect_single_column_swap(A.tensors[p]->as<Tensor>(),
                                               B.tensors[p]->as<Tensor>());
    if (swap == 0) continue;
    RecombEntry R = A;
    R.tensors[p] = make_antisymm(A.tensors[p]->as<Tensor>());
    return R;
  }
  return std::nullopt;
}

}  // namespace

ExprPtr antisymm_recombine(const ExprPtr& expr) {
  if (!expr->is<Sum>()) return expr;
  auto const& sum = expr->as<Sum>();

  std::vector<RecombEntry> entries;
  entries.reserve(sum.size());
  for (auto&& s : sum) {
    RecombEntry e;
    if (!parse_recomb_entry(s, e)) return expr;  // unrecognized — bail safely
    entries.push_back(std::move(e));
  }

  // Iterate to a fixed point: each pass merges one recombinable pair. A
  // second pass over the once-recombined terms catches doubly-antisymm
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
// Burnside-enumeration Kramers tracer (Phase 2.12 prototype).
// See spinor.hpp for the algorithm summary.
// =====================================================================

namespace {

/// Collect the Tensor factors of a Product (or wrap a bare Tensor).
container::svector<const Tensor*> product_tensors(const Product& p) {
  container::svector<const Tensor*> out;
  for (auto&& f : p) {
    if (f->is<Tensor>()) out.push_back(&f->as<Tensor>());
  }
  return out;
}

/// True iff indices @p p and @p q both sit in the same antisymmetric
/// bra OR ket group of Tensor @p t (i.e. t is Antisymm and {p,q} ⊆ bra
/// or {p,q} ⊆ ket).
bool same_antisymm_group(const Tensor& t, const Index& p, const Index& q) {
  if (t.symmetry() != Symmetry::Antisymm) return false;
  auto in = [](const auto& range, const Index& x) {
    return std::find(range.begin(), range.end(), x) != range.end();
  };
  const bool p_bra = in(t.bra(), p), q_bra = in(t.bra(), q);
  const bool p_ket = in(t.ket(), p), q_ket = in(t.ket(), q);
  return (p_bra && q_bra) || (p_ket && q_ket);
}

/// True iff index @p p appears anywhere in Tensor @p t.
bool tensor_has_index(const Tensor& t, const Index& p) {
  for (auto&& idx : t.const_braket())
    if (idx == p) return true;
  return false;
}

}  // namespace

ExprPtr kramers_trace_burnside(const ExprPtr& expr) {
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& t : *expr) s->append(kramers_trace_burnside(t));
    ExprPtr r = s;
    expand(r);
    rapid_simplify(r);
    canonicalize(r);
    rapid_simplify(r);
    return r;
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

  // --- 1. Find transposition generators of the index-permutation
  //        symmetry group. (p,q) qualifies iff in every tensor that
  //        contains both they share an antisymmetric bra/ket group, they
  //        co-occur in at least one tensor, and the sign is +1.
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
          // appears in only one slot of this tensor — a transposition
          // would move it to a different tensor position; not a symmetry
          ok = false;
          break;
        }
      }
      if (ok && n_common > 0 && (n_common % 2 == 0)) {
        generators.emplace_back(a, b);
      }
    }
  }

  // --- 2. Build the permutation group as a BFS closure over generators.
  //        Each element is a permutation of {0..n-1} (svector<size_t>).
  auto apply_transposition =
      [](container::svector<std::size_t> perm, std::size_t a,
         std::size_t b) -> container::svector<std::size_t> {
    std::swap(perm[a], perm[b]);
    return perm;
  };
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
          auto ng = apply_transposition(g, a, b);
          if (seen.insert(ng).second) next.push_back(ng);
        }
      }
      frontier = std::move(next);
    }
  }

  // --- 3. Burnside orbit enumeration over the 2^n bit-string configs.
  //        A permutation acts: (perm·cfg) has bit perm[k] = bit k of cfg.
  auto apply_perm = [](std::uint64_t cfg,
                       const container::svector<std::size_t>& perm) {
    std::uint64_t out = 0;
    for (std::size_t k = 0; k < perm.size(); ++k)
      if ((cfg >> k) & 1u) out |= (std::uint64_t{1} << perm[k]);
    return out;
  };
  const std::uint64_t n_configs = std::uint64_t{1} << n;
  container::svector<bool> visited(n_configs, false);
  // each rep: (config, orbit_size)
  container::svector<std::pair<std::uint64_t, std::size_t>> reps;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    if (visited[cfg]) continue;
    container::set<std::uint64_t> orbit;
    for (auto const& g : group) orbit.insert(apply_perm(cfg, g));
    for (auto o : orbit) visited[o] = true;
    // representative = lex-min of the orbit
    const std::uint64_t rep = *orbit.begin();
    reps.emplace_back(rep, orbit.size());
  }

  // --- 4. Per orbit rep: substitute Kramers labels, expand_antisymm,
  //        canonicalize, scale by orbit size.
  auto result = std::make_shared<Sum>();
  for (auto const& [cfg, orbit_size] : reps) {
    container::map<Index, Index> repl;
    for (std::size_t k = 0; k < n; ++k) {
      const bool is_dn = (cfg >> k) & 1u;
      repl[indices[k]] =
          is_dn ? make_kramers_dn(indices[k]) : make_kramers_up(indices[k]);
    }
    auto term = append_kramers(expr, repl);
    term = expand_antisymm(term, /*skip_spinsymm*/ false);
    expand(term);
    rapid_simplify(term);
    canonicalize(term);
    rapid_simplify(term);
    if (orbit_size != 1) {
      term =
          ex<Constant>(rational{static_cast<std::intmax_t>(orbit_size)}) * term;
      expand(term);
      rapid_simplify(term);
    }
    result->append(term);
  }

  // --- 5. Final cross-term canonicalize.
  ExprPtr r = result;
  expand(r);
  rapid_simplify(r);
  canonicalize(r);
  rapid_simplify(r);
  return r;
}

}  // namespace sequant::mbpt
