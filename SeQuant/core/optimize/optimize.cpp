#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/optimize/cost_model.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/optimize/sum.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {

namespace {

index_to_extent_t default_idx_to_size() {
  return [](Index const& ix) { return ix.space().approximate_size(); };
}

/// Diagnostic (env SEQUANT_FACTORIZER_DEBUG): for the chosen factorization
/// \p result of a single term, log each intermediate's result footprint AS THE
/// COST MODEL SIZES IT (idx_to_extent + inner_pow), plus the peak (the value
/// the DenseSpaceTime objective minimizes), so one can see why the factorizer
/// accepted a given intermediate -- e.g. an under-sized multi-composite tensor.
/// Footprints in mega-elements; outer{...} lists each free outer index extent,
/// inner{Np:e,...} lists each CSV/PNO composite's proto-index count N and the
/// extent e the model assigns it.
void log_chosen_factorization(ExprPtr const& result,
                              OptimizeOptions const& opts) {
  if (!opts.idx_to_extent || !result) return;
  auto describe = [&](ExprPtr const& e) -> std::pair<double, std::string> {
    auto g = get_unique_indices(e);
    std::vector<Index> all;
    all.insert(all.end(), g.bra.begin(), g.bra.end());
    all.insert(all.end(), g.ket.begin(), g.ket.end());
    all.insert(all.end(), g.aux.begin(), g.aux.end());
    auto const tot = tot_indices(all);
    double const fp = opt::detail::inner_aware_volume(tot, opts.idx_to_extent,
                                                      opts.inner_pow);
    std::string s = "outer{";
    for (auto const& o : tot.outer)
      s += std::to_string(opts.idx_to_extent(o)) + ",";
    s += "} inner{";
    for (auto const& i : tot.inner)
      s += std::to_string(i.proto_indices().size()) +
           "p:" + std::to_string(opts.idx_to_extent(i)) + ",";
    s += "}";
    return {fp, std::move(s)};
  };
  double peak = 0.;
  std::string peak_s;
  std::function<void(ExprPtr const&)> rec = [&](ExprPtr const& e) {
    if (!e->is<Product>()) return;
    auto const desc = describe(e);
    std::clog << "[FACTORIZE-NODE] fp=" << desc.first / 1e6 << "Me "
              << desc.second << "\n";
    if (desc.first > peak) {
      peak = desc.first;
      peak_s = desc.second;
    }
    for (auto const& f : e->as<Product>().factors()) rec(f);
  };
  std::clog << "[FACTORIZE] -------- chosen tree --------\n";
  rec(result);
  std::clog << "[FACTORIZE] PEAK fp=" << peak / 1e6 << "Me " << peak_s << "\n";
  std::clog.flush();
}

/// Optimize a Product that contains only Tensor and scalar factors.
/// Guards every access to OptimizeOptions::term_batch_axes: optimize_impl
/// optimizes a Sum's summands in parallel (sequant::for_each) and each summand
/// inserts into, re-keys and reads that one shared unordered_map; an unlocked
/// insert racing a rehash is a segfault (seen on a rank of an 8-rank run).
static std::mutex term_batch_axes_mutex;

ExprPtr opt_pure_product(Product const& prod, OptimizeOptions const& opts) {
  bool const subnet_cse = opts.CSE.subnet;
  // Build the cost knobs field-by-field from OptimizeOptions / its BatchPolicy.
  // Batching config (both role predicates, batch_target_size, inner_pow,
  // batch_persistent_only) travels on CostParams.
  CostParams cost;
  cost.is_volatile_leaf = opts.batch_policy.is_volatile_leaf;
  cost.volatile_weight = opts.volatile_weight;
  cost.footprint_weight = opts.footprint_weight;
  cost.peak_flops_tolerance = opts.peak_flops_tolerance;
  cost.roofline = opts.roofline;
  cost.accumulation_factor = opts.batch_policy.accumulation_factor;
  cost.peak_threshold = opts.batch_policy.peak_threshold;
  cost.prune_outer_products = opts.prune_outer_products;
  cost.batch_spectator_indices = opts.batch_policy.batch_spectator_indices;
  cost.is_batchable_contracted_index =
      opts.batch_policy.is_batchable_contracted_index;
  cost.is_batchable_external_index =
      opts.batch_policy.is_batchable_external_index;
  cost.batch_target_size = opts.batch_policy.batch_target_size;
  cost.inner_pow = opts.inner_pow;
  cost.batch_persistent_only = opts.batch_policy.persistent_only;
  // Filled by either batched arm below (both pass &node_axes as out_axes when
  // term_batch_axes is set); every other objective leaves it empty, so the
  // term_batch_axes insertion at the end is then a no-op-shaped empty-vector
  // entry (harmless: the binarizer only consumes entries for summands a
  // batched objective annotated).
  container::vector<NodeBatchAnnotation> node_axes;
  auto run = [&]() -> ExprPtr {
    if (opts.objective_function == ObjectiveFunction::DenseFLOPs)
      return opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSize)
      return opt::single_term_opt<ObjectiveFunction::DenseSize>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTime)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTime>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseTimeSpace)
      return opt::single_term_opt<ObjectiveFunction::DenseTimeSpace>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTimeBatched)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTimeBatched>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.term_batch_axes ? &node_axes : nullptr);
    SEQUANT_ASSERT(opts.objective_function ==
                   ObjectiveFunction::DenseTimeSpaceBatched);
    return opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
        prod, opts.idx_to_extent, subnet_cse, cost,
        opts.term_batch_axes ? &node_axes : nullptr);
  };
  ExprPtr result = run();
  if (opts.term_batch_axes) {
    // optimize_impl optimizes a Sum's summands with sequant::for_each
    // (std::execution::par_unseq), so opt_pure_product runs concurrently across
    // summands -- and this insert into the shared term_batch_axes map is not
    // thread-safe (std::unordered_map: concurrent inserts race even on distinct
    // keys -- a rehash tears the structure). Serialize just the insert; the
    // heavy DP above stays parallel. Without this the map is corrupted and the
    // downstream whole-Sum re-key reads a wrong-sized node_batch_axes, tripping
    // binarize's node_counter == size assertion (a nondeterministic, thread-
    // count-dependent SIGABRT, absent under a sequential par_unseq fallback
    // such as libc++).
    std::lock_guard<std::mutex> lock(term_batch_axes_mutex);
    (*opts.term_batch_axes)[result.get()] = std::move(node_axes);
  }
  if (std::getenv("SEQUANT_FACTORIZER_DEBUG"))
    log_chosen_factorization(result, opts);
  return result;
}

/// Deliberately non-identifier label prefix used to stand in for non-Tensor,
/// non-scalar factors during single-term optimization. Chosen so that no
/// user-defined tensor label can collide with it.
inline constexpr std::wstring_view placeholder_label_prefix = L"@__opt_";

/// The non_tensors slot a placeholder tensor label (see
/// placeholder_label_prefix) stands in for. The prefix is internal; anything
/// carrying it must have been emitted by opt_mixed_product with a
/// pure-decimal suffix, so any deviation is a programming error.
std::size_t placeholder_index(std::wstring_view label) {
  SEQUANT_ASSERT(label.starts_with(placeholder_label_prefix));
  auto suffix_view = label.substr(placeholder_label_prefix.size());
  SEQUANT_ASSERT(!suffix_view.empty());
  std::size_t suffix = 0;
  for (wchar_t c : suffix_view) {
    SEQUANT_ASSERT(c >= L'0' && c <= L'9');
    suffix = suffix * 10 + static_cast<std::size_t>(c - L'0');
  }
  return suffix;
}

/// Number of contraction (DP) nodes binarize builds for \p e, i.e. how many
/// BinarizationOptions::node_batch_axes entries it consumes, mirroring
/// impl::binarize: a Product consumes its factors' counts (in factor order)
/// plus one per tensor x tensor left-fold step (#non-scalar factors - 1); a
/// Sum factor is binarized on a private counter (consumes none); a Re/Im
/// wrapper factor shares the counter only when every factor is a scalar (its
/// inner nodes are then the product's DP nodes).
std::size_t binarize_dp_node_count(ExprPtr const& e) {
  if (e->is<RealPart>())
    return binarize_dp_node_count(e->as<RealPart>().inner());
  if (e->is<ImagPart>())
    return binarize_dp_node_count(e->as<ImagPart>().inner());
  if (!e->is<Product>()) return 0;
  auto const& prod = e->as<Product>();
  bool const wrapper_shares_counter = ranges::all_of(
      prod.factors(), [](ExprPtr const& f) { return f->is_scalar(); });
  std::size_t n = 0;
  std::size_t n_tensor = 0;
  for (auto const& f : prod.factors()) {
    if (f->is<Sum>()) {
      // private counter, no entries
    } else if (f->is<RealPart>() || f->is<ImagPart>()) {
      if (wrapper_shares_counter) n += binarize_dp_node_count(f);
    } else {
      n += binarize_dp_node_count(f);
    }
    if (!f->is_scalar()) ++n_tensor;
  }
  return n + (n_tensor > 1 ? n_tensor - 1 : 0);
}

/// Optimize a Product that contains some non-Tensor, non-scalar factors by
/// substituting placeholder tensors with target indices, optimizing the
/// resulting tensor-only product, then swapping the originals back in.
ExprPtr optimize_impl(ExprPtr const& expr, OptimizeOptions const& opts,
                      bool reorder, bool parallel_outer);

ExprPtr opt_mixed_product(Product const& prod, OptimizeOptions const& opts) {
  container::svector<ExprPtr> non_tensors(prod.size());
  container::svector<ExprPtr> new_factors;
  new_factors.reserve(prod.size());

  for (std::size_t i = 0; i < prod.size(); ++i) {
    auto&& f = prod.factor(i);
    if (f->is<Tensor>() || f->is_scalar()) {
      new_factors.emplace_back(f);
    } else {
      // A non-tensor factor (a Sum of products, e.g. the flavor bracket a
      // CSV transform wraps around a projected leaf, sum_flavors g.C.C; or a
      // nested product) is opaque to the outer contraction order, but its
      // own contraction order matters just as much: put back as written it
      // evaluates in its authored left-to-right order. Measured on DCH
      // cc-pVDZ PNS-CCD (2026-09-05): a projection bracket whose external-
      // pair C came first materialized an n_occ^4 n_v n_csv intermediate
      // (3.3 GB each, 32 GB of them cached) where the optimal order, which
      // the DF cost model had assumed, peaks at n_occ^2 n_v n_csv.
      non_tensors[i] = optimize_impl(f, opts, /*reorder=*/false,
                                     /*parallel_outer=*/false);
      auto target_idxs = get_unique_indices(f);
      new_factors.emplace_back(ex<Tensor>(
          std::wstring(placeholder_label_prefix) + std::to_wstring(i),
          bra(target_idxs.bra), ket(target_idxs.ket), aux(target_idxs.aux)));
    }
  }

  auto result = opt_pure_product(
      Product{prod.scalar(), new_factors, Product::Flatten::No}, opts);

  // Per-node batch annotations (opts.term_batch_axes): opt_pure_product keyed
  // the outer network's entries -- one per DP node over the placeholders, in
  // binarize's left-first post-order -- on `result`, and each nested Product
  // factor's own optimization above keyed its entries on non_tensors[i].
  // binarize consumes one shared counter in post-order over the whole tree,
  // a nested Product factor's contraction nodes included (only a Sum factor
  // gets a private counter and no entries), so splice each nested product's
  // entries in at its placeholder's position and re-key the merged list on
  // `result`. Without this the outer entries land on the brackets' inner
  // nodes: e.g. the DF driver (g C C)(K) . (g C C)(K), contracted over the
  // aux index K at the root with each bracket keeping K open, had the root's
  // contracted-K mark stamped on the first bracket's inner node, whose result
  // still carries K -- the batched runtime then accumulated per-batch
  // partials of unequal K extent (Kramers-union PNS-CCD, 2026-09-09).
  if (opts.term_batch_axes) {
    std::lock_guard<std::mutex> lock(term_batch_axes_mutex);
    container::vector<NodeBatchAnnotation> outer;
    if (auto it = opts.term_batch_axes->find(result.get());
        it != opts.term_batch_axes->end())
      outer = std::move(it->second);
    container::vector<NodeBatchAnnotation> merged;
    std::size_t next_outer = 0;
    std::function<void(ExprPtr const&)> walk = [&](ExprPtr const& e) {
      if (e->is<Product>()) {
        std::size_t n_tensor = 0;
        for (auto const& f : e->as<Product>().factors()) {
          walk(f);
          if (!f->is_scalar()) ++n_tensor;
        }
        for (std::size_t k = 1; k < n_tensor; ++k, ++next_outer)
          merged.push_back(next_outer < outer.size() ? outer[next_outer]
                                                     : NodeBatchAnnotation{});
        return;
      }
      if (!e->is<Tensor>()) return;
      auto const label = e->as<Tensor>().label();
      if (!label.starts_with(placeholder_label_prefix)) return;
      auto const& inner = non_tensors[placeholder_index(label)];
      SEQUANT_ASSERT(inner);
      // the nested product's own entries, in its post-order; a Sum bracket
      // contributes none (private counter in binarize). A count mismatch
      // (an inner optimization path that did not record) degrades to
      // unannotated inner nodes rather than misaligning the outer ones.
      std::size_t const need = binarize_dp_node_count(inner);
      auto it = opts.term_batch_axes->find(inner.get());
      if (it != opts.term_batch_axes->end() && it->second.size() == need)
        merged.insert(merged.end(), it->second.begin(), it->second.end());
      else
        merged.insert(merged.end(), need, NodeBatchAnnotation{});
    };
    walk(result);
    (*opts.term_batch_axes)[result.get()] = std::move(merged);
  }

  auto replacer = [&non_tensors](ExprPtr& out) {
    if (!out->is<Tensor>()) return;
    auto label = out->as<Tensor>().label();
    if (!label.starts_with(placeholder_label_prefix)) return;
    auto const suffix = placeholder_index(label);
    SEQUANT_ASSERT(suffix < non_tensors.size() && non_tensors[suffix]);
    out = non_tensors[suffix].clone();
  };

  result->visit(replacer, /* atoms_only = */ true);
  return result;
}

/// Recursive workhorse. \p parallel_outer controls whether the (single)
/// outermost Sum's summands are processed in parallel; nested recursive
/// calls always run sequentially to avoid `sequant::for_each` oversubscription.
ExprPtr optimize_impl(ExprPtr const& expr, OptimizeOptions const& opts,
                      bool reorder, bool parallel_outer) {
  // Re/Im wrappers are transparent to optimization: optimize the wrapped
  // expression and re-wrap. Without this the wrapper is returned untouched
  // and its inner product evaluates in naive left-to-right order (measured:
  // 14.4 GB vs 1.7 GB peak RSS on a Kramers-CSV MP2 energy whose TRS fold
  // wrapped three terms).
  // A wrapper at the summand root: its inner contraction nodes are the
  // summand's DP nodes (binarize shares the node counter with it), so its
  // batch axes are re-keyed under the wrapper the caller keys on.
  auto rekey_axes = [&opts](ExprPtr const& inner, ExprPtr const& wrapper) {
    if (!opts.term_batch_axes) return;
    std::lock_guard<std::mutex> lock(term_batch_axes_mutex);
    auto it = opts.term_batch_axes->find(inner.get());
    if (it != opts.term_batch_axes->end())
      (*opts.term_batch_axes)[wrapper.get()] = it->second;
  };
  if (expr->is<RealPart>()) {
    auto inner = optimize_impl(expr->as<RealPart>().inner(), opts,
                               /*reorder=*/false, /*parallel_outer=*/false);
    auto wrapped = ex<RealPart>(inner);
    rekey_axes(inner, wrapped);
    return wrapped;
  }
  if (expr->is<ImagPart>()) {
    auto inner = optimize_impl(expr->as<ImagPart>().inner(), opts,
                               /*reorder=*/false, /*parallel_outer=*/false);
    auto wrapped = ex<ImagPart>(inner);
    rekey_axes(inner, wrapped);
    return wrapped;
  }
  if (expr->is<Product>()) {
    auto const& prod_in = expr->as<Product>();
    // Re/Im wrapper factors are transparent too (the conjugate-pair fold
    // emits `2 Re[A]`): RealPart::is_scalar() would otherwise let the
    // wrapper pass through opt_pure_product as an opaque scalar with A left
    // in its naive left-to-right order. Optimize each wrapper's inner first.
    auto const has_wrapper = ranges::any_of(prod_in, [](auto&& x) {
      return x->template is<RealPart>() || x->template is<ImagPart>();
    });
    Product::factors_type factors;
    container::svector<ExprPtr> inners;  // optimized wrapper inners, in order
    if (has_wrapper) {
      for (auto const& f : prod_in) {
        if (f->is<RealPart>()) {
          inners.push_back(
              optimize_impl(f->as<RealPart>().inner(), opts, false, false));
          factors.push_back(ex<RealPart>(inners.back()));
        } else if (f->is<ImagPart>()) {
          inners.push_back(
              optimize_impl(f->as<ImagPart>().inner(), opts, false, false));
          factors.push_back(ex<ImagPart>(inners.back()));
        } else
          factors.push_back(f);
      }
    }
    Product const prod_rewrapped =
        has_wrapper ? Product{prod_in.scalar(), factors, Product::Flatten::No}
                    : Product{};
    auto const& prod = has_wrapper ? prod_rewrapped : prod_in;
    bool pure = ranges::all_of(prod, [](auto&& x) {
      return x->template is<Tensor>() || x->is_scalar();
    });
    auto result =
        pure ? opt_pure_product(prod, opts) : opt_mixed_product(prod, opts);
    // Wrappers whose siblings are all scalars (the fold's `2 Re[A]`): the
    // product has no contraction nodes of its own, so the summand's DP nodes
    // are exactly the wrappers' inner nodes, in factor order (binarize
    // shares its node counter with such wrappers). Re-key their batch axes
    // under the summand pointer the caller keys on.
    if (has_wrapper && opts.term_batch_axes) {
      const bool scalar_siblings =
          ranges::all_of(prod_in, [](auto&& x) { return x->is_scalar(); });
      if (scalar_siblings) {
        std::lock_guard<std::mutex> lock(term_batch_axes_mutex);
        container::vector<NodeBatchAnnotation> axes;
        for (auto const& inner : inners) {
          auto it = opts.term_batch_axes->find(inner.get());
          if (it != opts.term_batch_axes->end())
            axes.insert(axes.end(), it->second.begin(), it->second.end());
        }
        (*opts.term_batch_axes)[result.get()] = std::move(axes);
      }
    }
    return result;
  }

  if (expr->is<Sum>()) {
    auto const& in_sum = expr->as<Sum>();
    Sum::summands_type new_smands(in_sum.size());

    // Every summand is optimized on a private clone, taken here, sequentially,
    // before the (possibly parallel) loop below. Summands routinely share
    // subexpression objects -- tensors reused by expand(), or a whole nested
    // Sum factor (the flavor bracket a CSV transform wraps around a projected
    // leaf) reused across the terms it appears in -- and Index/Expr memoize
    // labels and hashes lazily in unsynchronized mutable members. Optimizing
    // (and, since opt_mixed_product also optimizes nested Sum factors, walking
    // and canonicalizing) a shared object from several threads races on those
    // caches and can yield a run-to-run different tree; on a distributed
    // evaluation that is a deadlock, since every rank must build the same
    // tree (DCH PNS-MP1 on 8 ranks, 2026-09-05: two ranks built a different
    // residual tree and the run hung in iteration 1). The clones make
    // invariant (1) below hold by construction; the input is never touched
    // concurrently.
    Sum::summands_type private_smands;
    private_smands.reserve(in_sum.size());
    for (auto const& s : in_sum.summands())
      private_smands.push_back(s->clone());

    auto do_term = [&](std::size_t i) {
      new_smands[i] = optimize_impl(private_smands[i], opts,
                                    /*reorder=*/false,
                                    /*parallel_outer=*/false);
    };

    // Thread-safety of the parallel branch rests on two invariants; do NOT
    // break them without re-auditing:
    //   1. Each task writes a distinct, pre-allocated new_smands[i] slot, and
    //      the work below (single_term_opt ->
    //      TensorNetwork::canonicalize_slots) operates on per-task *clones* of
    //      the input tensors. The lazily populated `mutable` caches on
    //      Expr/Index (hash_value_, label_, ...) are unsynchronized, so
    //      concurrent work must never read/write them on a shared (non-cloned)
    //      node. Index comparison touches only immutable members, so building
    //      index sets over shared indices is safe.
    //   2. The binarize() pass below DOES read Index::label() (a lazy cache
    //      write) on the optimized summands, so it is run *sequentially, after*
    //      for_each() has joined -- never inside do_term().
    // The default Context and cardinal_tensor_labels must also be configured
    // before entering here (their writes are unsynchronized unless
    // SEQUANT_CONTEXT_MANIPULATION_THREADSAFE); optimize() only reads them.
    if (parallel_outer && in_sum.size() > 1) {
      auto indices = ranges::views::iota(std::size_t{0}, in_sum.size());
      sequant::for_each(indices, do_term);
    } else {
      for (std::size_t i = 0; i < in_sum.size(); ++i) do_term(i);
    }

    Sum new_sum(std::move(new_smands), Sum::move_only_tag{});

    // Re-key the per-summand batch annotations onto the final reassembled Sum.
    // opt_pure_product keyed each summand's node_batch_axes (one entry per
    // contraction node, left-first post-order) by that optimized summand's
    // Product pointer. But the caller binarizes the whole reassembled Sum in
    // one call and looks the annotation up by the final Sum pointer -- and
    // under reorder, opt::reorder's clone-on-append (Sum::append clones) gives
    // the final summands new pointers while new_sum (which still holds the
    // keyed pointers) is destroyed on return. So gather the per-summand vectors
    // in the final summand order into one whole-tree vector -- binarize walks
    // the Sum-tree in that same order, one entry per contraction node, so the
    // flat node_batch_axes stays aligned with its node counter -- and store it
    // under the final Sum pointer, dropping the now-unreachable per-summand
    // entries. Without this, every batch annotation is silently lost and
    // over-budget intermediates materialize whole. `order` is a list of
    // clusters, each a list of positions into new_sum, flattened in emission
    // order (identity for the no-reorder path); it must match how the final Sum
    // orders its summands.
    // The keys are the optimized summands' addresses, taken here, before the
    // no-reorder path below moves new_sum into its result (the ExprPtrs keep
    // their pointees, so the keys stay valid; new_sum's summand list does not).
    container::vector<Expr const*> smand_keys;
    smand_keys.reserve(new_sum.size());
    for (auto const& s : new_sum.summands()) smand_keys.push_back(s.get());
    auto rekey_onto =
        [&](ExprPtr const& result,
            container::vector<container::vector<std::size_t>> const& order) {
          if (!opts.term_batch_axes) return;
          // This re-key runs for every Sum, including a nested one (a
          // flavor-sum bracket of a Kramers-split residual) whose enclosing
          // summand is being optimized in parallel with its siblings, which
          // insert into the same map: the find/erase/insert below must hold
          // the map's mutex like every other access, or a concurrent insert
          // (rehash) tears the map and entries go missing -- binarize's
          // node_counter == node_batch_axes.size() then fails at random.
          std::lock_guard<std::mutex> lock(term_batch_axes_mutex);
          container::vector<NodeBatchAnnotation> combined;
          // SEQUANT_BATCH_AXES_DEBUG=1: one line per summand (entries found
          // or missing), to align with binarize's per-summand node counts
          static const bool debug = std::getenv("SEQUANT_BATCH_AXES_DEBUG");
          for (auto const& clstr : order)
            for (auto p : clstr) {
              auto it = opts.term_batch_axes->find(smand_keys.at(p));
              if (debug) {
                // new_sum may already be moved-from here; the keyed pointees
                // stay alive (see above)
                Expr const& s = *smand_keys.at(p);
                std::string spelling = toUtf8(s.to_latex());
                std::cerr << "[batch-axes] optimizer summand " << p << ": "
                          << (it == opts.term_batch_axes->end()
                                  ? std::string("NO ENTRY")
                                  : std::to_string(it->second.size()) +
                                        " entries")
                          << " type="
                          << (s.is<Product>() ? "Product"
                              : s.is<Sum>()   ? "Sum"
                                              : "other")
                          << " | " << spelling.substr(0, 160) << "\n";
              }
              if (it == opts.term_batch_axes->end()) continue;
              combined.insert(combined.end(),
                              std::make_move_iterator(it->second.begin()),
                              std::make_move_iterator(it->second.end()));
              opts.term_batch_axes->erase(it);
            }
          (*opts.term_batch_axes)[result.get()] = std::move(combined);
        };

    if (!reorder) {
      container::vector<container::vector<std::size_t>> identity;
      identity.reserve(new_sum.size());
      for (std::size_t i = 0; i < new_sum.size(); ++i) identity.push_back({i});
      auto result = ex<Sum>(std::move(new_sum));
      rekey_onto(result, identity);
      return result;
    }

    // Binarize once per optimized summand and hand the nodes to reorder()
    // so they aren't re-built inside clusters(). NOTE: this runs sequentially
    // by design -- see invariant (2) above.
    container::vector<FullBinaryNode<EvalExpr>> nodes;
    nodes.reserve(new_sum.size());
    // per-summand binarize for ordering only; positional head doesn't escape.
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    for (auto const& s : new_sum.summands()) nodes.push_back(binarize(s));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    // Same (new_sum, nodes) opt::reorder consumes, so the flattened cluster
    // order equals the final summand order the reordered Sum emits.
    auto const order = opt::clusters(new_sum, nodes);
    auto result = ex<Sum>(opt::reorder(new_sum, nodes));
    rekey_onto(result, order);
    return result;
  }

  return expr->clone();
}

}  // namespace

ExprPtr optimize(ExprPtr const& expr, OptimizeOptions opts) {
  if (!opts.idx_to_extent) opts.idx_to_extent = default_idx_to_size();
  return optimize_impl(expr, opts, opts.reorder == ReorderSum::Reorder,
                       /*parallel_outer=*/true);
}

OptimizeResult optimize_result(ExprPtr const& expr, OptimizeOptions opts) {
  if (!opts.idx_to_extent) opts.idx_to_extent = default_idx_to_size();
  OptimizeResult res;
  res.expr = optimize_impl(expr, opts, opts.reorder == ReorderSum::Reorder,
                           /*parallel_outer=*/true);
  return res;
}

ResultExpr& optimize(ResultExpr& expr, OptimizeOptions opts) {
  expr.expression() = optimize(expr.expression(), std::move(opts));
  return expr;
}

ResultExpr& optimize(ResultExpr&& expr, OptimizeOptions opts) {
  return optimize(expr, std::move(opts));
}

// backwards compatibility overloads

namespace {
inline OptimizeOptions compatibility_opts(bool reorder_sum) {
  return OptimizeOptions{
      .reorder = reorder_sum ? ReorderSum::Reorder : ReorderSum::NoReorder,
      .inner_pow = {}};
}
}  // namespace

ExprPtr optimize(ExprPtr const& expr, bool reorder_sum) {
  return optimize(expr, compatibility_opts(reorder_sum));
}

ResultExpr& optimize(ResultExpr& expr, bool reorder_sum) {
  return optimize(expr, compatibility_opts(reorder_sum));
}

ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum) {
  return optimize(std::move(expr), compatibility_opts(reorder_sum));
}

}  // namespace sequant
