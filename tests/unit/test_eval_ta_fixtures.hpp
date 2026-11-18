//
// Shared prelude of the TiledArray-backend unit tests (test_eval_ta*.cpp):
// the includes, the random-array fixtures and the array comparison helpers.
//
// The tests are spread over several translation units because a single one
// needed more compiler memory than a CI runner has. Everything lives in a
// named namespace (not an anonymous one) with `inline` functions so that no
// translation unit warns about helpers it does not use.
//

#ifndef SEQUANT_TESTS_UNIT_TEST_EVAL_TA_FIXTURES_HPP
#define SEQUANT_TESTS_UNIT_TEST_EVAL_TA_FIXTURES_HPP

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"
#include "tensor_power_fixtures.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/backends/tiledarray/array_ops.hpp>
#include <SeQuant/core/eval/backends/tiledarray/eval_context.hpp>
#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <cstddef>
#include <utility>

#include <TiledArray/tensor/arena_kernels.h>
#include <TiledArray/tensor/arena_tensor.h>
#include <tiledarray.h>
#include <boost/regex.hpp>

#include <data_info.hpp>
#include <eval_utils.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/intersperse.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace TiledArray {
template <typename>
constexpr bool is_tnsr_expr_v{};

template <typename Arg, auto... Args>
constexpr bool is_tnsr_expr_v<expressions::TsrExpr<Arg, Args...>> = true;

template <typename T>
concept tnsr_expr = is_tnsr_expr_v<T>;

template <typename T>
concept array = TA::detail::array_tos<T> || TA::detail::array_tot<T>;

}  // namespace TiledArray

namespace sequant::tests::eval_ta {

///
/// \brief Represents the outer indices and the inner indices of a nested
/// tensor.
///
/// \note The nested tensor is a concept that generalizes the sequant::Tensor
/// with and without proto indices. sequant::Tensors with proto indices have
/// outer and inner indices, whereas, those without proto indices only have
/// outer indices.
///
struct NestedTensorIndices {
  sequant::container::svector<sequant::Index> outer, inner;

  explicit NestedTensorIndices(sequant::Tensor const& tnsr) {
    using ranges::views::join;
    using ranges::views::transform;
    using namespace sequant;

    auto append_unique = [](auto& cont, auto const& el) {
      if (!ranges::contains(cont, el)) cont.emplace_back(el);
    };

    for (Index const& ix : tnsr.const_braket_indices()) {
      append_unique(ix.has_proto_indices() ? inner : outer, ix);
    }

    for (Index const& ix :
         tnsr.const_braket_indices() | transform(&Index::proto_indices) | join)
      append_unique(outer, ix);

    for (auto&& ix : tnsr.aux()) {
      SEQUANT_ASSERT(!ix.has_proto_indices() &&
                     "Aux indices with proto indices not supported");
      outer.emplace_back(ix);
    }
  }

  [[nodiscard]] auto outer_inner() const noexcept {
    return ranges::views::concat(outer, inner);
  }
};

inline auto to_ta_node(sequant::FullBinaryNode<sequant::EvalExpr> node) {
  using namespace sequant;
  // The array-mode list of every node is its canon_indices() (outer modes
  // first, then the proto-carrying inner ones), as produced by the core
  // canonicalizer; EvalExprTA derives its annotation from that list. No
  // per-node re-derivation from the tensor's slots.
  return transform_node(node, [](auto&& val) { return EvalExprTA(val); });
}

inline auto eval_node(sequant::ExprPtr const& expr) {
  // sequant::binarize(ExprPtr) is deprecated for caller-visible head
  // construction; this helper exists for legacy test sections that don't
  // depend on the head's bra/ket layout.
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return to_ta_node(sequant::binarize(expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

inline auto eval_node(sequant::ResultExpr const& res) {
  return to_ta_node(sequant::binarize(res));
}

inline auto tensor_to_key(sequant::Tensor const& tnsr) {
  // One random tensor per tensor LABEL and index STRUCTURE (spaces, spins,
  // proto pattern), independent of the index names: that is the DAG's own
  // value identity (a CSE-shared leaf is one array), so forest references
  // and the ordered executor agree. The PAO-like `m` space is normalized too.
  static auto const idx_rgx = boost::wregex{L"([iaxm])([↑↓])?(_?\\d+)"};
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i"   ? L"o"
            : mo[1].str() == L"a" ? L"v"
            : mo[1].str() == L"m" ? L"m"
                                  : L"x") +
           mo[2].str();
  };

  NestedTensorIndices oixs{tnsr};
  if (oixs.inner.empty()) {
    auto const tnsr_deparsed =
        sequant::serialize(tnsr.clone(), {.annot_symm = false});
    return boost::regex_replace(tnsr_deparsed, idx_rgx, formatter);
  } else {
    using ranges::views::intersperse;
    using ranges::views::join;
    using ranges::views::transform;
    using namespace sequant;

    auto ix_lbl = [&formatter](Index const& ix) -> std::wstring {
      std::wstring lbl(ix.label().data());
      return boost::regex_replace(lbl, idx_rgx, formatter);
    };

    auto ixs_lbl = [&ix_lbl](auto const& ixs) -> std::wstring {
      return ixs | transform(ix_lbl) | intersperse(L",") | join |
             ranges::to<std::wstring>;
    };

    std::wstring result(tnsr.label());
    result += L"{" + ixs_lbl(oixs.outer) + L";" + ixs_lbl(oixs.inner) + L"}";
    return result;
  }
}

inline auto tensor_to_key(std::wstring_view spec) {
  return tensor_to_key(sequant::deserialize<sequant::ExprPtr>(
                           spec, {.def_perm_symm = sequant::Symmetry::Nonsymm})
                           ->as<sequant::Tensor>());
}

template <typename NumericT>
auto random_tensor(TA::Range const& rng) {
  TA::Tensor<NumericT> result{rng};
  std::generate(result.begin(), result.end(),
                TA::detail::MakeRandom<NumericT>::generate_value);
  return result;
}

/// The inner tile of every nested (ToT) array in these tests: the arena-pinned
/// view a CSV/PNO array stores (the inner tile SeQuant-eval-ta instantiates its
/// nested Result for), not an owning TA::Tensor<T>. An arena cell is a view
/// into its outer tile's slab, so a nested tile is built in one pass through
/// make_tot_tile() instead of by assigning inner tensors into an outer tile.
template <typename NumericT>
using tot_inner_t = TA::ArenaTensor<NumericT>;
template <typename NumericT>
using tot_tile_t = TA::Tensor<tot_inner_t<NumericT>>;

/// A nested tile over \p outer_rng whose inner cell at outer ordinal \c o spans
/// \c inner_rng(o) (a TA::Range; zero volume leaves the cell empty) and holds
/// \c fill(o, k) at its ordinal \c k. An arena cell's range is zero-based, so
/// only the extents of \c inner_rng(o) carry over.
template <typename NumericT, typename InnerRangeFn, typename FillFn>
tot_tile_t<NumericT> make_tot_tile(TA::Range const& outer_rng,
                                   InnerRangeFn&& inner_rng, FillFn&& fill) {
  using inner_range_t = typename tot_inner_t<NumericT>::range_type;
  return TA::detail::make_nested_tile<tot_tile_t<NumericT>>(
      outer_rng,
      [&](auto const& idx) {
        // extent() is a view into the range, which must outlive it; the arena
        // range type is built from an extent vector, not from a TA::Range.
        TA::Range const r = inner_rng(outer_rng.ordinal(idx));
        auto const ext = r.extent();
        return inner_range_t(std::vector<std::size_t>(ext.begin(), ext.end()));
      },
      [&](auto& cell, auto const& idx) {
        auto const o = outer_rng.ordinal(idx);
        for (std::size_t k = 0; k < cell.size(); ++k)
          cell.data()[k] = fill(o, k);
      });
}

// note: all the inner tensors (elements of the outer tensor)
//       have the same @c inner_rng
template <typename NumericT>
auto random_tensor_of_tensor(TA::Range const& outer_rng,
                             TA::Range const& inner_rng) {
  return make_tot_tile<NumericT>(
      outer_rng, [&inner_rng](std::size_t) { return inner_rng; },
      [](std::size_t, std::size_t) {
        return TA::detail::MakeRandom<NumericT>::generate_value();
      });
}

template <typename NumericT = double, typename TAPolicyT = TA::SparsePolicy>
class rand_tensor_yield {
  TA::World& world;
  size_t nocc_;
  size_t nvirt_;
  size_t naux_;
  // max tile size along each mode; the default (~0) makes a single tile per
  // mode (the original behavior). Set smaller to produce multi-tile arrays.
  size_t max_tile_ = ~size_t{0};
  mutable sequant::container::map<std::wstring, sequant::ResultPtr>
      label_to_er_;
  /// the TA annotation (layout) of the array stored under a tensor key
  mutable sequant::container::map<std::wstring, std::string> layout_of_;

 public:
  /// Produce arrays whose modes are tiled in blocks of at most \p n.
  /// \p n must be positive (0 would make the tiling loop in make_tr1 spin).
  void set_max_tile(size_t n) {
    REQUIRE(n > 0);
    max_tile_ = n;
  }

  /// Forget the cached tensor for label \p lbl so the next yield regenerates a
  /// fresh (different) random value -- models an amplitude tensor changing
  /// between CC iterations (a "volatile" leaf), while every other leaf stays
  /// fixed. Used to test cross-iteration cache correctness.
  void reset_label(std::wstring const& lbl) {
    label_to_er_.erase(lbl);
    layout_of_.erase(lbl);
  }

  /// Serve \p r for every occurrence of tensor \p tnsr (by tensor_to_key);
  /// \p r must be in the tensor's slot layout. For tests that need leaf data
  /// they choose: same-shaped random arrays may hold identical values.
  void set_leaf(sequant::Tensor const& tnsr, sequant::ResultPtr r) {
    label_to_er_.insert_or_assign(tensor_to_key(tnsr), std::move(r));
  }

  using array_type = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;

  /// Serve f(slot ordinals of \p tnsr) for every occurrence of the leaf
  /// \p tnsr, in either orientation: a leaf node of a bra-ket symmetric
  /// tensor may hold either. Occupied (i) modes have nocc elements, all others
  /// nvirt.
  template <typename F>
  void pin(sequant::Tensor const& tnsr, F f) {
    using namespace sequant;
    auto const nb = tnsr.bra_rank();
    auto pin_as = [&](Tensor const& as, bool swapped) {
      std::vector<TA::TiledRange1> tr1s;
      for (Index const& ix : as.const_braket_indices())
        tr1s.push_back(
            tr1_for(L"", ix.space().base_key() == L"i" ? nocc_ : nvirt_));
      array_type arr(world, TA::TiledRange(tr1s.begin(), tr1s.end()));
      arr.init_tiles([f, nb, swapped](TA::Range const& r) {
        TA::Tensor<NumericT> tile(r);
        for (auto const& idx : r) {
          std::vector<std::size_t> x(idx.begin(), idx.end());
          if (swapped) std::rotate(x.begin(), x.end() - nb, x.end());
          tile(idx) = f(x);
        }
        return tile;
      });
      // f may refer to the caller's locals, and nothing else waits for the
      // tiles of a leaf that no evaluation reads
      world.gop.fence();
      set_leaf(as, eval_result<ResultTensorTA<array_type>>(arr));
    };
    pin_as(tnsr, false);
    pin_as(
        Tensor(tnsr.label(), bra(tnsr.ket()), ket(tnsr.bra()), aux(tnsr.aux()),
               TensorSymmetries{.perm = tnsr.symmetry(),
                                .braket = tnsr.braket_symmetry(),
                                .column = tnsr.column_symmetry()}),
        true);
  }

  using array_tot_type = TA::DistArray<tot_tile_t<NumericT>, TAPolicyT>;
  using numeric_type = NumericT;

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt}, naux_{nvirt * 2} {}

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt, size_t naux)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt}, naux_{naux} {}

  /// Backend array-ops (zero destination + axis chunking) for a batched eval
  /// driven by this yield, built from THIS yield's OWN per-space tranges (i ->
  /// nocc, a/m -> nvirt, x -> naux, tiled in blocks of max_tile_). Nothing is
  /// read out of the eval DAG; the leaf source owns the tilings, exactly as
  /// mpqc's registries do in production. \tparam ToTArray the nested array type
  /// when a scatter destination is nested (defaults to flat).
  /// Per-space tile size override (SEQUANT_UT_TILE_<space base key>, e.g.
  /// SEQUANT_UT_TILE_i=1 SEQUANT_UT_TILE_x=168); falls back to max_tile_.
  size_t tile_for(std::wstring const& bk) const {
    std::string key = "SEQUANT_UT_TILE_" + sequant::toUtf8(bk);
    if (char const* e = std::getenv(key.c_str()); e && std::atol(e) > 0)
      return static_cast<size_t>(std::atol(e));
    return max_tile_;
  }
  TA::TiledRange1 tr1_for(std::wstring const& bk, size_t e) const {
    size_t const t = tile_for(bk);
    sequant::container::svector<size_t> b;
    for (size_t x = 0; x < e; x += t) b.push_back(x);
    b.push_back(e);
    return TA::TiledRange1(b.begin(), b.end());
  }

  template <typename ToTArray = array_type>
  sequant::BackendArrayOps array_ops() const {
    auto ibr = sequant::get_default_context().index_basis_registry();
    std::map<std::wstring, TA::TiledRange1> m;
    auto const bk = [&](wchar_t const* l) {
      return std::wstring(ibr->retrieve_space(l).base_key());
    };
    m[bk(L"i")] = tr1_for(bk(L"i"), nocc_);
    m[bk(L"a")] = tr1_for(bk(L"a"), nvirt_);
    m[bk(L"x")] = tr1_for(bk(L"x"), naux_);
    m[bk(L"m")] = tr1_for(bk(L"m"), nvirt_);
    return sequant::make_ta_array_ops<array_type, ToTArray>(std::move(m),
                                                            world);
  }

  sequant::ResultPtr operator()(sequant::Variable const& var) const {
    using result_t = sequant::ResultScalar<NumericT>;

    auto make_var = []() {
      return sequant::eval_result<result_t>(
          TA::detail::MakeRandom<NumericT>::generate_value());
    };

    return label_to_er_.try_emplace(std::wstring{var.label()}, make_var())
        .first->second;
  }

  sequant::ResultPtr operator()(sequant::Power const& pw) const {
    using result_t = sequant::ResultScalar<NumericT>;

    // evaluate base
    NumericT base_val;
    if (pw.base()->template is<sequant::Constant>()) {
      base_val = pw.base()
                     ->template as<sequant::Constant>()
                     .template value<NumericT>();
    } else {
      SEQUANT_ASSERT(pw.base()->template is<sequant::Variable>());
      auto base_result = (*this)(pw.base()->template as<sequant::Variable>());
      base_val = base_result->template get<NumericT>();
    }

    auto exp_val = static_cast<double>(pw.exponent());
    return sequant::eval_result<result_t>(
        static_cast<NumericT>(std::pow(base_val, exp_val)));
  }

  sequant::ResultPtr operator()(
      sequant::meta::can_evaluate auto const& node) const {
    using namespace sequant;
    if (node->is_tensor())
      return (*this)(node->as_tensor(), node->canon_indices(),
                     std::string(node->annot()));

    if (node->is_variable()) return (*this)(node->as_variable());

    if (node->is_power()) return (*this)(node->as_power());

    SEQUANT_ASSERT(node->is_constant());

    using result_t = ResultScalar<NumericT>;

    auto d = (node->as_constant()).template value<NumericT>();
    return eval_result<result_t>(d);
  }

  /// A tensor leaf by itself: laid out per its own canonical node.
  sequant::ResultPtr operator()(sequant::Tensor const& tnsr) const {
    sequant::EvalExprTA const canon{tnsr};
    return (*this)(canon.as_tensor(), canon.canon_indices(), canon.annot());
  }

  /// Slot layout of a tensor: outer modes = the proto-free bra/ket/aux
  /// indices in slot order followed by the composites' proto atoms (each
  /// once, by label, in slot order); inner modes = the composites in slot
  /// order. The one physical layout every occurrence of a tensor is
  /// generated in; a node's canon_indices() is a permutation of it.
  struct SlotLayout {
    sequant::container::svector<sequant::Index> outer, inner;
  };
  static SlotLayout slot_layout(sequant::Tensor const& tnsr) {
    using namespace sequant;
    SlotLayout out;
    auto has = [](auto const& cont, Index const& ix) {
      return ranges::any_of(cont, [&ix](Index const& x) {
        return x.full_label() == ix.full_label();
      });
    };
    for (Index const& ix : tnsr.const_braket_indices())
      if (ix.has_proto_indices()) {
        if (!has(out.inner, ix)) out.inner.push_back(ix);
      } else if (!has(out.outer, ix)) {
        out.outer.push_back(ix);
      }
    for (Index const& ix : tnsr.const_braket_indices())
      for (Index const& pr : ix.proto_indices())
        if (!has(out.outer, pr)) out.outer.push_back(pr);
    for (Index const& ix : tnsr.aux()) {
      SEQUANT_ASSERT(!ix.has_proto_indices() &&
                     "Aux indices with proto indices not supported");
      if (!has(out.outer, ix)) out.outer.push_back(ix);
    }
    return out;
  }
  /// TA annotation naming modes by slot id: "s0,s1,...;c0,..."
  static std::string slot_annot(SlotLayout const& sl) {
    std::string a;
    for (std::size_t k = 0; k < sl.outer.size(); ++k)
      a += (k ? "," : "") + std::string("s") + std::to_string(k);
    for (std::size_t k = 0; k < sl.inner.size(); ++k)
      a += (k ? "," : ";") + std::string("c") + std::to_string(k);
    return a;
  }
  /// The requested layout (\p canon, the node's canon_indices()) in slot ids.
  static std::string canon_slot_annot(
      SlotLayout const& sl, sequant::EvalExpr::index_vector const& canon) {
    using namespace sequant;
    auto pos = [](auto const& cont,
                  Index const& ix) -> std::optional<std::size_t> {
      for (std::size_t k = 0; k < cont.size(); ++k)
        if (cont[k].full_label() == ix.full_label()) return k;
      return std::nullopt;
    };
    std::string outer, inner;
    for (Index const& ix : canon) {
      if (ix.has_proto_indices()) {
        auto k = pos(sl.inner, ix);
        SEQUANT_ASSERT(k && "canonical composite not among the tensor's slots");
        inner +=
            (inner.empty() ? ";" : ",") + std::string("c") + std::to_string(*k);
      } else {
        auto k = pos(sl.outer, ix);
        SEQUANT_ASSERT(k && "canonical index not among the tensor's slots");
        outer +=
            (outer.empty() ? "" : ",") + std::string("s") + std::to_string(*k);
      }
    }
    return outer + inner;
  }

  /// The leaf array for \p tnsr laid out per \p canon (the requesting node's
  /// canon_indices()). One random array per tensor (keyed by tensor_to_key)
  /// generated in slot layout; a request is served as a permutation of it to
  /// the requested layout (cached per layout), so every occurrence sees the
  /// same tensor.
  sequant::ResultPtr operator()(sequant::Tensor const& tnsr,
                                sequant::EvalExpr::index_vector const& canon,
                                std::string const& /*annot*/) const {
    using namespace ranges::views;
    using namespace sequant;

    std::wstring const label = tensor_to_key(tnsr);
    SlotLayout const sl = slot_layout(tnsr);
    std::string const requested = canon_slot_annot(sl, canon);
    std::string const stored = slot_annot(sl);
    auto const serve = [&](ResultPtr const& base) -> ResultPtr {
      if (std::getenv("SEQUANT_UT_YIELD_DIAG"))
        std::cerr << "[yield] key=" << sequant::toUtf8(label)
                  << " stored=" << stored << " requested=" << requested
                  << " base_layout=" << base->layout_desc() << std::endl;
      if (requested == stored) return base;
      std::wstring const key =
          label + L"|" + std::wstring(requested.begin(), requested.end());
      if (auto&& f2 = label_to_er_.find(key); f2 != label_to_er_.end())
        return f2->second;
      auto permuted = base->permute(std::array<std::any, 2>{stored, requested});
      label_to_er_.emplace(key, permuted);
      return permuted;
    };
    if (auto&& found = label_to_er_.find(label); found != label_to_er_.end())
      return serve(found->second);

    // generate in slot layout
    struct {
      container::svector<Index> outer, inner;
    } nested{sl.outer, sl.inner};

    ResultPtr result{nullptr};
    auto ibr = get_default_context().index_basis_registry();

    auto make_extents = [this, &ibr](auto&& ixs) -> container::svector<size_t> {
      return ixs | transform([this, &ibr](auto const& ix) -> size_t {
               SEQUANT_ASSERT(ix.space() == ibr->retrieve_space(L"i") ||
                              ix.space() == ibr->retrieve_space(L"a") ||
                              ix.space() == ibr->retrieve_space(L"x") ||
                              ix.space() == ibr->retrieve_space(L"m"));
               return ix.space() == ibr->retrieve_space(L"i")   ? nocc_
                      : ix.space() == ibr->retrieve_space(L"a") ? nvirt_
                      : ix.space() == ibr->retrieve_space(L"x")
                          ? naux_
                          : nvirt_;  // m (flat
                                     // PAO-like
                                     // leg)
             }) |
             ranges::to<container::svector<size_t>>;
    };

    auto const outer_extent = make_extents(nested.outer);

    auto const outer_tr = [&outer_extent, &nested, this]() {
      container::vector<TA::TiledRange1> tr1s;
      tr1s.reserve(outer_extent.size());
      std::size_t k = 0;
      for (auto const& ix : nested.outer) {
        tr1s.emplace_back(
            tr1_for(std::wstring(ix.space().base_key()), outer_extent[k]));
        ++k;
      }
      return TA::TiledRange(tr1s.begin(), tr1s.end());
    }();

    auto const outer_r = TA::Range(outer_extent);

    if (nested.inner.empty()) {
      // regular tensor
      using ArrayT = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
      ArrayT array{world, outer_tr};
      // SEQUANT_UT_YIELD_SPARSITY=<p>: under a sparse policy, leave a
      // deterministic fraction p of the tiles ZERO (real block sparsity), so
      // sparse-shape estimation and truncation are exercised.
      static double const sparsity = [] {
        char const* e = std::getenv("SEQUANT_UT_YIELD_SPARSITY");
        return e ? std::atof(e) : 0.0;
      }();
      std::size_t ord = 0;
      for (auto it = array.begin(); it != array.end(); ++it, ++ord) {
        if (!array.is_local(it.index())) continue;
        bool const zero =
            sparsity > 0.0 &&
            (static_cast<double>((ord * 2654435761u) % 1000u) / 1000.0 <
             sparsity);
        if (zero)
          *it = world.taskq.add(
              [](TA::Range const& r) {
                return TA::Tensor<NumericT>(r, NumericT{0});
              },
              it.make_range());
        else
          *it = world.taskq.add(random_tensor<NumericT>, it.make_range());
      }
      if constexpr (!std::is_same_v<TAPolicyT, TA::DensePolicy>)
        if (sparsity > 0.0) array.truncate();
      result = eval_result<ResultTensorTA<ArrayT>>(array);
    } else {
      // tensor of tensor
      using ArrayT = TA::DistArray<tot_tile_t<NumericT>, TAPolicyT>;

      auto const inner_extent = make_extents(nested.inner);
      auto const inner_r = TA::Range(inner_extent);

      auto make_tile = [inner_r](TA::Range const& orng) {
        return random_tensor_of_tensor<NumericT>(orng, inner_r);
      };

      ArrayT array{world, outer_tr};

      for (auto it = array.begin(); it != array.end(); ++it)
        if (array.is_local(it.index()))
          *it = world.taskq.add(make_tile, it.make_range());

      result = eval_result<ResultTensorOfTensorTA<ArrayT>>(array);
    }

    auto success = label_to_er_.emplace(label, result);
    SEQUANT_ASSERT(success.second && "couldn't store ResultPtr!");
    SEQUANT_ASSERT(success.first->second);
    return serve(success.first->second);
  }

  ///
  /// \param label eg.
  ///  - 't_vvoo', 'f_ov' for generic tensor key strings
  ///  - 't{a1,a2;i1,i2}', 'f{i1;a1}' supported by sequant::parse_expr
  ///
  /// \return ResultPtr
  ///
  /// \note The ResultPtr should already exist in the cache otherwise throws.
  ///       This overload is only intended to access already existing ERPtrs
  ///       from the cache. To create a new cache entry use the
  ///       operator()(Tensor const&) overload.
  ///
  sequant::ResultPtr operator()(std::wstring_view label) const {
    auto&& found = label_to_er_.find(label.data());
    if (found == label_to_er_.end())
      found = label_to_er_.find(tensor_to_key(label));
    if (found == label_to_er_.end())
      throw sequant::Exception{"attempted access of non-existent ResultPtr!"};
    return found->second;
  }
};

enum struct ErrorTol : int { Loose = 1000, Normal = 100, Tight = 2 };

using enum ErrorTol;

template <ErrorTol Tol, std::floating_point T>
constexpr bool approx_equal(T val1, T val2) {
  constexpr auto margin =
      static_cast<int>(Tol) * std::numeric_limits<T>::epsilon();
  return (val1 - val2) == Catch::Approx(0.).margin(margin);
}

template <ErrorTol Tol = Normal, TA::tnsr_expr ArrExpr>
bool equal_tarrays(ArrExpr const& arr1, ArrExpr const& arr2) {
  typename ArrExpr::array_type diff;
  diff(arr1.annotation()) = arr1 - arr2;
  return approx_equal<Tol>(TA::norm2(diff), 0.);
}

template <ErrorTol Tol = Normal, TA::array Array>
bool equal_tarrays(Array arr1, Array arr2, std::string const& annot1,
                   std::string const& annot2) {
  return equal_tarrays<Tol>(arr1(annot1), arr2(annot2));
}

template <ErrorTol Tol = Normal, TA::array Array>
bool equal_tarrays(Array const& arr1,  //
                   Array const& arr2,  //
                   std::string const& annot) {
  return equal_tarrays<Tol>(arr1, arr2, annot, annot);
}

template <ErrorTol Tol = Normal, TA::detail::array_tos Array>
bool equal_tarrays(Array const& arr1, Array const& arr2) {
  return equal_tarrays<Tol>(arr1, arr2,
                            TA::detail::dummy_annotation(rank(arr1)),
                            TA::detail::dummy_annotation(rank(arr2)));
}

}  // namespace sequant::tests::eval_ta

#endif  // SEQUANT_TESTS_UNIT_TEST_EVAL_TA_FIXTURES_HPP
