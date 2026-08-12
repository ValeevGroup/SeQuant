#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/backends/tiledarray/eval_context.hpp>
#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/schedule_dump.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <tiledarray.h>
#include <boost/regex.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/intersperse.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <string>
#include <vector>

// Force compile-instantiation of the complex tensor-of-tensors Result so its
// adjoint() override (`result(annot) = arr(annot).conj()`, relying on TA's
// recursive conj for nested tiles) is type-checked. No TA eval test constructs
// a complex ToT adjoint, and Result::adjoint() is private (reachable only
// through the EvalOp::Adjoint IR node); the ta_tot_conj_complex test below
// runtime-checks the underlying TA conj while this instantiation compile-checks
// the override.
template class sequant::ResultTensorOfTensorTA<
    TA::DistArray<TA::Tensor<TA::Tensor<std::complex<double>>>>>;

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

namespace {

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

auto to_ta_node(sequant::FullBinaryNode<sequant::EvalExpr> node) {
  using namespace sequant;
  return transform_node(node, [](auto&& val) {
    if (val.is_tensor()) {
      return EvalExprTA(*val.op_type(), val.result_type(), val.expr(),
                        NestedTensorIndices(val.as_tensor()).outer_inner() |
                            ranges::to<EvalExpr::index_vector>(),
                        val.canon_phase(), val.hash_value(),
                        val.copy_connectivity_graph());
    } else
      return EvalExprTA(val);
  });
}

auto eval_node(sequant::ExprPtr const& expr) {
  // sequant::binarize(ExprPtr) is deprecated for caller-visible head
  // construction; this helper exists for legacy test sections that don't
  // depend on the head's bra/ket layout.
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return to_ta_node(sequant::binarize(expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

auto eval_node(sequant::ResultExpr const& res) {
  return to_ta_node(sequant::binarize(res));
}

auto tensor_to_key(sequant::Tensor const& tnsr) {
  static auto const idx_rgx = boost::wregex{L"([iax])([↑↓])?(_?\\d+)"};
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i"   ? L"o"
            : mo[1].str() == L"a" ? L"v"
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

auto tensor_to_key(std::wstring_view spec) {
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

// note: all the inner tensors (elements of the outer tensor)
//       have the same @c inner_rng
template <typename NumericT>
auto random_tensor_of_tensor(TA::Range const& outer_rng,
                             TA::Range const& inner_rng) {
  TA::Tensor<TA::Tensor<NumericT>> result{outer_rng};

  std::generate(result.begin(), result.end(),
                [inner_rng]() { return random_tensor<NumericT>(inner_rng); });

  return result;
}

template <typename NumericT = double, typename TAPolicyT = TA::DensePolicy>
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

 public:
  /// Produce arrays whose modes are tiled in blocks of at most \p n.
  /// \p n must be positive (0 would make the tiling loop in make_tr1 spin).
  void set_max_tile(size_t n) {
    REQUIRE(n > 0);
    max_tile_ = n;
  }

  using array_type = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
  using array_tot_type =
      TA::DistArray<TA::Tensor<TA::Tensor<NumericT>>, TAPolicyT>;
  using numeric_type = NumericT;

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt}, naux_{nvirt * 2} {}

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt, size_t naux)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt}, naux_{naux} {}

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
    if (node->is_tensor()) return (*this)(node->as_tensor());

    if (node->is_variable()) return (*this)(node->as_variable());

    if (node->is_power()) return (*this)(node->as_power());

    SEQUANT_ASSERT(node->is_constant());

    using result_t = ResultScalar<NumericT>;

    auto d = (node->as_constant()).template value<NumericT>();
    return eval_result<result_t>(d);
  }

  sequant::ResultPtr operator()(sequant::Tensor const& tnsr) const {
    using namespace ranges::views;
    using namespace sequant;

    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_er_.find(label); found != label_to_er_.end()) {
      //      std::cout << "label = [" << sequant::to_string(label)
      //                << "] FOUND in cache. Returning.." << std::endl;
      return found->second;
    }

    ResultPtr result{nullptr};
    auto isr = get_default_context().index_space_registry();

    auto make_extents = [this, &isr](auto&& ixs) -> container::svector<size_t> {
      return ixs | transform([this, &isr](auto const& ix) -> size_t {
               SEQUANT_ASSERT(ix.space() == isr->retrieve(L"i") ||
                              ix.space() == isr->retrieve(L"a") ||
                              ix.space() == isr->retrieve(L"x") ||
                              ix.space() == isr->retrieve(L"m"));
               return ix.space() == isr->retrieve(L"i")   ? nocc_
                      : ix.space() == isr->retrieve(L"a") ? nvirt_
                      : ix.space() == isr->retrieve(L"x") ? naux_
                                                          : nvirt_;  // m (flat
                                                                     // PAO-like
                                                                     // leg)
             }) |
             ranges::to<container::svector<size_t>>;
    };

    NestedTensorIndices nested{tnsr};

    auto const outer_extent = make_extents(nested.outer);

    auto const outer_tr = [&outer_extent, this]() {
      auto make_tr1 = [this](size_t e) {
        container::svector<size_t> b;
        for (size_t x = 0; x < e; x += max_tile_) b.push_back(x);
        b.push_back(e);
        return TA::TiledRange1(b.begin(), b.end());
      };
      container::vector<TA::TiledRange1> tr1s;
      tr1s.reserve(outer_extent.size());
      for (auto e : outer_extent) tr1s.emplace_back(make_tr1(e));
      return TA::TiledRange(tr1s.begin(), tr1s.end());
    }();

    auto const outer_r = TA::Range(outer_extent);

    if (nested.inner.empty()) {
      // regular tensor
      using ArrayT = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
      ArrayT array{world, outer_tr};
      for (auto it = array.begin(); it != array.end(); ++it)
        if (array.is_local(it.index()))
          *it = world.taskq.add(random_tensor<NumericT>, it.make_range());
      result = eval_result<ResultTensorTA<ArrayT>>(array);
    } else {
      // tensor of tensor
      using ArrayT = TA::DistArray<TA::Tensor<TA::Tensor<NumericT>>, TAPolicyT>;

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
    //    std::cout << "label = [" << sequant::to_string(label)
    //              << "] NotFound in cache. Creating.." << std::endl;
    SEQUANT_ASSERT(success.first->second);
    return success.first->second;
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
      throw std::runtime_error{"attempted access of non-existent ResultPtr!"};
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

}  // namespace

// Regression tests for tensor-of-tensor (ToT) ops on arrays whose inner tiles
// are all empty (e.g. a fully screened CSV residual). tot_inner_rank() reads
// the rank from a populated inner tile, so it returns 0 here; downstream ToT
// ops must not turn that 0 into a degenerate (non-ToT) annotation. Reproduces
// the mpqc4 PAO-CSV abort fixed in add_inplace and column_symmetrize_ta.
TEST_CASE("tot_all_empty_inner", "[eval][tot]") {
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultTensorOfTensorTA;
  using sequant::detail::tot_inner_rank;
  using ToTArray = TA::DistArray<TA::Tensor<TA::Tensor<double>>>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;

  auto& world = TA::get_default_world();

  // Build a ToT array; `empty_inner` leaves every inner tensor empty.
  auto build = [&world](unsigned outer_rank, bool empty_inner,
                        TA::TiledRange1 outer_tr1 =
                            TA::TiledRange1{0, 3}) -> ToTArray {
    std::vector<TA::TiledRange1> tr1s(outer_rank, outer_tr1);
    TA::TiledRange outer_tr(tr1s.begin(), tr1s.end());
    std::vector<std::size_t> inner_ext(outer_rank, 4);
    TA::Range inner_r(inner_ext);

    auto tile_fn = [inner_r, empty_inner](TA::Range const& orng) {
      TA::Tensor<TA::Tensor<double>> t{orng};
      if (!empty_inner)
        for (auto& inner : t) {
          inner = TA::Tensor<double>{inner_r};
          std::fill(inner.begin(), inner.end(), 1.0);
        }
      return t;
    };

    ToTArray arr{world, outer_tr};
    for (auto it = arr.begin(); it != arr.end(); ++it)
      if (arr.is_local(it.index()))
        *it = world.taskq.add(tile_fn, it.make_range());
    world.gop.fence();
    return arr;
  };

  SECTION("add_inplace: empty += populated adopts the populated inner rank") {
    auto empty = build(1, /*empty_inner=*/true);
    auto full = build(1, /*empty_inner=*/false);
    REQUIRE(tot_inner_rank(empty) == 0);
    REQUIRE(tot_inner_rank(full) == 1);

    auto r_empty = eval_result<ResultToT>(empty);
    auto r_full = eval_result<ResultToT>(full);
    REQUIRE_NOTHROW(r_empty->add_inplace(*r_full));
    auto const& res = r_empty->get<ToTArray>();
    REQUIRE(tot_inner_rank(res) == 1);
    // the populated operand's data (all 1.0) must actually have landed.
    for (auto it = res.begin(); it != res.end(); ++it)
      if (res.is_local(it.index()))
        for (auto const& inner : it->get()) {
          REQUIRE_FALSE(inner.empty());
          for (auto const& x : inner) REQUIRE(x == 1.0);
        }
  }

  SECTION("add_inplace: populated += empty is a no-op") {
    auto full = build(1, false);
    auto empty = build(1, true);
    auto r_full = eval_result<ResultToT>(full);
    auto r_empty = eval_result<ResultToT>(empty);
    REQUIRE_NOTHROW(r_full->add_inplace(*r_empty));
    REQUIRE(tot_inner_rank(r_full->get<ToTArray>()) == 1);
  }

  SECTION("add_inplace: empty += empty is an identity no-op") {
    auto e1 = build(1, true);
    auto e2 = build(1, true);
    auto r1 = eval_result<ResultToT>(e1);
    auto r2 = eval_result<ResultToT>(e2);
    REQUIRE_NOTHROW(r1->add_inplace(*r2));
    REQUIRE(tot_inner_rank(r1->get<ToTArray>()) == 0);
  }

  SECTION("symmetrize: all-empty ToT is the identity") {
    auto empty = build(1, true);
    REQUIRE(tot_inner_rank(empty) == 0);
    auto r_empty = eval_result<ResultToT>(empty);
    ResultPtr sym;
    REQUIRE_NOTHROW(sym = r_empty->symmetrize());
    REQUIRE(tot_inner_rank(sym->get<ToTArray>()) == 0);
  }

  SECTION("slice_array_over_mode: all-empty ToT slices to a zero ToT") {
    auto empty = build(1, true);
    REQUIRE(tot_inner_rank(empty) == 0);
    ToTArray sliced;
    REQUIRE_NOTHROW(sliced = sequant::slice_array_over_mode(empty, 0, 0, 1));
    REQUIRE(tot_inner_rank(sliced) == 0);
    REQUIRE(sliced.trange().dim(0).tile_extent() == 1);
  }

  SECTION(
      "slice_array_over_mode: rebasing matches block() on a multi-tile cut") {
    // Slice across two tiles at a non-zero base; hand-built trange must match
    // what the populated path gets from TA's block().
    TA::TiledRange1 const tr1{0, 3, 6, 9};
    auto empty = build(1, /*empty_inner=*/true, tr1);
    auto full = build(1, /*empty_inner=*/false, tr1);
    auto const s_empty = sequant::slice_array_over_mode(empty, 0, 1, 3);
    auto const s_full = sequant::slice_array_over_mode(full, 0, 1, 3);
    REQUIRE(tot_inner_rank(s_empty) == 0);
    REQUIRE(s_empty.trange() == s_full.trange());
  }
}

TEST_CASE("eval_with_tiledarray", "[eval]") {
  // Reproducer for the mpqc4 cck real-field NaN regression. The eval-graph
  // head's bra/ket split is purely positional: each external index ends up in
  // whichever slot (bra or ket) it occupied in its source tensor, summed
  // across the term. Two orientations of a tensor that are equivalent under
  // bra<->ket-swap (Symm braket_symmetry) produce different head bra_rank,
  // which breaks downstream code (e.g. mpqc's jacobi_update) that assumes a
  // conventional 2:2 (vir,vir;occ,occ) layout for a CCSD T2 residual head.
  // Bug is independent of scalar Field — fires under both Conjugate and Symm
  // whenever the canonical orientation puts an external on the "wrong" side.
  SECTION("eval-graph head bra/ket split is positional, not external-aware") {
    using sequant::deserialize;
    using sequant::EvalExprTA;
    using sequant::ExprPtr;

    // The test expressions below put an internal index (e.g. a_3) in the bra
    // of two factors (g.bra and t.bra), which violates the default-context
    // strict bra↔ket-symmetry policy (each internal must appear at most once
    // per side under Conjugate). Disable the policy for this scope; we are
    // probing the eval-graph head's bra/ket layout, not the canonicalizer's
    // covariance assumptions.
    auto ctx_resetter = sequant::set_scoped_default_context(
        sequant::Context{sequant::get_default_context()}.set(
            sequant::AssertStrictBraKetSymmetry::No));

    auto report = [](sequant::ExprPtr const& e, std::string const& label) {
      auto node = eval_node(e);
      auto const& head = node->as_tensor();
      std::wstring head_str = sequant::to_latex(head);
      std::string head_str8{head_str.begin(), head_str.end()};
      INFO(label + " head: " + head_str8 +
           "  bra_rank=" + std::to_string(head.bra_rank()) +
           "  ket_rank=" + std::to_string(head.ket_rank()));
      return std::make_pair(head.bra_rank(), head.ket_rank());
    };

    // Representative SF R2 residual term (h2o-cck-2-631g-pvdz) with externals
    // {a_1, a_2, i_1, i_2}. For a CCSD T2 the head ought to be I{a_1,a_2; i_1,
    // i_2} (vir,vir;occ,occ) so downstream Jacobi-style updates index orbital
    // energies correctly.

    // (A) g written with the occ external i_2 in its bra slot: regardless of
    // braket_symmetry, the head ends up I{i_2,a_1,a_2; i_1} — head bra/ket is
    // assigned by external SLOT, not by space.
    auto expr_bra_external_symm = deserialize(
        L"2 g{i_2,a_3;i_3,i_4}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_bra_external_symm);
    auto [br_symm, kr_symm] =
        report(expr_bra_external_symm, "g{i_2,a_3;...}:N-S-S");

    auto expr_bra_external_conj = deserialize(
        L"2 g{i_2,a_3;i_3,i_4}:N-C-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_bra_external_conj);
    auto [br_conj, kr_conj] =
        report(expr_bra_external_conj, "g{i_2,a_3;...}:N-C-S");

    // (B) Same expression with g's bra/ket pre-swapped (mathematically
    // equivalent under Symm braket_symmetry). Braket-symmetric tensors now
    // canonicalize by bra/ket *color* (index spaces), so g canonicalizes to the
    // same orientation as in (A) regardless of how it is written: the head
    // comes out the same 3:1 as (A). Pre-swapping no longer steers the head
    // layout — use the ResultExpr API (C) to pin it.
    auto expr_ket_external_swap = deserialize(
        L"2 g{i_3,i_4;i_2,a_3}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_ket_external_swap);
    auto [br_swap, kr_swap] =
        report(expr_ket_external_swap, "g{i_3,i_4;i_2,a_3}:N-S-S (pre-swap)");

    // (C) Caller-supplied head layout via ResultExpr: the right shape of the
    // public API for this. The caller writes the LHS with the bra/ket layout
    // it wants the head to have; binarize(ResultExpr) at eval_expr.hpp:435
    // overwrites the eval-tree root's tensor with res.result_as_tensor(),
    // making the IR's positional choice irrelevant.
    auto res_explicit_layout = sequant::deserialize<sequant::ResultExpr>(
        L"R2{a_1,a_2;i_1,i_2}:N-N-S = "
        L"2 g{i_2,a_3;i_3,i_4}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    auto node_explicit = eval_node(res_explicit_layout);
    auto const& head_explicit = node_explicit->as_tensor();
    std::wstring head_explicit_str = sequant::to_latex(head_explicit);
    std::string head_explicit_str8{head_explicit_str.begin(),
                                   head_explicit_str.end()};
    INFO("ResultExpr R2{a_1,a_2;i_1,i_2} head: " + head_explicit_str8 +
         "  bra_rank=" + std::to_string(head_explicit.bra_rank()) +
         "  ket_rank=" + std::to_string(head_explicit.ket_rank()));

    // Document the current behavior:
    // (A) same 3:1 split regardless of Symm vs Conjugate — head bra/ket is
    // assigned positionally (by external slot), independent of scalar Field.
    CHECK(br_symm == 3);
    CHECK(kr_symm == 1);
    CHECK(br_conj == 3);
    CHECK(kr_conj == 1);
    // (B) the pre-swapped form canonicalizes (by bra/ket color) to the same
    // orientation as (A), so the head is the same 3:1 — pre-swapping is no
    // longer a lever on head layout.
    CHECK(br_swap == 3);
    CHECK(kr_swap == 1);
    // (C) ResultExpr API gives the caller exact control; head matches the LHS
    // verbatim no matter how the RHS factors are oriented internally.
    CHECK(head_explicit.bra_rank() == 2);
    CHECK(head_explicit.ket_rank() == 2);
    CHECK(head_explicit.bra().at(0).label() == L"a_1");
    CHECK(head_explicit.bra().at(1).label() == L"a_2");
    CHECK(head_explicit.ket().at(0).label() == L"i_1");
    CHECK(head_explicit.ket().at(1).label() == L"i_2");
  }

  SECTION("real") {
    using ranges::views::transform;
    using sequant::EvalExprTA;
    using sequant::evaluate;
    using sequant::evaluate_antisymm;
    using sequant::evaluate_symm;

    using TA::TArrayD;

    auto parse_antisymm = [](auto const& xpr) {
      return sequant::deserialize<sequant::ExprPtr>(
          xpr, {.def_perm_symm = sequant::Symmetry::Antisymm});
    };

    auto& world = TA::get_default_world();

    const size_t nocc = 2, nvirt = 20;
    auto yield_ =
        rand_tensor_yield<double, TA::DensePolicy>{world, nocc, nvirt};
    auto yield = [&yield_](std::wstring_view lbl) -> TA::TArrayD const& {
      return yield_(lbl)->get<TA::TArrayD>();
    };

    auto yield_d = [&yield_](std::wstring_view lbl) ->
        typename TA::TArrayD::numeric_type {
          return yield_(lbl)->get<typename TA::TArrayD::numeric_type>();
        };

    auto eval = [&yield_](sequant::ExprPtr const& expr,
                          std::string const& target_labels) {
      return evaluate(eval_node(expr), target_labels, yield_)
          ->get<TA::TArrayD>();
    };

    auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                               std::string const& target_labels) {
      return evaluate_symm(eval_node(expr), target_labels, yield_)
          ->get<TA::TArrayD>();
    };

    auto eval_antisymm = [&yield_](sequant::ExprPtr const& expr,
                                   std::string const& target_labels) {
      return evaluate_antisymm(eval_node(expr), target_labels, yield_)
          ->get<TA::TArrayD>();
    };

    auto eval_biorthogonal_nns_project = [&yield_](
                                             sequant::ExprPtr const& expr,
                                             std::string const& target_labels) {
      auto result = evaluate(eval_node(expr), target_labels, yield_);
      return sequant::mbpt::biorthogonal_nns_project(
          result->get<TA::TArrayD>(), eval_node(expr)->as_tensor().bra_rank());
    };

    SECTION("summation") {
      auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayD{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum1_eval, sum1_man));

      auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");
      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("product") {
      auto expr1 =
          parse_antisymm(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayD prod1_man{};
      prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(equal_tarrays(prod1_eval, prod1_man));

      auto expr2 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayD{};
      prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      REQUIRE(equal_tarrays(prod2_eval, prod2_man));

      auto expr3 = sequant::deserialize<sequant::ExprPtr>(
          L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayD{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");
      REQUIRE(equal_tarrays(prod3_eval, prod3_man));

      auto expr4 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayD{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(equal_tarrays(prod4_eval, prod4_man));
    }

    SECTION("sum and product") {
      auto expr1 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
      auto eval1 = eval(expr1, "a_1,a_2,i_1,i_2");

      auto man1 = TArrayD{};
      man1("a1,a2,i1,i2") = -1.0 / 4 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4") +
                            1.0 / 16 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a1,a2;i3,i4}")("a1,a2,i3,i4") *
                                yield(L"t{a3,a4;i1,i2}")("a3,a4,i1,i2");
      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayD{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("variable at leaves") {
      auto expr2 =
          parse_antisymm(L"(α * 2 * t_{a1}^{i1} * β) + (3/2 * f_{i1}^{a1})");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          yield_d(L"α") * 2 * yield(L"t{a1;i1}")("a1,i1") * yield_d(L"β") +
          1.5 * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("power at leaves") {
      using sequant::Constant;
      using sequant::ex;
      using sequant::Power;
      using sequant::rational;
      using sequant::Variable;

      {
        auto pw = ex<Power>(ex<Constant>(2), rational{1, 2});
        auto t = parse_antisymm(L"t_{a1}^{i1}");
        auto expr = pw * t;

        auto eval1 = eval(expr, "i_1,a_1");

        auto man1 = TArrayD{};
        man1("i1,a1") = std::sqrt(2.0) * yield(L"t{a1;i1}")("a1,i1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Variable>(L"α"), rational{2, 1});
        auto f = parse_antisymm(L"f_{i1}^{a1}");
        auto expr = pw * f;
        auto eval1 = eval(expr, "i_1,a_1");

        auto alpha_val = yield_d(L"α");
        auto man1 = TArrayD{};
        man1("i1,a1") = (alpha_val * alpha_val) * yield(L"f{i1;a1}")("i1,a1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Constant>(2), rational{1, 2});
        auto t = parse_antisymm(L"t_{a1}^{i1}");
        auto alpha = ex<Variable>(L"α");
        auto f = parse_antisymm(L"f_{i1}^{a1}");

        auto expr = pw * t + alpha * f;

        auto eval1 = eval(expr, "i_1,a_1");

        auto man1 = TArrayD{};
        man1("i1,a1") = std::sqrt(2.0) * yield(L"t{a1;i1}")("a1,i1") +
                        yield_d(L"α") * yield(L"f{i1;a1}")("i1,a1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Constant>(rational{1, 3}), rational{2, 1});
        auto g = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
        auto expr = pw * g;

        auto eval1 = eval(expr, "i_1,i_2,a_1,a_2");

        auto man1 = TArrayD{};
        man1("i1,i2,a1,a2") =
            std::pow(1.0 / 3.0, 2.0) * yield(L"g{i1,i2;a1,a2}")("i1,i2,a1,a2");
        REQUIRE(equal_tarrays(eval1, man1));
      }
    }

    SECTION("Antisymmetrization") {
      auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = 0.25 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // odd-ranked tensor
      auto expr2 = parse_antisymm(L"g_{i1, i2, i3}^{a1, a2}");
      auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2}");

      auto man2 = TArrayD{};
      man2("0,1,2,3,4") =
          arr2("0,1,2,3,4") - arr2("1,0,2,3,4") + arr2("1,2,0,3,4") -
          arr2("2,1,0,3,4") + arr2("2,0,1,3,4") - arr2("0,2,1,3,4") -
          arr2("0,1,2,4,3") + arr2("1,0,2,4,3") - arr2("1,2,0,4,3") +
          arr2("2,1,0,4,3") - arr2("2,0,1,4,3") + arr2("0,2,1,4,3");

      REQUIRE(equal_tarrays(eval2, man2));

      auto expr3 = parse_antisymm(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayD{};
      man3("0,1") = arr3("0,1") - arr3("1,0");
      man3("0,1") = 0.5 * man3("0,1");

      REQUIRE(equal_tarrays(eval3, man3));
    }

    SECTION("Symmetrization") {
      auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = parse_antisymm(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayD man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");
      man2("0,1,2,3,4,5") = (1.0 / 6.0) * man2("0,1,2,3,4,5");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Biorthogonal Cleanup") {
      // low-rank residuals: skip nns
      auto expr1 = parse_antisymm(L"R_{a1, a2}^{i1, i2}");
      auto eval1 = eval_biorthogonal_nns_project(expr1, "a_1,a_2,i_1,i_2");
      auto const& arr1 = yield(L"R{a1,a2;i1,i2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // for rank 3 residual, nns applies:
      // result = NNS_P * sum_of_ket_permutations
      auto expr2 = parse_antisymm(L"R_{a1, a2, a3}^{i1, i2, i3}");
      auto eval2 =
          eval_biorthogonal_nns_project(expr2, "a_1,a_2,a_3,i_1,i_2,i_3");
      auto const& arr2 = yield(L"R{a1,a2,a3;i1,i2,i3}");

      auto man2 = TArrayD{};
      man2("0,1,2,3,4,5") =
          arr2("0,1,2,3,4,5") -
          (1.0 / 5.0) *
              (arr2("0,1,2,3,5,4") + arr2("0,1,2,4,3,5") + arr2("0,1,2,4,5,3") +
               arr2("0,1,2,5,3,4") + arr2("0,1,2,5,4,3"));

      REQUIRE(equal_tarrays(eval2, man2));

      // for rank 4 residual, nns applies:
      // result = NNS_P * sum_of_ket_permutations
      auto expr3 = parse_antisymm(L"R_{a1, a2, a3, a4}^{i1, i2, i3, i4}");
      auto eval3 = eval_biorthogonal_nns_project(
          expr3, "a_1,a_2,a_3,a_4,i_1,i_2,i_3,i_4");
      auto const& arr3 = yield(L"R{a1,a2,a3,a4;i1,i2,i3,i4}");

      auto man3 = TArrayD{};
      man3("0,1,2,3,4,5,6,7") = 1.0 * arr3("0,1,2,3,4,5,6,7") +
                                -4.0 / 14.0 * arr3("0,1,2,3,4,5,7,6") +
                                -4.0 / 14.0 * arr3("0,1,2,3,4,6,5,7") +
                                -1.0 / 14.0 * arr3("0,1,2,3,4,6,7,5") +
                                -1.0 / 14.0 * arr3("0,1,2,3,4,7,5,6") +
                                -4.0 / 14.0 * arr3("0,1,2,3,4,7,6,5") +
                                -4.0 / 14.0 * arr3("0,1,2,3,5,4,6,7") +
                                2.0 / 14.0 * arr3("0,1,2,3,5,4,7,6") +
                                -1.0 / 14.0 * arr3("0,1,2,3,5,6,4,7") +
                                2.0 / 14.0 * arr3("0,1,2,3,5,6,7,4") +
                                2.0 / 14.0 * arr3("0,1,2,3,5,7,4,6") +
                                -1.0 / 14.0 * arr3("0,1,2,3,5,7,6,4") +
                                -1.0 / 14.0 * arr3("0,1,2,3,6,4,5,7") +
                                2.0 / 14.0 * arr3("0,1,2,3,6,4,7,5") +
                                -4.0 / 14.0 * arr3("0,1,2,3,6,5,4,7") +
                                -1.0 / 14.0 * arr3("0,1,2,3,6,5,7,4") +
                                2.0 / 14.0 * arr3("0,1,2,3,6,7,4,5") +
                                2.0 / 14.0 * arr3("0,1,2,3,6,7,5,4") +
                                2.0 / 14.0 * arr3("0,1,2,3,7,4,5,6") +
                                -1.0 / 14.0 * arr3("0,1,2,3,7,4,6,5") +
                                -1.0 / 14.0 * arr3("0,1,2,3,7,5,4,6") +
                                -4.0 / 14.0 * arr3("0,1,2,3,7,5,6,4") +
                                2.0 / 14.0 * arr3("0,1,2,3,7,6,4,5") +
                                2.0 / 14.0 * arr3("0,1,2,3,7,6,5,4");

      REQUIRE(equal_tarrays<Loose>(eval3, man3));
    }

    SECTION("Others") {
      using namespace std::string_literals;
      auto expr1 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

      auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                       ->get<TA::TArrayD>();

      auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
        return eval_node(x);
      }) | ranges::to_vector;

      auto eval2 =
          evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TA::TArrayD>();

      REQUIRE(equal_tarrays(eval1, eval2));
    }

    SECTION("non-covariant indices") {
      using sequant::deserialize;
      using sequant::EvalExprTA;
      using sequant::evaluate;

      using TA::TArrayD;
      auto& world = TA::get_default_world();
      const size_t nocc = 2, nvirt = 4, naux = 12;

      auto yield_ =
          rand_tensor_yield<double, TA::DensePolicy>{world, nocc, nvirt, naux};
      auto yield = [&yield_](std::wstring_view lbl) -> TA::TArrayD const& {
        return yield_(lbl)->get<TA::TArrayD>();
      };

      auto eval = [&yield_](sequant::ExprPtr const& expr,
                            std::string const& target_labels) {
        return evaluate(eval_node(expr), target_labels, yield_)
            ->get<TA::TArrayD>();
      };

      auto expr1 = deserialize(
          L"((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2})(X{a3;;x2} X{;a4;x2})");
      auto eval1 = eval(expr1, "a_1,a_2,a_3,a_4");
      auto man1 = [&]() {
        auto X1 = yield(L"X{a1;;x1}");
        REQUIRE(X1.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X1.trange().elements_range().extent(1) == naux);
        auto X2 = yield(L"X{;a2;x1}");
        REQUIRE(X2.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X2.trange().elements_range().extent(1) == naux);
        auto X3 = yield(L"X{a3;;x2}");
        REQUIRE(X3.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X3.trange().elements_range().extent(1) == naux);
        auto X4 = yield(L"X{;a4;x2}");
        REQUIRE(X4.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X4.trange().elements_range().extent(1) == naux);
        auto Y = yield(L"Y{;;x1,x2}");
        REQUIRE(Y.trange().elements_range().extent(0) == naux);
        REQUIRE(Y.trange().elements_range().extent(1) == naux);
        auto X12 = TA::einsum("ax,bx->abx", X1, X2);
        REQUIRE(X12.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X12.trange().elements_range().extent(1) == nvirt);
        REQUIRE(X12.trange().elements_range().extent(2) == naux);
        auto X12Y = TA::einsum("abx,xy->aby", X12, Y);
        auto X34 = TA::einsum("cy,dy->cdy", X3, X4);
        return TA::einsum("aby,cdy->abcd", X12Y, X34);
      }();
      REQUIRE(equal_tarrays(eval1, man1, "a1,a2,a3,a4"));

      // cluster-specific RDM: γ{a2;a1;i1,i2} = t{i1,i2;a1,a3} T2{a2,a3;i1,i2}
      // i1,i2 form standard bra-ket pairs across factors but are external
      // (in aux of result); binarize(ResultExpr) must keep them uncontracted
      {
        auto res = deserialize<sequant::ResultExpr>(
            L"GAM{a2;a1;i1,i2} = t{i1,i2;a1,a3} T2{a2,a3;i1,i2}");
        auto node = eval_node(res);
        auto eval_rdm = evaluate(node, std::string("a_2,a_1,i_1,i_2"), yield_)
                            ->get<TA::TArrayD>();
        auto man_rdm = [&]() {
          auto t = yield(L"t{i1,i2;a1,a3}");
          auto T2 = yield(L"T2{a2,a3;i1,i2}");
          return TA::einsum("ijab,cbij->caij", t, T2);
        }();
        REQUIRE(equal_tarrays(eval_rdm, man_rdm, "a2,a1,i1,i2"));
      }

      {  // multiple bra or ket indices require AssertStrictBraKetSymmetry::No
        auto ctx_resetter = sequant::set_scoped_default_context(
            sequant::Context{sequant::get_default_context()}.set(
                sequant::AssertStrictBraKetSymmetry::No));

        // hyperindex i1 in ket slots of 3 tensors
        auto expr2 = deserialize(L"T{a1;i1} T{a2;i1} T{a3;i1}");
        auto eval2 = eval(expr2, "a_1,a_2,a_3");
        auto man2 = [&]() {
          auto T1 = yield(L"T{a1;i1}");
          auto T2 = yield(L"T{a2;i1}");
          auto T3 = yield(L"T{a3;i1}");
          auto T12 = TA::einsum("ai,bi->abi", T1, T2);
          return TA::einsum("abi,ci->abc", T12, T3);
        }();
        REQUIRE(equal_tarrays(eval2, man2, "a1,a2,a3"));

        // hyperindex a1 in bra slots of 3 tensors
        auto expr3 = deserialize(L"T{a1;i1} T{a1;i2} T{a1;i3}");
        auto eval3 = eval(expr3, "i_1,i_2,i_3");
        auto man3 = [&]() {
          auto T1 = yield(L"T{a1;i1}");
          auto T2 = yield(L"T{a1;i2}");
          auto T3 = yield(L"T{a1;i3}");
          auto T12 = TA::einsum("ai,aj->aij", T1, T2);
          return TA::einsum("aij,ak->ijk", T12, T3);
        }();
        REQUIRE(equal_tarrays(eval3, man3, "i1,i2,i3"));
      }

      {  // high-order aux hyperindex carried into the result: aux index x1 is
         // shared by 3 tensors AND external (in the result aux slot). This is
         // the >2-tensor named-hyperedge case the TNv3 canonicalize_slots fix
         // unlocks; TA evaluates the resulting batched-over-x1 contraction
         // natively (einsum keeps an index common to both operands and result).
         //   R{;;x1} = A{;a1;x1} B{a1;a2;x1} C{a2;;x1}
         //   R[x] = sum_{a1,a2} A[a1,x] B[a1,a2,x] C[a2,x]
        auto res = deserialize<sequant::ResultExpr>(
            L"R{;;x1} = A{;a1;x1} B{a1;a2;x1} C{a2;;x1}");
        auto evalR = evaluate(eval_node(res), std::string("x_1"), yield_)
                         ->get<TA::TArrayD>();
        auto manR = [&]() {
          auto A = yield(L"A{;a1;x1}");    // [a1, x]
          auto B = yield(L"B{a1;a2;x1}");  // [a1, a2, x]
          auto C = yield(L"C{a2;;x1}");    // [a2, x]
          auto AB = TA::einsum("ax,abx->bx", A, B);
          return TA::einsum("bx,bx->x", AB, C);
        }();
        REQUIRE(equal_tarrays(evalR, manR, "x1"));
      }
    }  // multiple bra or ket indices
  }

  SECTION("complex") {
    using TArrayC = TA::DistArray<TA::Tensor<std::complex<double>>>;

    const size_t nocc = 2, nvirt = 20;
    auto& world = TA::get_default_world();

    auto yield_ = rand_tensor_yield<std::complex<double>, TA::DensePolicy>{
        world, nocc, nvirt};

    auto yield = [&yield_](std::wstring_view lbl) -> TArrayC const& {
      return yield_(lbl)->get<TArrayC>();
    };

    auto eval = [&yield_](sequant::ExprPtr const& expr,
                          std::string const& target_labels) {
      return evaluate(eval_node(expr), target_labels, yield_)->get<TArrayC>();
    };

    auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                               std::string const& target_labels) {
      return evaluate_symm(eval_node(expr), target_labels, yield_)
          ->get<TArrayC>();
    };

    auto eval_antisymm = [&yield_](sequant::ExprPtr const& expr,
                                   std::string const& target_labels) {
      return evaluate_antisymm(eval_node(expr), target_labels, yield_)
          ->get<TArrayC>();
    };

    using namespace sequant;
    using namespace std::string_literals;

    SECTION("summation") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayC{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum1_eval, sum1_man));

      auto expr2 =
          deserialize<sequant::ExprPtr>(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayC{};
      sum2_man("i1,a1") =
          std::complex<double>{2} * yield(L"t{a1;i1}")("a1,i1") +
          std::complex<double>{1.5} * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("product") {
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayC prod1_man{};
      prod1_man("i4,a1,a4,i1") = std::complex<double>{1 / 2.0} *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(equal_tarrays(prod1_eval, prod1_man));

      auto expr2 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayC{};
      prod2_man("a1,a2,i1,i2") = std::complex<double>{-1 / 4.0} *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      REQUIRE(equal_tarrays(prod2_eval, prod2_man));

      auto expr3 = sequant::deserialize<sequant::ExprPtr>(
          L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayC{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");

      REQUIRE(equal_tarrays(prod3_eval, prod3_man));

      auto expr4 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayC{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(equal_tarrays(prod4_eval, prod4_man));
    }

    SECTION("sum and product") {
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
      auto eval1 = eval(expr1, "a_1,a_2,i_1,i_2");

      auto man1 = TArrayC{};
      man1("a1,a2,i1,i2") = std::complex<double>{-1.0 / 4} *
                                yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4") +
                            std::complex<double>{1.0 / 16} *
                                yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a1,a2;i3,i4}")("a1,a2,i3,i4") *
                                yield(L"t{a3,a4;i1,i2}")("a3,a4,i1,i2");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayC{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Antisymmetrization") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = std::complex<double>{0.25} * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // odd-ranked tensor
      auto expr2 = deserialize<sequant::ExprPtr>(L"g_{i1, i2, i3}^{a1, a2}");
      auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2}");

      auto man2 = TArrayC{};
      man2("0,1,2,3,4") =
          arr2("0,1,2,3,4") - arr2("1,0,2,3,4") + arr2("1,2,0,3,4") -
          arr2("2,1,0,3,4") + arr2("2,0,1,3,4") - arr2("0,2,1,3,4") -
          arr2("0,1,2,4,3") + arr2("1,0,2,4,3") - arr2("1,2,0,4,3") +
          arr2("2,1,0,4,3") - arr2("2,0,1,4,3") + arr2("0,2,1,4,3");

      REQUIRE(equal_tarrays(eval2, man2));

      auto expr3 = deserialize<sequant::ExprPtr>(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayC{};
      man3("0,1") = arr3("0,1") - arr3("1,0");
      man3("0,1") = std::complex<double>{0.5} * man3("0,1");

      REQUIRE(equal_tarrays(eval3, man3));
    }

    SECTION("Symmetrization") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = deserialize<sequant::ExprPtr>(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayC man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");
      man2("0,1,2,3,4,5") = (1.0 / 6.0) * man2("0,1,2,3,4,5");

      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Others") {
      using namespace std::string_literals;
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

      auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                       ->get<TArrayC>();

      auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
        return eval_node(x);
      }) | ranges::to_vector;

      auto eval2 = evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TArrayC>();

      REQUIRE(equal_tarrays(eval1, eval2));
    }
  }

  SECTION("tot") {
    using namespace sequant;

    //
    // eg: approx_equal("i,j;a,b", arr1, arr2)
    // - arr1 and arr2 are DistArrays with equal TiledRange and matching Range
    // for
    //   the inner tensors at the corresponding tile positions.
    // - 'i', 'j', 'a', and 'b' are dummy indices that annotate the modes of
    // outer
    //   and inner tensors. Why? Because TA::norm2 function is not supported for
    //   tensor-of-tensor tiles
    //
    auto approx_equal = [](std::string const& annot, auto const& lhs,
                           auto const& rhs) -> bool {
      return Catch::Approx(lhs(annot).dot(lhs(annot))) ==
             rhs(annot).dot(rhs(annot));
    };

    auto& world = TA::get_default_world();

    size_t const nocc = 2;
    size_t const nvirt = 3;

    rand_tensor_yield<int> yield{world, nocc, nvirt};

    using ArrayT = typename decltype(yield)::array_type;
    using ArrayToT = typename decltype(yield)::array_tot_type;
    using NumericT = typename decltype(yield)::numeric_type;

    SECTION("T_times_ToT_to_ToT") {
      constexpr std::wstring_view expr_str =
          L"3"
          L" * "
          L"f{i3;i1}"
          L" * "
          L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}";
      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));
      std::string const target_layout{"i_1,i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
      auto result = evaluate(node, target_layout, yield)->get<ArrayToT>();
      ArrayToT ref;
      {
        auto const& lhs = yield(L"f{i3;i1}")->get<ArrayT>();
        auto const& rhs =
            yield(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}")->get<ArrayToT>();
        ref = TA::einsum(lhs("i_3,i_1"), rhs("i_2,i_3;a_3i_2i_3,a_4i_2i_3"),
                         target_layout);
        ref(target_layout) = 3 * ref(target_layout);
      }
      REQUIRE(approx_equal("i,j,k;a,b", result, ref));
    }

    SECTION("ToT_times_ToT_to_ToT") {
      constexpr std::wstring_view expr_str =
          L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}"
          L" * "
          L"s{a2<i1,i2>;a4<i2,i3>}";

      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));
      std::string const target_layout{"i_2,i_1;a_1i_1i_2,a_2i_1i_2"};

      auto result = evaluate(node, target_layout, yield)->get<ArrayToT>();

      ArrayToT ref;
      {
        auto const& lhs =
            yield(L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}")->get<ArrayToT>();
        auto const& rhs = yield(L"s{a2<i1,i2>;a4<i2,i3>}")->get<ArrayToT>();
        ref = TA::einsum(lhs("i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"),
                         rhs("i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"), target_layout);
      }
      REQUIRE(approx_equal("i,j;a,b", result, ref));
    }

    SECTION("ToT_times_ToT_to_Scalar") {
      constexpr std::wstring_view expr_str =
          L"I{a1<i1,i2>,a2<i1,i2>;i1,i2}"
          L" * "
          L"g{i1,i2;a2<i1,i2>,a1<i1,i2>}";
      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));

      auto result = evaluate(node, yield)->get<NumericT>();

      NumericT ref;
      {
        auto const& lhs =
            yield(L"I{a1<i1,i2>,a2<i1,i2>;i1,i2}")->get<ArrayToT>();
        auto const& rhs =
            yield(L"g{i1,i2;a2<i1,i2>,a1<i1,i2>}")->get<ArrayToT>();
        ref = TA::dot(lhs("i_1,i_2;a_1i_1i_2,a_2i_1i_2"),
                      rhs("i_1,i_2;a_2i_1i_2,a_1i_1i_2"));
      }
      REQUIRE(result == Catch::Approx(ref));
    }

    SECTION("symmetrize") {
      // A double-valued yield is used because the 1/n! prefactor is fractional.
      rand_tensor_yield<double> dyield{world, nocc, nvirt};
      using DArrayToT = typename decltype(dyield)::array_tot_type;

      // n=2 (logical rank 4): outer = occupied (i1,i2), inner = virtual
      // (a1,a2). symmetrize() must produce S(i,j;a,b) = 1/2! * (R(i,j;a,b) +
      // R(j,i;b,a))
      // -- outer and inner modes permute in lockstep.
      {
        auto const Rnode = eval_node(
            deserialize<sequant::ExprPtr>(L"R{a1<i1,i2>,a2<i1,i2>;i1,i2}"));
        auto const Rres = dyield(Rnode);
        auto const& R = Rres->get<DArrayToT>();

        auto const symm = Rres->symmetrize()->get<DArrayToT>();

        std::string const annot{"i_1,i_2;a_1i_1i_2,a_2i_1i_2"};
        DArrayToT ref;
        ref(annot) = R(annot) + R("i_2,i_1;a_2i_1i_2,a_1i_1i_2");
        ref(annot) = 0.5 * ref(annot);

        DArrayToT diff;
        diff(annot) = symm(annot) - ref(annot);
        // Norm-squared of the difference; an absolute margin is needed because
        // Approx(0.0) is a pure relative test. The reference sums the same
        // permutations in a different order than the production helper, so a
        // tiny rounding residual (~1e-30) is expected; a real mismatch is O(1).
        REQUIRE(diff(annot).dot(diff(annot)) ==
                Catch::Approx(0.0).margin(1e-12));
      }

      // n=3 (logical rank 6): exercises the 6-permutation / 1/3! path. The 6
      // particle permutations permute outer (i1,i2,i3) and inner (a1,a2,a3) in
      // lockstep; each inner index keeps its proto-suffix "i_1i_2i_3".
      {
        auto const Rnode = eval_node(deserialize<sequant::ExprPtr>(
            L"R{a1<i1,i2,i3>,a2<i1,i2,i3>,a3<i1,i2,i3>;i1,i2,i3}"));
        auto const Rres = dyield(Rnode);
        auto const& R = Rres->get<DArrayToT>();

        auto const symm = Rres->symmetrize()->get<DArrayToT>();

        std::string const annot{
            "i_1,i_2,i_3;a_1i_1i_2i_3,a_2i_1i_2i_3,a_3i_1i_2i_3"};
        DArrayToT ref;
        ref(annot) = R("i_1,i_2,i_3;a_1i_1i_2i_3,a_2i_1i_2i_3,a_3i_1i_2i_3") +
                     R("i_1,i_3,i_2;a_1i_1i_2i_3,a_3i_1i_2i_3,a_2i_1i_2i_3") +
                     R("i_3,i_1,i_2;a_3i_1i_2i_3,a_1i_1i_2i_3,a_2i_1i_2i_3") +
                     R("i_3,i_2,i_1;a_3i_1i_2i_3,a_2i_1i_2i_3,a_1i_1i_2i_3") +
                     R("i_2,i_3,i_1;a_2i_1i_2i_3,a_3i_1i_2i_3,a_1i_1i_2i_3") +
                     R("i_2,i_1,i_3;a_2i_1i_2i_3,a_1i_1i_2i_3,a_3i_1i_2i_3");
        ref(annot) = (1.0 / 6.0) * ref(annot);

        DArrayToT diff;
        diff(annot) = symm(annot) - ref(annot);
        // Norm-squared of the difference; an absolute margin is needed because
        // Approx(0.0) is a pure relative test. The reference sums the same
        // permutations in a different order than the production helper, so a
        // tiny rounding residual (~1e-30) is expected; a real mismatch is O(1).
        REQUIRE(diff(annot).dot(diff(annot)) ==
                Catch::Approx(0.0).margin(1e-12));
      }
    }
  }
}

TEST_CASE("eval_custom_evaluator", "[eval]") {
  using sequant::evaluate;
  using sequant::ResultPtr;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 20;
  auto yield_ = rand_tensor_yield<double, TA::DensePolicy>{world, nocc, nvirt};

  // a multi-product expression: several non-leaf nodes in the eval tree.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}",
      {.def_perm_symm = sequant::Symmetry::Antisymm});
  std::string const target = "a_1,a_2,i_1,i_2";
  auto const node = eval_node(expr);

  // standard-scheme reference (no custom evaluator)
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  SECTION("declining evaluator defers to the standard scheme") {
    // A custom evaluator that always returns null is consulted at every
    // non-leaf node and the standard scheme produces the result.
    int consulted = 0;
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(
        [&consulted](node_t const&, cache_t&) -> ResultPtr {
          ++consulted;
          return nullptr;
        });
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(equal_tarrays(res, ref));
    REQUIRE(consulted > 1);  // multiple non-leaf nodes, all declined
  }

  SECTION("intercepting evaluator takes over a subtree") {
    // A custom evaluator that takes over the first (root) node it is consulted
    // on -- here by re-evaluating that subtree via the standard scheme on a
    // scratch cache. The result must still match, and the non-null return must
    // short-circuit the recursion, so the evaluator fires exactly once.
    int consulted = 0;
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(
        [&consulted, &yield_](node_t const& n, cache_t&) -> ResultPtr {
          ++consulted;
          auto scratch = cache_t::empty();
          return evaluate(n, yield_, scratch);
        });
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(equal_tarrays(res, ref));
    REQUIRE(consulted == 1);
  }
}

TEST_CASE("eval_batch_axis", "[eval]") {
  using sequant::batch_axis;
  using sequant::contracted_indices;

  auto node_of = [](std::wstring_view xpr) {
    return eval_node(sequant::deserialize<sequant::ExprPtr>(xpr));
  };

  SECTION("single contracted index") {
    // R_{a1}^{i1,i3} * f_{i3}^{i2} sums over i3.
    auto const node = node_of(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 1);
    auto const mode = batch_axis(node);
    REQUIRE(mode.has_value());
    REQUIRE(mode.value() == c.front());
  }

  SECTION("two contracted indices") {
    // g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4} sums over a1,a2.
    auto const node = node_of(L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 2);
    auto const mode = batch_axis(node);
    REQUIRE(mode.has_value());
    REQUIRE(ranges::contains(c, mode.value()));
  }

  SECTION("leaf has no contracted index") {
    auto const node = node_of(L"f_{i1}^{a1}");
    REQUIRE(contracted_indices(node).empty());
    REQUIRE_FALSE(batch_axis(node).has_value());
  }

  SECTION("a sum is not a contraction") {
    auto const node = node_of(L"f_{i1}^{a1} + t_{a1}^{i1}");
    REQUIRE(contracted_indices(node).empty());
    REQUIRE_FALSE(batch_axis(node).has_value());
  }

  SECTION("predicate scopes the batch mode") {
    // contracts a1,a2 (unoccupied)
    auto const node = node_of(L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 2);

    // accept exactly one contracted index -> batch_axis returns it
    auto const only_first = [&c](sequant::Index const& ix) {
      return ix == c.front();
    };
    REQUIRE(batch_axis(node, only_first) == c.front());

    // accept none -> nullopt
    REQUIRE_FALSE(batch_axis(node, [](sequant::Index const&) {
                    return false;
                  }).has_value());

    // scope batching to a specific IndexSpace
    auto const unocc = c.front().space();
    auto const in_unocc = [&unocc](sequant::Index const& ix) {
      return ix.space() == unocc;
    };
    REQUIRE(batch_axis(node, in_unocc).has_value());

    // a node whose only contracted index is in a different (occupied) space
    auto const node_occ = node_of(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");  // sums i3
    REQUIRE_FALSE(batch_axis(node_occ, in_unocc).has_value());
  }
}

TEST_CASE("eval_slice_array_over_mode", "[eval]") {
  using sequant::slice_array_over_mode;
  auto& world = TA::get_default_world();

  // arr: a(2 tiles), b(3 tiles), c(1 tile); bb shares b's TiledRange1.
  TA::TArrayD arr(world, TA::TiledRange{{0, 2, 4}, {0, 3, 6, 9}, {0, 5}});
  TA::TArrayD bb(world, TA::TiledRange{{0, 3, 6, 9}, {0, 7}});
  arr.fill_random();
  bb.fill_random();
  world.gop.fence();

  SECTION("trange of the sliced mode") {
    auto const s = slice_array_over_mode(arr, 1, 1, 3);  // b-tiles [1,3)
    REQUIRE(s.trange().dim(1).tile_extent() == 2);
    REQUIRE(s.trange().dim(0).tile_extent() ==
            arr.trange().dim(0).tile_extent());
    REQUIRE(s.trange().dim(2).tile_extent() ==
            arr.trange().dim(2).tile_extent());
  }

  SECTION(
      "blocked contraction over the sliced mode reconstructs the full one") {
    // full contraction over b
    TA::TArrayD full;
    full("a,c,d") = arr("a,b,c") * bb("b,d");

    // split b's 3 tiles into [0,1) and [1,3), contract each, sum
    auto const a0 = slice_array_over_mode(arr, 1, 0, 1);
    auto const a1 = slice_array_over_mode(arr, 1, 1, 3);
    auto const b0 = slice_array_over_mode(bb, 0, 0, 1);
    auto const b1 = slice_array_over_mode(bb, 0, 1, 3);
    TA::TArrayD p0, p1, summed;
    p0("a,c,d") = a0("a,b,c") * b0("b,d");
    p1("a,c,d") = a1("a,b,c") * b1("b,d");
    summed("a,c,d") = p0("a,c,d") + p1("a,c,d");

    REQUIRE(equal_tarrays(summed, full));
  }

  SECTION("Result::slice_mode takes tile-aligned element bounds") {
    // mode 1 (b) tiles {0,3,6,9}; element range [3,9) (tile-aligned, as
    // mode_batches produces) corresponds to tiles [1,3).
    sequant::ResultPtr const r =
        sequant::eval_result<sequant::ResultTensorTA<TA::TArrayD>>(arr);
    auto const via_result = r->slice_mode(1, 3, 9)->get<TA::TArrayD>();
    auto const direct = slice_array_over_mode(arr, 1, 1, 3);
    REQUIRE(equal_tarrays(via_result, direct));
  }

  SECTION("Result::mode_batches partitions a mode into element ranges") {
    using batches_t = sequant::container::svector<std::pair<size_t, size_t>>;
    sequant::ResultPtr const r =
        sequant::eval_result<sequant::ResultTensorTA<TA::TArrayD>>(arr);
    // mode 1 (b) has 3 tiles of 3 elements each (extent 9). target_batch_size
    // is an UPPER BOUND: each batch is the largest whole-tile group whose total
    // size does not exceed the target, with a floor of one tile (so a target
    // below the tile size still yields one tile per batch). A batch must never
    // exceed the target except via that one-tile floor.
    // (extra parens: compare as a single bool so Catch2 needn't stringify
    // pairs) target >= extent -> a single batch (whole mode).
    REQUIRE((r->mode_batches(1, 100) == batches_t{{0, 9}}));
    // target == 2 tiles -> two-tile batches (6 <= 6).
    REQUIRE((r->mode_batches(1, 6) == batches_t{{0, 6}, {6, 9}}));
    // target just above the tile size (but below 2 tiles) -> ONE tile per
    // batch: a 2-tile batch (6) would exceed the target. Regression: the old
    // `acc >= target` rule rounded UP to a 2-tile batch, so any target a hair
    // above the tile size doubled the realized batch -- the aux_target_size
    // 236->243 crash (236-wide K tiles, 243 target -> 472-wide batch).
    REQUIRE((r->mode_batches(1, 5) == batches_t{{0, 3}, {3, 6}, {6, 9}}));
    REQUIRE((r->mode_batches(1, 4) == batches_t{{0, 3}, {3, 6}, {6, 9}}));
    // target == tile size -> one tile per batch.
    REQUIRE((r->mode_batches(1, 3) == batches_t{{0, 3}, {3, 6}, {6, 9}}));
    // target below the tile size -> still one tile per batch (the floor).
    REQUIRE((r->mode_batches(1, 1) == batches_t{{0, 3}, {3, 6}, {6, 9}}));
  }
}

// Numeric-invariance for Result::write_into_slice on the TA backend: the union
// of disjoint, tile-aligned blocks scattered into a pre-sized destination must
// reconstruct the whole array EXACTLY. write_into_slice() is the inverse of
// slice_array_over_mode() (which GATHERS a block out): here we gather disjoint
// blocks of a whole array R and scatter each back into a fresh, pre-sized
// destination, then require the reassembled destination == R elementwise.
TEST_CASE("eval_write_into_slice", "[eval]") {
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultTensorOfTensorTA;
  using sequant::ResultTensorTA;
  using sequant::slice_array_over_mode;
  auto& world = TA::get_default_world();

  SECTION("flat: disjoint tile-aligned blocks reassemble exactly") {
    // mode 0 has 3 multi-element tiles (size 3 each; occ_tile_size>1 analog),
    // mode 1 a single tile. Split mode 0 into element ranges [0,3) and [3,9)
    // (tiles [0,1) and [1,3)): disjoint, tile-aligned, and tile-SPANNING.
    TA::TArrayD R(world, TA::TiledRange{{0, 3, 6, 9}, {0, 5}});
    R.fill_random();
    world.gop.fence();

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // elements [0,3)
    auto const b1 = slice_array_over_mode(R, 0, 1, 3);  // elements [3,9)

    // destination pre-sized to R's full shape, then filled block by block.
    TA::TArrayD dest(world, R.trange());
    dest.fill_local(0.0);
    world.gop.fence();

    auto rdest = eval_result<ResultTensorTA<TA::TArrayD>>(dest);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TArrayD>>(b0), 0, 0,
                            3);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TArrayD>>(b1), 0, 3,
                            9);

    REQUIRE(equal_tarrays<Tight>(rdest->get<TA::TArrayD>(), R));
  }

  SECTION("flat: nonzero element lobound (frozen-core offset) is preserved") {
    // mode 0 has element lobound 2 (a frozen-core-like offset) and two size-3
    // tiles; split into [2,5) and [5,8). The block bounds honor the lobound.
    TA::TArrayD R(world, TA::TiledRange{{2, 5, 8}, {0, 4}});
    R.fill_random();
    world.gop.fence();

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // elements [2,5)
    auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // elements [5,8)
    // the gathered blocks keep the source element lobound
    REQUIRE(b0.trange().dim(0).elements_range().first == 2);

    TA::TArrayD dest(world, R.trange());
    dest.fill_local(0.0);
    world.gop.fence();

    auto rdest = eval_result<ResultTensorTA<TA::TArrayD>>(dest);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TArrayD>>(b0), 0, 2,
                            5);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TArrayD>>(b1), 0, 5,
                            8);

    auto const& out = rdest->get<TA::TArrayD>();
    REQUIRE(out.trange().dim(0).elements_range().first == 2);  // lobound kept
    REQUIRE(equal_tarrays<Tight>(out, R));
  }

  SECTION("tot: disjoint tile-aligned blocks reassemble exactly") {
    using ToTArray = TA::DistArray<TA::Tensor<TA::Tensor<double>>>;
    using ResultToT = ResultTensorOfTensorTA<ToTArray>;

    // outer mode 0: two size-3 tiles (multi-element, tile-spanning split);
    // outer mode 1: one size-2 tile. Inner tensors are 2x2. Inner values are
    // position-dependent (derived from the outer coordinate) so a mis-scattered
    // block -- wrong offset -- would change the reassembled norm.
    TA::TiledRange const outer_tr{{0, 3, 6}, {0, 2}};
    TA::Range const inner_r(std::array<std::size_t, 2>{2, 2});

    auto build = [&world, inner_r](TA::TiledRange const& otr,
                                   bool zero) -> ToTArray {
      auto tile_fn = [inner_r, zero](TA::Range const& orng) {
        TA::Tensor<TA::Tensor<double>> t{orng};
        std::size_t o = 0;
        for (auto const& coord : orng) {
          auto& inner = t[o++];
          inner = TA::Tensor<double>{inner_r};
          double base = 0.0;
          if (!zero)
            for (auto c : coord)
              base = base * 37.0 + static_cast<double>(c + 1);
          std::size_t k = 0;
          for (auto& x : inner) x = zero ? 0.0 : base * 100.0 + (++k);
        }
        return t;
      };
      ToTArray arr{world, otr};
      for (auto it = arr.begin(); it != arr.end(); ++it)
        if (arr.is_local(it.index()))
          *it = world.taskq.add(tile_fn, it.make_range());
      world.gop.fence();
      return arr;
    };

    ToTArray R = build(outer_tr, /*zero=*/false);

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // outer elements [0,3)
    auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // outer elements [3,6)

    // destination pre-sized to R's full shape with well-formed (zero) inners.
    ToTArray dest = build(outer_tr, /*zero=*/true);

    auto rdest = eval_result<ResultToT>(dest);
    rdest->write_into_slice(*eval_result<ResultToT>(b0), 0, 0, 3);
    rdest->write_into_slice(*eval_result<ResultToT>(b1), 0, 3, 6);

    // ToT elementwise invariance: norm of the difference (dot of the diff with
    // itself) must vanish. TA::norm2 does not support ToT tiles, so use dot.
    auto const& out = rdest->get<ToTArray>();
    ToTArray diff;
    diff("i,j;a,b") = out("i,j;a,b") - R("i,j;a,b");
    REQUIRE(Catch::Approx(diff("i,j;a,b").dot(diff("i,j;a,b"))) == 0.0);
    // and the reassembled result is nonzero (guards against a trivial all-zero
    // pass, e.g. if both blocks silently no-op'd).
    REQUIRE(out("i,j;a,b").dot(out("i,j;a,b")) > 0.0);
  }
}

TEST_CASE("eval_batched_custom_evaluator", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // multi-tile arrays so batching over a contracted index actually engages:
  // unoccupied extent 12 in tiles of <=4 -> 3 tiles.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // contracts a1,a2 (unoccupied) -> batch mode is an unoccupied index (3 tiles)
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto node = eval_node(expr);

  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node);
  REQUIRE(ax.has_value());
  node->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});

  // Reference first, so yield_'s (random) leaf arrays are generated and cached;
  // the batched evaluator below copies yield_ and thus reuses the same arrays.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Batched evaluation must reproduce the reference for any target batch size,
  // since sum_K = sum_{batches} sum_{K in batch}. The batch mode is unoccupied
  // (extent 12, tiles of 4 -> 3 tiles). target_batch_size is in *elements*:
  // 100 -> 1 batch (no-op), 8 -> 2 batches ([0,8),[8,12)), 4 -> 3 batches, and
  // 1 -> 3 batches (each tile its own batch).
  for (std::size_t target_batch_size :
       {std::size_t{100}, std::size_t{8}, std::size_t{4}, std::size_t{1}}) {
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        }));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    // batched summation reorders the contraction, so allow a looser FP margin
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }
}

TEST_CASE("eval_batched_custom_evaluator persistence gate", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // Contracts a1,a2 (unoccupied, 3 tiles) -> batchable over an unoccupied mode.
  // The subtree contains a "t" leaf, which we treat as volatile.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto node = eval_node(expr);
  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node);
  REQUIRE(ax.has_value());
  node->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Volatile-leaf predicate: the amplitude "t".
  auto is_volatile_t = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  // (1) Gate ON (persistent_only=true): the subtree contains a volatile "t"
  // leaf, so batching is DECLINED -- the spy scope-guard (invoked only once a
  // node passes every gate and yields >1 batch) is never called, yet the
  // standard scheme still gives the correct result.
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
        sequant::accept_any_index{}, spy, is_volatile_t,
        /*persistent_only=*/true));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE_FALSE(batched);  // gated volatile subtree -> not batched
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }

  // (2) No gate (default never_volatile): the SAME node DOES batch -- confirms
  // (1)'s decline is due to the volatility gate, not the index/tiling setup.
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
        sequant::accept_any_index{}, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(batched);  // no gate -> batched as before
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }
}

TEST_CASE("eval_batched_custom_evaluator_tot", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Multi-tile the occupied space so the contracted occupied index i3 spans
  // more than one tile and batching over it actually engages: extent 8 in
  // tiles of <=4 -> 2 tiles. The inner (virtual) space is left single-tiled.
  rand_tensor_yield<double> yield{world, /*nocc=*/8, /*nvirt=*/3};
  yield.set_max_tile(4);

  using ArrayToT = typename decltype(yield)::array_tot_type;

  // ToT * ToT -> ToT (the same expression as the "tot" section above). The
  // contracted indices are the occupied i3 (an *outer* mode of both leaves)
  // and the inner virtual a4. Scoping the batch mode to the occupied space
  // selects i3, so the batched partials slice ToT leaves over i3
  // (slice_array_over_mode) and sum ToT partials (add_inplace) -- the two
  // annotation-free ToT array operations that must emit an "outer;inner"
  // annotation rather than a flat one (else DistArray's is_tot_index() trips).
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"I{a4<i2,i3>,a1<i1,i2>;i1,i2} * s{a2<i1,i2>;a4<i2,i3>}");
  std::string const target = "i_2,i_1;a_1i_1i_2,a_2i_1i_2";
  auto node = eval_node(expr);
  // Reference first (non-batched), so yield's random leaf arrays are generated
  // and cached; the batched evaluator reuses the same arrays.
  auto const ref = evaluate(node, target, yield)->get<ArrayToT>();

  auto const occ =
      sequant::get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node, accept_occ);
  REQUIRE(ax.has_value());
  node->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});

  // TA::norm2 is unsupported for tensor-of-tensor tiles, so compare via the
  // self-dot of each array (a scalar norm^2); reordering the contraction over
  // i3 must not change it.
  auto self_dot = [](auto const& arr) {
    return arr("i,j;a,b").dot(arr("i,j;a,b"));
  };
  auto const ref_dot = self_dot(ref);

  for (std::size_t target_batch_size :
       {std::size_t{100}, std::size_t{4}, std::size_t{1}}) {
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_occ));
    auto const res = evaluate(node, target, yield, cache)->get<ArrayToT>();
    REQUIRE(self_dot(res) == Catch::Approx(ref_dot));
  }
}

TEST_CASE("ta_tot_conj_complex", "[eval]") {
  // Exercises TiledArray's .conj() on a complex tensor-of-tensors through the
  // expression engine — the capability SeQuant's ToT adjoint() relies on after
  // teaching TA to recurse conj into nested tiles. With genuinely complex data
  // a missing/incorrect inner conj would flip the wrong imaginary sign.
  using namespace sequant;
  auto& world = TA::get_default_world();
  size_t const nocc = 2, nvirt = 3;
  rand_tensor_yield<std::complex<double>, TA::DensePolicy> yield{world, nocc,
                                                                 nvirt};
  using ArrayToT = typename decltype(yield)::array_tot_type;

  std::string const annot{"i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
  auto const t = deserialize<sequant::ExprPtr>(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}");
  auto const& src = yield(t->as<sequant::Tensor>())->get<ArrayToT>();

  ArrayToT conjd;
  conjd(annot) = src(annot).conj();
  ArrayToT::wait_for_lazy_cleanup(world);

  // elementwise: conjd inner == conj(src inner); identical tiling/layout (no
  // permutation in the assignment) lets us walk local tiles in lockstep.
  auto it_s = src.begin();
  auto it_c = conjd.begin();
  for (; it_s != src.end(); ++it_s, ++it_c) {
    auto const& souter = it_s->get();
    auto const& couter = it_c->get();
    REQUIRE(souter.size() == couter.size());
    for (std::size_t o = 0; o < souter.size(); ++o) {
      auto const& sinner = souter[o];
      auto const& cinner = couter[o];
      if (sinner.empty()) continue;
      for (std::size_t k = 0; k < sinner.size(); ++k) {
        CHECK(cinner[k].real() == Catch::Approx(sinner[k].real()));
        CHECK(cinner[k].imag() == Catch::Approx(-sinner[k].imag()));
      }
    }
  }
}

TEST_CASE("eval_batched_scratch", "[eval]") {
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using sequant::Index;

  // W-analog of the PNO-CCSD PPL intermediate: two canonically-equal internal
  // siblings, both carrying the auxiliary batch mode x_1 (free at the
  // children, contracted at the root; an aux-aux edge, like the DF index K).
  // Every orbital contraction pairs a bra with a ket.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (g{a_3;i_2;x_1} * h{i_4;a_3})");
  auto const node = eval_node(expr);
  REQUIRE_FALSE(node.leaf());
  REQUIRE_FALSE(node.left().leaf());
  REQUIRE_FALSE(node.right().leaf());
  // structural precondition: the two children ARE canonically equal
  {
    auto cm = sequant::cache_manager(ranges::views::single(node), 2);
    REQUIRE(cm.max_life(node.left()) == 2);
  }
  auto const x1 = [&] {  // the contracted (batch) mode
    auto modes = sequant::contracted_indices(node);
    REQUIRE(modes.size() == 1);
    return modes[0];
  }();

  auto real = cache_t::empty();

  SECTION("registers consistent repeats with in-pass counts") {
    std::vector<std::pair<node_t const*, Index>> const members{{&node, x1}};
    auto bs = sequant::detail::make_batched_scratch(members, real);
    REQUIRE(bs.cache.exists(node.left()));
    REQUIRE(bs.cache.max_life(node.left()) == 2);
    REQUIRE_FALSE(bs.cache.exists(node));  // member roots are not registered
    REQUIRE(bs.seeds.empty());
  }

  SECTION("signature-inconsistent subnodes are not registered") {
    // the same subtree appears under two members, but the second member's
    // mode (an index the shared subnode does not carry) gives it signature
    // 'absent' while the first gives a position -> inconsistent -> unshared
    auto const expr2 = sequant::deserialize<sequant::ExprPtr>(
        L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * p{i_5;i_6;x_1}");
    auto const node2 = eval_node(expr2);
    auto const bogus_axis = Index(L"i_9");
    std::vector<std::pair<node_t const*, Index>> const members{
        {&node, x1}, {&node2, bogus_axis}};
    auto bs = sequant::detail::make_batched_scratch(members, real);
    REQUIRE_FALSE(bs.cache.exists(node.left()));
  }

  SECTION("descends through inconsistently-sliced re-encounters") {
    // M1 (mode x_1) = X * D2 with X = (D * u): D's two occurrences (inside X
    // and as the root's sibling D2) are visited with the same signature, so D
    // alone would be registered. M2 (bogus mode) re-encounters X with
    // signature 'absent' -- inconsistent, unshared. The walk must descend
    // through that re-encounter: under it D's signature is also 'absent', so
    // sharing D would serve M2's (per-occurrence) evaluation of X a wrongly
    // sliced value. D must end up unregistered.
    auto const expr_m1 = sequant::deserialize<sequant::ExprPtr>(
        L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * u{i_5;i_3}) * "
        L"(g{a_3;i_2;x_1} * h{i_4;a_3})");
    auto const m1 = eval_node(expr_m1);
    auto const expr_m2 = sequant::deserialize<sequant::ExprPtr>(
        L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * u{i_5;i_3}) * p{i_6;i_7;x_1}");
    auto const m2 = eval_node(expr_m2);
    // structural preconditions: m1 = X * D2, X = D * u, D2 == D == m2's X
    // child canonically
    sequant::TreeNodeEqualityComparator<node_t> const eq;
    REQUIRE_FALSE(m1.left().leaf());
    REQUIRE_FALSE(m1.left().left().leaf());
    auto const& D = m1.left().left();
    REQUIRE(eq(D, m1.right()));
    REQUIRE(eq(m1.left(), m2.left()));

    // positive control: M1 alone registers D (count 2, consistent signature)
    {
      std::vector<std::pair<node_t const*, Index>> const members{{&m1, x1}};
      auto bs = sequant::detail::make_batched_scratch(members, real);
      REQUIRE(bs.cache.exists(D));
      REQUIRE(bs.cache.max_life(D) == 2);
    }
    // adding M2 makes X inconsistent AND must expose D's 'absent' signature
    // beneath X's second occurrence
    {
      auto const bogus_axis = Index(L"i_9");
      std::vector<std::pair<node_t const*, Index>> const members{
          {&m1, x1}, {&m2, bogus_axis}};
      auto bs = sequant::detail::make_batched_scratch(members, real);
      REQUIRE_FALSE(bs.cache.exists(m1.left()));  // X: inconsistent
      REQUIRE_FALSE(bs.cache.exists(D));  // D: inconsistent via pruned branch
    }
  }
}

// BLOCKED, hidden by default -- see .superpowers/sdd/oamb-a0-note.md sections
// 10.4 and 11. This test batched an UNANNOTATED node via the runtime's depth-0
// heuristic fallback, which has been removed (annotations are now
// authoritative). Re-pointing it at an explicit annotation is NOT
// behaviour-preserving: cache_manager vetoes caching for a node whose own
// batched_here() carries a sliced batchable mode, and the heuristic never set
// batched_here() -- so the same mode batched by the two routes gives different
// CSE. The correct expectations depend on Phase B replacing that veto with
// per-context (per-slice) caching. Re-enable and re-derive the counts then;
// do not "fix" the numbers against veto behaviour that is about to be deleted.

// Task 3 of the whole-scope batched DAG execution design
// (doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md):
// EQUIVALENCE. A two-root forest whose roots SHARE an aux-carrying composite
// S = g*h (the DF gC analogue) and both CONTRACT the aux batch mode x_1 at
// their root (aux-aux edge, like the DF index K). The whole-scope executor
// drives ONE x_1 loop for the whole forest and must reproduce the unbatched
// forest-descent result to within FP noise (the batched summation reorders the
// contraction). Uses REAL TA data at a small tractable size: the DryRun backend
// is zero-data (see test_scope_executor.cpp's note -- it cannot witness a
// dropped/double-counted root) and the water-20 forest is intractable at real
// tensor sizes, so this exercises the identical sharing structure small.
TEST_CASE("evaluate_whole_scope matches forest descent over one aux loop",
          "[eval][scope-executor]") {
  using sequant::evaluate;
  using sequant::eval::build_scope_schedule;
  using sequant::eval::build_value_node_map;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_whole_scope;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // aux multi-tiled (naux 12 in tiles of 4 -> 3 tiles) so x_1 batches.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  auto const aux =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };

  // Two summable roots sharing S = g*h (carries aux x_1, free at the child,
  // contracted at each root):
  //   F1 = (g{a2;i1;x1} * h{i3;a2}) * (p{a3;i2;x1} * q{i4;a3})
  //   F2 = (g{a2;i1;x1} * h{i3;a2}) * (r{a3;i2;x1} * w{i4;a3})
  // both contract x_1 (aux-aux) -> results over {i1,i2,i3,i4}.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (r{a_3;i_2;x_1} * w{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  sequant::Index x1;
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    x1 = *ax;
    nd->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
    nd->set_batch_order_aware(true);
  }
  REQUIRE(x1.space() == aux);

  // Reference: unbatched forest descent (ground truth). Without a custom
  // evaluator batched_here is ignored, so this is the exact F1 + F2.
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // Scope schedule: x_1 is the single realized batch loop.
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = compute_dag_boulevard(forest, cm, block_of);
  auto sched = build_scope_schedule<std::wstring>(
      rich, {std::wstring(x1.space().base_key())});
  REQUIRE(sched.root.children.size() == 1);
  REQUIRE(sched.root.children.front().kind ==
          sequant::BatchModeType::Contracted);

  // value_id -> node bridge: the aux scope must home the shared composite S.
  auto const vmap = build_value_node_map(forest);
  std::size_t n_homed = 0;
  for (auto vid : sched.root.children.front().homed_values)
    if (vmap.count(rich.cells[vid].hash)) ++n_homed;
  REQUIRE(n_homed > 0);

  // Whole-scope evaluation over the single aux loop.
  auto ws_cache = sequant::CacheManager<node_t>::empty();
  auto const target_batch = [](sequant::Index const&) -> std::size_t {
    return 4;
  };
  auto const got = evaluate_whole_scope(forest, sched, rich, target, yield_,
                                        ws_cache, target_batch)
                       ->get<TArrayD>();

  // Agreement to within FP noise (batched summation reorders the contraction).
  // Both `got` and `ref` are permuted to the `target` layout, so `target` is
  // their common annotation.
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

// SP3 Task 2 of the ordered-scope batched-eval design (the sequel to the
// whole-scope batched DAG execution design; see ordered_schedule.hpp's own
// doc comments for the SP2 OrderedSchedule IR and ordered_executor.hpp's
// detail::run_ordered_contracted_block for the executor this test drives):
// EQUIVALENCE for a realized Contracted loop block (AccumulateSum reduction).
// Same sharing structure as the "evaluate_whole_scope matches forest descent
// over one aux loop" test above (a shared aux-carrying composite S = g*h,
// both roots contracting the aux batch mode x_1 AT THEIR OWN ROOT, so each
// root itself is the block's AccumulateSum output -- "a reduction value
// consumed by a value built after the loop", per the SP3 brief, realized here
// as the shared per-root combine step that runs after evaluate_ordered_
// schedule's loop-block walk closes), but driven through the SP1/SP2/SP3
// pipeline (analyze_legality -> build_ordered_schedule -> evaluate_ordered_
// schedule) instead of the whole-scope ScopeSchedule/walk_scope path. Real TA
// data at a small tractable size (the DryRun backend is zero-data and cannot
// witness a dropped/mis-accumulated batch).
TEST_CASE(
    "evaluate_ordered_schedule matches forest descent over one Contracted "
    "loop block",
    "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::build_value_node_map;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using sequant::eval::ScopeBlock;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // aux multi-tiled (naux 12 in tiles of 4 -> 3 tiles) so x_1 batches.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  auto const aux =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");

  // Two summable roots sharing S = g*h (carries aux x_1, free at the child,
  // contracted at each root -- an AccumulateSum output of the x_1 block):
  //   F1 = (g{a2;i1;x1} * h{i3;a2}) * (p{a3;i2;x1} * q{i4;a3})
  //   F2 = (g{a2;i1;x1} * h{i3;a2}) * (r{a3;i2;x1} * w{i4;a3})
  // both contract x_1 (aux-aux) -> results over {i1,i2,i3,i4}. Same forest as
  // the whole-scope aux-only equivalence test above.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (r{a_3;i_2;x_1} * w{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  // Reference: unbatched forest descent (ground truth), taken BEFORE
  // batched_here() is stamped below -- forest descent ignores it regardless
  // (no custom evaluator installed), but this keeps the reference call
  // manifestly independent of the annotation.
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // SP1's compute_dag_boulevard (peak_profile.hpp) builds each occurrence's
  // enclosing-loop context (OccurrenceRec::ectx) off batched_here() during its
  // descent -- the SAME annotation the runtime batched evaluator consults, and
  // the one build_ordered_schedule's own fixture (test_ordered_schedule.cpp's
  // W/{Κ} forest) stamps on the node that REALIZES the loop. Without it no
  // occurrence has an enclosing x_1 loop, so classify_axis's LoopLocal check
  // (bound lockstep to an enclosing loop of the SAME index) can never find
  // one and every x_1-carrying value falls to LoopCarried (AccumulateScatter)
  // instead of Reduction/LoopLocal -- exactly the wrong classification this
  // test's own precondition check caught. Each root itself is the node that
  // contracts x_1 (batch_axis finds it as the shared index between the root's
  // two factors), so batched_here() goes there, mirroring the whole-scope
  // aux-only equivalence test above.
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  sequant::Index x1;
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    x1 = *ax;
    nd->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
  }
  REQUIRE(x1.space() == aux);

  // SP1/SP2: x_1 (the aux space) is the only batchable Contracted axis; no
  // External axis at all -- the natural fixture for a Contracted-only Task 2.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"x"});

  // Precondition this test exists to exercise: a realized Contracted loop
  // block at the root, with at least one AccumulateSum output (the shared
  // reduction each root's own build resolves to) -- the shape Task 1's own
  // test explicitly does NOT produce.
  ScopeBlock const* x_block = nullptr;
  for (auto const& step : ordered.root.steps)
    if (auto const* b = std::get_if<ScopeBlock>(&step.value)) x_block = b;
  REQUIRE(x_block != nullptr);
  REQUIRE(x_block->kind == sequant::BatchModeType::Contracted);
  REQUIRE(!x_block->outputs.empty());
  for (auto const& [vid, kind] : x_block->outputs)
    REQUIRE(kind == sequant::eval::OutputKind::AccumulateSum);

  // value_id -> node bridge: the x_1 block must home the shared composite S
  // as a plain BuildStep (LoopLocal), exactly like the whole-scope test's
  // n_homed > 0 check.
  auto const vmap = build_value_node_map(forest);
  std::size_t n_build_ids = 0;
  for (auto const& step : x_block->steps)
    if (std::holds_alternative<sequant::eval::BuildStep>(step.value))
      ++n_build_ids;
  REQUIRE(n_build_ids > 0);

  // evaluate_ordered_schedule over the single x_1 loop block.
  auto cache = sequant::CacheManager<node_t>::empty();
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  // Agreement to within FP noise (the batched accumulation reorders the
  // x_1 contraction vs the unbatched reference).
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

// Task 4 of the whole-scope batched DAG execution design
// (doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md):
// NESTED-scope EQUIVALENCE. A two-root forest with a NESTED scope tree --
// aux x_1 (CONTRACTED, OUTER) over occ i_1 (EXTERNAL/spectator, INNER) --
// exercising the recursive walk (one batch loop per level) with accumulate at
// the outer (contracted) exit and scatter at the inner (external) exit. The
// roots SHARE an aux-only composite S = g*h (carries x_1, invariant to i_1),
// homed at the OUTER x level: the executor must build it once per x-block and
// reuse it across every inner i-block, while still reproducing the unbatched
// forest-descent result to FP noise (the batched summation reorders both the
// x-contraction and the i-scatter). Real TA data at a small tractable size,
// exactly as the Task-3 aux-only equivalence test (the DryRun backend is
// zero-data; a dropped/double-counted slice needs real arithmetic to witness).
TEST_CASE("evaluate_whole_scope matches forest descent over nested aux+occ",
          "[eval][scope-executor]") {
  using sequant::evaluate;
  using sequant::eval::build_scope_schedule;
  using sequant::eval::build_value_node_map;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_whole_scope;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // occ multi-tiled (nocc 8 in tiles of 4 -> 2 tiles) so i_1 batches (external
  // scatter), and aux multi-tiled (naux 12 in tiles of 4 -> 3 tiles) so x_1
  // batches (contracted accumulate).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 8, 6, 12};
  yield_.set_max_tile(4);

  // Two summable roots sharing the aux-only composite S = g*h (carries x_1,
  // invariant to i_1 -> homed at the OUTER x level), each with its own occ+aux
  // subproduct P_k = u_k*w_k (carries both x_1 and i_1,i_2 -> homed at the
  // INNER i level):
  //   F1 = (g{a2;a3;x1} * h{a3;a2}) * (u1{i1;a4;x1} * w1{a4;i2})
  //   F2 = (g{a2;a3;x1} * h{a3;a2}) * (u2{i1;a4;x1} * w2{a4;i2})
  // both contract x_1 (aux-aux) -> results over {i_1, i_2}; i_1 is the batched
  // external (spectator) mode, i_2 an un-batched external.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;a_3;x_1} * h{a_3;a_2}) * (u1{i_1;a_4;x_1} * w1{a_4;i_2})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;a_3;x_1} * h{a_3;a_2}) * (u2{i_1;a_4;x_1} * w2{a_4;i_2})");
  std::string const target = "i_1,i_2";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  sequant::Index const x1{L"x_1"};
  sequant::Index const i1{L"i_1"};
  for (auto& nd : forest) {
    nd->set_batched_here({{x1, sequant::BatchModeType::Contracted},
                          {i1, sequant::BatchModeType::External}});
    nd->set_batch_order_aware(true);
  }

  // Reference: unbatched forest descent (ground truth). Without a custom
  // evaluator batched_here is ignored, so this is the exact F1 + F2.
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // Scope schedule: aux x_1 outer (contracted), occ i_1 inner (external).
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = compute_dag_boulevard(forest, cm, block_of);
  auto sched = build_scope_schedule<std::wstring>(rich, {L"x", L"i"});
  REQUIRE(sched.root.children.size() == 1);
  auto const& xscope = sched.root.children.front();
  REQUIRE(xscope.mode.space().base_key() == L"x");
  REQUIRE(xscope.kind == sequant::BatchModeType::Contracted);
  REQUIRE(xscope.children.size() == 1);
  auto const& iscope = xscope.children.front();
  REQUIRE(iscope.mode.space().base_key() == L"i");
  REQUIRE(iscope.kind == sequant::BatchModeType::External);
  REQUIRE(iscope.children.empty());

  // value_id -> node bridge: the OUTER x scope must home the shared aux-only
  // composite S (carries x_1, not i_1).
  auto const vmap = build_value_node_map(forest);
  std::size_t n_x_homed = 0;
  for (auto vid : xscope.homed_values)
    if (vmap.count(rich.cells[vid].hash)) ++n_x_homed;
  REQUIRE(n_x_homed > 0);

  // Whole-scope evaluation over the nested aux+occ loop nest.
  auto ws_cache = sequant::CacheManager<node_t>::empty();
  auto const target_batch = [](sequant::Index const&) -> std::size_t {
    return 4;
  };
  auto const got = evaluate_whole_scope(forest, sched, rich, target, yield_,
                                        ws_cache, target_batch)
                       ->get<TArrayD>();

  // Agreement to within FP noise (the nested batched loop reorders the
  // x-contraction and assembles i_1 by disjoint-slice scatter).
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

// Task 6 of the whole-scope batched DAG execution design
// (doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md):
// COEXISTENCE. sequant::evaluate(forest, policy, layout, leaf, cache, ...) --
// the driver-entry overload added by scope_executor.hpp -- must route to
// eval::evaluate_whole_scope when policy.whole_scope_execution is true, and
// to the UNMODIFIED sequant::evaluate(Nodes const&, layout, ...) forest
// descent (BYTE-IDENTICAL, no schedule built) when it is false. Same aux-only
// two-root forest as the Task-3 equivalence test above (S = g*h shared,
// contracts x_1 at each root).
//
// Discriminator: without a custom evaluator, plain forest descent ignores
// batched_here entirely, so a NAIVE numeric-equivalence-to-reference check
// alone cannot tell "routed to whole-scope" apart from "silently fell through
// to forest descent" -- both would land close to the unbatched reference. The
// flag-ON branch is instead cross-checked against an INDEPENDENT direct call
// to eval::evaluate_whole_scope (same schedule/target, fresh cache): agreement
// there, given whole-scope's batched summation reorders arithmetic vs the
// unbatched reference (see the Task-3 test above), is real evidence the
// dispatcher actually took the whole-scope path. The flag-OFF branch is
// checked against the reference at a tolerance (1e-13) far tighter than the
// flag-ON branch's (1e-10 -- a genuine batched reordering): it is (by
// construction) the very same function call, so any disagreement above the
// machine-epsilon noise floor of TA's own tile-wise reduction order (not
// guaranteed deterministic across independently scheduled World tasks, even
// for the identical, unreordered contraction) would be a routing bug.
TEST_CASE(
    "sequant::evaluate(forest, policy, ...) routes to whole-scope vs forest "
    "descent by BatchPolicy::whole_scope_execution",
    "[eval][scope-executor]") {
  using sequant::evaluate;
  using sequant::eval::build_scope_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_whole_scope;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  auto const aux =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };

  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (r{a_3;i_2;x_1} * w{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  sequant::Index x1;
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    x1 = *ax;
    nd->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
    nd->set_batch_order_aware(true);
  }
  REQUIRE(x1.space() == aux);

  // Reference: today's unbatched forest descent (fresh cache).
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  auto const target_batch = [](sequant::Index const&) -> std::size_t {
    return 4;
  };

  // Independent direct whole-scope call (fresh cache), for the flag-ON
  // cross-check.
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = compute_dag_boulevard(forest, cm, block_of);
  auto sched = build_scope_schedule<std::wstring>(
      rich, {std::wstring(x1.space().base_key())});
  REQUIRE(sched.root.children.size() == 1);
  auto ws_ref_cache = sequant::CacheManager<node_t>::empty();
  auto const direct_ws =
      evaluate_whole_scope(forest, sched, rich, target, yield_, ws_ref_cache,
                           target_batch)
          ->get<TArrayD>();

  // ---- flag ON: dispatcher must agree with the independent whole-scope
  // call (both realize the identical batched-summation algorithm; agreement
  // only up to FP noise, same convention as "make_evaluator BatchPolicy
  // adapter" above) and stay within the Task-3 tolerance of the unbatched
  // reference.
  {
    sequant::BatchPolicy policy_on;
    policy_on.batch_target_size = target_batch;
    policy_on.whole_scope_execution = true;

    auto cache_on = sequant::CacheManager<node_t>::empty();
    auto const got_on = evaluate(forest, policy_on, target, yield_, cache_on,
                                 {std::wstring(x1.space().base_key())})
                            ->get<TArrayD>();

    CHECK(equal_tarrays<Loose>(got_on, direct_ws));

    TArrayD diff;
    diff(target) = got_on(target) - ref(target);
    double const rel = TA::norm2(diff) / TA::norm2(ref);
    INFO("relative L2 diff vs unbatched reference = " << rel);
    CHECK(rel < 1e-10);
  }

  // ---- flag OFF: dispatcher must equal today's forest descent -- the same
  // underlying sequant::evaluate(Nodes const&, layout, ...) call, no schedule
  // built, no reordering. Compared by RELATIVE L2 norm rather than
  // equal_tarrays<Tight>'s fixed ABSOLUTE margin (~2 ULP): TA's own tile-wise
  // reduction order is not guaranteed deterministic across two independently
  // scheduled World tasks even for the identical, unreordered contraction (a
  // pre-existing TA property, unrelated to this driver), so a fixed absolute
  // margin can spuriously fail on an O(10)-scale result even when the
  // relative agreement is at the machine-epsilon noise floor. The bound here
  // (1e-13) is far tighter than the flag-ON branch's 1e-10 -- which reflects
  // a GENUINE reordering (batched summation) -- so the two remain clearly
  // distinguishable.
  {
    sequant::BatchPolicy const policy_off;  // whole_scope_execution = false
    auto cache_off = sequant::CacheManager<node_t>::empty();
    auto const got_off =
        evaluate(forest, policy_off, target, yield_, cache_off)->get<TArrayD>();
    TArrayD diff;
    diff(target) = got_off(target) - ref(target);
    double const rel_off = TA::norm2(diff) / TA::norm2(ref);
    INFO("relative L2 diff (flag OFF vs reference) = " << rel_off);
    CHECK(rel_off < 1e-13);
  }
}

TEST_CASE("eval_batched_custom_evaluator dedups within-batch repeats",
          "[.][blocked-on-per-context-caching]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12, 12};
  yield_.set_max_tile(4);

  // W-analog: root contracts the aux index x_1; the two children are
  // canonically equal
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (g{a_3;i_2;x_1} * h{i_4;a_3})");
  std::string const target = "i_1,i_3,i_2,i_4";
  auto node = eval_node(expr);
  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node);
  REQUIRE(ax.has_value());
  node->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  std::map<std::wstring, int> n_yield;
  auto counting_yield = [&yield_, &n_yield](node_t const& n) {
    if (n->is_tensor()) ++n_yield[std::wstring(n->as_tensor().label())];
    return yield_(n);
  };

  // x_1 is auxiliary: extent 12, tiles of 4; target 4 elements -> 3 batches
  int const n_b = 3;
  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; }));
  auto const res =
      evaluate(node, target, counting_yield, cache)->get<TArrayD>();
  REQUIRE(equal_tarrays<Loose>(res, ref));

  // with within-batch dedup: per batch the left child evaluates (g and h
  // yielded once each), the right child is a scratch cache hit; plus one g
  // from the evaluator's mode_batches probe of the x_1-carrying leaf
  CHECK(n_yield[L"g"] == n_b + 1);
  CHECK(n_yield[L"h"] == n_b);
}

// BLOCKED, hidden by default -- see .superpowers/sdd/oamb-a0-note.md sections
// 10.4 and 11. This test batched an UNANNOTATED node via the runtime's depth-0
// heuristic fallback, which has been removed (annotations are now
// authoritative). Re-pointing it at an explicit annotation is NOT
// behaviour-preserving: cache_manager vetoes caching for a node whose own
// batched_here() carries a sliced batchable mode, and the heuristic never set
// batched_here() -- so the same mode batched by the two routes gives different
// CSE. The correct expectations depend on Phase B replacing that veto with
// per-context (per-slice) caching. Re-enable and re-derive the counts then;
// do not "fix" the numbers against veto behaviour that is about to be deleted.
TEST_CASE("eval_batched_custom_evaluator group replay",
          "[.][blocked-on-per-context-caching]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12, 12};
  yield_.set_max_tile(4);

  // Two persistent finals sharing the sub-intermediate S = g*h (carries the
  // auxiliary batch mode x_1):
  //   F1 = S * S'   (canonically-equal siblings; contracts x_1, aux-aux)
  //   F2 = S * p    (contracts x_1, aux-aux)
  // Volatile heads (label "t") make F1 and F2 persistent. Every orbital
  // contraction pairs a bra with a ket.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * (g{a_3;i_2;x_1} * h{i_4;a_3}))"
      L" * t{i_1,i_2;i_3,i_9}");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * p{i_5;i_6;x_1}) * t{i_1;i_3}");
  std::string const tgt1 = "i_4,i_9";
  std::string const tgt2 = "i_5,i_6";
  auto n1 = eval_node(t1);
  auto n2 = eval_node(t2);
  // Annotate explicitly: the runtime batches STRICTLY on the optimizer's
  // annotations (no heuristic fallback), so each final must state the mode it
  // wants batched. batch_axis() picks exactly what the old depth-0 heuristic
  // picked, so this reproduces the behaviour these assertions were written for.
  for (auto* nd : {&n1, &n2}) {
    auto const ax = sequant::batch_axis(*nd);
    REQUIRE(ax.has_value());
    (*nd)->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
  }

  auto is_volatile_t = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  // references (unbatched, uncached, uncounted)
  auto const ref1 = evaluate(n1, tgt1, yield_)->get<TArrayD>();
  auto const ref2 = evaluate(n2, tgt2, yield_)->get<TArrayD>();

  // real cache over the two-term forest; F1 and F2 must classify persistent
  auto cache = sequant::cache_manager(std::vector{n1, n2}, is_volatile_t);
  REQUIRE(cache.persistent(n1.left()));
  REQUIRE(cache.persistent(n2.left()));

  std::map<std::wstring, int> n_yield;
  auto counting_yield = [&yield_, &n_yield](node_t const& n) {
    if (n->is_tensor()) ++n_yield[std::wstring(n->as_tensor().label())];
    return yield_(n);
  };
  // persistent_only=true: gate volatile term heads so the group replay streams
  // only the persistent finals F1/F2 (the across-the-board default would also
  // batch volatile triggers, changing the yield counts asserted below).
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      sequant::accept_any_index{}, sequant::make_no_scope_guard{},
      is_volatile_t,
      /*persistent_only=*/true));

  // evaluating term 1 triggers at F1 and must prebuild F2 in the same passes
  auto const res1 = evaluate(n1, tgt1, counting_yield, cache)->get<TArrayD>();
  REQUIRE(cache.alive(n2.left()));  // F2 prebuilt by the group replay

  auto const res2 = evaluate(n2, tgt2, counting_yield, cache)->get<TArrayD>();
  REQUIRE(equal_tarrays<Loose>(res1, ref1));
  REQUIRE(equal_tarrays<Loose>(res2, ref2));

  // S evaluated once per batch (n_b = 3), shared by F1 (twice) and F2 (once):
  // g: 3 (S evals) + 1 (trigger probe) + 1 (F2 candidacy probe) = 5
  // h: 3; p: 3 (sliced, once per batch); t: 2 (one per term head)
  CHECK(n_yield[L"g"] == 5);
  CHECK(n_yield[L"h"] == 3);
  CHECK(n_yield[L"p"] == 3);
  CHECK(n_yield[L"t"] == 2);
}

// BLOCKED, hidden by default -- see .superpowers/sdd/oamb-a0-note.md sections
// 10.4 and 11. This test batched an UNANNOTATED node via the runtime's depth-0
// heuristic fallback, which has been removed (annotations are now
// authoritative). Re-pointing it at an explicit annotation is NOT
// behaviour-preserving: cache_manager vetoes caching for a node whose own
// batched_here() carries a sliced batchable mode, and the heuristic never set
// batched_here() -- so the same mode batched by the two routes gives different
// CSE. The correct expectations depend on Phase B replacing that veto with
// per-context (per-slice) caching. Re-enable and re-derive the counts then;
// do not "fix" the numbers against veto behaviour that is about to be deleted.
TEST_CASE("eval_batched_custom_evaluator group replay layers nested finals",
          "[.][blocked-on-per-context-caching]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12, 12};
  yield_.set_max_tile(4);

  // F_in = (g*h)*p (contracts the aux index x_1; persistent) nests inside
  // F_out, which contracts its own aux mode x_2:  F_out = (F_in * r) * q,
  // with r and q carrying x_2. Triggering at F_out must build F_in in an
  // inner layer first, then seed its full value into F_out's pass (F_in
  // carries no x_2). Every orbital contraction pairs a bra with a ket.
  auto const t_out = sequant::deserialize<sequant::ExprPtr>(
      L"((((g{a_2;i_1;x_1} * h{i_3;a_2}) * p{i_5;i_6;x_1}) * r{i_6;i_7;x_2})"
      L" * q{i_7;i_8;x_2}) * t{i_1;i_3,i_9}");
  auto const t_in = sequant::deserialize<sequant::ExprPtr>(
      L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * p{i_5;i_6;x_1}) * t{i_1;i_3,i_7}");
  std::string const tgt_out = "i_5,i_8,i_9";
  std::string const tgt_in = "i_5,i_6,i_7";
  auto n_out = eval_node(t_out);
  auto n_in = eval_node(t_in);

  auto is_volatile_t = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  auto const ref_out = evaluate(n_out, tgt_out, yield_)->get<TArrayD>();
  auto const ref_in = evaluate(n_in, tgt_in, yield_)->get<TArrayD>();

  auto cache = sequant::cache_manager(std::vector{n_out, n_in}, is_volatile_t);
  // structural preconditions: both finals persistent; F_in nests in F_out
  REQUIRE(cache.persistent(n_out.left()));
  REQUIRE(cache.persistent(n_in.left()));
  {
    sequant::TreeNodeEqualityComparator<node_t> const eq;
    REQUIRE(eq(n_out.left().left().left(), n_in.left()));
  }

  std::map<std::wstring, int> n_yield;
  auto counting_yield = [&yield_, &n_yield](node_t const& n) {
    if (n->is_tensor()) ++n_yield[std::wstring(n->as_tensor().label())];
    return yield_(n);
  };
  // batch only over auxiliary indices (as a DF-batched application would)
  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Annotate explicitly: the runtime batches STRICTLY on the optimizer's
  // annotations (no heuristic fallback), so each final must state the mode it
  // wants batched. batch_axis() picks exactly what the old depth-0 heuristic
  // picked, so this reproduces the behaviour these assertions were written for.
  for (auto* nd : {&n_out, &n_in}) {
    auto const ax = sequant::batch_axis(*nd, accept_aux);
    REQUIRE(ax.has_value());
    (*nd)->set_batched_here({{*ax, sequant::BatchModeType::Contracted}});
  }

  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, sequant::make_no_scope_guard{}, is_volatile_t));

  auto const res_out =
      evaluate(n_out, tgt_out, counting_yield, cache)->get<TArrayD>();
  REQUIRE(cache.alive(n_in.left()));  // inner layer built and stored

  auto const res_in =
      evaluate(n_in, tgt_in, counting_yield, cache)->get<TArrayD>();
  REQUIRE(equal_tarrays<Loose>(res_out, ref_out));
  REQUIRE(equal_tarrays<Loose>(res_in, ref_in));

  // layer 1 (F_in over x_1, 3 batches): g,h,p once per batch; layer 2 (F_out
  // over x_2): F_in SEEDED (g,h,p untouched), r,q sliced once per batch.
  // Probes: +1 r (trigger, x_2-carrying leaf), +1 g (F_in candidacy, x_1
  // leaf). If seeding failed and F_in were re-derived per batch, g/h/p
  // would be 6+.
  CHECK(n_yield[L"g"] == 4);
  CHECK(n_yield[L"h"] == 3);
  CHECK(n_yield[L"p"] == 3);
  CHECK(n_yield[L"r"] == 4);
  CHECK(n_yield[L"q"] == 3);
  CHECK(n_yield[L"t"] == 2);
}

TEST_CASE("eval_batched_custom_evaluator nests inner mode", "[eval]") {
  // Task 4.2 exactness gate: two batchable modes annotated at DIFFERENT nodes
  // of one tree must batch by nesting -- the outer node slices mode A at the
  // top and, WITHIN each A batch, the evaluator re-enters on the per-batch
  // scratch to slice mode B at an inner node (`for A-batch: for B-batch:
  // replay`). The inner slice must compose on top of the outer one; the nested
  // sum equals the unbatched contraction because `sum_K = sum_b sum_{K in b}`
  // applies per mode.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles) so both
  // batch modes slice into 3 batches each and depth-2 nesting genuinely
  // engages.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // Root (outer) contracts the aux mode x_1; the inner product (u*v) contracts
  // a DIFFERENT aux mode x_2. Crucially, the inner leaf u carries BOTH x_1 and
  // x_2, so the inner re-derivation MUST slice x_1 to the outer batch (compose
  // le_g) -- reusing the original leaf evaluator would leave u at the full x_1
  // extent and the x_1 contraction above would then mismatch. Every product is
  // written fully binary so binarize preserves the nesting.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((u{i_1;i_2;x_1,x_2} * v{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);  // mutable: batch modes are annotated below

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Annotate the root with its aux mode (x_1) and the unique descendant that
  // contracts an aux mode (the inner u*v, x_2). Locating the inner node by its
  // contracted aux mode keeps this robust to binarize's operand ordering.
  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_batched_here({{*root_axis, sequant::BatchModeType::Contracted}});

  node_t* inner = nullptr;
  std::optional<sequant::Index> inner_axis;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        inner = &n;
        inner_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(inner_axis.has_value());
  REQUIRE(*inner_axis != *root_axis);  // the two nested modes differ
  (*inner)->set_batched_here(
      {{*inner_axis, sequant::BatchModeType::Contracted}});

  // Reference: plain (unbatched) evaluation. Computed first so yield_'s random
  // leaf arrays are generated and cached; the batched evaluator reuses them.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Spy scope-guard: records the batch count each time an evaluator FIRES (i.e.
  // picks a mode partitioning into >1 batch). Without the re-entrant scratch
  // only the root fires (one record); with nesting the inner evaluator fires
  // once per outer batch -- this vector is what proves depth-2 nesting engaged.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, spy, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness gate: the nested batched result equals the unbatched one.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Depth-2 nesting engaged: the root fires once over x_1 (3 batches) and, per
  // x_1 batch, the inner evaluator fires over x_2 (3 batches) -- 1 + 3 = 4
  // firings, each partitioning into 3 batches. Before the re-entrant scratch,
  // only the root fired (guard_calls == {3}); size > 1 is the RED/GREEN gate.
  REQUIRE(guard_calls.size() > 1);
  CHECK(guard_calls.size() == 4);
  for (auto n : guard_calls) CHECK(n == 3);
}

TEST_CASE("eval_batched_custom_evaluator nests two modes on one node",
          "[eval]") {
  // Finding N2 regression gate: the single-term-opt DP can stamp MORE THAN ONE
  // batch mode on a SINGLE eval node (aprime is a bitmask;
  // reconstruct_batched_modes pushes one Index per set bit into that node's
  // batched_here()). The runtime must slice ALL of them by nesting -- `for
  // K-batch: for mu1-batch: replay` at the SAME node -- otherwise the modeled
  // peak (which priced BOTH modes sliced) is a lie. Here ONE product node
  // carries two aux batch modes x_1 and x_2, both present on both leaves; the
  // evaluator must slice x_1 at depth 0 and, WITHIN each x_1 batch, re-enter on
  // the SAME node to slice x_2 at depth 1. Before the sliceability-aware pick
  // the runtime re-picked x_1 (the first annotated mode) on the re-entry, found
  // it already fully sliced (one batch), and DECLINED -- x_2 was never sliced
  // (max_depth stayed 1).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // ONE product contracting i_1 (occ) plus BOTH aux modes x_1 and x_2; each
  // leaf carries x_1 and x_2, so slicing either aux mode slices both leaves.
  // Result free indices are i_5, i_6.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{i_1;i_5;x_1,x_2} * h{i_6;i_1;x_1,x_2})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);  // mutable: both batch modes annotated below

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Stamp BOTH aux contracted indices on the single node (the two-set-bit
  // aprime case). contracted_indices lists them; keep the aux ones.
  sequant::container::svector<sequant::Index> two_axes;
  for (sequant::Index const& ix : sequant::contracted_indices(node))
    if (accept_aux(ix)) two_axes.push_back(ix);
  REQUIRE(two_axes.size() == 2);
  REQUIRE(two_axes[0] != two_axes[1]);
  sequant::container::svector<std::pair<sequant::Index, sequant::BatchModeType>>
      two_axes_typed;
  for (sequant::Index const& ix : two_axes)
    two_axes_typed.push_back({ix, sequant::BatchModeType::Contracted});
  node->set_batched_here(two_axes_typed);

  // Reference: plain (unbatched) evaluation; also generates yield_'s random
  // leaf arrays that the batched evaluator reuses.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Distinct per-mode target sizes so the two nested levels realize DIFFERENT
  // batch counts. Batches snap to whole tiles of 4 and the target is an upper
  // bound: x_1 -> target 4 over extent 12 = 3 single-tile batches; x_2 ->
  // target 8 = 2 batches (two 4-tiles fit, then the last). A depth-2 product of
  // 3*2 = 6 (rather than 3*3 or 2*2) proves the outer and inner levels sliced
  // DISTINCT modes with distinct partitions.
  auto const& x1 = two_axes[0];
  auto target_batch_size = [x1](sequant::Index const& ix) -> std::size_t {
    return ix == x1 ? std::size_t{4} : std::size_t{8};
  };

  // Tracking guard (as in the compose test): records the live-guard stack so
  // depth-2 simultaneity and the per-instant product are observable.
  struct GuardState {
    std::vector<std::size_t> live;
    std::size_t max_depth = 0;
    std::vector<std::size_t> counts;
    std::vector<std::size_t> products_at_depth2;
  } state;
  struct TrackingGuard {
    GuardState* st;
    TrackingGuard(GuardState* s, std::size_t n) : st(s) {
      st->counts.push_back(n);
      st->live.push_back(n);
      st->max_depth = std::max(st->max_depth, st->live.size());
      if (st->live.size() == 2)
        st->products_at_depth2.push_back(st->live[0] * st->live[1]);
    }
    TrackingGuard(TrackingGuard const&) = delete;
    TrackingGuard& operator=(TrackingGuard const&) = delete;
    ~TrackingGuard() { st->live.pop_back(); }
  };
  auto make_tracking_guard = [&state](std::size_t n) {
    return TrackingGuard(&state, n);
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, target_batch_size, accept_aux, make_tracking_guard,
      sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: nesting two modes on one node changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Depth 2 reached on ONE node: x_1 sliced at depth 0, x_2 at depth 1. Before
  // the fix the node re-picked x_1 and declined -> max_depth would stay 1.
  REQUIRE(state.max_depth == 2);
  // The outer level fires once over x_1 (3 batches); per outer batch the inner
  // level fires once over x_2 (2 batches) -> 1 + 3 = 4 firings, counts a
  // permutation of {3, 2, 2, 2}.
  REQUIRE(state.counts.size() == 4);
  CHECK(std::count(state.counts.begin(), state.counts.end(), 3u) == 1);
  CHECK(std::count(state.counts.begin(), state.counts.end(), 2u) == 3);
  // Every depth-2 instant multiplies 3 (outer x_1) by 2 (inner x_2) = 6, i.e.
  // two DISTINCT modes with DISTINCT partitions were live together.
  REQUIRE(state.products_at_depth2.size() == 3);
  for (auto const p : state.products_at_depth2) CHECK(p == 6);
  CHECK(state.live.empty());
}

TEST_CASE("eval_batched_custom_evaluator hoists loop-invariant descendant",
          "[eval]") {
  // Task 4 (B2) hoisting gate: a network batched over an OUTER mode (x_1) whose
  // subtree contains an intermediate INVARIANT to that mode (I2 = a*b contracts
  // x_2 and carries no x_1) must build the invariant ONCE at its emitted scope
  // level and reuse it across the x_1 batches -- NOT rebuild it per batch. The
  // root R contracts x_1 (the outer batch loop, depth 0); I2 is a descendant
  // with scope_level == -1 (invariant to the whole nest), so it hoists to the
  // real cache and the scope chain serves it to every x_1 batch body.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles), so
  // x_1 slices into 3 batches.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // R = (((g{i_1;i_2;x_2} * h{i_1;i_3;x_2}) * w{i_2;i_5;x_1}) *
  // p{i_3;i_6;x_1}), a left-deep chain (same shape as the two-mode-nesting
  // test) whose INNERMOST pair is the x_1-invariant intermediate:
  // - I2 = g*h contracts i_1 and x_2 -> {i_2;i_3}: carries NO aux, so it is
  //   invariant to the outer x_1 loop.
  // - M = I2*w contracts i_2 -> {i_3;i_5;x_1}: carries x_1 (not hoistable).
  // - R = M*p contracts i_3 and x_1 -> {i_5;i_6}: the x_1 batch trigger.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((g{i_1;i_2;x_2} * h{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Root batches over its contracted aux mode x_1 (outer loop). Its own
  // scope_level is irrelevant (the trigger is never hoisted); leave it unset.
  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_batched_here({{*root_axis, sequant::BatchModeType::Contracted}});

  // Locate I2 = the unique NON-root node that contracts an aux mode (x_2).
  node_t* i2 = nullptr;
  std::optional<sequant::Index> i2_axis;
  auto find_i2 = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        i2 = &n;
        i2_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_i2(find_i2, node);
  REQUIRE(i2 != nullptr);
  REQUIRE(i2_axis.has_value());
  REQUIRE(*i2_axis != *root_axis);  // I2 contracts a DIFFERENT aux mode
  // Annotate I2: batched over x_2, order-aware (emitted), and EMPTY residency
  // (empty sliced_modes) -- it is invariant to the whole enclosing nest, so
  // per-level placement hoists it above every batch loop to the real/term
  // (root) cache and builds it once. This drives the post-hoist_invariants
  // residency signals (batch_order_aware + sliced_modes), not a per-node scalar
  // placement level.
  (*i2)->set_batched_here({{*i2_axis, sequant::BatchModeType::Contracted}});
  (*i2)->set_batch_order_aware(true);

  // Reference: plain (unbatched) evaluation; also fills yield_'s leaf cache.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Counting leaf evaluator: counts how many times leaf `g` (which appears ONLY
  // inside I2) is requested. I2 built once => g requested once; I2 rebuilt per
  // x_1 (and per inner x_2) batch => requested many times. This is the
  // build-count proxy for I2.
  int g_evals = 0;
  auto counting_yield = [&yield_,
                         &g_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && leaf->as_tensor().label() == L"g") ++g_evals;
    return yield_(leaf);
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
  auto const res =
      evaluate(node, target, counting_yield, cache)->get<TArrayD>();

  // Exactness gate: hoisting changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Hoisting gate (RED/GREEN): the x_1-invariant intermediate I2 is built ONCE,
  // not once per x_1 batch. Before the hoist path, I2 was rebuilt (and itself
  // x_2-batched) inside every x_1 batch, so leaf `g` was evaluated many times.
  CHECK(g_evals == 1);
}

TEST_CASE(
    "eval_batched_custom_evaluator hoists to an intermediate contracted-mode "
    "level",
    "[eval]") {
  // Task 3 review Finding 2 witness: the runtime residency is sliced_modes(),
  // the all-batched-modes cross-occurrence meet -- which carries a node's
  // enclosing CONTRACTED (aux) modes directly, since a node variant to an outer
  // aux loop carries that aux free on a result slot. The pre-existing hoist
  // test only exercises the EMPTY-residency case (whole-nest invariant ->
  // root). This test exercises the NON-EMPTY case: a node P carries an OUTER
  // contracted mode (x_1) but not the INNER one (x_2) its containing trigger
  // batches over, so sliced_modes(P) == {x_1} and its home level is the OUTER
  // (x_1) scratch, not the term/root cache and not T's own inner (x_2) scratch.
  // P must therefore be built ONCE PER OUTER (x_1) BLOCK and reused across
  // every inner (x_2) batch within that block: rebuilt when x_1 advances,
  // untouched while only x_2 advances.
  //
  // Network (left-deep chain; verified by direct inspection of the binarized
  // tree -- canon_indices()/batch_axis() are LOCAL to a node's own top-level
  // pairing, not recomputed transitively, so the grouping below is exactly
  // what deserialize+binarize produce for this expression):
  //   D2 = q{i_1;a_9;x_1} * r{a_9;i_2}         -> {i_1;i_2;x_1}: carries the
  //        OUTER mode x_1, does NOT carry the INNER mode x_2 -- the witness
  //        node (non-empty sliced_modes = {x_1}).
  //   Q  = D2 * s{;;x_1,x_2}                    -> {i_1;i_2;x_2}: a drop-in
  //        structural replacement for a plain 3-index leaf (same free-index
  //        role as `g` in the "hoists loop-invariant descendant" test), so Q
  //        contracts x_1 internally and carries x_2 onward -- the rest of
  //        the network is unaffected by what is inside Q.
  //   T  = Q * h{i_1;i_3;x_2}                    -> {i_2;i_3;x_1}: contracts
  //        x_2 (T is the INNER trigger; D2 is T's grandchild).
  //   M  = T * w{i_2;i_5;x_1}                    -> {i_3;i_5;x_1}
  //   R  = M * p{i_3;i_6;x_1}                    -> {i_5;i_6}: contracts x_1
  //        (R is the OUTER trigger).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ/virt single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 blocks)
  // so both the outer (x_1) and inner (x_2) loops genuinely partition into
  // more than one block each.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((((q{i_1;a_9;x_1} * r{a_9;i_2}) * s{;;x_1,x_2}) * h{i_3;i_1;x_2}) * "
      L"w{i_2;i_5;x_1}) * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Root's own aux mode (x_1, the OUTER trigger).
  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_batched_here({{*root_axis, sequant::BatchModeType::Contracted}});

  // T = the unique non-root node contracting an aux mode different from
  // root's (x_2, the INNER trigger). Located via batch_axis (robust to
  // binarize's operand ordering), matching the existing hoist test's pattern.
  node_t* t = nullptr;
  std::optional<sequant::Index> t_axis;
  auto find_t = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        t = &n;
        t_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_t(find_t, node);
  REQUIRE(t != nullptr);
  REQUIRE(t_axis.has_value());
  REQUIRE(*t_axis != *root_axis);
  (*t)->set_batched_here({{*t_axis, sequant::BatchModeType::Contracted}});

  // D2 = T's grandchild (T.left().left()), the q*r product. Reached by
  // direct structural navigation (confirmed empirically for this exact
  // expression, since batch_axis cannot locate D2 -- it contracts no aux
  // mode itself); guarded by a structural REQUIRE so a future
  // canonicalization change fails loudly here instead of silently testing
  // the wrong node.
  node_t& d2 = node.left().left().left().left();
  REQUIRE(sequant::index_position(d2, *root_axis).has_value());  // carries x_1
  REQUIRE_FALSE(sequant::index_position(d2, *t_axis).has_value());  // not x_2
  REQUIRE(d2.left().leaf());
  REQUIRE(d2.right().leaf());
  REQUIRE(d2.left()->as_tensor().label() == L"q");
  REQUIRE(d2.right()->as_tensor().label() == L"r");

  // Annotate D2: order-aware, with NON-EMPTY sliced_modes = {x_1} (the signal
  // under test -- the all-batched-modes meet carries the enclosing aux mode
  // x_1 that D2 holds free on its result) and no External batched_here stamp
  // (so it is not demoted). Its residency is {x_1} only -- present at T's
  // firing (ectx = [x_1]) but NOT at root's own (ectx = []) -- so it is
  // correctly skipped at the outer pre-loop placement call and collected at
  // T's own (inner) placement call, homed at the OUTER (x_1) scratch level,
  // not the term/root cache and not T's own (x_2) scratch.
  d2->set_batch_order_aware(true);
  d2->set_sliced_modes({*root_axis});

  // Reference: plain (unbatched) evaluation; also fills yield_'s leaf cache.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Counting leaf evaluator: `r` carries no aux index at all, so it is never
  // touched by the sliceability PROBE for either x_1 or x_2 (which only ever
  // requests a leaf that actually carries the probed mode) -- it is touched
  // ONLY when D2 is materialized. This isolates D2's build count, mirroring
  // the h_evals technique used for the external-carrier witness below.
  int r_evals = 0;
  auto counting_yield = [&yield_,
                         &r_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && leaf->as_tensor().label() == L"r") ++r_evals;
    return yield_(leaf);
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
  auto const res =
      evaluate(node, target, counting_yield, cache)->get<TArrayD>();

  // Exactness gate: placing D2 at the x_1 scratch (not term-global, not
  // rebuilt per x_2 sub-batch) changes no numerics -- it is genuinely
  // resliced to the correct x_1 block each time x_1 advances.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // Placement gate: x_1 and x_2 both partition the 12-element aux extent into
  // 3 blocks of 4 with this target size, so a CORRECTLY-homed D2 is rebuilt
  // once per outer block: r_evals == 3. (Empirically checked, not shipped
  // here, since the fetch-time slice-on-use safety net keeps VALUES correct
  // either way -- this is a footprint/avoidable-work gate, not an exactness
  // one, matching the C60 story where a misplaced node is a peak/time bug,
  // not a wrong-answer bug: with d2->set_sliced_modes({*root_axis})
  // commented out, D2's residency is wrongly empty, it is misclassified as a
  // whole-nest invariant and hoisted to the term/root cache instead --
  // r_evals == 1, i.e. held at the FULL x_1 extent for the entire run
  // instead of resliced per block; with batch_order_aware(true) also removed
  // (no hoisting at all -- the pre-existing behavior before this residency
  // component existed), D2 is rebuilt inline once per (x_1, x_2) pair --
  // r_evals == 9.)
  CHECK(r_evals == 3);
}

TEST_CASE("batched_eval_external_axis_scatter", "[eval][batched-external]") {
  // Task 5 exactness gate: a batch mode stamped EXTERNAL on a node is a
  // external index that survives FREE onto the node's result (a Hadamard /
  // batched-matmul mode present on every operand and the result, contracted at
  // no node). The runtime must SCATTER the per-block partials into DISJOINT
  // slices of a pre-sized result (write_into_slice), NOT accumulate them
  // (add_inplace, correct only for a contracted mode). Because the blocks tile
  // the mode without overlap, the scattered result equals the unbatched one
  // EXACTLY -- a memory schedule, never an approximation.
  //
  // The external here is an AUXILIARY index x_1 rather than an occupied one
  // only because canonicalize forbids a NON-auxiliary index shared among > 2
  // tensor slots (no well-defined bra/ket slot type); a high-order aux
  // hyperindex carried into the result IS supported. The runtime scatter
  // mechanism under test is mode-agnostic -- occ vs aux matters only at the
  // DP/cost level (covered by the [dryrun-occ-*] gates). Same rationale as
  // "eval_forest_over_external_occ nests intra-term aux batching".
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles)
  // so the external mode genuinely partitions into > 1 block.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/4,
                                                    /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  // A flat, fully binary tree whose result carries aux x_1 as a pure external:
  //   inner = g*h  -> {a_1, a_3, x_1}   (contracts a_2)
  //   root  = inner*p -> {a_1, a_4, x_1} (contracts a_3)
  // x_1 is present on g, h, p AND the result, contracted at no node; a_1, a_4
  // are ordinary external modes. The ResultExpr LHS pins x_1 EXTERNAL (the
  // ExprPtr binarize would instead infer it contracted); stamping x_1 External
  // on the root and the inner node exercises the scatter branch at both levels.
  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4;x_1})");
  auto node = eval_node(res_expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The external mode x_1, taken from the node's free/result outer modes so
  // the Index object matches the parsed one exactly (a plain outer mode, no
  // protos).
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on result

  // Stamp External on the root AND the inner node (the mode appears free at
  // both levels).
  node->set_batched_here({{mode, sequant::BatchModeType::External}});
  node_t* inner = nullptr;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && sequant::index_position(n, mode).has_value()) inner = &n;
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(sequant::index_position(*inner, mode).has_value());
  (*inner)->set_batched_here({{mode, sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation (batched_here are ignored without a
  // custom evaluator). Computed first so yield_'s random leaves are generated
  // and cached, then reused by the batched run below.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Spy scope-guard: records the block count each time the evaluator FIRES
  // (picks a mode partitioning into > 1 block). A silently-ignored External
  // stamp would never fire -- this is the RED/GREEN witness that the ON path
  // actually blocked.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  // Two block sizes, both strictly below x_1's extent (12): 4 -> 3 blocks,
  // 8 -> 2 blocks. Exactness must hold for every partition.
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    guard_calls.clear();
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_aux, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

    // Exactness: scattering disjoint blocks reconstructs the whole result.
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // The ON path blocked: the root fired at least once, partitioning into
    // > 1 block (block size strictly less than x_1's extent).
    REQUIRE_FALSE(guard_calls.empty());
    for (auto const n : guard_calls) CHECK(n > 1);
    std::size_t const expected_root_blocks = target_batch_size == 4 ? 3 : 2;
    CHECK(guard_calls.front() == expected_root_blocks);
  }
}

TEST_CASE(
    "batched external-carrying intermediate hoists to its seed home; the "
    "result stays exact via slice-on-use (Phase 4b-3: no demotion veto)",
    "[eval][batched-external]") {
  // Phase 4b-3 T1 re-baseline. This test used to pin the `has_demoted_external`
  // demotion veto: in an external SCATTER over a mode e, a DESCENDANT
  // intermediate M that itself CARRIES e free onto its result (an External
  // `batched_here()` stamp absent from `sliced_modes()` -- a meet-demoted
  // carrier, the cross-pair giants in the C60 story) was NEVER hoisted, so it
  // was rebuilt SLICED once per e-block (h requested nblocks times). The veto
  // is now DELETED: placement is purely router-override-or-seed, so with an
  // EMPTY router (no MPQC pre-pass here) M is HOISTED to its seed home. Its
  // residency (`sliced_modes`) is empty, so the seed home is the chain root:
  // M is built exactly ONCE (h requested once), not nblocks times.
  //
  // Crucially this changes PLACEMENT only, never the RESULT: M is cached FULL
  // at the root, and when R consumes it inside the e-scatter the batched
  // Enter-stage slice-on-use slices the ancestor-scope value to the current
  // e-block (the `[eval][slice-on-use]` witness below), so the scatter still
  // reassembles the whole result exactly (assertion 1, unchanged). The move
  // that trimmed a demoted giant's peak is now the remat router's job (T2);
  // here, with no router, the seed placement is what the runtime uses.
  //
  // The sibling INV (order-aware, NO External stamp) was already hoisted under
  // the veto and still is: both siblings now hoist, so h_evals and v_evals are
  // both 1. The test remains a real witness -- it proves the demoted carrier
  // now hoists (was nblocks, now 1) AND that the numerical result is untouched.
  //
  // Network:
  //   R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1})
  //                      * (u{a_3;a_5} * v{a_5;a_4}))
  //  - M = g*h contracts a_2 -> {a_1;a_3;x_1}: CARRIES the external x_1.
  //  - INV = u*v contracts a_5 -> {a_3;a_4}: a plain sibling of M, carrying
  //    NO aux mode at all (a drop-in structural replacement for the leaf `p`
  //    the earlier version of this test used, so R's own shape is unchanged).
  //  - R = M*INV contracts a_3 -> {a_1;a_4;x_1}: the x_1 scatter trigger.
  // M is stamped External on x_1 + order_aware(true), with sliced_modes left
  // empty (default). INV is stamped order_aware(true) only -- no batched_here
  // entries at all. Both are collected + hoisted to the term/root cache (built
  // once, regardless of e-blocks).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles)
  // so the external mode genuinely partitions into > 1 block.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/4,
                                                    /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1})"
      L" * (u{a_3;a_5} * v{a_5;a_4}))");
  auto node = eval_node(res_expr);
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The external mode x_1 (a plain outer aux mode, no protos), taken from the
  // node's free/result modes so the Index matches the parsed one exactly.
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on R

  // Stamp External on the root R (the scatter trigger over x_1).
  node->set_batched_here({{mode, sequant::BatchModeType::External}});

  // Locate M (carries x_1 free) and INV (the sibling that does not) -- the
  // only two non-root, non-leaf nodes in this network.
  node_t* m = nullptr;
  node_t* inv = nullptr;
  auto find_mn = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (sequant::index_position(n, mode).has_value())
        m = &n;
      else
        inv = &n;
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_mn(find_mn, node);
  REQUIRE(m != nullptr);
  REQUIRE(inv != nullptr);
  REQUIRE(m != inv);
  REQUIRE(sequant::index_position(*m, mode).has_value());  // M carries x_1
  REQUIRE_FALSE(
      sequant::index_position(*inv, mode).has_value());  // INV does not

  // M: External stamp, order-aware, EMPTY sliced_modes (default) -- a
  // meet-demoted external carrier. With the veto deleted, M is collected +
  // hoisted to its seed home (the chain root, since its residency is empty),
  // built exactly once; slice-on-use slices the hoisted-full M per e-block.
  (*m)->set_batched_here({{mode, sequant::BatchModeType::External}});
  (*m)->set_batch_order_aware(true);

  // INV: order-aware, no batched_here entries at all (no External stamp) ->
  // collected + hoisted to the term/root cache, built exactly once (as before
  // -- INV was never demoted).
  (*inv)->set_batch_order_aware(true);

  // Reference: plain unbatched evaluation (fills yield_'s leaf cache too).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Counting leaf evaluator. We count leaf h -- the RIGHT operand of M = g*h.
  // find_leaf_carrying(., x_1) returns the LEFTMOST carrier (g), so the
  // sliceability PROBE (mode_batches on the carrier leaf) only ever requests g;
  // h is requested purely when M is actually MATERIALIZED (the g*h
  // contraction). Counting h therefore isolates M's build count from probe
  // noise: hoisted-full => M built once at full e (h requested once, the RED
  // bug); rebuilt sliced per e-block => M built once per block (h requested
  // nblocks times, each at block-e extent, the fix). Leaf v (INV's right
  // operand) is never touched by any sliceability probe (INV carries no aux
  // mode at all), so it isolates INV's build count the same way.
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    std::size_t const nblocks = target_batch_size == 4 ? 3 : 2;
    int h_evals = 0;
    int v_evals = 0;
    auto counting_yield = [&yield_, &h_evals,
                           &v_evals](node_t const& leaf) -> sequant::ResultPtr {
      if (leaf->is_tensor()) {
        if (leaf->as_tensor().label() == L"h") ++h_evals;
        if (leaf->as_tensor().label() == L"v") ++v_evals;
      }
      return yield_(leaf);
    };

    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        counting_yield,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
    auto const res =
        evaluate(node, target, counting_yield, cache)->get<TArrayD>();

    // 1. Slice-exactness: scattering disjoint blocks reconstructs the whole
    //    result exactly (a memory schedule, never an approximation).
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // 2. M (the meet-demoted external carrier) is now HOISTED to its seed
    //    home and built exactly ONCE (h requested once), regardless of
    //    nblocks -- the veto that used to force nblocks sliced rebuilds is
    //    deleted. The result is still exact (assertion 1) because slice-on-use
    //    slices the hoisted-full M to each e-block at the consumer.
    CHECK(h_evals == 1);

    // 3. INV (no External stamp) is hoisted once to the term/root cache too,
    //    regardless of nblocks -- unchanged by the veto deletion.
    CHECK(v_evals == 1);
  }
}

TEST_CASE("batched cached intermediate is sliced to the batch block on use",
          "[eval][slice-on-use]") {
  // Task 3 (slice-on-use) RED/GREEN witness for the LOAD-BEARING fix. A shared
  // intermediate M that carries an external mode e can be cached FULL at an
  // OUTER scope (correct where it is used unbatched / peak-acceptable there)
  // and then consumed inside a NESTED e-external loop -- the C60 summand-46
  // `s.C` pattern (design spec 2026-07-27 Sec 1.1: `I(mu~; a<i1,i2>)` stored
  // full at scope -1, served WHOLE to the external use, contracting the
  // full 8.7 GB into a 2930 GB giant). The fix: the Enter-stage slice-on-use
  // slices a value fetched from an ANCESTOR scope (access_at hops > 0) to the
  // current block for the loops the fetch crossed, so M is realized at BLOCK-e
  // extent at the consumer, never full.
  //
  // Network: R{a_1;a_4;x_1} = (g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4}.
  //  - M = g*h contracts a_2 -> {a_1;a_3;x_1}: carries the external x_1.
  //  - R = M*p contracts a_3 -> {a_1;a_4;x_1}: consumes M under an x_1 block.
  // M is pre-stored FULL in an OUTER cache; an INNER child cache carries a
  // batch_context of one x_1 block. Evaluating R on the inner cache fetches M
  // from the outer scope (hops == 1) and MUST slice it to the block.
  using sequant::evaluate;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using hasher_t = sequant::TreeNodeHasher<node_t>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/4,
                                                    /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4})");
  auto node = eval_node(res_expr);  // R
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  // The external mode x_1 (a plain outer aux mode, no protos).
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (ix.space() == aux_space && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on R

  // Locate M: the unique non-root internal node that carries x_1 free on its
  // result (g*h).
  node_t* m = nullptr;
  auto find_m = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && sequant::index_position(n, mode).has_value()) m = &n;
    self(self, n.left());
    self(self, n.right());
  };
  find_m(find_m, node);
  REQUIRE(m != nullptr);
  REQUIRE(sequant::index_position(*m, mode).has_value());  // M carries x_1
  // Keep the canonical phase trivial so the pre-stored value convention is
  // unambiguous (store V, Enter-stage hit applies phase == 1).
  REQUIRE((*m)->canon_phase() == 1);

  // The x_1 block under test: the first tile, [0, 4) of the extent-12 mode
  // (tile-aligned, as slice_mode requires).
  std::size_t const blk_lo = 0, blk_hi = 4;

  // Reference (full x_1): plain unbatched evaluation. Also fills yield_'s leaf
  // cache so the later runs reuse the same random leaves.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);
  auto const ref_vol = ref.trange().elements_range().volume();

  // Independent BLOCK reference: R computed with every x_1-carrying leaf sliced
  // to [0, 4) (the old le_g-style leaf slicing), i.e. the correct block-x_1
  // result -- built WITHOUT the slice-on-use path.
  auto block_ref_leaf = [&](node_t const& ln) -> sequant::ResultPtr {
    sequant::ResultPtr r = yield_(ln);
    if (auto const p = sequant::index_position(ln, mode))
      return r->slice_mode(*p, blk_lo, blk_hi);
    return r;
  };
  auto const block_ref = evaluate(node, target, block_ref_leaf)->get<TArrayD>();
  auto const block_vol = block_ref.trange().elements_range().volume();
  REQUIRE(block_vol > 0);
  REQUIRE(block_vol * 3 == ref_vol);  // one of three x_1 tiles

  // --- Slice-on-use run. ---
  // OUTER cache holds M FULL; a counting leaf evaluator proves M's children
  // (g, h) are NOT re-evaluated during the run (M is served from cache, sliced
  // on use -- not rebuilt).
  int gh_evals = 0;
  auto counting_yield = [&yield_,
                         &gh_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && (leaf->as_tensor().label() == L"g" ||
                              leaf->as_tensor().label() == L"h"))
      ++gh_evals;
    return yield_(leaf);
  };

  auto const m_full = evaluate(*m, (*m)->annot(), counting_yield);
  REQUIRE(m_full->get<TArrayD>().trange().elements_range().volume() ==
          block_ref.trange().elements_range().volume() * 3);  // full x_1

  std::unordered_map<node_t, std::size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(*m, 100);  // generous life; M stays alive across reads
  auto outer = cache_t(std::move(outer_counts));
  (void)outer.store(*m, m_full);
  REQUIRE(outer.alive(*m));

  auto inner = cache_t::empty();
  inner.set_parent(&outer);
  inner.set_batch_context({{mode, {blk_lo, blk_hi}}});

  gh_evals = 0;  // count only the slice-on-use run
  auto const res =
      evaluate(node, target, counting_yield, inner)->get<TArrayD>();

  // (a) SIZE PROBE: M is realized at BLOCK-x_1 extent at the consumer, so R is
  //     block-sized -- NOT the full-x_1 giant. (RED: served whole -> full-x_1
  //     product -> res_vol == ref_vol.)
  auto const res_vol = res.trange().elements_range().volume();
  REQUIRE(res_vol == block_vol);
  REQUIRE(res_vol * 3 == ref_vol);

  // M's children were never re-evaluated: M was fetched from the outer cache
  // and sliced on use, not rebuilt from g, h.
  REQUIRE(gh_evals == 0);

  // (b) EXACTNESS: the sliced-on-use result equals the independent block
  //     reference exactly (a memory schedule, never an approximation).
  TArrayD diff;
  diff("i,j,k") = block_ref("i,j,k") - res("i,j,k");
  REQUIRE(TA::norm2(diff) < 1e-12);
}

TEST_CASE("batched_scratch_no_seed_external", "[eval][batched-external]") {
  // Task 6 (Part A) RED/GREEN witness for the SEED DECISION in
  // make_batched_scratch. When an EXTERNAL mode is batched, a
  // persistent, alive intermediate that is INVARIANT under the batch's
  // *contracted* member mode but CARRIES that batched External mode must NOT be
  // seeded from the real cache: its full, unsliced-external value is wrong
  // under the external slice, so it must be recomputed sliced. Before this fix
  // the per-node signature and the `!e.sig && persistent && alive` seed rule
  // knew only the contracted member mode, so such a node (contracted signature
  // == absent) was wrongly seeded FULL -- the "seed full i" defect and the root
  // of the original w8-occbatch abort. After the fix the batch's External modes
  // enter the signature: a node carrying one is non-invariant and excluded.
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using sequant::evaluate;
  using sequant::Index;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/4,
                                                    /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  // P = g*h carries the aux external x_1 (the External mode); the volatile
  // head t makes P a PERSISTENT final. P does not carry the contracted batch
  // mode played by `mu` below.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * t{i_1;i_3}");
  auto n = eval_node(expr);
  REQUIRE_FALSE(n.leaf());
  auto const& P = n.left();  // g*h, carries x_1
  REQUIRE_FALSE(P.leaf());

  // the External mode carried by P (a plain aux outer mode)
  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  Index x_ext;
  for (auto const& ix : P->canon_indices())
    if (ix.space() == aux_space && !ix.has_proto_indices()) {
      x_ext = ix;
      break;
    }
  REQUIRE(x_ext.nonnull());
  REQUIRE(sequant::index_position(P, x_ext).has_value());  // P carries it

  // real cache: P classifies persistent (volatile head t); store a value so it
  // is ALIVE (the state that makes a node a seed candidate). Built with n's
  // batched_here() still empty (no External stamp anywhere yet), so the
  // gated cache_manager()'s own internal stamp_lifetime_masks() call (Task 2:
  // it stamps every caller's forest, not just build_dryrun_cache's) sees no
  // External mode and leaves every mask all-full -- the mask veto stays inert
  // and P is persistent here exactly as before that change.
  auto is_volatile_t = [](node_t const& k) {
    return k.leaf() && k->is_tensor() && k->as_tensor().label() == L"t";
  };
  auto real = sequant::cache_manager(std::vector{n}, is_volatile_t);
  REQUIRE(real.persistent(P));
  real.store(P, evaluate(P, P->annot(), yield_));
  REQUIRE(real.alive(P));

  // NOW stamp x_ext External on the member root, as the optimizer would; this
  // is the batch's External-mode set that make_batched_scratch partitions
  // out below. This is a DIFFERENT mechanism from the cache's lifetime mask
  // (make_batched_scratch reads batched_here() directly, never the mask), and
  // is applied only after `real` is built so it does not also feed the
  // cache's own veto -- keeping this test isolated to the seed-decision guard
  // under test.
  n->set_batched_here({{x_ext, sequant::BatchModeType::External}});

  // member mode = a contracted mode P does NOT carry (mu-analog): its per-node
  // signature over P is 'absent', so on the contracted mode alone P looks
  // seedable. The External mode is what must veto the seed.
  Index const mu(L"i_9");  // absent from the tree -> null contracted signature
  REQUIRE_FALSE(sequant::index_position(P, mu).has_value());
  std::vector<std::pair<node_t const*, Index>> const members{{&n, mu}};

  sequant::TreeNodeEqualityComparator<node_t> const eq;
  auto seeds_P = [&](auto const& bs) {
    return std::any_of(bs.seeds.begin(), bs.seeds.end(),
                       [&](node_t const* s) { return eq(*s, P); });
  };

  // POST-FIX: P carries the batched External mode, so it is NOT seeded (it will
  // be recomputed under the external slice).
  auto const bs = sequant::detail::make_batched_scratch(members, real);
  REQUIRE_FALSE(seeds_P(bs));

  // Control: with NO External mode stamped (a purely contracted batch), the
  // SAME invariant persistent P IS seeded -- proving the exclusion is driven by
  // the External stamp and that the contracted-only path is byte-identical.
  n->set_batched_here({});
  auto const bs2 = sequant::detail::make_batched_scratch(members, real);
  REQUIRE(seeds_P(bs2));
}

TEST_CASE("batched_scratch_tot_presize_scatter", "[eval][batched-external]") {
  // Task 6 (Part B): the ToT ResultTensorOfTensorTA::pre_sized_zeros_over_mode
  // must produce a destination that the ToT scatter primitives
  // (write_into_slice -> write_array_into_mode) reassemble EXACTLY. This is the
  // ToT analog of the flat pre-size Task 5 added; CSV/PNO-CCSD residuals carry
  // ToT tiles, so the external-mode scatter needs a ToT pre-size. Here we drive
  // the exact runtime sequence: pre-size from the FIRST block partial (widening
  // its OUTER mode to the carrier's FULL tiling), then write_into_slice
  // every disjoint block. The reassembled ToT must equal the original.
  using sequant::eval_result;
  using sequant::ResultTensorOfTensorTA;
  using sequant::slice_array_over_mode;
  using ToTArray = TA::DistArray<TA::Tensor<TA::Tensor<double>>>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;
  auto& world = TA::get_default_world();

  // outer mode 0 (the external mode): two size-3 tiles (tile-spanning
  // split); outer mode 1: one size-2 tile. Inner tensors are 2x2 with
  // position-dependent values so a mis-scattered block changes the norm.
  TA::TiledRange const outer_tr{{0, 3, 6}, {0, 2}};
  TA::Range const inner_r(std::array<std::size_t, 2>{2, 2});

  auto build = [&world, inner_r](TA::TiledRange const& otr,
                                 bool zero) -> ToTArray {
    auto tile_fn = [inner_r, zero](TA::Range const& orng) {
      TA::Tensor<TA::Tensor<double>> t{orng};
      std::size_t o = 0;
      for (auto const& coord : orng) {
        auto& inner = t[o++];
        inner = TA::Tensor<double>{inner_r};
        double base = 0.0;
        if (!zero)
          for (auto c : coord) base = base * 37.0 + static_cast<double>(c + 1);
        std::size_t k = 0;
        for (auto& x : inner) x = zero ? 0.0 : base * 100.0 + (++k);
      }
      return t;
    };
    ToTArray arr{world, otr};
    for (auto it = arr.begin(); it != arr.end(); ++it)
      if (arr.is_local(it.index()))
        *it = world.taskq.add(tile_fn, it.make_range());
    world.gop.fence();
    return arr;
  };

  ToTArray R = build(outer_tr, /*zero=*/false);

  // gather the disjoint, tile-aligned blocks over outer mode 0 (the external
  // mode), exactly as the scatter branch's per-block evaluate would produce.
  auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // outer elements [0,3)
  auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // outer elements [3,6)

  // PRE-SIZE from the FIRST block partial: widen its OUTER mode 0 (sliced to
  // [0,3)) to the carrier's full mode-0 tiling. R itself is a ToT carrier of
  // the full external mode (axis_src.is<this_type>() branch); axis_src_mode =
  // 0.
  sequant::ResultPtr const first = eval_result<ResultToT>(b0);
  sequant::ResultPtr const carrier = eval_result<ResultToT>(R);
  auto rdest = first->pre_sized_zeros_over_mode(/*mode=*/0, *carrier,
                                                /*axis_src_mode=*/0);
  // the pre-sized destination spans the FULL external mode (mode 0: {0,3,6})
  REQUIRE(rdest->get<ToTArray>().trange().dim(0).tile_extent() == 2);
  REQUIRE(rdest->get<ToTArray>().trange().dim(1).tile_extent() == 1);

  // SCATTER every block into its disjoint slice of the pre-sized destination.
  rdest->write_into_slice(*eval_result<ResultToT>(b0), 0, 0, 3);
  rdest->write_into_slice(*eval_result<ResultToT>(b1), 0, 3, 6);

  auto const& out = rdest->get<ToTArray>();
  ToTArray diff;
  diff("i,j;a,b") = out("i,j;a,b") - R("i,j;a,b");
  REQUIRE(Catch::Approx(diff("i,j;a,b").dot(diff("i,j;a,b"))) == 0.0);
  // guard against a trivial all-zero pass (e.g. a silently no-op'd scatter).
  REQUIRE(out("i,j;a,b").dot(out("i,j;a,b")) > 0.0);
}

TEST_CASE("batched_eval_external_proto_occ_scatter",
          "[eval][batched-external]") {
  // D2.2 regression: the external-mode scatter locates and slices an external
  // occupied index that a CSV/PNO-CCk giant carries ONLY as a protoindex of a
  // composite PNO leg `a<i,j>` -- never written as a plain top-level slot on
  // the giant or on the operands that build it. The giant here is the
  // particle-particle-ladder-shaped intermediate
  //   W{a1<i,j>,a2<i,j>} = (g{m1;m2} * C{a1<i,j>;m2}) * C{a2<i,j>;m1}
  // whose operand legs carry the occ i_1,i_2 only inside the composite `a`.
  //
  // WHY NO PROTO-AWARE LOCATOR IS NEEDED: occupied indices are GUARANTEED to be
  // OUTER modes of a nested (tensor-of-tensor) DistArray -- only the CSV/PNO
  // modes (the composite's base index `a`) are inner. Canonicalization promotes
  // such a proto occ to a plain, top-level, NON-proto entry of canon_indices()
  // (annot() = "i_1,i_2;a_1i_1i_2,a_2i_1i_2"), so index_position() already
  // finds it at its physical outer mode and the scatter slices it as an
  // ordinary outer mode -- there is no "proto-only" blind spot. D1 (the DP) now
  // emits BatchModeType::External for exactly these occ; THIS scatter is what
  // realizes the DP's modeled memory bound at runtime. The test pins that
  // end-to-end: the scatter fires over the occ, the assembled result
  // reconstructs the unbatched reference EXACTLY, and the per-block footprint
  // scales ~block/extent. It also guards the 2026-07-11 defect-1 over-throw (a
  // plain PAO leg alongside composite proto legs, or a leaf that does not carry
  // the occ at all): locating the occ returns nullopt there instead of
  // asserting.
  using sequant::evaluate;
  using sequant::Index;
  using sequant::index_position;
  using sequant::make_batched_custom_evaluator;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using ToTArray = TA::DistArray<TA::Tensor<TA::Tensor<double>>>;

  auto& world = TA::get_default_world();
  // occupied multi-tiled (12 in tiles of 4 -> 3 tiles) so the proto occ mode
  // genuinely partitions into > 1 block; virtual/pao single-tiled.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/12,
                                                    /*nvirt=*/4, /*naux=*/8};
  yield_.set_max_tile(4);

  // g is a flat (mu~) operand carrying NO occ; the two C legs are composite
  // (a<i_1,i_2>) carrying the occ only as protos.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{m_1;m_2} * C{m_2;a1<i_1,i_2>}) * C{a2<i_1,i_2>;m_1}");
  auto node = eval_node(expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const occ =
      sequant::get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](Index const& ix) {
    return ix.space() == occ && !ix.has_proto_indices();
  };

  // The external occ mode, taken from the giant's own (promoted) outer canon
  // indices so the Index matches exactly.
  Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  // The proto occ IS locatable: promoted to a plain outer canon index, so
  // index_position finds it at its physical outer mode -- no blind spot.
  REQUIRE(index_position(node, mode).has_value());

  // Mixed-leaf guard (2026-07-11 defect-1): the C leaf carries a plain PAO
  // outer leg (m_1/m_2) ALONGSIDE the composite a<i,j> whose proto is the occ;
  // the g leaf carries NO occ at all. Locating the occ must resolve on the C
  // leaf and return nullopt on the g leaf -- never over-throw on a leaf that
  // does not carry the target.
  node_t const* c_leaf = nullptr;
  node_t const* g_leaf = nullptr;
  auto find_leaves = [&](auto&& self, node_t const& n) -> void {
    if (n.leaf()) {
      if (n->is_tensor()) {
        auto const lbl = n->as_tensor().label();
        if (lbl == L"C" && !c_leaf) c_leaf = &n;
        if (lbl == L"g" && !g_leaf) g_leaf = &n;
      }
      return;
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_leaves(find_leaves, node);
  REQUIRE(c_leaf != nullptr);
  REQUIRE(g_leaf != nullptr);
  // mixed leaf carries the occ as an outer mode (no over-throw)...
  CHECK(index_position(*c_leaf, mode).has_value());
  // ...and find_leaf_carrying descends past the g leaf to it, returning a
  // physical outer mode (the mode fed to slice_mode / mode_batches).
  auto const carrier = sequant::find_leaf_carrying(node, mode);
  REQUIRE(carrier.has_value());
  CHECK(index_position(carrier->first, mode) == carrier->second);
  // a leaf that does not carry the occ resolves to nullopt (no over-throw).
  CHECK_FALSE(index_position(*g_leaf, mode).has_value());

  // Stamp External on every node whose result carries the occ (the root and the
  // inner g*C product), as the optimizer would for a forest-level external
  // mode.
  node->set_batched_here({{mode, sequant::BatchModeType::External}});
  auto stamp_carriers = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && index_position(n, mode).has_value())
      n->set_batched_here({{mode, sequant::BatchModeType::External}});
    self(self, n.left());
    self(self, n.right());
  };
  stamp_carriers(stamp_carriers, node);

  // Reference: plain unbatched evaluation (batched_here ignored without a
  // custom evaluator). Computed first so yield_'s random leaves are cached and
  // reused.
  auto const ref = evaluate(node, target, yield_)->get<ToTArray>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Footprint scaling: the occ outer mode spans the full extent (12) in the
  // unbatched result, and a single tile-aligned block spans only the block
  // width (4) -- the ~block/extent shrink the scatter buys per block.
  auto const dest_mode = *index_position(node, mode);
  REQUIRE(ref.trange().dim(dest_mode).extent() == 12);

  // Spy scope-guard: records the block count each time the scatter fires.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, [](Index const&) -> std::size_t { return 4; }, accept_occ, spy,
      sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<ToTArray>();

  // Exactness: scattering the disjoint occ blocks reconstructs the whole ToT.
  ToTArray diff;
  diff(target) = ref(target) - res(target);
  REQUIRE(diff(target).dot(diff(target)) == Catch::Approx(0.0).margin(1e-20));
  // guard against a trivial all-zero pass.
  REQUIRE(res(target).dot(res(target)) > 0.0);

  // The scatter genuinely blocked the proto occ: fired, > 1 block (occ extent
  // 12, block width 4 -> 3 blocks at the root).
  REQUIRE_FALSE(guard_calls.empty());
  for (auto const n : guard_calls) CHECK(n > 1);
  CHECK(guard_calls.front() == 3);
}

TEST_CASE("batched_eval_external_two_occ", "[eval][batched-external]") {
  // Task 9 (rank-2 multi-mode): batched_eval_external_axis_scatter stamps a
  // SINGLE external mode External; this test stamps TWO DISTINCT occupied
  // indices External on the SAME (only) product node and proves the runtime
  // nests both scatter loops as a PRODUCT of block counts, not just one
  // mode's worth. Structurally this is the External/scatter analog of
  // "eval_batched_custom_evaluator nests two modes on one node" (which does
  // the same two-modes-per-node nesting for BatchModeType::Contracted).
  //
  // R{i_1;i_2} = g{i_1;a_1} * h{a_1;i_2}: i_1 is free on g and the result
  // only (not shared with h); i_2 is free on h and the result only (not
  // shared with g); a_1 is contracted between g and h. Because i_1 and i_2
  // each appear on only ONE operand (as ordinary result indices), plain
  // ExprPtr binarization already keeps them free -- unlike the aux
  // hyperindex in batched_eval_external_axis_scatter (shared by > 2 tensor
  // slots there), no ResultExpr pinning is needed to keep them external.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occupied multi-tiled (12 in tiles of 4 -> 3 tiles) so BOTH occ
  // external modes i_1, i_2 partition into > 1 block; virtual single-tiled
  // (4, the contracted mode a_1).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/12,
                                                    /*nvirt=*/4};
  yield_.set_max_tile(4);

  auto const expr =
      sequant::deserialize<sequant::ExprPtr>(L"(g{i_1;a_1} * h{a_1;i_2})");
  std::string const target = "i_1,i_2";
  auto node = eval_node(expr);  // mutable: batch modes stamped below

  auto const occ =
      sequant::get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  // The two distinct occupied result indices, taken from the node's own
  // canonical (free/result) indices -- both must be free on the node.
  sequant::container::svector<sequant::Index> occ_axes;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) occ_axes.push_back(ix);
  REQUIRE(occ_axes.size() == 2);
  REQUIRE(occ_axes[0] != occ_axes[1]);
  REQUIRE(sequant::index_position(node, occ_axes[0]).has_value());
  REQUIRE(sequant::index_position(node, occ_axes[1]).has_value());

  // Stamp BOTH occupied indices External on the single product node.
  node->set_batched_here({{occ_axes[0], sequant::BatchModeType::External},
                          {occ_axes[1], sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation (the OFF path -- batched_here are
  // ignored without a custom evaluator). Computed first so yield_'s random
  // leaf arrays are generated and cached, then reused by the batched run.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Distinct per-mode target sizes (as in "nests two modes on one node") so
  // the two nested levels realize DIFFERENT batch counts and their product
  // is unambiguous: occ_axes[0] -> target 4 over extent 12 = 3 blocks;
  // occ_axes[1] -> target 8 = 2 blocks (two tiles fused, then the
  // remainder). Both target sizes are strictly below the mode extent (12),
  // i.e. the ON configuration for both modes.
  auto const& axis0 = occ_axes[0];
  auto target_batch_size = [axis0](sequant::Index const& ix) -> std::size_t {
    return ix == axis0 ? std::size_t{4} : std::size_t{8};
  };

  // Firing-witness: records the block count each time the scatter evaluator
  // fires, and the product of the two live (outer, inner) counts whenever
  // both levels are simultaneously live -- the RED/GREEN proof that the two
  // External loops genuinely NEST as a PRODUCT rather than only one mode
  // firing or the second mode being silently skipped.
  struct GuardState {
    std::vector<std::size_t> live;
    std::vector<std::size_t> counts;
    std::vector<std::size_t> products_at_depth2;
    std::size_t max_depth = 0;
  } state;
  struct TrackingGuard {
    GuardState* st;
    TrackingGuard(GuardState* s, std::size_t n) : st(s) {
      st->counts.push_back(n);
      st->live.push_back(n);
      st->max_depth = std::max(st->max_depth, st->live.size());
      if (st->live.size() == 2)
        st->products_at_depth2.push_back(st->live[0] * st->live[1]);
    }
    TrackingGuard(TrackingGuard const&) = delete;
    TrackingGuard& operator=(TrackingGuard const&) = delete;
    ~TrackingGuard() { st->live.pop_back(); }
  };
  auto make_tracking_guard = [&state](std::size_t n) {
    return TrackingGuard(&state, n);
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, target_batch_size, accept_occ, make_tracking_guard,
      sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: scattering disjoint blocks over both modes reconstructs the
  // whole result.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // Depth-2 nesting engaged: occ_axes[0] sliced at the outer level (3
  // blocks), occ_axes[1] at the inner level (2 blocks) WITHIN each outer
  // block.
  REQUIRE(state.max_depth == 2);
  REQUIRE(state.counts.size() == 4);  // 1 outer firing + 3 inner firings
  CHECK(std::count(state.counts.begin(), state.counts.end(), 3u) == 1);
  CHECK(std::count(state.counts.begin(), state.counts.end(), 2u) == 3);
  // Every depth-2 instant multiplies 3 (outer) by 2 (inner) = 6: BOTH
  // external loops were live together with DISTINCT partitions, i.e. the
  // PRODUCT of block counts, not just one mode's worth.
  REQUIRE(state.products_at_depth2.size() == 3);
  for (auto const p : state.products_at_depth2) CHECK(p == 6);
  CHECK(state.live.empty());
}

TEST_CASE("batched_eval_external_hadamard", "[eval][batched-external]") {
  // Task 9 (genuine Hadamard external): guards the "slice EVERY carrying
  // operand to the same block" contract. i_3 is shared ELEMENTWISE (not
  // summed) between BOTH operands of ONE product node and survives to the
  // result -- a true Hadamard mode, not merely a free index that happens to
  // pass through a single operand untouched (that weaker shape is what the
  // free-only fixtures, e.g. batched_scratch_no_seed_external's g/h/t tree,
  // already exercise). A bug that sliced only ONE of the two operands
  // carrying i_3 (leaving the other at its full extent) would make the
  // per-block elementwise product over i_3 mismatched/wrong for THIS
  // fixture, whereas a free-only fixture (where only one operand ever
  // carries the mode) cannot distinguish "slice every carrying operand" from
  // "slice the first carrying operand" -- there is only ever one.
  //
  // R{i_3;a_1,a_2} = (g{i_3;a_1} * h{i_3;a_2}): i_3 appears on g, on h, AND
  // on the result, contracted at no node (a pure elementwise/outer product
  // over i_3 combined with a_1 from g and a_2 from h). Because i_3 is
  // repeated across exactly 2 operand tensors (the is_valid cap is > 2), it
  // is syntactically legal directly on the occupied space; plain binarize
  // would default to CONTRACTING a repeated-and-otherwise-unlisted index, so
  // the ResultExpr LHS pins i_3 into R's free/output indices (the same
  // technique batched_eval_external_axis_scatter uses for its aux
  // external).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occupied multi-tiled (12 in tiles of 4 -> 3 tiles) so the Hadamard
  // external i_3 genuinely partitions into > 1 block; virtual single-tiled
  // (4).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/12,
                                                    /*nvirt=*/4};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{;a_1,a_2;i_3} = (g{;a_1;i_3} * h{;a_2;i_3})");
  auto node = eval_node(res_expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const occ =
      sequant::get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  // The Hadamard external i_3, taken from the node's own free/result outer
  // modes.
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on result

  // Confirm i_3 is genuinely carried by BOTH operands (not just one) -- the
  // structural property that makes this fixture "genuinely Hadamard" rather
  // than free-only.
  REQUIRE(sequant::index_position(node.left(), mode).has_value());
  REQUIRE(sequant::index_position(node.right(), mode).has_value());

  node->set_batched_here({{mode, sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation (the OFF path).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Firing-witness spy, as in batched_eval_external_axis_scatter: records the
  // block count each time the evaluator fires. A silently-ignored External
  // stamp (or a stamp that degraded to a no-op because only one operand got
  // sliced and the runtime bailed) would never fire this.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  // Two block sizes, both strictly below i_3's extent (12): 4 -> 3 blocks,
  // 8 -> 2 blocks (the ON configuration for both).
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    guard_calls.clear();
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_occ, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

    // Exactness: if the runtime had sliced only ONE of the two Hadamard
    // operands, the per-block elementwise product over i_3 would combine
    // mismatched slices and this would fail (or the two mismatched ranges
    // would not even align into a well-formed contraction).
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // The ON path blocked: the evaluator fired at least once, partitioning
    // into > 1 block.
    REQUIRE_FALSE(guard_calls.empty());
    for (auto const n : guard_calls) CHECK(n > 1);
    std::size_t const expected_blocks = target_batch_size == 4 ? 3 : 2;
    CHECK(guard_calls.front() == expected_blocks);
  }
}

TEST_CASE("batched_eval_external_nested_contracted",
          "[eval][batched-external]") {
  // Task 9 (bonus, carried from Task 5's review): an External mode on an
  // OUTER node nested with a Contracted mode on an INNER node -- the mirror
  // image of "eval_batched_custom_evaluator nests inner mode" (which nests
  // Contracted-over-Contracted). Here the SCATTER branch (Task 5) at the
  // root must re-enter and let a plain accumulate (Contracted) fire WITHIN
  // each external block, i.e. "for external-block: for contracted-block:
  // accumulate, then scatter the block". This exercises the opposite
  // composition order from batched_eval_external_axis_scatter (which nests
  // External outside External at two DIFFERENT nodes, never mixing kinds).
  //
  // inner = u{i_1;i_2;x_1,x_2} * v{i_1;i_3;x_2}: contracts i_1 and x_2
  // (both present on both operands, absent from the term's declared free
  // set), leaving free x_1, i_2, i_3.
  // root = inner * p{i_3;i_6;x_1}: contracts i_3 (present on inner and p,
  // summed), leaving free x_1, i_2, i_6. x_1 is present on BOTH root
  // operands (inner and p) and is declared part of R's free indices, so it
  // is EXTERNAL (scattered) at the root -- a plain ExprPtr would instead
  // default to contracting a doubly-repeated, undeclared index, so (as in
  // batched_eval_external_axis_scatter) the ResultExpr LHS pins x_1 free.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ/virt single-tiled (4, unused as batch modes here); aux multi-tiled
  // (12 in tiles of 4 -> 3 tiles) so both x_1 (outer, External) and x_2
  // (inner, Contracted) partition into > 1 block.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{i_2;i_6;x_1} = ((u{i_2;i_1;x_1,x_2} * v{i_1;i_3;x_2})"
      L" * p{i_3;i_6;x_1})");
  auto node = eval_node(res_expr);  // mutable: batch modes stamped below
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The External mode x_1, taken from the ROOT's own free/result outer
  // modes (a plain outer mode, no protos).
  sequant::Index x1_axis;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      x1_axis = ix;
      break;
    }
  REQUIRE(x1_axis.nonnull());
  REQUIRE(sequant::index_position(node, x1_axis).has_value());  // free on R

  // The inner node (u*v) carries a DIFFERENT aux mode (x_2), contracted
  // there -- locate it the same way "nests inner mode" does, by its
  // contracted aux index rather than tree position (robust to binarize's
  // operand ordering).
  node_t* inner = nullptr;
  std::optional<sequant::Index> x2_axis;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        inner = &n;
        x2_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(x2_axis.has_value());
  REQUIRE(*x2_axis != x1_axis);  // the two nested modes differ

  // Stamp the root's external mode External, the inner's contracted mode
  // Contracted.
  node->set_batched_here({{x1_axis, sequant::BatchModeType::External}});
  (*inner)->set_batched_here({{*x2_axis, sequant::BatchModeType::Contracted}});

  // Reference: plain unbatched evaluation (the OFF path). Computed first so
  // yield_'s random leaf arrays are generated and cached, then reused below.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Firing-witness spy: records the block count each time EITHER level
  // fires (the root's scatter or the inner's accumulate). Distinct target
  // sizes for x_1 (outer, External) vs x_2 (inner, Contracted) make the two
  // levels' partitions unambiguous: x_1 -> target 4 over extent 12 = 3
  // blocks; x_2 -> target 8 = 2 blocks.
  auto target_batch_size = [x1_axis](sequant::Index const& ix) -> std::size_t {
    return ix == x1_axis ? std::size_t{4} : std::size_t{8};
  };
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, target_batch_size, accept_aux, spy, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: External-scatter-outside-Contracted-accumulate composes to
  // the same result as the unbatched contraction.
  TArrayD diff;
  diff("i,j,k") = ref("i,j,k") - res("i,j,k");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // Depth-2 nesting engaged: the root fires once over x_1 (3 blocks) and,
  // per x_1 block, the inner fires once over x_2 (2 blocks) -- 1 + 3 = 4
  // firings, counts a permutation of {3, 2, 2, 2}. Before Task 5/6 this
  // composition (External outside Contracted) did not exist, so this is the
  // RED/GREEN witness that the scatter branch's per-block reinstall lets an
  // inner Contracted mode nest exactly as it does under an outer Contracted
  // mode.
  REQUIRE(guard_calls.size() == 4);
  CHECK(std::count(guard_calls.begin(), guard_calls.end(), 3u) == 1);
  CHECK(std::count(guard_calls.begin(), guard_calls.end(), 2u) == 3);
}

TEST_CASE(
    "eval_batched_custom_evaluator nested scope guards compose "
    "multiplicatively",
    "[eval]") {
  // Task 4.3: prove that screening relaxation composes over nested batch
  // levels. The re-entrant inner evaluator (Task 4.2) is built by threading
  // make_scope_guard (along with accept/is_volatile/persistent_only/
  // target_batch_size) INTO the nested make_batched_custom_evaluator call
  // unchanged -- so the inner level constructs its own guard via
  // make_scope_guard(inner_batches), not a no-op. The outer scope_guard is a
  // local RAII variable held for the outer level's entire batch loop, which
  // includes the per-batch evaluate() calls that trigger the inner re-entry;
  // so when the inner guard is constructed, the outer guard is STILL ALIVE.
  // A backend guard that relaxes block-sparse screening scaled by its own
  // level's batch count therefore composes MULTIPLICATIVELY across nesting:
  // net relaxation = outer_batches * inner_batches, matching "a contribution
  // significant over the full product of batch modes must not be screened
  // away in an individual (outer-cell, inner-cell) cell."
  //
  // This test cannot exercise REAL block-sparse screening -- the TA eval
  // tests here use dense TensorD, which has no SparseShape to relax. Instead
  // it proves the STRUCTURAL composition: a custom ScopeGuardFactory whose
  // RAII guard records, on construction/destruction, the batch count it was
  // built with against a shared "currently alive" stack. If both guards are
  // ever alive at once with the stack holding {outer_n, inner_n}, that is
  // exactly the multiplicative-composition invariant a real backend guard
  // would exploit. Numeric validation against an actual screening threshold
  // is deferred to the Phase 6 end-to-end MPQC run, where a real
  // TiledArray-SparseShape-backed guard exists to relax.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Same depth-2 nesting shape as the "nests inner mode" test above: occ
  // single-tiled (4), aux multi-tiled (12 in tiles of 4 -> 3 tiles) so both
  // batch modes realize 3 batches each.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((u{i_1;i_2;x_1,x_2} * v{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_space_registry()->retrieve(L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_batched_here({{*root_axis, sequant::BatchModeType::Contracted}});

  node_t* inner = nullptr;
  std::optional<sequant::Index> inner_axis;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        inner = &n;
        inner_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(inner_axis.has_value());
  REQUIRE(*inner_axis != *root_axis);
  (*inner)->set_batched_here(
      {{*inner_axis, sequant::BatchModeType::Contracted}});

  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Shared recording state: a stack of currently-alive guards' batch counts
  // (pushed on construction, popped on destruction). Whenever a guard is
  // constructed while another is already alive (stack depth becomes 2), the
  // product of the two counts is recorded -- that product is the net
  // relaxation factor a composing backend guard would apply at that instant.
  struct GuardState {
    std::vector<std::size_t> live;
    std::size_t max_depth = 0;
    std::vector<std::size_t> products_at_depth2;
  } state;

  struct TrackingGuard {
    GuardState* st;
    TrackingGuard(GuardState* s, std::size_t n) : st(s) {
      st->live.push_back(n);
      st->max_depth = std::max(st->max_depth, st->live.size());
      if (st->live.size() == 2)
        st->products_at_depth2.push_back(st->live[0] * st->live[1]);
    }
    TrackingGuard(TrackingGuard const&) = delete;
    TrackingGuard& operator=(TrackingGuard const&) = delete;
    ~TrackingGuard() { st->live.pop_back(); }
  };
  auto make_tracking_guard = [&state](std::size_t n) {
    return TrackingGuard(&state, n);
  };

  auto cache = cache_t::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, make_tracking_guard, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: nesting the guard-instrumented evaluator changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // The real deliverable: both guards were alive simultaneously (depth 2
  // reached) at least once, and every such simultaneity records the product
  // of outer_batches (3, over x_1) and inner_batches (3, over x_2) -- 9, the
  // net relaxation a composing backend would apply. This happens once per
  // outer batch (3 outer batches), matching the "nests inner mode" test's
  // guard_calls count of 4 total firings (1 outer + 3 inner).
  REQUIRE(state.max_depth == 2);
  REQUIRE(state.products_at_depth2.size() == 3);
  for (auto const p : state.products_at_depth2) CHECK(p == 9);
  // All guards are popped by the end: no leaked / mismatched RAII scope.
  CHECK(state.live.empty());
}

TEST_CASE("make_evaluator BatchPolicy adapter", "[eval]") {
  // Strong equivalence: make_evaluator(policy, yielder) must produce the same
  // numerical result as a hand-built make_batched_custom_evaluator with the
  // matching Tensor->EvalNode volatile lift, on the same expression.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using sequant::make_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Multi-tile unoccupied space so batching actually engages (3 tiles of <=4).
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // Contracts a1,a2 (unoccupied) -> batch mode is an unoccupied index (3
  // tiles). The subtree contains a "t" leaf, which the policy marks volatile.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto const node = eval_node(expr);

  // Reference (no batching).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Build a BatchPolicy: batch over any index at size 4; "t" is volatile;
  // persistent_only=true turns the volatility gate ON (default is across the
  // board).
  sequant::BatchPolicy policy{
      .is_batchable_contracted_index =
          [](sequant::Index const&) { return true; },
      .batch_target_size = [](sequant::Index const&) -> std::size_t {
        return 4;
      },
      .is_volatile_leaf =
          [](sequant::Tensor const& t) { return t.label() == L"t"; },
      .persistent_only = true};

  // (1) make_evaluator with volatile "t" gate (persistent_only=true): batching
  // should be DECLINED (same as the persistence gate test above for
  // make_batched_custom_evaluator).
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_evaluator(policy, yield_, spy));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE_FALSE(batched);  // volatile subtree -> not batched
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }

  // (2) Strong equivalence: make_evaluator (no volatile gate: empty predicate)
  // == hand-built make_batched_custom_evaluator, same result as reference.
  {
    sequant::BatchPolicy policy_nv{
        .is_batchable_contracted_index = policy.is_batchable_contracted_index,
        .batch_target_size = policy.batch_target_size,
        .is_volatile_leaf = {}  // no volatile gate
    };

    // make_evaluator path
    auto cache_me = cache_t::empty();
    cache_me.set_custom_evaluator(make_evaluator(policy_nv, yield_));
    auto const res_me =
        evaluate(node, target, yield_, cache_me)->get<TArrayD>();

    // hand-built path
    auto hand_nv = [](node_t const&) -> bool { return false; };
    auto cache_hb = cache_t::empty();
    cache_hb.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, policy_nv.batch_target_size, policy_nv.is_batchable_index(),
        sequant::make_no_scope_guard{}, hand_nv));
    auto const res_hb =
        evaluate(node, target, yield_, cache_hb)->get<TArrayD>();

    // Both must match the reference (batched summation -> loose FP tolerance).
    REQUIRE(equal_tarrays<Loose>(res_me, ref));
    REQUIRE(equal_tarrays<Loose>(res_hb, ref));
    // And agree with each other to the same loose tolerance: both come from
    // the same batched-summation algorithm, whose accumulation order is
    // thread-non-deterministic, so the two independent evaluations agree only
    // up to FP noise (a Tight/exact compare here flakes by a few ULPs).
    REQUIRE(equal_tarrays<Loose>(res_me, res_hb));
  }
}

TEST_CASE("shape_spike_tot_general_product", "[shape-spike]") {
  // De-risk spike: can (A(la) * B(ra)).set_shape(s) for a ToT general product
  // (a) evaluate without throwing, (b) honor the imposed SparseShape (zeroed
  // tiles stay zero), and (c) produce the same inner data as TA::einsum for
  // the surviving tiles?
  //
  // Uses the same ToT*ToT->ToT annotation as the ToT_times_ToT_to_ToT section
  // of eval_with_tiledarray, but with SparsePolicy and a multi-tile outer
  // TiledRange so there are enough outer tiles to zero one.
  //
  // Annotation (same structure as the existing ToT*ToT test):
  //   lhs: "i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"
  //   rhs: "i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"
  //   result: "i_2,i_1;a_1i_1i_2,a_2i_1i_2"
  // Contraction: outer i_3, inner a_4.  Fused outer (Hadamard): i_1, i_2.
  // Result outer: (i_2, i_1), each with 2 tiles -> 4 outer tiles total.

  using ToTArray =
      TA::DistArray<TA::Tensor<TA::Tensor<double>>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  // 2 tiles per occ mode; 1 tile per virt mode.
  // occ extent = 4, tile size = 2; virt extent = 3.
  TA::TiledRange1 const tr1_occ{0, 2, 4};
  std::size_t const virt_ext = 3;

  // Outer TiledRange for lhs and rhs: 3D (i1, i2, i3) all occ.
  TA::TiledRange const outer_tr{tr1_occ, tr1_occ, tr1_occ};

  // Inner tensor Range:
  //   lhs inner = (a4, a1), each virt_ext = (3, 3)
  //   rhs inner = (a2, a4), each virt_ext = (3, 3)
  TA::Range const inner_r{virt_ext, virt_ext};

  // Build a sparse ToT array with all outer tiles populated (norm > threshold).
  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    // "Dense" SparseShape: every tile has scaled norm 1.0.
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        // Capture ir by value: the lambda must not capture it by reference
        // since the taskq may run after this scope.
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto lhs = make_tot(outer_tr, inner_r);
  auto rhs = make_tot(outer_tr, inner_r);

  // Matching the existing ToT*ToT test annotations.
  std::string const la{"i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"};
  std::string const ra{"i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"};
  std::string const ca{"i_2,i_1;a_1i_1i_2,a_2i_1i_2"};

  // Baseline: einsum.
  ToTArray C0 = TA::einsum(lhs(la), rhs(ra), ca);
  world.gop.fence();
  REQUIRE_FALSE(C0.is_initialized() == false);

  // Result outer TiledRange: 2D (i_2, i_1), both with {0,2,4} -> 2x2=4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // Construct SparseShape that zeros outer tile (0,0) of the result.
  // tile_norms uses do_not_scale=true: the values are treated as already-scaled
  // (so a value of 0.0 is exactly zero and >0 is nonzero).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};
  // s must outlive the assignment below (set_shape stores &s internally).

  // Standard-layer evaluation with imposed shape.
  ToTArray C1;
  REQUIRE_NOTHROW(C1(ca) = (lhs(la) * rhs(ra)).set_shape(s));
  world.gop.fence();

  // Outcome (1): shape constraint honored -- tile (0,0) must be zero/absent.
  REQUIRE(C1.is_zero({0, 0}));

  // Outcome (2): surviving outer tiles must match the einsum baseline.
  // Walk C1's local tiles; for each nonzero outer tile, compare inner data
  // element-wise with C0.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    auto const& outer0 = C0.find(idx).get();  // TA::Tensor<TA::Tensor<double>>
    auto const& outer1 = it->get();
    REQUIRE(outer0.range() == outer1.range());

    for (std::size_t o = 0; o < outer0.size(); ++o) {
      auto const& inner0 = outer0[o];
      auto const& inner1 = outer1[o];
      REQUIRE(inner0.range() == inner1.range());
      for (std::size_t k = 0; k < inner0.size(); ++k) {
        CHECK(inner0[k] == Catch::Approx(inner1[k]));
      }
    }
  }
}

TEST_CASE("shape_spike_T_times_ToT_general_product", "[shape-spike]") {
  // Case A: T x ToT -> ToT (mixed operand). Models a dressed DF integral g
  // (flat T) contracted with a PNO coefficient C (ToT) to yield a ToT result.
  // Confirms (T_op(la) * ToT_op(ra)).set_shape(s) for this mixed general
  // product (a) evaluates, (b) honors the imposed SparseShape, (c) matches the
  // TA::einsum baseline on surviving tiles.
  //
  // Annotation mirrors the existing T_times_ToT_to_ToT section:
  //   T_op (flat):  "i_3,i_1"
  //   ToT_op (ToT): "i_2,i_3;a_3i_2i_3,a_4i_2i_3"
  //   result (ToT): "i_1,i_2;a_3i_2i_3,a_4i_2i_3"
  // Shared/contracted outer mode: i_3.  Free outer: i_1 (from T), i_2 (from
  // ToT).  Result outer (i_1, i_2), each with 2 tiles -> 4 outer tiles.

  using FlatArray = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using ToTArray =
      TA::DistArray<TA::Tensor<TA::Tensor<double>>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  TA::TiledRange1 const tr1_occ{0, 2, 4};  // 2 tiles, occ extent 4
  std::size_t const virt_ext = 3;

  // T_op outer: (i_3, i_1), both occ.
  TA::TiledRange const t_outer_tr{tr1_occ, tr1_occ};
  // ToT_op outer: (i_2, i_3), both occ; inner (a_3, a_4) each virt.
  TA::TiledRange const tot_outer_tr{tr1_occ, tr1_occ};
  TA::Range const inner_r{virt_ext, virt_ext};

  // Build a flat sparse T array with all tiles populated.
  auto make_flat = [&](TA::TiledRange const& otr) -> FlatArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    FlatArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        *it = world.taskq.add(
            [](TA::Range const& rng) { return random_tensor<double>(rng); },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  // Build a sparse ToT array with all outer tiles populated.
  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto T_op = make_flat(t_outer_tr);
  auto ToT_op = make_tot(tot_outer_tr, inner_r);

  std::string const la{"i_3,i_1"};
  std::string const ra{"i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
  std::string const ca{"i_1,i_2;a_3i_2i_3,a_4i_2i_3"};

  // Baseline: einsum.
  ToTArray C0 = TA::einsum(T_op(la), ToT_op(ra), ca);
  world.gop.fence();

  // Result outer TiledRange: 2D (i_1, i_2) -> 2x2 = 4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // SparseShape zeroing outer tile (0,0).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};

  // Standard-layer evaluation with imposed shape.
  ToTArray C1;
  REQUIRE_NOTHROW(C1(ca) = (T_op(la) * ToT_op(ra)).set_shape(s));
  world.gop.fence();

  // (1) shape constraint honored.
  REQUIRE(C1.is_zero({0, 0}));

  // (2) surviving tiles match einsum baseline.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    auto const& outer0 = C0.find(idx).get();
    auto const& outer1 = it->get();
    REQUIRE(outer0.range() == outer1.range());

    for (std::size_t o = 0; o < outer0.size(); ++o) {
      auto const& inner0 = outer0[o];
      auto const& inner1 = outer1[o];
      REQUIRE(inner0.range() == inner1.range());
      for (std::size_t k = 0; k < inner0.size(); ++k) {
        CHECK(inner0[k] == Catch::Approx(inner1[k]));
      }
    }
  }
}

TEST_CASE("shape_spike_ToT_inner_contraction_to_flat_T", "[shape-spike]") {
  // Case B: ToT x ToT inner-contraction -> flat T (the dot_inner / DeNest::True
  // path; see result.hpp:581 TA::einsum<TA::DeNest::True>). Two ToT operands
  // share their outer (occ) modes (Hadamard) and fully contract their inner
  // (composite) modes, denesting to a plain tensor-of-scalars result.
  //
  // The standard-expression-layer equivalent of einsum<DeNest::True> is the
  // .dot_inner() expression (TA einsum/tiledarray.h:603):
  //   C(c) = A(a + inner.a).dot_inner(B(b + inner.b));
  // DotInnerExpr derives from Expr, so it exposes set_shape(); the spike here
  // is whether that override is honored for the denested flat-T result.
  //
  // Annotation:
  //   lhs (ToT): "i_1,i_2;a_1i_1i_2,a_2i_1i_2"
  //   rhs (ToT): "i_1,i_2;a_1i_1i_2,a_2i_1i_2"
  //   result (flat T): "i_1,i_2"   (inner a_1,a_2 dotted away)
  // Outer i_1,i_2 are Hadamard (survive); result outer (i_1,i_2) -> 2x2 tiles.

  using FlatArray = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using ToTArray =
      TA::DistArray<TA::Tensor<TA::Tensor<double>>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  TA::TiledRange1 const tr1_occ{0, 2, 4};  // 2 tiles, occ extent 4
  std::size_t const virt_ext = 3;

  // Both operands: outer (i_1, i_2) occ; inner (a_1, a_2) virt.
  TA::TiledRange const outer_tr{tr1_occ, tr1_occ};
  TA::Range const inner_r{virt_ext, virt_ext};

  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto lhs = make_tot(outer_tr, inner_r);
  auto rhs = make_tot(outer_tr, inner_r);

  // Outer annotations (a/b) and inner annotations (inner.a/inner.b) for
  // .dot_inner(); the result is the bare outer annotation (no inner part).
  std::string const a_outer{"i_1,i_2"};
  std::string const a_inner{";a_1i_1i_2,a_2i_1i_2"};
  std::string const b_outer{"i_1,i_2"};
  std::string const b_inner{";a_1i_1i_2,a_2i_1i_2"};
  std::string const ca{"i_1,i_2"};

  // Baseline: einsum<DeNest::True>. The full ToT annotation carries the inner
  // part; the result annotation is outer-only (denested).
  FlatArray C0 = TA::einsum<TA::DeNest::True>(lhs(a_outer + a_inner),
                                              rhs(b_outer + b_inner), ca);
  world.gop.fence();

  // Result outer TiledRange: 2D (i_1, i_2) -> 2x2 = 4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // SparseShape zeroing outer tile (0,0).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};

  // Standard-layer denesting product via .dot_inner(), with imposed shape.
  FlatArray C1;
  REQUIRE_NOTHROW(C1(ca) = (lhs(a_outer + a_inner))
                               .dot_inner(rhs(b_outer + b_inner))
                               .set_shape(s));
  world.gop.fence();

  // (1) shape constraint honored.
  REQUIRE(C1.is_zero({0, 0}));

  // (2) surviving flat tiles match the einsum<DeNest::True> baseline.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    // Copy the tiles out BY VALUE: C0.find(idx) and *it are temporary Futures
    // whose payload get() only references, so binding `auto const&` would
    // dangle once those temporaries die at the semicolon (ASan
    // stack-use-after-scope on the tile's Range). A TA::Tensor copy is a cheap
    // shallow handle that keeps its own Range alive for the loop body.
    TA::Tensor<double> const t0 = C0.find(idx).get();
    TA::Tensor<double> const t1 = it->get();
    REQUIRE(t0.range() == t1.range());
    for (std::size_t k = 0; k < t0.size(); ++k) {
      CHECK(t0[k] == Catch::Approx(t1[k]));
    }
  }
}

TEST_CASE("shape_provider_general_product", "[shape-provider]") {
  // Task 2: a provider returning a REAL SparseShape for a ToT*ToT->ToT general
  // product (the `(lhs(la) * rhs(ra)).set_shape(s)` emission form) must
  // constrain the eval result: the zeroed outer tiles are is_zero and the
  // surviving tiles equal the unshaped (einsum) baseline. Also covers the
  // no-op (full-ones shape) and nullopt (decline) cases.
  //
  // Graph: ToT * ToT -> ToT, single internal Product node. SparsePolicy so a
  // SparseShape is meaningful; multi-tiled occ so >=1 outer tile can be zeroed.
  using sequant::CacheManager;
  using sequant::evaluate;
  using sequant::TAEvalContext;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ extent 4 in tiles of <=2 -> 2 outer tiles per occ mode.
  rand_tensor_yield<double, TA::SparsePolicy> yield{world, /*nocc=*/4,
                                                    /*nvirt=*/3};
  yield.set_max_tile(2);

  using ToTArray = typename decltype(yield)::array_tot_type;

  constexpr std::wstring_view expr_str =
      L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}"
      L" * "
      L"s{a2<i1,i2>;a4<i2,i3>}";
  auto const node = eval_node(sequant::deserialize<sequant::ExprPtr>(expr_str));
  std::string const target{"i_2,i_1;a_1i_1i_2,a_2i_1i_2"};

  // Unshaped reference (no hook).
  auto const ref = evaluate(node, target, yield)->get<ToTArray>();
  TA::TiledRange const res_tr = ref.trange();

  auto make_ctx = [](auto shape_fn) {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [shape_fn](
            sequant::FullBinaryNode<sequant::EvalExprTA> const&,
            TA::TiledRange const& tr) -> std::optional<TA::SparseShape<float>> {
      return shape_fn(tr);
    };
    return ctx;
  };

  SECTION("real shape: zeroed tile is_zero, survivors match baseline") {
    // Zero outer tile (0,0); keep the rest.
    auto ctx = make_ctx([](TA::TiledRange const& tr) {
      TA::Tensor<float> norms{tr.tiles_range(), 1.0f};
      norms[{0, 0}] = 0.0f;
      return TA::SparseShape<float>{norms, tr, /*do_not_scale=*/true};
    });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();

    // (1) zeroed tile gone.
    REQUIRE(res.is_zero({0, 0}));
    // (2) survivors equal the unshaped baseline.
    for (auto it = res.begin(); it != res.end(); ++it) {
      if (!res.is_local(it.index()) || res.is_zero(it.index())) continue;
      auto const& idx = it.index();
      REQUIRE_FALSE(ref.is_zero(idx));
      auto const& o0 = ref.find(idx).get();
      auto const& o1 = it->get();
      REQUIRE(o0.range() == o1.range());
      for (std::size_t o = 0; o < o0.size(); ++o) {
        auto const& in0 = o0[o];
        auto const& in1 = o1[o];
        REQUIRE(in0.range() == in1.range());
        for (std::size_t k = 0; k < in0.size(); ++k)
          CHECK(in0[k] == Catch::Approx(in1[k]));
      }
    }
  }

  SECTION("no-op full-ones shape: result equals unshaped") {
    auto ctx = make_ctx([](TA::TiledRange const& tr) {
      return TA::SparseShape<float>{1.0f, tr};
    });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();
    std::string const annot{"i,j;a,b"};
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }

  SECTION("nullopt: hook declines, falls through to unshaped prod()") {
    auto ctx = make_ctx(
        [](TA::TiledRange const&) -> std::optional<TA::SparseShape<float>> {
          return std::nullopt;
        });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();
    std::string const annot{"i,j;a,b"};
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }

  SECTION("no hook installed => unshaped behavior unchanged") {
    auto const res = evaluate(node, target, yield)->get<ToTArray>();
    std::string const annot{"i,j;a,b"};
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }
}

TEST_CASE("shape_provider_denest_to_flat", "[shape-provider]") {
  // Task 2: the DeNest::True path (ToT * ToT -> flat T), emitted as
  // `lhs(la).dot_inner(rhs(ra)).set_shape(s)`. A ResultExpr pins the flat head
  // (bare occ indices, no protos), so node->tot() is false while both operands
  // are ToT => de_nest == True at the binary site. A provider returning a real
  // shape must zero the designated outer tile and leave survivors matching the
  // unshaped (einsum<DeNest::True>) baseline.
  using sequant::CacheManager;
  using sequant::evaluate;
  using sequant::TAEvalContext;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield{world, /*nocc=*/4,
                                                    /*nvirt=*/3};
  yield.set_max_tile(2);

  using FlatArray = typename decltype(yield)::array_type;

  // Result head D{i1,i2} (flat): both operands carry composite virtuals
  // a1<i1,i2>,a2<i1,i2> that fully contract (inner), leaving bare i1,i2.
  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"D{i1,i2} = "
      L"R{a1<i1,i2>,a2<i1,i2>;i1,i2} * g{i1,i2;a1<i1,i2>,a2<i1,i2>}");
  auto const node = eval_node(res_expr);
  // Confirm this really is the denest (DeNest::True) path.
  REQUIRE(node.left()->tot());
  REQUIRE(node.right()->tot());
  REQUIRE_FALSE(node->tot());

  std::string const target{"i_1,i_2"};
  auto const ref = evaluate(node, target, yield)->get<FlatArray>();
  TA::TiledRange const res_tr = ref.trange();

  SECTION("real shape on the denest result: zeroed tile + survivors match") {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [](sequant::FullBinaryNode<sequant::EvalExprTA> const&,
           TA::TiledRange const& tr) -> std::optional<TA::SparseShape<float>> {
      TA::Tensor<float> norms{tr.tiles_range(), 1.0f};
      norms[{0, 0}] = 0.0f;
      return TA::SparseShape<float>{norms, tr, /*do_not_scale=*/true};
    };
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<FlatArray>();

    REQUIRE(res.is_zero({0, 0}));
    for (auto it = res.begin(); it != res.end(); ++it) {
      if (!res.is_local(it.index()) || res.is_zero(it.index())) continue;
      auto const& idx = it.index();
      REQUIRE_FALSE(ref.is_zero(idx));
      auto const& t0 = ref.find(idx).get();
      auto const& t1 = it->get();
      REQUIRE(t0.range() == t1.range());
      for (std::size_t k = 0; k < t0.size(); ++k)
        CHECK(t0[k] == Catch::Approx(t1[k]));
    }
  }

  SECTION("nullopt on the denest result: unchanged") {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [](sequant::FullBinaryNode<sequant::EvalExprTA> const&,
           TA::TiledRange const&) -> std::optional<TA::SparseShape<float>> {
      return std::nullopt;
    };
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<FlatArray>();
    REQUIRE(equal_tarrays(res, ref, "i,j"));
  }
}

// ---------------------------------------------------------------------------
// Node-level external placement correctness reproducer (order_aware_recompute).
//
// Drives the FULL optimizer -> binarize -> stamp_lifetime_masks -> batched TA
// evaluate pipeline (unlike the eval_batched_* cases, which stamp batch modes
// by hand), so the External BatchModeType stamps come from the optimizer's
// node-level placement, not a fixture. Runs the SAME term twice --
// order_aware_recompute = false (root-level "forest seed") and = true
// (node-level placement) -- and compares each batched result against the plain
// unbatched reference (ground truth). Hidden ([.]) diagnostic.
// ---------------------------------------------------------------------------
namespace {

// Copy the batch-annotation fields from a plain EvalExpr binarized tree
// onto the structurally-identical EvalExprTA tree (to_ta_node's tensor-branch
// reconstructor drops them). Both trees are built by the same binarize/
// transform_node post-order, so a lockstep walk aligns node-for-node.
void copy_batch_annotations(
    sequant::FullBinaryNode<sequant::EvalExpr> const& from,
    sequant::FullBinaryNode<sequant::EvalExprTA>& to) {
  const_cast<sequant::EvalExprTA&>(*to).set_batched_here(from->batched_here());
  const_cast<sequant::EvalExprTA&>(*to).set_sliced_modes(from->sliced_modes());
  const_cast<sequant::EvalExprTA&>(*to).set_batch_order_aware(
      from->batch_order_aware());
  const_cast<sequant::EvalExprTA&>(*to).set_batch_effective_count(
      from->batch_effective_count());
  if (!from.leaf()) {
    copy_batch_annotations(from.left(), to.left());
    copy_batch_annotations(from.right(), to.right());
  }
}

// Dump per-node annotations of a TA eval tree (result indices, batched_here
// modes+kind, sliced_modes, order_aware) for eyeballing placement differences.
void dump_annotations(sequant::FullBinaryNode<sequant::EvalExprTA> const& n,
                      int depth = 0) {
  using sequant::BatchModeType;
  std::string pad(2 * depth, ' ');
  std::string res;
  for (auto const& ix : n->canon_indices()) {
    res += sequant::toUtf8(ix.full_label());
    res += ' ';
  }
  std::string bh;
  for (auto const& [ix, knd] : n->batched_here()) {
    bh += sequant::toUtf8(ix.full_label());
    bh += (knd == BatchModeType::External ? ":ext " : ":con ");
  }
  std::string sm;
  for (auto const& ix : n->sliced_modes()) {
    sm += sequant::toUtf8(ix.full_label());
    sm += ' ';
  }
  std::string ident;
  try {
    ident = sequant::toUtf8(sequant::io::serialization::to_string(n->expr()));
  } catch (...) {
    ident = "?";
  }
  std::cerr << pad << (n.leaf() ? "LEAF " : "NODE ") << "res=[" << res
            << "] batched=[" << bh << "] sliced=[" << sm
            << "] oa=" << (n->batch_order_aware() ? 1 : 0)
            << " | expr=" << ident << "\n";
  if (!n.leaf()) {
    dump_annotations(n.left(), depth + 1);
    dump_annotations(n.right(), depth + 1);
  }
}

}  // namespace

TEST_CASE("node_level_external_placement_correctness",
          "[.][eval][batched-external][node-placement]") {
  using namespace sequant;
  using TA::TArrayD;
  using node_t = FullBinaryNode<EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();

  auto isr = get_default_context().index_space_registry();
  auto occ = isr->retrieve(L"i");
  auto virt = isr->retrieve(L"a");
  auto aux = isr->retrieve(L"x");

  // Actual arrays: occ multi-tiled (40 in tiles of 4 -> 10 blocks) so the
  // external occ mode partitions into > 1 batch, and large vs virt/aux (4) so
  // an i_2-carrying intermediate (~i*i) is the peak and slicing i_2 lowers it,
  // forcing node-level placement to adopt i_2 as External.
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/12,
                                                    /*nvirt=*/4, /*naux=*/4};
  yield_.set_max_tile(4);

  // Term with an external occ spectator (i_2, free on the result and carried on
  // internal nodes) plus a subtree invariant to it. Non-left-deep so the right
  // branch (w*p) carries NO i_2 (invariant), the left branch (g*h) carries i_2.
  // Flat 4-tensor product (no inner parens: those would parse as a nested
  // Product of sub-products, which optimize() treats as < 3 factors and leaves
  // unfactorized => no batch annotations). The optimizer factorizes it. i_2,i_6
  // are occ externals (free on the result, carried across internal nodes);
  // x_1,x_2 are contracted aux.
  std::string const expr_str =
      "g{i_2,a_1;x_1} * h{a_1,a_2;x_1} * w{a_2,a_3;x_2} * p{a_3,i_6;x_2}";
  std::string const target = "i_2,i_6";

  auto const is_occ = [occ](Index const& ix) { return ix.space() == occ; };
  auto const is_aux = [aux](Index const& ix) { return ix.space() == aux; };

  // idx_to_extent inflated so the modeled peak far exceeds peak_threshold and
  // the external-emit gate fires on BOTH the legacy (order_aware=false) and the
  // node-level (order_aware=true) path.
  std::function<std::size_t(Index const&)> idxsz =
      [occ, virt, aux](Index const& ix) -> std::size_t {
    if (ix.space() == occ) return 12;
    if (ix.space() == virt) return 4;
    if (ix.space() == aux) return 4;
    return 1;
  };
  std::function<std::size_t(Index const&)> bts =
      [occ](Index const& ix) -> std::size_t {
    return ix.space() == occ ? std::size_t{4} : std::size_t{60};
  };

  auto run = [&](bool order_aware) -> TArrayD {
    auto const expr =
        deserialize<ExprPtr>(std::wstring(expr_str.begin(), expr_str.end()));

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();

    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
    opts.reorder = ReorderSum::NoReorder;
    opts.idx_to_extent = idxsz;
    opts.inner_pow = [](Index const&, std::size_t) -> double { return 1.0; };
    opts.batch_policy.is_batchable_contracted_index = is_aux;
    opts.batch_policy.is_batchable_external_index = is_occ;
    opts.batch_policy.batch_target_size = bts;
    opts.batch_policy.batch_spectator_indices = true;
    opts.batch_policy.order_aware_recompute = order_aware;
    // The sweep var now drives node-level EMISSION too (decoupled from the cost
    // model): order_aware=false => root-seed, true => node-level placement.
    // This test's point is that BOTH emissions give the correct energy.
    opts.batch_policy.node_level_placement = order_aware;
    opts.batch_policy.peak_threshold = 1.0e3;  // tiny budget => force external
    opts.term_batch_axes = axes_map;

    auto optimized = optimize(expr, opts);
    REQUIRE(optimized);
    auto it = axes_map->find(optimized.get());
    REQUIRE(it != axes_map->end());
    auto const& node_axes = it->second;
    std::cerr << "[diag] node_axes.size()=" << node_axes.size() << "\n";
    for (std::size_t j = 0; j < node_axes.size(); ++j) {
      std::cerr << "  node_axes[" << j << "] oa=" << node_axes[j].order_aware
                << " naxes=" << node_axes[j].axes.size() << " :";
      for (auto const& [ix, knd] : node_axes[j].axes)
        std::cerr << " " << toUtf8(ix.full_label())
                  << (knd == BatchModeType::External ? ":ext" : ":con");
      std::cerr << "\n";
    }

    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto enode = binarize(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    auto ta_node = to_ta_node(enode);
    copy_batch_annotations(enode, ta_node);

    // Populate sliced_modes from the External stamps (the real pipeline does
    // this inside CacheManager; the manual empty cache below does not).
    std::array<node_t, 1> forest{ta_node};
    stamp_lifetime_masks(forest);
    ta_node = forest[0];

    std::cerr << "\n=== order_aware_recompute = "
              << (order_aware ? "true" : "false") << " ===\n";
    dump_annotations(ta_node);
    std::cerr << "SCHEDULE_IR_JSON "
              << sequant::eval::schedule_ir_json(ta_node,
                                                 order_aware ? "oa" : "root")
              << "\n";

    auto accept = opts.batch_policy.is_batchable_index();
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, bts, accept, make_no_scope_guard{}, never_volatile{}));
    return evaluate(ta_node, target, yield_, cache)->get<TArrayD>();
  };

  // Ground truth: plain unbatched evaluation (no custom evaluator, batch
  // stamps ignored). Also fills yield_'s random leaf cache first.
  auto const ref_expr =
      deserialize<ExprPtr>(std::wstring(expr_str.begin(), expr_str.end()));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto ref_node = to_ta_node(binarize(ref_expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  auto const ref = evaluate(ref_node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);

  auto const res_off = run(false);
  auto const res_on = run(true);

  TArrayD d_off, d_on;
  d_off("i,j") = ref("i,j") - res_off("i,j");
  d_on("i,j") = ref("i,j") - res_on("i,j");
  double const err_off = TA::norm2(d_off);
  double const err_on = TA::norm2(d_on);
  std::cerr << "\n[node-placement] |ref| = " << TA::norm2(ref)
            << "  err(order_aware=false) = " << err_off
            << "  err(order_aware=true) = " << err_on << "\n";

  CHECK(err_off < 1e-10);
  CHECK(err_on < 1e-10);
}
