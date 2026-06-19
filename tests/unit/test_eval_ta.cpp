#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
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
                              ix.space() == isr->retrieve(L"x"));
               return ix.space() == isr->retrieve(L"i")   ? nocc_
                      : ix.space() == isr->retrieve(L"a") ? nvirt_
                                                          : naux_;
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
    auto const axis = batch_axis(node);
    REQUIRE(axis.has_value());
    REQUIRE(axis.value() == c.front());
  }

  SECTION("two contracted indices") {
    // g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4} sums over a1,a2.
    auto const node = node_of(L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 2);
    auto const axis = batch_axis(node);
    REQUIRE(axis.has_value());
    REQUIRE(ranges::contains(c, axis.value()));
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

  SECTION("predicate scopes the batch axis") {
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
    // mode 1 (b) has 3 tiles of 3 elements each (extent 9).
    // (extra parens: compare as a single bool so Catch2 needn't stringify
    // pairs) target larger than the extent -> a single batch (caller declines).
    REQUIRE((r->mode_batches(1, 100) == batches_t{{0, 9}}));
    // target 4: tile0 (3<4) + tile1 reaches 6>=4 -> [0,6); remainder [6,9).
    REQUIRE((r->mode_batches(1, 4) == batches_t{{0, 6}, {6, 9}}));
    // target 1: every tile is its own batch.
    REQUIRE((r->mode_batches(1, 1) == batches_t{{0, 3}, {3, 6}, {6, 9}}));
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

  // contracts a1,a2 (unoccupied) -> batch axis is an unoccupied index (3 tiles)
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto const node = eval_node(expr);

  // Reference first, so yield_'s (random) leaf arrays are generated and cached;
  // the batched evaluator below copies yield_ and thus reuses the same arrays.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Batched evaluation must reproduce the reference for any target batch size,
  // since sum_K = sum_{batches} sum_{K in batch}. The batch axis is unoccupied
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

  // Contracts a1,a2 (unoccupied, 3 tiles) -> batchable over an unoccupied axis.
  // The subtree contains a "t" leaf, which we treat as volatile.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto const node = eval_node(expr);
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Volatile-leaf predicate: the amplitude "t".
  auto is_volatile_t = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  // (1) Gate ON: the subtree contains a volatile "t" leaf, so batching is
  // DECLINED -- the spy scope-guard (invoked only once a node passes every gate
  // and yields >1 batch) is never called, yet the standard scheme still gives
  // the correct result.
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
        sequant::accept_any_index{}, spy, is_volatile_t));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE_FALSE(batched);  // volatile subtree -> not batched
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
  // and the inner virtual a4. Scoping the batch axis to the occupied space
  // selects i3, so the batched partials slice ToT leaves over i3
  // (slice_array_over_mode) and sum ToT partials (add_inplace) -- the two
  // annotation-free ToT array operations that must emit an "outer;inner"
  // annotation rather than a flat one (else DistArray's is_tot_index() trips).
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"I{a4<i2,i3>,a1<i1,i2>;i1,i2} * s{a2<i1,i2>;a4<i2,i3>}");
  std::string const target = "i_2,i_1;a_1i_1i_2,a_2i_1i_2";
  auto const node = eval_node(expr);

  // Reference first (non-batched), so yield's random leaf arrays are generated
  // and cached; the batched evaluator reuses the same arrays.
  auto const ref = evaluate(node, target, yield)->get<ArrayToT>();

  auto const occ =
      sequant::get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

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
  // siblings, both carrying the auxiliary batch axis x_1 (free at the
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
  auto const x1 = [&] {  // the contracted (batch) axis
    auto axes = sequant::contracted_indices(node);
    REQUIRE(axes.size() == 1);
    return axes[0];
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
    // axis (an index the shared subnode does not carry) gives it signature
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
    // M1 (axis x_1) = X * D2 with X = (D * u): D's two occurrences (inside X
    // and as the root's sibling D2) are visited with the same signature, so D
    // alone would be registered. M2 (bogus axis) re-encounters X with
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

TEST_CASE("eval_batched_custom_evaluator dedups within-batch repeats",
          "[eval]") {
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
  auto const node = eval_node(expr);
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

TEST_CASE("eval_batched_custom_evaluator group replay", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12, 12};
  yield_.set_max_tile(4);

  // Two persistent finals sharing the sub-intermediate S = g*h (carries the
  // auxiliary batch axis x_1):
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
  auto const n1 = eval_node(t1);
  auto const n2 = eval_node(t2);

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
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      sequant::accept_any_index{}, sequant::make_no_scope_guard{},
      is_volatile_t));

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

TEST_CASE("eval_batched_custom_evaluator group replay layers nested finals",
          "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TA::TArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, 4, 12, 12};
  yield_.set_max_tile(4);

  // F_in = (g*h)*p (contracts the aux index x_1; persistent) nests inside
  // F_out, which contracts its own aux axis x_2:  F_out = (F_in * r) * q,
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
  auto const n_out = eval_node(t_out);
  auto const n_in = eval_node(t_in);

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
