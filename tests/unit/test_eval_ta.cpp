#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <tiledarray.h>
#include <boost/regex.hpp>

#include <cstdlib>
#include <string>
#include <vector>

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

    for (auto&& ix : tnsr.aux()) {
      SEQUANT_ASSERT(!ix.has_proto_indices() &&
                     "Aux indices with proto indices not supported");
      outer.emplace_back(ix);
    }

    auto append_unique = [](auto& cont, auto const& el) {
      if (!ranges::contains(cont, el)) cont.emplace_back(el);
    };

    for (Index const& ix : tnsr.const_braket_indices()) {
      append_unique(ix.has_proto_indices() ? inner : outer, ix);
    }

    for (Index const& ix :
         tnsr.const_braket_indices() | transform(&Index::proto_indices) | join)
      append_unique(outer, ix);
  }

  [[nodiscard]] auto outer_inner() const noexcept {
    return ranges::views::concat(outer, inner);
  }
};

auto eval_node(sequant::ExprPtr const& expr) {
  using namespace sequant;
  auto node = binarize(expr);
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

auto tensor_to_key(sequant::Tensor const& tnsr) {
  static auto const idx_rgx = boost::wregex{L"([ia])([↑↓])?(_?\\d+)"};
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i" ? L"o" : L"v") + mo[2].str();
  };

  NestedTensorIndices oixs{tnsr};
  if (oixs.inner.empty()) {
    auto const tnsr_deparsed = sequant::deparse(tnsr.clone(), false);
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
  return tensor_to_key(sequant::parse_expr(spec, sequant::Symmetry::Nonsymm)
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
  mutable sequant::container::map<std::wstring, sequant::ResultPtr>
      label_to_er_;

 public:
  using array_type = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
  using array_tot_type =
      TA::DistArray<TA::Tensor<TA::Tensor<NumericT>>, TAPolicyT>;
  using numeric_type = NumericT;

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt} {}

  sequant::ResultPtr operator()(sequant::Variable const& var) const {
    using result_t = sequant::ResultScalar<NumericT>;

    auto make_var = []() {
      return sequant::eval_result<result_t>(
          TA::detail::MakeRandom<NumericT>::generate_value());
    };

    return label_to_er_.try_emplace(std::wstring{var.label()}, make_var())
        .first->second;
  }

  sequant::ResultPtr operator()(
      sequant::meta::can_evaluate auto const& node) const {
    using namespace sequant;
    if (node->is_tensor()) return (*this)(node->as_tensor());

    if (node->is_variable()) return (*this)(node->as_variable());

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
                              ix.space() == isr->retrieve(L"a"));
               return ix.space() == isr->retrieve(L"i") ? nocc_ : nvirt_;
             }) |
             ranges::to<container::svector<size_t>>;
    };

    NestedTensorIndices nested{tnsr};

    auto const outer_extent = make_extents(nested.outer);

    auto const outer_tr = [&outer_extent]() {
      container::vector<TA::TiledRange1> tr1s;
      tr1s.reserve(outer_extent.size());
      for (auto e : outer_extent) tr1s.emplace_back(TA::TiledRange1(0, e));
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
}  // namespace

TEST_CASE("eval_with_tiledarray", "[eval]") {
  SECTION("real") {
    using ranges::views::transform;
    using sequant::EvalExprTA;
    using sequant::evaluate;
    using sequant::evaluate_antisymm;
    using sequant::evaluate_symm;

    using TA::TArrayD;

    auto parse_antisymm = [](auto const& xpr) {
      return parse_expr(xpr, sequant::Symmetry::Antisymm);
    };

    // tnsr is assumed to be single-tiled
    auto norm = [](TArrayD const& tnsr) { return TA::norm2(tnsr); };

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
      return sequant::biorthogonal_nns_project(
          result->get<TA::TArrayD>(), eval_node(expr)->as_tensor().bra_rank());
    };

    SECTION("summation") {
      auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayD{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(norm(sum1_man) == Catch::Approx(norm(sum1_eval)));

      auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(norm(sum2_man) == Catch::Approx(norm(sum2_eval)));
    }

    SECTION("product") {
      auto expr1 =
          parse_antisymm(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayD prod1_man{};
      prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(norm(prod1_man) == Catch::Approx(norm(prod1_eval)));

      auto expr2 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayD{};
      prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      REQUIRE(norm(prod2_man) == Catch::Approx(norm(prod2_eval)));

      auto expr3 = sequant::parse_expr(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayD{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");
      REQUIRE(norm(prod3_man) == Catch::Approx(norm(prod3_eval)));

      auto expr4 = sequant::parse_expr(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayD{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(norm(prod4_man) == Catch::Approx(norm(prod4_eval)));
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

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      auto expr2 = sequant::parse_expr(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayD{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(norm(man2) == Catch::Approx(norm(eval2)));
    }

    SECTION("variable at leaves") {
      auto expr2 =
          parse_antisymm(L"(α * 2 * t_{a1}^{i1} * β) + (3/2 * f_{i1}^{a1})");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          yield_d(L"α") * 2 * yield(L"t{a1;i1}")("a1,i1") * yield_d(L"β") +
          1.5 * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(norm(sum2_man) == Catch::Approx(norm(sum2_eval)));
    }

    SECTION("Antisymmetrization") {
      auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      TArrayD zero1;
      zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");

      // https://github.com/catchorg/Catch2/issues/1444
      REQUIRE(norm(zero1) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));

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
      TArrayD zero2;
      zero2("0,1,2,3,4") = man2("0,1,2,3,4") - eval2("0,1,2,3,4");
      REQUIRE(norm(zero2) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));

      auto expr3 = parse_antisymm(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayD{};
      man3("0,1") = arr3("0,1") - arr3("1,0");

      TArrayD zero3;
      zero3("0,1") = man3("0,1") - eval3("0,1");
      REQUIRE(norm(zero3) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));
    }

    SECTION("Symmetrization") {
      auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      auto expr2 = parse_antisymm(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayD man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");

      REQUIRE(norm(man2) == Catch::Approx(norm(eval2)));
    }

    SECTION("Biorthogonal Cleanup") {
      // low-rank residuals: skip nns
      auto expr1 = parse_antisymm(L"R_{a1, a2}^{i1, i2}");
      auto eval1 = eval_biorthogonal_nns_project(expr1, "a_1,a_2,i_1,i_2");
      auto const& arr1 = yield(L"R{a1,a2;i1,i2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3");

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));
      TArrayD zero1;
      zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");
      REQUIRE(norm(zero1) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));

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

      REQUIRE(norm(man2) == Catch::Approx(norm(eval2)));
      TArrayD zero2;
      zero2("0,1,2,3,4,5") = man2("0,1,2,3,4,5") - eval2("0,1,2,3,4,5");
      REQUIRE(norm(zero1) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));

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

      REQUIRE(norm(man3) == Catch::Approx(norm(eval3)));
      TArrayD zero3;
      zero3("0,1,2,3,4,5,6,7") =
          man3("0,1,2,3,4,5,6,7") - eval3("0,1,2,3,4,5,6,7");
      REQUIRE(norm(zero3) ==
              Catch::Approx(0).margin(1000 *
                                      std::numeric_limits<double>::epsilon()));
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

      REQUIRE(norm(eval1) == Catch::Approx(norm(eval2)));
    }
  }

  SECTION("complex") {
    using TArrayC = TA::DistArray<TA::Tensor<std::complex<double>>>;
    auto norm = [](TArrayC const& tnsr) { return TA::norm2(tnsr); };

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
      auto expr1 = parse_expr(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayC{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      // todo:
      REQUIRE(norm(sum1_man) == Catch::Approx(norm(sum1_eval)));

      auto expr2 = parse_expr(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayC{};
      sum2_man("i1,a1") =
          std::complex<double>{2} * yield(L"t{a1;i1}")("a1,i1") +
          std::complex<double>{1.5} * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(norm(sum2_man) == Catch::Approx(norm(sum2_eval)));
    }

    SECTION("product") {
      auto expr1 = parse_expr(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayC prod1_man{};
      prod1_man("i4,a1,a4,i1") = std::complex<double>{1 / 2.0} *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(norm(prod1_man) == Catch::Approx(norm(prod1_eval)));

      auto expr2 = parse_expr(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayC{};
      prod2_man("a1,a2,i1,i2") = std::complex<double>{-1 / 4.0} *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      REQUIRE(norm(prod2_man) == Catch::Approx(norm(prod2_eval)));

      auto expr3 = sequant::parse_expr(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayC{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");
      REQUIRE(norm(prod3_man) == Catch::Approx(norm(prod3_eval)));

      auto expr4 = sequant::parse_expr(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayC{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(norm(prod4_man) == Catch::Approx(norm(prod4_eval)));
    }

    SECTION("sum and product") {
      auto expr1 = parse_expr(
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

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      auto expr2 = sequant::parse_expr(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayC{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(norm(man2) == Catch::Approx(norm(eval2)));
    }

    SECTION("Antisymmetrization") {
      auto expr1 = parse_expr(L"1/2 * g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = std::complex<double>{0.5} * man1("0,1,2,3");

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      TArrayC zero1;
      zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");

      // todo: Catch::Approx(0.0) == 0 fails. probably update catch2 version
      // REQUIRE(Approx(norm(zero1)) == 0);

      // odd-ranked tensor
      auto expr2 = parse_expr(L"g_{i1, i2, i3}^{a1, a2}");
      auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2}");

      auto man2 = TArrayC{};
      man2("0,1,2,3,4") =
          arr2("0,1,2,3,4") - arr2("1,0,2,3,4") + arr2("1,2,0,3,4") -
          arr2("2,1,0,3,4") + arr2("2,0,1,3,4") - arr2("0,2,1,3,4") -
          arr2("0,1,2,4,3") + arr2("1,0,2,4,3") - arr2("1,2,0,4,3") +
          arr2("2,1,0,4,3") - arr2("2,0,1,4,3") + arr2("0,2,1,4,3");
      TArrayC zero2;
      zero2("0,1,2,3,4") = man2("0,1,2,3,4") - eval2("0,1,2,3,4");
      REQUIRE(norm(zero2) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));

      auto expr3 = parse_expr(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayC{};
      man3("0,1") = arr3("0,1") - arr3("1,0");

      TArrayC zero3;
      zero3("0,1") = man3("0,1") - eval3("0,1");
      REQUIRE(norm(zero3) == Catch::Approx(0).margin(
                                 100 * std::numeric_limits<double>::epsilon()));
    }

    SECTION("Symmetrization") {
      auto expr1 = parse_expr(L"1/2 * g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(norm(man1) == Catch::Approx(norm(eval1)));

      auto expr2 = parse_expr(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayC man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");

      REQUIRE(norm(man2) == Catch::Approx(norm(eval2)));
    }

    SECTION("Others") {
      using namespace std::string_literals;
      auto expr1 = parse_expr(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

      auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                       ->get<TArrayC>();

      auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
        return eval_node(x);
      }) | ranges::to_vector;

      auto eval2 = evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TArrayC>();

      REQUIRE(norm(eval1) == Catch::Approx(norm(eval2)));
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
      auto const node = eval_node(parse_expr(expr_str));
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

      auto const node = eval_node(parse_expr(expr_str));
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
      auto const node = eval_node(parse_expr(expr_str));

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
  }
}
