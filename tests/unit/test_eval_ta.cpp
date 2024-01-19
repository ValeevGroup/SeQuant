#include "catch.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval.hpp>

#include <tiledarray.h>
#include <boost/regex.hpp>

#include <cstdlib>
#include <string>
#include <vector>

namespace {
auto eval_node(sequant::ExprPtr const& expr) {
  return sequant::eval_node<sequant::EvalExprTA>(expr);
}

auto tensor_to_key(sequant::Tensor const& tnsr) {
  static auto const idx_rgx = boost::wregex{L"([ia])([↑↓])?(_?\\d+)"};
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i" ? L"o" : L"v") + mo[2].str();
  };

  auto const tnsr_deparsed = sequant::deparse_expr(tnsr.clone(), false);
  return boost::regex_replace(tnsr_deparsed, idx_rgx, formatter);
}

auto tensor_to_key(std::wstring_view spec) {
  return tensor_to_key(sequant::parse_expr(spec, sequant::Symmetry::nonsymm)
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
  mutable sequant::container::map<std::wstring, sequant::ERPtr> label_to_er_;

 public:
  using array_type = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
  using array_tot_type =
      TA::DistArray<TA::Tensor<TA::Tensor<NumericT>>, TAPolicyT>;
  using numeric_type = NumericT;

  rand_tensor_yield(TA::World& world_, size_t nocc, size_t nvirt)
      : world{world_}, nocc_{nocc}, nvirt_{nvirt} {}

  sequant::ERPtr operator()(sequant::Variable const& var) const {
    using result_t = sequant::EvalScalar<NumericT>;

    auto make_var = []() {
      return sequant::eval_result<result_t>(
          TA::detail::MakeRandom<NumericT>::generate_value());
    };

    return label_to_er_.try_emplace(std::wstring{var.label()}, make_var())
        .first->second;
  }

  template <typename T, typename = std::enable_if_t<sequant::IsEvaluable<T>>>
  sequant::ERPtr operator()(T const& node) const {
    using namespace sequant;
    if (node->is_tensor()) return (*this)(node->as_tensor());

    if (node->is_variable()) return (*this)(node->as_variable());

    assert(node->is_constant());

    using result_t = EvalScalar<NumericT>;

    auto d = (node->as_constant()).template value<NumericT>();
    return eval_result<result_t>(d);
  }

  sequant::ERPtr operator()(sequant::Tensor const& tnsr) const {
    using namespace ranges::views;
    using namespace sequant;

    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_er_.find(label); found != label_to_er_.end()) {
      //      std::cout << "label = [" << sequant::to_string(label)
      //                << "] FOUND in cache. Returning.." << std::endl;
      return found->second;
    }

    ERPtr result{nullptr};

    auto make_extents = [this](auto&& ixs) -> container::svector<size_t> {
      return ixs | transform([this](auto const& ix) -> size_t {
               assert(ix.space() == IndexSpace::active_occupied ||
                      ix.space() == IndexSpace::active_unoccupied);
               return ix.space() == IndexSpace::active_occupied ? nocc_
                                                                : nvirt_;
             }) |
             ranges::to<container::svector<size_t>>;
    };

    NestedTensorIndices nested{tnsr};

    auto const outer_extent = make_extents(nested.outer);
    auto const outer_tr = TA::TiledRange{
        outer_extent | transform([](auto e) { return TA::TiledRange1(0, e); }) |
        ranges::to_vector};
    auto const outer_r = TA::Range(outer_extent);

    if (nested.inner.empty()) {
      // regular tensor
      using ArrayT = TA::DistArray<TA::Tensor<NumericT>, TAPolicyT>;
      ArrayT array{world, outer_tr};
      for (auto it = array.begin(); it != array.end(); ++it)
        if (array.is_local(it.index()))
          *it = world.taskq.add(random_tensor<NumericT>, it.make_range());
      result = eval_result<EvalTensorTA<ArrayT>>(array);
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

      result = eval_result<EvalTensorOfTensorTA<ArrayT>>(array);
    }

    auto success = label_to_er_.emplace(label, result);
    assert(success.second && "couldn't store ERPtr!");
    //    std::cout << "label = [" << sequant::to_string(label)
    //              << "] NotFound in cache. Creating.." << std::endl;
    assert(success.first->second);
    return success.first->second;
  }

  ///
  /// \param label eg.
  ///  - 't_vvoo', 'f_ov' for generic tensor key strings
  ///  - 't{a1,a2;i1,i2}', 'f{i1;a1}' supported by sequant::parse_expr
  ///
  /// \return ERPtr
  ///
  /// \note The ERPtr should already exist in the cache otherwise throws.
  ///       This overload is only intended to access already existing ERPtrs
  ///       from the cache. To create a new cache entry use the
  ///       operator()(Tesnor const&) overload.
  ///
  sequant::ERPtr operator()(std::wstring_view label) const {
    auto&& found = label_to_er_.find(label.data());
    if (found == label_to_er_.end())
      found = label_to_er_.find(tensor_to_key(label));
    if (found == label_to_er_.end())
      throw std::runtime_error{"attempted access of non-existent ERPtr!"};
    return found->second;
  }
};

}  // namespace

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using sequant::EvalExprTA;
  using sequant::evaluate;
  using sequant::evaluate_antisymm;
  using sequant::evaluate_symm;

  using TA::TArrayD;

  auto parse_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, sequant::Symmetry::antisymm);
  };

  // tnsr is assumed to be single-tiled
  auto norm = [](TArrayD const& tnsr) { return TA::norm2(tnsr); };

  auto& world = TA::get_default_world();

  const size_t nocc = 2, nvirt = 20;
  auto yield_ = rand_tensor_yield<double, TA::DensePolicy>{world, nocc, nvirt};
  auto yield = [&yield_](std::wstring_view lbl) -> TA::TArrayD const& {
    return yield_(lbl)->get<TA::TArrayD>();
  };

  auto yield_d = [&yield_](std::wstring_view lbl) ->
      typename TA::TArrayD::numeric_type {
        return yield_(lbl)->get<typename TA::TArrayD::numeric_type>();
      };

  auto eval = [&yield_](sequant::ExprPtr const& expr,
                        std::string const& target_labels) {
    return evaluate(eval_node(expr), target_labels, yield_)->get<TA::TArrayD>();
  };

  auto eval_symm =
      [&yield_](sequant::ExprPtr const& expr, std::string const& target_labels,
                sequant::container::svector<std::array<size_t, 3>> const&
                    groups = {}) {
        return evaluate_symm(eval_node(expr), target_labels, groups, yield_)
            ->get<TA::TArrayD>();
      };

  auto eval_antisymm =
      [&yield_](sequant::ExprPtr const& expr, std::string const& target_labels,
                sequant::container::svector<std::array<size_t, 3>> const&
                    groups = {}) {
        return evaluate_antisymm(eval_node(expr), target_labels, groups, yield_)
            ->get<TA::TArrayD>();
      };

  SECTION("summation") {
    auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

    auto sum1_eval = eval(expr1, "i_1,a_1");

    auto sum1_man = TArrayD{};
    sum1_man("i1,a1") =
        yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

    auto sum2_eval = eval(expr2, "i_1,a_1");

    auto sum2_man = TArrayD{};
    sum2_man("i1,a1") =
        2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 = parse_antisymm(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                               yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                               yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, i4}");
    auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                               yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                               yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                               yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
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

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("variable at leaves") {
    auto expr2 =
        parse_antisymm(L"(α * 2 * t_{a1}^{i1} * β) + (3/2 * f_{i1}^{a1})");

    auto sum2_eval = eval(expr2, "i_1,a_1");

    auto sum2_man = TArrayD{};
    sum2_man("i1,a1") =
        yield_d(L"α") * 2 * yield(L"t{a1;i1}")("a1,i1") * yield_d(L"β") +
        1.5 * yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("Antisymmetrization") {
    auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
    auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

    auto man1 = TArrayD{};
    man1("0,1,2,3") =
        arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));

    TArrayD zero1;
    zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");

    // todo: Approx(0.0) == 0 fails. probably update catch2 version
    // REQUIRE(Approx(norm(zero1)) == 0);

    // partial antisymmetrization

    // g("0,1,2,3,4,5")
    // suppose [3, 4,]
    //         [0, 1,]
    // form particle antisymmetric pair of index ranges
    auto expr2 = parse_antisymm(L"g_{i1,i2,i3}^{a1,a2,a3}");
    auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3", {{0, 3, 2}});
    auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
    auto man2 = TArrayD{};
    man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") - arr2("0,1,2,4,3,5") +
                          arr2("1,0,2,4,3,5") - arr2("1,0,2,3,4,5");

    REQUIRE(norm(man2) == Approx(norm(eval2)));

    TArrayD zero2;
    zero2("0,1,2,3,4,5") = man2("0,1,2,3,4,5") - eval2("0,1,2,3,4,5");
    // todo: might fail. update catch2 version
    //    REQUIRE(Approx(norm(zero2)) == 0);
  }

  SECTION("Symmetrization") {
    auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
    auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

    auto man1 = TArrayD{};
    man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));

    TArrayD zero1;
    zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");
    // todo: update catch2
    // REQUIRE(Approx(norm(zero1)) == 0);

    // partial symmetrization

    // g("0,1,2,3,4,5,6,7")
    //    |-| |-| |-| |-|
    //     ^-------^
    //         ^-------^
    // suppose [0,1] and [4,5] form a particle symmetric pair of index ranges
    // and, [2,3] and [6,7] form another particle symmetric pair of index ranges

    auto expr2 = parse_antisymm(L"g_{i1,i2,i3,i4}^{a1,a2,a3,a4}");
    //                                0  1  2  3    4  5  6  7

    auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,i_4,a_1,a_2,a_3,a_4",
                           {{0, 4, 2}, {2, 6, 2}});
    auto const& arr2 = yield(L"g{i1,i2,i3,i4;a1,a2,a3,a4}");
    TArrayD man2;
    man2("i,j,k,l,a,b,c,d") = arr2("i,j,k,l,a,b,c,d") +
                              arr2("i,j,l,k,a,b,d,c") + arr2("j,i,k,l,b,a,c,d");

    REQUIRE(norm(man2) == Approx(norm(eval2)));

    TArrayD zero2;
    zero2("i,j,k,l,a,b,c,d") =
        man2("i,j,k,l,a,b,c,d") - eval2("i,j,k,l,a,b,c,d");
    // todo: update catch2
    // REQUIRE(Approx(norm(zero2)) == 0);
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

    REQUIRE(norm(eval1) == Approx(norm(eval2)));
  }
}

TEST_CASE("TEST_EVAL_USING_TA_COMPLEX", "[eval]") {
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

  auto eval_symm =
      [&yield_](sequant::ExprPtr const& expr, std::string const& target_labels,
                sequant::container::svector<std::array<size_t, 3>> const&
                    groups = {}) {
        return evaluate_symm(eval_node(expr), target_labels, groups, yield_)
            ->get<TArrayC>();
      };

  auto eval_antisymm =
      [&yield_](sequant::ExprPtr const& expr, std::string const& target_labels,
                sequant::container::svector<std::array<size_t, 3>> const&
                    groups = {}) {
        return evaluate_antisymm(eval_node(expr), target_labels, groups, yield_)
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
    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_expr(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

    auto sum2_eval = eval(expr2, "i_1,a_1");

    auto sum2_man = TArrayC{};
    sum2_man("i1,a1") = std::complex<double>{2} * yield(L"t{a1;i1}")("a1,i1") +
                        std::complex<double>{1.5} * yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 = parse_expr(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

    TArrayC prod1_man{};
    prod1_man("i4,a1,a4,i1") = std::complex<double>{1 / 2.0} *
                               yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                               yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_expr(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, i4}");
    auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

    auto prod2_man = TArrayC{};
    prod2_man("a1,a2,i1,i2") = std::complex<double>{-1 / 4.0} *
                               yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                               yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                               yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
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

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Antisymmetrization") {
    auto expr1 = parse_expr(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
    auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

    auto man1 = TArrayC{};
    man1("0,1,2,3") =
        arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

    man1("0,1,2,3") = std::complex<double>{0.5} * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));

    TArrayC zero1;
    zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");

    // todo: Approx(0.0) == 0 fails. probably update catch2 version
    // REQUIRE(Approx(norm(zero1)) == 0);

    // partial antisymmetrization

    // g("0,1,2,3,4,5")
    // suppose [3, 4,]
    //         [0, 1,]
    // form particle antisymmetric pair of index ranges
    auto expr2 = parse_expr(L"g_{i1,i2,i3}^{a1,a2,a3}");
    auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3", {{0, 3, 2}});
    auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
    auto man2 = TArrayC{};
    man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") - arr2("0,1,2,4,3,5") +
                          arr2("1,0,2,4,3,5") - arr2("1,0,2,3,4,5");

    REQUIRE(norm(man2) == Approx(norm(eval2)));

    TArrayC zero2;
    zero2("0,1,2,3,4,5") = man2("0,1,2,3,4,5") - eval2("0,1,2,3,4,5");
    // todo: might fail. update catch2 version
    //    REQUIRE(Approx(norm(zero2)) == 0);
  }

  SECTION("Symmetrization") {
    auto expr1 = parse_expr(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
    auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

    auto man1 = TArrayC{};
    man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));

    TArrayC zero1;
    zero1("0,1,2,3") = man1("0,1,2,3") - eval1("0,1,2,3");
    // todo: update catch2
    // REQUIRE(Approx(norm(zero1)) == 0);

    // partial symmetrization

    // g("0,1,2,3,4,5,6,7")
    //    |-| |-| |-| |-|
    //     ^-------^
    //         ^-------^
    // suppose [0,1] and [4,5] form a particle symmetric pair of index ranges
    // and, [2,3] and [6,7] form another particle symmetric pair of index ranges

    auto expr2 = parse_expr(L"g_{i1,i2,i3,i4}^{a1,a2,a3,a4}");
    //                                0  1  2  3    4  5  6  7

    auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,i_4,a_1,a_2,a_3,a_4",
                           {{0, 4, 2}, {2, 6, 2}});
    auto const& arr2 = yield(L"g{i1,i2,i3,i4;a1,a2,a3,a4}");
    TArrayC man2;
    man2("i,j,k,l,a,b,c,d") = arr2("i,j,k,l,a,b,c,d") +
                              arr2("i,j,l,k,a,b,d,c") + arr2("j,i,k,l,b,a,c,d");

    REQUIRE(norm(man2) == Approx(norm(eval2)));

    TArrayC zero2;
    zero2("i,j,k,l,a,b,c,d") =
        man2("i,j,k,l,a,b,c,d") - eval2("i,j,k,l,a,b,c,d");
    // todo: update catch2
    // REQUIRE(Approx(norm(zero2)) == 0);

    auto expr3 = parse_expr(L"g_{i1,i2,i3}^{a1,a2,a3}");
    auto eval3 = eval_symm(expr3, "i_1,i_2,i_3,a_1,a_2,a_3", {{0, 3, 3}});
    auto const& arr3 = yield(L"g{i1,i2,i3;a1,a2,a3}");
    TArrayC man3;
    man3("i,j,k,a,b,c") = arr3("i,j,k,a,b,c") + arr3("i,k,j,a,c,b") +
                          arr3("j,i,k,b,a,c") + arr3("j,k,i,b,c,a") +
                          arr3("k,i,j,c,a,b") + arr3("k,j,i,c,b,a");
    REQUIRE(norm(man3) == Approx(norm(eval3)));
    TArrayC zero3;
    zero3("i,j,k,a,b,c") = eval3("i,j,k,a,b,c") - man3("i,j,k,a,b,c");
  }

  SECTION("Others") {
    using namespace std::string_literals;
    auto expr1 = parse_expr(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

    auto eval1 =
        evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)->get<TArrayC>();

    auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
      return eval_node(x);
    }) | ranges::to_vector;

    auto eval2 = evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TArrayC>();

    REQUIRE(norm(eval1) == Approx(norm(eval2)));
  }
}

TEST_CASE("TEST_EVAL_USING_TA_TOT", "[eval_tot]") {
  using namespace sequant;

  //
  // eg: approx_equal("i,j;a,b", arr1, arr2)
  // - arr1 and arr2 are DistArrays with equal TiledRange and matching Range for
  //   the inner tensors at the corresponding tile positions.
  // - 'i', 'j', 'a', and 'b' are dummy indices that annotate the modes of outer
  //   and inner tensors. Why? Because TA::norm2 function is not supported for
  //   tensor-of-tensor tiles
  //
  auto approx_equal = [](std::string const& annot, auto const& lhs,
                         auto const& rhs) -> bool {
    return Approx(lhs(annot).dot(lhs(annot))) == rhs(annot).dot(rhs(annot));
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
      auto const& rhs = yield(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}")->get<ArrayToT>();
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
      auto const& lhs = yield(L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}")->get<ArrayToT>();
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
      auto const& lhs = yield(L"I{a1<i1,i2>,a2<i1,i2>;i1,i2}")->get<ArrayToT>();
      auto const& rhs = yield(L"g{i1,i2;a2<i1,i2>,a1<i1,i2>}")->get<ArrayToT>();
      ref = TA::dot(lhs("i_1,i_2;a_1i_1i_2,a_2i_1i_2"),
                    rhs("i_1,i_2;a_2i_1i_2,a_1i_1i_2"));
    }
    REQUIRE(result == Approx(ref));
  }
}
