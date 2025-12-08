#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse.hpp>

#include <btas/btas.h>
#include <btas/tensor_func.h>
#include <boost/regex.hpp>

#include <string>
#include <vector>

namespace {

auto eval_node(sequant::ExprPtr const& expr) {
  using namespace sequant;
  return binarize<EvalExprBTAS>(expr);
}

static auto const idx_rgx = boost::wregex{L"([ia])([↑↓])?_?(\\d+)"};
auto tensor_to_key(sequant::Tensor const& tnsr) {
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i" ? L"o" : L"v") + mo[2].str();
  };

  auto const tnsr_deparsed = sequant::deparse(tnsr.clone(), false);
  return boost::regex_replace(tnsr_deparsed, idx_rgx, formatter);
}

[[maybe_unused]] auto tensor_to_key(std::wstring_view spec) {
  return tensor_to_key(sequant::parse_expr(spec, sequant::Symmetry::Nonsymm)
                           ->as<sequant::Tensor>());
}

template <typename Tensor_t>
class rand_tensor_yield {
 private:
  size_t const nocc_;
  size_t const nvirt_;
  mutable std::map<std::wstring, sequant::ResultPtr> label_to_tnsr_;

 public:
  rand_tensor_yield(size_t noccupied, size_t nvirtual)
      : nocc_{noccupied}, nvirt_{nvirtual} {}

  [[nodiscard]] Tensor_t make_rand_tensor(sequant::Tensor const& tnsr) const {
    using ranges::views::repeat_n;
    using ranges::views::transform;
    using sequant::IndexSpace;
    auto isr = sequant::get_default_context().index_space_registry();

    SEQUANT_ASSERT(
        ranges::all_of(tnsr.const_braket_indices(),
                       [&isr](auto const& idx) {
                         return idx.space() == isr->retrieve(L"i") ||
                                idx.space() == isr->retrieve(L"a");
                       }) &&
        "Unsupported IndexSpace type found while generating tensor.");

    auto rng = btas::Range{
        tnsr.const_braket_indices() | transform([this, &isr](auto const& idx) {
          return idx.space() == isr->retrieve(L"i") ? nocc_ : nvirt_;
        }) |
        ranges::to_vector};

    auto result = Tensor_t{rng};
    result.generate(
        []() { return static_cast<double>(std::rand()) / RAND_MAX; });
    return result;
  }

  sequant::ResultPtr operator()(sequant::Tensor const& tnsr) const {
    using namespace sequant;
    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_tnsr_.find(label);
        found != label_to_tnsr_.end()) {
      //      std::wcout << "label = [" << label << "] FOUND in cache.
      //      Returning.."
      //                 << std::endl;
      return found->second;
    }
    auto t = make_rand_tensor(tnsr);
    auto&& success = label_to_tnsr_.emplace(
        label, eval_result<ResultTensorBTAS<Tensor_t>>(std::move(t)));
    SEQUANT_ASSERT(success.second && "couldn't store tensor!");
    //    std::wcout << "label = [" << label << "] NotFound in cache.
    //    Creating.."
    //               << std::endl;
    return success.first->second;
  }

  sequant::ResultPtr operator()(
      sequant::meta::can_evaluate auto const& node) const {
    using namespace sequant;
    if (node->result_type() == sequant::ResultType::Tensor) {
      SEQUANT_ASSERT(node->expr()->template is<sequant::Tensor>());
      return (*this)(node->expr()->template as<sequant::Tensor>());
    }

    using result_t = ResultScalar<double>;

    SEQUANT_ASSERT(node->expr()->template is<sequant::Constant>());
    auto d = node->as_constant().template value<double>();
    return eval_result<result_t>(d);
  }

  ///
  /// \param label eg. t{v,v;o,o}, f{o;v}
  /// \return const ref to Tensor_t type tensor
  /// \note The tensor should be already present in the yielder cache
  ///       otherwise throws assertion error. To avoid that use the other
  ///       overload of operator() that takes sequant::Tensor const&
  sequant::ResultPtr operator()(std::wstring_view label) const {
    auto&& found = label_to_tnsr_.find(label.data());
    if (found == label_to_tnsr_.end()) {
      SEQUANT_ASSERT(false && "attempted access of non-existent tensor!");
    }
    return found->second;
  }
};

using namespace sequant;

template <
    typename Iterable,
    std::enable_if_t<
        std::is_convertible_v<sequant::meta::range_value_t<Iterable>, Index> &&
            !sequant::meta::is_statically_castable_v<
                Iterable const&, std::wstring>  // prefer the ctor taking the
                                                // std::wstring
        ,
        bool> = true>
container::svector<long> tidxs(Iterable const& indices) noexcept {
  return sequant::EvalExprBTAS::index_hash(indices) |
         ranges::to<container::svector<long>>;
}

container::svector<long> tidxs(Tensor const& tnsr) noexcept {
  return sequant::EvalExprBTAS::index_hash(tnsr.const_braket_indices()) |
         ranges::to<container::svector<long>>;
}

container::svector<long> tidxs(
    ExprPtr expr, std::initializer_list<size_t> tnsr_coords) noexcept {
  auto tnsr_p = expr;
  for (auto i : tnsr_coords) tnsr_p = tnsr_p->at(i);
  SEQUANT_ASSERT(tnsr_p->is<Tensor>());
  return tidxs(tnsr_p->as<Tensor>());
}

auto terse_index = [](std::wstring const& spec) {
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return mo[1].str() + mo[2].str() + L"_" + mo[3].str();
  };
  return boost::regex_replace(spec, idx_rgx, formatter);
};

container::svector<long> tidxs(std::wstring const& csv) noexcept {
  using ranges::views::all;
  using ranges::views::split;
  using ranges::views::transform;
  auto const detersed = terse_index(csv);
  return tidxs(detersed | split(L',') |
               transform([](auto&& v) { return ranges::to<std::wstring>(v); }));
}

}  // namespace

TEST_CASE("eval_with_btas", "[eval_btas]") {
  using ranges::views::transform;
  using namespace sequant;
  using namespace sequant;

  using BTensorD = btas::Tensor<double>;

  auto norm = [](BTensorD const& tnsr) { return btas::norm(tnsr); };

  std::srand(2023);
  const size_t nocc = 2, nvirt = 20;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};
  auto yield = [&yield_](std::wstring_view lbl) -> BTensorD const& {
    return yield_(lbl)->get<BTensorD>();
  };

  auto eval = [&yield_](sequant::ExprPtr const& expr,
                        container::svector<long> const& target_labels) {
    return evaluate(eval_node(expr), target_labels, yield_)->get<BTensorD>();
  };

  auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                             container::svector<long> const& target_labels) {
    auto cache = sequant::CacheManager::empty();
    return evaluate_symm(eval_node(expr), target_labels, yield_)
        ->get<BTensorD>();
  };

  auto eval_antisymm = [&yield_](
                           sequant::ExprPtr const& expr,
                           container::svector<long> const& target_labels) {
    return evaluate_antisymm(eval_node(expr), target_labels, yield_)
        ->get<BTensorD>();
  };

  auto eval_biorthogonal_nns_project =
      [&yield_](sequant::ExprPtr const& expr,
                container::svector<long> const& target_labels) {
        return evaluate_biorthogonal_nns_project(eval_node(expr), target_labels,
                                                 yield_)
            ->get<BTensorD>();
      };

  auto parse_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, sequant::Symmetry::Antisymm);
  };

  SECTION("Summation") {
    auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");
    auto const tidx1 = tidxs(expr1, {0});
    auto eval1 = eval(expr1, tidx1);

    auto man1 = yield(L"t{v;o}");
    man1 += BTensorD{btas::permute(yield(L"f{o;v}"), {1, 0})};

    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));

    auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");
    auto const tidx2 = tidxs(expr2, {0, 0});
    auto eval2 = eval(expr2, tidx2);

    auto man2 = yield(L"t{v;o}");
    btas::scal(2.0, man2);
    auto temp = BTensorD{btas::permute(yield(L"f{o;v}"), {1, 0})};
    btas::scal(1.5, temp);
    man2 += temp;

    REQUIRE(norm(eval2) == Catch::Approx(norm(man2)));
  }

  SECTION("Product") {
    auto expr1 =
        parse_antisymm(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{ i1, i2}");
    auto const tidx1 = tidxs(L"i1,i4,a1,a4");
    auto eval1 = eval(expr1, tidx1);

    // mnemonics
    // ===
    // i looks like 1
    // a looks like 7

    BTensorD man1;
    auto const& g = yield(L"g{o,o;v,v}");
    auto const& t2 = yield(L"t{v,v;o,o}");
    btas::contract(0.5, g, {12, 14, 72, 74}, t2, {71, 72, 11, 12}, 0.0, man1,
                   {11, 14, 71, 74});
    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));

    auto expr2 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a1,a3}^{i3,i4} * t_{a2,a4}^{ i1, i2}");
    auto tidx2 = tidxs(L"i1,i2,a1,a2");
    auto eval2 = eval(expr2, tidx2);

    BTensorD man2, temp;
    btas::contract(1.0, g, {13, 14, 73, 74}, t2, {71, 73, 13, 14}, 0.0, temp,
                   {71, 74});
    btas::contract(-0.25, temp, {71, 74}, t2, {72, 74, 11, 12}, 0.0, man2,
                   {11, 12, 71, 72});
    REQUIRE(norm(eval2) == Catch::Approx(norm(man2)));
  }

  SECTION("Summation and Product") {
    auto expr1 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
    auto tidx1 = tidxs(L"i1,i2,a1,a2");
    auto eval1 = eval(expr1, tidx1);

    auto const& g = yield(L"g{o,o;v,v}");
    auto const& t2 = yield(L"t{v,v;o,o}");
    BTensorD temp1, man1;
    btas::contract(1.0, g, {13, 14, 73, 74}, t2, {71, 73, 13, 14}, 0.0, temp1,
                   {71, 74});
    btas::contract((-1 / 4.0), temp1, {71, 74}, t2, {72, 74, 11, 12}, 0.0, man1,
                   {11, 12, 71, 72});
    temp1.clear();
    btas::contract(1.0, g, {13, 14, 73, 74}, t2, {73, 74, 11, 12}, 0.0, temp1,
                   {11, 12, 13, 14});
    BTensorD temp2;
    btas::contract((1 / 16.0), temp1, {11, 12, 13, 14}, t2, {71, 72, 13, 14},
                   0.0, temp2, {11, 12, 71, 72});
    man1 += temp2;
    temp1.clear();
    temp2.clear();
    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));
  }

  SECTION("Antisymmetrization") {
    using btas::permute;

    auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto tidx1 = tidxs(L"i_1,i_2,a_1,a_2");
    auto eval1 = eval_antisymm(expr1, tidx1);

    auto const& g = yield(L"g{o,o;v,v}");
    BTensorD man1{g.range()}, temp{g.range()};
    man1.fill(0);
    temp.fill(0);

    man1 += BTensorD{permute(g, {0, 1, 2, 3})};

    temp += BTensorD{permute(g, {1, 0, 2, 3})};
    btas::scal(-1.0, temp);
    man1 += temp;

    temp.fill(0);
    temp += BTensorD{permute(g, {0, 1, 3, 2})};
    btas::scal(-1.0, temp);
    man1 += temp;

    temp.clear();
    man1 += BTensorD{permute(g, {1, 0, 3, 2})};

    btas::scal(0.5, man1);

    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));

    auto expr2 = parse_antisymm(L"R_{a1,a2}^{i1}");
    auto tidx2 = tidxs(L"a_1,a_2,i_1");
    auto eval2 = eval_antisymm(expr2, tidx2);

    auto const& r = yield(L"R{v,v;o}");
    BTensorD man2{r.range()}, temp2{r.range()};
    man2.fill(0.0);
    temp2.fill(0.0);

    man2 += BTensorD{permute(r, {0, 1, 2})};

    temp2 += BTensorD{permute(r, {1, 0, 2})};
    btas::scal(-1.0, temp2);
    man2 += temp2;
    temp2.clear();

    REQUIRE(norm(eval2) == Catch::Approx(norm(man2)));
  }

  SECTION("Symmetrization") {
    using btas::permute;

    auto expr1 = parse_antisymm(L"1/2 * g_{i1, i2}^{a1, a2}");
    auto tidx1 = tidxs(L"i_1,i_2,a_1,a_2");
    auto eval1 = eval_symm(expr1, tidx1);

    auto const& g = yield(L"g{o,o;v,v}");

    BTensorD man1{g.range()};
    man1.fill(0);

    man1 += BTensorD{permute(g, {0, 1, 2, 3})};
    man1 += BTensorD{permute(g, {1, 0, 3, 2})};
    btas::scal(0.5, man1);

    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));
  }

  SECTION("Biorthogonal Cleanup") {
    using btas::permute;
    // low-rank residuals: skip cleanup
    auto expr1 = parse_antisymm(L"R_{a1, a2}^{i1, i2}");
    auto tidx1 = tidxs(L"a_1,a_2,i_1,i_2");
    auto eval1 = eval_biorthogonal_nns_project(expr1, tidx1);
    auto const& r1 =
        yield(L"R{v,v;o,o}");  // Assuming v = virtual, o = occupied

    BTensorD man1{r1.range()};
    man1.fill(0);
    man1 = r1;

    REQUIRE(norm(eval1) == Catch::Approx(norm(man1)));

    BTensorD zero1{r1.range()};
    zero1 = man1 - eval1;
    REQUIRE(norm(zero1) == Catch::Approx(0).margin(
                               100 * std::numeric_limits<double>::epsilon()));

    // high-rank residuals: cleanup applies:
    // result = identity - (1/ket_rank!) * sum_of_ket_permutations
    auto expr2 = parse_antisymm(L"R_{a1, a2, a3}^{i1, i2, i3}");
    auto tidx2 = tidxs(L"a_1,a_2,a_3,i_1,i_2,i_3");
    auto eval2 = eval_biorthogonal_nns_project(expr2, tidx2);
    auto const& r2 = yield(L"R{v,v,v;o,o,o}");

    BTensorD man2{r2.range()};
    man2.fill(0);
    man2 = r2;

    BTensorD perm_sum{r2.range()};
    perm_sum.fill(0);

    perm_sum += r2;
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 3, 5, 4})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 4, 3, 5})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 4, 5, 3})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 5, 3, 4})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 5, 4, 3})};

    btas::scal(1.0 / 6.0, perm_sum);
    man2 -= perm_sum;
    REQUIRE(norm(eval2) == Catch::Approx(norm(man2)));

    BTensorD zero2{r2.range()};
    zero2 = man2 - eval2;
    REQUIRE(norm(zero2) == Catch::Approx(0).margin(
                               100 * std::numeric_limits<double>::epsilon()));
  }

  SECTION("Others") {
    auto expr1 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

    auto tidx1 = tidxs(L"i1,i2,a1,a2");

    auto const eval1 =
        evaluate(eval_node(expr1), tidx1, yield_)->get<BTensorD>();

    auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
      return eval_node(x);
    }) | ranges::to_vector;

    auto const eval2 = evaluate(nodes1, tidx1, yield_)->get<BTensorD>();

    REQUIRE(norm(eval1) == Catch::Approx(norm(eval1)));
  }
}
