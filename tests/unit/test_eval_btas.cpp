#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/btas/eval_expr.hpp>
#include <SeQuant/core/eval/backends/btas/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>  // mbpt::Spin

#include <btas/btas.h>
#include <btas/tensor_func.h>

#include <boost/regex.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/split.hpp>
#include <range/v3/view/transform.hpp>

#include <cmath>
#include <complex>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace {

auto eval_node(sequant::ExprPtr const& expr) {
  using namespace sequant;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return binarize<EvalExprBTAS>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

static auto const idx_rgx = boost::wregex{L"([ia])([↑↓])?_?(\\d+)"};
auto tensor_to_key(sequant::Tensor const& tnsr) {
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i" ? L"o" : L"v") + mo[2].str();
  };

  auto const tnsr_deparsed =
      sequant::serialize(tnsr.clone(), {.annot_symm = false});
  return boost::regex_replace(tnsr_deparsed, idx_rgx, formatter);
}

[[maybe_unused]] auto tensor_to_key(std::wstring_view spec) {
  return tensor_to_key(
      sequant::deserialize<sequant::ExprPtr>(
          spec, {.def_perm_symm = sequant::Symmetry::Nonsymm,
                 .def_braket_symm = sequant::Hermiticity::NonHermitian})
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
    using numeric_type = typename Tensor_t::numeric_type;
    result.generate([]() -> numeric_type {
      auto const re = static_cast<double>(std::rand()) / RAND_MAX;
      if constexpr (sequant::meta::is_complex_v<numeric_type>) {
        // genuinely complex data so a missing conjugation is observable
        auto const im = static_cast<double>(std::rand()) / RAND_MAX;
        return numeric_type(re, im);
      } else {
        return static_cast<numeric_type>(re);
      }
    });
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
    auto&& found = label_to_tnsr_.find(std::wstring{label});
    if (found == label_to_tnsr_.end()) {
      SEQUANT_ASSERT(false && "attempted access of non-existent tensor!");
    }
    return found->second;
  }
};

// Like rand_tensor_yield but aux-aware: lays out tensor axes in full
// bra-ket-aux order (const_braketaux_indices) and sizes indices by space,
// including the batching space z. Used to evaluate expressions that carry an
// auxiliary/batching hyperindex (Laplace-MP2 / THC style) end to end.
template <typename Tensor_t>
class aux_rand_tensor_yield {
 private:
  size_t const nocc_, nvirt_, nz_;
  mutable std::map<std::wstring, sequant::ResultPtr> cache_;

  static std::wstring key(sequant::Tensor const& t) {
    // normalize away specific index ordinals within a space so a canonicalized
    // leaf and the original deserialized tensor map to the same entry
    return boost::regex_replace(
        sequant::serialize(t.clone(), {.annot_symm = false}),
        boost::wregex{L"([iaz])[↑↓]?_?\\d+"}, L"$1");
  }

  size_t extent(sequant::Index const& idx) const {
    auto isr = sequant::get_default_context().index_space_registry();
    if (idx.space() == isr->retrieve(L"i")) return nocc_;
    if (idx.space() == isr->retrieve(L"a")) return nvirt_;
    if (idx.space() == isr->retrieve(L"z")) return nz_;
    SEQUANT_ASSERT(false && "aux_rand_tensor_yield: unsupported IndexSpace");
    return 0;
  }

 public:
  aux_rand_tensor_yield(size_t nocc, size_t nvirt, size_t nz)
      : nocc_{nocc}, nvirt_{nvirt}, nz_{nz} {}

  sequant::ResultPtr operator()(sequant::Tensor const& tnsr) const {
    using namespace sequant;
    auto const k = key(tnsr);
    if (auto it = cache_.find(k); it != cache_.end()) return it->second;
    auto rng = btas::Range{tnsr.const_braketaux_indices() |
                           ranges::views::transform(
                               [this](auto const& ix) { return extent(ix); }) |
                           ranges::to_vector};
    Tensor_t t{rng};
    t.generate([]() { return static_cast<double>(std::rand()) / RAND_MAX; });
    return cache_
        .emplace(k, eval_result<ResultTensorBTAS<Tensor_t>>(std::move(t)))
        .first->second;
  }

  sequant::ResultPtr operator()(
      sequant::meta::can_evaluate auto const& node) const {
    using namespace sequant;
    if (node->result_type() == ResultType::Tensor)
      return (*this)(node->expr()->template as<Tensor>());
    SEQUANT_ASSERT(node->expr()->template is<Constant>());
    return eval_result<ResultScalar<double>>(
        node->as_constant().template value<double>());
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

/// Leaf yielder backed by pre-built arrays, so that a symmetry relation
/// between two spellings of the same tensor (d{a;i} = -conj(d{i;a}) for an
/// anti-Hermitian d) holds exactly in the data. rand_tensor_yield draws the
/// two spellings independently, which no sign convention could then satisfy.
template <typename Tensor_t>
class pinned_tensor_yield {
 private:
  std::map<std::wstring, sequant::ResultPtr> label_to_tnsr_;

 public:
  void put(sequant::Tensor const& tnsr, Tensor_t data) {
    [[maybe_unused]] auto success = label_to_tnsr_.emplace(
        tensor_to_key(tnsr),
        sequant::eval_result<sequant::ResultTensorBTAS<Tensor_t>>(
            std::move(data)));
    SEQUANT_ASSERT(success.second && "couldn't store tensor!");
  }

  sequant::ResultPtr operator()(sequant::Tensor const& tnsr) const {
    auto found = label_to_tnsr_.find(tensor_to_key(tnsr));
    if (found == label_to_tnsr_.end())
      throw sequant::Exception("pinned_tensor_yield: no data for leaf " +
                               sequant::toUtf8(tensor_to_key(tnsr)));
    return found->second;
  }

  sequant::ResultPtr operator()(
      sequant::meta::can_evaluate auto const& node) const {
    using namespace sequant;
    if (node->result_type() == ResultType::Tensor) {
      SEQUANT_ASSERT(node->expr()->template is<Tensor>());
      return (*this)(node->expr()->template as<Tensor>());
    }
    SEQUANT_ASSERT(node->expr()->template is<Constant>());
    using numeric_type = typename Tensor_t::numeric_type;
    return eval_result<ResultScalar<numeric_type>>(
        node->as_constant().template value<numeric_type>());
  }
};

}  // namespace

TEST_CASE("eval_with_btas", "[eval_btas]") {
  using ranges::views::transform;
  using namespace sequant;
  using namespace sequant;

  using BTensorD = btas::Tensor<double>;

  auto norm = [](BTensorD const& tnsr) {
    return std::sqrt(btas::dotc(tnsr, tnsr));
  };

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
        auto result = evaluate(eval_node(expr), target_labels, yield_);
        return mbpt::biorthogonal_nns_project(
            result->get<BTensorD>(), eval_node(expr)->as_tensor().bra_rank());
      };

  auto parse_antisymm = [](auto const& xpr) {
    return deserialize<sequant::ExprPtr>(
        xpr, {.def_perm_symm = sequant::Symmetry::Antisymm,
              .def_braket_symm = sequant::Hermiticity::NonHermitian});
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

    auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
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

    btas::scal(0.25, man1);

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
    btas::scal(0.5, man2);
    temp2.clear();

    REQUIRE(norm(eval2) == Catch::Approx(norm(man2)));
  }

  SECTION("Symmetrization") {
    using btas::permute;

    auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
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

    perm_sum += BTensorD{permute(r2, {0, 1, 2, 3, 5, 4})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 4, 3, 5})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 4, 5, 3})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 5, 3, 4})};
    perm_sum += BTensorD{permute(r2, {0, 1, 2, 5, 4, 3})};

    btas::scal(1.0 / 5.0, perm_sum);
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

    REQUIRE(norm(eval1) == Catch::Approx(norm(eval2)));

    BTensorD zero2{eval2.range()};
    zero2 = eval1 - eval2;
    REQUIRE(norm(zero2) == Catch::Approx(0).margin(
                               100 * std::numeric_limits<double>::epsilon()));
  }
}

// Task 2 (multiroot-single-dag-eval): a Sum node marked accumulate_in_place()
// (Task 1's binarize()-applied mark on the left-accumulator chain of a folded
// N-ary Sum) must evaluate via Result::add_inplace() into the left operand's
// own buffer rather than the allocating Result::sum() -- EXCEPT when its left
// child is a leaf: a leaf's ResultPtr comes straight out of the caller's
// leaf_evaluator (rand_tensor_yield here memoizes by label, standing in for a
// production evaluator that reuses/caches AO integrals or amplitudes across
// calls), so this engine cannot know whether mutating it would corrupt some
// OTHER, unrelated read. This was found empirically: an earlier version of
// this fix (no leaf exclusion) newly broke the pre-existing "eval_with_btas /
// Summation" test, which reads a leaf back out of the SAME yield_ after
// evaluating a marked Sum that used it as the chain seed -- the leaf came
// back already mutated. So for the chain (((t1+t2)+t3)+t4), all 3 Sum nodes
// are STILL marked (Task 1's static property, unaffected), but only the two
// OUTER ones (whose left child is itself a Sum result, never a leaf) actually
// evaluate in place; the innermost (t1+t2), whose left child t1 is a leaf,
// falls back to the allocating sum() -- one bounded extra allocation for the
// whole chain, not per term.
//
// Verified two ways: (1) the eval trace's per-op mode tally (a real counter,
// via Logger::instance().eval.stream when eval.level > 0) confirms exactly
// which 2 of the 3 marked Sum nodes actually ran the in-place path
// (SumInplace) versus the allocating one (Sum) -- this observes the executed
// branch directly, so it does not depend on inferring anything from pointer
// values. Additionally, the leaf-exclusion safety fix itself is directly
// exercised: every leaf's buffer is byte-identical before and after the
// marked evaluation (nothing was corrupted), and the marked result's pointer
// aliases none of the 4 leaves (no accidental leaf mutation slipped through).
// (2) Numeric equality against an unmarked reference (the SAME chain with
// every mark forced off, reproducing the pre-Task-2 always-allocate path).
TEST_CASE("eval_sum_accumulate_in_place_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  auto norm = [](BTensorD const& tnsr) {
    return std::sqrt(btas::dotc(tnsr, tnsr));
  };

  // Counts non-overlapping occurrences of `needle` in `haystack` -- used to
  // tally "| Sum |" vs "| SumInplace |" trace lines. The two are
  // distinguished unambiguously: "| SumInplace |" never matches "| Sum |"
  // (the character right after "Sum" differs, 'I' vs ' '), so no collision.
  auto count_substr = [](std::string const& haystack,
                         std::string const& needle) {
    std::size_t n = 0, pos = 0;
    while ((pos = haystack.find(needle, pos)) != std::string::npos) {
      ++n;
      pos += needle.size();
    }
    return n;
  };

  std::srand(2024);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  auto sum = deserialize(L"t1{i1;a1} + t2{i1;a1} + t3{i1;a1} + t4{i1;a1}");
  REQUIRE(sum->is<Sum>());
  auto const& summands = sum->as<Sum>().summands();
  REQUIRE(summands.size() == 4);

  // Pre-populate (and capture) every leaf's buffer BEFORE either tree is
  // evaluated, so both evaluations below see the SAME leaf inputs and each
  // leaf's pristine value is known ahead of time for the safety check below.
  container::svector<ResultPtr> leaves;
  container::svector<BTensorD> leaves_pristine;
  for (auto const& s : summands) {
    leaves.push_back(yield_(s->as<Tensor>()));
    leaves_pristine.push_back(leaves.back()->get<BTensorD>());
  }

  auto marked = eval_node(sum);
  auto unmarked = eval_node(sum);

  // Sanity on the marking itself (Task 1): the whole 3-Sum chain is marked.
  std::size_t total_sum = 0, inplace = 0;
  marked.visit([&](auto const& n) {
    if (!n->is_sum()) return;
    ++total_sum;
    if (n->accumulate_in_place()) ++inplace;
  });
  REQUIRE(total_sum == 3);
  REQUIRE(inplace == 3);

  // Force the SAME chain unmarked, reproducing the pre-Task-2 always-allocate
  // sum() path, to serve as the reference. The chain is left-leaning
  // (fold_left_to_node): unmarked = ((t1+t2)+t3)+t4, so the 3 Sum nodes are
  // unmarked, unmarked.left(), and unmarked.left().left().
  unmarked->set_accumulate_in_place(false);
  unmarked.left()->set_accumulate_in_place(false);
  unmarked.left().left()->set_accumulate_in_place(false);
  REQUIRE_FALSE(unmarked->accumulate_in_place());
  REQUIRE_FALSE(unmarked.left()->accumulate_in_place());
  REQUIRE_FALSE(unmarked.left().left()->accumulate_in_place());

  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;

  // Evaluate the unmarked reference FIRST: sum() never mutates its operands,
  // so this leaves every leaf's buffer untouched for the marked run below.
  std::ostringstream trace_unmarked;
  logger.eval.level = 1;
  logger.eval.stream = &trace_unmarked;
  ResultPtr const result_unmarked = evaluate<Trace::On>(unmarked, yield_);
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(count_substr(trace_unmarked.str(), "| Sum |") == 3);
  REQUIRE(count_substr(trace_unmarked.str(), "| SumInplace |") == 0);

  // Now the marked run: 2 of the 3 marked Sum nodes (the two whose left
  // child is itself a Sum result, never a leaf) accumulate in place; the
  // innermost (t1+t2), leaf-seeded, falls back to the allocating sum().
  std::ostringstream trace_marked;
  logger.eval.level = 1;
  logger.eval.stream = &trace_marked;
  ResultPtr const result_marked = evaluate<Trace::On>(marked, yield_);
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(count_substr(trace_marked.str(), "| SumInplace |") == 2);
  REQUIRE(count_substr(trace_marked.str(), "| Sum |") == 1);

  // (1a) Safety: no leaf was mutated as a side effect -- every leaf's buffer
  // is still byte-identical to its pristine value, and the marked result
  // does not alias any of them.
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    BTensorD const diff = leaves[i]->get<BTensorD>() - leaves_pristine[i];
    REQUIRE(norm(diff) == Catch::Approx(0).margin(
                              100 * std::numeric_limits<double>::epsilon()));
    REQUIRE(result_marked.get() != leaves[i].get());
  }
  REQUIRE(result_marked.get() != result_unmarked.get());

  // (2) Numeric correctness: marked (in-place) equals unmarked (allocating).
  BTensorD zero =
      result_marked->get<BTensorD>() - result_unmarked->get<BTensorD>();
  REQUIRE(norm(zero) == Catch::Approx(0).margin(
                            100 * std::numeric_limits<double>::epsilon()));
}

TEST_CASE("eval_adjoint_complex_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorC = btas::Tensor<std::complex<double>>;

  std::srand(2023);
  const size_t nocc = 2, nvirt = 5;
  auto yield_ = rand_tensor_yield<BTensorC>{nocc, nvirt};

  // A non-Hermitian tensor's adjoint() swaps bra/ket and sets the adjointed
  // state (spelled with a trailing '⁺'); binarize lowers that to a leaf whose
  // CanonTransform conjugate-transposes the operand on retrieval. With
  // genuinely complex data the conjugation is observable (a missing conj
  // would leave imaginary parts unflipped; a pure transpose would still pass
  // a norm-only check).
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor t_adj = t;
  REQUIRE(t_adj.adjoint() == 1);
  REQUIRE(t_adj.label() == L"t");
  REQUIRE(t_adj.adjointed());

  auto node = eval_node(ex<Tensor>(t_adj));
  REQUIRE(node.leaf());
  REQUIRE(node->canon_transform().conj);
  REQUIRE(node->canon_transform().braket_swap);

  // the provider is asked for the bare array 't{a_1;i_1}': shape
  // [nvirt, nocc], indexed (a, i)
  auto const& src = yield_(node->as_tensor())->get<BTensorC>();
  REQUIRE_FALSE(node->as_tensor().adjointed());
  REQUIRE(src.extent(0) == nvirt);
  REQUIRE(src.extent(1) == nocc);

  // the transform is applied on retrieval and permutes nothing
  // (apply_canon_transform passes the node's own annot twice), so the served
  // buffer's modes keep their own labels: the requested layout is the stored
  // array's and the value is elementwise conj(src)
  auto const adj = evaluate(node, tidxs(L"a_1,i_1"), yield_)->get<BTensorC>();
  REQUIRE(adj.extent(0) == nvirt);
  REQUIRE(adj.extent(1) == nocc);

  for (size_t a = 0; a < nvirt; ++a)
    for (size_t i = 0; i < nocc; ++i) {
      auto const expected = std::conj(src(a, i));
      auto const got = adj(a, i);
      CHECK(got.real() == Catch::Approx(expected.real()).margin(1e-12));
      CHECK(got.imag() == Catch::Approx(expected.imag()).margin(1e-12));
    }
}

// A real-field odd-parity Hermitian tensor is antisymmetric under the whole
// bra<->ket exchange, p{a;i} = -p{i;a}. At the eval boundary a flat leaf keeps
// such a signed orientation as written: the leaf canonicalizer does not trade
// a bra<->ket exchange for a sign (fold_signed_braket is false), so the
// provider is asked for the spelling as written and the leaf's transform
// carries no phase. The two orientations are therefore distinct leaves, each
// with phase +1, and a sum that uses both evaluates to the denoted value, with
// and without a cache. (The orientations differ in space because the pinned
// yielder keys leaves by label and slot spaces only.)
TEST_CASE("eval_signed_leaf_phase_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  // an Index whose space carries a real field (the default field is Complex)
  auto ridx = [](std::wstring_view label) {
    Index i(label);
    IndexSpace sp = i.space();
    sp.field(Field::Real);
    return Index(label, sp);
  };
  auto p = [&ridx](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(
        L"p", bra{ridx(b)}, ket{ridx(k)},
        TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                         .conjugation_parity = ConjugationParity::Odd});
  };
  auto t = [&ridx](std::wstring_view lbl, std::wstring_view b,
                   std::wstring_view k) {
    return ex<Tensor>(lbl, bra{ridx(b)}, ket{ridx(k)});
  };
  REQUIRE(p(L"i_1", L"a_1")->as<Tensor>().braket_symmetry() ==
          BraKetSymmetry::Antisymm);

  const std::size_t nocc = 3, nvirt = 4;
  std::srand(2025);
  auto rnd = []() { return static_cast<double>(std::rand()) / RAND_MAX - 0.5; };
  // the swapped spelling reads the same array with the slots exchanged, at
  // the sign the trait asserts: p{a_1;i_1}(a, i) = -p{i_1;a_1}(i, a)
  BTensorD P{btas::Range{nocc, nvirt}};
  P.generate(rnd);
  BTensorD Pswapped{btas::Range{nvirt, nocc}};
  for (std::size_t i = 0; i < nocc; ++i)
    for (std::size_t a = 0; a < nvirt; ++a) Pswapped(a, i) = -P(i, a);
  BTensorD U{btas::Range{nvirt, nocc}};
  U.generate(rnd);
  BTensorD V{btas::Range{nvirt, nocc}};
  V.generate(rnd);

  pinned_tensor_yield<BTensorD> yield;
  yield.put(p(L"i_1", L"a_1")->as<Tensor>(), P);
  yield.put(p(L"a_1", L"i_1")->as<Tensor>(), Pswapped);
  yield.put(t(L"u", L"a_1", L"i_3")->as<Tensor>(), U);
  yield.put(t(L"v", L"a_1", L"i_3")->as<Tensor>(), V);

  // r{i_1;i_3} = p{i_1;a_1} u{a_1;i_3} + p{a_1;i_1} v{a_1;i_3}
  auto e = p(L"i_1", L"a_1") * t(L"u", L"a_1", L"i_3") +
           p(L"a_1", L"i_1") * t(L"v", L"a_1", L"i_3");
  auto node = eval_node(e);

  // the two p leaves keep their spellings: distinct hashes, no phase
  using NodeT = decltype(node);
  container::svector<NodeT const*> p_leaves;
  auto collect = [&p_leaves](auto& self, NodeT const& nd) -> void {
    if (nd.leaf()) {
      if (nd->is_tensor() && nd->as_tensor().label() == L"p")
        p_leaves.push_back(&nd);
      return;
    }
    self(self, nd.left());
    self(self, nd.right());
  };
  collect(collect, node);
  REQUIRE(p_leaves.size() == 2);
  REQUIRE((*p_leaves[0])->hash_value() != (*p_leaves[1])->hash_value());
  REQUIRE((*p_leaves[0])->canon_phase() == 1);
  REQUIRE((*p_leaves[1])->canon_phase() == 1);
  REQUIRE((*p_leaves[0])->as_tensor().bra()[0].label() == L"i_1");
  REQUIRE((*p_leaves[1])->as_tensor().bra()[0].label() == L"a_1");
  REQUIRE(node->canon_phase() == 1);

  BTensorD ref{btas::Range{nocc, nocc}};
  ref.fill(0.0);
  for (std::size_t i1 = 0; i1 < nocc; ++i1)
    for (std::size_t i3 = 0; i3 < nocc; ++i3)
      for (std::size_t a = 0; a < nvirt; ++a)
        ref(i1, i3) += P(i1, a) * U(a, i3) + Pswapped(a, i1) * V(a, i3);

  // the root is evaluated in its own canonical layout
  auto const& root_ix = node->canon_indices();
  REQUIRE(root_ix.size() == 2);
  const bool i1_first = root_ix[0].label() == L"i_1";
  auto check = [&](BTensorD const& got) {
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t i3 = 0; i3 < nocc; ++i3) {
        const double g = i1_first ? got(i1, i3) : got(i3, i1);
        CHECK(g == Catch::Approx(ref(i1, i3)).margin(1e-12));
      }
  };

  SECTION("without a cache") {
    check(evaluate(node, node->annot(), yield)->get<BTensorD>());
  }

  SECTION("through a shared cache") {
    auto cache =
        cache_manager(std::array{node}, [](auto const&) { return false; });
    check(evaluate(node, node->annot(), yield, cache)->get<BTensorD>());
  }
}

// A flat leaf whose written slot order is not the block-canonical one is
// stored under the canonical spelling and carries the canonicalization's sign
// in its CanonTransform. That sign reaches the value: a leaf's transform is
// applied once on the way out of the retrieval, so the engine hands up the
// written spelling's value, phase and all, and the two spellings share one
// provider array and one cache slot while denoting values that differ by the
// sign. (t{a_1,i_3;i_1,i_2} is Antisymm with symmetric columns, so sorting its
// bra costs one transposition.)
TEST_CASE("eval_leaf_phase_reaches_the_value_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  const std::size_t nocc = 2, nvirt = 3;
  std::srand(7);
  auto rnd = [](std::vector<std::size_t> const& extents) {
    BTensorD r{btas::Range{extents}};
    r.generate(
        []() { return static_cast<double>(std::rand()) / RAND_MAX - 0.5; });
    return r;
  };
  auto isr = get_default_context().index_space_registry();
  auto extents_of = [&isr, nocc, nvirt](Tensor const& tn) {
    std::vector<std::size_t> e;
    for (auto const& ix : tn.const_braket_indices())
      e.push_back(ix.space() == isr->retrieve(L"i") ? nocc : nvirt);
    return e;
  };

  auto t = [](std::wstring_view b0, std::wstring_view b1) {
    return ex<Tensor>(L"t", bra{Index{b0}, Index{b1}},
                      ket{Index{L"i_1"}, Index{L"i_2"}}, Symmetry::Antisymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  auto t_as_written = t(L"a_1", L"i_3");
  auto t_canon = t(L"i_3", L"a_1");
  auto v = ex<Tensor>(L"v", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"});
  auto u = ex<Tensor>(L"u", bra{L"i_1", L"i_2"}, ket{L"i_3", L"a_1"});

  auto leaf = eval_node(t_as_written);
  auto leaf_canon = eval_node(t_canon);
  REQUIRE(leaf.leaf());
  // one slot, two phases: the stored spelling is the canonical one
  REQUIRE(leaf->hash_value() == leaf_canon->hash_value());
  REQUIRE(leaf->canon_phase() == -1);
  REQUIRE(leaf_canon->canon_phase() == 1);
  REQUIRE(leaf->as_tensor().bra()[0].label() == L"i_3");
  REQUIRE(leaf->as_tensor().bra()[1].label() == L"a_1");
  REQUIRE(leaf->as_tensor().ket()[0].label() == L"i_1");
  REQUIRE(leaf->as_tensor().ket()[1].label() == L"i_2");
  REQUIRE(leaf->annot() == leaf_canon->annot());

  // T is the provider's array for the stored spelling t{i_3,a_1;i_1,i_2},
  // laid out (i_3, a_1, i_1, i_2)
  BTensorD T = rnd(extents_of(leaf->as_tensor()));
  BTensorD V = rnd({nocc, nocc, nvirt, nocc});  // v{i_1,i_2;a_1,i_3}
  BTensorD U = rnd({nocc, nocc, nocc, nvirt});  // u{i_1,i_2;i_3,a_1}
  pinned_tensor_yield<BTensorD> yield;
  yield.put(leaf->as_tensor(), T);
  yield.put(v->as<Tensor>(), V);
  yield.put(u->as<Tensor>(), U);

  auto scalar_of = [](ResultPtr const& res) {
    REQUIRE(res->is<ResultScalar<double>>());
    return res->as<ResultScalar<double>>().value();
  };
  // each network is checked with no cache and then twice through one shared
  // CacheManager, so that the second read of a slot both spellings share
  // would surface a conversion applied the wrong number of times
  auto check_scalar = [&scalar_of](auto const& node, auto const& leaf_yield,
                                   double ref) {
    auto check = [&scalar_of, ref](ResultPtr const& res) {
      CHECK(scalar_of(res) == Catch::Approx(ref).margin(1e-12));
    };
    check(evaluate(node, node->annot(), leaf_yield));
    auto cache =
        cache_manager(std::array{node}, [](auto const&) { return false; });
    check(evaluate(node, node->annot(), leaf_yield, cache));
    check(evaluate(node, node->annot(), leaf_yield, cache));
  };

  SECTION("the leaf alone is the stored array times its phase") {
    auto const got = evaluate(leaf, leaf->annot(), yield)->get<BTensorD>();
    REQUIRE(got.range() == T.range());
    for (std::size_t k = 0; k < nocc; ++k)
      for (std::size_t a = 0; a < nvirt; ++a)
        for (std::size_t i = 0; i < nocc; ++i)
          for (std::size_t j = 0; j < nocc; ++j)
            CHECK(got(k, a, i, j) ==
                  Catch::Approx(-T(k, a, i, j)).margin(1e-12));
  }

  // Σ over all four slots; the as-written spelling contributes -T
  double ref_v = 0., ref_u = 0.;
  for (std::size_t k = 0; k < nocc; ++k)
    for (std::size_t a = 0; a < nvirt; ++a)
      for (std::size_t i = 0; i < nocc; ++i)
        for (std::size_t j = 0; j < nocc; ++j) {
          ref_v += -T(k, a, i, j) * V(i, j, a, k);
          ref_u += T(k, a, i, j) * U(i, j, k, a);
        }

  SECTION("a closed product carries the leaf's phase") {
    auto node = eval_node(t_as_written * v);
    REQUIRE(node->is_scalar());
    check_scalar(node, yield, ref_v);
  }

  SECTION("both spellings in one sum, one leaf slot, two phases") {
    auto node = eval_node(t_as_written * v + t_canon * u);
    REQUIRE(node->is_scalar());
    REQUIRE(node.left().left()->hash_value() ==
            node.right().left()->hash_value());
    check_scalar(node, yield, ref_v + ref_u);
  }
}

// A product whose last combination is scalar-valued (scalar * scalar) has
// opaque children: nothing is flattened, so each child hands up the value it
// denotes, and the node's transform must carry what the slot hash leaves out
// of its children -- their phases and the conjugation a hoisted prefix
// strips. Two spellings that share the slot are read through one cache, so
// the second read surfaces a transform the node failed to record.
TEST_CASE("eval_scalar_product_node_carries_its_children_transform_btas",
          "[eval_btas]") {
  using namespace sequant;
  using C = std::complex<double>;

  SECTION("a uniformly conjugated scalar pair hoists onto the bare slot") {
    auto x = ex<Variable>(L"x"), y = ex<Variable>(L"y");
    auto xs = x->clone(), ys = y->clone();
    xs->as<Variable>().conjugate();
    ys->as<Variable>().conjugate();
    auto n_conj = eval_node(xs * ys);
    auto n_bare = eval_node(x * y);
    REQUIRE(n_conj->is_scalar());
    REQUIRE(n_conj->hash_value() == n_bare->hash_value());
    REQUIRE(n_conj->canon_transform() == CanonTransform{.conj = true});
    REQUIRE(n_bare->canon_transform().trivial());

    auto yield = [](auto const& n) -> ResultPtr {
      REQUIRE(n->is_variable());
      return eval_result<ResultScalar<C>>(
          n->as_variable().label() == L"x" ? C(1, 2) : C(3, -1));
    };
    auto scalar_of = [](ResultPtr const& res) {
      REQUIRE(res->is<ResultScalar<C>>());
      return res->as<ResultScalar<C>>().value();
    };
    auto close_to = [](C got, C ref) {
      CHECK(got.real() == Catch::Approx(ref.real()).margin(1e-12));
      CHECK(got.imag() == Catch::Approx(ref.imag()).margin(1e-12));
    };
    C const xy = C(1, 2) * C(3, -1);
    auto cache = cache_manager(std::array{n_conj, n_bare});
    close_to(scalar_of(evaluate(n_conj, n_conj->annot(), yield, cache)),
             std::conj(xy));
    // the bare spelling reads the slot the conjugated one filled
    close_to(scalar_of(evaluate(n_bare, n_bare->annot(), yield, cache)), xy);
    close_to(scalar_of(evaluate(n_bare, n_bare->annot(), yield)), xy);
  }

  SECTION("a scalar-valued product child's phase reaches the node") {
    using BTensorD = btas::Tensor<double>;
    const std::size_t nocc = 2, nvirt = 3;
    std::srand(13);
    auto rnd = [](std::vector<std::size_t> const& extents) {
      BTensorD r{btas::Range{extents}};
      r.generate(
          []() { return static_cast<double>(std::rand()) / RAND_MAX - 0.5; });
      return r;
    };
    auto t = [](std::wstring_view b0, std::wstring_view b1) {
      return ex<Tensor>(L"t", bra{Index{b0}, Index{b1}},
                        ket{Index{L"i_1"}, Index{L"i_2"}}, Symmetry::Antisymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    };
    auto v = ex<Tensor>(L"v", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"});
    auto y = ex<Variable>(L"y");
    // (t v) y, with t in the two slot orders the antisymmetry relates
    auto p_w = eval_node(ex<Product>(ExprPtrList{t(L"a_1", L"i_3"), v, y},
                                     Product::Flatten::No));
    auto p_c = eval_node(ex<Product>(ExprPtrList{t(L"i_3", L"a_1"), v, y},
                                     Product::Flatten::No));
    REQUIRE(p_w->is_scalar());
    REQUIRE(p_w.left()->is_scalar());
    // the inner products share a slot with opposite phases ...
    REQUIRE(p_w.left()->hash_value() == p_c.left()->hash_value());
    REQUIRE(p_w.left()->canon_phase() == -1);
    REQUIRE(p_c.left()->canon_phase() == 1);
    // ... and so do the outer nodes, whose transform carries that phase
    REQUIRE(p_w->hash_value() == p_c->hash_value());
    REQUIRE(p_w->canon_phase() == -1);
    REQUIRE(p_c->canon_phase() == 1);

    // T is the provider's array for the stored spelling t{i_3,a_1;i_1,i_2}
    BTensorD T = rnd({nocc, nvirt, nocc, nocc});
    BTensorD V = rnd({nocc, nocc, nvirt, nocc});  // v{i_1,i_2;a_1,i_3}
    pinned_tensor_yield<BTensorD> tensors;
    tensors.put(p_c.left().left()->as_tensor(), T);
    tensors.put(v->as<Tensor>(), V);
    auto yield = [&tensors](auto const& n) -> ResultPtr {
      if (n->is_variable()) return eval_result<ResultScalar<double>>(2.0);
      return tensors(n);
    };
    double ref_c = 0.;
    for (std::size_t k = 0; k < nocc; ++k)
      for (std::size_t a = 0; a < nvirt; ++a)
        for (std::size_t i = 0; i < nocc; ++i)
          for (std::size_t j = 0; j < nocc; ++j)
            ref_c += T(k, a, i, j) * V(i, j, a, k);
    ref_c *= 2.0;
    auto scalar_of = [](ResultPtr const& res) {
      REQUIRE(res->is<ResultScalar<double>>());
      return res->as<ResultScalar<double>>().value();
    };
    auto cache = cache_manager(std::array{p_w, p_c});
    CHECK(scalar_of(evaluate(p_w, p_w->annot(), yield, cache)) ==
          Catch::Approx(-ref_c).margin(1e-12));
    // the canonical spelling reads the slot the as-written one filled
    CHECK(scalar_of(evaluate(p_c, p_c->annot(), yield, cache)) ==
          Catch::Approx(ref_c).margin(1e-12));
    CHECK(scalar_of(evaluate(p_c, p_c->annot(), yield)) ==
          Catch::Approx(ref_c).margin(1e-12));
  }
}

// The default simplify() folds conjugate pairs of scalar-valued summands
// only, so a tensor-valued sum with a conjugate pair reaches the evaluator as
// the sum it was written as, whose adjointed summand is a network of leaves
// stored bare under the adjoint channel.
TEST_CASE("eval_tensor_valued_conjugate_pair_btas", "[eval_btas]") {
  using namespace sequant;
  using C = std::complex<double>;
  using BTensorC = btas::Tensor<C>;

  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);
  REQUIRE(Index{L"a_1"}.space().field() == Field::Complex);

  const std::size_t nocc = 2, nvirt = 3;
  std::srand(17);
  auto rnd = [](std::vector<std::size_t> const& extents) {
    BTensorC r{btas::Range{extents}};
    r.generate([]() {
      return C(static_cast<double>(std::rand()) / RAND_MAX - 0.5,
               static_cast<double>(std::rand()) / RAND_MAX - 0.5);
    });
    return r;
  };

  auto const u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"a_2"});
  auto const v = ex<Tensor>(L"v", bra{L"a_2"}, ket{L"i_1"});
  auto sum = u * v + conjugate(u * v);
  simplify(sum);
  REQUIRE(sum->is<Sum>());
  REQUIRE(sum->as<Sum>().summands().size() == 2);

  auto node = eval_node(sum);
  REQUIRE(node->op_type() == EvalOp::Sum);
  REQUIRE(node->is_tensor());

  BTensorC U = rnd({nvirt, nvirt});  // u{a_1;a_2}: U(a1, a2)
  BTensorC V = rnd({nvirt, nocc});   // v{a_2;i_1}: V(a2, i)
  pinned_tensor_yield<BTensorC> yield;
  yield.put(u->as<Tensor>(), U);
  yield.put(v->as<Tensor>(), V);

  auto const got =
      evaluate(node, tidxs(std::vector<Index>{Index{L"a_1"}, Index{L"i_1"}}),
               yield)
          ->get<BTensorC>();
  REQUIRE(got.rank() == 2);
  REQUIRE(got.extent(0) == nvirt);
  REQUIRE(got.extent(1) == nocc);
  for (std::size_t a = 0; a < nvirt; ++a)
    for (std::size_t i = 0; i < nocc; ++i) {
      C x{0., 0.};
      for (std::size_t b = 0; b < nvirt; ++b) x += U(a, b) * V(b, i);
      C const ref = x + std::conj(x);
      CHECK(got(a, i).real() == Catch::Approx(ref.real()).margin(1e-12));
      CHECK(got(a, i).imag() == Catch::Approx(0.).margin(1e-12));
    }
}

// The cell table holds a value in its canonical orientation while every
// reader wants the node's own, and the conversion is the node's whole
// CanonTransform: a node that carries a hoisted elementwise conjugation
// (t꙳ g꙳ over a real basis, whose two factors are both hoistable, so the
// intermediate carries {conj} and shares the unconjugated t g slot) must come
// back out of the table conjugated. Two roots sharing that intermediate make
// it a cell of its own rather than a transient of one production tree, which
// is what puts the store side's conversion on the path.
TEST_CASE("eval_ordered_conj_node_is_converted_whole_btas",
          "[eval_btas][ordered]") {
  using namespace sequant;
  using C = std::complex<double>;
  using BTensorC = btas::Tensor<C>;

  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  const std::size_t nocc = 2, nvirt = 3;
  std::srand(11);
  auto rnd = [](std::vector<std::size_t> const& extents) {
    BTensorC r{btas::Range{extents}};
    r.generate([]() {
      return C(static_cast<double>(std::rand()) / RAND_MAX - 0.5,
               static_cast<double>(std::rand()) / RAND_MAX - 0.5);
    });
    return r;
  };
  // an Index whose space carries a real field (the default field is Complex)
  auto ridx = [](std::wstring_view label) {
    Index i(label);
    IndexSpace sp = i.space();
    sp.field(Field::Real);
    return Index(label, sp);
  };
  // parity None keeps the '꙳' over a real basis; the default Even consumes it
  auto starred = [&ridx](std::wstring_view lbl, std::wstring_view b,
                         std::wstring_view k) {
    Tensor t(lbl, bra{ridx(b)}, ket{ridx(k)},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.kconjugate() == 1);
    REQUIRE(t.kconjugated());
    return ex<Tensor>(std::move(t));
  };
  auto bare = [&ridx](std::wstring_view lbl, std::wstring_view b,
                      std::wstring_view k) {
    return ex<Tensor>(lbl, bra{ridx(b)}, ket{ridx(k)});
  };

  auto const t_ = starred(L"t", L"a_1", L"i_1");
  auto const g_ = starred(L"g", L"i_1", L"a_2");
  auto const v_ = bare(L"v", L"a_2", L"a_1");
  auto const u_ = bare(L"u", L"a_2", L"a_1");

  auto root1 = eval_node(t_->clone() * g_->clone() * v_->clone());
  auto root2 = eval_node(t_->clone() * g_->clone() * u_->clone());
  using NodeT = decltype(root1);
  REQUIRE(root1->is_scalar());
  REQUIRE(root2->is_scalar());
  // the two-factor prefix hoists the conj; v and u break the three-factor one
  REQUIRE_FALSE(root1->canon_transform().conj);
  REQUIRE_FALSE(root2->canon_transform().conj);
  auto const& imed = root1.left();
  REQUIRE(imed->is_tensor());
  REQUIRE(imed->canon_transform() == CanonTransform{.conj = true});
  // one cell, shared by both roots
  REQUIRE(imed->hash_value() == root2.left()->hash_value());

  // the stored (bare) spellings the provider is asked for
  auto const& tl = imed.left();
  auto const& gl = imed.right();
  REQUIRE(tl.leaf());
  REQUIRE(gl.leaf());
  REQUIRE(tl->canon_transform() == CanonTransform{.conj = true});
  REQUIRE(gl->canon_transform() == CanonTransform{.conj = true});
  REQUIRE_FALSE(tl->as_tensor().kconjugated());
  REQUIRE_FALSE(gl->as_tensor().kconjugated());
  REQUIRE(tl->as_tensor().bra()[0].label() == L"a_1");
  REQUIRE(tl->as_tensor().ket()[0].label() == L"i_1");
  REQUIRE(gl->as_tensor().bra()[0].label() == L"i_1");
  REQUIRE(gl->as_tensor().ket()[0].label() == L"a_2");

  BTensorC T = rnd({nvirt, nocc});   // t{a_1;i_1}: T(a1, i)
  BTensorC G = rnd({nocc, nvirt});   // g{i_1;a_2}: G(i, a2)
  BTensorC V = rnd({nvirt, nvirt});  // v{a_2;a_1}: V(a2, a1)
  BTensorC U = rnd({nvirt, nvirt});  // u{a_2;a_1}: U(a2, a1)
  pinned_tensor_yield<BTensorC> yield;
  yield.put(tl->as_tensor(), T);
  yield.put(gl->as_tensor(), G);
  yield.put(v_->as<Tensor>(), V);
  yield.put(u_->as<Tensor>(), U);

  // A(a1, a2) = Σ_i conj(T(a1, i)) conj(G(i, a2)); the roots close it
  C ref1{0., 0.}, ref2{0., 0.};
  for (std::size_t a1 = 0; a1 < nvirt; ++a1)
    for (std::size_t a2 = 0; a2 < nvirt; ++a2) {
      C a{0., 0.};
      for (std::size_t i = 0; i < nocc; ++i)
        a += std::conj(T(a1, i)) * std::conj(G(i, a2));
      ref1 += a * V(a2, a1);
      ref2 += a * U(a2, a1);
    }

  auto scalar_of = [](ResultPtr const& res) {
    REQUIRE(res->is<ResultScalar<C>>());
    return res->as<ResultScalar<C>>().value();
  };
  auto close_to = [](C got, C ref) {
    CHECK(got.real() == Catch::Approx(ref.real()).margin(1e-12));
    CHECK(got.imag() == Catch::Approx(ref.imag()).margin(1e-12));
  };

  SECTION("the tree-walking engine agrees with the hand reference") {
    close_to(scalar_of(evaluate(root1, root1->annot(), yield)), ref1);
    close_to(scalar_of(evaluate(root2, root2->annot(), yield)), ref2);
  }

  SECTION("the ordered executor reads the shared conj node back conjugated") {
    BatchPolicy const policy;
    eval::dryrun::SizeRegime const regime;
    eval::dryrun::CostModel const cm{regime};
    auto const block_of = [](Index const&) -> std::size_t { return 1; };
    std::function<std::size_t(Index const&)> const target =
        [](Index const&) -> std::size_t { return 1; };

    container::svector<NodeT> const roots{root1, root2};
    auto const rich = eval::compute_dag_boulevard(roots, cm, block_of);
    auto const legality = eval::analyze_legality(rich, roots, policy);
    auto const ordered =
        eval::build_ordered_schedule(rich, legality, policy, {});
    auto cache = CacheManager<NodeT>::empty();
    // one schedule over both roots: the forest form sums them, so the
    // reference is ref1 + ref2
    auto const got = eval::evaluate_ordered_schedule(
        roots, ordered, rich, EvalExprBTAS::annot_t{}, yield, cache, target);
    close_to(scalar_of(got), ref1 + ref2);
  }
}

// A product that contracts every index of both operands is a bilinear dot,
// sum_k L[k] R[k], with no conjugation: SeQuant spells a conjugated operand
// explicitly (an adjointed or K-conjugated state, lowered to a leaf whose
// CanonTransform conjugates), so the backend must not conjugate on its own.
// BTAS's
// btas::dot is dotc (BLAS zdotc, the first operand conjugated), which is only
// observable on complex data.
TEST_CASE("eval_dot_complex_btas", "[eval_btas]") {
  using namespace sequant;
  using C = std::complex<double>;
  using BTensorC = btas::Tensor<C>;

  // arbitrary distinct annotation labels; x is a batch axis
  const long k = 3, x = 9;
  const std::size_t nk = 4, nx = 3;

  // deterministic fill with nonzero imaginary parts, distinct per element
  auto make = [](std::initializer_list<std::size_t> ext, double phase) {
    BTensorC t{btas::Range{container::svector<std::size_t>{ext}}};
    double v = 0.0;
    t.generate([&v, phase]() {
      v += 1.0;
      return C{std::cos(v * phase) * v / 7.0, std::sin(v * phase) * v / 5.0};
    });
    return t;
  };
  auto annots = [](container::svector<long> l, container::svector<long> r,
                   container::svector<long> c) {
    return std::array<std::any, 3>{std::move(l), std::move(r), std::move(c)};
  };
  auto prod = [](BTensorC const& L, BTensorC const& R,
                 std::array<std::any, 3> const& ann) {
    ResultPtr lr = eval_result<ResultTensorBTAS<BTensorC>>(L);
    ResultPtr rr = eval_result<ResultTensorBTAS<BTensorC>>(R);
    return lr->prod(*rr, ann, DeNest::False);
  };

  SECTION("full contraction: L[k] R[k] -> scalar") {
    auto L = make({nk}, 0.7), R = make({nk}, 1.3);
    auto res = prod(L, R, annots({k}, {k}, {}));
    REQUIRE(res->is<ResultScalar<C>>());
    C ref{0.0, 0.0};
    for (std::size_t ik = 0; ik < nk; ++ik) ref += L(ik) * R(ik);
    auto const got = res->as<ResultScalar<C>>().value();
    CHECK(got.real() == Catch::Approx(ref.real()).margin(1e-12));
    CHECK(got.imag() == Catch::Approx(ref.imag()).margin(1e-12));
  }

  SECTION("full contraction of the non-batch indices: L[k,x] R[k,x] -> C[x]") {
    auto L = make({nk, nx}, 0.7), R = make({nk, nx}, 1.3);
    auto res = prod(L, R, annots({k, x}, {k, x}, {x}));
    REQUIRE(res->is<ResultTensorBTAS<BTensorC>>());
    auto const& got = res->get<BTensorC>();
    REQUIRE(got.rank() == 1);
    REQUIRE(got.extent(0) == nx);
    for (std::size_t ix = 0; ix < nx; ++ix) {
      C ref{0.0, 0.0};
      for (std::size_t ik = 0; ik < nk; ++ik) ref += L(ik, ix) * R(ik, ix);
      CHECK(got(ix).real() == Catch::Approx(ref.real()).margin(1e-12));
      CHECK(got(ix).imag() == Catch::Approx(ref.imag()).margin(1e-12));
    }
  }
}

// A high-order hyperindex is an index shared among MORE than two tensor slots.
// When such an index is also *external* (named -- it survives into the result),
// TensorNetworkV3::canonicalize_slots used to assert that every named-index
// edge connects at most two vertices, so binarizing an expression carrying a
// batching/auxiliary index across many factors threw. This reproduces the
// Laplace-transform MP2 energy denominator, where the batching index z1 rides
// in the aux slot of every factor and of the result E{;;z1}.
//
// NB This only exercises tree construction (what the reported failure did):
// binarize<EvalExprBTAS>(deserialize<ResultExpr>(...)). Actually *evaluating*
// the tree is a separate matter -- every node is a batched (Hadamard)
// contraction over z1, which the BTAS backend's btas::contract (an index in
// both operands is always summed, never batched) does not support.
TEST_CASE("binarize_highorder_aux_hyperindex", "[eval_btas][hyperindex]") {
  using namespace sequant;

  // register the OBS AO space (μ) and the batching space (z) on top of the
  // standard single-reference spaces
  auto isr = mbpt::make_sr_spaces();
  mbpt::add_ao_spaces(isr, mbpt::Spin::any);  // μ (OBS AO)
  mbpt::add_batching_spaces(isr);             // z (batching)
  auto ctx_resetter =
      set_scoped_default_context(Context{get_default_context()}.set(isr));

  const std::wstring e_a_expr =
      L"E{;;z1} = w{;;z1} * g{μ5,μ6;μ7,μ8;z1} * c{μ1;i1;z1} * ć{i1;μ5;z1} * "
      L"d{μ7;μ3;z1} * c{μ2;i2;z1} * ć{i2;μ6;z1} * d{μ8;μ4;z1} * "
      L"g{μ3,μ4;μ1,μ2;z1}"
      L" - w{;;z1} * g{μ5,μ6;μ7,μ8;z1} * c{μ1;i1;z1} * ć{i1;μ5;z1} * "
      L"d{μ7;μ3;z1} * c{μ2;i2;z1} * ć{i2;μ6;z1} * d{μ8;μ4;z1} * "
      L"g{μ4,μ3;μ1,μ2;z1}";

  auto res = deserialize<ResultExpr>(
      e_a_expr, {.def_braket_symm = sequant::Hermiticity::NonHermitian});

  // the regression: this used to throw
  //   SEQUANT_ASSERT(edge_it->vertex_count() == 2) in canonicalize_slots
  REQUIRE_NOTHROW(binarize<EvalExprBTAS>(res));
  auto node = binarize<EvalExprBTAS>(res);

  const auto z1 = Index{L"z_1"};
  auto has_aux_z1 = [&z1](Tensor const& t) {
    return ranges::any_of(t.aux(), [&z1](Index const& ix) { return ix == z1; });
  };

  // the result carries z1 in its (only) aux slot ...
  REQUIRE(node->is_tensor());
  CHECK(node->as_tensor().aux_rank() == 1);
  CHECK(has_aux_z1(node->as_tensor()));

  // ... and z1 rides in the aux slot of every tensor node of the tree (leaves
  // and contraction intermediates alike): the batching index is never
  // contracted away.
  std::size_t tensor_nodes = 0;
  node.visit(
      [&](auto const& n) {
        if (n->is_tensor()) {
          ++tensor_nodes;
          CHECK(has_aux_z1(n->as_tensor()));
        }
      },
      TreeTraversal::PreOrder);
  // 9 leaves + 8 contraction intermediates per summand, plus the sum head and
  // the (-1) scaling node -- comfortably more than a handful of tensor nodes.
  CHECK(tensor_nodes > 10);
}

// Directly exercises ResultTensorBTAS::prod on a contraction that carries a
// batch (Hadamard) axis -- a label shared by both operands AND the result. A
// plain btas::contract cannot express one (it always sums an index common to
// both operands), so prod() detects the batch axis and slices over it. Every
// slice shape that batched_contract() must handle is covered: a genuine
// contraction, a dot (result-remainder empty), and a scalar*tensor (one
// operand-remainder empty).
TEST_CASE("btas_batched_contract", "[eval_btas][hyperindex]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  // arbitrary distinct annotation labels; x is the batch axis
  const long a = 1, b = 2, k = 3, x = 9;
  const std::size_t na = 2, nb = 4, nk = 3, nx = 5;

  auto make = [](std::initializer_list<std::size_t> ext) {
    BTensorD t{btas::Range{container::svector<std::size_t>{ext}}};
    // deterministic, non-trivial fill: distinct value per element
    double v = 0.0;
    t.generate([&v]() { return (v += 1.0) / 7.0; });
    return t;
  };
  auto annots = [](container::svector<long> l, container::svector<long> r,
                   container::svector<long> c) {
    return std::array<std::any, 3>{std::move(l), std::move(r), std::move(c)};
  };
  auto prod = [](BTensorD const& L, BTensorD const& R,
                 std::array<std::any, 3> const& ann) {
    ResultPtr lr = eval_result<ResultTensorBTAS<BTensorD>>(L);
    ResultPtr rr = eval_result<ResultTensorBTAS<BTensorD>>(R);
    return lr->prod(*rr, ann, DeNest::False)->get<BTensorD>();
  };

  SECTION("contraction + batch: L[a,k,x] R[k,b,x] -> C[a,b,x]") {
    auto L = make({na, nk, nx}), R = make({nk, nb, nx});
    auto C = prod(L, R, annots({a, k, x}, {k, b, x}, {a, b, x}));
    REQUIRE(C.rank() == 3);
    CHECK(C.extent(0) == na);
    CHECK(C.extent(1) == nb);
    CHECK(C.extent(2) == nx);
    for (std::size_t ia = 0; ia < na; ++ia)
      for (std::size_t ib = 0; ib < nb; ++ib)
        for (std::size_t ix = 0; ix < nx; ++ix) {
          double ref = 0.0;
          for (std::size_t ik = 0; ik < nk; ++ik)
            ref += L(ia, ik, ix) * R(ik, ib, ix);
          CHECK(C(ia, ib, ix) == Catch::Approx(ref));
        }
  }

  SECTION("dot + batch (empty result remainder): L[k,x] R[k,x] -> C[x]") {
    auto L = make({nk, nx}), R = make({nk, nx});
    auto C = prod(L, R, annots({k, x}, {k, x}, {x}));
    REQUIRE(C.rank() == 1);
    CHECK(C.extent(0) == nx);
    for (std::size_t ix = 0; ix < nx; ++ix) {
      double ref = 0.0;
      for (std::size_t ik = 0; ik < nk; ++ik) ref += L(ik, ix) * R(ik, ix);
      CHECK(C(ix) == Catch::Approx(ref));
    }
  }

  SECTION("scalar*tensor + batch (empty operand remainder): L[x] R[b,x]") {
    auto L = make({nx}), R = make({nb, nx});
    auto C = prod(L, R, annots({x}, {b, x}, {b, x}));
    REQUIRE(C.rank() == 2);
    CHECK(C.extent(0) == nb);
    CHECK(C.extent(1) == nx);
    for (std::size_t ib = 0; ib < nb; ++ib)
      for (std::size_t ix = 0; ix < nx; ++ix)
        CHECK(C(ib, ix) == Catch::Approx(L(ix) * R(ib, ix)));
  }

  SECTION(
      "result axis order differs from [batch,rest]: L[a,x] R[b,x] -> "
      "C[x,a,b]") {
    auto L = make({na, nx}), R = make({nb, nx});
    auto C = prod(L, R, annots({a, x}, {b, x}, {x, a, b}));
    REQUIRE(C.rank() == 3);
    CHECK(C.extent(0) == nx);
    CHECK(C.extent(1) == na);
    CHECK(C.extent(2) == nb);
    for (std::size_t ix = 0; ix < nx; ++ix)
      for (std::size_t ia = 0; ia < na; ++ia)
        for (std::size_t ib = 0; ib < nb; ++ib)
          CHECK(C(ix, ia, ib) == Catch::Approx(L(ia, ix) * R(ib, ix)));
  }
}

// End-to-end through binarize + evaluate: a batched-over-aux contraction whose
// batching index z1 is shared by 3 tensors (the >2-tensor case the TNv3 fix
// unlocks) and carried into the result. Every contraction node evaluates via
// the batched prod() path, and the answer must match a hand-computed reference.
//   R{;;z1} = A{;a1;z1} * B{a1;a2;z1} * C{a2;;z1}
//   R[z] = sum_{a1,a2} A[a1,z] * B[a1,a2,z] * C[a2,z]
TEST_CASE("eval_btas_batched_over_aux", "[eval_btas][hyperindex]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  auto isr = mbpt::make_sr_spaces();
  mbpt::add_batching_spaces(isr);  // z
  auto ctx_resetter =
      set_scoped_default_context(Context{get_default_context()}.set(isr));

  // Nonsymm so bra<->ket orientation is never folded by canonicalization: the
  // leaf identities (hence the yielder cache keys) stay put.
  const io::serialization::DeserializationOptions opts{
      .def_perm_symm = Symmetry::Nonsymm,
      .def_braket_symm = BraKetSymmetry::Nonsymm};

  std::srand(42);
  const size_t nocc = 2, nvirt = 3, nz = 5;
  aux_rand_tensor_yield<BTensorD> yield_{nocc, nvirt, nz};

  auto res = deserialize<ResultExpr>(
      L"R{;;z1} = A{;a1;z1} * B{a1;a2;z1} * C{a2;;z1}", opts);
  auto node = binarize<EvalExprBTAS>(res);

  auto got = evaluate(node, node->annot(), yield_)->get<BTensorD>();
  REQUIRE(got.rank() == 1);
  REQUIRE(got.extent(0) == nz);

  auto leaf = [&](std::wstring_view s) -> BTensorD const& {
    return yield_(deserialize<ExprPtr>(s, opts)->as<Tensor>())->get<BTensorD>();
  };
  auto const& A = leaf(L"A{;a1;z1}");    // [a1, z]
  auto const& B = leaf(L"B{a1;a2;z1}");  // [a1, a2, z]
  auto const& C = leaf(L"C{a2;;z1}");    // [a2, z]

  for (size_t z = 0; z < nz; ++z) {
    double ref = 0.0;
    for (size_t a1 = 0; a1 < nvirt; ++a1)
      for (size_t a2 = 0; a2 < nvirt; ++a2)
        ref += A(a1, z) * B(a1, a2, z) * C(a2, z);
    CHECK(got(z) == Catch::Approx(ref));
  }
}

TEST_CASE("eval_signed_network_btas", "[eval_btas]") {
  using namespace sequant;
  using C = std::complex<double>;
  using BTensorC = btas::Tensor<C>;

  // A leaf keeps its as-written orientation, and the states ride the leaf's
  // retrieval transform over the bare array: '⁺' as {conj, braket_swap} (the
  // adjoint of the stored array -- the bundles are exchanged in the leaf's
  // spelling, the modes keep their own labels), a '꙳' over a real basis as
  // {conj} alone (a pure elementwise conjugation, identity layout). A stated
  // tensor never carries a sign; signs live in scalars. Symbolic
  // canonicalization never exchanges a Conjugate tensor's bundles, so a
  // network keeps its value through it.
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  std::srand(2024);
  const size_t nocc = 2, nvirt = 3;
  auto rnd = []() {
    return C(static_cast<double>(std::rand()) / RAND_MAX,
             static_cast<double>(std::rand()) / RAND_MAX);
  };

  // d is anti-Hermitian: the two spellings are related exactly by
  // d{a_1;i_1} = -conj(d{i_1;a_1}), which the yielder below pins in the data
  BTensorC D{btas::Range{nocc, nvirt}};
  D.generate(rnd);
  BTensorC Dswapped{btas::Range{nvirt, nocc}};
  for (size_t i = 0; i < nocc; ++i)
    for (size_t a = 0; a < nvirt; ++a) Dswapped(a, i) = -std::conj(D(i, a));
  BTensorC U{btas::Range{nvirt, nocc}};  // generic u{a_1;i_2}
  U.generate(rnd);
  BTensorC V{btas::Range{nocc, nvirt}};  // generic v{i_1;a_1}
  V.generate(rnd);

  auto d = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(
        L"d", bra{b}, ket{k},
        TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
  };
  auto u2 = [] { return ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_2"}); };
  auto v = [] { return ex<Tensor>(L"v", bra{L"i_1"}, ket{L"a_1"}); };

  pinned_tensor_yield<BTensorC> yield;
  yield.put(d(L"i_1", L"a_1")->as<Tensor>(), D);
  yield.put(d(L"a_1", L"i_1")->as<Tensor>(), Dswapped);
  yield.put(u2()->as<Tensor>(), U);
  yield.put(v()->as<Tensor>(), V);

  // The engine returns the root's value in the requested layout; no node
  // here has a canonicalization phase, since a leaf keeps its orientation
  // and the networks are spelled in their canonical form.
  auto eval_open = [&yield](ExprPtr const& expr,
                            container::svector<long> const& layout) {
    auto node = eval_node(expr);
    REQUIRE(node->canon_phase() == 1);
    return evaluate(node, layout, yield)->get<BTensorC>();
  };

  SECTION("a contraction keeps its value through canonicalization") {
    // r{i_1;i_2} = d{i_1;a_1} u{a_1;i_2}; an open contraction, so the check
    // runs through btas::contract.
    BTensorC ref{btas::Range{nocc, nocc}};
    ref.fill(C{0., 0.});
    for (size_t i = 0; i < nocc; ++i)
      for (size_t j = 0; j < nocc; ++j)
        for (size_t a = 0; a < nvirt; ++a) ref(i, j) += D(i, a) * U(a, j);

    auto e = d(L"i_1", L"a_1") * u2();
    auto const layout = tidxs(L"i_1,i_2");
    auto const as_written = eval_open(e, layout);
    auto const canonicalized = eval_open(canonicalize(e->clone()), layout);
    for (size_t i = 0; i < nocc; ++i)
      for (size_t j = 0; j < nocc; ++j) {
        CHECK(as_written(i, j).real() ==
              Catch::Approx(ref(i, j).real()).margin(1e-12));
        CHECK(as_written(i, j).imag() ==
              Catch::Approx(ref(i, j).imag()).margin(1e-12));
        CHECK(canonicalized(i, j).real() ==
              Catch::Approx(ref(i, j).real()).margin(1e-12));
        CHECK(canonicalized(i, j).imag() ==
              Catch::Approx(ref(i, j).imag()).margin(1e-12));
      }
  }

  SECTION("a sum keeps its value through canonicalization") {
    auto e = d(L"i_1", L"a_1") + v();
    auto const got = eval_open(canonicalize(e->clone()), tidxs(L"i_1,a_1"));
    for (size_t i = 0; i < nocc; ++i)
      for (size_t a = 0; a < nvirt; ++a) {
        auto const expected = D(i, a) + V(i, a);
        CHECK(got(i, a).real() == Catch::Approx(expected.real()).margin(1e-12));
        CHECK(got(i, a).imag() == Catch::Approx(expected.imag()).margin(1e-12));
      }
  }

  // t{a_1;i_1} and its partner in a closed network; T(a, i) is t's array
  BTensorC T{btas::Range{nvirt, nocc}};
  T.generate(rnd);
  BTensorC W{btas::Range{nocc, nvirt}};  // W(i, a)
  W.generate(rnd);
  // Σ_{a,i} conj(T(a, i)) W(i, a)
  C ref_conj{0., 0.};
  for (size_t a = 0; a < nvirt; ++a)
    for (size_t i = 0; i < nocc; ++i) ref_conj += std::conj(T(a, i)) * W(i, a);
  // the root is checked against the reference without a cache, then twice
  // through one shared cache, so that a cache-slot collision between the
  // bare leaf and its stated spelling (t and t꙳ over a real basis, or t and
  // t⁺ -- each pair shares one slot) would surface on the second read
  auto check_closed = [&ref_conj](auto const& node, auto const& leaf_yield) {
    REQUIRE(node->is_scalar());
    REQUIRE(node->canon_phase() == 1);
    auto check = [&ref_conj](ResultPtr const& res) {
      REQUIRE(res->template is<ResultScalar<C>>());
      auto const got = res->template as<ResultScalar<C>>().value();
      CHECK(got.real() == Catch::Approx(ref_conj.real()).margin(1e-12));
      CHECK(got.imag() == Catch::Approx(ref_conj.imag()).margin(1e-12));
    };
    check(evaluate(node, node->annot(), leaf_yield));
    auto cache =
        cache_manager(std::array{node}, [](auto const&) { return false; });
    check(evaluate(node, node->annot(), leaf_yield, cache));
    check(evaluate(node, node->annot(), leaf_yield, cache));
  };

  SECTION("a K-conjugated leaf over a real basis is the conjugate array") {
    // over a real basis with complex data, t꙳{a_1;i_1} w{i_1;a_1} evaluates
    // to Σ conj(T(a, i)) W(i, a): the K-conjugate of t is the elementwise
    // conjugate of its array, served by a {conj} transform over the bare
    // leaf, which permutes nothing
    auto ridx = [](std::wstring_view label) {
      Index i(label);
      IndexSpace sp = i.space();
      sp.field(Field::Real);
      return Index(label, sp);
    };
    // parity None keeps the '꙳'; under the default Even it normalizes away
    Tensor t(L"t", bra{ridx(L"a_1")}, ket{ridx(L"i_1")},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor w(L"w", bra{ridx(L"i_1")}, ket{ridx(L"a_1")});
    REQUIRE(t.base_field() == Field::Real);
    pinned_tensor_yield<BTensorC> ryield;
    ryield.put(t, T);
    ryield.put(w, W);

    // the conjugate of a matrix element over a real basis is the K-conjugate
    // with the slots in place
    auto tk = conjugate(ex<Tensor>(t));
    REQUIRE(tk->is<Tensor>());
    REQUIRE(tk->as<Tensor>().kconjugated());
    REQUIRE_FALSE(tk->as<Tensor>().adjointed());
    REQUIRE(tk->as<Tensor>().bra()[0].label() == L"a_1");

    auto node = eval_node(tk * ex<Tensor>(w));
    REQUIRE(node->op_type() == EvalOp::Product);
    auto const& kn = node.left();
    REQUIRE(kn.leaf());
    REQUIRE(kn->canon_transform() == CanonTransform{.conj = true});
    REQUIRE_FALSE(kn->as_tensor().kconjugated());

    check_closed(node, ryield);
  }

  SECTION("an adjointed leaf over a complex basis rides the transform") {
    // over a complex basis conjugate(t{a_1;i_1}) is t⁺{i_1;a_1}, a leaf whose
    // transform conjugate-transposes: t⁺{i_1;a_1} w'{a_1;i_1} = Σ conj(T(a, i))
    // W'(a, i), with W'(a, i) = W(i, a) so that the reference is the same
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    Tensor w(L"w", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.base_field() == Field::Complex);
    BTensorC Wt{btas::Range{nvirt, nocc}};
    for (size_t a = 0; a < nvirt; ++a)
      for (size_t i = 0; i < nocc; ++i) Wt(a, i) = W(i, a);
    pinned_tensor_yield<BTensorC> cyield;
    cyield.put(t, T);
    cyield.put(w, Wt);

    auto tadj = adjoint(ex<Tensor>(t));
    REQUIRE(tadj->is<Tensor>());
    REQUIRE(tadj->as<Tensor>().adjointed());
    REQUIRE_FALSE(tadj->as<Tensor>().kconjugated());
    REQUIRE(tadj->as<Tensor>().bra()[0].label() == L"i_1");
    REQUIRE(*conjugate(ex<Tensor>(t)) == *tadj);

    auto node = eval_node(tadj * ex<Tensor>(w));
    REQUIRE(node->op_type() == EvalOp::Product);
    auto const& an = node.left();
    REQUIRE(an.leaf());
    REQUIRE_FALSE(an->as_tensor().adjointed());
    REQUIRE(an->as_tensor().bra()[0].label() == L"a_1");
    REQUIRE(an->canon_transform().conj);
    REQUIRE(an->canon_transform().braket_swap);

    check_closed(node, cyield);
  }
}
