#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"
#include "tensor_power_fixtures.hpp"

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/btas/array_ops.hpp>
#include <SeQuant/core/eval/backends/btas/eval_expr.hpp>
#include <SeQuant/core/eval/backends/btas/result.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
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
  return tensor_to_key(sequant::deserialize<sequant::ExprPtr>(
                           spec, {.def_perm_symm = sequant::Symmetry::Nonsymm})
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

  /// Serve \p r for every occurrence of tensor \p tnsr (by tensor_to_key);
  /// \p r must be in the tensor's slot layout.
  void set_leaf(sequant::Tensor const& tnsr, sequant::ResultPtr r) {
    label_to_tnsr_.insert_or_assign(tensor_to_key(tnsr), std::move(r));
  }

  [[nodiscard]] Tensor_t make_rand_tensor(sequant::Tensor const& tnsr) const {
    using ranges::views::repeat_n;
    using ranges::views::transform;
    using sequant::IndexSpace;
    auto ibr = sequant::get_default_context().index_basis_registry();

    SEQUANT_ASSERT(
        ranges::all_of(tnsr.const_braket_indices(),
                       [&ibr](auto const& idx) {
                         return idx.space() == ibr->retrieve_space(L"i") ||
                                idx.space() == ibr->retrieve_space(L"a");
                       }) &&
        "Unsupported IndexSpace type found while generating tensor.");

    auto rng = btas::Range{
        tnsr.const_braket_indices() | transform([this, &ibr](auto const& idx) {
          return idx.space() == ibr->retrieve_space(L"i") ? nocc_ : nvirt_;
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

  /// array ops sizing the occupied and the virtual space as this yield's
  /// tensors
  [[nodiscard]] sequant::BackendArrayOps array_ops() const {
    auto ibr = sequant::get_default_context().index_basis_registry();
    auto const bk = [&ibr](wchar_t const* l) {
      return std::wstring(ibr->retrieve_space(l).base_key());
    };
    return sequant::make_btas_array_ops<Tensor_t>(
        {{bk(L"i"), nocc_}, {bk(L"a"), nvirt_}});
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
    auto ibr = sequant::get_default_context().index_basis_registry();
    if (idx.space() == ibr->retrieve_space(L"i")) return nocc_;
    if (idx.space() == ibr->retrieve_space(L"a")) return nvirt_;
    if (idx.space() == ibr->retrieve_space(L"z")) return nz_;
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
        xpr, {.def_perm_symm = sequant::Symmetry::Antisymm});
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
    // result = 5/6 * identity - 1/6 * sum of the other ket permutations
    // (the idempotent NNS projector row)
    auto expr2 = parse_antisymm(L"R_{a1, a2, a3}^{i1, i2, i3}");
    auto tidx2 = tidxs(L"a_1,a_2,a_3,i_1,i_2,i_3");
    auto eval2 = eval_biorthogonal_nns_project(expr2, tidx2);
    auto const& r2 = yield(L"R{v,v,v;o,o,o}");

    BTensorD man2{r2.range()};
    man2.fill(0);
    man2 = r2;
    btas::scal(5.0 / 6.0, man2);

    BTensorD perm_sum{r2.range()};
    perm_sum.fill(0);

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

    REQUIRE(norm(eval1) == Catch::Approx(norm(eval2)));

    BTensorD zero2{eval2.range()};
    zero2 = eval1 - eval2;
    REQUIRE(norm(zero2) == Catch::Approx(0).margin(
                               100 * std::numeric_limits<double>::epsilon()));
  }
}

// A Sum node marked accumulate_in_place()
// (binarize()'s mark on the left-accumulator chain of a folded
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
// are STILL marked (a static property, unaffected), but only the two
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
// every mark forced off, reproducing the always-allocate path).
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

  // Sanity on the marking itself: the whole 3-Sum chain is marked.
  std::size_t total_sum = 0, inplace = 0;
  marked.visit([&](auto const& n) {
    if (!n->is_sum()) return;
    ++total_sum;
    if (n->accumulate_in_place()) ++inplace;
  });
  REQUIRE(total_sum == 3);
  REQUIRE(inplace == 3);

  // Force the SAME chain unmarked, reproducing the always-allocate
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

  // A Nonsymm-braket tensor's adjoint() marks the label with '⁺' and swaps
  // bra/ket; binarize lowers that to an EvalOp::Adjoint node, and evaluating it
  // must conjugate-transpose the operand. With genuinely complex data the
  // conjugation is observable (a missing conj would leave imaginary parts
  // unflipped — a pure transpose would still pass a norm-only check).
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor t_adj = t;
  t_adj.adjoint();
  REQUIRE(t_adj.label() == L"t⁺");

  auto node = eval_node(ex<Tensor>(t_adj));
  REQUIRE(node->op_type() == EvalOp::Adjoint);

  // operand tensor (bare 't{a_1;i_1}'): shape [nvirt, nocc], indexed (a, i)
  auto const& src = yield_(node.left()->as_tensor())->get<BTensorC>();
  REQUIRE(src.extent(0) == nvirt);
  REQUIRE(src.extent(1) == nocc);

  // adjoint result: shape [nocc, nvirt], indexed (i, a) == conj(src(a, i))
  auto const adj = evaluate(node, tidxs(t_adj), yield_)->get<BTensorC>();
  REQUIRE(adj.extent(0) == nocc);
  REQUIRE(adj.extent(1) == nvirt);

  for (size_t i = 0; i < nocc; ++i)
    for (size_t a = 0; a < nvirt; ++a) {
      auto const expected = std::conj(src(a, i));
      auto const got = adj(i, a);
      CHECK(got.real() == Catch::Approx(expected.real()).margin(1e-12));
      CHECK(got.imag() == Catch::Approx(expected.imag()).margin(1e-12));
    }
}

TEST_CASE("eval_tensor_power_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2024);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  // the harness fills tensors with values in [0, 1], so every exponent applies
  const TensorSymmetries dsym{.perm = Symmetry::Symm,
                              .braket = BraKetSymmetry::Symm,
                              .column = ColumnSymmetry::Symm};
  auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym);
  auto M = ex<Tensor>(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      TensorSymmetries{.perm = Symmetry::Antisymm});
  auto const target = tidxs(D->as<Tensor>());

  auto eval_at = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), target, yield_)->get<BTensorD>();
  };
  // every element of got == f(element of ref), in lockstep (same layout)
  auto check_elementwise = [](BTensorD const& got, BTensorD const& ref,
                              auto const& f) {
    REQUIRE(got.range() == ref.range());
    auto g = got.begin();
    for (auto r = ref.begin(); r != ref.end(); ++r, ++g)
      CHECK(*g == Catch::Approx(f(*r)).margin(1e-12));
  };

  auto node = eval_node(pow(D, rational{-1}));
  REQUIRE(node->op_type() == EvalOp::Power);
  auto const src = yield_(node.left()->as_tensor())->get<BTensorD>();
  auto const d = eval_at(D);

  check_elementwise(eval_at(pow(D, rational{-1})), d,
                    [](double x) { return 1 / x; });
  check_elementwise(eval_at(pow(D, rational{0})), d,
                    [](double) { return 1.0; });
  check_elementwise(eval_at(pow(D, rational{2})), d,
                    [](double x) { return x * x; });
  check_elementwise(eval_at(pow(D, rational{1, 2})), d,
                    [](double x) { return std::sqrt(x); });

  // ℓ = M ⊙ D^-1, a Hadamard product
  auto const m = eval_at(M);
  ResultExpr ell(Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                        TensorSymmetries{.perm = Symmetry::Antisymm}),
                 M * pow(D, rational{-1}));
  auto const got_ell =
      evaluate(binarize<EvalExprBTAS>(ell), target, yield_)->get<BTensorD>();
  {
    BTensorD expected = m;
    auto di = d.begin();
    for (auto& x : expected) x /= *di++;
    check_elementwise(got_ell, expected, [](double x) { return x; });
  }

  // E = L · (M ⊙ D^-1): every index appears on three tensors
  {
    auto ctx = get_default_context();
    ctx.set(AssertStrictBraKetSymmetry::No);
    auto _ = set_scoped_default_context(ctx);
    auto L = ex<Tensor>(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                        TensorSymmetries{.perm = Symmetry::Antisymm});
    ResultExpr E(Variable(L"E"), L * M * pow(D, rational{-1}));
    auto const got_E =
        evaluate(binarize<EvalExprBTAS>(E), container::svector<long>{}, yield_)
            ->get<double>();
    auto const l = eval_at(L);
    double expected = 0;
    auto mi = m.begin();
    auto di = d.begin();
    for (auto li = l.begin(); li != l.end(); ++li, ++mi, ++di)
      expected += *li * *mi / *di;
    CHECK(got_E == Catch::Approx(expected).margin(1e-10));
  }

  // M ⊙ D^-1 + D^-1, with D^-1 as either summand, evaluates correctly and
  // leaves the cached D untouched
  for (bool const power_first : {false, true}) {
    auto const hadamard = M * pow(D, rational{-1});
    auto const power = pow(D, rational{-1});
    ResultExpr sum(Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                          TensorSymmetries{.perm = Symmetry::Nonsymm}),
                   power_first ? power + hadamard : hadamard + power);
    auto const got_sum =
        evaluate(binarize<EvalExprBTAS>(sum), target, yield_)->get<BTensorD>();
    BTensorD expected = m;
    auto di = d.begin();
    for (auto& x : expected) {
      x = (x + 1) / *di;
      ++di;
    }
    check_elementwise(got_sum, expected, [](double x) { return x; });
    check_elementwise(yield_(node.left()->as_tensor())->get<BTensorD>(), src,
                      [](double x) { return x; });
  }

  // an antisymmetric base spelled out of canonical order: its leaf
  // canonicalizes with phase -1 and the power node carries that phase raised
  // to the exponent; in the base's own layout pow(A, e) is the element-wise
  // power of A as A evaluates there, the negative of its canonical spelling
  // (eval_leaf_reordered_with_a_sign_btas pins that sign).
  {
    auto A = ex<Tensor>(L"A", bra{L"a_1", L"i_1"}, ket{L"i_2", L"a_2"},
                        TensorSymmetries{.perm = Symmetry::Antisymm});
    REQUIRE(eval_node(pow(A, rational{2})).left()->canon_phase() == -1);
    REQUIRE(eval_node(pow(A, rational{2}))->canon_phase() == 1);
    REQUIRE(eval_node(pow(A, rational{-1}))->canon_phase() == -1);
    auto const a_target = tidxs(A->as<Tensor>());
    auto const at_a = [&](ExprPtr const& expr) {
      return evaluate(eval_node(expr), a_target, yield_)->get<BTensorD>();
    };
    auto const a = at_a(A);
    check_elementwise(at_a(pow(A, rational{2})), a,
                      [](double x) { return x * x; });
    check_elementwise(at_a(pow(A, rational{-1})), a,
                      [](double x) { return 1 / x; });
  }
}

TEST_CASE("eval_broadcast_sum_btas", "[eval_btas][broadcast]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2024);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};
  // a Broadcast node shapes the modes it adds through the array ops
  auto const aops = yield_.array_ops();
  auto eval_aops = [&](auto const& node, auto const& layout) {
    auto cache = CacheManager<std::remove_cvref_t<decltype(node)>>::empty();
    cache.set_array_ops(&aops);
    return evaluate(node, layout, yield_, cache);
  };

  auto eps_a = [](wchar_t const* a) { return ex<Tensor>(L"ε", bra{}, ket{a}); };
  auto eps_i = [](wchar_t const* i) { return ex<Tensor>(L"ε", bra{i}, ket{}); };
  // every ket-ε occurrence shares one yielded array by key, likewise every
  // bra-ε occurrence
  auto const ea = yield_(eps_a(L"a_1")->as<Tensor>())->get<BTensorD>();
  auto const ei = yield_(eps_i(L"i_1")->as<Tensor>())->get<BTensorD>();
  REQUIRE(ea.size() == nvirt);
  REQUIRE(ei.size() == nocc);

  const TensorSymmetries dsym{.perm = Symmetry::Symm,
                              .braket = BraKetSymmetry::Symm,
                              .column = ColumnSymmetry::Symm};
  Tensor const head(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym);
  ResultExpr res(head,
                 eps_a(L"a_1") + eps_a(L"a_2") - eps_i(L"i_1") - eps_i(L"i_2"));
  auto const got =
      eval_aops(binarize<EvalExprBTAS>(res), tidxs(head))->get<BTensorD>();
  REQUIRE(got.range() == btas::Range{nocc, nocc, nvirt, nvirt});
  for (std::size_t i1 = 0; i1 < nocc; ++i1)
    for (std::size_t i2 = 0; i2 < nocc; ++i2)
      for (std::size_t a1 = 0; a1 < nvirt; ++a1)
        for (std::size_t a2 = 0; a2 < nvirt; ++a2)
          CHECK(got(i1, i2, a1, a2) ==
                Catch::Approx(ea(a1) + ea(a2) - ei(i1) - ei(i2)).margin(1e-12));

  // distinct leaves and a nonsymmetric head: the expected tensor changes under
  // any mode swap, so this pins the result layout
  {
    auto const e1 = ex<Tensor>(L"ε1", bra{}, ket{L"a_1"});
    auto const e2 = ex<Tensor>(L"ε2", bra{}, ket{L"a_2"});
    auto const e3 = ex<Tensor>(L"ε3", bra{L"i_1"}, ket{});
    auto const e4 = ex<Tensor>(L"ε4", bra{L"i_2"}, ket{});
    auto vec = [&](ExprPtr const& e) {
      return yield_(e->as<Tensor>())->get<BTensorD>();
    };
    auto const v1 = vec(e1), v2 = vec(e2), v3 = vec(e3), v4 = vec(e4);
    const TensorSymmetries nsym{.perm = Symmetry::Nonsymm,
                                .braket = BraKetSymmetry::Nonsymm,
                                .column = ColumnSymmetry::Nonsymm};
    Tensor const nhead(L"D'", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, nsym);
    ResultExpr nres(nhead, e1 + e2 - e3 - e4);
    auto const ngot =
        eval_aops(binarize<EvalExprBTAS>(nres), tidxs(nhead))->get<BTensorD>();
    REQUIRE(ngot.range() == btas::Range{nocc, nocc, nvirt, nvirt});
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t i2 = 0; i2 < nocc; ++i2)
        for (std::size_t a1 = 0; a1 < nvirt; ++a1)
          for (std::size_t a2 = 0; a2 < nvirt; ++a2)
            CHECK(
                ngot(i1, i2, a1, a2) ==
                Catch::Approx(v1(a1) + v2(a2) - v3(i1) - v4(i2)).margin(1e-12));
  }

  // a constant summand is broadcast over every head mode
  {
    Tensor const chead(L"D", bra{L"i_1"}, ket{L"a_1"});
    ResultExpr cres(
        chead, eps_a(L"a_1") - eps_i(L"i_1") + ex<Constant>(rational{3, 2}));
    auto const cnode = binarize<EvalExprBTAS>(cres);
    REQUIRE(cnode.right()->is_broadcast());
    REQUIRE(cnode.right().left()->is_constant());
    auto const cgot = eval_aops(cnode, tidxs(chead))->get<BTensorD>();
    REQUIRE(cgot.range() == btas::Range{nocc, nvirt});
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t a1 = 0; a1 < nvirt; ++a1)
        CHECK(cgot(i1, a1) ==
              Catch::Approx(ea(a1) - ei(i1) + 1.5).margin(1e-12));
  }
}

TEST_CASE("eval_inlined_system_btas", "[eval_btas][inline]") {
  using namespace sequant;
  namespace tp = sequant::tests::tensor_power;
  namespace data = tp::data;
  using BTensorD = btas::Tensor<double>;
  using System = container::svector<ResultExpr>;
  using Flags = container::svector<bool>;

  // three occupied orbitals: with two, an antisymmetric pair leaves E{i_2;i_3}
  // diagonal
  const size_t nocc = 3, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};
  // a Broadcast node shapes the modes it adds through the array ops
  auto const aops = yield_.array_ops();
  auto eval_aops = [&](auto const& node, auto const& layout) {
    auto cache = CacheManager<std::remove_cvref_t<decltype(node)>>::empty();
    cache.set_array_ops(&aops);
    return evaluate(node, layout, yield_, cache);
  };
  // serves f(slot ordinals of \p tnsr) for every occurrence of the leaf
  // \p tnsr, in either orientation: a leaf node of a bra-ket symmetric tensor
  // may hold either
  auto pin = [&](ExprPtr const& tnsr, auto f) {
    Tensor const& t = tnsr->as<Tensor>();
    auto const nb = t.bra_rank();
    auto pin_as = [&](Tensor const& as, bool swapped) {
      std::vector<std::size_t> extents;
      for (Index const& ix : as.const_braket_indices())
        extents.push_back(ix.space().base_key() == L"i" ? nocc : nvirt);
      BTensorD arr{btas::Range{extents}};
      auto it = arr.begin();
      for (auto const& idx : arr.range()) {
        std::vector<std::size_t> x(idx.begin(), idx.end());
        if (swapped) std::rotate(x.begin(), x.end() - nb, x.end());
        *it++ = f(x);
      }
      yield_.set_leaf(as, eval_result<ResultTensorBTAS<BTensorD>>(arr));
    };
    pin_as(t, false);
    pin_as(Tensor(t.label(), bra(t.ket()), ket(t.bra()), aux(t.aux()),
                  TensorSymmetries{.perm = t.symmetry(),
                                   .braket = t.braket_symmetry(),
                                   .column = t.column_symmetry()}),
           true);
  };
  using Ordinals = std::vector<std::size_t>;
  pin(tp::antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
      [](Ordinals const& x) { return data::M(x[0], x[1], x[2], x[3]); });
  pin(tp::contractor(),
      [](Ordinals const& x) { return data::L(x[0], x[1], x[2], x[3]); });

  auto const e_layout = tidxs(tp::head(L"E"));
  // the value of the one root of a system
  auto eval_system = [&](System const& sys) {
    auto const res = binarize<EvalExprBTAS>(sys);
    REQUIRE(res.roots().size() == 1);
    return eval_aops(res.roots()[0].get(), e_layout)->get<BTensorD>();
  };
  // E{i_2;i_3} of consumer_E()
  auto check_E = [&](BTensorD const& got) {
    REQUIRE(got.range() == btas::Range{nocc, nocc});
    for (std::size_t i2 = 0; i2 < nocc; ++i2)
      for (std::size_t i3 = 0; i3 < nocc; ++i3)
        CHECK(got(i2, i3) == Catch::Approx(data::contract_with_ell(
                                               data::L, i3, i2, nocc, nvirt))
                                 .margin(1e-12));
  };

  SECTION("inlined E equals the driver-chained E and the direct form") {
    pin(tp::power_base(),
        [](Ordinals const& x) { return data::D(x[0], x[1], x[2], x[3]); });
    System const sys{tp::ell_definition(), tp::consumer_E()};
    REQUIRE(binarize(sys).consumed == Flags{true, false});
    auto const inlined = eval_system(sys);
    check_E(inlined);

    // ℓ evaluated by itself, then served as the leaf E reads, in either
    // orientation (see pin)
    for (Tensor const& use :
         {tp::relabelled(tp::ell_definition().result_as_tensor(),
                         bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
          tp::relabelled(tp::ell_definition().result_as_tensor(),
                         bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"})})
      yield_.set_leaf(use,
                      evaluate(binarize<EvalExprBTAS>(tp::ell_definition()),
                               tidxs(use), yield_));
    auto same = [&inlined](BTensorD const& got) {
      REQUIRE(got.range() == inlined.range());
      auto g = got.begin();
      for (auto r = inlined.begin(); r != inlined.end(); ++r, ++g)
        CHECK(*g == Catch::Approx(*r).margin(1e-12));
    };
    same(evaluate(binarize<EvalExprBTAS>(tp::consumer_E()), e_layout, yield_)
             ->get<BTensorD>());
    // the direct form contracts L with M and D in one network, which the
    // strict bra-ket rule rejects
    auto _ = set_scoped_default_context(tp::lax_braket_context());
    same(evaluate(binarize<EvalExprBTAS>(tp::hadamard_power_term()), e_layout,
                  yield_)
             ->get<BTensorD>());
  }

  SECTION("three levels: ε into D into ℓ into E") {
    pin(ex<Tensor>(L"ε", bra{}, ket{L"a_1"}),
        [](Ordinals const& x) { return data::eps_virt(x[0]); });
    pin(ex<Tensor>(L"ε", bra{L"i_1"}, ket{}),
        [](Ordinals const& x) { return data::eps_occ(x[0]); });
    System const sys{tp::D_broadcast(), tp::ell_definition(), tp::consumer_E()};
    REQUIRE(binarize(sys).consumed == Flags{true, true, false});
    check_E(eval_system(sys));
  }
}

// A product that contracts every index of both operands is a bilinear dot,
// sum_k L[k] R[k], with no conjugation: SeQuant spells a conjugated operand
// explicitly (a Conjugate value modifier, an EvalOp::Adjoint node), so the
// backend must not conjugate on its own. BTAS's btas::dot is dotc (BLAS zdotc,
// the first operand conjugated), which is only observable on complex data.
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

  SECTION("batch-full-contraction") {
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
  auto ibr = mbpt::make_sr_spaces();
  mbpt::add_ao_basis(ibr, mbpt::Spin::any);  // μ (OBS AO)
  mbpt::add_batching_spaces(ibr);            // z (batching)
  auto ctx_resetter =
      set_scoped_default_context(Context{get_default_context()}.set(ibr));

  const std::wstring e_a_expr =
      L"E{;;z1} = w{;;z1} * g{μ5,μ6;μ7,μ8;z1} * c{μ1;i1;z1} * ć{i1;μ5;z1} * "
      L"d{μ7;μ3;z1} * c{μ2;i2;z1} * ć{i2;μ6;z1} * d{μ8;μ4;z1} * "
      L"g{μ3,μ4;μ1,μ2;z1}"
      L" - w{;;z1} * g{μ5,μ6;μ7,μ8;z1} * c{μ1;i1;z1} * ć{i1;μ5;z1} * "
      L"d{μ7;μ3;z1} * c{μ2;i2;z1} * ć{i2;μ6;z1} * d{μ8;μ4;z1} * "
      L"g{μ4,μ3;μ1,μ2;z1}";

  auto res = deserialize<ResultExpr>(e_a_expr);

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

  SECTION("batch-result-axis-order") {
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

  auto ibr = mbpt::make_sr_spaces();
  mbpt::add_batching_spaces(ibr);  // z
  auto ctx_resetter =
      set_scoped_default_context(Context{get_default_context()}.set(ibr));

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

TEST_CASE("result_size_in_bytes_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  BTensorD t{btas::Range{container::svector<std::size_t>{2, 3}}};
  ResultPtr r = eval_result<ResultTensorBTAS<BTensorD>>(std::move(t));
  REQUIRE(r->size_in_bytes() == 6 * sizeof(double));
}

TEST_CASE("eval_btas_annotation_tells_basis_instances_apart",
          "[eval_btas][basis]") {
  using namespace sequant;
  const Index a1(L"a_1");
  const std::vector<Index> idxs{a1, a1.replace_basis_instance(1),
                                a1.replace_basis_instance(2)};
  const auto annot =
      sequant::EvalExprBTAS::index_hash(idxs) | ranges::to<std::vector<long>>;
  CHECK(annot[0] != annot[1]);
  CHECK(annot[0] != annot[2]);
  CHECK(annot[1] != annot[2]);
}

TEST_CASE("eval_leaf_reordered_with_a_sign_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2026);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  // an antisymmetric leaf whose bra holds two spaces: the canonical order of
  // the bundle is i_3, a_1, so the leaf as written is its negative
  const TensorSymmetries antisymm{.perm = Symmetry::Antisymm};
  auto written =
      ex<Tensor>(L"t", bra{L"a_1", L"i_3"}, ket{L"i_1", L"i_2"}, antisymm);
  auto sorted =
      ex<Tensor>(L"t", bra{L"i_3", L"a_1"}, ket{L"i_1", L"i_2"}, antisymm);
  // one layout for both: a_1, i_3, i_1, i_2
  auto const target = tidxs(written->as<Tensor>());
  auto eval_at = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), target, yield_)->get<BTensorD>();
  };

  auto const t_written = eval_at(written);
  auto const t_sorted = eval_at(sorted);
  REQUIRE(t_written.range() == t_sorted.range());
  double norm = 0;
  auto s = t_sorted.begin();
  for (auto w = t_written.begin(); w != t_written.end(); ++w, ++s) {
    norm += *s * *s;
    CHECK(*w == Catch::Approx(-*s).margin(1e-12));
  }
  REQUIRE(norm > 0.1);

  // and under a contraction
  auto W = ex<Tensor>(L"W", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"},
                      TensorSymmetries{.perm = Symmetry::Nonsymm});
  auto scalar_of = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), yield_)->get<double>();
  };
  auto const e_written = scalar_of(written * W);
  auto const e_sorted = scalar_of(sorted * W);
  REQUIRE(std::abs(e_sorted) > 1e-3);
  CHECK(e_written == Catch::Approx(-e_sorted).margin(1e-12));
}

TEST_CASE("eval_cached_intermediate_of_reordered_leaf_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2026);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  // the two spellings of eval_leaf_reordered_with_a_sign_btas: the leaf as
  // written is the negative of the sorted one, and their nodes hash alike
  const TensorSymmetries antisymm{.perm = Symmetry::Antisymm};
  const TensorSymmetries nonsymm{.perm = Symmetry::Nonsymm};
  auto written =
      ex<Tensor>(L"t", bra{L"a_1", L"i_3"}, ket{L"i_1", L"i_2"}, antisymm);
  auto sorted =
      ex<Tensor>(L"t", bra{L"i_3", L"a_1"}, ket{L"i_1", L"i_2"}, antisymm);
  auto scalar_of = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), yield_)->get<double>();
  };
  // evaluates both expressions with one cache, so that an intermediate of the
  // first is the second's if their nodes hash alike
  auto with_shared_cache = [&](ExprPtr const& first, ExprPtr const& second) {
    std::vector nodes{eval_node(first), eval_node(second)};
    auto cache = cache_manager(nodes);
    auto const v1 = evaluate(nodes[0], yield_, cache)->get<double>();
    auto const v2 = evaluate(nodes[1], yield_, cache)->get<double>();
    return std::pair{v1, v2};
  };

  SECTION("a product") {
    auto W =
        ex<Tensor>(L"W", bra{L"i_1", L"i_2"}, ket{L"i_3", L"a_2"}, nonsymm);
    auto V = ex<Tensor>(L"V", bra{L"a_2"}, ket{L"a_1"}, nonsymm);
    // the nodes of the two t W have one hash and opposite phases
    auto const n_sorted = eval_node(sorted * W * V);
    auto const n_written = eval_node(written * W * V);
    REQUIRE(n_sorted.left()->hash_value() == n_written.left()->hash_value());
    REQUIRE(n_sorted.left()->canon_phase() == -n_written.left()->canon_phase());

    auto const e = scalar_of(sorted * W * V);
    REQUIRE(std::abs(e) > 1e-3);
    for (bool sorted_first : {true, false}) {
      auto const [e1, e2] =
          sorted_first ? with_shared_cache(sorted * W * V, written * W * V)
                       : with_shared_cache(written * W * V, sorted * W * V);
      CHECK(e1 == Catch::Approx(sorted_first ? e : -e).margin(1e-12));
      CHECK(e2 == Catch::Approx(sorted_first ? -e : e).margin(1e-12));
    }
  }

  SECTION("a sum") {
    auto Y =
        ex<Tensor>(L"Y", bra{L"a_1", L"i_3"}, ket{L"i_1", L"i_2"}, nonsymm);
    auto Z =
        ex<Tensor>(L"Z", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"}, nonsymm);
    // the two sums differ in the sign of a summand, hence in value
    REQUIRE(eval_node((sorted + Y) * Z).left()->hash_value() !=
            eval_node((written + Y) * Z).left()->hash_value());

    auto const t = scalar_of(sorted * Z);
    auto const y = scalar_of(Y * Z);
    REQUIRE(std::abs(t) > 1e-3);
    for (bool sorted_first : {true, false}) {
      auto const [e1, e2] =
          sorted_first ? with_shared_cache((sorted + Y) * Z, (written + Y) * Z)
                       : with_shared_cache((written + Y) * Z, (sorted + Y) * Z);
      CHECK(e1 == Catch::Approx(sorted_first ? y + t : y - t).margin(1e-12));
      CHECK(e2 == Catch::Approx(sorted_first ? y - t : y + t).margin(1e-12));
    }
  }
}

TEST_CASE("eval_ordered_reordered_leaf_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2026);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  // the two spellings of eval_leaf_reordered_with_a_sign_btas as roots of one
  // schedule: the leaf and t W are one value each, read with opposite signs
  const TensorSymmetries antisymm{.perm = Symmetry::Antisymm};
  const TensorSymmetries nonsymm{.perm = Symmetry::Nonsymm};
  auto written =
      ex<Tensor>(L"t", bra{L"a_1", L"i_3"}, ket{L"i_1", L"i_2"}, antisymm);
  auto sorted =
      ex<Tensor>(L"t", bra{L"i_3", L"a_1"}, ket{L"i_1", L"i_2"}, antisymm);
  auto W = ex<Tensor>(L"W", bra{L"i_1", L"i_2"}, ket{L"i_3", L"a_2"}, nonsymm);
  auto V = ex<Tensor>(L"V", bra{L"a_2"}, ket{L"a_1"}, nonsymm);

  auto const e = evaluate(eval_node(sorted * W * V), yield_)->get<double>();
  REQUIRE(std::abs(e) > 1e-3);

  BatchPolicy policy;
  policy.scheduler = BatchScheduler::ordered;
  for (bool sorted_first : {true, false}) {
    std::vector nodes{eval_node((sorted_first ? sorted : written) * W * V),
                      eval_node((sorted_first ? written : sorted) * W * V)};
    REQUIRE(nodes[0].left()->hash_value() == nodes[1].left()->hash_value());
    auto cache = cache_manager(nodes);
    auto const results = evaluate_each(nodes, policy, yield_, cache);
    REQUIRE(results.size() == 2);
    CHECK(results[0]->get<double>() ==
          Catch::Approx(sorted_first ? e : -e).margin(1e-12));
    CHECK(results[1]->get<double>() ==
          Catch::Approx(sorted_first ? -e : e).margin(1e-12));
  }
}

TEST_CASE("eval_adjoint_of_reordered_leaf_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;

  std::srand(2026);
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};

  // the adjoints of the two spellings of eval_leaf_reordered_with_a_sign_btas
  // differ by the sign their operands do; their nodes hash alike but are not
  // one value (their layouts, as written, differ)
  auto adjoint_of = [](std::initializer_list<std::wstring> bra_labels) {
    Tensor t(L"t", bra(bra_labels), ket{L"i_1", L"i_2"}, Symmetry::Antisymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    t.adjoint();
    return ex<Tensor>(std::move(t));
  };
  auto written = adjoint_of({L"a_1", L"i_3"});
  auto sorted = adjoint_of({L"i_3", L"a_1"});
  auto W = ex<Tensor>(L"W", bra{L"a_1", L"i_3"}, ket{L"i_1", L"i_2"},
                      TensorSymmetries{.perm = Symmetry::Nonsymm});
  REQUIRE(eval_node(written)->op_type() == EvalOp::Adjoint);

  auto scalar_of = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), yield_)->get<double>();
  };
  auto const e = scalar_of(sorted * W);
  REQUIRE(std::abs(e) > 1e-3);
  CHECK(scalar_of(written * W) == Catch::Approx(-e).margin(1e-12));

  // each spelling twice, so that every node is a value of the cache or of the
  // schedule
  for (bool ordered : {false, true})
    for (bool sorted_first : {true, false}) {
      CAPTURE(ordered, sorted_first);
      std::vector nodes{eval_node(sorted_first ? sorted : written),
                        eval_node(sorted_first ? written : sorted),
                        eval_node(sorted_first ? sorted : written),
                        eval_node(sorted_first ? written : sorted)};
      REQUIRE(nodes[0]->hash_value() == nodes[1]->hash_value());
      auto cache = cache_manager(nodes);
      BatchPolicy policy;
      if (ordered) policy.scheduler = BatchScheduler::ordered;
      auto const results = evaluate_each(nodes, policy, yield_, cache);
      for (int k : {0, 1}) {
        // i_1, i_2, i_3, a_1 in one and i_1, i_2, a_1, i_3 in the other
        auto const& s =
            results[sorted_first ? k * 2 : k * 2 + 1]->get<BTensorD>();
        auto const& w =
            results[sorted_first ? k * 2 + 1 : k * 2]->get<BTensorD>();
        REQUIRE(s.extent(2) == nocc);
        REQUIRE(w.extent(2) == nvirt);
        for (size_t i = 0; i < nocc; ++i)
          for (size_t j = 0; j < nocc; ++j)
            for (size_t l = 0; l < nocc; ++l)
              for (size_t c = 0; c < nvirt; ++c)
                CHECK(w(i, j, c, l) ==
                      Catch::Approx(-s(i, j, l, c)).margin(1e-12));
      }
    }
}

TEST_CASE("eval_inlined_equation_phase_btas", "[eval_btas][inline]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;
  using Labels = std::initializer_list<std::wstring>;

  std::srand(2026);
  const size_t nocc = 3, nvirt = 4;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};
  const TensorSymmetries antisymm{.perm = Symmetry::Antisymm};
  const TensorSymmetries nonsymm{.perm = Symmetry::Nonsymm};

  // the scalar roots of a system, alone, with one cache and as one schedule,
  // against \p expected
  auto check = [&](container::svector<ResultExpr> const& sys,
                   std::array<double, 2> const& expected) {
    auto const res = binarize<EvalExprBTAS>(sys);
    std::vector<EvalNode<EvalExprBTAS>> nodes;
    for (auto const& r : res.roots())
      if (r.get()->is_scalar()) nodes.push_back(r.get());
    REQUIRE(nodes.size() == 2);
    // the two consumers are one canonical value, of opposite phases
    REQUIRE(nodes[0].left()->hash_value() == nodes[1].left()->hash_value());
    REQUIRE(nodes[0].left()->canon_phase() == -nodes[1].left()->canon_phase());
    auto cache = cache_manager(nodes);
    BatchPolicy ordered;
    ordered.scheduler = BatchScheduler::ordered;
    auto ordered_cache = cache_manager(nodes);
    auto const scheduled = evaluate_each(nodes, ordered, yield_, ordered_cache);
    for (std::size_t k : {0, 1}) {
      CAPTURE(k);
      CHECK(evaluate(nodes[k], yield_)->get<double>() ==
            Catch::Approx(expected[k]).margin(1e-10));
      CHECK(evaluate(nodes[k], yield_, cache)->get<double>() ==
            Catch::Approx(expected[k]).margin(1e-10));
      CHECK(scheduled[k]->get<double>() ==
            Catch::Approx(expected[k]).margin(1e-10));
    }
  };
  auto scalar_of = [&](ExprPtr const& expr) {
    return evaluate(eval_node(expr), yield_)->get<double>();
  };

  SECTION("a tensor equation read with its head's labels permuted") {
    // I is antisymmetric, so I{a_2,a_1;...} is its negative: the inlined
    // subtrees have one hash and opposite phases, and the head, as written,
    // carries that sign into the network of the consumer
    auto antisymmetric = [&](ExprPtr const& t, std::size_t n01,
                             std::size_t n23) {
      btas::Range const range{n01, n01, n23, n23};
      BTensorD r{range}, a{range};
      for (auto& x : r) x = double(std::rand()) / RAND_MAX - 0.5;
      for (auto const& x : range)
        a(x[0], x[1], x[2], x[3]) =
            r(x[0], x[1], x[2], x[3]) - r(x[1], x[0], x[2], x[3]) -
            r(x[0], x[1], x[3], x[2]) + r(x[1], x[0], x[3], x[2]);
      yield_.set_leaf(t->as<Tensor>(),
                      eval_result<ResultTensorBTAS<BTensorD>>(a));
    };
    auto I = [&](Labels b) {
      return ex<Tensor>(L"I", bra(b), ket{L"i_1", L"i_2"}, antisymm);
    };
    auto u =
        ex<Tensor>(L"u", bra{L"a_1", L"a_2"}, ket{L"a_3", L"a_4"}, antisymm);
    auto v =
        ex<Tensor>(L"v", bra{L"a_3", L"a_4"}, ket{L"i_1", L"i_2"}, antisymm);
    antisymmetric(u, nvirt, nvirt);
    antisymmetric(v, nvirt, nocc);
    auto W =
        ex<Tensor>(L"W", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_5"}, nonsymm);
    auto Z = ex<Tensor>(L"Z", bra{L"a_5"}, ket{L"a_2"}, nonsymm);

    auto const e = scalar_of(u * v * W * Z);
    REQUIRE(std::abs(e) > 1e-3);
    Labels const a12{L"a_1", L"a_2"}, a21{L"a_2", L"a_1"};
    for (bool as_defined_first : {true, false}) {
      CAPTURE(as_defined_first);
      Labels const first = as_defined_first ? a12 : a21;
      Labels const second = as_defined_first ? a21 : a12;
      check({ResultExpr(I({L"a_1", L"a_2"})->as<Tensor>(), u * v),
             ResultExpr(Variable(L"R1"), I(first) * W * Z),
             ResultExpr(Variable(L"R2"), I(second) * W * Z)},
            {as_defined_first ? e : -e, as_defined_first ? -e : e});
    }
  }

  SECTION("scalar equations that differ by a sign") {
    // the two t W have one hash and opposite phases (see
    // eval_cached_intermediate_of_reordered_leaf_btas), and a scalar is no
    // factor of the network of its consumer
    auto t = [&](Labels b) {
      return ex<Tensor>(L"t", bra(b), ket{L"i_1", L"i_2"}, antisymm);
    };
    auto W =
        ex<Tensor>(L"W", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"}, nonsymm);
    auto X = ex<Tensor>(L"X", bra{L"a_2"}, ket{L"i_4"}, nonsymm);
    auto Y = ex<Tensor>(L"Y", bra{L"i_4"}, ket{L"a_3"}, nonsymm);
    auto Z = ex<Tensor>(L"Z", bra{L"a_3"}, ket{L"a_2"}, nonsymm);
    auto E = [](std::wstring label) { return ex<Variable>(label); };

    auto const e = scalar_of(t({L"a_1", L"i_3"}) * W) * scalar_of(X * Y * Z);
    REQUIRE(std::abs(e) > 1e-3);
    for (bool written_first : {true, false}) {
      CAPTURE(written_first);
      check({ResultExpr(Variable(L"E1"), t({L"a_1", L"i_3"}) * W),
             ResultExpr(Variable(L"E2"), t({L"i_3", L"a_1"}) * W),
             ResultExpr(Variable(L"R1"),
                        E(written_first ? L"E1" : L"E2") * X * Y * Z),
             ResultExpr(Variable(L"R2"),
                        E(written_first ? L"E2" : L"E1") * X * Y * Z)},
            {written_first ? e : -e, written_first ? -e : e});
    }
  }
}

TEST_CASE("eval_slot_alignment_btas", "[eval_btas]") {
  using namespace sequant;
  using BTensorD = btas::Tensor<double>;
  using Labels = std::initializer_list<std::wstring>;

  std::srand(2026);
  const size_t nocc = 3, nvirt = 4;
  auto yield_ = rand_tensor_yield<BTensorD>{nocc, nvirt};
  // no tensor has a symmetry: the two orders of a_1, a_2 are different values
  auto tensor = [](std::wstring const& label, Labels b, Labels k) {
    return ex<Tensor>(label, bra(b), ket(k),
                      TensorSymmetries{.perm = Symmetry::Nonsymm});
  };
  Labels const a12{L"a_1", L"a_2"}, a21{L"a_2", L"a_1"}, i12{L"i_1", L"i_2"};
  auto g = [&](Labels b) { return tensor(L"g", b, i12); };
  auto h = [&](Labels b) { return tensor(L"h", b, i12); };
  auto const W = tensor(L"W", i12, a12);
  auto const X = tensor(L"X", {L"a_3"}, {L"i_3"});
  auto const Y = tensor(L"Y", {L"i_3"}, {L"a_3"});
  // (f W) (X Y): f W is an intermediate of the tree
  auto consumer = [&](ExprPtr const& f) {
    return eval_node(ex<Product>(
        ExprPtrList{ex<Product>(ExprPtrList{f, W}, Product::Flatten::No),
                    X * Y},
        Product::Flatten::No));
  };
  // two trees that differ in value evaluate with one cache as each does alone
  auto check_nodes = [&](EvalNode<EvalExprBTAS> const& n1,
                         EvalNode<EvalExprBTAS> const& n2) {
    std::vector nodes{n1, n2};
    std::array const alone{evaluate(nodes[0], yield_)->get<double>(),
                           evaluate(nodes[1], yield_)->get<double>()};
    REQUIRE(std::abs(alone[0] - alone[1]) > 1e-3);
    auto cache = cache_manager(nodes);
    for (std::size_t k : {0, 1})
      CHECK(evaluate(nodes[k], yield_, cache)->get<double>() ==
            Catch::Approx(alone[k]).margin(1e-10));
  };
  auto check = [&](ExprPtr const& f1, ExprPtr const& f2) {
    check_nodes(consumer(f1), consumer(f2));
  };

  SECTION("summands aligned differently") {
    check(g(a12) + h(a12), g(a12) + h(a21));
    check(g(a12) + h(a12) + g(a21), g(a12) + h(a21) + g(a21));
  }

  SECTION("factors that differ in the labels of their slots") {
    check(g(a12) + h(a12), g(a21) + h(a21));
    check(pow(g(a12), rational{3}), pow(g(a21), rational{3}));
  }

  SECTION("two factors that are results of operations") {
    // (p + q) (r + s) and (r + s) (p + q), contracted alike
    auto m = [&](std::wstring const& label, Labels b, Labels k) {
      return tensor(label, b, k);
    };
    auto pq = [&](Labels b, Labels k) { return m(L"p", b, k) + m(L"q", b, k); };
    auto rs = [&](Labels b, Labels k) { return m(L"r", b, k) + m(L"s", b, k); };
    auto pw = [&](std::wstring const& label, Labels b, Labels k, int e) {
      return pow(m(label, b, k), rational{e});
    };
    auto const Z = tensor(L"Z", {L"a_5"}, {L"a_1"});
    Labels const a1{L"a_1"}, a4{L"a_4"}, a5{L"a_5"};
    auto two = [&](ExprPtr const& f1, ExprPtr const& f2) {
      return eval_node(ex<Product>(
          ExprPtrList{ex<Product>(ExprPtrList{f1, f2}, Product::Flatten::No),
                      Z},
          Product::Flatten::No));
    };
    check_nodes(two(pq(a1, a4), rs(a4, a5)), two(rs(a1, a4), pq(a4, a5)));
    check_nodes(two(pw(L"p", a1, a4, 2), pw(L"q", a4, a5, 3)),
                two(pw(L"q", a1, a4, 3), pw(L"p", a4, a5, 2)));
  }

  SECTION("bra and ket of a factor exchanged by a bra-ket symmetry") {
    // f and e are symmetric matrices: the sums f e + q of the two trees are
    // one value, laid out a_1, a_2 in one and a_2, a_1 in the other, so a_1 is
    // the bra of one and the ket of the other
    auto symmetric = [&](std::wstring const& label, Labels b, Labels k) {
      auto t = ex<Tensor>(label, bra(b), ket(k),
                          TensorSymmetries{.perm = Symmetry::Nonsymm,
                                           .braket = BraKetSymmetry::Symm});
      BTensorD r{btas::Range{nvirt, nvirt}}, m{btas::Range{nvirt, nvirt}};
      for (auto& x : r) x = double(std::rand()) / RAND_MAX - 0.5;
      for (std::size_t c = 0; c < nvirt; ++c)
        for (std::size_t d = 0; d < nvirt; ++d) m(c, d) = r(c, d) + r(d, c);
      yield_.set_leaf(t->as<Tensor>(),
                      eval_result<ResultTensorBTAS<BTensorD>>(m));
      return t;
    };
    Labels const a1{L"a_1"}, a2{L"a_2"}, a6{L"a_6"};
    auto const Z = tensor(L"Z", a2, a1);
    auto with_Z = [&](ExprPtr const& f) {
      return eval_node(ex<Product>(
          ExprPtrList{ex<Product>(ExprPtrList{f, Z}, Product::Flatten::No),
                      X * Y},
          Product::Flatten::No));
    };
    check_nodes(with_Z(symmetric(L"f", a1, a6) * symmetric(L"e", a6, a2) +
                       tensor(L"q", a1, a2)),
                with_Z(symmetric(L"f", a6, a2) * symmetric(L"e", a1, a6) +
                       tensor(L"q", a2, a1)));
  }
}
