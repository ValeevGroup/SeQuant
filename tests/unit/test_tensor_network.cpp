//
// Created by Eduard Valeyev on 3/23/18.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick_graph.hpp>

#include <algorithm>
#include <cmath>
#include <codecvt>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <SeQuant/core/timer.hpp>
#include <range/v3/all.hpp>

std::string to_utf8(const std::wstring& wstr) {
  using convert_type = std::codecvt_utf8<wchar_t>;
  std::wstring_convert<convert_type, wchar_t> converter;
  return converter.to_bytes(wstr);
}

template <typename Container>
std::vector<sequant::ExprPtr> to_tensors(const Container& cont) {
  std::vector<sequant::ExprPtr> tensors;

  std::transform(cont.begin(), cont.end(), std::back_inserter(tensors),
                 [](const auto& tensor) {
                   auto casted =
                       std::dynamic_pointer_cast<sequant::Expr>(tensor);
                   REQUIRE(casted != nullptr);
                   return casted;
                 });
  return tensors;
}

template <typename Container>
sequant::ExprPtr to_product(const Container& container) {
  return sequant::ex<sequant::Product>(to_tensors(container));
}

namespace sequant {
class TensorNetworkAccessor {
 public:
  auto get_canonical_bliss_graph(
      sequant::TensorNetwork tn,
      const sequant::TensorNetwork::named_indices_t* named_indices = nullptr) {
    tn.canonicalize_graph(named_indices ? *named_indices : tn.ext_indices_);
    tn.init_edges();
    auto graph = tn.create_graph(named_indices);
    return std::make_pair(std::move(graph.bliss_graph), graph.vertex_labels);
  }
};
}  // namespace sequant

TEST_CASE("TensorNetwork", "[elements]") {
  using namespace sequant;
  using namespace sequant::mbpt::tensor;

  sequant::set_default_context(Context(
      mbpt::make_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, SPBasis::spinorbital));

  SECTION("Edges") {
    using Vertex = TensorNetwork::Vertex;
    using Edge = TensorNetwork::Edge;
    using Origin = TensorNetwork::Origin;

    Vertex v1(Origin::Bra, 0, 1, Symmetry::antisymm);
    Vertex v2(Origin::Bra, 0, 0, Symmetry::antisymm);
    Vertex v3(Origin::Ket, 1, 0, Symmetry::symm);
    Vertex v4(Origin::Ket, 1, 3, Symmetry::symm);
    Vertex v5(Origin::Bra, 3, 0, Symmetry::nonsymm);
    Vertex v6(Origin::Bra, 3, 2, Symmetry::nonsymm);
    Vertex v7(Origin::Ket, 3, 1, Symmetry::nonsymm);
    Vertex v8(Origin::Ket, 5, 0, Symmetry::symm);

    const Index dummy(L"a_1");

    Edge e1(v1, dummy);
    e1.connect_to(v4);
    Edge e2(v2, dummy);
    e2.connect_to(v3);
    Edge e3(v3, dummy);
    e3.connect_to(v5);
    Edge e4(v4, dummy);
    e4.connect_to(v6);

    Edge e5(v8, dummy);
    e5.connect_to(v6);
    Edge e6(v8, dummy);
    e6.connect_to(v7);

    Edge e7(v4, dummy);

    // Due to tensor symmetries, these edges are considered equal
    REQUIRE(e1 == e2);
    REQUIRE(!(e1 < e2));
    REQUIRE(!(e2 < e1));

    // Smallest terminal index wins
    REQUIRE(!(e1 == e3));
    REQUIRE(e1 < e3);
    REQUIRE(!(e3 < e1));

    // For non-symmetric tensors the connection slot is taken into account
    REQUIRE(!(e3 == e4));
    REQUIRE(e3 < e4);
    REQUIRE(!(e4 < e3));

    REQUIRE(!(e5 == e6));
    REQUIRE(e5 < e6);
    REQUIRE(!(e6 < e5));

    // Unconnected edges always come before fully connected ones
    REQUIRE(!(e7 == e1));
    REQUIRE(e7 < e1);
    REQUIRE(!(e1 < e7));
  }

  SECTION("constructors") {
    {  // with Tensors
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t2 =
          ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));

      auto t1_x_t2_p_t2 = t1 * (t2 + t2);  // can only use a flat tensor product
      REQUIRE_THROWS_AS(TensorNetwork(*t1_x_t2_p_t2), std::logic_error);
    }

    {  // with NormalOperators
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_2"}, V);
      auto t2 = ex<FNOperator>(WstrList{L"i_2"}, WstrList{L"i_1"}, V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));
    }

    {  // with Tensors and NormalOperators
      auto tmp = A(-2) * H_(2) * T_(2) * T_(2);
      REQUIRE_NOTHROW(TensorNetwork(tmp->as<Product>().factors()));
    }

  }  // SECTION("constructors")

  SECTION("accessors") {
    {
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t2 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_3"}, V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));
      TensorNetwork tn(*t1_x_t2);

      // edges
      auto edges = tn.edges();
      REQUIRE(edges.size() == 3);

      // ext indices
      auto ext_indices = tn.ext_indices();
      REQUIRE(ext_indices.size() == 2);

      // tensors
      auto tensors = tn.tensors();
      REQUIRE(size(tensors) == 2);
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[0]));
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[1]));
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[0]) == *t1);
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[1]) == *t2);
    }
  }  // SECTION("accessors")

  SECTION("canonicalizer") {
    SECTION("no externals") {
      Index::reset_tmp_index();
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t2 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_2"}, V);
      auto t1_x_t2 = t1 * t2;
      TensorNetwork tn(*t1_x_t2);
      tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

      REQUIRE(size(tn.tensors()) == 2);
      REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
      REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
      //        std::wcout <<
      //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
      //        std::endl; std::wcout <<
      //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
      //        std::endl;
      REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
              L"{F^{{i_2}}_{{i_1}}}");
      REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
              L"{\\tilde{a}^{{i_1}}_{{i_2}}}");
    }

    SECTION("with externals") {
      Index::reset_tmp_index();
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_2"}, WstrList{L"i_17"}, WstrList{});
      auto t2 = ex<FNOperator>(WstrList{L"i_2"}, WstrList{L"i_3"}, V);
      auto t1_x_t2 = t1 * t2;

      // with all external named indices
      SECTION("implicit") {
        TensorNetwork tn(*t1_x_t2);
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

        REQUIRE(size(tn.tensors()) == 2);
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
        // std::wcout <<
        // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
        // std::endl; std::wcout <<
        // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
        // std::endl;
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_1}}_{{i_3}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                L"{F^{{i_{17}}}_{{i_1}}}");
      }

      // with explicit named indices
      SECTION("explicit") {
        Index::reset_tmp_index();
        TensorNetwork tn(*t1_x_t2);

        using named_indices_t = TensorNetwork::named_indices_t;
        named_indices_t indices{Index{L"i_17"}};
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                        &indices);

        REQUIRE(size(tn.tensors()) == 2);
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
        //        std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]))
        //        << std::endl; std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]))
        //        << std::endl;
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_2}}_{{i_1}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                L"{F^{{i_{17}}}_{{i_2}}}");
      }
    }

    SECTION("particle non-conserving") {
      const auto input1 = parse_expr(L"P{;a1,a3}");
      const auto input2 = parse_expr(L"P{a1,a3;}");
      const std::wstring expected1 = L"{{P^{{a_1}{a_3}}_{}}}";
      const std::wstring expected2 = L"{{P^{}_{{a_1}{a_3}}}}";

      for (int variant : {1, 2}) {
        for (bool fast : {true, false}) {
          TensorNetwork tn(
              std::vector<ExprPtr>{variant == 1 ? input1 : input2});
          tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), fast);
          REQUIRE(tn.tensors().size() == 1);
          auto result = ex<Product>(to_tensors(tn.tensors()));
          REQUIRE(to_latex(result) == (variant == 1 ? expected1 : expected2));
        }
      }
    }

    SECTION("non-symmetric") {
      const auto input =
          parse_expr(L"A{i7,i3;i9,i12}:A I1{i7,i3;;x5}:N I2{;i9,i12;x5}:N")
              .as<Product>()
              .factors();
      const std::wstring expected =
          L"A{i_1,i_2;i_3,i_4}:A * I1{i_1,i_2;;x_1}:N * I2{;i_3,i_4;x_1}:N";

      for (bool fast : {true, false}) {
        TensorNetwork tn(input);
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), fast);
        const auto result = ex<Product>(to_tensors(tn.tensors()));
        REQUIRE(deparse_expr(result) == expected);
      }
    }

    SECTION("particle-1,2-symmetry") {
      const std::vector<std::pair<std::wstring, std::wstring>> pairs = {
          {L"S{i_1,i_2,i_3;a_1,a_2,a_3}:N * f{i_4;i_2}:N * "
           L"t{a_1,a_2,a_3;i_4,i_3,i_1}:N",
           L"S{i_1,i_2,i_3;a_1,a_2,a_3}:N * f{i_4;i_1}:N * "
           L"t{a_1,a_2,a_3;i_2,i_3,i_4}:N"},
          {L"Γ{o_2,o_4;o_1,o_3}:N * g{i_1,o_1;o_2,e_1}:N * "
           L"t{o_3,e_1;o_4,i_1}:N",
           L"Γ{o_2,o_4;o_1,o_3}:N * g{i_1,o_3;o_4,e_1}:N * "
           L"t{o_1,e_1;o_2,i_1}:N"}};
      for (const auto& pair : pairs) {
        const auto first = parse_expr(pair.first).as<Product>().factors();
        const auto second = parse_expr(pair.second).as<Product>().factors();

        TensorNetworkAccessor accessor;
        auto [first_graph, first_labels] =
            accessor.get_canonical_bliss_graph(TensorNetwork(first));
        auto [second_graph, second_labels] =
            accessor.get_canonical_bliss_graph(TensorNetwork(second));
        if (first_graph->cmp(*second_graph) != 0) {
          std::wstringstream stream;
          stream << "First graph:\n";
          first_graph->write_dot(stream, first_labels, true);
          stream << "Second graph:\n";
          second_graph->write_dot(stream, second_labels, true);
          stream << "Wick graph:\n";
          auto [wick_graph, labels, d1, d2] =
              WickGraph(first).make_bliss_graph();
          wick_graph->write_dot(stream, labels, true);

          FAIL(to_utf8(stream.str()));
        }

        TensorNetwork tn1(first);
        TensorNetwork tn2(second);

        tn1.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);
        tn2.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

        REQUIRE(tn1.tensors().size() == tn2.tensors().size());
        for (std::size_t i = 0; i < tn1.tensors().size(); ++i) {
          auto t1 = std::dynamic_pointer_cast<Expr>(tn1.tensors()[i]);
          auto t2 = std::dynamic_pointer_cast<Expr>(tn2.tensors()[i]);
          REQUIRE(t1);
          REQUIRE(t2);
          REQUIRE(to_latex(t1) == to_latex(t2));
        }
      }
    }

    SECTION("miscellaneous") {
      const std::vector<std::pair<std::wstring, std::wstring>> inputs = {
          {L"g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A",
           L"g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A"},
          {L"g{a_1,i_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A",
           L"-1 g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A"},

          {L"g{i_1,a_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N",
           L"g{i_1,a_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N"},
          {L"g{a_1,i_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N",
           L"g{i_1,a_1;i_2,i_3}:N * I{i_3,i_2;i_1,a_1}:N"},
      };

      for (const auto& [input, expected] : inputs) {
        const auto input_tensors = parse_expr(input).as<Product>().factors();

        TensorNetwork tn(input_tensors);
        ExprPtr factor = tn.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), true);

        ExprPtr prod = to_product(tn.tensors());
        if (factor) {
          prod = ex<Product>(
              prod.as<Product>().scale(factor.as<Constant>().value()));
        }

        REQUIRE(deparse_expr(prod) == expected);
      }
    }

#ifndef SEQUANT_SKIP_LONG_TESTS
    SECTION("Exhaustive SRCC example") {
      // Note: the exact canonical form written here is implementation-defined
      // and doesn't actually matter What does, is that all equivalent ways of
      // writing it down, canonicalizes to the same exact form
      const Product expectedExpr =
          parse_expr(
              L"A{i1,i2;a1,a2} g{i3,i4;a3,a4} t{a1,a3;i3,i4} t{a2,a4;i1,i2}",
              Symmetry::antisymm)
              .as<Product>();

      const auto expected = expectedExpr.factors();

      TensorNetworkAccessor accessor;
      const auto [canonical_graph, canonical_graph_labels] =
          accessor.get_canonical_bliss_graph(TensorNetwork(expected));

      std::wcout << "Canonical graph:\n";
      canonical_graph->write_dot(std::wcout, canonical_graph_labels);
      std::wcout << std::endl;

      std::vector<Index> indices;
      for (std::size_t i = 0; i < expected.size(); ++i) {
        const Tensor& tensor = expected[i].as<Tensor>();
        for (const Index& idx : tensor.indices()) {
          if (std::find(indices.begin(), indices.end(), idx) == indices.end()) {
            indices.push_back(idx);
          }
        }
      }
      std::sort(indices.begin(), indices.end());

      const auto original_indices = indices;

      // Make sure to clone all expressions in order to not accidentally modify
      // the ones in expected (even though they are const... the pointer-like
      // semantics of expressions messes with const semantics)
      std::remove_const_t<decltype(expected)> factors;
      for (const auto& factor : expected) {
        factors.push_back(factor.clone());
      }
      std::sort(factors.begin(), factors.end());

      const auto is_occ = [](const Index& idx) {
        return idx.space() == Index(L"i_1").space();
      };

      // Iterate over all tensor permutations and all permutations of possible
      // index name swaps
      REQUIRE(std::is_sorted(factors.begin(), factors.end()));
      REQUIRE(std::is_sorted(indices.begin(), indices.end()));
      REQUIRE(std::is_partitioned(indices.begin(), indices.end(), is_occ));
      REQUIRE(std::partition_point(indices.begin(), indices.end(), is_occ) ==
              indices.begin() + 4);
      std::size_t total_variations = 0;
      do {
        do {
          do {
            total_variations++;

            // Compute index replacements
            container::map<Index, Index> idxrepl;
            for (std::size_t i = 0; i < indices.size(); ++i) {
              REQUIRE(original_indices[i].space() == indices[i].space());

              idxrepl.insert(
                  std::make_pair(original_indices.at(i), indices.at(i)));
            }

            // Apply index replacements to a copy of the current tensor
            // permutation
            auto copy = factors;
            for (ExprPtr& expr : copy) {
              expr.as<Tensor>().transform_indices(idxrepl);
              reset_tags(expr.as<Tensor>());
            }

            TensorNetwork tn(copy);

            // At the heart of our canonicalization lies the fact that we can
            // always create the uniquely defined canonical graph for a given
            // network
            const auto [current_graph, current_graph_labels] =
                accessor.get_canonical_bliss_graph(tn);
            if (current_graph->cmp(*canonical_graph) != 0) {
              std::wcout << "Canonical graph for "
                         << deparse_expr(ex<Product>(copy)) << ":\n";
              current_graph->write_dot(std::wcout, current_graph_labels);
              std::wcout << std::endl;
            }
            REQUIRE(current_graph->cmp(*canonical_graph) == 0);

            tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                            false);

            std::vector<ExprPtr> actual;
            std::transform(tn.tensors().begin(), tn.tensors().end(),
                           std::back_inserter(actual), [](const auto& t) {
                             assert(std::dynamic_pointer_cast<Expr>(t));
                             return std::dynamic_pointer_cast<Expr>(t);
                           });

            // The canonical graph must not change due to the other
            // canonicalization steps we perform
            REQUIRE(accessor.get_canonical_bliss_graph(TensorNetwork(actual))
                        .first->cmp(*canonical_graph) == 0);

            REQUIRE(actual.size() == expected.size());

            if (!std::equal(expected.begin(), expected.end(), actual.begin())) {
              std::wostringstream sstream;
              sstream
                  << "Expected all tensors to be equal (actual == expected), "
                     "but got:\n";
              for (std::size_t i = 0; i < expected.size(); ++i) {
                std::wstring equality =
                    actual[i] == expected[i] ? L" == " : L" != ";

                sstream << deparse_expr(actual[i]) << equality
                        << deparse_expr(expected[i]) << "\n";
              }
              sstream << "\nInput was " << deparse_expr(ex<Product>(factors))
                      << "\n";
              FAIL(to_utf8(sstream.str()));
            }
          } while (std::next_permutation(indices.begin() + 4, indices.end()));
        } while (std::next_permutation(indices.begin(), indices.begin() + 4));
      } while (std::next_permutation(factors.begin(), factors.end()));

      // 4! (tensors) * 4! (internal indices) * 4! (external indices)
      REQUIRE(total_variations == 24 * 24 * 24);
    }
#endif

    SECTION("idempotency") {
      const std::vector<std::wstring> inputs = {
          L"F{i1;i8} g{i8,i9;i1,i7}",
          L"A{i7,i3;i9,i12}:A I1{i7,i3;;x5}:N I2{;i9,i12;x5}:N",
          L"f{i4;i1}:N t{a1,a2,a3;i2,i3,i4}:N S{i1,i2,i3;a1,a2,a3}:N",
          L"P{a1,a3;} k{i8;i2}",
          L"L{x6;;x2} P{;a1,a3}",
      };

      for (const std::wstring& current : inputs) {
        auto factors1 = parse_expr(current).as<Product>().factors();
        auto factors2 = parse_expr(current).as<Product>().factors();

        TensorNetwork reference_tn(factors1);
        reference_tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                                  false);

        TensorNetwork check_tn(factors2);
        check_tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                              false);

        REQUIRE(to_latex(to_product(reference_tn.tensors())) ==
                to_latex(to_product(check_tn.tensors())));

        for (bool fast : {true, false, true, true, false, false, true}) {
          reference_tn.canonicalize(
              TensorCanonicalizer::cardinal_tensor_labels(), fast);

          REQUIRE(to_latex(to_product(reference_tn.tensors())) ==
                  to_latex(to_product(check_tn.tensors())));
        }
      }
    }
  }  // SECTION("canonicalizer")

  SECTION("misc1") {
    if (false) {
      Index::reset_tmp_index();
      // TN1 from manuscript
      auto g =
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},
                     WstrList{}, Symmetry::antisymm);
      auto ta =
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_3"}, WstrList{L"i_1", L"i_2"},
                     WstrList{}, Symmetry::antisymm);
      auto tb =
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_4"}, WstrList{L"i_3", L"i_4"},
                     WstrList{}, Symmetry::antisymm);

      auto tmp = g * ta * tb;
      // std::wcout << "TN1 = " << to_latex(tmp) << std::endl;
      TensorNetwork tn(tmp->as<Product>().factors());

      // make graph
      // N.B. treat all indices as dummy so that the automorphism ignores the
      using named_indices_t = TensorNetwork::named_indices_t;
      named_indices_t indices{};
      REQUIRE_NOTHROW(tn.create_graph(&indices));
      TensorNetwork::Graph graph = tn.create_graph(&indices);

      // create dot
      {
        std::basic_ostringstream<wchar_t> oss;
        REQUIRE_NOTHROW(graph.bliss_graph->write_dot(oss, graph.vertex_labels));
        // std::wcout << "oss.str() = " << std::endl << oss.str() << std::endl;
      }

      bliss::Stats stats;
      graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

      std::vector<std::vector<unsigned int>> aut_generators;
      auto save_aut = [&aut_generators](const unsigned int n,
                                        const unsigned int* aut) {
        aut_generators.emplace_back(aut, aut + n);
      };
      graph.bliss_graph->find_automorphisms(
          stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);
      CHECK(aut_generators.size() ==
            2);  // there are 2 generators, i1<->i2, i3<->i4

      std::basic_ostringstream<wchar_t> oss;
      bliss::print_auts(aut_generators, oss, graph.vertex_labels);
      CHECK(oss.str() == L"({i_3},{i_4})\n({i_1},{i_2})\n");
      // std::wcout << oss.str() << std::endl;
    }

    // profile canonicalization for synthetic tests in
    // DOI 10.1016/j.cpc.2018.02.014
    if (false) {
      for (auto testcase : {0, 1, 2, 3}) {
        // - testcase=0,2 are "equivalent" and correspond to the "frustrated"
        //   case in Section 5.3 of DOI 10.1016/j.cpc.2018.02.014
        // - testcase=1 corresponds to the "frustrated" case in Section 5.4 of
        //   DOI 10.1016/j.cpc.2018.02.014
        // - testcase=3 corresponds to the "No symmetry dummy"
        //   case in Section 5.1 of DOI 10.1016/j.cpc.2018.02.014
        if (testcase == 0)
          std::wcout
              << "canonicalizing network with 1 totally-symmetric tensor with "
                 "N indices and 1 asymmetric tensor with N indices\n";
        else if (testcase == 3)
          std::wcout << "canonicalizing network with 1 asymmetric tensor with "
                        "N indices and 1 asymmetric tensor with N indices\n";
        else
          std::wcout << "canonicalizing network with n equivalent asymmetric "
                        "tensors with N/n indices each and 1 asymmetric tensor "
                        "with N indices\n";

        std::wcout << "N,n,min_time,geommean_time,max_time\n";

        for (auto N :
             {1, 2, 4, 8, 16, 32, 64, 128, 256}) {  // total number of indices

          int n;
          switch (testcase) {
            case 0:
              n = 1;
              break;
            case 1:
              n = N / 2;
              break;
            case 2:
              n = N;
              break;
            case 3:
              n = 1;
              break;
            default:
              abort();
          }
          if (n == 0 || n > N) continue;

          auto ctx_resetter = set_scoped_default_context(
              (static_cast<std::size_t>(N) > Index::min_tmp_index())
                  ? Context(get_default_context())
                        .set_first_dummy_index_ordinal(N + 1)
                  : get_default_context());

          // make list of indices
          std::vector<Index> indices;
          for (auto i = 0; i != N; ++i) {
            std::wostringstream oss;
            oss << "i_" << i;
            indices.emplace_back(oss.str());
          }
          std::random_device rd;

          // randomly sample connectivity between bra and ket tensors
          const auto S = 10;  // how many samples to take

          auto product_time =
              1.;  // product of all times, need to get geometric mean
          auto min_time =
              std::numeric_limits<double>::max();  // total time for all samples
          auto max_time =
              std::numeric_limits<double>::min();  // total time for all samples
          for (auto s = 0; s != S; ++s) {
            // make tensors of independently (and randomly) permuted
            // contravariant and covariant indices
            auto contravariant_indices = indices;
            auto covariant_indices = indices;

            std::shuffle(contravariant_indices.begin(),
                         contravariant_indices.end(), std::mt19937{rd()});
            std::shuffle(covariant_indices.begin(), covariant_indices.end(),
                         std::mt19937{rd()});

            auto utensors =
                covariant_indices | ranges::views::chunk(N / n) |
                ranges::views::transform([&](const auto& idxs) {
                  return ex<Tensor>(
                      L"u", idxs, std::vector<Index>{}, std::vector<Index>{},
                      (testcase == 3
                           ? Symmetry::nonsymm
                           : ((n == 1) ? Symmetry::symm : Symmetry::nonsymm)));
                }) |
                ranges::to_vector;
            CHECK(utensors.size() == static_cast<std::size_t>(n));
            auto dtensors = contravariant_indices | ranges::views::chunk(N) |
                            ranges::views::transform([&](const auto& idxs) {
                              return ex<Tensor>(L"d", std::vector<Index>{},
                                                std::vector<Index>{}, idxs,
                                                Symmetry::nonsymm);
                            }) |
                            ranges::to_vector;
            CHECK(dtensors.size() == 1);

            ExprPtr expr;
            for (auto g = 0; g != n; ++g) {
              if (g == 0)
                expr = utensors[0] * dtensors[0];
              else
                expr = expr * utensors[g];
            }

            TensorNetwork tn(expr->as<Product>().factors());

            // produce misc data for publication
            if (false && s == 0) {
              std::wcout << "N=" << N << " n=" << n << " expr:\n"
                         << expr->to_latex() << std::endl;

              // make graph
              REQUIRE_NOTHROW(tn.create_graph());
              TensorNetwork::Graph graph = tn.create_graph();

              // create dot
              std::basic_ostringstream<wchar_t> oss;
              REQUIRE_NOTHROW(
                  graph.bliss_graph->write_dot(oss, graph.vertex_labels));
              std::wcout << "bliss graph:" << std::endl
                         << oss.str() << std::endl;
            }

            sequant::TimerPool<> timer;
            timer.start();
            tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                            false);
            timer.stop();
            const auto elapsed_seconds = timer.read();
            product_time *= elapsed_seconds;
            min_time = std::min(min_time, elapsed_seconds);
            max_time = std::max(max_time, elapsed_seconds);
          }

          const auto geommean_time = std::pow(product_time, 1. / S);
          std::wcout << N << "," << n << "," << min_time << "," << geommean_time
                     << "," << max_time << "\n";
        }
      }
    }

  }  // SECTION("misc1")

}  // TEST_CASE("Tensor")
