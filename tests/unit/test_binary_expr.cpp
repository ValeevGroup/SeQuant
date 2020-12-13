#include <SeQuant/domain/utils/binary_expr.hpp>
#include <iostream>

#include "catch.hpp"

TEST_CASE("TEST BINARY_EXPR_CLASS", "[binary_expr]") {
  using sequant::utils::binary_expr;
  using sequant::utils::make_binary_expr;
  using sequant::utils::detail::binary_expr_internal;
  using sequant::utils::detail::binary_expr_leaf;

  SECTION("binary_expr constructors") {
    using namespace std::string_literals;
    REQUIRE_NOTHROW(binary_expr_internal<size_t>(0, nullptr, nullptr));

    auto node1 = binary_expr_leaf(0);
    REQUIRE(node1.leaf());

    auto node2 = make_binary_expr(0, make_binary_expr(1), make_binary_expr(2));

    REQUIRE(node2->left() != nullptr);
    REQUIRE(node2->right() != nullptr);
    REQUIRE_FALSE(node2->leaf());

    REQUIRE_NOTHROW([&node2] { auto node3 = std::move(node2); });

    REQUIRE_NOTHROW(make_binary_expr(1, 2, 3));

    auto node3 = make_binary_expr("P"s, "L"s, "R"s);
    REQUIRE(node3->data() == "P"s);
    REQUIRE(node3->left()->data() == "L"s);
    REQUIRE(node3->right()->data() == "R"s);
  }

  SECTION("methods") {
    using data_type = size_t;
    auto node1 = make_binary_expr(0, make_binary_expr(1), make_binary_expr(2));
    REQUIRE(node1->data() == 0);
    REQUIRE(node1->left()->data() == 1);
    REQUIRE(node1->right()->data() == 2);
  }

  SECTION("visitor methods") {
    using namespace std::literals::string_literals;
    auto node1 = make_binary_expr("foo"s,  //
                                  "bar"s,  //
                                  "bazz"s);
    // node1 tree:
    //
    //       foo
    //      /   \
    //   bar    bazz
    //

    struct {
      std::vector<std::string> words;
      void operator()(const std::string& str) { words.emplace_back(str); }

    } visitor;

    sequant::utils::visit_inorder_binary_expr<std::string>(node1, std::ref(visitor));
    REQUIRE(visitor.words == std::vector<std::string>{"bar", "foo", "bazz"});

    visitor.words.clear();
    sequant::utils::visit_preorder_binary_expr<std::string>(node1, std::ref(visitor));
    REQUIRE(visitor.words == std::vector<std::string>{"foo", "bar", "bazz"});

    visitor.words.clear();
    sequant::utils::visit_postorder_binary_expr<std::string>(node1, std::ref(visitor));
    REQUIRE(visitor.words == std::vector<std::string>{"bar", "bazz", "foo"});
  }

  SECTION("evaluator") {
    auto node = make_binary_expr(2, make_binary_expr(3, 5, 7),
                                 make_binary_expr(11, 13, 17));
    // node:
    //
    //             2
    //          /     \
    //         3       11
    //       /  \      / \
    //      5    7   13  17

    //
    // evaluator
    //
    struct {
      double operator()(int x) { return x; }

      // sums left and right result then divides by this node's data
      double operator()(typename binary_expr<int>::node_ptr const& n, double l,
                        double r) {
        return (l + r) / n->data();
      }

    } evaluator;

    //
    // expected
    //                 7.4/(1.1 * 2)
    //              /       \
    //            /          \
    //          /             \
    //         4.0            3/1.1
    //       /  \            / \
    //      5.0    7.0   13.0  17.0
    //

    REQUIRE(sequant::utils::evaluate_binary_expr<int, double>(
                node, evaluator) == Approx(7.4 / (1.1 * 2)));
  }

  SECTION("transformation") {
    using namespace std::literals::string_literals;
    auto node1 = make_binary_expr("foo"s, "bar"s, "bazz"s);
    // node1 tree:
    //
    //       foo
    //      /   \
    //   bar    bazz
    //

    struct {
      size_t operator()(const std::string& str) { return str.size(); }
    } transformer;

    auto node1_trans =
        sequant::utils::transform_binary_expr<std::string, size_t>(node1,
                                                                   transformer);

    // node1_trans tree:
    //
    //        3
    //      /   \
    //     3     4

    struct {
      std::vector<size_t> lengths;
      void operator()(size_t x) { lengths.emplace_back(x); }
    } visitor;

    sequant::utils::visit_inorder_binary_expr<size_t>(node1_trans, visitor);
    REQUIRE(visitor.lengths == std::vector<size_t>{3, 3, 4});
  }
}
