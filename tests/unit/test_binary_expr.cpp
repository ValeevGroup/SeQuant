#include <SeQuant/domain/utils/binary_expr.hpp>
#include <iostream>
#include <sstream>

#include "catch.hpp"

using sequant::utils::binary_expr;

template <typename T, typename F>
void visit_preorder_binary_expr(typename binary_expr<T>::node_ptr const& node,
                                F&& visitor) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "visitor signature not matched");
  if (!node) return;  // if nullptr passed explicitly

  if (node->leaf()) {
    visitor(node);
    return;
  }

  // visit this node
  visitor(node);

  // visit left node
  visit_preorder_binary_expr<T, F>(node->left(), std::forward<F>(visitor));

  // visit right node
  visit_preorder_binary_expr<T, F>(node->right(), std::forward<F>(visitor));
}

template <typename T, typename F>
void visit_postorder_binary_expr(typename binary_expr<T>::node_ptr const& node,
                                 F&& visitor) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "visitor signature not matched");

  if (node->leaf()) {
    visitor(node);
    return;
  }
  // visit left node
  visit_postorder_binary_expr<T, F>(node->left(), std::forward<F>(visitor));

  // visit right node
  visit_postorder_binary_expr<T, F>(node->right(), std::forward<F>(visitor));

  // visit this node
  visitor(node);
}

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
    auto node1 = make_binary_expr(0, make_binary_expr(1), make_binary_expr(2));
    REQUIRE(node1->data() == 0);
    REQUIRE(node1->left()->data() == 1);
    REQUIRE(node1->right()->data() == 2);
  }

  SECTION("visitor methods") {
    using namespace std::literals::string_literals;

    using node_type = binary_expr<std::string>;

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
      void operator()(node_type::node_ptr const& node) {
        words.emplace_back(node->data());
      }

    } visitor;

    visitor.words.clear();
    sequant::utils::visit_inorder_binary_expr<std::string>(node1, visitor);
    REQUIRE(visitor.words == std::vector<std::string>{"bar", "foo", "bazz"});

    visitor.words.clear();
    visit_preorder_binary_expr<std::string>(node1, visitor);
    REQUIRE(visitor.words == std::vector<std::string>{"foo", "bar", "bazz"});

    visitor.words.clear();
    visit_postorder_binary_expr<std::string>(node1, visitor);
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
      double operator()(typename binary_expr<int>::node_ptr const& node) {
        return node->data();
      }

      // sums left and right result then divides by this node's data
      double operator()(typename binary_expr<int>::node_ptr const& node,
                        double leval, double reval) {
        return (leval + reval) / node->data();
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

    REQUIRE(sequant::utils::evaluate_binary_expr<int>(node, evaluator) ==
            Approx(7.4 / (1.1 * 2)));
  }

  SECTION("digraph") {
    const auto& node = make_binary_expr(1, 2, 3);
    std::stringstream os;
    sequant::utils::digraph_binary_expr<int>(
        os, node,
        [](binary_expr<int>::node_ptr const& x) { return x->data(); });

    auto expected =
        "digraph binary_expr {\n"
        "node0[label=1];\n"
        "node0 -> node1;\n"
        "node1[label=2];\n"
        "node0 -> node2;\n"
        "node2[label=3];\n"
        "}\n";

    REQUIRE(os.str() == expected);
  }
}
