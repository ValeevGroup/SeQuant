#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/binary_node.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <utility>

#include <range/v3/all.hpp>

TEST_CASE("TEST BINARY_NODE", "[FullBinaryNode]") {
  using ranges::views::iota;
  using ranges::views::take;
  using sequant::FullBinaryNode;
  using sequant::TreeTraversal;

  SECTION("construction") {
    REQUIRE_NOTHROW(FullBinaryNode{0});
    REQUIRE_NOTHROW(FullBinaryNode{'a', 'b', 'c'});
    REQUIRE_NOTHROW(
        FullBinaryNode{'a', FullBinaryNode{'b'}, FullBinaryNode{'c'}});
  }

  SECTION("copy ctor and assign") {
    auto const n1 = FullBinaryNode{'a', 'b', 'c'};
    auto const n2{n1};
    REQUIRE(&n1 != &n2);
    REQUIRE(n1 == n2);

    auto const n3 = n1;
    REQUIRE(&n1 != &n3);
    REQUIRE(n1 == n3);
  }

  SECTION("move ctor and assign") {
    auto n1{FullBinaryNode{1, 2, 3}};
    auto n2 = std::move(n1);
    REQUIRE(n2 == FullBinaryNode{1, 2, 3});
  }

  SECTION("derefence") {
    auto const n1 = FullBinaryNode{100};

    REQUIRE_NOTHROW(*n1);
    REQUIRE(*n1 == 100);

    struct dummy {
      void dummy_fun() const {}
    };

    auto const n2 = FullBinaryNode{dummy{}};

    REQUIRE_NOTHROW(n2->dummy_fun());
  }

  SECTION("internal node") {
    auto const n = FullBinaryNode{3, 2, 5};

    REQUIRE_FALSE(n.leaf());
    REQUIRE(*n == 3);
    REQUIRE(*n.left() == 2);
    REQUIRE(*n.right() == 5);
  }

  SECTION("leaf node") {
    auto const n = FullBinaryNode{'n'};
    REQUIRE(*n == 'n');
    REQUIRE(n.leaf());

    REQUIRE_THROWS(n.left());
    REQUIRE_THROWS(n.right());
  }

  SECTION("advanced construction") {
    auto constexpr leaves = std::array{3, 2, 4};

    struct adder {
      int operator()(int x) const { return x; }
      int operator()(int x, int y) const { return x + y; }
    };

    auto const node = FullBinaryNode<int>{leaves, adder{}};

    REQUIRE(*node == 9);
    REQUIRE(*node.left() == 5);
    REQUIRE(*node.right() == 4);
    REQUIRE(*node.left().left() == 3);
    REQUIRE(*node.left().right() == 2);

    auto const leaves2 = ranges::views::all(leaves);
    auto const node2 = FullBinaryNode<int>{leaves2, adder{}};
  }

  SECTION("reconnect") {
    auto subtree1 = FullBinaryNode<int>{3, 2, 1};
    auto subtree2 = FullBinaryNode<int>{7, subtree1, FullBinaryNode<int>(4)};
    auto tree = FullBinaryNode<int>{13, subtree2, FullBinaryNode<int>(5)};

    REQUIRE(*tree == 13);
    REQUIRE(*tree.left() == 7);
    REQUIRE(*tree.left().left() == 3);
    REQUIRE(*tree.left().left().left() == 2);
    REQUIRE(*tree.left().left().right() == 1);
    REQUIRE(*tree.left().right() == 4);
    REQUIRE(*tree.right() == 5);

    tree.left() = std::move(tree.left().left());

    REQUIRE(*tree == 13);
    REQUIRE(*tree.left() == 3);
    REQUIRE(*tree.left().left() == 2);
    REQUIRE(*tree.left().right() == 1);
    REQUIRE(*tree.right() == 5);
    REQUIRE(!tree.left().leaf());
    REQUIRE(tree.left().left().leaf());
    REQUIRE(tree.left().right().leaf());
    REQUIRE(tree.right().leaf());
  }

  SECTION("digraph generation") {
    auto take_nums = [](size_t count, int from = 1) {
      return iota(from) | take(count);
    };

    struct make_sum {
      int operator()(int x) const { return x; }
      int operator()(int x, int y) const { return x + y; }
    };

    auto ms = make_sum{};
    auto const node1 = FullBinaryNode<int>{take_nums(1), ms};
    auto const node2 = FullBinaryNode<int>{take_nums(2), ms};
    auto const node3 = FullBinaryNode<int>{take_nums(3), ms};
    auto const node4 = FullBinaryNode<int>{
        6, FullBinaryNode<int>{1}, FullBinaryNode<int>{take_nums(2, 2), ms}};

    auto label_gen_str = [](auto const& n) { return std::to_string(*n); };
    auto label_gen_wstr = [](auto const& n) { return std::to_wstring(*n); };

    REQUIRE(node1.digraph<std::string>(label_gen_str, "node1") ==
            std::string{"digraph node1{\n"
                        "node0[label=1];\n"
                        "}"});

    REQUIRE(node2.digraph<std::string>(label_gen_str) ==
            std::string{"digraph {\n"
                        "node0[label=3];\n"
                        "node1[label=1];\n"
                        "node2[label=2];\n"
                        "node0 -> node1;\n"
                        "node0 -> node2;\n"
                        "}"});

    REQUIRE(node3.digraph<std::wstring>(label_gen_wstr, L"node3") ==
            std::wstring{L"digraph node3{\n"
                         "node0[label=6];\n"
                         "node1[label=3];\n"
                         "node2[label=1];\n"
                         "node3[label=2];\n"
                         "node1 -> node2;\n"
                         "node1 -> node3;\n"
                         "node4[label=3];\n"
                         "node0 -> node1;\n"
                         "node0 -> node4;\n"
                         "}"});

    REQUIRE(node4.digraph<std::wstring>(label_gen_wstr, L"node4") ==
            std::wstring{L"digraph node4{\n"
                         "node0[label=6];\n"
                         "node1[label=1];\n"
                         "node2[label=5];\n"
                         "node3[label=2];\n"
                         "node4[label=3];\n"
                         "node2 -> node3;\n"
                         "node2 -> node4;\n"
                         "node0 -> node1;\n"
                         "node0 -> node2;\n"
                         "}"});

    /*
       tree from node1:
                    1

       tree from node2:
                    3
                  /   \
                1      2

       tree from node3:
                    6
                  /   \
                3       3
              /   \
            1      2

       tree from node4:
                    6
                  /   \
                1      5
                      /  \
                    2     3
    */
  }

  SECTION("visitor") {
    using node_t = FullBinaryNode<std::string>;
    /*
                 C
               /   \
             A       B
    */
    auto const node = node_t{"C", "A", "B"};
    REQUIRE(*node == "C");
    REQUIRE(*node.left() == "A");
    REQUIRE(*node.right() == "B");

    std::string str1;
    std::string str2;
    auto visitor = [&str1](node_t const& n) { str1 += *n; };
    SECTION("single order") {
      auto extended_visitor = [&str2](node_t const& n, TreeTraversal) {
        str2 += *n;
      };

      node.visit(visitor, TreeTraversal::PostOrder);
      node.visit(extended_visitor, TreeTraversal::PostOrder);
      REQUIRE(str1 == "ABC");
      REQUIRE(str1 == str2);

      str1.clear();
      str2.clear();
      node.visit(visitor, TreeTraversal::PreOrder);
      node.visit(extended_visitor, TreeTraversal::PreOrder);
      REQUIRE(str1 == "CAB");
      REQUIRE(str1 == str2);

      str1.clear();
      str2.clear();
      node.visit(visitor, TreeTraversal::InOrder);
      node.visit(extended_visitor, TreeTraversal::InOrder);
      REQUIRE(str1 == "ACB");
      REQUIRE(str1 == str2);

      str1.clear();
      str2.clear();
      node.visit_leaf(visitor);
      node.visit_leaf(extended_visitor);
      REQUIRE(str1 == "AB");
      REQUIRE(str1 == str2);

      str1.clear();
      str2.clear();
      node.visit_internal(visitor);
      node.visit_internal(extended_visitor);
      REQUIRE(str1 == "C");
      REQUIRE(str1 == str2);

      auto not_vowel = [&str1](node_t const& n) -> bool {
        bool yn = std::string{"AEIOU"}.find(*n) == std::string::npos;
        if (yn) str1 += *n;
        return yn;
      };
      auto extended_not_vowel = [&str2](node_t const& n,
                                        TreeTraversal) -> bool {
        bool yn = std::string{"AEIOU"}.find(*n) == std::string::npos;
        if (yn) str2 += *n;
        return yn;
      };

      auto alphabet = node_t{"W", node_t{"E", "A", "B"}, node_t{"Z", "X", "O"}};

      str1.clear();
      str2.clear();
      alphabet.visit(not_vowel);
      alphabet.visit(extended_not_vowel);
      REQUIRE(str1 == "WZX");
      REQUIRE(str1 == str2);
    }
    SECTION("multi order") {
      auto extended_visitor = [&str2](node_t const& n, TreeTraversal context) {
        if (context == TreeTraversal::PreOrder) {
          str2 += "Pre";
        } else if (context == TreeTraversal::PostOrder) {
          str2 += "Post";
        } else if (context == TreeTraversal::InOrder) {
          str2 += "In";
        } else if (context == TreeTraversal::Any) {
          str2 += "Any";
        } else {
          str2 += "Invalid";
        }

        str2 += *n;
      };

      node.visit(extended_visitor, TreeTraversal::PreAndPostOrder);
      REQUIRE(str2 == "PreCAnyAAnyBPostC");
      str2.clear();

      node.visit(extended_visitor, TreeTraversal::PreAndInOrder);
      REQUIRE(str2 == "PreCAnyAInCAnyB");
      str2.clear();

      node.visit(extended_visitor, TreeTraversal::PostAndInOrder);
      REQUIRE(str2 == "AnyAInCAnyBPostC");
      str2.clear();

      node.visit(extended_visitor, TreeTraversal::Any);
      REQUIRE(str2 == "PreCAnyAInCAnyBPostC");
      str2.clear();
    }
  }
}
