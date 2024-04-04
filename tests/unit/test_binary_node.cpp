#include "catch.hpp"

#include <SeQuant/core/binary_node.hpp>

#include <array>
#include <string>
#include <cstddef>
#include <utility>

#include <range/v3/all.hpp>

TEST_CASE("TEST BINARY_NODE", "[FullBinaryNode]") {
  using ranges::views::iota;
  using ranges::views::take;
  using sequant::FullBinaryNode;

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
    auto const node = node_t{"C", "A", "B"};
    REQUIRE(*node == "C");
    REQUIRE(*node.left() == "A");
    REQUIRE(*node.right() == "B");

    std::string str;
    auto visitor = [&str](node_t const& n) { str += *n; };

    node.visit(visitor);
    REQUIRE(str == "ABC");

    str.clear();
    node.visit(visitor, sequant::PreOrder{});
    REQUIRE(str == "CAB");

    str.clear();
    node.visit(visitor, sequant::InOrder{});
    REQUIRE(str == "ACB");

    str.clear();
    node.visit_leaf(visitor);
    REQUIRE(str == "AB");

    str.clear();
    node.visit_internal(visitor);
    REQUIRE(str == "C");

    auto not_vowel = [&str](node_t const& n) -> bool {
      bool yn = std::string{"AEIOU"}.find(*n) == std::string::npos;
      if (yn) str += *n;
      return yn;
    };

    auto alphabet = node_t{"W", node_t{"E", "A", "B"}, node_t{"Z", "X", "O"}};

    str.clear();
    alphabet.visit(not_vowel);
    REQUIRE(str == "WZX");
  }
}
