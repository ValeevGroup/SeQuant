//
// created by Bimal Gaudel on 01/19/2020
//

#include "catch.hpp"

#include <iostream>
#include <memory>
#include <SeQuant/domain/evaluate/path_tree.hpp>

TEST_CASE("PathTree is computed", "[path_tree]") {
  using sequant::evaluate::PathTree;

  REQUIRE_NOTHROW(PathTree());

  auto ptree0 = std::make_shared<PathTree>(0);
  REQUIRE(ptree0->get_label() == 0);
  REQUIRE(ptree0->get_children().size() == 0);

  REQUIRE(ptree0->is_leaf() == true);

  auto ptree1 = std::make_shared<PathTree>(1);
  auto ptree2 = std::make_shared<PathTree>(2);
  auto ptree3 = std::make_shared<PathTree>(3);

  ptree0->add_child(ptree1);
  ptree0->add_child(ptree2);
  ptree0->add_child(ptree3);

  REQUIRE(ptree0->is_leaf() == false);
  REQUIRE(ptree1->is_leaf() == true);
  REQUIRE(ptree2->is_leaf() == true);
  REQUIRE(ptree3->is_leaf() == true);

  REQUIRE(ptree0->get_children().size() == 3);

  for (auto i = 0; i < 3; ++i) {
    REQUIRE(ptree0->get_children()[i]->get_label() == i + 1);
  }
}
