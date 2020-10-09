#include <SeQuant/domain/factorize/eval_sequence.hpp>

#include <set>
#include <sstream>

#include "catch.hpp"

TEST_CASE("TEST_EVAL_SEQUENCE", "[eval_sequence]") {
  using namespace sequant::factorize;
  auto init_rt_vec = [](size_t n) {
    std::vector<rooted_tree> vec;
    vec.reserve(n);
    for (size_t ii = 0; ii < n; ++ii) vec.emplace_back(rooted_tree{ii});
    return vec;
  };

  SECTION("rooted_tree") {
    REQUIRE_NOTHROW(rooted_tree{0});
    REQUIRE(rooted_tree{0}.children.empty());
    REQUIRE(rooted_tree{0}.label == 0);

    auto t0 = rooted_tree{0};
    REQUIRE_NOTHROW(t0.children.push_back(rooted_tree{0}));
  }

  SECTION("enumerate_eval_sequence") {
    //
    // Return a lambda expression that can be passed to the
    // enumerate_eval_sequence function.
    // The lambda will insert the string representations of
    // eval sequence trees generated by the enumerate_eval_sequence function
    // to to_container
    //
    auto reap_seqs_to = [](std::set<std::wstring>& to_container) {
      return [&to_container](const rooted_tree& t) {
        std::wostringstream oss;
        oss << t;
        oss.flush();
        auto result = to_container.insert(oss.str());
        // assert the same sequence was not inserted previously
        assert(result.second);
      };
    };

    std::set<std::wstring> result1;

    enumerate_eval_sequence(init_rt_vec(3), reap_seqs_to(result1));
    REQUIRE(result1 == std::set<std::wstring>{
                           L"(0 1 2)",  //
                           L"(0 2 1)",  //
                           L"(1 2 0)"   //
                       });

    std::set<std::wstring> result2;
    enumerate_eval_sequence(init_rt_vec(4), reap_seqs_to(result2));
    REQUIRE(result2 == std::set<std::wstring>{
                           L"(0 1 2 3)",    //
                           L"(0 1 3 2)",    //
                           L"(0 2 1 3)",    //
                           L"(0 2 3 1)",    //
                           L"(0 3 1 2)",    //
                           L"(0 3 2 1)",    //
                           L"(1 2 0 3)",    //
                           L"(1 2 3 0)",    //
                           L"(0 3 (1 2))",  //
                           L"(1 3 0 2)",    //
                           L"(1 3 2 0)",    //
                           L"(0 2 (1 3))",  //
                           L"(2 3 0 1)",    //
                           L"(2 3 1 0)",    //
                           L"(0 1 (2 3))",  //
                       });

    //
    // In mathematics, the double factorial or semifactorial of a number n,
    // denoted by n!!, is the product of all the integers from 1 up to n that
    // have the same parity (odd or even) as n. -- wikipedia.org
    //
    auto double_factorial = [](size_t x) {
      if (x == 0) return 1;
      long long n = x;
      auto result = 1;
      if (n % 2 == 0) {
        // even
        for (; n >= 2; n -= 2) result *= n;
      } else {
        // odd
        for (; n >= 1; n -= 2) result *= n;
      }
      return result;
    };

    auto count_num_eval_seqs = [](const std::vector<rooted_tree>& init) {
      //
      // counts the number of trees enumerated by the enumerate_eval_sequence
      // function for a given initial vector of rooted trees init
      //
      size_t count = 0;
      enumerate_eval_sequence(init, [&count](const rooted_tree&) { ++count; });
      return count;
    };

    // the number of ways to evaluate a tensor product depends on the number
    // of factors it has -- obviously
    // the growth follows the odd factorial sequence
    // ie.
    // # factors in a product | # ways to evaluate product
    // -----------------------|---------------------------
    //                     2  | 1     = 1!!
    //                     3  | 3     = 3!!
    //                     4  | 15    = 5!!
    //                     5  | 105   = 7!!
    //                     6  | 945   = 9!!
    //
    // and so on.

    auto iter = 0;
    size_t numFacs = 2, oddFacInp = 1;
    // increasing the iter boundary exponentially increases comp time
    while (iter < 7) {
      REQUIRE(count_num_eval_seqs(init_rt_vec(numFacs)) ==
              double_factorial(oddFacInp));
      ++numFacs;
      oddFacInp += 2;
      ++iter;
    }

  }  // SECTION
}