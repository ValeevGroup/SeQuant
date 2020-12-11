#include "catch.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>

#include <iostream>
#include <range/v3/view.hpp>
#include <regex>

auto validate_regex = [](std::wstring_view target, std::wregex rgx) -> bool {
  return std::regex_match(target.data(), rgx);
};

TEST_CASE("TEST_EXPR_PARSE", "[expr_parse]") {
  using namespace sequant::utils::detail;
  using sequant::container::svector;

  SECTION("index") {
    const auto& index_regex =
        std::wregex(expr_rgx_pat.at("indices"), std::regex_constants::nosubs);

    for (std::wstring_view idx : {L"(i_1, i_2, a_1, a_2)",        //
                                  L"(\ni_1,\n i_2,\n a_1, a_2)",  //
                                  L"(\ta_100\t,\n i_198\n)",      //
                                  L"(i_1)",                       //
                                  L"(a_1)",                       //
                                  L"(i1, i2, a1, a2)"})
      REQUIRE(validate_regex(prune_space(idx), index_regex));
  }

  SECTION("braket") {
    const auto& bra_regex = std::wregex(expr_rgx_pat.at("bra"));
    const auto& ket_regex = std::wregex(expr_rgx_pat.at("ket"));

    for (const auto& x : {L"_(i1, i2)", L"_ (i1)"})
      REQUIRE(validate_regex(prune_space(x), bra_regex));

    for (const auto& x : {L"^(i1, i2)", L"^ (i1)"})
      REQUIRE(validate_regex(prune_space(x), ket_regex));
  }

  SECTION("tensor") {
    using index_cont = sequant::container::svector<std::wstring>;
    const auto& tensor_regex = std::wregex(expr_rgx_pat.at("tensor"));

    for (std::wstring_view tnsr : {
             L"t _(i1, i2)^(a1, a2)",   //
             L"g _(i1, i2)^(a1, a2)",   //
             L"I10 ^(a1, a2)_(i1, i2)"  //
         })
      REQUIRE(validate_regex(prune_space(tnsr), tensor_regex));
  }

  SECTION("decimal") {
    const auto& decimal_regex = std::wregex(expr_rgx_pat.at("decimal"));
    for (std::wstring_view num : {L"11.23",  //
                                  L"0.1",    //
                                  L".1",     //
                                  L"1.",     //
                                  L"1.100",  //
                                  L"1100",   //
                                  L"-0.1",   //
                                  L"+0.1"})
      REQUIRE(validate_regex(prune_space(num), decimal_regex));

    REQUIRE(11.23 == Approx(to_decimal(L"11.23")));
    REQUIRE(0.1 == Approx(to_decimal(L"0.1")));
    REQUIRE(.1 == Approx(to_decimal(L".1")));
    REQUIRE(1. == Approx(to_decimal(L"1.")));
    REQUIRE(1.100 == Approx(to_decimal(L"1.100")));
    REQUIRE(1100 == Approx(to_decimal(L"1100")));
    REQUIRE(-.1 == Approx(to_decimal(L"-.1")));
    REQUIRE(+.1 == Approx(to_decimal(L"+.1")));
  }

  SECTION("fraction") {
    const auto& fraction_regex = std::wregex(expr_rgx_pat.at("fraction"));

    for (std::wstring_view frac : {
             L"1/2",      //
             L"1.1/2.2",  //
             L"1",        //
             L"1.",       //
             L".1/.20"    //
         })
      REQUIRE(validate_regex(frac, fraction_regex));

    const auto& num_val_map =
        std::initializer_list<std::pair<std::wstring, double>>{
            {L"1/2", 1.0 / 2.0},            //
            {L"1.1/2.2", 1.1 / 2.2},        //
            {L".1/.20", 0.1 / 0.2},         //
            {L"-.1/-.20", -0.1 / -0.20},    //
            {L"-1.1/-1.20", -1.1 / -1.20},  //
            {L"1.", 1},                     //
            {L"1", 1}                       //
        };

    for (const auto& [str, num] : num_val_map) {
      std::wsmatch match_obj;
      std::regex_match(str, match_obj, fraction_regex);
      REQUIRE(as_fraction(match_obj[1].str(), match_obj[2].str()) ==
              Approx(num));
    }
  }

  SECTION("product") {
    // HELPER LAMBDA
    auto all_tokens = [](std::wstring_view target, const std::wregex& rgx) {
      return target | ranges::views::tokenize(rgx) |
             ranges::to<std::vector<std::wstring>>;
    };
    // HELPER LAMBDA END

    std::wregex product_regex{expr_rgx_pat.at("product_term")};
    const std::wstring p1 = L"t_(i1)^(a1) * g^(a2)_(i2)";
    const std::wstring p2 =
        L"-1/2*t_(i1, i2)^(a1, a2) * t_(i1)^(a1) * g^(a2)_(i2)";
    const std::wstring p3 = L"- t_(i1)^(a1) * g^(a2)_(i2)";
    const std::wstring p4 = L"- t_(i1)^(a1)";

    for (const auto& x : {p1, p2, p3, p4}) {
      auto target = prune_space(x);
      REQUIRE(validate_regex(target, product_regex));
    }

    const std::wstring& p5 =
        prune_space(L"t_(i1)^(a1) * f^(a2)_(i2) + g_(i1, i2)^(a1, a2)");
    const std::wstring& p6 = prune_space(
        L"t_(i1)^(a1) * f^(a2)_(i2) "
        L"  + 1.4/2.3 * g_(i1, i3)^(a1, a3) "
        L"* t_(a2, a3)^(i2, i3)");
    for (const auto& [str, vec] : std::initializer_list<
             std::pair<std::wstring, std::vector<std::wstring>>>{
             {p5, {prune_space(L"t_(i1)^(a1) * f^(a2)_(i2)")}},
             {p6,
              {prune_space(L"t_(i1)^(a1) * f^(a2)_(i2) "),
               prune_space(L"  + 1.4/2.3 * g_(i1, i3)^(a1, a3) * t_(a2, "
                           L"a3)^(i2, i3)")}}}) {
      REQUIRE(all_tokens(str, product_regex) == vec);
    }
  }

  SECTION("term") {
    const auto& term_regex = std::wregex(expr_rgx_pat.at("term"));

    auto single_terms = {
        L"t_(i1, i2)^(a1, a2)",
        L"g^(i1, i2)_(a1, a2)",
        L"g123^(i1, i2)_(a1, a2)",
        L"t_(i1)^(a1) * g^(a2)_(i2)",
        L"-1/2*t_(i1, i2)^(a1, a2) * t_(i1)^(a1) * g^(a2)_(i2)",
        L"- t_(i1)^(a1) * g^(a2)_(i2)",
        L"- t_(i1)^(a1)"};

    for (const auto& t : single_terms) {
      auto target = prune_space(t);
      REQUIRE(validate_regex(target, term_regex));
    }

    auto terms = {L"t_(i1)^(a1) + g_(i1)^(a1)",
                  L"t_(i1)^(a1) * f_(a2)^(i2) + 1/2*g_(i1)^(a1)",
                  L"t_(i1)^(a1) + g_(i1)^(a1) + 0.5/-0.4 * f_(a2)^(i2)",
                  L"t_(i1)^(a1) + g_(i1)^(a1) * f_(a2)^(i2)"};

    std::vector<std::vector<std::wstring>> terms_extract{
        {L"t_(i1)^(a1)", L" + g_(i1)^(a1)"},
        {L"t_(i1)^(a1) * f_(a2)^(i2) ", L"+ 1/2*g_(i1)^(a1)"},
        {L"t_(i1)^(a1) ", L"+ g_(i1)^(a1) ", L"+ 0.5/-0.4 * f_(a2)^(i2)"},
        {L"t_(i1)^(a1) ", L"+ g_(i1)^(a1) * f_(a2)^(i2)"}};
    for (const auto& [t, x] : ranges::views::zip(terms, terms_extract)) {
      auto target = prune_space(t);

      auto to_match = x | ranges::views::transform(prune_space) |
                      ranges::to<std::vector<std::wstring>>;

      auto xtrxt = target | ranges::views::tokenize(term_regex) |
                   ranges::to<std::vector<std::wstring>>;
      REQUIRE(to_match == xtrxt);
    }
  }

  SECTION("sum") {
    auto sum_regex = std::wregex(expr_rgx_pat.at("sum"));
    auto terms = {L"t_(i1, i2)^(a1, a2)",
                  L"g^(i1, i2)_(a1, a2)",
                  L"g123^(i1, i2)_(a1, a2)",
                  L"t_(i1)^(a1) * g^(a2)_(i2)",
                  L"-1/2*t_(i1, i2)^(a1, a2) * t_(i1)^(a1) * g^(a2)_(i2)",
                  L"- t_(i1)^(a1) * g^(a2)_(i2)",
                  L"- t_(i1)^(a1)"};

    for (const auto& t : terms) {
      auto target = prune_space(t);
      REQUIRE_FALSE(validate_regex(target, sum_regex));
    }

    terms = {L"t_(i1)^(a1) + g_(i1)^(a1)",
             L"t_(i1)^(a1) * f_(a2)^(i2) + 1/2*g_(i1)^(a1)",
             L"t_(i1)^(a1) + g_(i1)^(a1) + 0.5/-0.4 * f_(a2)^(i2)",
             L"t_(i1)^(a1) + g_(i1)^(a1) * f_(a2)^(i2)"};

    for (const auto& t : terms) {
      auto target = prune_space(t);
      REQUIRE(validate_regex(target, sum_regex));
    }
  }
}

TEST_CASE("TEST_MAKE_EXPR_BY_PARSE", "[expr_parse]") {
  using namespace std::string_literals;
  using namespace sequant;
  using utils::parse_expr;

  using index_list = container::svector<Index>;

  auto A = ex<Tensor>(Tensor{L"A", index_list{L"i_1", L"i_2"},
                             index_list{L"a_1", L"a_2"}, Symmetry::antisymm});

  auto g1 = ex<Tensor>(Tensor{L"g", index_list{L"i_3", L"i_4"},
                              index_list{L"a_3", L"a_4"}, Symmetry::antisymm});

  auto t1 = ex<Tensor>(
      Tensor{L"t", {L"a_3", L"a_4"}, {L"i_1", L"i_2"}, Symmetry::antisymm});

  auto t2 = ex<Tensor>(
      Tensor{L"t", {L"a_1", L"a_2"}, {L"i_3", L"i_4"}, Symmetry::antisymm});

  auto t3 = ex<Tensor>(Tensor{L"t", {L"a_1"}, {L"i_3"}, Symmetry::antisymm});

  auto t4 = ex<Tensor>(
      Tensor{L"t", {L"a_2", L"a_3"}, {L"i_1", L"i_2"}, Symmetry::antisymm});

  auto f1 = ex<Tensor>(Tensor{L"f", {L"i_3"}, {L"a_3"}, Symmetry::antisymm});

  auto prod1 = ex<Product>(Product{1.0 / 16, {A, g1, t1, t2}});

  auto prod2 = ex<Product>(Product{1.0 / 2, {A, f1, t3, t4}});

  auto sum1 = ex<Sum>(Sum{prod1, prod2});
  auto sum2 =
      ex<Sum>(Sum{prod1, ex<Tensor>(Tensor{L"g", index_list{L"a_1, a_2"},
                                           index_list{L"i_1", L"i_2"},
                                           Symmetry::antisymm})});
  const auto& str_g1_1 =
      L"g"
      "_(i3, i4)"
      "^(a_3, a_4)";

  const auto& str_g1_2 =
      L"g"
      "^(a_3, a_4)"
      "_(i3, i4)";

  const auto& str_A = L"A_(i1, i2)^(a1, a2)";
  const auto& str_t1 = L"t^(i1, i2)_(a3, a4)";
  const auto& str_t2 = L"t_(a1, a2)^(i3, i4)";
  const auto& str_t3 = L"t_(a_1)^(i_3)";
  const auto& str_t4 = L"t_(a2, a3)^(i1, i2)";
  const auto& str_f1 = L"f_(i3)^(a3)";

  for (const auto& [s, t] :
       container::map<std::wstring_view, ExprPtr>{{str_g1_1, g1},
                                                  {str_g1_2, g1},
                                                  {str_t1, t1},
                                                  {str_t2, t2},
                                                  {str_t3, t3},
                                                  {str_t4, t4},
                                                  {str_f1, f1}}) {
    auto parsed = parse_expr(s, Symmetry::antisymm);
    REQUIRE(parsed);
    REQUIRE(parsed->is<Tensor>());
    REQUIRE(*parsed == *t);
  }

  const auto& str_prod1 = L"  1/16 * "s + str_A + L" * "s + str_g1_1 + L" * "s +
                          str_t1 + L" * "s + str_t2;

  const auto& str_prod2 = L" 1/2 *"s + str_A + L" * "s + str_f1 + L" * "s +
                          str_t3 + L" * "s + str_t4;

  for (const auto& [s, e] : container::map<std::wstring_view, ExprPtr>{
           {str_prod1, prod1}, {str_prod2, prod2}}) {
    auto parsed = parse_expr(s, Symmetry::antisymm);
    REQUIRE(parsed);
    REQUIRE(parsed->is<Product>());
    // std::wcout << "Orig:   " << e->to_latex() << "\n"
    //            << "Parsed: " << parsed->to_latex() << "\n\n";
    // REQUIRE(*parsed == *e);
  }
  // std::wcout.flush();

  const auto& str_sum1 = str_prod1 + L" + "s + str_prod2;
  const auto& str_sum2 = str_prod1 + L" + g_(a1, a2)^(i1, i2) ";

  for (const auto& [s, e] : container::map<std::wstring_view, ExprPtr>{
           {str_sum1, sum1}, {str_sum2, sum2}}) {
    auto parsed = parse_expr(s, Symmetry::antisymm);
    REQUIRE(parsed);
    REQUIRE(parsed->is<Sum>());
    // REQUIRE(*parsed == *e);
  }

  auto parsed = parse_expr(str_sum2, Symmetry::antisymm);
  REQUIRE(parsed);
  REQUIRE(parsed->at(1)->is<Tensor>());
}
