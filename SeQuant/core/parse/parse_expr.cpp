#include "SeQuant/core/parse_expr.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse/regex_sequant.hpp>
#include <SeQuant/core/parse/rpn.hpp>
#include <SeQuant/core/parse/token_sequant.hpp>
#include <boost/regex.hpp>

namespace sequant {

namespace {

ExprPtr times(ExprPtr const& left, ExprPtr const& right) {
  using Flatten = Product::Flatten;

  if (left.is<Product>() && right.is<Product>()) {
    auto f = left.as<Product>().scalar() * right.as<Product>().scalar();
    auto result =
        Product{f, ranges::begin(*left), ranges::end(*left), Flatten::No};
    auto new_right =
        Product{1, ranges::begin(*right), ranges::end(*right), Flatten::No};
    if (ranges::distance(new_right.factors()) == 1)
      result.append(1, ranges::front(new_right.factors()), Flatten::No);
    else
      result.append(1, ex<Product>(std::move(new_right)), Flatten::No);
    return ex<Product>(std::move(result));
  }

  if (left.is<Constant>() && right.is<Product>()) {
    auto result = Product{};
    result.scale(left.as<Constant>().value());
    result.scale(right.as<Product>().scalar());
    auto new_right =
        ex<Product>(1, ranges::begin(*right), ranges::end(*right), Flatten::No);
    result.append(1, new_right, Flatten::No);
    return ex<Product>(std::move(result));
  }

  if (left.is<Product>()) {
    auto new_left = Product(left->as<Product>().scalar(), ranges::begin(*left),
                            ranges::end(*left), Flatten::No);
    new_left.append(1, right, Flatten::No);
    return ex<Product>(std::move(new_left));
  }

  if (right.is<Product>()) {
    auto new_right =
        Product(right->as<Product>().scalar(), ranges::begin(*right),
                ranges::end(*right), Flatten::No);
    new_right.prepend(1, left, Flatten::No);
    return ex<Product>(std::move(new_right));
  }

  if (left->is_atom() && right->is_atom()) return left * right;

  return ex<Product>(1, ExprPtrList{left, right}, Flatten::No);
}

}  // namespace

namespace parse {
auto throw_invalid_expr = [](std::string const& arg = "Invalid expression") {
  throw std::runtime_error(arg);
};

sequant::rational to_fraction(boost::wsmatch const& match) {
  auto const num = match[1].str();
  auto const denom = match[2].str();
  return {num.empty() ? 1 : std::stol(num),
          denom.empty() ? 1 : std::stol(denom)};
}

std::wstring pure_index_string(boost::wssub_match const& mo) {
  static const boost::wregex rgx{regex_patterns::pure_index_capture()};
  boost::wsmatch match;
  boost::regex_match(mo.begin(), mo.end(), match, rgx);
  return match.str(1) + L"_" + match.str(2);
}

std::vector<Index> parse_pure_indices(boost::wssub_match const& mo) {
  auto result = std::vector<Index>{};
  auto rgx = boost::wregex{regex_patterns::pure_index()};
  auto end = boost::wsregex_iterator{};
  for (auto iter = boost::wsregex_iterator{mo.begin(), mo.end(), rgx};
       iter != end; ++iter) {
    result.emplace_back(pure_index_string((*iter)[0]));
  }
  return result;
}

container::svector<Index> parse_indices(boost::wssub_match const& mo) {
  auto result = container::svector<Index>{};
  auto rgx = boost::wregex{regex_patterns::index_capture()};
  auto end = boost::wsregex_iterator{};
  for (auto iter = boost::wsregex_iterator{mo.begin(), mo.end(), rgx};
       iter != end; ++iter) {
    if ((*iter)[2].matched) {
      auto proto_indices = parse_pure_indices((*iter)[2]);
      auto ispace = proto_indices.begin()->space();

      // all the proto-indices of an index belong to the same IndexSpace
      assert(ranges::all_of(proto_indices, [ispace](auto const& idx) {
        return idx.space() == ispace;
      }));

      result.emplace_back(pure_index_string((*iter)[1]), ispace,
                          std::move(proto_indices));
    } else {
      result.emplace_back(pure_index_string((*iter)[1]));
    }
  }
  return result;
}

std::unique_ptr<Token> to_operand_tensor(boost::wsmatch const& match,
                                         Symmetry s) {
  if (match[4].matched) {
    auto sm_char = match[4].str();
    if (sm_char == L"S")
      s = Symmetry::symm;
    else if (sm_char == L"A")
      s = Symmetry::antisymm;
    else if (sm_char == L"N")
      s = Symmetry::nonsymm;
  }
  return token<OperandTensor>(match[1].str(), parse_indices(match[2]),
                              parse_indices(match[3]), s);
}

std::unique_ptr<Token> to_operand_variable(boost::wsmatch const& match) {
  return token<OperandVariable>(match[1].str(), match[2].matched);
}

}  // namespace parse

namespace deparse {
std::wstring deparse_pure_index(Index const& idx) {
  return idx.label() |
         ranges::views::filter([](wchar_t c) { return c != L'_'; }) |
         ranges::to<std::wstring>;
}

std::wstring deparse_index(Index const& idx) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;

  std::wstring pure = deparse_pure_index(idx);
  if (idx.has_proto_indices()) {
    auto proto = idx.proto_indices() | transform(deparse_pure_index) |
                 intersperse(L",") | join | ranges::to<std::wstring>;
    return pure + L"<" + proto + L">";
  } else
    return pure;
}

template <typename IterableOfIndices>
std::wstring deparse_indices(IterableOfIndices const& indices) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;

  return indices | transform(deparse_index) | intersperse(L",") | join |
         ranges::to<std::wstring>;
}

std::wstring deparse_expr(Tensor const& tnsr, bool annot_sym) {
  std::wstring result = tnsr.label().data();
  result += L"{";
  result += deparse_indices(tnsr.bra());
  result += L";";
  result += deparse_indices(tnsr.ket());
  result += L"}";
  auto tsym = tnsr.symmetry();
  result += !annot_sym                   ? L""
            : tsym == Symmetry::antisymm ? L":A"
            : tsym == Symmetry::symm     ? L":S"
                                         : L":N";
  return result;
}

std::wstring deparse_expr(Constant const& c) {
  auto val = c.value<Complex<sequant::rational>>();
  return val.im == 0
             ? to_wstring(val.re)
             : L"(" + to_wstring(val.re) + L"," + to_wstring(val.im) + L")";
}

std::wstring deparse_expr(Variable const& v) {
  return v.label().data() +
         (v.conjugated() ? std::wstring{L"^*"} : std::wstring{});
}

std::wstring deparse_expr(Product const& prod, bool annot_sym) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;

  std::wstring result;
  result += prod.scalar() == 1 ? L""
            : prod.scalar() == -1
                ? L"-"
                : deparse_expr(Constant{prod.scalar()}) + L" * ";

  result += prod.factors() | transform([annot_sym](ExprPtr const& f) {
              return f.is<Sum>() || f.is<Product>()
                         ? L"(" + deparse_expr(f, annot_sym) + L")"
                         : deparse_expr(f, annot_sym);
            }) |
            intersperse(L" * ") | join | ranges::to<std::wstring>;

  return result;
}

std::wstring deparse_expr(Sum const& sum, bool annot_sym) {
  auto str = deparse_expr(sum.summand(0), annot_sym);
  for (auto&& xpr : ranges::views::tail(sum.summands())) {
    auto to_app = deparse_expr(xpr, annot_sym);
    if (to_app.front() == L'-') {
      str += L" - " +
             std::wstring{std::next(std::begin(to_app)), std::end(to_app)};
    } else {
      str += L" + " + to_app;
    }
  }
  return str;
}

}  // namespace deparse

ExprPtr parse_expr(std::wstring_view raw_expr, Symmetry symmetry) {
  using namespace parse;
  using pattern = regex_patterns;

  auto prev_locale = std::locale{};
  std::locale::global(std::locale::classic());

  static auto const tensor_exp =
      boost::wregex{pattern::tensor_expanded().data()};
  static auto const tensor_terse =
      boost::wregex{pattern::tensor_terse().data()};
  /// matches sequant::Variable
  static auto const variable = boost::wregex{pattern::sequant_variable()};
  static auto const fraction = boost::wregex{pattern::abs_real_frac().data()};
  static auto const times = boost::wregex{L"\\*"};
  static auto const plus = boost::wregex{L"\\+"};
  static auto const minus = boost::wregex{L"\\-"};
  static auto const paren_left = boost::wregex{L"\\("};
  static auto const paren_right = boost::wregex{L"\\)"};
  static auto const space = boost::wregex{L"\\s+?"};

  std::wstring const expr = raw_expr.data();

  boost::wsmatch match;
  auto iter = expr.begin();
  auto validate = [&iter, &match, &expr](boost::wregex const& rgx) -> bool {
    auto yes = boost::regex_search(iter, expr.end(), match, rgx);

    if (yes && iter == match.begin()->first) {
      iter = match.begin()->second;
      return true;
    }
    return false;
  };

  auto rpn = parse::ReversePolishNotation{};

  std::unique_ptr<Token> last_token = token<TokenDummy>();

  auto add_times_if_needed = [&rpn, &last_token]() {
    if (last_token->is<RightParenthesis>() || last_token->is<Operand>())
      rpn.add_operator(token<OperatorTimes>());
  };

  while (iter != expr.end()) {
    if (validate(space)) {
      // do nothing
      // iterator is now at the non-space character
    } else if (validate(fraction)) {
      add_times_if_needed();
      rpn.add_operand(token<OperandConstant>(to_fraction(match)));
      last_token = token<OperandConstant>(0);  // dummy
    } else if (validate(times)) {
      if (last_token->is<Operator>()) throw_invalid_expr();
      rpn.add_operator(token<OperatorTimes>());
      last_token = token<OperatorTimes>();
    } else if (validate(plus)) {
      if (last_token->is<parse::Operator>()) throw_invalid_expr();
      if (last_token->is<TokenDummy>() || last_token->is<LeftParenthesis>()) {
        rpn.add_operator(token<OperatorPlusUnary>());
        last_token = token<OperatorPlusUnary>();
      } else {
        rpn.add_operator(token<OperatorPlus>());
        last_token = token<OperatorPlus>();
      }
    } else if (validate(minus)) {
      if (last_token->is<Operator>()) throw_invalid_expr();
      if (last_token->is<TokenDummy>() || last_token->is<LeftParenthesis>()) {
        rpn.add_operator(token<OperatorMinusUnary>());
        last_token = token<OperatorMinusUnary>();
      } else {
        rpn.add_operator(token<OperatorMinus>());
        last_token = token<OperatorMinus>();
      }
    } else if (validate(paren_left)) {
      add_times_if_needed();
      rpn.add_left_parenthesis();
      last_token = token<LeftParenthesis>();
    } else if (validate(paren_right)) {
      if (last_token->is<Operator>()) throw_invalid_expr();
      if (!rpn.add_right_parenthesis()) throw_invalid_expr();
      last_token = token<RightParenthesis>();
    } else if (validate(tensor_terse) || validate(tensor_exp)) {
      add_times_if_needed();

      rpn.add_operand(to_operand_tensor(match, symmetry));
      last_token =
          token<OperandTensor>(L"Dummy", IndexList{L"i_1"}, IndexList{L"a_1"});
    } else if (validate(variable)) {
      add_times_if_needed();
      rpn.add_operand(to_operand_variable(match));
      last_token = token<OperandVariable>(L"0");  // dummy
    } else {
      throw_invalid_expr();
    }
  }

  // RPN stack built. done parsing
  std::locale::global(prev_locale);

  if (!rpn.close()) throw_invalid_expr();

  if (rpn.tokens().empty()) return nullptr;

  auto result = container::vector<ExprPtr>{};
  for (auto const& t : rpn.tokens()) {
    if (t->is<OperandTensor>())
      result.push_back(t->as<OperandTensor>().clone());
    else if (t->is<OperandConstant>())
      result.push_back(t->as<OperandConstant>().clone());
    else if (t->is<OperandVariable>())
      result.push_back(t->as<OperandVariable>().clone());
    else if (t->is<OperatorPlusUnary>()) {
      if (result.empty()) throw_invalid_expr();
      continue;
    } else if (t->is<OperatorMinusUnary>()) {
      if (result.empty()) throw_invalid_expr();
      auto liter = result.end() - 1;
      *liter = ex<Constant>(-1) * (*liter);
    } else {
      if (result.empty()) throw_invalid_expr();

      auto rhs_operand = result[result.size() - 1];
      result.pop_back();
      auto lhs_operand = result[result.size() - 1];
      result.pop_back();

      if (t->is<OperatorPlus>())
        result.push_back(lhs_operand + rhs_operand);
      else if (t->is<OperatorMinus>())
        result.push_back(lhs_operand - rhs_operand);
      else if (t->is<OperatorTimes>()) {
        result.push_back(sequant::times(lhs_operand, rhs_operand));
      } else
        assert(false && "Unknown token type");
    }
  }

  if (result.size() != 1) throw_invalid_expr();
  return result[0];
}

std::wstring deparse_expr(ExprPtr expr, bool annot_sym) {
  using namespace deparse;
  if (expr->is<Tensor>())
    return deparse_expr(expr->as<Tensor>(), annot_sym);
  else if (expr->is<Sum>())
    return deparse_expr(expr->as<Sum>(), annot_sym);
  else if (expr->is<Product>())
    return deparse_expr(expr->as<Product>(), annot_sym);
  else if (expr->is<Constant>())
    return deparse_expr(expr->as<Constant>());
  else if (expr->is<Variable>())
    return deparse_expr(expr->as<Variable>());
  else
    throw std::runtime_error("Unsupported expr type for deparse!");
}

}  // namespace sequant
