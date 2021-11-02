#include "SeQuant/core/parse_expr.hpp"

#include <SeQuant/core/parse/regex_sequant.hpp>
#include <SeQuant/core/parse/rpn.hpp>
#include <SeQuant/core/parse/token_sequant.hpp>
#include <boost/regex.hpp>

namespace sequant::parse {
std::wstring prune_space(std::wstring const& raw) {
  return boost::regex_replace(raw, boost::wregex(L"[[:space:]]+"), L"");
}

double to_decimal(std::wstring_view numstr) {
  double num;
  std::wistringstream wiss{numstr.data()};
  wiss >> num;
  return wiss ? num : 0;
}

double to_fraction(boost::wsmatch const& match) {
  auto const num = match[1].str();
  auto const denom = match[2].str();
  return (num.empty() ? 1.0 : to_decimal(num)) /
         (denom.empty() ? 1.0 : to_decimal(denom));
}

container::vector<Index> to_indices(std::wstring_view raw_csv_indices) {
  using ranges::views::transform;
  using ranges::views::split;
    std::wostringstream oss{};
    wchar_t last_char = L' ';
    for (auto x: raw_csv_indices){
      if (std::iswalpha(last_char) && std::iswdigit(x))
        oss << L'_';
      oss << x;
      last_char = x;
    }
    auto with_uscores = oss.str();

    return split(with_uscores,L',')
           | transform([](auto&& idx_rng){
             return Index{ranges::to<std::wstring>(idx_rng)}; })
           | ranges::to<container::vector<Index>>;
}

std::unique_ptr<Token> to_operand_tensor(boost::wsmatch const& match, Symmetry s){
  if (match[4].matched){
    auto sm_char = match[4].str();
    if (sm_char == L"S") s = Symmetry::symm;
    else if (sm_char == L"A") s = Symmetry::antisymm;
    else if (sm_char == L"N") s = Symmetry::nonsymm;
  }
  return token<OperandTensor>(match[1].str(),
                              to_indices(match[2].str()),
                              to_indices(match[3].str()), s);
}
}

namespace sequant {

ExprPtr parse_expr(std::wstring_view raw_expr, Symmetry symmetry){
  using namespace parse;
  using pattern = regex_patterns;

  auto throw_invalid_expr = []() {
    throw std::runtime_error("Invalid expression!");
  };

  static auto const tensor_exp = boost::wregex{pattern::tensor_expanded().data()};
  static auto const tensor_terse = boost::wregex{pattern::tensor_terse().data()};
  static auto const fraction = boost::wregex{pattern::fraction().data()};
  static auto const times = boost::wregex{L"\\*"};
  static auto const plus  = boost::wregex{L"\\+"};
  static auto const minus = boost::wregex{L"\\-"};
  static auto const paren_left = boost::wregex{L"\\("};
  static auto const paren_right = boost::wregex{L"\\)"};

  auto const expr = parse::prune_space(raw_expr.data());

  boost::wsmatch match;
  auto iter = expr.begin();
  auto validate = [&iter, &match, &expr](boost::wregex const& rgx) -> bool {
    boost::regex_search(iter, expr.end(), match, rgx);

    if (match.empty() || iter != match.begin()->first) return false;
    else {
      iter = match.begin()->second;
      return true;
    }
  };

  auto rpn = parse::ReversePolishNotation{};

  std::unique_ptr<Token> last_token = token<TokenDummy>();

  auto add_times_if_needed = [&rpn, &last_token]() {
    if (last_token->is<RightParenthesis>()
        || last_token->is<Operand>())
      rpn.add_operator(token<OperatorTimes>());
  };

  while (iter != expr.end()) {
    if (validate(tensor_terse)||validate(tensor_exp)) {
      add_times_if_needed();

      rpn.add_operand(to_operand_tensor(match, symmetry));
      last_token = token<OperandTensor>(L"Dummy",
                                        IndexList{L"i_1"},
                                        IndexList{L"a_1"});
    } else if (validate(fraction)) {
      add_times_if_needed();
      rpn.add_operand(token<OperandConstant>(to_fraction(match)));
      last_token = token<OperandConstant>(0); // dummy
    } else if (validate(times)) {
      if (last_token->is<Operator>()) throw_invalid_expr();
      rpn.add_operator(token<OperatorTimes>());
      last_token = token<OperatorTimes>();
    } else if (validate(plus)) {
      if (last_token->is<parse::Operator>()) throw_invalid_expr();
      if (last_token->is<TokenDummy>() || last_token->is<LeftParenthesis>()) {
        rpn.add_operator(token<OperatorPlusUnary>());
        last_token = token<OperatorPlusUnary>();
      }
      else {
        rpn.add_operator(token<OperatorPlus>());
        last_token = token<OperatorPlus>();
      }
    } else if (validate(minus)) {
      if (last_token->is<Operator>()) throw_invalid_expr();
      if (last_token->is<TokenDummy>() || last_token->is<LeftParenthesis>()) {
        rpn.add_operator(token<OperatorMinusUnary>());
        last_token = token<OperatorMinusUnary>();
      }
      else {
        rpn.add_operator(token<OperatorMinus>());
        last_token = token<OperatorMinus>();
      }
    } else if (validate(paren_left)) {
      add_times_if_needed();
      rpn.add_left_parenthesis();
      last_token = token<LeftParenthesis>();
    } else if (validate(paren_right)){
      if (last_token->is<Operator>()) throw_invalid_expr();
      if (!rpn.add_right_parenthesis())
        throw_invalid_expr();
      last_token = token<RightParenthesis>();
    } else {
      throw_invalid_expr();
    }
  }

  if (!rpn.close()) throw_invalid_expr();

  if (rpn.tokens().empty())
    return nullptr;

  auto result = container::vector<ExprPtr>{};
  for (auto const& t: rpn.tokens()){
    if (t->is<OperandTensor>())
      result.push_back(t->as<OperandTensor>().clone());
    else if (t->is<OperandConstant>())
      result.push_back(t->as<OperandConstant>().clone());
    else if (t->is<OperatorPlusUnary>()) {
      if (result.empty()) throw_invalid_expr();
      continue;
    }
    else if (t->is<OperatorMinusUnary>()){
      if (result.empty()) throw_invalid_expr();
      auto liter = result.end() - 1;
      *liter = ex<Constant>(-1)*(*liter);
    }
    else {
      if (result.empty()) throw_invalid_expr();

      auto rhs_operand = result[result.size()-1];
      result.pop_back();
      auto lhs_operand = result[result.size()-1];
      result.pop_back();

      if (t->is<OperatorPlus>())
        result.push_back(lhs_operand + rhs_operand);
      else if (t->is<OperatorMinus>())
        result.push_back(lhs_operand - rhs_operand);
      else if (t->is<OperatorTimes>()) {
        auto prod_ptr = lhs_operand->is<Product>() ? lhs_operand:
                                                   ex<Product>(1, ExprPtrList{lhs_operand});
        auto& prod = prod_ptr->as<Product>();
        auto append_prod = [&prod](ExprPtr lrhs){
          auto& p = lrhs->as<Product>();
          prod.scale(p.scalar());
          p.scale(1./p.scalar());
          if (p.size() == 1)
            prod.append(1.,p.factor(0));
          else prod.append(lrhs);
        };

        if (rhs_operand->is<Product>())
          append_prod(rhs_operand);
        else prod.append(1., rhs_operand);

        result.push_back(prod_ptr);
      }
      else assert(false && "Unknown token type");
    }
  }

  if (result.size() != 1) throw_invalid_expr();
  return result[0];
}

ExprPtr parse_expr_asymm(std::wstring_view raw) {
  return parse_expr(raw, Symmetry::antisymm);
}

}  // namespace sequant
