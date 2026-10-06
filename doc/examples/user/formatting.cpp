#include <SeQuant/core/io/format.hpp>

#include <SeQuant/core/expr.hpp>

#include <format>
#include <iostream>
#include <sstream>

int main() {
  // start-snippet-1
  const auto expr = sequant::ex<sequant::Variable>(L"x");
  std::cout << expr << '\n';  // {x}

  const auto latex = std::format("{}", expr);                       // {x}
  const auto short_latex = std::format("{:l}", expr);               // {x}
  const auto named_latex = std::format("{:latex}", expr);           // {x}
  const auto serialized = std::format("{:s}", expr);                // x
  const auto named_serialized = std::format("{:serialize}", expr);  // x
  // end-snippet-1

  // start-snippet-2
  std::wostringstream stream;
  stream << expr;
  const auto wide_latex = std::format(L"{:l}", expr);        // L"{x}"
  const auto wide_serialized = std::format(L"{:s}", expr);   // L"x"
  const auto empty = std::format("{}", sequant::ExprPtr{});  // NULL
  // end-snippet-2

  return latex == "{x}" && short_latex == latex && named_latex == latex &&
                 serialized == "x" && named_serialized == serialized &&
                 stream.str() == L"{x}" && wide_latex == L"{x}" &&
                 wide_serialized == L"x" && empty == "NULL"
             ? 0
             : 1;
}
