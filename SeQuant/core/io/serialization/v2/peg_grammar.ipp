using namespace std::literals;
constexpr static std::string_view peg_serialization_grammar = R"(

Integer <- < [0-9]+ >
Float <- < Integer? '.' Integer >
Real <- Float / Integer
Imaginary <- Real 'i'
Complex <- Imaginary / Real
Constant <- Complex

# This is a very liberal rule that simply tries to exclude the most important
# non-ID characters (in ASCII range) and then allows everything else
Identifier <- ( [^\u0000-\u002f\u003a-\u0040\u005b-\u0060\u007b-\u00a0] / '_' )+

Variable <- < Identifier ('^*')? >

Index <- < Identifier >

Indices <- ','* Index (','+ Index)*

IndexGroups <- ';'* Indices (';'+ Indices)* ';'*

PredefinedSymm <- 'A' / 'S' / 'N' / 'bkS' / 'bkC' / 'bkN' / 'pS' / 'pN'

Symmetry <- PredefinedSymm

Tensor <- Identifier ('^*')? '[' IndexGroups ']' ( ':' Symmetry (',' Symmetry)* )?

Statistic <- 'F' / 'B'

NormalOrderedOperator <- Identifier '{' IndexGroups '}' (':' Statistic)?

Function <- '~' Identifier '(' [^)]* ')'


Nullary <- '(' Expression ')' / Function / Constant / Tensor / NormalOrderedOperator / Variable

UnaryOperator <- '+' / '-'

Unary <- UnaryOperator Nullary

BinaryOperator <- '+' / '-' / '*' / '/'

ExprAtom <- Unary / Nullary

InfixExpression <- ExprAtom ( BinaryOperator ExprAtom )* {
  precedence
    L - +
    L * /
}

Expression <- ExprAtom Expression / InfixExpression Expression?

Result <- Tensor / Variable

ResultExpression <- Result '=' Expression

%whitespace  <-  [ \t\r\n]*

)"sv;
