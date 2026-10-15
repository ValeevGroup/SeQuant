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
# U+207A/U+207B (superscript +/-) are excluded so that they are available as
# the Unicode creator/annihilator marks of an action-tagged index.
Identifier <- < ( [^\u0000-/:-@[-`{-\u00A0⁺⁻] / '_' )+ >

Variable <- Identifier Conjugation?

Conjugation <- '^*'

# --- Indices ---------------------------------------------------------------

# An index may carry an explicit, otherwise unnamed, basis instance  i#N
# (N may be negative) and/or a proto-index domain  i<p1,p2> . A proto index may
# itself carry a basis instance  p#N .
BasisInstance <- '#' < '-'? [0-9]+ >
ProtoIndex <- Identifier BasisInstance?
ProtoList <- ProtoIndex (',' ProtoIndex)*
IndexDomain <- '<' ProtoList? '>'
Index <- Identifier BasisInstance? IndexDomain?

# A slot inside a bra/ket/aux group may be left blank (positional empty index).
# The empty alternative deliberately matches nothing so that the position is
# preserved in the produced value list.
EmptyIndex <- ''
IndexItem <- Index / EmptyIndex
Indices <- IndexItem (',' IndexItem)*
IndexGroups <- Indices (';' Indices)*

# --- Symmetry --------------------------------------------------------------

# Longer (prefixed) tokens must precede their single-letter prefixes.
PredefinedSymm <- 'bkS' / 'bkC' / 'bkN' / 'pS' / 'pN' / 'A' / 'S' / 'N'
Symmetry <- < PredefinedSymm >
SymmetrySpec <- ':' Symmetry ('-' Symmetry)*

Tensor <- Identifier '[' IndexGroups ']' SymmetrySpec?

Statistic <- 'F' / 'B'
NormalOrderedOperator <- Identifier '{' IndexGroups '}' (':' Statistic)?

# --- Power -----------------------------------------------------------------

PowerBase <- '(' Expression ')' / Constant / Variable
Power <- PowerBase '^' '(' Real ')' Conjugation?

# --- Functions (extension hook) --------------------------------------------

# ~name(arg, arg, ...) dispatches to a registered handler. An argument is
# either an action-tagged index (creator '+'/'⁺', annihilator '-'/'⁻')
# or an arbitrary sub-expression. The optional ':' Statistic selects Fermi/Bose.
CreatorMark <- '+' / '⁺'
AnnihilatorMark <- '-' / '⁻'
ActionIndex <- Index ( CreatorMark / AnnihilatorMark )
FunctionArg <- ActionIndex / Expression
ArgList <- FunctionArg (',' FunctionArg)*
Function <- '~' Identifier '(' ArgList? ')' (':' Statistic)?

# --- Expression ------------------------------------------------------------

Nullary <- Function / Power / '(' Expression ')' / Tensor / NormalOrderedOperator / Constant / Variable

# A leading unary sign binds tighter than a product. Juxtaposition (implicit
# multiplication) binds like an explicit '*'; because the juxtaposed factor is
# a bare Nullary (no sign), a '+'/'-' can never be swallowed as a factor and so
# always reaches the Sum level below.
Sign <- < [-+] >
Factor <- Sign* Nullary

MulOp <- < [*/] >
ProductTail <- MulOp Factor / Nullary
Product <- Factor ProductTail*

AddOp <- < [-+] >
Sum <- Product ( AddOp Product )*

Expression <- Sum

Result <- Tensor / Variable

ResultExpression <- Result '=' Expression

%whitespace  <-  [ \t\r\n]*

)"sv;
