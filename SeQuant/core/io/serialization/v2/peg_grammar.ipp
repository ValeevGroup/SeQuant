using namespace std::literals;
constexpr static std::string_view peg_serialization_grammar = R"(

Real <- < [0-9]* '.' [0-9]+ / [0-9]+ >
# the 'i' must follow the number directly and not start a label, so that e.g.
# '2 i' and '2ix' are products with variables i and ix
Imaginary <- < Real 'i' !IdentifierChar >
Complex <- Imaginary / Real
Constant <- Complex

# A bare identifier excludes the ASCII characters that have a meaning in the
# grammar (and whitespace and control characters) and allows everything else;
# serialize.cpp's is_bare_label() must agree with it. Any other label is written
# in backticks, with '\' escaping a '`' or '\' in it.
IdentifierChar <- [^\u0000-&(-/:-@[-^`{-\u00A0]
BareIdentifier <- < IdentifierChar+ >
QuotedIdentifier <- < '`' ( '\\' . / [^`\\] )* '`' >
Identifier <- QuotedIdentifier / BareIdentifier

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

# A slot inside a bra/ket/aux group may hold a blank (null) index, spelled '_'
# or left empty. The empty alternative deliberately matches nothing so that the
# position is preserved in the produced value list; a group that is empty
# altogether has no slots, so a lone blank must be spelled '_'.
BlankIndex <- '_' &[,;\]}]
EmptyIndex <- ''
IndexItem <- BlankIndex / Index / EmptyIndex
Indices <- IndexItem (',' IndexItem)*
IndexGroups <- Indices (';' Indices)*

# --- Symmetry --------------------------------------------------------------

# Longer (prefixed) tokens must precede their single-letter prefixes.
PredefinedSymm <- 'bkS' / 'bkC' / 'bkN' / 'pS' / 'pN' / 'A' / 'S' / 'N'
Symmetry <- < PredefinedSymm >
SymmetrySpec <- ':' Symmetry ('-' Symmetry)*

Tensor <- Identifier '[' IndexGroups ']' SymmetrySpec?

# The label selects the statistics (and the vacuum) of the operator.
NormalOrderedOperator <- Identifier '{' IndexGroups '}'

# --- Power -----------------------------------------------------------------

PowerBase <- '(' Expression ')' / Constant / Variable
Exponent <- < [-+]? [0-9]+ ('/' [0-9]+)? >
# The second form is the one v1 writes for a conjugated power.
Power <- PowerBase '^' '(' Exponent ')' Conjugation? / '(' Power ')' Conjugation

# --- Functions (extension hook) --------------------------------------------

# ~name(...) dispatches to a registered handler, which receives the text
# between the parentheses as written and parses it itself. Only brackets are
# matched here, so that the call ends at the right ')'.
RawText <- ( '(' RawText ')' / '[' RawText ']' / '{' RawText '}' / QuotedIdentifier / ![()[\]{}`] . )*
RawArguments <- < RawText >
Function <- '~' Identifier '(' RawArguments ')'

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
