//
// Created by Robert Adam on 2023-09-22
//

#ifndef SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP
#define SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/v1/ast.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <boost/variant.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/find.hpp>

#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <variant>

namespace sequant::io::serialization::v1::transform {

using DefaultSymmetries =
    std::tuple<Symmetry, std::variant<BraKetSymmetry, Hermiticity>,
               ColumnSymmetry>;

template <typename AST, typename PositionCache, typename Iterator>
std::tuple<std::size_t, std::size_t> get_pos(const AST &ast,
                                             const PositionCache &cache,
                                             const Iterator &begin) {
  const auto range = cache.position_of(ast);

  return {std::distance(begin, range.begin()),
          std::distance(range.begin(), range.end())};
}

template <typename PositionCache, typename Iterator>
Index to_index(const io::serialization::v1::ast::Index &index,
               const PositionCache &position_cache, const Iterator &begin) {
  container::vector<Index> protoIndices;
  protoIndices.reserve(index.protoLabels.size());

  for (const io::serialization::v1::ast::IndexLabel &current :
       index.protoLabels) {
    try {
      std::wstring label = current.label + L"_" + std::to_wstring(current.id);
      IndexSpace space =
          get_default_context().index_space_registry()->retrieve(label);
      protoIndices.push_back(Index(std::move(label), std::move(space)));
    } catch (const IndexSpace::bad_key &) {
      auto [offset, length] = get_pos(current, position_cache, begin);
      throw SerializationError(offset, length,
                               "Unknown index space '" + toUtf8(current.label) +
                                   "' in proto index specification");
    } catch (const Exception &e) {
      auto [offset, length] = get_pos(current, position_cache, begin);
      throw SerializationError(offset, length,
                               "Invalid index '" + toUtf8(current.label) + "_" +
                                   std::to_string(current.id) + ": " +
                                   e.what());
    }
  }

  try {
    IndexSpace space = get_default_context().index_space_registry()->retrieve(
        index.label.label);
    return Index(std::move(space), index.label.id, std::move(protoIndices));
  } catch (const IndexSpace::bad_key &e) {
    auto [offset, length] = get_pos(index.label, position_cache, begin);
    throw SerializationError(offset, length,
                             "Unknown index space '" +
                                 toUtf8(index.label.label) +
                                 "' in index specification");
  } catch (const Exception &e) {
    auto [offset, length] = get_pos(index.label, position_cache, begin);
    throw SerializationError(offset, length,
                             "Invalid index '" + toUtf8(index.label.label) +
                                 "_" + std::to_string(index.label.id) + ": " +
                                 e.what());
  }
}

template <typename PositionCache, typename Iterator>
std::tuple<container::vector<Index>, container::vector<Index>,
           container::vector<Index>>
make_indices(const io::serialization::v1::ast::IndexGroups &groups,
             const PositionCache &position_cache, const Iterator &begin) {
  container::vector<Index> braIndices;
  container::vector<Index> ketIndices;
  container::vector<Index> auxiliaries;
  container::vector<Index> auxIndices;

  static_assert(std::is_same_v<decltype(groups.bra), decltype(groups.ket)>,
                "Types for bra and ket indices must be equal for pointer "
                "aliasing to work");

  const auto *bra = &groups.bra;
  const auto *ket = &groups.ket;
  if (groups.reverse_bra_ket) {
    bra = &groups.ket;
    ket = &groups.bra;
  }

  braIndices.reserve(bra->size());
  ketIndices.reserve(ket->size());

  for (const io::serialization::v1::ast::Index &current : *bra) {
    braIndices.push_back(to_index(current, position_cache, begin));
  }
  for (const io::serialization::v1::ast::Index &current : *ket) {
    ketIndices.push_back(to_index(current, position_cache, begin));
  }
  for (const io::serialization::v1::ast::Index &current : groups.auxiliaries) {
    auxiliaries.push_back(to_index(current, position_cache, begin));
  }

  return {std::move(braIndices), std::move(ketIndices), std::move(auxiliaries)};
}

template <typename Iterator>
Symmetry to_perm_symmetry(char c, std::size_t offset, const Iterator &,
                          Symmetry default_symmetry) {
  if (c == io::serialization::v1::ast::SymmetrySpec::unspecified) {
    return default_symmetry;
  }

  switch (c) {
    case 'A':
    case 'a':
      return Symmetry::Antisymm;
    case 'S':
    case 's':
      return Symmetry::Symm;
    case 'N':
    case 'n':
      return Symmetry::Nonsymm;
  }

  throw SerializationError(
      offset, 1, std::string("Invalid symmetry specifier '") + c + "'");
}

template <typename Iterator>
std::variant<BraKetSymmetry, Hermiticity> to_braket_symmetry(
    char c, std::size_t offset, const Iterator &,
    std::variant<BraKetSymmetry, Hermiticity> default_symmetry) {
  if (c == io::serialization::v1::ast::SymmetrySpec::unspecified) {
    return default_symmetry;
  }

  // The letter is either a pinned exchange symmetry or a hermiticity trait.
  // 'C' / 'S' / 'N' are concrete BraKetSymmetry::{Conjugate, Symm, Nonsymm}
  // values, from which Tensor back-fills the traits that derive them; 'H'
  // (Hermitian) and 'A' (AntiHermitian) name the field-agnostic trait itself
  // and leave the exchange symmetry to be derived from it, the parity and the
  // indices' basis at Tensor construction. The two are kept apart: reading a
  // pin as a trait instead would resolve 'C' against Real-field indices to
  // Symm and lose the AntiHermitian preimage of 'N'. The serializer spells a
  // definite hermiticity with its trait letter (see serialize_symm), so a pin
  // letter reaches here only from hand-written or older input.
  switch (c) {
    case 'C':
    case 'c':
      return BraKetSymmetry::Conjugate;
    case 'S':
    case 's':
      return BraKetSymmetry::Symm;
    case 'N':
    case 'n':
      return BraKetSymmetry::Nonsymm;
    case 'H':
    case 'h':
      return Hermiticity::Hermitian;
    case 'A':
    case 'a':
      return Hermiticity::AntiHermitian;
  }

  throw SerializationError(
      offset, 1,
      std::string("Invalid BraKet symmetry / Hermiticity specifier '") + c +
          "'");
}

template <typename Iterator>
std::optional<ConjugationParity> to_conjugation_parity(char c,
                                                       std::size_t offset,
                                                       const Iterator &) {
  if (c == io::serialization::v1::ast::SymmetrySpec::unspecified) {
    return std::nullopt;
  }

  switch (c) {
    case 'E':
    case 'e':
      return ConjugationParity::Even;
    case 'O':
    case 'o':
      return ConjugationParity::Odd;
    case 'N':
    case 'n':
      return ConjugationParity::None;
  }

  throw SerializationError(
      offset, 1,
      std::string("Invalid conjugation parity specifier '") + c + "'");
}

template <typename Iterator>
ColumnSymmetry to_column_symmetry(char c, std::size_t offset, const Iterator &,
                                  ColumnSymmetry default_symmetry) {
  if (c == io::serialization::v1::ast::SymmetrySpec::unspecified) {
    return default_symmetry;
  }

  switch (c) {
    case 'S':
    case 's':
      return ColumnSymmetry::Symm;
    case 'N':
    case 'n':
      return ColumnSymmetry::Nonsymm;
  }

  throw SerializationError(
      offset, 1,
      std::string("Invalid particle symmetry specifier '") + c + "'");
}

template <typename PositionCache, typename Iterator>
Constant to_constant(const io::serialization::v1::ast::Number &number,
                     const PositionCache &, const Iterator &) {
  const ::sequant::rational magnitude =
      (static_cast<std::int64_t>(number.numerator) == number.numerator &&
       static_cast<std::int64_t>(number.denominator) == number.denominator)
          // Integer fraction
          ? ::sequant::rational(static_cast<std::int64_t>(number.numerator),
                                static_cast<std::int64_t>(number.denominator))
          // Construct from floating point value
          : ::sequant::rational(number.numerator / number.denominator);

  // an imaginary literal is the same magnitude on the imaginary axis
  if (number.imaginary) {
    return Constant(Constant::scalar_type{::sequant::rational{0}, magnitude});
  }
  return Constant(magnitude);
}

template <typename PositionCache, typename Iterator>
std::tuple<Symmetry, std::variant<BraKetSymmetry, Hermiticity>, ColumnSymmetry,
           std::optional<ConjugationParity>>
to_symmetries(
    const boost::optional<io::serialization::v1::ast::SymmetrySpec> &symm_spec,
    const DefaultSymmetries &default_symms, const PositionCache &cache,
    const Iterator &begin) {
  if (!symm_spec.has_value()) {
    return {std::get<0>(default_symms), std::get<1>(default_symms),
            std::get<2>(default_symms), std::nullopt};
  }

  const ast::SymmetrySpec &spec = symm_spec.get();

  auto [offset, length] = get_pos(spec, cache, begin);

  // Note: symmetry specifications are a separator (colon or dash) followed by
  // an uppercase letter each (no whitespace allowed in-between)
  Symmetry perm_symm = to_perm_symmetry(spec.perm_symm, offset + 1, begin,
                                        std::get<0>(default_symms));
  std::variant<BraKetSymmetry, Hermiticity> braket_symm = to_braket_symmetry(
      spec.braket_symm, offset + 3, begin, std::get<1>(default_symms));
  ColumnSymmetry column_symm = to_column_symmetry(
      spec.column_symm, offset + 5, begin, std::get<2>(default_symms));
  // the fourth letter is optional and has no Context-level default: absent
  // means "let Tensor::resolve_symmetries derive it", not "Even"
  std::optional<ConjugationParity> parity =
      to_conjugation_parity(spec.conjugation_parity, offset + 7, begin);

  return {perm_symm, braket_symm, column_symm, parity};
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const io::serialization::v1::ast::Product &product,
                    const PositionCache &position_cache, const Iterator &begin,
                    const DefaultSymmetries &default_symms);
template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const io::serialization::v1::ast::Sum &sum,
                    const PositionCache &position_cache, const Iterator &begin,
                    const DefaultSymmetries &default_symms);

/// the vacuum of a tilde-labelled NormalOperator: NormalOperator::label()
/// is the tilde label for every non-Physical vacuum
inline Vacuum tilde_vacuum(Statistics s) {
  const auto vac = get_default_context(s).vacuum();
  return vac == Vacuum::Physical ? Vacuum::SingleProduct : vac;
}

template <typename PositionCache, typename Iterator>
struct Transformer {
  std::reference_wrapper<const PositionCache> position_cache;
  std::reference_wrapper<const Iterator> begin;
  std::reference_wrapper<const DefaultSymmetries> default_symms;

  ExprPtr operator()(const io::serialization::v1::ast::Product &product) const {
    return ast_to_expr<PositionCache>(product, position_cache.get(),
                                      begin.get(), default_symms.get());
  }

  ExprPtr operator()(const io::serialization::v1::ast::Sum &sum) const {
    return ast_to_expr<PositionCache>(sum, position_cache.get(), begin.get(),
                                      default_symms.get());
  }

  /// reports @p message at the source range @p node came from
  template <typename AST>
  [[noreturn]] void throw_at(const AST &node, std::string message) const {
    auto [offset, length] = get_pos(node, position_cache.get(), begin.get());
    throw SerializationError(offset, length, std::move(message));
  }

  /// splits the trailing state marks off the name of @p node ; the marks may
  /// come in either order, at most one of each
  /// @return the bare name, whether it was adjointed, whether it was
  ///         K-conjugated
  template <typename AST>
  std::tuple<std::wstring, bool, bool> split_marks(const AST &node) const {
    std::wstring name = node.name;
    try {
      const auto [adjointed, kconjugated] = split_state_marks(name);
      return {std::move(name), adjointed, kconjugated};
    } catch (const Exception &e) {
      throw_at(node, e.what());
    }
  }

  /// refuses a state on a name that admits none
  template <typename AST>
  void refuse_marks(const AST &node, bool adjointed, bool kconjugated) const {
    if (adjointed || kconjugated)
      throw_at(node, "an operator name carries no adjoint or conjugation mark");
  }

  ExprPtr operator()(const io::serialization::v1::ast::Tensor &tensor) const {
    auto [braIndices, ketIndices, auxiliaries] =
        make_indices(tensor.indices, position_cache.get(), begin.get());

    auto [perm_symm, braket_symm, column_symm, parity] =
        to_symmetries(tensor.symmetry, default_symms.get(),
                      position_cache.get(), begin.get());

    // the two core states are spelled as trailing marks of the name; they are
    // split off here and applied after construction rather than left to the
    // Tensor constructor's own mark adoption, which cannot hand back the sign
    // their normalization can carry
    auto [name, adjointed, kconjugated] = split_marks(tensor);

    // create NormalOperator or Tensor
    decltype(ranges::begin(FNOperator::labels())) fit;
    if ((fit = ranges::find(FNOperator::labels(), name)) !=
        ranges::end(FNOperator::labels())) {
      // an operator-valued tensor carries neither state: its bra<->ket swap
      // exchanges creators and annihilators
      refuse_marks(tensor, adjointed, kconjugated);
      SEQUANT_ASSERT(ranges::size(auxiliaries) == 0);
      SEQUANT_ASSERT(!tensor.symmetry.has_value() ||
                     ((tensor.symmetry.value().perm_symm ==
                           ast::SymmetrySpec::unspecified ||
                       tensor.symmetry.value().perm_symm == 'A') &&
                      (tensor.symmetry.value().column_symm ==
                           ast::SymmetrySpec::unspecified ||
                       tensor.symmetry.value().column_symm == 'S')));
      Vacuum vac = fit == ranges::begin(FNOperator::labels())
                       ? Vacuum::Physical
                       : tilde_vacuum(Statistics::FermiDirac);
      return ex<FNOperator>(cre(std::move(ketIndices)),
                            ann(std::move(braIndices)), vac);
    }
    decltype(ranges::begin(BNOperator::labels())) bit;
    if ((bit = ranges::find(BNOperator::labels(), name)) !=
        ranges::end(BNOperator::labels())) {
      refuse_marks(tensor, adjointed, kconjugated);
      SEQUANT_ASSERT(ranges::size(auxiliaries) == 0);
      SEQUANT_ASSERT(!tensor.symmetry.has_value() ||
                     ((tensor.symmetry.value().perm_symm ==
                           ast::SymmetrySpec::unspecified ||
                       tensor.symmetry.value().perm_symm == 'S') &&
                      (tensor.symmetry.value().column_symm ==
                           ast::SymmetrySpec::unspecified ||
                       tensor.symmetry.value().column_symm == 'S')));
      Vacuum vac = bit == ranges::begin(BNOperator::labels())
                       ? Vacuum::Physical
                       : tilde_vacuum(Statistics::BoseEinstein);
      return ex<BNOperator>(cre(std::move(ketIndices)),
                            ann(std::move(braIndices)), vac);
    }

    // whether the input spelled a symmetry out, as opposed to inheriting it
    // from the Context/DeserializationOptions defaults. Only a *spelled-out*
    // value may contradict a symmetry that is fixed by definition -- a
    // defaulted one is silently replaced by the correct value below, since the
    // ambient default says nothing about this particular tensor.
    const auto specified = [&tensor](char ast::SymmetrySpec::*field) {
      return tensor.symmetry.has_value() &&
             tensor.symmetry.value().*field != ast::SymmetrySpec::unspecified;
    };
    const bool perm_symm_specified = specified(&ast::SymmetrySpec::perm_symm);
    const bool braket_symm_specified =
        specified(&ast::SymmetrySpec::braket_symm);
    const bool column_symm_specified =
        specified(&ast::SymmetrySpec::column_symm);

    // (Anti)symmetry in bra and ket implies column symmetry, and the Tensor
    // ctor rejects an *explicit* ColumnSymmetry that contradicts it. Apply the
    // implication here for a column symmetry that merely came from the
    // defaults, so that only a genuinely contradicting explicit spec reaches
    // (and is rejected by) the ctor.
    if (!column_symm_specified &&
        (perm_symm == Symmetry::Symm || perm_symm == Symmetry::Antisymm))
      column_symm = ColumnSymmetry::Symm;

    // Force the defining symmetries of the reserved (anti)symmetrization
    // operators; see sequant::{anti,}symmetrizer_symmetries.
    const bool is_reserved_symmetrizer =
        name == reserved::antisymm_label() || name == reserved::symm_label();
    // Â antisymmetrizes within bra and within ket, Ŝ only across the
    // {bra,ket} particle columns (i.e. it is perm-Nonsymm). Supply the
    // defining value only when none was spelled out, so that a contradicting
    // explicit spec reaches the Tensor ctor and is rejected there rather than
    // silently overwritten here.
    if (is_reserved_symmetrizer && !perm_symm_specified)
      perm_symm = name == reserved::antisymm_label() ? Symmetry::Antisymm
                                                     : Symmetry::Nonsymm;
    // (anti)symmetrization operators act on indistinguishable particles, hence
    // are always column symmetric; supply that rather than passing the
    // Context's column default through, which the Tensor ctor would reject as
    // a contradicting *explicit* request
    if (is_reserved_symmetrizer && !column_symm_specified)
      column_symm = ColumnSymmetry::Symm;
    // likewise, force braket-Nonsymm rather than inheriting the Context's
    // default Hermiticity, which could derive a non-Nonsymm braket and make a
    // plain "Ŝ{...}"/"Â{...}" fail to construct; an explicit braket spec is
    // left untouched so that the Tensor ctor still rejects it
    if (is_reserved_symmetrizer && !braket_symm_specified)
      braket_symm = BraKetSymmetry::Nonsymm;

    // the reserved metric and Kronecker tensors are Hermitian by definition;
    // force that rather than inheriting the Context's default Hermiticity, so
    // that a deserialized s/δ equals the one make_overlap()/make_kronecker()
    // builds (they participate in the tensor hash, so a mismatch would keep
    // otherwise-equal terms from merging)
    if ((name == reserved::overlap_label() ||
         name == reserved::kronecker_label()) &&
        !braket_symm_specified)
      braket_symm = Hermiticity::Hermitian;

    // a reference density's symmetries are fixed by its label and rank (see
    // density::symmetries()); supply them where none was spelled out, and let
    // the Tensor ctor reject a spelled-out one that contradicts them. An
    // aux-only tensor is a layout representation, not a density.
    const bool is_density =
        ranges::contains(reserved::density_labels(), name) &&
        !(braIndices.empty() && ketIndices.empty());
    TensorSymmetries syms = is_density
                                ? density::symmetries(name, braIndices.size())
                                : TensorSymmetries{};
    if (!is_density || perm_symm_specified) syms.perm = perm_symm;
    if (!is_density || column_symm_specified) syms.column = column_symm;
    // the braket spec is a BraKetSymmetry or a Hermiticity; the parity, when
    // spelled out, is passed alongside it, and otherwise left unset so
    // Tensor::resolve_symmetries derives or defaults it the same way a
    // programmatic construction would
    if (!is_density || braket_symm_specified)
      std::visit(
          [&syms](auto symm) {
            using SymmType = std::decay_t<decltype(symm)>;
            if constexpr (std::is_same_v<SymmType, BraKetSymmetry>) {
              syms.braket = symm;
              syms.hermiticity = std::nullopt;
            } else {
              static_assert(std::is_same_v<SymmType, Hermiticity>);
              syms.hermiticity = symm;
            }
          },
          braket_symm);
    if (parity.has_value()) syms.conjugation_parity = *parity;

    ExprPtr t;
    try {
      t = ex<Tensor>(name, bra(std::move(braIndices)),
                     ket(std::move(ketIndices)), aux(std::move(auxiliaries)),
                     syms);
    } catch (const Exception &e) {
      // a symmetry the tensor cannot have (an exchange symmetry the traits do
      // not derive over the indices' basis, a contradicting column symmetry,
      // a density without its defining symmetries): report it at the tensor
      auto [offset, length] =
          get_pos(tensor, position_cache.get(), begin.get());
      throw SerializationError(offset, length, e.what());
    }
    if (adjointed || kconjugated) {
      const auto sign =
          t->template as<Tensor>().set_states(adjointed, kconjugated);
      // the normalization can consume a sign (an anti-Hermitian or
      // odd-parity tensor), which only a scalar factor can carry
      if (sign != 1) return ex<Product>(sign, ExprPtrList{std::move(t)});
    }
    return t;
  }

  ExprPtr operator()(
      const io::serialization::v1::ast::Variable &variable) const {
    // a Variable has the one conjugated state, spelled by a trailing `꙳`
    auto [name, adjointed, kconjugated] = split_marks(variable);
    if (adjointed)
      throw_at(variable, "a variable name carries no adjoint mark");
    ExprPtr var = ex<Variable>(std::move(name));
    if (kconjugated) var->as<Variable>().conjugate();
    return var;
  }

  ExprPtr operator()(const io::serialization::v1::ast::Number &number) const {
    return ex<Constant>(to_constant(number, position_cache.get(), begin.get()));
  }

  ExprPtr operator()(
      const io::serialization::v1::ast::RealImagPart &part) const {
    ExprPtr inner = ast_to_expr<PositionCache>(
        part.inner, position_cache.get(), begin.get(), default_symms.get());
    if (!inner) throw_at(part, "Re[]/Im[] wraps no expression");
    // the smart builders apply the eager composition rules, so what comes
    // back can be the inner expression itself (`Re[Re[x]]`), a constant
    // (`Re[3]`) or a scaled wrapper (`Re[1/2 x]`), exactly as a programmatic
    // real_part()/imaginary_part() call would give
    return part.imaginary ? imaginary_part(std::move(inner))
                          : real_part(std::move(inner));
  }

  ExprPtr operator()(const io::serialization::v1::ast::Power &power) const {
    // build base from Number or Variable
    ExprPtr base = boost::apply_visitor(*this, power.base);

    // exponent must be a real rational, reject otherwise
    if (power.exponent.imaginary) {
      auto [offset, length] = get_pos(power, position_cache.get(), begin.get());
      throw SerializationError(offset, length,
                               "Power exponent must be a real number");
    }
    if (static_cast<std::int64_t>(power.exponent.numerator) !=
            power.exponent.numerator ||
        static_cast<std::int64_t>(power.exponent.denominator) !=
            power.exponent.denominator) {
      auto [offset, length] = get_pos(power, position_cache.get(), begin.get());
      throw SerializationError(
          offset, length,
          "Power exponent must be a rational number (integer numerator and "
          "denominator)");
    }
    sequant::rational exponent(
        static_cast<std::int64_t>(power.exponent.numerator),
        static_cast<std::int64_t>(power.exponent.denominator));

    auto pw = ex<Power>(std::move(base), std::move(exponent));
    if (power.conjugated) {
      pw->as<Power>().conjugate();
    }
    return pw;
  }
};

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const io::serialization::v1::ast::NullaryValue &value,
                    const PositionCache &position_cache, const Iterator &begin,
                    DefaultSymmetries default_symms) {
  return boost::apply_visitor(
      Transformer<PositionCache, Iterator>{
          std::ref(position_cache), std::ref(begin), std::ref(default_symms)},
      value);
}

template <typename T, typename... Ts>
bool holds_alternative(const boost::variant<Ts...> &v) noexcept {
  return boost::get<T>(&v) != nullptr;
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const io::serialization::v1::ast::Product &product,
                    const PositionCache &position_cache, const Iterator &begin,
                    const DefaultSymmetries &default_symms) {
  if (product.factors.empty()) {
    SEQUANT_ABORT("Product factors must not be empty");
  }

  if (product.factors.size() == 1) {
    return ast_to_expr(product.factors.front(), position_cache, begin,
                       default_symms);
  }

  std::vector<ExprPtr> factors;
  factors.reserve(product.factors.size());
  Constant prefactor(1);

  // We perform constant folding
  for (const io::serialization::v1::ast::NullaryValue &value :
       product.factors) {
    if (holds_alternative<io::serialization::v1::ast::Number>(value)) {
      prefactor *=
          to_constant(boost::get<io::serialization::v1::ast::Number>(value),
                      position_cache, begin);
    } else {
      ExprPtr factor = ast_to_expr(value, position_cache, begin, default_symms);
      // a group that collapses to a constant joins the prefactor too: `(1 +
      // 2i)`, the spelling a composite scalar is emitted with, lands back on
      // the Product's scalar it was serialized from
      if (factor && factor->is<Constant>()) {
        prefactor *= factor->as<Constant>();
      } else {
        factors.push_back(std::move(factor));
      }
    }
  }

  if (factors.empty()) {
    // Only constants
    return ex<Constant>(std::move(prefactor));
  }

  if (factors.size() == 1 && prefactor.value() == 1) {
    return factors.front();
  }

  return ex<Product>(prefactor.value(), std::move(factors),
                     Product::Flatten::No);
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const io::serialization::v1::ast::Sum &sum,
                    const PositionCache &position_cache, const Iterator &begin,
                    const DefaultSymmetries &default_symms) {
  if (sum.summands.empty()) {
    return {};
  }
  if (sum.summands.size() == 1) {
    return ast_to_expr(sum.summands.front(), position_cache, begin,
                       default_symms);
  }

  std::vector<ExprPtr> summands;
  summands.reserve(sum.summands.size());
  std::transform(
      sum.summands.begin(), sum.summands.end(), std::back_inserter(summands),
      [&](const io::serialization::v1::ast::Product &product) {
        return ast_to_expr(product, position_cache, begin, default_symms);
      });

  ExprPtr folded = ex<Sum>(std::move(summands));
  // Sum::append adds constants up and drops zeros, so a written sum can
  // collapse; a collapsed one is its remaining term, not a one-summand Sum.
  // That is what makes `1 + 2i` -- the spelling of a composite scalar -- the
  // Constant it was serialized from
  const auto &folded_summands = folded->as<Sum>().summands();
  if (folded_summands.empty()) {
    return ex<Constant>(0);
  }
  if (folded_summands.size() == 1) {
    return folded_summands.front();
  }
  return folded;
}

template <typename PositionCache, typename Iterator>
ResultExpr ast_to_result(const io::serialization::v1::ast::ResultExpr &result,
                         const PositionCache &position_cache,
                         const Iterator &begin,
                         DefaultSymmetries default_symms) {
  ExprPtr lhs = std::visit(
      Transformer<PositionCache, Iterator>{
          std::ref(position_cache), std::ref(begin), std::ref(default_symms)},
      result.lhs);
  ExprPtr rhs = ast_to_expr(result.rhs, position_cache, begin, default_symms);

  if (lhs.is<Tensor>()) {
    return {std::move(lhs.as<Tensor>()), std::move(rhs)};
  } else if (lhs.is<Variable>()) {
    return {std::move(lhs.as<Variable>()), std::move(rhs)};
  } else {
    auto [offset, length] = get_pos(result.lhs, position_cache, begin);
    throw SerializationError(
        offset, length, "LHS of a ResultExpr must be a Tensor or a Variable");
  }
}

}  // namespace sequant::io::serialization::v1::transform

#endif  // SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP
