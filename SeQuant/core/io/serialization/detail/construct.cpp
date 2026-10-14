#include <SeQuant/core/context.hpp>
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/detail/construct.hpp>
#include <SeQuant/core/reserved.hpp>

#include <range/v3/algorithm/contains.hpp>

namespace sequant::io::serialization::detail {

DefaultSymmetries resolve_default_symmetries(
    const DeserializationOptions &options) {
  // unspecified symmetries default to the active Context's; the Hermiticity is
  // resolved against each tensor's base_field downstream, as in the Tensor ctor
  const Context &ctx = get_default_context();
  DefaultSymmetries symms{ctx.deserialization_symmetry(),
                          ctx.deserialization_hermiticity(),
                          ctx.deserialization_column_symmetry()};

  if (options.def_perm_symm.has_value()) {
    std::get<0>(symms) = options.def_perm_symm.value();
  }
  if (options.def_braket_symm.has_value()) {
    std::get<1>(symms) = options.def_braket_symm.value();
  }
  if (options.def_col_symm.has_value()) {
    std::get<2>(symms) = options.def_col_symm.value();
  }
  if (std::get<0>(symms) != Symmetry::Nonsymm) {
    // antisymmetry/symmetric bra and ket imply particle symmetry
    std::get<2>(symms) = ColumnSymmetry::Symm;
  }

  return symms;
}

Vacuum tilde_vacuum(Statistics s) {
  const auto vac = get_default_context(s).vacuum();
  return vac == Vacuum::Physical ? Vacuum::SingleProduct : vac;
}

ExprPtr make_tensor(const std::wstring &name, container::vector<Index> bra,
                    container::vector<Index> ket, container::vector<Index> aux,
                    Symmetry perm,
                    std::variant<BraKetSymmetry, Hermiticity> braket,
                    ColumnSymmetry column, bool perm_specified,
                    bool braket_specified, bool column_specified) {
  // (Anti)symmetry in bra and ket implies column symmetry, and the Tensor
  // ctor rejects an *explicit* ColumnSymmetry that contradicts it. Apply the
  // implication here for a column symmetry that merely came from the defaults,
  // so that only a genuinely contradicting explicit spec reaches (and is
  // rejected by) the ctor.
  if (!column_specified &&
      (perm == Symmetry::Symm || perm == Symmetry::Antisymm))
    column = ColumnSymmetry::Symm;

  // Force the defining symmetries of the reserved (anti)symmetrization
  // operators; see sequant::{anti,}symmetrizer_symmetries.
  const bool is_reserved_symmetrizer =
      name == reserved::antisymm_label() || name == reserved::symm_label();
  // Â antisymmetrizes within bra and within ket, Ŝ only across the {bra,ket}
  // particle columns (i.e. it is perm-Nonsymm). Supply the defining value only
  // when none was spelled out, so that a contradicting explicit spec reaches
  // the Tensor ctor and is rejected there rather than silently overwritten.
  if (is_reserved_symmetrizer && !perm_specified)
    perm = name == reserved::antisymm_label() ? Symmetry::Antisymm
                                              : Symmetry::Nonsymm;
  // (anti)symmetrization operators act on indistinguishable particles, hence
  // are always column symmetric; supply that rather than passing the Context's
  // column default through, which the Tensor ctor would reject as a
  // contradicting *explicit* request
  if (is_reserved_symmetrizer && !column_specified)
    column = ColumnSymmetry::Symm;
  // likewise, force braket-Nonsymm rather than inheriting the Context's default
  // Hermiticity, which could derive a non-Nonsymm braket and make a plain
  // "Ŝ{...}"/"Â{...}" fail to construct; an explicit braket spec is left
  // untouched so that the Tensor ctor still rejects it
  if (is_reserved_symmetrizer && !braket_specified)
    braket = BraKetSymmetry::Nonsymm;

  // the reserved metric and Kronecker tensors are Hermitian by definition;
  // force that rather than inheriting the Context's default Hermiticity, so
  // that a deserialized s/δ equals the one make_overlap()/make_kronecker()
  // builds (they participate in the tensor hash, so a mismatch would keep
  // otherwise-equal terms from merging)
  if ((name == reserved::overlap_label() ||
       name == reserved::kronecker_label()) &&
      !braket_specified)
    braket = Hermiticity::Hermitian;

  // a reference density's symmetries are fixed by its label and rank (see
  // density::symmetries()); supply them where none was spelled out, and let the
  // Tensor ctor reject a spelled-out one that contradicts them. An aux-only
  // tensor is a layout representation, not a density.
  if (ranges::contains(reserved::density_labels(), name) &&
      !(bra.empty() && ket.empty())) {
    auto syms = density::symmetries(name, bra.size());
    if (perm_specified) syms.perm = perm;
    if (column_specified) syms.column = column;
    if (braket_specified)
      std::visit(
          [&syms](auto symm) {
            if constexpr (std::is_same_v<decltype(symm), BraKetSymmetry>) {
              syms.braket = symm;
              syms.hermiticity = std::nullopt;
            } else
              syms.hermiticity = symm;
          },
          braket);
    return ex<Tensor>(name, ::sequant::bra(std::move(bra)),
                      ::sequant::ket(std::move(ket)),
                      ::sequant::aux(std::move(aux)), syms);
  }

  // Dispatch to correct Tensor constructor (taking either BraKetSymmetry or
  // Hermiticity)
  return std::visit(
      [&](auto symm) {
        return ex<Tensor>(name, ::sequant::bra(std::move(bra)),
                          ::sequant::ket(std::move(ket)),
                          ::sequant::aux(std::move(aux)), perm, symm, column);
      },
      braket);
}

}  // namespace sequant::io::serialization::detail
