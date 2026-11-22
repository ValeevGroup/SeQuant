#include <SeQuant/domain/mbpt/models/hf.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <string>
#include <string_view>
#include <utility>

namespace sequant::mbpt::hf {

namespace {

constexpr std::wstring_view rotation_label = L"κₒ";

[[noreturn]] void reject(std::string_view function, const std::string& text) {
  throw Exception("mbpt::hf::" + std::string(function) + ": " + text);
}

/// the reference-occupied space and its complement, after the checks every
/// function here makes
std::pair<IndexSpace, IndexSpace> spaces(std::string_view function) {
  if (!get_default_mbpt_context().op_registry()->contains(
          std::wstring(rotation_label)))
    reject(function, "the operator label κₒ is not registered");
  if (get_default_context().vacuum() != Vacuum::SingleProduct ||
      !detail::reference_is_vacuum())
    reject(function,
           "the vacuum must be SingleProduct and the reference the vacuum");
  if (get_default_mbpt_context().csv() == CSV::Yes)
    reject(function,
           "orbital rotations with pair-specific orbitals are not supported");
  // the reference being the vacuum, the vacuum-unoccupied space is the
  // complement of the reference-occupied one
  const auto ibr = get_default_context().index_basis_registry();
  const auto occupied = ibr->reference_occupied_space(Spin::any);
  try {
    return {occupied, ibr->vacuum_unoccupied_space(Spin::any)};
  } catch (const Exception&) {
    reject(function, "the unoccupied space is not registered");
  }
}

ExprPtr rotated(std::size_t order, const std::pair<IndexSpace, IndexSpace>& s) {
  const auto K = ex<op_t>(
      [] { return rotation_label; },
      [s] {
        return OpMaker<Statistics::FermiDirac>(std::wstring(rotation_label),
                                               ncre(1), nann(1), cre(s.second),
                                               ann(s.first))();
      },
      [s](qnc_t& qns) {
        qns = combine(generic_excitation_qns(1, 1, s.second, s.first), qns);
      });
  return simplify(op::ref_av(lst(op::H(), K - adjoint(K), order)));
}

}  // namespace

ExprPtr rotated_energy(std::size_t order) {
  return rotated(order, spaces("rotated_energy"));
}

}  // namespace sequant::mbpt::hf
