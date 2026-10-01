//
// Amplitude tensors and the check that they carry their basis grants.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/basis_grants.hpp>
#include <SeQuant/domain/mbpt/context.hpp>

#include <string>

namespace sequant::mbpt {

namespace {
std::string instance_text(const IndexBasis::optional_instance& instance) {
  return instance ? std::to_string(*instance) : std::string("none");
}
}  // namespace

bool is_amplitude_tensor(const AbstractTensor& t, const OpRegistry& reg) {
  const std::wstring label(strip_adjoint_label(t._label()));
  if (!reg.contains(label)) return false;
  const auto cls = reg.to_class(label);
  return cls == OpClass::Ex || cls == OpClass::Deex;
}

void assert_amplitudes_carry_granted_basis(const ExprPtr& expr) {
  const auto& reg = *get_default_mbpt_context().op_registry();
  container::set<std::wstring> granted, ungranted;
  auto check_slot = [&reg](const AbstractTensor& t, const std::wstring& label,
                           const Index& idx) {
    const auto grant = reg.basis_grant(label, idx.space());
    const auto& instance = idx.basis().basis_instance();
    if (instance == grant || (!grant && !idx.has_proto_indices())) return;
    throw Exception("mbpt::assert_amplitudes_carry_granted_basis: slot " +
                    toUtf8(idx.full_label()) + " of " +
                    toUtf8(std::wstring(t._label())) +
                    " carries basis instance " + instance_text(instance) +
                    " but the grant of " + toUtf8(label) + " on its space is " +
                    instance_text(grant));
  };
  expr->visit(
      [&](const ExprPtr& x) {
        if (!x->is<AbstractTensor>()) return;
        const auto& t = x->as<AbstractTensor>();
        if (!is_amplitude_tensor(t, reg)) return;
        const std::wstring label(strip_adjoint_label(t._label()));
        (reg.has_basis_grants(label) ? granted : ungranted).insert(label);
        for (const Index& idx : t._bra()) check_slot(t, label, idx);
        for (const Index& idx : t._ket()) check_slot(t, label, idx);
      },
      /* atoms_only = */ true);
  if (!granted.empty() && !ungranted.empty())
    throw Exception("mbpt::assert_amplitudes_carry_granted_basis: amplitude " +
                    toUtf8(*ungranted.begin()) + " has no basis grant while " +
                    toUtf8(*granted.begin()) +
                    " has one; grant every amplitude of the expression");
}

}  // namespace sequant::mbpt
