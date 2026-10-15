#ifndef SEQUANT_DOMAIN_MBPT_SERIALIZATION_HPP
#define SEQUANT_DOMAIN_MBPT_SERIALIZATION_HPP

namespace sequant::mbpt {

/// Registers the v2 serialization functions for mbpt::Operator<qns_t, S>:
/// `~fmbop(label, ranks...)` for fermionic and `~bmbop(label, ranks...)` for
/// bosonic operators. The ranks are those of the factory that makes the
/// operator: `K` for a general operator (e.g. `~fmbop(g, 2)`) and the
/// quasiparticle ranks `nₚ, nₕ` for an (de)excitation operator, written as a
/// single `K` if they are equal (e.g. `~fmbop(t, 2)`, `~fmbop(R, 2, 1)`). A
/// trailing `batch=[...]` gives the batch ordinals. The label carries the
/// adjoint mark and the perturbation order, if any (e.g. `~fmbop(t¹⁺, 2)`).
///
/// An operator is rebuilt from this through the operator registry of the
/// default mbpt::Context, as the standard factories (mbpt::t(), mbpt::R(),
/// mbpt::A(), ...) build it, so the registry must define its label when
/// deserializing. Serializing an operator that the standard factories cannot
/// reproduce (e.g. one acting on custom index spaces) throws a
/// io::serialization::SerializationError.
///
/// Registering again is harmless.
void register_serialization_functions();

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SERIALIZATION_HPP
