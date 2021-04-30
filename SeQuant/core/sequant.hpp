#ifndef SEQUANT_SEQUANT_H
#define SEQUANT_SEQUANT_H

#include "space.hpp"
#include "attr.hpp"
#include "index.hpp"

namespace sequant {

/// Specifies second quantization context, such as vacuum choice, whether index
/// spaces are orthonormal, sizes of index spaces, etc.
class SeQuant {
 public:
  SeQuant() = default;
  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  explicit SeQuant(Vacuum vac,
                   IndexSpaceMetric m,
                   BraKetSymmetry bks,
                   Spinbasis spb = sequant::Spinbasis::spin_orbit)
      : vacuum_(vac), metric_(m), braket_symmetry_(bks),spin_basis_(spb) {}
  ~SeQuant() = default;

  Vacuum vacuum() const { return vacuum_; }
  IndexSpaceMetric metric() const { return metric_; }
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }
  Spinbasis spin_basis() const { return spin_basis_;}
  /// @return the IndexRegistry object
  std::shared_ptr<IndexRegistry> index_registry() const;

 private:
  Vacuum vacuum_ = Vacuum::SingleProduct;
  IndexSpaceMetric metric_ = IndexSpaceMetric::Unit;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::conjugate;
  Spinbasis spin_basis_ = sequant::Spinbasis::spin_orbit;
};

const SeQuant &get_default_context();
void set_default_context(const SeQuant &ctx);
void reset_default_context();

}  // namespace sequant

#endif