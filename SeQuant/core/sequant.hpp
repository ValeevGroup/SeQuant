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
  /// @param spb single-particle basis (spin-free or spin-dependent)
  explicit SeQuant(Vacuum vac,
                   IndexSpaceMetric m,
                   BraKetSymmetry bks,
                   SPBasis spb = sequant::SPBasis::spinorbital)
      : vacuum_(vac), metric_(m), braket_symmetry_(bks), spbasis_(spb) {}
  ~SeQuant() = default;

  Vacuum vacuum() const { return vacuum_; }
  IndexSpaceMetric metric() const { return metric_; }
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }
  SPBasis spbasis() const { return spbasis_;}
  /// @return the IndexRegistry object
  std::shared_ptr<IndexRegistry> index_registry() const;

 private:
  Vacuum vacuum_ = Vacuum::SingleProduct;
  IndexSpaceMetric metric_ = IndexSpaceMetric::Unit;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::conjugate;
  SPBasis spbasis_ = sequant::SPBasis::spinorbital;
};

const SeQuant &get_default_context();
void set_default_context(const SeQuant &ctx);
void reset_default_context();

}  // namespace sequant

#endif