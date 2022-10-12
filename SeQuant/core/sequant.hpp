#ifndef SEQUANT_SEQUANT_H
#define SEQUANT_SEQUANT_H

#include "attr.hpp"
#include "index.hpp"
#include "space.hpp"

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
  explicit SeQuant(Vacuum vac, IndexSpaceMetric m = IndexSpaceMetric::Unit,
                   BraKetSymmetry bks = BraKetSymmetry::conjugate,
                   SPBasis spb = SPBasis::spinorbital)
      : vacuum_(vac), metric_(m), braket_symmetry_(bks), spbasis_(spb) {}
  ~SeQuant() = default;

  Vacuum vacuum() const { return vacuum_; }
  IndexSpaceMetric metric() const { return metric_; }
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }
  SPBasis spbasis() const { return spbasis_; }
  /// @return the IndexRegistry object
  std::shared_ptr<IndexRegistry> index_registry() const;

 private:
  Vacuum vacuum_ = Vacuum::Physical;
  IndexSpaceMetric metric_ = IndexSpaceMetric::Unit;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::conjugate;
  SPBasis spbasis_ = sequant::SPBasis::spinorbital;
};

/// SeQuant object equality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are equal
/// \warning does not compare index registries
bool operator==(const SeQuant& ctx1, const SeQuant& ctx2);

/// SeQuant object inequality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are not equal
/// \warning does not compare index registries
bool operator!=(const SeQuant& ctx1, const SeQuant& ctx2);

const SeQuant& get_default_context();
void set_default_context(const SeQuant& ctx);
void reset_default_context();

namespace detail {
struct ContextResetter {
  ContextResetter() noexcept;
  ~ContextResetter() noexcept;

 private:
  SeQuant previous_ctx_;
};
}  // namespace detail

detail::ContextResetter set_scoped_default_context(const SeQuant& ctx);

}  // namespace sequant

#endif
