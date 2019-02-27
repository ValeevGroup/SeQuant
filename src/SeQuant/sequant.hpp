#ifndef SEQUANT_SEQUANT_H
#define SEQUANT_SEQUANT_H

#include "space.hpp"
#include "attr.hpp"

namespace sequant {

/// Specifies second quantization context, such as vacuum choice, whether index
/// spaces are orthonormal, etc.
class SeQuant {
 public:
  /// @param vac the vacuum type
  explicit SeQuant(Vacuum vac = Vacuum::SingleProduct,
                   IndexSpaceMetric m = IndexSpaceMetric::Unit)
      : vacuum_(vac), metric_(m) {}
  ~SeQuant() = default;

  Vacuum vacuum() const { return vacuum_; }
  IndexSpaceMetric metric() const { return metric_; }

 private:
  Vacuum vacuum_;
  IndexSpaceMetric metric_;
};

const SeQuant &get_default_context();
void set_default_context(const SeQuant &ctx);
void reset_default_context();

}  // namespace sequant

#endif