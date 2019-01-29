#ifndef SEQUANT2_SEQUANT_H
#define SEQUANT2_SEQUANT_H

#include "space.hpp"
#include "vacuum.hpp"

namespace sequant2 {

/// Specifies second quantization context, such as vacuum choice, whether index
/// spaces are orthonormal, etc.
class SeQuant2 {
 public:
  /// @param vac the vacuum type
  explicit SeQuant2(Vacuum vac = Vacuum::SingleProduct,
                    IndexSpaceMetric m = IndexSpaceMetric::Unit)
      : vacuum_(vac), metric_(m) {}
  ~SeQuant2() = default;

  Vacuum vacuum() const { return vacuum_; }
  IndexSpaceMetric metric() const { return metric_; }

 private:
  Vacuum vacuum_;
  IndexSpaceMetric metric_;
};

const SeQuant2 &get_default_context();
void set_default_context(const SeQuant2 &ctx);
void reset_default_context();

}  // namespace sequant2

#endif