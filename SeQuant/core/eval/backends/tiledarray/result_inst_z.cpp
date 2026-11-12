#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/backends/tiledarray/tiledarray.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY

#include <complex>

namespace sequant {

template class SEQUANT_EVAL_TA_TEMPLATE_DEF_API ResultTensorTA<
    TA::DistArray<TA::Tensor<std::complex<double>>, TA::SparsePolicy>>;

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY
