#include <SeQuant/core/eval/backends/tiledarray/result.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY

#include <tiledarray.h>

#include <complex>

namespace sequant {

template class SEQUANT_EVAL_TA_TEMPLATE_DEF_API
    ResultTensorOfTensorTA<TA::DistArray<
        TA::Tensor<TA::ArenaTensor<std::complex<double>>>, TA::SparsePolicy>>;

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY
