#include <SeQuant/core/eval/backends/tiledarray/result.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY

#include <tiledarray.h>

namespace sequant {

template class ResultTensorOfTensorTA<
    TA::DistArray<TA::Tensor<TA::Tensor<double>>, TA::SparsePolicy>>;

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY
