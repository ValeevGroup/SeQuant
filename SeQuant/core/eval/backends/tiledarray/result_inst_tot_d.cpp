#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/backends/tiledarray/tiledarray.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY

namespace sequant {

template class SEQUANT_EVAL_TA_TEMPLATE_DEF_API ResultTensorOfTensorTA<
    TA::DistArray<TA::Tensor<TA::ArenaTensor<double>>, TA::SparsePolicy>>;

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY
