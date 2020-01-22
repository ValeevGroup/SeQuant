#ifndef SEQUANT_EVAL_FWD
#define SEQUANT_EVAL_FWD

#include <SeQuant/core/container.hpp>
#include <boost/functional/hash_fwd.hpp>

#include <cstddef>
#include <complex>
#include <memory>
#include <string>

namespace sequant::evaluate {
class EvalTensor;

template <typename DataTensorType>
class ExprEvaluator;

using EvTensorPtr = std::shared_ptr<EvalTensor>;

using label_container_type = container::svector<std::wstring, 4>;

using constant_type = std::complex<double>;

using hash_type = std::size_t;

using btas_index_container = container::svector<hash_type>;

template <typename DataTensorType>
using hash_to_dtensor_map =
    container::map<hash_type, std::shared_ptr<DataTensorType>>;

}  // namespace sequant::evaluate

#endif  // ifndef SEQUANT_EVAL_FWD
