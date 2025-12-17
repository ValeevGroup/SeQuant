#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/swap.hpp>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/detail/variant.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/interface.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace sequant::mbpt {

namespace detail {

Index make_index_with_spincase(const Index& idx, mbpt::Spin s) {
  // sanity check: make sure have only one spin label
  SEQUANT_ASSERT(!(idx.label().find(L'↑') != std::wstring::npos &&
                   idx.label().find(L'↓') != std::wstring::npos));

  // to preserve rest of bits first unset spin bit, then set them to the desired
  // state
  auto qns = mbpt::spinannotation_remove(idx.space().qns()).unIon(s);

  IndexSpace space;
  // try looking up space in registry
  const auto label = mbpt::spinannotation_replacе(idx.space().base_key(), s);
  if (auto isr = get_default_context().index_space_registry()) {
    auto* space_ptr = isr->retrieve_ptr(label);
    if (space_ptr && space_ptr->type() == idx.space().type() &&
        space_ptr->qns() == qns) {
      space = *space_ptr;
    }
  }
  // if space not found, construct
  if (!space) {
    space = IndexSpace{label, idx.space().type(), qns,
                       // N.B. assume size does not depend on spin
                       idx.space().approximate_size()};
  }
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_spincase(pidx, s);
  return Index{space, idx.ordinal(), protoindices};
}

// The argument really should be non-const but const semantics are broken
// for the ExprPtr type so we are required to make this const in order
// to be able to use this function everywhere we want to.
void reset_idx_tags(const ExprPtr& expr) {
  expr->visit(
      [](ExprPtr& current) {
        if (current.is<AbstractTensor>()) {
          current.as<AbstractTensor>()._reset_tags();
        }
      },
      true);
}

template <typename Container, typename TraceFunction, typename... Args>
[[nodiscard]] Container wrap_trace(const ResultExpr& expr,
                                   TraceFunction&& tracer, Args&&... args) {
  bool searchForNonEquivalentResults = expr.symmetry() != Symmetry::Nonsymm;
  searchForNonEquivalentResults &=
      expr.bra().size() > 1 || expr.ket().size() > 1;
  const bool brasSameSpace = std::all_of(
      expr.bra().begin(), expr.bra().end(),
      [&](const Index& idx) { return idx.space() == expr.bra()[0].space(); });
  const bool ketsSameSpace = std::all_of(
      expr.ket().begin(), expr.ket().end(),
      [&](const Index& idx) { return idx.space() == expr.ket()[0].space(); });
  searchForNonEquivalentResults &= !brasSameSpace && !ketsSameSpace;

  if (!searchForNonEquivalentResults) {
    ResultExpr traced = expr.clone();
    traced.expression() =
        tracer(traced.expression(),
               traced.index_particle_grouping<container::svector<Index>>(),
               std::forward<Args>(args)...);

    traced.set_symmetry(Symmetry::Nonsymm);

    return {std::move(traced)};
  }

  SEQUANT_ASSERT(expr.symmetry() == Symmetry::Antisymm ||
                 expr.symmetry() == Symmetry::Symm);

  // TODO: Do we have to track the sign?
  const bool permuteBra = expr.bra().size() >= expr.ket().size();
  auto permIndices = permuteBra ? expr.bra() : expr.ket();
  const std::size_t unchangedSize =
      permuteBra ? expr.ket().size() : expr.bra().size();

  [[maybe_unused]] auto get_phase = [](auto container) {
    reset_ts_swap_counter<Index>();
    bubble_sort(container.begin(), container.end(), std::less<Index>{});
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  reset_ts_swap_counter<Index>();
  bubble_sort(permIndices.begin(), permIndices.end(), std::less<Index>{});
  const int initialSign = ts_swap_counter_is_even<Index>() ? 1 : -1;
  const auto originalIndices = permIndices;

  container::svector<container::set<std::pair<IndexSpace, IndexSpace>>>
      idxPairings;

  Container resultSet;

  // For next_permutation to work in this context, permIndices must be sorted
  SEQUANT_ASSERT(std::is_sorted(permIndices.begin(), permIndices.end()));

  int sign = initialSign;
  do {
    const int currentSign = sign;
    // std::next_permutation creates one lexicographical permutation after the
    // other, which should imply that the phase should alternate between
    // iterations.
    sign *= -1;
    SEQUANT_ASSERT(currentSign == get_phase(permIndices) * initialSign);

    container::set<std::pair<IndexSpace, IndexSpace>> currentPairing;

    for (std::size_t i = 0; i < unchangedSize; ++i) {
      if (permuteBra) {
        currentPairing.insert(
            std::make_pair(permIndices[i].space(), expr.ket()[i].space()));
      } else {
        currentPairing.insert(
            std::make_pair(expr.bra()[i].space(), permIndices[i].space()));
      }
    }

    for (std::size_t i = unchangedSize; i < permIndices.size(); ++i) {
      currentPairing.insert(
          std::make_pair(permIndices[i].space(), IndexSpace::null));
    }

    if (std::find(idxPairings.begin(), idxPairings.end(), currentPairing) !=
        idxPairings.end()) {
      continue;
    }

    // Found a new index pairing

    ExprPtr expression = expr.expression().clone();

    expression *= ex<Constant>(currentSign);
    expression = simplify(expression);

    ResultExpr result = [&]() {
      SEQUANT_ASSERT(expr.has_label());
      if (permuteBra) {
        return ResultExpr(bra(permIndices), ket(expr.ket()), aux(expr.aux()),
                          expr.symmetry(), expr.braket_symmetry(),
                          expr.column_symmetry(), expr.label(),
                          std::move(expression));
      } else {
        return ResultExpr(bra(expr.bra()), ket(permIndices), aux(expr.aux()),
                          expr.symmetry(), expr.braket_symmetry(),
                          expr.column_symmetry(), expr.label(),
                          std::move(expression));
      }
    }();

    result.expression() =
        tracer(result.expression(),
               result.index_particle_grouping<container::svector<Index>>(),
               std::forward<Args>(args)...);

    result.set_symmetry(Symmetry::Nonsymm);

    resultSet.push_back(std::move(result));
  } while (std::next_permutation(permIndices.begin(), permIndices.end()));

  return resultSet;
}

}  // namespace detail

Index make_spinalpha(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::alpha);
};

Index make_spinbeta(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::beta);
};

Index make_spinfree(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::any);
};

ExprPtr swap_bra_ket(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Lambda for tensor
  auto tensor_swap = [](const Tensor& tensor) {
    return ex<Tensor>(tensor.label(), bra(tensor.ket().value()),
                      ket(tensor.bra().value()), tensor.symmetry(),
                      tensor.braket_symmetry(), tensor.column_symmetry());
  };

  // Lambda for product
  auto product_swap = [&tensor_swap](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        result->append(1, tensor_swap(term->as<Tensor>()),
                       Product::Flatten::No);
      } else if (term->is<Variable>() || term->is<Constant>()) {
        result->append(1, term);
      } else {
        throw std::runtime_error("Invalid Expr type in product_swap: " +
                                 term->type_name());
      }
    }
    return result;
  };

  if (expr->is<Tensor>())
    return tensor_swap(expr->as<Tensor>());
  else if (expr->is<Product>())
    return product_swap(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      result->append(swap_bra_ket(term));
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in swap_bra_ket: " +
                             expr->type_name());
  }
}

ExprPtr append_spin(const ExprPtr& expr,
                    const container::map<Index, Index>& index_replacements) {
  auto add_spin_to_tensor = [&index_replacements](const Tensor& tensor) {
    auto spin_tensor = std::make_shared<Tensor>(tensor);
    spin_tensor->transform_indices(index_replacements);
    return spin_tensor;
  };

  auto add_spin_to_product = [&add_spin_to_tensor](const Product& product) {
    auto spin_product = std::make_shared<Product>();
    spin_product->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        spin_product->append(1, add_spin_to_tensor(term->as<Tensor>()));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        spin_product->append(1, term);
      } else {
        throw std::runtime_error(
            "Invalid Expr type in append_spin::add_spin_to_product: " +
            term->type_name());
      }
    }
    return spin_product;
  };

  if (expr->is<Tensor>()) {
    return add_spin_to_tensor(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return add_spin_to_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    auto spin_expr = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      spin_expr->append(append_spin(summand, index_replacements));
    }
    return spin_expr;
  } else if (expr->is<Constant>() || expr->is<Variable>()) {
    return expr;
  }

  throw std::runtime_error("Unsupported Expr type in append_spin");
}

ExprPtr remove_spin(const ExprPtr& expr) {
  auto remove_spin_from_tensor = [](const Tensor& tensor) {
    container::svector<Index> b(tensor.bra().begin(), tensor.bra().end());
    container::svector<Index> k(tensor.ket().begin(), tensor.ket().end());
    {
      for (auto&& idx : ranges::views::concat(b, k)) {
        idx = make_spinfree(idx);
      }
    }
    return ex<Tensor>(tensor.label(), bra(std::move(b)), ket(std::move(k)),
                      tensor.aux(), tensor.symmetry(),
                      tensor.braket_symmetry());
  };

  auto remove_spin_from_product =
      [&remove_spin_from_tensor](const Product& product) {
        auto result = std::make_shared<Product>();
        result->scale(product.scalar());
        for (auto&& term : product) {
          if (term->is<Tensor>()) {
            result->append(1, remove_spin_from_tensor(term->as<Tensor>()));
          } else if (term->is<Constant>() || term->is<Variable>()) {
            result->append(1, term);
          } else {
            throw std::runtime_error(
                "Invalid Expr type in remove_spin::remove_spin_from_product: " +
                term->type_name());
          }
        }
        return result;
      };

  if (expr->is<Tensor>()) {
    return remove_spin_from_tensor(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return remove_spin_from_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(remove_spin(summand));
    }
    return result;
  } else if (expr->is<Constant>() || expr->is<Variable>()) {
    return expr;
  } else {
    throw std::runtime_error("Invalid Expr type in remove_spin: " +
                             expr->type_name());
  }
}

bool ms_conserving_columns(const AbstractTensor& tensor) {
  for (const auto& [bra, ket] :
       ranges::zip_view(tensor._bra(), tensor._ket())) {
    if (bra.nonnull() && ket.nonnull()) {
      const auto bra_ms = mbpt::to_spin(bra.space().qns());
      ;
      const auto ket_ms = mbpt::to_spin(ket.space().qns());
      ;
      if (bra_ms != ket_ms) return false;
    }
  }
  return true;
}

bool ms_uniform_tensor(const AbstractTensor& tensor) {
  auto braket = tensor._braket();
  SEQUANT_ASSERT(ranges::empty(braket) == false);
  std::optional<mbpt::Spin> ms;
  return ranges::all_of(braket, [&ms](const auto& idx) {
    if (idx.nonnull()) {
      const auto idx_ms = mbpt::to_spin(idx.space().qns());
      if (ms.has_value())
        return idx_ms == *ms;
      else {
        ms = idx_ms;
        return true;
      }
    } else
      return true;
  });
}

bool can_expand(const AbstractTensor& tensor) {
  SEQUANT_ASSERT(tensor._bra_rank() == tensor._ket_rank() &&
                 "can_expand(Tensor) failed.");
  if (tensor._bra_rank() != tensor._ket_rank()) return false;

  // indices must have specific spin
  [[maybe_unused]] auto all_have_spin =
      ranges::all_of(tensor._braket(), [](const auto& idx) {
        auto idx_spin = mbpt::to_spin(idx.space().qns());
        return idx_spin == mbpt::Spin::alpha || idx_spin == mbpt::Spin::beta;
      });
  SEQUANT_ASSERT(ranges::all_of(tensor._braket(), [](const auto& idx) {
    auto idx_spin = mbpt::to_spin(idx.space().qns());
    return idx_spin == mbpt::Spin::alpha || idx_spin == mbpt::Spin::beta;
  }));

  // count alpha indices in bra
  auto is_alpha = [](const Index& idx) {
    return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha;
  };

  // count alpha indices in bra
  auto a_bra = ranges::count_if(tensor._bra(), is_alpha);

  // count alpha indices in ket
  auto a_ket = ranges::count_if(tensor._ket(), is_alpha);

  return a_bra == a_ket;
}

ExprPtr expand_antisymm(const Tensor& tensor, bool skip_spinsymm) {
  SEQUANT_ASSERT(tensor.bra_rank() == tensor.ket_rank());
  // Return non-symmetric tensor if rank is 1
  if (tensor.bra_rank() <= 1) {
    Tensor new_tensor(tensor.label(), tensor.bra(), tensor.ket(), tensor.aux(),
                      Symmetry::Nonsymm, tensor.braket_symmetry(),
                      tensor.column_symmetry());
    return std::make_shared<Tensor>(new_tensor);
  }

  // If all indices have the same spin label,
  // return the antisymm tensor
  if (skip_spinsymm && ms_uniform_tensor(tensor)) {
    return std::make_shared<Tensor>(tensor);
  }

  SEQUANT_ASSERT(tensor.bra_rank() > 1 && tensor.ket_rank() > 1);

  auto get_phase = [](const Tensor& t) {
    container::svector<Index> bra(t.bra().begin(), t.bra().end());
    container::svector<Index> ket(t.ket().begin(), t.ket().end());
    reset_ts_swap_counter<Index>();
    bubble_sort(std::begin(bra), std::end(bra));
    bubble_sort(std::begin(ket), std::end(ket));
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  // Generate a sum of asymmetric tensors if the input tensor is antisymmetric
  // and greater than one body otherwise, return the tensor
  if (tensor.symmetry() == Symmetry::Antisymm) {
    const auto prefactor = get_phase(tensor);
    container::set<Index> bra_list(tensor.bra().begin(), tensor.bra().end());
    container::set<Index> ket_list(tensor.ket().begin(), tensor.ket().end());
    auto expr_sum = std::make_shared<Sum>();
    do {
      // N.B. must copy
      auto new_tensor =
          Tensor(tensor.label(), bra(bra_list), ket(ket_list), tensor.aux(),
                 Symmetry::Nonsymm, tensor.braket_symmetry(),
                 tensor.column_symmetry());

      if (ms_conserving_columns(new_tensor)) {
        auto new_tensor_product = std::make_shared<Product>();
        new_tensor_product->append(get_phase(new_tensor),
                                   ex<Tensor>(new_tensor));
        new_tensor_product->scale(prefactor);
        expr_sum->append(new_tensor_product);
      }
    } while (std::next_permutation(bra_list.begin(), bra_list.end()));

    return expr_sum;
  } else {
    return std::make_shared<Tensor>(tensor);
  }
}

ExprPtr expand_antisymm(const ExprPtr& expr, bool skip_spinsymm) {
  if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else if (expr->is<Tensor>())
    return expand_antisymm(expr->as<Tensor>(), skip_spinsymm);

  // Product lambda
  auto expand_product = [&skip_spinsymm](const Product& expr) {
    Product temp{};
    temp.scale(expr.scalar());
    for (auto&& term : expr) {
      if (term->is<Tensor>()) {
        temp.append(1, expand_antisymm(term->as<Tensor>(), skip_spinsymm),
                    Product::Flatten::No);
      } else if (term->is<Variable>() || term->is<Constant>()) {
        temp.append(1, term, Product::Flatten::No);
      } else {
        throw std::runtime_error(
            "Invalid Expr type in expand_antisymm::expand_product: " +
            term->type_name());
      }
    }
    ExprPtr result = std::make_shared<Product>(temp);
    rapid_simplify(result);
    return result;
  };

  if (expr->is<Product>())
    return expand_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      result->append(expand_antisymm(term, skip_spinsymm));
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in expand_antisymm: " +
                             expr->type_name());
  }
}

container::svector<container::map<Index, Index>> A_maps(const Tensor& A) {
  SEQUANT_ASSERT(A.label() == L"A");

  container::svector<std::size_t> bra_indices(A.bra_rank());
  container::svector<std::size_t> ket_indices(A.ket_rank());
  std::iota(bra_indices.begin(), bra_indices.end(), 0);
  std::iota(ket_indices.begin(), ket_indices.end(), 0);

  container::svector<container::map<Index, Index>> result;

  do {
    do {
      container::map<Index, Index> current_replacements;

      for (std::size_t i = 0; i < bra_indices.size(); ++i) {
        current_replacements.emplace(A.bra()[i], A.bra()[bra_indices[i]]);
      }
      for (std::size_t i = 0; i < ket_indices.size(); ++i) {
        current_replacements.emplace(A.ket()[i], A.ket()[ket_indices[i]]);
      }

      result.push_back(std::move(current_replacements));
    } while (std::next_permutation(bra_indices.begin(), bra_indices.end()));
  } while (std::next_permutation(ket_indices.begin(), ket_indices.end()));

  return result;
}

ExprPtr expand_A_op(const ProductPtr& product) {
  bool has_A_operator = false;

  // Check A and build replacement map
  container::svector<container::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto A = term->as<Tensor>();
      if (A.label() == L"A" && A.bra_rank() <= 1 && A.ket_rank() <= 1) {
        return remove_tensor(product, L"A");
      } else if ((A.label() == L"A")) {
        has_A_operator = true;
        map_list = A_maps(A);
        break;
      }
    }
  }

  if (!has_A_operator) return product;

  auto new_result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    // Get phase of the transformation
    int phase;
    {
      container::svector<Index> transformed_list;
      for (const auto& [key, val] : map) transformed_list.push_back(val);

      reset_ts_swap_counter<Index>();
      bubble_sort(std::begin(transformed_list), std::end(transformed_list));
      phase = ts_swap_counter_is_even<Index>() ? 1 : -1;
    }

    ProductPtr new_product = std::make_shared<Product>();
    new_product->scale(product->scalar());
    auto temp_product = remove_tensor(product, L"A");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product->append(1, ex<Tensor>(new_tensor));
      } else {
        new_product->append(1, term);
      }
    }
    new_product->scale(phase);
    new_result->append(new_product);
  }  // map_list

  detail::reset_idx_tags(new_result);

  return new_result;
}

ExprPtr symmetrize_expr(const ProductPtr& product) {
  auto result = std::make_shared<Sum>();

  // Drops canonical-order assumption; handles arbitrary sequence and variables.
  const auto& factors = product->factors();
  auto it = ranges::find_if(factors, [](const ExprPtr& factor) {
    return factor->is<Tensor>() && factor->as<Tensor>().label() == L"A";
  });
  if (it == ranges::end(factors)) return product;
  const auto& A_tensor = (*it)->as<Tensor>();
  SEQUANT_ASSERT(A_tensor.label() == L"A");

  auto A_is_nconserving = A_tensor.bra_rank() == A_tensor.ket_rank();

  if (A_is_nconserving && A_tensor.bra_rank() == 1)
    return remove_tensor(product, L"A");

  SEQUANT_ASSERT(A_tensor.rank() > 1);

  auto S = Tensor{};
  if (A_is_nconserving) {
    S = Tensor(L"S", A_tensor.bra(), A_tensor.ket(), A_tensor.aux(),
               Symmetry::Nonsymm);
  } else {  // A is N-nonconserving
    auto n = std::min(A_tensor.bra_rank(), A_tensor.ket_rank());
    container::svector<Index> bra_list(A_tensor.bra().begin(),
                                       A_tensor.bra().begin() + n);
    container::svector<Index> ket_list(A_tensor.ket().begin(),
                                       A_tensor.ket().begin() + n);
    S = Tensor(L"S", bra(std::move(bra_list)), ket(std::move(ket_list)),
               A_tensor.aux(), Symmetry::Nonsymm);
  }

  // Generate replacement maps from a list of Index type (could be a bra or a
  // ket)
  // Uses a permuted list of int to generate permutations
  // TODO factor out for reuse
  auto maps_from_list = [](const container::svector<Index>& list) {
    container::svector<int> int_list(list.size());
    std::iota(int_list.begin(), int_list.end(), 0);
    container::svector<container::map<Index, Index>> result;
    do {
      container::map<Index, Index> map;
      auto list_ptr = list.begin();
      for (auto&& i : int_list) {
        map.emplace(*list_ptr, list[i]);
        list_ptr++;
      }
      result.push_back(map);
    } while (std::next_permutation(int_list.begin(), int_list.end()));
    SEQUANT_ASSERT(result.size() ==
                   boost::numeric_cast<size_t>(factorial(list.size())));
    return result;
  };

  // Get phase relative to the canonical order
  // TODO factor out for reuse
  auto get_phase = [](const container::map<Index, Index>& map) {
    container::svector<Index> idx_list;
    for (const auto& [key, val] : map) idx_list.push_back(val);
    reset_ts_swap_counter<Index>();
    bubble_sort(std::begin(idx_list), std::end(idx_list));
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  container::svector<container::map<Index, Index>> maps;
  // CASE 1: n_bra = n_ket on all tensors
  if (A_is_nconserving) {
    maps = maps_from_list(A_tensor.bra());
  } else {
    SEQUANT_ASSERT(A_tensor.bra_rank() != A_tensor.ket_rank());
    maps = A_tensor.bra_rank() > A_tensor.ket_rank()
               ? maps_from_list(A_tensor.bra())
               : maps_from_list(A_tensor.ket());
  }
  SEQUANT_ASSERT(!maps.empty());
  for (auto&& map : maps) {
    Product new_product{};
    new_product.scale(product->scalar());
    new_product.append(get_phase(map), ex<Tensor>(S));
    auto temp_product = remove_tensor(product, L"A");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product.append(1, ex<Tensor>(new_tensor));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        new_product.append(1, term);
      } else {
        throw std::runtime_error("Invalid Expr type in symmetrize_expr: " +
                                 term->type_name());
      }
    }
    result->append(ex<Product>(new_product));
  }  // map
  return result;
}

ExprPtr symmetrize_expr(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  if (expr->is<Product>())
    return symmetrize_expr(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(symmetrize_expr(summand));
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in symmetrize_expr: " +
                             expr->type_name());
  }
}

ExprPtr expand_A_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  if (expr->is<Product>())
    return expand_A_op(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(expand_A_op(summand));
    }
    return result;
  }

  throw std::runtime_error("Invalid Expr type in expand_A_op: " +
                           expr->type_name());
}

container::svector<container::map<Index, Index>> P_maps(const Tensor& P) {
  SEQUANT_ASSERT(P.label() == L"P");

  // Return pair-wise replacements
  // P_ij -> {{i,j},{j,i}}
  // P_ijkl \equiv P_ij P_kl -> {{i,j},{j,i},{k,l},{l,k}}
  // P_ij^ab \equiv P_ij P^ab -> {{i,j},{j,i},{a,b},{b,a}}
  SEQUANT_ASSERT(P.bra_rank() % 2 == 0 && P.ket_rank() % 2 == 0);
  container::map<Index, Index> idx_rep;
  auto indices = P.const_braket_indices();
  for (auto it = indices.begin(); it != indices.end(); ranges::advance(it, 2)) {
    auto& idx1 = *it;
    auto it_next = it;
    ++it_next;
    SEQUANT_ASSERT(it_next != indices.end());
    auto& idx2 = *it_next;
    idx_rep.emplace(idx1, idx2);
    idx_rep.emplace(idx2, idx1);
  }

  SEQUANT_ASSERT(idx_rep.size() == (P.bra_net_rank() + P.ket_net_rank()));
  return container::svector<container::map<Index, Index>>{idx_rep};
}

ExprPtr expand_P_op(const ProductPtr& product) {
  bool has_P_operator = false;

  // Check P and build a replacement map
  // Assuming a product can have multiple P operators
  container::svector<container::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      const auto& P = term->as<Tensor>();
      if ((P.label() == L"P") && (P.bra_rank() > 1 || (P.ket_rank() > 1))) {
        has_P_operator = true;
        auto map = P_maps(P);
        map_list.insert(map_list.end(), map.begin(), map.end());
      } else if ((P.label() == L"P") &&
                 (P.bra_rank() == 1 && (P.ket_rank() == 1))) {
        return remove_tensor(product, L"P");
      }
    }
  }

  if (!has_P_operator) return product;

  auto result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    ProductPtr new_product = std::make_shared<Product>();
    new_product->scale(product->scalar());
    auto temp_product = remove_tensor(product, L"P");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_tensor.reset_tags();
        new_product->append(1, ex<Tensor>(new_tensor));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        new_product->append(1, term);
      } else {
        throw std::runtime_error("Invalid Expr type in expand_P_op: " +
                                 term->type_name());
      }
    }
    result->append(new_product);
  }  // map_list

  return result;
}

ExprPtr expand_P_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;
  else if (expr->is<Product>())
    return expand_P_op(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& summand : *expr) {
      result->append(expand_P_op(summand));
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in expand_P_op: " +
                             expr->type_name());
  }
}

container::svector<container::map<Index, Index>> S_replacement_maps(
    const Tensor& S) {
  SEQUANT_ASSERT(S.label() == L"S");
  SEQUANT_ASSERT(S.bra_rank() > 1);
  SEQUANT_ASSERT(S.bra().size() == S.ket().size());
  container::svector<int> int_list(S.bra().size());
  std::iota(std::begin(int_list), std::end(int_list), 0);

  container::svector<container::map<Index, Index>> maps;
  do {
    container::map<Index, Index> map;
    auto S_bra_ptr = S.bra().begin();
    auto S_ket_ptr = S.ket().begin();
    for (auto&& i : int_list) {
      map.emplace(*S_bra_ptr, S.bra()[i]);
      ++S_bra_ptr;
      map.emplace(*S_ket_ptr, S.ket()[i]);
      ++S_ket_ptr;
    }
    maps.push_back(map);
  } while (std::next_permutation(int_list.begin(), int_list.end()));

  return maps;
}

ExprPtr S_maps(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  auto result = std::make_shared<Sum>();

  // Check if S operator is present
  if (!has_tensor(expr, L"S")) return expr;

  detail::reset_idx_tags(expr);

  // Lambda for applying S on products
  auto expand_S_product = [](const ProductPtr& product) -> ExprPtr {
    // check if S is present
    if (!has_tensor(product, L"S")) return product;

    container::svector<container::map<Index, Index>> maps;
    // supports arbitrary sequence and variables
    for (auto&& factor : product->factors()) {
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"S") {
        maps = S_replacement_maps(factor->as<Tensor>());
        break;
      }
    }
    SEQUANT_ASSERT(!maps.empty());
    Sum sum{};
    for (auto&& map : maps) {
      ProductPtr new_product = std::make_shared<Product>();
      new_product->scale(product->scalar());
      auto temp_product = remove_tensor(product, L"S").as_shared_ptr<Product>();
      for (auto&& term : temp_product) {
        if (term->is<Tensor>()) {
          auto new_tensor = term->as<Tensor>();
          new_tensor.transform_indices(map);
          new_product->append(1, ex<Tensor>(new_tensor));
        } else {
          new_product->append(1, term);
        }
      }
      sum.append(new_product);
    }
    ExprPtr result = std::make_shared<Sum>(sum);
    return result;
  };

  if (expr->is<Product>()) {
    result->append(expand_S_product(expr.as_shared_ptr<Product>()));
  } else if (expr->is<Sum>()) {
    for (auto&& term : *expr) {
      if (term->is<Product>()) {
        result->append(expand_S_product(term.as_shared_ptr<Product>()));
      } else if (term->is<Tensor>() || term->is<Constant>() ||
                 expr->is<Variable>()) {
        result->append(term);
      }
    }
  }

  detail::reset_idx_tags(result);
  return result;
}

ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs) {
  if (!expr->is<Sum>()) return expr;
  if (ext_idxs.size() <= 2) return expr;  // always skip R1 and R2

  // hash filtering logic for R > 2
  container::map<std::size_t, container::vector<ExprPtr>> largest_coeff_terms;

  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;

    auto product = term.as_shared_ptr<Product>();
    auto scalar = product->scalar();

    sequant::TensorNetwork tn(*product);
    auto hash =
        tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
            .hash_value();

    auto it = largest_coeff_terms.find(hash);
    if (it == largest_coeff_terms.end()) {
      largest_coeff_terms[hash] = {term};
    } else {
      if (!it->second.empty()) {
        auto existing_scalar = it->second[0]->as<Product>().scalar();
        auto existing_abs = abs(existing_scalar);
        auto current_abs = abs(scalar);

        if (current_abs > existing_abs) {
          it->second.clear();
          it->second.push_back(term);
        } else if (current_abs == existing_abs) {
          it->second.push_back(term);
        }
      }
    }
  }

  Sum filtered;
  for (const auto& [_, terms] : largest_coeff_terms) {
    for (const auto& t : terms) {
      filtered.append(t);
    }
  }
  auto result = ex<Sum>(filtered);

  return result;
}

ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool full_expansion) {
  // Symmetrize and expression
  // Partially expand the antisymmetrizer and write it in terms of S operator.
  // See symmetrize_expr(expr) function for implementation details. We want an
  // expression with non-symmetric tensors, hence we are partially expanding the
  // antisymmetrizer (A) and fully expanding the anti-symmetric tensors to
  // non-symmetric.
  // full_expansion: it fully expands the antisymmetrizer directly (can be used
  // for v2 eqs, however it is not an optimized way).
  auto partially_or_fully_expand = [&full_expansion](const ExprPtr& expr) {
    auto temp = expr;
    if (has_tensor(temp, L"A")) {
      if (full_expansion) {
        temp = expand_A_op(temp);
      } else {
        temp = symmetrize_expr(temp);
      }
    }
    temp = expand_antisymm(temp);
    rapid_simplify(temp);
    return temp;
  };
  auto expr = partially_or_fully_expand(expression);

  // Index tags are cleaned prior to calling the fast canonicalizer
  detail::reset_idx_tags(expr);  // This call is REQUIRED
  expand(expr);                  // This call is REQUIRED
  simplify(expr);  // full simplify to combine terms before count_cycles

  // Lambda for spin-tracing a product term
  // For closed-shell case, a spin-traced result is a product term scaled by
  // 2^{n_cycles}, where n_cycles are counted by the lambda function described
  // above. For every product term, the bra indices on all tensors are merged
  // into a single list, so are the ket indices. External indices are
  // substituted with either one of the index (because the two vectors should be
  // permutations of each other to count cycles). All tensors must be nonsymm.
  auto trace_product = [&ext_index_groups](const Product& product) {
    // Remove S if present in a product
    Product temp_product{};
    temp_product.scale(product.scalar());
    if (product.factor(0).is<Tensor>() &&
        product.factor(0)->as<Tensor>().label() == L"S") {
      for (auto&& term : product.factors()) {
        if (term->is<Tensor>() && term->as<Tensor>().label() != L"S") {
          temp_product.append(1, term, Product::Flatten::No);
        } else if (!term->is<Tensor>()) {
          temp_product.append(1, term, Product::Flatten::No);
        }
      }
    } else {
      temp_product = product;
    }

    const bool collect_symmetrizer_indices = ext_index_groups.empty();

    auto get_ket_indices = [&](const Product& prod) {
      container::svector<Index> ket_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>() && (collect_symmetrizer_indices ||
                                (t->as<Tensor>().label() != L"A" &&
                                 t->as<Tensor>().label() != L"S"))) {
          const Tensor& tensor = t->as<Tensor>();
          ket_idx.insert(ket_idx.end(), tensor.ket().begin(),
                         tensor.ket().end());
        }
      }
      return ket_idx;
    };
    auto product_kets = get_ket_indices(temp_product);

    auto get_bra_indices = [&](const Product& prod) {
      container::svector<Index> bra_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>() && (collect_symmetrizer_indices ||
                                (t->as<Tensor>().label() != L"A" &&
                                 t->as<Tensor>().label() != L"S"))) {
          const Tensor& tensor = t->as<Tensor>();
          bra_idx.insert(bra_idx.end(), tensor.bra().begin(),
                         tensor.bra().end());
        }
      }
      return bra_idx;
    };
    auto product_bras = get_bra_indices(temp_product);

    auto substitute_ext_idx = [&product_bras, &product_kets](
                                  const container::svector<Index>& idx_pair) {
      SEQUANT_ASSERT(idx_pair.size() == 2);
      const auto& what = idx_pair[0];
      const auto& with = idx_pair[1];
      std::replace(product_bras.begin(), product_bras.end(), what, with);
      std::replace(product_kets.begin(), product_kets.end(), what, with);
    };

    // Substitute indices from external index list
    ranges::for_each(ext_index_groups, substitute_ext_idx);

    auto n_cycles = count_cycles(product_kets, product_bras);

    auto result = std::make_shared<Product>(product);
    result->scale(pow2(n_cycles));
    return result;
  };

  if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else if (expr->is<Tensor>())
    return trace_product(
        (ex<Constant>(1) * expr)->as<Product>());  // expand_all(expr);
  else if (expr->is<Product>())
    return trace_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      if (summand->is<Product>()) {
        result->append(trace_product(summand->as<Product>()));
      } else if (summand->is<Tensor>()) {
        result->append(
            trace_product((ex<Constant>(1) * summand)->as<Product>()));
      } else {
        SEQUANT_ASSERT(summand->is<Constant>() || summand->is<Variable>());
        result->append(summand);
      }
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in closed_shell_spintrace: " +
                             expr->type_name());
  }
}

container::svector<ResultExpr> closed_shell_spintrace(const ResultExpr& expr,
                                                      bool full_expansion) {
  using TraceFunction =
      ExprPtr (*)(const ExprPtr&,
                  const container::svector<container::svector<Index>>&, bool);

  return detail::wrap_trace<container::svector<ResultExpr>>(
      expr, static_cast<TraceFunction>(&closed_shell_spintrace),
      full_expansion);
}

ExprPtr closed_shell_CC_spintrace_v1(ExprPtr const& expr,
                                     ClosedShellCCSpintraceOptions options) {
  SEQUANT_ASSERT(options.method == BiorthogonalizationMethod::V1);
  using ranges::views::transform;

  auto const ext_idxs = external_indices(expr);
  auto st_expr = options.naive_spintrace
                     ? spintrace(expr, ext_idxs)
                     : closed_shell_spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Remove S operator
    for (auto& term : *st_expr) {
      if (term->is<Product>())
        term = remove_tensor(term.as_shared_ptr<Product>(), L"S");
    }

    // Biorthogonal transformation
    st_expr = biorthogonal_transform(st_expr, ext_idxs);

    auto bixs = ext_idxs | transform([](auto&& vec) { return vec[1]; });
    auto kixs = ext_idxs | transform([](auto&& vec) { return vec[0]; });
    st_expr =
        ex<Tensor>(Tensor{L"S", bra(std::move(bixs)), ket(std::move(kixs))}) *
        st_expr;
  }

  simplify(st_expr);

  return st_expr;
}

ExprPtr closed_shell_CC_spintrace_v2(ExprPtr const& expr,
                                     ClosedShellCCSpintraceOptions options) {
  SEQUANT_ASSERT(options.method == BiorthogonalizationMethod::V2);
  using ranges::views::transform;

  auto const ext_idxs = external_indices(expr);
  auto st_expr = options.naive_spintrace
                     ? spintrace(expr, ext_idxs)
                     : closed_shell_spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Remove S operator to apply biorthogonal transformation
    for (auto& term : *st_expr) {
      if (term->is<Product>())
        term = remove_tensor(term.as_shared_ptr<Product>(), L"S");
    }
    st_expr = biorthogonal_transform(st_expr, ext_idxs);

    // adding S in order to expand it and have all the raw equations
    auto bixs = ext_idxs | transform([](auto&& vec) { return vec[1]; });
    auto kixs = ext_idxs | transform([](auto&& vec) { return vec[0]; });
    if (bixs.size() > 1) {
      st_expr =
          ex<Tensor>(Tensor{L"S", bra(std::move(bixs)), ket(std::move(kixs))}) *
          st_expr;
    }
    simplify(st_expr);
    // expanding S after spintracing and biorthogonalization, to avoid dealing
    // with large number of terms
    st_expr = S_maps(st_expr);
    // canonicalizer must be called before hash-filter to combine terms
    canonicalize(st_expr);

    // apply hash filter method to get unique set of terms
    st_expr = WK_biorthogonalization_filter(st_expr, ext_idxs);
    // add S tensor again
    st_expr =
        ex<Tensor>(Tensor{L"S", bra(std::move(bixs)), ket(std::move(kixs))}) *
        st_expr;

    rational combined_factor;
    if (ext_idxs.size() <= 2) {
      combined_factor = rational(1, factorial(ext_idxs.size()));
    } else {
      auto fact_n = factorial(ext_idxs.size());
      combined_factor =
          rational(1, fact_n - 1);  // this is (1/fact_n) * (fact_n/(fact_n-1))
    }
    st_expr = ex<Constant>(combined_factor) * st_expr;
  }

  simplify(st_expr);
  // std::wcout << "final eqs after symm: "
  //            << sequant::to_latex_align(
  //                   sequant::ex<sequant::Sum>(
  //                       sequant::opt::reorder(st_expr->as<sequant::Sum>())),
  //                   0, 4)
  //            << std::endl;

  return st_expr;
}

ExprPtr closed_shell_CC_spintrace(ExprPtr const& expr,
                                  ClosedShellCCSpintraceOptions options) {
  switch (options.method) {
    case BiorthogonalizationMethod::V1:
      return closed_shell_CC_spintrace_v1(expr, options);
    case BiorthogonalizationMethod::V2:
      return closed_shell_CC_spintrace_v2(expr, options);
    default:
      SEQUANT_ASSERT(false && "unreachable code reached");
      abort();
  }
}

/// Collect all indices from an expression
container::set<Index, Index::LabelCompare> index_list(const ExprPtr& expr) {
  container::set<Index, Index::LabelCompare> grand_idxlist;
  if (expr->is<Tensor>()) {
    ranges::for_each(expr->as<Tensor>().const_indices(),
                     [&grand_idxlist](const Index& idx) {
                       idx.reset_tag();
                       grand_idxlist.insert(idx);
                     });
  }

  for (const ExprPtr& subExpr : expr) {
    grand_idxlist.merge(index_list(subExpr));
  }

  return grand_idxlist;
}

Tensor swap_spin(const Tensor& t) {
  auto is_any_spin = [](const Index& i) {
    return mbpt::to_spin(i.space().qns()) == mbpt::Spin::any;
  };

  // Return tensor if there are no spin labels
  if (std::all_of(t.const_braket_indices().begin(),
                  t.const_braket_indices().end(), is_any_spin)) {
    return t;
  }

  // Return new index where the spin-label is flipped
  auto spin_flipped_idx = [](const Index& idx) {
    SEQUANT_ASSERT(mbpt::to_spin(idx.space().qns()) != mbpt::Spin::any);
    return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha
               ? make_spinbeta(idx)
               : make_spinalpha(idx);
  };

  container::svector<Index> b(t.rank()), k(t.rank());

  for (std::size_t i = 0; i < t.rank(); ++i) {
    b.at(i) = spin_flipped_idx(t.bra().at(i));
    k.at(i) = spin_flipped_idx(t.ket().at(i));
  }

  return {t.label(),    bra(std::move(b)),   ket(std::move(k)),  t.aux(),
          t.symmetry(), t.braket_symmetry(), t.column_symmetry()};
}

ExprPtr swap_spin(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  auto swap_tensor = [](const Tensor& t) { return ex<Tensor>(swap_spin(t)); };

  auto swap_product = [&swap_tensor](const Product& p) {
    Product result{};
    result.scale(p.scalar());
    for (auto& t : p) {
      if (t->is<Tensor>()) {
        result.append(1, swap_tensor(t->as<Tensor>()), Product::Flatten::No);
      } else if (t->is<Constant>() || t->is<Variable>()) {
        result.append(1, t, Product::Flatten::No);
      } else {
        throw std::runtime_error(
            "Invalid Expr type in swap_spin::swap_product: " + t->type_name());
      }
    }
    return ex<Product>(result);
  };

  if (expr->is<Tensor>())
    return swap_tensor(expr->as<Tensor>());
  else if (expr->is<Product>())
    return swap_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    Sum result;
    for (auto& term : *expr) {
      result.append(swap_spin(term));
    }
    return ex<Sum>(result);
  } else {
    throw std::runtime_error("Invalid Expr type in swap_spin: " +
                             expr->type_name());
  }
}

ExprPtr merge_tensors(const Tensor& O1, const Tensor& O2) {
  SEQUANT_ASSERT(O1.label() == O2.label());
  SEQUANT_ASSERT(O1.symmetry() == O2.symmetry());
  auto b = ranges::views::concat(O1.bra(), O2.bra());
  auto k = ranges::views::concat(O1.ket(), O2.ket());
  auto a = ranges::views::concat(O1.aux(), O2.aux());
  return ex<Tensor>(Tensor(O1.label(), bra(b), ket(k), aux(a), O1.symmetry()));
}

std::vector<ExprPtr> open_shell_A_op(const Tensor& A) {
  SEQUANT_ASSERT(A.label() == L"A");
  SEQUANT_ASSERT(A.bra_rank() == A.ket_rank());
  auto rank = A.bra_rank();

  std::vector<ExprPtr> result(rank + 1);
  result.at(0) = ex<Constant>(1);
  result.at(rank) = ex<Constant>(1);

  for (std::size_t i = 1; i < rank; ++i) {
    auto spin_bra = A.bra();
    auto spin_ket = A.ket();
    std::transform(spin_bra.begin(), spin_bra.end() - i, spin_bra.begin(),
                   make_spinalpha);
    std::transform(spin_ket.begin(), spin_ket.end() - i, spin_ket.begin(),
                   make_spinalpha);
    std::transform(spin_bra.end() - i, spin_bra.end(), spin_bra.end() - i,
                   make_spinbeta);
    std::transform(spin_ket.end() - i, spin_ket.end(), spin_ket.end() - i,
                   make_spinbeta);
    ranges::for_each(spin_bra, [](const Index& i) { i.reset_tag(); });
    ranges::for_each(spin_ket, [](const Index& i) { i.reset_tag(); });
    result.at(i) = ex<Tensor>(
        Tensor(L"A", spin_bra, spin_ket, A.aux(), Symmetry::Antisymm));
    // std::wcout << to_latex(result.at(i)) << " ";
  }
  // std::wcout << "\n" << std::endl;
  return result;
}

std::vector<ExprPtr> open_shell_P_op_vector(const Tensor& A) {
  SEQUANT_ASSERT(A.label() == L"A");

  // N+1 spin-cases for corresponding residual
  std::vector<ExprPtr> result_vector(A.bra_rank() + 1);

  // List of indices
  const auto rank = A.bra_rank();
  container::svector<int> idx(rank);
  std::iota(idx.begin(), idx.end(), 0);

  // Anti-symmetrizer is preserved for all identical spin cases,
  // So return a constant
  result_vector.at(0) = ex<Constant>(1);     // all alpha
  result_vector.at(rank) = ex<Constant>(1);  // all beta

  // This loop generates all the remaining spin cases
  for (std::size_t i = 1; i < rank; ++i) {
    container::svector<int> alpha_spin(idx.begin(), idx.end() - i);
    container::svector<int> beta_spin(idx.end() - i, idx.end());

    container::svector<Tensor> P_bra_list, P_ket_list;
    for (auto& j : alpha_spin) {
      for (auto& k : beta_spin) {
        if (!alpha_spin.empty() && !beta_spin.empty()) {
          P_bra_list.emplace_back(Tensor(
              L"P", bra{A.bra().at(j), A.bra().at(k)}, ket{}, Symmetry::Symm));
          P_ket_list.emplace_back(Tensor(
              L"P", bra{}, ket{A.ket().at(j), A.ket().at(k)}, Symmetry::Symm));
        }
      }
    }

    // The P4 terms
    if (alpha_spin.size() > 1 && beta_spin.size() > 1) {
      for (std::size_t a = 0; a != alpha_spin.size() - 1; ++a) {
        auto i1 = alpha_spin[a];
        for (std::size_t b = a + 1; b != alpha_spin.size(); ++b) {
          auto i2 = alpha_spin[b];
          for (std::size_t c = 0; c != beta_spin.size() - 1; ++c) {
            auto i3 = beta_spin[c];
            for (std::size_t d = c + 1; d != beta_spin.size(); ++d) {
              auto i4 = beta_spin[d];
              P_bra_list.emplace_back(
                  Tensor(L"P",
                         bra{A.bra().at(i1), A.bra().at(i3), A.bra().at(i2),
                             A.bra().at(i4)},
                         ket{}, Symmetry::Symm));
              P_ket_list.emplace_back(
                  Tensor(L"P", bra{},
                         ket{A.ket().at(i1), A.ket().at(i3), A.ket().at(i2),
                             A.ket().at(i4)},
                         Symmetry::Symm));
            }
          }
        }
      }
    }

    Sum bra_permutations{};
    bra_permutations.append(ex<Constant>(1));
    Sum ket_permutations{};
    ket_permutations.append(ex<Constant>(1));

    for (auto& p : P_bra_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      bra_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    for (auto& p : P_ket_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      ket_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    ExprPtr spin_case_result =
        ex<Sum>(bra_permutations) * ex<Sum>(ket_permutations);
    expand(spin_case_result);

    // Merge P operators if it encounters alpha_spin product of operators
    for (auto& term : *spin_case_result) {
      if (term->is<Product>()) {
        auto P = term->as<Product>();
        if (P.factors().size() == 2) {
          auto P1 = P.factor(0)->as<Tensor>();
          auto P2 = P.factor(1)->as<Tensor>();
          term = merge_tensors(P1, P2);
        }
      }
    }
    result_vector.at(i) = spin_case_result;
  }
  return result_vector;
}

std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    std::optional<int> target_spin_case) {
  if (expr->is<Constant>() || expr->is<Variable>()) {
    return std::vector<ExprPtr>{expr};
  }

  // Grand index list contains both internal and external indices
  container::set<Index, Index::LabelCompare> grand_idxlist = index_list(expr);

  container::set<Index> ext_idxlist;
  for (auto&& idxgrp : ext_index_groups) {
    for (auto&& idx : idxgrp) {
      idx.reset_tag();
      ext_idxlist.insert(idx);
    }
  }

  container::set<Index> int_idxlist;
  for (auto&& gidx : grand_idxlist) {
    if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
      int_idxlist.insert(gidx);
    }
  }

  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> int_index_groups;
  for (auto&& i : int_idxlist) {
    int_index_groups.emplace_back(IndexGroup(1, i));
  }

  SEQUANT_ASSERT(grand_idxlist.size() ==
                 int_idxlist.size() + ext_idxlist.size());

  // make a spin-specific index, orientation is given by spin_bit: 0 =
  // spin-down/beta, 1 = spin-up/alpha
  auto make_spinspecific = [](const Index& idx, const long int& spin_bit) {
    return spin_bit == 0 ? make_spinalpha(idx) : make_spinbeta(idx);
  };

  // Generate index replacement maps
  auto spin_cases = [&make_spinspecific](
                        const container::svector<IndexGroup>& idx_group) {
    const auto ncases = pow2(idx_group.size());
    container::svector<container::map<Index, Index>> all_replacements(ncases);

    for (uint64_t i = 0; i != ncases; ++i) {
      container::map<Index, Index> idx_rep;
      for (size_t idxg = 0; idxg != idx_group.size(); ++idxg) {
        auto spin_bit = (i << (64 - idxg - 1)) >> 63;
        SEQUANT_ASSERT((spin_bit == 0) || (spin_bit == 1));
        for (auto& idx : idx_group[idxg]) {
          auto spin_idx = make_spinspecific(idx, spin_bit);
          idx_rep.emplace(idx, spin_idx);
        }
      }
      all_replacements[i] = idx_rep;
    }
    return all_replacements;
  };

  // External index replacement maps
  auto ext_spin_cases =
      [&make_spinspecific](const container::svector<IndexGroup>& idx_group) {
        container::svector<container::map<Index, Index>> all_replacements;

        // container::svector<int> spins(idx_group.size(), 0);
        for (std::size_t i = 0; i <= idx_group.size(); ++i) {
          container::svector<int> spins(idx_group.size(), 0);
          std::fill(spins.end() - i, spins.end(), 1);

          container::map<Index, Index> idx_rep;
          for (std::size_t j = 0; j != idx_group.size(); ++j) {
            for (auto& idx : idx_group[j]) {
              auto spin_idx = make_spinspecific(idx, spins[j]);
              idx_rep.emplace(idx, spin_idx);
            }
          }
          all_replacements.push_back(idx_rep);
        }
        return all_replacements;
      };

  // Internal and external index replacements are independent
  auto i_rep = spin_cases(int_index_groups);
  auto e_rep = ext_spin_cases(ext_index_groups);

  // For a single spin case, keep only the relevant spin case
  // PS: all alpha indexing start at 0
  if (target_spin_case) {
    auto external_replacement_map = e_rep.at(*target_spin_case);
    e_rep.clear();
    e_rep.push_back(external_replacement_map);
  }

  // Expand 'A' operator and 'antisymm' tensors
  auto expanded_expr = expand_A_op(expr);
  detail::reset_idx_tags(expanded_expr);
  expand(expanded_expr);
  simplify(expanded_expr);

  std::vector<ExprPtr> result{};

  // return true if a product is spin-symmetric
  auto spin_symm_product = [](const Product& product) {
    container::svector<Index> cBra, cKet;  // concat Bra and concat Ket
    for (auto& term : product) {
      if (term->is<Tensor>()) {
        auto tnsr = term->as<Tensor>();
        cBra.insert(cBra.end(), tnsr.bra().begin(), tnsr.bra().end());
        cKet.insert(cKet.end(), tnsr.ket().begin(), tnsr.ket().end());
      } else if (term->is<Product>() || term->is<Sum>()) {
        throw std::runtime_error(
            "Nested Product and Sum not supported in spin_symm_product");
      }
    }
    SEQUANT_ASSERT(cKet.size() == cBra.size());

    auto i_ket = cKet.begin();
    for (auto& b : cBra) {
      if (b.space().qns() != i_ket->space().qns()) return false;
      ++i_ket;
    }
    return true;
  };

  //
  // SPIN-TRACING algorithm begins here
  //

  // Loop over external index replacement maps
  for (auto& e : e_rep) {
    // Add spin labels to external indices
    auto spin_expr = append_spin(expanded_expr, e);
    detail::reset_idx_tags(spin_expr);
    Sum e_result{};

    // Loop over internal index replacement maps
    for (auto& i : i_rep) {
      // Add spin labels to internal indices, expand antisymmetric tensors
      auto spin_expr_i = append_spin(spin_expr, i);
      spin_expr_i = expand_antisymm(spin_expr_i, true);
      expand(spin_expr_i);
      detail::reset_idx_tags(spin_expr_i);
      Sum i_result{};

      if (spin_expr_i->is<Tensor>() || spin_expr_i->is<Constant>() ||
          spin_expr_i->is<Variable>()) {
        e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Product>()) {
        if (spin_symm_product(spin_expr_i->as<Product>()))
          e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Sum>()) {
        for (auto& pr : *spin_expr_i) {
          if (pr->is<Product>()) {
            if (spin_symm_product(pr->as<Product>())) i_result.append(pr);
          } else if (pr->is<Tensor>()) {
            if (ms_conserving_columns(pr->as<Tensor>())) i_result.append(pr);
          } else if (pr->is<Constant>() || pr->is<Variable>()) {
            i_result.append(pr);
          } else
            throw("Unknown ExprPtr type.");
        }
        e_result.append(std::make_shared<Sum>(i_result));
      }

    }  // loop over internal indices
    result.push_back(std::make_shared<Sum>(e_result));
  }  // loop over external indices

  if (target_spin_case) {
    SEQUANT_ASSERT(result.size() == 1 &&
                   "Spin-specific case must return one expression.");
  }

  // Canonicalize and simplify all expressions
  for (auto& expression : result) {
    detail::reset_idx_tags(expression);
    canonicalize(expression);
    rapid_simplify(expression);
  }
  return result;
}

std::vector<ExprPtr> open_shell_CC_spintrace(const ExprPtr& expr) {
  Tensor A = expr->at(0)->at(0)->as<Tensor>();
  SEQUANT_ASSERT(A.label() == L"A");
  size_t const i = A.rank();
  auto P_vec = open_shell_P_op_vector(A);
  auto A_vec = open_shell_A_op(A);
  SEQUANT_ASSERT(P_vec.size() == i + 1);
  std::vector<Sum> concat_terms(i + 1);
  [[maybe_unused]] size_t n_spin_orbital_term = 0;
  for (auto& product_term : *expr) {
    auto term = remove_tensor(product_term.as_shared_ptr<Product>(), L"A");
    std::vector<ExprPtr> os_st(i + 1);

    // Apply the P operators on the product term without the A,
    // Expand the P operators and spin-trace the expression
    // Then apply A operator, canonicalize and remove A operator
    for (std::size_t s = 0; s != os_st.size(); ++s) {
      os_st.at(s) = P_vec.at(s) * term;
      expand(os_st.at(s));
      os_st.at(s) = expand_P_op(os_st.at(s));
      os_st.at(s) =
          open_shell_spintrace(os_st.at(s), external_indices(A), s).at(0);
      if (i > 2) {
        os_st.at(s) = A_vec.at(s) * os_st.at(s);
        simplify(os_st.at(s));
        os_st.at(s) = remove_tensor(os_st.at(s), L"A");
      }
    }

    for (size_t j = 0; j != os_st.size(); ++j) {
      concat_terms.at(j).append(os_st.at(j));
    }
    ++n_spin_orbital_term;
  }

  // Combine spin-traced terms for the current residual
  std::vector<ExprPtr> expr_vec;
  for (auto& spin_case : concat_terms) {
    auto ptr = sequant::ex<Sum>(spin_case);
    expr_vec.push_back(ptr);
  }

  return expr_vec;
}

ExprPtr spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool spinfree_index_spaces) {
  // Escape immediately if expression is a constant
  if (expression->is<Constant>() || expression->is<Variable>()) {
    return expression;
  }

#ifdef SEQUANT_ASSERT_ENABLED
  // Verify that the number of external indices matches the number of indices in
  // ext_index_groups, UNLESS user overrode external definitions in default
  // context
  const auto& copts = get_default_context().canonicalization_options();
  if (!copts.has_value() || !copts->named_indices) {
    auto count_indices = [](const auto& range) {
      auto sizes = range | ranges::views::transform(
                               [](const auto& list) { return list.size(); });
      return std::accumulate(sizes.begin(), sizes.end(), 0);
    };
    auto determined_externals =
        external_indices<container::svector<container::svector<Index>>>(
            expression);

    SEQUANT_ASSERT(count_indices(ext_index_groups) ==
                   count_indices(determined_externals));
  }
#endif

  // This function must be used for tensors with spin-specific indices only. If
  // the spin-symmetry is conserved: the tensor is expanded; else: zero is
  // returned.
  auto spintrace_tensor = [](const Tensor& tensor) {
    return can_expand(tensor) ? expand_antisymm(tensor) : ex<Constant>(0);
  };

  // This function is used to spin-trace a product terms with spin-specific
  // indices. It checks if all tensors can be expanded and spintraces individual
  // tensors by call to the spin_trace_tensor lambda.
  auto spintrace_product =
      [&spintrace_tensor](const ProductPtr& product) -> ExprPtr {
    ProductPtr spin_product = std::make_shared<Product>();

    // Check if all tensors in this product can be expanded
    // If NOT all tensors can be expanded, return zero
    for (const ExprPtr& expr : product->factors()) {
      if (expr.is<Tensor>()) {
        if (!can_expand(expr.as<Tensor>())) {
          return ex<Constant>(0);
        }
      } else if (expr.is<Sum>() || expr.is<Product>()) {
        throw std::runtime_error(
            "Nested sums/products not supported in spin_trace_product");
      }
    }

    spin_product->scale(product->scalar());
    for (const ExprPtr& expr : product->factors()) {
      if (expr.is<Tensor>()) {
        spin_product->append(1, spintrace_tensor(expr.as<Tensor>()));
      } else if (expr.is<Variable>() || expr.is<Constant>()) {
        spin_product->append(1, expr.clone());
      } else {
        // Would need some sort of recursion but it is not clear how that would
        // interact with other code in here yet so prefer to error instead.
        throw std::runtime_error(
            "spin_trace_product: Nested products or sums inside of a Product "
            "not supported (yet)");
      }
    }

    ExprPtr result = spin_product;
    expand(result);
    rapid_simplify(result);
    return result;
  };

  // Most important lambda of this function
  auto trace_product = [&ext_index_groups, &spintrace_tensor,
                        &spintrace_product,
                        spinfree_index_spaces](const ProductPtr& product) {
    ExprPtr expr = product->clone();
    // List of all indices in the expression
    container::set<Index, Index::LabelCompare> grand_idxlist = index_list(expr);

    // List of external indices, i.e. indices that are not summed over Einstein
    // style (indices that are not repeated in an expression)
    container::set<Index> ext_idxlist;
    for (auto&& idxgrp : ext_index_groups) {
      for (auto&& idx : idxgrp) {
        idx.reset_tag();
        ext_idxlist.insert(idx);
      }
    }

    // List of internal indices, i.e. indices that are contracted over
    container::set<Index> int_idxlist;
    for (auto&& gidx : grand_idxlist) {
      if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
        int_idxlist.insert(gidx);
      }
    }

    // EFV: generate the grand list of index groups by concatenating list of
    // external index groups with the groups of internal indices (each
    // internal index = 1 group)
    // TODO some internal indices can be a priori placed in the same group, if
    // they refer to the same particle of a spin-free non-antisymmetrized Tensor
    //      so visit all Tensors in the expression and locate such groups of
    //      internal indices before placing the rest into separate groups
    using IndexGroup = container::svector<Index>;
    container::svector<IndexGroup> index_groups;
    for (auto&& i : int_idxlist) index_groups.emplace_back(IndexGroup(1, i));
    index_groups.insert(index_groups.end(), ext_index_groups.begin(),
                        ext_index_groups.end());

    // EFV: for each spincase (loop over integer from 0 to 2^n-1, n=#of index
    // groups)
    SEQUANT_ASSERT(index_groups.size() <= 64);
    const uint64_t nspincases = pow2(index_groups.size());

    auto result = std::make_shared<Sum>();
    for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
         ++spincase_bitstr) {
      // EFV:  assign spin to each index group => make a replacement list
      container::map<Index, Index> index_replacements;

      uint64_t index_group_count = 0;
      for (auto&& index_group : index_groups) {
        auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
        SEQUANT_ASSERT(spin_bit == 0 || spin_bit == 1);

        for (auto&& index : index_group) {
          index_replacements.emplace(index, spin_bit == 0
                                                ? make_spinalpha(index)
                                                : make_spinbeta(index));
        }
        ++index_group_count;
      }

      // Append spin labels to indices in the expression
      auto spin_expr = append_spin(expr, index_replacements);
      rapid_simplify(spin_expr);  // This call is required for Tensor case

      // NB: There are temporaries in the following code to enable
      // printing intermediate expressions.
      if (spin_expr->is<Tensor>()) {
        auto st_expr = spintrace_tensor(spin_expr->as<Tensor>());
        result->append(spinfree_index_spaces ? remove_spin(st_expr) : st_expr);
      } else if (spin_expr->is<Product>()) {
        auto st_expr = spintrace_product(spin_expr.as_shared_ptr<Product>());
        if (!st_expr->is<Constant>() || st_expr->as<Constant>().value() != 0) {
          result->append(spinfree_index_spaces ? remove_spin(st_expr)
                                               : st_expr);
        }
      } else if (spin_expr->is<Sum>()) {
        for (auto&& summand : *spin_expr) {
          Sum st_expr{};
          if (summand->is<Tensor>())
            st_expr.append(spintrace_tensor(summand->as<Tensor>()));
          else if (summand->is<Product>())
            st_expr.append(spintrace_product(summand.as_shared_ptr<Product>()));
          else {
            st_expr.append(summand);
          }
          auto st_expr_ptr = ex<Sum>(st_expr);
          result->append(spinfree_index_spaces ? remove_spin(st_expr_ptr)
                                               : st_expr_ptr);
        }
      } else {
        result->append(expr);
      }
    }  // Permutation FOR loop
    return result;
  };

  // Expand antisymmetrizer operator (A) if present in the expression
  ExprPtr expr = expression;
  if (has_tensor(expr, L"A")) expr = expand_A_op(expr);

  if (expr->is<Tensor>()) expr = ex<Constant>(1) * expr;

  ExprPtr result;
  if (expr->is<Product>()) {
    result = trace_product(expr.as_shared_ptr<Product>());
  } else if (expr->is<Sum>()) {
    auto result_sum = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      if (term->is<Product>())
        result_sum->append(trace_product(term.as_shared_ptr<Product>()));
      else if (term->is<Tensor>()) {
        auto term_as_product = ex<Constant>(1) * term;
        result_sum->append(
            trace_product(term_as_product.as_shared_ptr<Product>()));
      } else
        result_sum->append(term);
      result = result_sum;
    }
    return result;
  } else {
    throw std::runtime_error("Invalid Expr type in spintrace: " +
                             expr->type_name());
  }

  detail::reset_idx_tags(result);
  return result;
}  // ExprPtr spintrace

container::svector<ResultExpr> spintrace(const ResultExpr& expr,
                                         bool spinfree_index_spaces) {
  using TraceFunction =
      ExprPtr (*)(const ExprPtr&,
                  const container::svector<container::svector<Index>>&, bool);

  return detail::wrap_trace<container::svector<ResultExpr>>(
      expr, static_cast<TraceFunction>(&spintrace), spinfree_index_spaces);
}

ExprPtr factorize_S(const ExprPtr& expression,
                    std::initializer_list<IndexList> ext_index_groups,
                    const bool fast_method) {
  // Canonicalize the expression
  ExprPtr expr = expression;
  // canonicalize(expr);

  // If expression has S operator, do nothing and exit
  if (has_tensor(expr, L"S")) return expr;

  // If S operator is absent: generate from ext_index_groups
  Tensor S{};
  {
    container::svector<Index> bra_list, ket_list;

    // Fill bras and kets
    ranges::for_each(ext_index_groups, [&](const IndexList& idx_pair) {
      auto it = idx_pair.begin();
      bra_list.push_back(*it);
      it++;
      ket_list.push_back(*it);
    });
    SEQUANT_ASSERT(bra_list.size() == ket_list.size());
    S = Tensor(L"S", bra(std::move(bra_list)), ket(std::move(ket_list)),
               Symmetry::Nonsymm);
  }

  // For any order CC residual equation:
  // Generate a list of permutation indices
  // Erase the canonical entry
  auto replacement_maps = S_replacement_maps(S);
  replacement_maps.erase(replacement_maps.begin());

  // Lambda function for index replacement in tensor
  auto transform_tensor =
      [](const Tensor& tensor,
         const container::map<Index, Index>& replacement_map) {
        auto result = std::make_shared<Tensor>(tensor);
        result->transform_indices(replacement_map);
        result->reset_tags();
        return result;
      };

  Sum result_sum{};
  ///////////////////////////////////////////////
  ///            Fast approach                ///
  ///////////////////////////////////////////////
  // This method hashes terms for faster run times

  if (fast_method) {
    // summands_hash_list sorted container of hash values of canonicalized
    // summands summands_hash_map unsorted map of (hash_val, summand) pairs
    // container::set<size_t> summands_hash_list;
    container::svector<size_t> summands_hash_list;
    std::unordered_map<size_t, ExprPtr> summands_hash_map;
    for (auto it = expr->begin(); it != expr->end(); ++it) {
      (*it)->canonicalize();
      auto hash = (*it)->hash_value();
      summands_hash_list.push_back(hash);
      summands_hash_map.emplace(hash, *it);
    }
    SEQUANT_ASSERT(summands_hash_list.size() == expr->size());
    SEQUANT_ASSERT(summands_hash_map.size() == expr->size());

    // Symmetrize every summand, assign its hash value to hash1
    // Check if hash1 exist in summands_hash_list
    // if(hash1 present in summands_hash_list) remove hash0, hash1
    // else continue
    [[maybe_unused]] int n_symm_terms = 0;
    auto symm_factor = factorial(S.bra_rank());
    for (auto it = expr->begin(); it != expr->end(); ++it) {
      // Exclude summand with value zero
      while ((*it)->hash_value() == ex<Constant>(0)->hash_value()) {
        ++it;
        if (it == expr->end()) break;
      }
      if (it == expr->end()) break;

      // Remove current hash_value from list and clone summand
      auto hash0 = (*it)->hash_value();
      summands_hash_list.erase(std::find(summands_hash_list.begin(),
                                         summands_hash_list.end(), hash0));
      auto new_product = (*it)->clone();
      new_product =
          ex<Constant>(rational{1, symm_factor}) * ex<Tensor>(S) * new_product;

      // CONTAINER OF HASH VALUES AND SYMMETRIZED TERMS
      // FOR GENERALIZED EXPRESSION WITH ARBITRARY S OPERATOR
      // Loop over replacement maps The entire code from here
      container::vector<size_t> hash1_list;
      for (auto&& replacement_map : replacement_maps) {
        size_t hash1 = 0;
        if ((*it)->is<Product>()) {
          // Clone *it, apply symmetrizer, store hash1 value
          auto product = (*it)->as<Product>();
          Product S_product{};
          S_product.scale(product.scalar());

          // Transform indices by action of S operator
          for (auto&& t : product) {
            if (t->is<Tensor>()) {
              S_product.append(
                  1, transform_tensor(t->as<Tensor>(), replacement_map),
                  Product::Flatten::No);
            } else if (t->is<Constant>() || t->is<Variable>()) {
              S_product.append(1, t, Product::Flatten::No);
            } else {
              throw std::runtime_error("Invalid Expr in factorize_S: " +
                                       t->type_name());
            }
          }
          auto new_product_expr = ex<Product>(S_product);
          new_product_expr->canonicalize();
          hash1 = new_product_expr->hash_value();

        } else if ((*it)->is<Tensor>()) {
          // Clone *it, apply symmetrizer, store hash value
          auto tensor = (*it)->as<Tensor>();

          // Transform indices by action of S operator
          auto new_tensor = transform_tensor(tensor, replacement_map);

          // Canonicalize the new tensor before computing hash value
          new_tensor->canonicalize();
          hash1 = new_tensor->hash_value();
        } else if ((*it)->is<Constant>() || (*it)->is<Variable>()) {
          hash1 = (*it)->hash_value();
        } else {
          throw std::runtime_error("Invalid Expr in factorize_S: " +
                                   (*it)->type_name());
        }
        hash1_list.push_back(hash1);
      }

      // bool symmetrizable = false;
      // auto hash1_found = [&](size_t h){ return summands_hash_list.find(h) !=
      // summands_hash_list.end();};
      auto hash1_found = [&summands_hash_list](size_t h) {
        return std::find(summands_hash_list.begin(), summands_hash_list.end(),
                         h) != summands_hash_list.end();
      };
      std::size_t n_hash_found = ranges::count_if(hash1_list, hash1_found);

      if (n_hash_found == hash1_list.size()) {
        // Prepend S operator
        // new_product = ex<Tensor>(S) * new_product;
        new_product = ex<Constant>(symm_factor) * new_product;
        ++n_symm_terms;

        // remove values from hash1_list from summands_hash_list
        ranges::for_each(hash1_list, [&](const size_t hash1) {
          // summands_hash_list.erase(hash1);
          summands_hash_list.erase(std::find(summands_hash_list.begin(),
                                             summands_hash_list.end(), hash1));

          auto term = summands_hash_map.find(hash1)->second;

          for (auto&& summand : *expr) {
            if (summand->hash_value() == hash1) summand = ex<Constant>(0);
          }
        });
      }
      result_sum.append(new_product);
    }
  } else {
    ///////////////////////////////////////////////
    ///            Lazy approach                ///
    ///////////////////////////////////////////////
    // This approach is slower because the hash values are computed on the fly.
    // Subsequently, this algorithm applies 'S' operator n^2 times

    [[maybe_unused]] int n_symm_terms = 0;

    // If a term was symmetrized, put the index in a list
    container::set<int> i_list;

    // The quick_method is a "faster" lazy approach that
    // symmetrizes *it instead of symmetrizing *find_it in the inside loop
    //    const auto tstart = std::chrono::high_resolution_clock::now();

    // Controls which term to symmetrize
    // true -> symmetrize the summand
    // false -> symmetrize lookup term
    // bool quick_method = true;

    // Loop over terms of expression (OUTER LOOP)
    for (auto it = expr->begin(); it != expr->end(); ++it) {
      // If *it is symmetrized, go to next
      while (std::find(i_list.begin(), i_list.end(),
                       std::distance(expr->begin(), it)) != i_list.end())
        ++it;
      // Clone the summand
      auto new_product = (*it)->clone();

      // hash value of summand
      container::vector<size_t> hash0_list;
      for (auto&& replacement_map : replacement_maps) {
        size_t hash0 = 0;
        if ((*it)->is<Product>()) {
          // Clone *it, apply symmetrizer, store hash value
          auto product = (*it)->as<Product>();
          Product S_product{};
          S_product.scale(product.scalar());

          // Transform indices by action of S operator
          for (auto&& t : product) {
            if (t->is<Tensor>()) {
              S_product.append(
                  1, transform_tensor(t->as<Tensor>(), replacement_map),
                  Product::Flatten::No);
            } else if (t->is<Constant>() || t->is<Variable>()) {
              S_product.append(1, t, Product::Flatten::No);
            } else {
              throw std::runtime_error("Invalid Expr type in factorize_S: " +
                                       t->type_name());
            }
          }
          auto new_product_expr = ex<Product>(S_product);
          new_product_expr->canonicalize();
          hash0 = new_product_expr->hash_value();
        } else if ((*it)->is<Tensor>()) {
          // Clone *it, apply symmetrizer, store hash value
          auto tensor = (*it)->as<Tensor>();

          // Transform indices by action of S operator
          auto new_tensor = transform_tensor(tensor, replacement_map);

          // Canonicalize the new tensor before computing hash value
          new_tensor->canonicalize();
          hash0 = new_tensor->hash_value();
        } else if ((*it)->is<Constant>() || (*it)->is<Variable>()) {
          hash0 = (*it)->hash_value();
        } else {
          throw std::runtime_error("Invalid Expr type in factorize_S: " +
                                   (*it)->type_name());
        }
        hash0_list.push_back(hash0);
      }

      for (auto&& hash0 : hash0_list) {
        std::size_t n_matches = 0;
        container::svector<size_t> idx_vec;
        for (auto find_it = it + 1; find_it != expr->end(); ++find_it) {
          auto idx = std::distance(expr->begin(), find_it);
          (*find_it)->canonicalize();

          if ((*find_it)->hash_value() == hash0) {
            ++n_matches;
            idx_vec.push_back(idx);
          }
        }
        if (n_matches == hash0_list.size()) {
          ++n_symm_terms;
          new_product = ex<Tensor>(S) * new_product;
          i_list.insert(idx_vec.begin(), idx_vec.end());
        }
      }

      // append product to running sum
      result_sum.append(new_product);
    }
    //    const auto tstop = std::chrono::high_resolution_clock::now();
    //    const auto time_elapsed =
    //        std::chrono::duration_cast<std::chrono::microseconds>(tstop -
    //        tstart);
    // std::cout << "Fast method: " << std::boolalpha << fast_method << "\n";
    // std::cout << "N symm terms found: " << n_symm_terms << "\n";
    // std::cout << "Time: " << time_elapsed.count() << " μs.\n";
  }

  ExprPtr result = std::make_shared<Sum>(result_sum);
  simplify(result);
  return result;
}

}  // namespace sequant::mbpt
