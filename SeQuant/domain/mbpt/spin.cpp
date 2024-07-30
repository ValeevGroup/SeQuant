#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr_algorithm.hpp>
#include <SeQuant/core/expr_operator.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/swap.hpp>

#ifdef SEQUANT_HAS_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/detail/variant.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/interface.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <algorithm>
#include <cassert>
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

namespace sequant {

namespace detail {

Index make_index_with_spincase(const Index& idx, mbpt::Spin s) {
  // sanity check: make sure have only one spin label
  assert(!(idx.label().find(L'↑') != std::wstring::npos &&
           idx.label().find(L'↓') != std::wstring::npos));

  // to preserve rest of bits first unset spin bit, then set them to the desired
  // state
  auto qns = mbpt::spinannotation_remove(idx.space().qns()).unIon(s);
  IndexSpace space{mbpt::spinannotation_replacе(idx.space().base_key(), s),
                   idx.space().type(), qns};
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_spincase(pidx, s);
  return Index{mbpt::spinannotation_replacе(idx.label(), s), space,
               protoindices};
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

ExprPtr transform_expr(const ExprPtr& expr,
                       const container::map<Index, Index>& index_replacements,
                       Constant::scalar_type scaling_factor) {
  if (expr->is<Constant>()) {
    return ex<Constant>(scaling_factor) * expr;
  }

  auto transform_tensor = [&index_replacements](const Tensor& tensor) {
    auto result = std::make_shared<Tensor>(tensor);
    result->transform_indices(index_replacements);
    result->reset_tags();
    return result;
  };

  auto transform_product = [&transform_tensor,
                            &scaling_factor](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        auto tensor = term->as<Tensor>();
        result->append(1, transform_tensor(tensor));
      }
    }
    result->scale(scaling_factor);
    return result;
  };

  if (expr->is<Tensor>()) {
    auto result =
        ex<Constant>(scaling_factor) * transform_tensor(expr->as<Tensor>());
    return result;
  } else if (expr->is<Product>()) {
    auto result = transform_product(expr->as<Product>());
    return result;
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& term : *expr) {
      result->append(transform_expr(term, index_replacements, scaling_factor));
    }
    return result;
  } else
    return nullptr;
}

ExprPtr swap_bra_ket(const ExprPtr& expr) {
  if (expr->is<Constant>()) return expr;

  // Lambda for tensor
  auto tensor_swap = [](const Tensor& tensor) {
    auto result = Tensor(tensor.label(), bra(tensor.ket().value()),
                         ket(tensor.bra().value()), tensor.symmetry(),
                         tensor.braket_symmetry(), tensor.particle_symmetry());
    return ex<Tensor>(result);
  };

  // Lambda for product
  auto product_swap = [&tensor_swap](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>())
        result->append(1, tensor_swap(term->as<Tensor>()),
                       Product::Flatten::No);
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
  } else
    return nullptr;
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
      if (term->is<Tensor>())
        spin_product->append(1, add_spin_to_tensor(term->as<Tensor>()));
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
  } else
    return nullptr;
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
    Tensor result(tensor.label(), bra(std::move(b)), ket(std::move(k)),
                  tensor.aux(), tensor.symmetry(), tensor.braket_symmetry());
    return std::make_shared<Tensor>(std::move(result));
  };

  auto remove_spin_from_product =
      [&remove_spin_from_tensor](const Product& product) {
        auto result = std::make_shared<Product>();
        result->scale(product.scalar());
        for (auto&& term : product) {
          if (term->is<Tensor>()) {
            result->append(1, remove_spin_from_tensor(term->as<Tensor>()));
          } else
            abort();
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
  } else
    return nullptr;
}

bool spin_symm_tensor(const Tensor& tensor) {
  assert(tensor.bra_rank() == tensor.ket_rank());
  for (std::size_t i = 0; i != tensor.rank(); ++i) {
    if (tensor.bra()[i].space().qns() != tensor.ket()[i].space().qns())
      return false;
  }
  return true;
}

bool same_spin_tensor(const Tensor& tensor) {
  auto braket = tensor.braket();
  auto spin_element = braket[0].space().qns();
  return std::all_of(braket.begin(), braket.end(),
                     [&spin_element](const auto& idx) {
                       return idx.space().qns() == spin_element;
                     });
}

bool can_expand(const Tensor& tensor) {
  assert(tensor.bra_rank() == tensor.ket_rank() &&
         "can_expand(Tensor) failed.");
  if (tensor.bra_rank() != tensor.ket_rank()) return false;

  // indices must have specific spin
  [[maybe_unused]] auto all_have_spin = std::all_of(
      tensor.const_braket().begin(), tensor.const_braket().end(),
      [](const auto& idx) {
        auto idx_spin = mbpt::to_spin(idx.space().qns());
        return idx_spin == mbpt::Spin::alpha || idx_spin == mbpt::Spin::beta;
      });
  assert(std::all_of(tensor.const_braket().begin(), tensor.const_braket().end(),
                     [](const auto& idx) {
                       auto idx_spin = mbpt::to_spin(idx.space().qns());
                       return idx_spin == mbpt::Spin::alpha ||
                              idx_spin == mbpt::Spin::beta;
                     }));

  // count alpha indices in bra
  auto is_alpha = [](const Index& idx) {
    return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha;
  };

  // count alpha indices in bra
  auto a_bra =
      std::count_if(tensor.bra().begin(), tensor.bra().end(), is_alpha);

  // count alpha indices in ket
  auto a_ket =
      std::count_if(tensor.ket().begin(), tensor.ket().end(), is_alpha);

  return a_bra == a_ket;
}

ExprPtr expand_antisymm(const Tensor& tensor, bool skip_spinsymm) {
  assert(tensor.bra_rank() == tensor.ket_rank());
  // Return non-symmetric tensor if rank is 1
  if (tensor.bra_rank() == 1) {
    Tensor new_tensor(tensor.label(), tensor.bra(), tensor.ket(), tensor.aux(),
                      Symmetry::nonsymm, tensor.braket_symmetry(),
                      tensor.particle_symmetry());
    return std::make_shared<Tensor>(new_tensor);
  }

  // If all indices have the same spin label,
  // return the antisymm tensor
  if (skip_spinsymm && same_spin_tensor(tensor)) {
    return std::make_shared<Tensor>(tensor);
  }

  assert(tensor.bra_rank() > 1 && tensor.ket_rank() > 1);

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
  if (tensor.symmetry() == Symmetry::antisymm) {
    const auto prefactor = get_phase(tensor);
    container::set<Index> bra_list(tensor.bra().begin(), tensor.bra().end());
    container::set<Index> ket_list(tensor.ket().begin(), tensor.ket().end());
    auto expr_sum = std::make_shared<Sum>();
    do {
      // N.B. must copy
      auto new_tensor = Tensor(tensor.label(), bra(bra_list), ket(ket_list),
                               tensor.aux(), Symmetry::nonsymm);

      if (spin_symm_tensor(new_tensor)) {
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
  if (expr->is<Constant>())
    return expr;
  else if (expr->is<Tensor>())
    return expand_antisymm(expr->as<Tensor>(), skip_spinsymm);

  // Product lambda
  auto expand_product = [&skip_spinsymm](const Product& expr) {
    Product temp{};
    temp.scale(expr.scalar());
    for (auto&& term : expr) {
      if (term->is<Tensor>())
        temp.append(1, expand_antisymm(term->as<Tensor>(), skip_spinsymm),
                    Product::Flatten::No);
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
  } else
    return nullptr;
}

bool has_tensor(const ExprPtr& expr, std::wstring label) {
  if (expr->is<Constant>()) return false;

  auto check_product = [&label](const Product& p) {
    return ranges::any_of(p.factors(), [&label](const auto& t) {
      return t->template is<Tensor>() && (t->template as<Tensor>()).label() == label;
    });
  };

  if (expr->is<Tensor>()) {
    return expr->as<Tensor>().label() == label;
  } else if (expr->is<Product>()) {
    return check_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    return ranges::any_of(
        *expr, [&label](const auto& term) { return has_tensor(term, label); });
  } else
    return false;
}

container::svector<container::map<Index, Index>> A_maps(const Tensor& A) {
  assert(A.label() == L"A");

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

ExprPtr remove_tensor(const Product& product, std::wstring label) {
  // filter out tensors with specified label
  auto new_product = std::make_shared<Product>();
  new_product->scale(product.scalar());
  for (auto&& term : product) {
    if (term->is<Tensor>()) {
      auto tensor = term->as<Tensor>();
      if (tensor.label() != label) new_product->append(1, ex<Tensor>(tensor));
    } else
      new_product->append(1, term);
  }
  return new_product;
}

ExprPtr remove_tensor(const ExprPtr& expr, std::wstring label) {
  if (expr->is<Sum>()) {
    Sum result{};
    for (auto& term : *expr) {
      result.append(remove_tensor(term, label));
    }
    return ex<Sum>(result);
  } else if (expr->is<Product>())
    return remove_tensor(expr->as<Product>(), label);
  else if (expr->is<Tensor>())
    return expr->as<Tensor>().label() == label ? ex<Constant>(1) : expr;
  else if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else
    return nullptr;
}

ExprPtr expand_A_op(const Product& product) {
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

  if (!has_A_operator) return std::make_shared<Product>(product);

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

    Product new_product{};
    new_product.scale(product.scalar());
    auto temp_product = remove_tensor(product, L"A");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product.append(1, ex<Tensor>(new_tensor));
      }
    }
    new_product.scale(phase);
    new_result->append(ex<Product>(new_product));
  }  // map_list

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>()) reset_tags(expr.as<Tensor>());
  };
  new_result->visit(reset_idx_tags, true);

  return new_result;
}

ExprPtr symmetrize_expr(const Product& product) {
  auto result = std::make_shared<Sum>();

  // Assumes canonical sequence of tensors in the product
  if (product.factor(0)->as<Tensor>().label() != L"A")
    return std::make_shared<Product>(product);

  // CHECK: A is present and >1 particle
  // GENERATE S tensor
  auto A_tensor = product.factor(0)->as<Tensor>();
  assert(A_tensor.label() == L"A");

  auto A_is_nconserving = A_tensor.bra_rank() == A_tensor.ket_rank();

  if (A_is_nconserving && A_tensor.bra_rank() == 1)
    return remove_tensor(product, L"A");

  assert(A_tensor.rank() > 1);

  auto S = Tensor{};
  if (A_is_nconserving) {
    S = Tensor(L"S", A_tensor.bra(), A_tensor.ket(), A_tensor.aux(),
               Symmetry::nonsymm);
  } else {  // A is N-nonconserving
    auto n = std::min(A_tensor.bra_rank(), A_tensor.ket_rank());
    container::svector<Index> bra_list(A_tensor.bra().begin(),
                                       A_tensor.bra().begin() + n);
    container::svector<Index> ket_list(A_tensor.ket().begin(),
                                       A_tensor.ket().begin() + n);
    S = Tensor(L"S", bra(std::move(bra_list)), ket(std::move(ket_list)),
               A_tensor.aux(), Symmetry::nonsymm);
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
    assert(result.size() ==
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
    assert(A_tensor.bra_rank() != A_tensor.ket_rank());
    maps = A_tensor.bra_rank() > A_tensor.ket_rank()
               ? maps_from_list(A_tensor.bra())
               : maps_from_list(A_tensor.ket());
  }
  assert(!maps.empty());
  for (auto&& map : maps) {
    Product new_product{};
    new_product.scale(product.scalar());
    new_product.append(get_phase(map), ex<Tensor>(S));
    auto temp_product = remove_tensor(product, L"A");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product.append(1, ex<Tensor>(new_tensor));
      }
    }
    result->append(ex<Product>(new_product));
  }  // map
  return result;
}

ExprPtr symmetrize_expr(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  if (expr->is<Product>())
    return symmetrize_expr(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(symmetrize_expr(summand));
    }
    return result;
  } else
    return nullptr;
}

ExprPtr expand_A_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  if (expr->is<Product>())
    return expand_A_op(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(expand_A_op(summand));
    }
    return result;
  } else
    return nullptr;
}

container::svector<container::map<Index, Index>> P_maps(const Tensor& P) {
  assert(P.label() == L"P");

  // Return pair-wise replacements
  // P_ij -> {{i,j},{j,i}}
  // P_ijkl \equiv P_ij P_kl -> {{i,j},{j,i},{k,l},{l,k}}
  // P_ij^ab \equiv P_ij P^ab -> {{i,j},{j,i},{a,b},{b,a}}
  assert(P.bra_rank() % 2 == 0 && P.ket_rank() % 2 == 0);
  container::map<Index, Index> idx_rep;
  for (std::size_t i = 0; i != P.const_braket().size(); i += 2) {
    idx_rep.emplace(P.const_braket().at(i), P.const_braket().at(i + 1));
    idx_rep.emplace(P.const_braket().at(i + 1), P.const_braket().at(i));
  }

  assert(idx_rep.size() == (P.bra_rank() + P.ket_rank()));
  return container::svector<container::map<Index, Index>>{idx_rep};
}

ExprPtr expand_P_op(const Product& product) {
  bool has_P_operator = false;

  // Check P and build a replacement map
  // Assuming a product can have multiple P operators
  container::svector<container::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto P = term->as<Tensor>();
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

  if (!has_P_operator) return std::make_shared<Product>(product);

  auto result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    Product new_product{};
    new_product.scale(product.scalar());
    auto temp_product = remove_tensor(product, L"P");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_tensor.reset_tags();
        new_product.append(1, ex<Tensor>(new_tensor));
      }
    }
    result->append(ex<Product>(new_product));
  }  // map_list

  return result;
}

ExprPtr expand_P_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>())
    return expr;
  else if (expr->is<Product>())
    return expand_P_op(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& summand : *expr) {
      result->append(expand_P_op(summand));
    }
    return result;
  } else
    return nullptr;
}

container::svector<container::map<Index, Index>> S_replacement_maps(
    const Tensor& S) {
  assert(S.label() == L"S");
  assert(S.bra_rank() > 1);
  assert(S.bra().size() == S.ket().size());
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
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  auto result = std::make_shared<Sum>();

  // Check if S operator is present
  if (!has_tensor(expr, L"S")) return expr;

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>()) expr->as<Tensor>().reset_tags();
  };
  expr->visit(reset_idx_tags);

  // Lambda for applying S on products
  auto expand_S_product = [](const Product& product) {
    // check if S is present
    if (!has_tensor(ex<Product>(product), L"S")) return ex<Product>(product);

    container::svector<container::map<Index, Index>> maps;
    if (product.factor(0)->as<Tensor>().label() == L"S")
      maps = S_replacement_maps(product.factor(0)->as<Tensor>());
    assert(!maps.empty());
    Sum sum{};
    for (auto&& map : maps) {
      Product new_product{};
      new_product.scale(product.scalar());
      auto temp_product = remove_tensor(product, L"S")->as<Product>();
      for (auto&& term : temp_product) {
        if (term->is<Tensor>()) {
          auto new_tensor = term->as<Tensor>();
          new_tensor.transform_indices(map);
          new_product.append(1, ex<Tensor>(new_tensor));
        }
      }
      sum.append(ex<Product>(new_product));
    }
    ExprPtr result = std::make_shared<Sum>(sum);
    return result;
  };

  if (expr->is<Product>()) {
    result->append(expand_S_product(expr->as<Product>()));
  } else if (expr->is<Sum>()) {
    for (auto&& term : *expr) {
      if (term->is<Product>()) {
        result->append(expand_S_product(term->as<Product>()));
      } else if (term->is<Tensor>() || term->is<Constant>()) {
        result->append(term);
      }
    }
  }

  result->visit(reset_idx_tags);
  return result;
}

ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<Index>>& ext_index_groups) {
  // Symmetrize and expression
  // Partially expand the antisymmetrizer and write it in terms of S operator.
  // See symmetrize_expr(expr) function for implementation details. We want an
  // expression with non-symmetric tensors, hence we are partially expanding the
  // antisymmetrizer (A) and fully expanding the anti-symmetric tensors to
  // non-symmetric.
  auto symm_and_expand = [](const ExprPtr& expr) {
    auto temp = expr;
    if (has_tensor(temp, L"A")) temp = symmetrize_expr(temp);
    temp = expand_antisymm(temp);
    rapid_simplify(temp);
    return temp;
  };
  auto expr = symm_and_expand(expression);

  // Index tags are cleaned prior to calling the fast canonicalizer
  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>()) expr->as<Tensor>().reset_tags();
  };

  // Cleanup index tags
  expr->visit(reset_idx_tags);  // This call is REQUIRED
  expand(expr);                 // This call is REQUIRED
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
    if (product.factor(0)->as<Tensor>().label() == L"S") {
      for (auto&& term : product.factors()) {
        if (term->is<Tensor>() && term->as<Tensor>().label() != L"S")
          temp_product.append(1, term, Product::Flatten::No);
      }
    } else {
      temp_product = product;
    }

    auto get_ket_indices = [](const Product& prod) {
      container::svector<Index> ket_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>())
          ranges::for_each(t->as<Tensor>().ket(), [&ket_idx](const Index& idx) {
            ket_idx.push_back(idx);
          });
      }
      return ket_idx;
    };
    auto product_kets = get_ket_indices(temp_product);

    auto get_bra_indices = [](const Product& prod) {
      container::svector<Index> bra_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>())
          ranges::for_each(t->as<Tensor>().bra(), [&bra_idx](const Index& idx) {
            bra_idx.push_back(idx);
          });
      }
      return bra_idx;
    };
    auto product_bras = get_bra_indices(temp_product);

    auto substitute_ext_idx = [&product_bras, &product_kets](
                                  const container::svector<Index>& idx_pair) {
      assert(idx_pair.size() == 2);
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

  if (expr->is<Constant>())
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
      } else  // summand->is<Constant>()
        result->append(summand);
    }
    return result;
  } else
    return nullptr;
}

container::svector<container::svector<Index>> external_indices(
    const ExprPtr& expr) {
  // Generate external index list from the projection manifold operator
  // (symmetrizer or antisymmetrizer)
  Tensor P{};
  for (const auto& prod : *expr) {
    if (prod->is<Product>()) {
      auto tensor = prod->as<Product>().factor(0)->as<Tensor>();
      if (tensor.label() == L"A" || tensor.label() == L"S") {
        P = tensor;
        break;
      }
    }
  }

  container::svector<container::svector<Index>> ext_index_groups;
  if (P) {  // if have the projection manifold operator
    assert(P.bra_rank() != 0 &&
           "Could not generate external index groups due to "
           "absence of (anti)symmetrizer (A or S) operator in expression.");
    assert(P.bra_rank() == P.ket_rank());
    ext_index_groups.resize(P.rank());
    for (std::size_t i = 0; i != P.rank(); ++i) {
      ext_index_groups[i] = container::svector<Index>{P.ket()[i], P.bra()[i]};
    }
  }
  return ext_index_groups;
}

container::svector<container::svector<Index>> external_indices(
    Tensor const& t) {
  using ranges::views::transform;
  using ranges::views::zip;

  assert(t.label() == L"S" || t.label() == L"A");
  assert(t.bra_rank() == t.ket_rank());
  return zip(t.ket(), t.bra()) | transform([](auto const& pair) {
           return container::svector<Index>{pair.first, pair.second};
         }) |
         ranges::to<container::svector<container::svector<Index>>>;
}

ExprPtr closed_shell_CC_spintrace(ExprPtr const& expr) {
  assert(expr->is<Sum>());
  using ranges::views::transform;

  auto const ext_idxs = external_indices(expr);
  auto st_expr = closed_shell_spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Remove S operator
    for (auto& term : *st_expr) {
      if (term->is<Product>()) term = remove_tensor(term->as<Product>(), L"S");
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

container::svector<ResultExpr> closed_shell_spintrace(ResultExpr expr) {
  bool searchForNonEquivalentResults = expr.symmetry() != Symmetry::nonsymm;
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
    expr.expression() = closed_shell_spintrace(
        expr.expression(),
        expr.index_particle_grouping<container::svector<Index>>());

    expr.set_symmetry(Symmetry::nonsymm);

    return {std::move(expr)};
  }

  assert(expr.symmetry() == Symmetry::antisymm ||
         expr.symmetry() == Symmetry::symm);

  // TODO: Do we have to track the sign?
  const bool permuteBra = expr.bra().size() >= expr.ket().size();
  auto permIndices = permuteBra ? expr.bra() : expr.ket();
  const std::size_t unchangedSize =
      permuteBra ? expr.ket().size() : expr.bra().size();

  auto get_phase = [](auto container) {
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

  container::svector<ResultExpr> resultSet;

  int sign = initialSign;
  do {
    const int currentSign = sign;
    // std::next_permutation creates one lexicographical permutation after the
    // other, which should imply that the phase should alternate between
    // iterations.
    sign *= -1;
    assert(currentSign == get_phase(permIndices) * initialSign);

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

    container::map<Index, Index> idxReplacements;

    for (std::size_t i = 0; i < permIndices.size(); ++i) {
      if (permIndices[i] != originalIndices[i]) {
        idxReplacements.insert({originalIndices[i], permIndices[i]});
        idxReplacements.insert({permIndices[i], originalIndices[i]});
      }
    }

    ExprPtr expression = expr.expression().clone();
    // expression->visit([&](ExprPtr &expr) { if (expr.is<Tensor>())
    // expr.as<Tensor>().transform_indices(idxReplacements); }, true);
    expression->visit(
        [](ExprPtr& expr) {
          if (expr.is<Tensor>()) expr.as<Tensor>().reset_tags();
        },
        true);

    expression *= ex<Constant>(currentSign);
    expression = simplify(expression);

    ResultExpr result = [&]() {
      assert(expr.has_label());
      if (permuteBra) {
        return ResultExpr(permIndices, expr.ket(), expr.aux(), expr.symmetry(),
                          expr.braket_symmetry(), expr.particle_symmetry(),
                          expr.label(), std::move(expression));
      } else {
        return ResultExpr(expr.bra(), permIndices, expr.aux(), expr.symmetry(),
                          expr.braket_symmetry(), expr.particle_symmetry(),
                          expr.label(), std::move(expression));
      }
    }();

    result.expression() = closed_shell_spintrace(
        result.expression(),
        result.index_particle_grouping<container::svector<Index>>());

    result.set_symmetry(Symmetry::nonsymm);

    resultSet.push_back(std::move(result));
  } while (std::next_permutation(permIndices.begin(), permIndices.end()));

  return resultSet;
}

ExprPtr closed_shell_CC_spintrace_rigorous(ExprPtr const& expr) {
  assert(expr->is<Sum>());
  using ranges::views::transform;

  auto const ext_idxs = external_indices(expr);
  auto st_expr = sequant::spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Remove S operator
    for (auto& term : *st_expr) {
      if (term->is<Product>()) term = remove_tensor(term->as<Product>(), L"S");
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

/// Collect all indices from an expression
auto index_list(const ExprPtr& expr) {
  container::set<Index, Index::LabelCompare> grand_idxlist;
  if (expr->is<Tensor>()) {
    ranges::for_each(expr->as<Tensor>().const_braket(),
                     [&grand_idxlist](const Index& idx) {
                       idx.reset_tag();
                       grand_idxlist.insert(idx);
                     });
  }

  return grand_idxlist;
}

Tensor swap_spin(const Tensor& t) {
  auto is_any_spin = [](const Index& i) {
    return mbpt::to_spin(i.space().qns()) == mbpt::Spin::any;
  };

  // Return tensor if there are no spin labels
  if (std::all_of(t.const_braket().begin(), t.const_braket().end(),
                  is_any_spin)) {
    return t;
  }

  // Return new index where the spin-label is flipped
  auto spin_flipped_idx = [](const Index& idx) {
    assert(mbpt::to_spin(idx.space().qns()) != mbpt::Spin::any);
    return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha
               ? make_spinbeta(idx)
               : make_spinalpha(idx);
  };

  container::svector<Index> b(t.rank()), k(t.rank());

  for (std::size_t i = 0; i < t.rank(); ++i) {
    b.at(i) = spin_flipped_idx(t.bra().at(i));
    k.at(i) = spin_flipped_idx(t.ket().at(i));
  }

  return {t.label(),    bra(std::move(b)),   ket(std::move(k)),    t.aux(),
          t.symmetry(), t.braket_symmetry(), t.particle_symmetry()};
}

ExprPtr swap_spin(const ExprPtr& expr) {
  if (expr->is<Constant>()) return expr;

  auto swap_tensor = [](const Tensor& t) { return ex<Tensor>(swap_spin(t)); };

  auto swap_product = [&swap_tensor](const Product& p) {
    Product result{};
    result.scale(p.scalar());
    for (auto& t : p) {
      assert(t->is<Tensor>());
      result.append(1, swap_tensor(t->as<Tensor>()), Product::Flatten::No);
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
  } else
    return nullptr;
}

ExprPtr merge_tensors(const Tensor& O1, const Tensor& O2) {
  assert(O1.label() == O2.label());
  assert(O1.symmetry() == O2.symmetry());
  auto b = ranges::views::concat(O1.bra(), O2.bra());
  auto k = ranges::views::concat(O1.ket(), O2.ket());
  auto a = ranges::views::concat(O1.aux(), O2.aux());
  return ex<Tensor>(Tensor(O1.label(), bra(b), ket(k), aux(a), O1.symmetry()));
}

std::vector<ExprPtr> open_shell_A_op(const Tensor& A) {
  assert(A.label() == L"A");
  assert(A.bra_rank() == A.ket_rank());
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
        Tensor(L"A", spin_bra, spin_ket, A.aux(), Symmetry::antisymm));
    // std::wcout << to_latex(result.at(i)) << " ";
  }
  // std::wcout << "\n" << std::endl;
  return result;
}

std::vector<ExprPtr> open_shell_P_op_vector(const Tensor& A) {
  assert(A.label() == L"A");

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
              L"P", bra{A.bra().at(j), A.bra().at(k)}, ket{}, Symmetry::nonsymm,
              BraKetSymmetry::nonsymm, ParticleSymmetry::nonsymm));
          P_ket_list.emplace_back(Tensor(
              L"P", bra{}, ket{A.ket().at(j), A.ket().at(k)}, Symmetry::nonsymm,
              BraKetSymmetry::nonsymm, ParticleSymmetry::nonsymm));
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
                         ket{}, Symmetry::nonsymm, BraKetSymmetry::nonsymm,
                         ParticleSymmetry::nonsymm));
              P_ket_list.emplace_back(
                  Tensor(L"P", bra{},
                         ket{A.ket().at(i1), A.ket().at(i3), A.ket().at(i2),
                             A.ket().at(i4)},
                         Symmetry::nonsymm, BraKetSymmetry::nonsymm,
                         ParticleSymmetry::nonsymm));
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
    const int single_spin_case) {
  if (expr->is<Constant>()) {
    return std::vector<ExprPtr>{expr};
  }

  // Grand index list contains both internal and external indices
  container::set<Index, Index::LabelCompare> grand_idxlist;
  auto collect_indices = [&grand_idxlist](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [&grand_idxlist](const Index& idx) {
                         idx.reset_tag();
                         grand_idxlist.insert(idx);
                       });
    }
  };
  expr->visit(collect_indices);

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

  assert(grand_idxlist.size() == int_idxlist.size() + ext_idxlist.size());

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
        assert((spin_bit == 0) || (spin_bit == 1));
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

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>()) expr->as<Tensor>().reset_tags();
  };

  // Internal and external index replacements are independent
  auto i_rep = spin_cases(int_index_groups);
  auto e_rep = ext_spin_cases(ext_index_groups);

  // For a single spin case, keep only the relevant spin case
  // PS: all alpha indexing start at 0
  if (single_spin_case) {
    auto external_replacement_map = e_rep.at(single_spin_case);
    e_rep.clear();
    e_rep.push_back(external_replacement_map);
  }

  // Expand 'A' operator and 'antisymm' tensors
  auto expanded_expr = expand_A_op(expr);
  expanded_expr->visit(reset_idx_tags);
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
      }
    }
    assert(cKet.size() == cBra.size());

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
    spin_expr->visit(reset_idx_tags);
    Sum e_result{};

    // Loop over internal index replacement maps
    for (auto& i : i_rep) {
      // Add spin labels to internal indices, expand antisymmetric tensors
      auto spin_expr_i = append_spin(spin_expr, i);
      spin_expr_i = expand_antisymm(spin_expr_i, true);
      expand(spin_expr_i);
      spin_expr_i->visit(reset_idx_tags);
      Sum i_result{};

      if (spin_expr_i->is<Tensor>()) {
        e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Product>()) {
        if (spin_symm_product(spin_expr_i->as<Product>()))
          e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Sum>()) {
        for (auto& pr : *spin_expr_i) {
          if (pr->is<Product>()) {
            if (spin_symm_product(pr->as<Product>())) i_result.append(pr);
          } else if (pr->is<Tensor>()) {
            if (spin_symm_tensor(pr->as<Tensor>())) i_result.append(pr);
          } else if (pr->is<Constant>()) {
            i_result.append(pr);
          } else
            throw("Unknown ExprPtr type.");
        }
        e_result.append(std::make_shared<Sum>(i_result));
      }

    }  // loop over internal indices
    result.push_back(std::make_shared<Sum>(e_result));
  }  // loop over external indices

  if (single_spin_case) {
    assert(result.size() == 1 &&
           "Spin-specific case must return one expression.");
  }

  // Canonicalize and simplify all expressions
  for (auto& expression : result) {
    expression->visit(reset_idx_tags);
    canonicalize(expression);
    rapid_simplify(expression);
  }
  return result;
}

std::vector<ExprPtr> open_shell_CC_spintrace(const ExprPtr& expr) {
  Tensor A = expr->at(0)->at(0)->as<Tensor>();
  assert(A.label() == L"A");
  size_t const i = A.rank();
  auto P_vec = open_shell_P_op_vector(A);
  auto A_vec = open_shell_A_op(A);
  assert(P_vec.size() == i + 1);
  std::vector<Sum> concat_terms(i + 1);
  [[maybe_unused]] size_t n_spin_orbital_term = 0;
  for (auto& product_term : *expr) {
    auto term = remove_tensor(product_term->as<Product>(), L"A");
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
    container::svector<container::svector<Index>> ext_index_groups,
    bool spinfree_index_spaces) {
  // Escape immediately if expression is a constant
  if (expression->is<Constant>()) {
    return expression;
  }

  // This function must be used for tensors with spin-specific indices only. If
  // the spin-symmetry is conserved: the tensor is expanded; else: zero is
  // returned.
  auto spin_trace_tensor = [](const Tensor& tensor) {
    return can_expand(tensor) ? expand_antisymm(tensor) : ex<Constant>(0);
  };

  // This function is used to spin-trace a product terms with spin-specific
  // indices. It checks if all tensors can be expanded and spintraces individual
  // tensors by call to the spin_trace_tensor lambda.
  auto spin_trace_product = [&spin_trace_tensor](const Product& product) {
    Product spin_product{};

    // Check if all tensors in this product can be expanded
    // If NOT all tensors can be expanded, return zero
    if (!std::all_of(product.factors().begin(), product.factors().end(),
                     [](const auto& t) {
                       return can_expand(t->template as<Tensor>());
                     })) {
      return ex<Constant>(0);
    }

    spin_product.scale(product.scalar());
    ranges::for_each(product.factors().begin(), product.factors().end(),
                     [&spin_trace_tensor, &spin_product](const auto& t) {
                       spin_product.append(
                           1, spin_trace_tensor(t->template as<Tensor>()));
                     });

    ExprPtr result = std::make_shared<Product>(spin_product);
    expand(result);
    rapid_simplify(result);
    return result;
  };

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>()) expr->as<Tensor>().reset_tags();
  };

  // Most important lambda of this function
  auto trace_product = [&ext_index_groups, &spin_trace_tensor,
                        &spin_trace_product,
                        spinfree_index_spaces](const Product& expression) {
    ExprPtr expr = std::make_shared<Product>(expression);

    // List of all indices in the expression
    container::set<Index, Index::LabelCompare> grand_idxlist;
    auto collect_indices = [&grand_idxlist](const ExprPtr& expr) {
      if (expr->is<Tensor>()) {
        ranges::for_each(expr->as<Tensor>().const_braket(),
                         [&grand_idxlist](const Index& idx) {
                           idx.reset_tag();
                           grand_idxlist.insert(idx);
                         });
      }
    };
    expr->visit(collect_indices);

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
    const uint64_t nspincases = pow2(index_groups.size());

    auto result = std::make_shared<Sum>();
    for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
         ++spincase_bitstr) {
      // EFV:  assign spin to each index group => make a replacement list
      container::map<Index, Index> index_replacements;

      uint64_t index_group_count = 0;
      for (auto&& index_group : index_groups) {
        auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
        assert(spin_bit == 0 || spin_bit == 1);

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
        auto st_expr = spin_trace_tensor(spin_expr->as<Tensor>());
        result->append(spinfree_index_spaces ? remove_spin(st_expr) : st_expr);
      } else if (spin_expr->is<Product>()) {
        auto st_expr = spin_trace_product(spin_expr->as<Product>());
        if (st_expr->size() != 0) {
          result->append(spinfree_index_spaces ? remove_spin(st_expr)
                                               : st_expr);
        }
      } else if (spin_expr->is<Sum>()) {
        for (auto&& summand : *spin_expr) {
          Sum st_expr{};
          if (summand->is<Tensor>())
            st_expr.append(spin_trace_tensor(summand->as<Tensor>()));
          else if (summand->is<Product>())
            st_expr.append(spin_trace_product(summand->as<Product>()));
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
    result = trace_product(expr->as<Product>());
  } else if ((expr->is<Sum>())) {
    auto result_sum = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      if (term->is<Product>())
        result_sum->append(trace_product(term->as<Product>()));
      else if (term->is<Tensor>()) {
        auto term_as_product = ex<Constant>(1) * term;
        result_sum->append(trace_product(term_as_product->as<Product>()));
      } else
        result_sum->append(term);
      result = result_sum;
    }
    return result;
  } else
    return nullptr;
  result->visit(reset_idx_tags);
  return result;
}  // ExprPtr spintrace

ResultExpr spintrace(ResultExpr expr, bool spinfree_index_spaces) {
  expr.expression() =
      spintrace(expr.expression(),
                expr.index_particle_grouping<container::svector<Index>>(),
                spinfree_index_spaces);

  expr.set_symmetry(Symmetry::nonsymm);

  return expr;
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
    assert(bra_list.size() == ket_list.size());
    S = Tensor(L"S", bra(std::move(bra_list)), ket(std::move(ket_list)),
               Symmetry::nonsymm);
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
    assert(summands_hash_list.size() == expr->size());
    assert(summands_hash_map.size() == expr->size());

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
      // FOR GENERALIZED EXPRESSION WITH ARBITRATY S OPERATOR
      // Loop over replacement maps The entire code from here
      container::vector<size_t> hash1_list;
      for (auto&& replacement_map : replacement_maps) {
        size_t hash1;
        if ((*it)->is<Product>()) {
          // Clone *it, apply symmetrizer, store hash1 value
          auto product = (*it)->as<Product>();
          Product S_product{};
          S_product.scale(product.scalar());

          // Transform indices by action of S operator
          for (auto&& t : product) {
            if (t->is<Tensor>())
              S_product.append(
                  1, transform_tensor(t->as<Tensor>(), replacement_map),
                  Product::Flatten::No);
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
        size_t hash0;
        if ((*it)->is<Product>()) {
          // Clone *it, apply symmetrizer, store hash value
          auto product = (*it)->as<Product>();
          Product S_product{};
          S_product.scale(product.scalar());

          // Transform indices by action of S operator
          for (auto&& t : product) {
            if (t->is<Tensor>())
              S_product.append(
                  1, transform_tensor(t->as<Tensor>(), replacement_map),
                  Product::Flatten::No);
          }
          auto new_product_expr = ex<Product>(S_product);
          new_product_expr->canonicalize();
          hash0 = new_product_expr->hash_value();
          hash0_list.push_back(hash0);
        } else if ((*it)->is<Tensor>()) {
          // Clone *it, apply symmetrizer, store hash value
          auto tensor = (*it)->as<Tensor>();

          // Transform indices by action of S operator
          auto new_tensor = transform_tensor(tensor, replacement_map);

          // Canonicalize the new tensor before computing hash value
          new_tensor->canonicalize();
          hash0 = new_tensor->hash_value();
          hash0_list.push_back(hash0);
        }
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

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups,
    const double threshold) {
  assert(!ext_index_groups.empty());
  const auto n_particles = ext_index_groups.size();

  using sequant::container::svector;

  // Coefficients
  container::svector<rational> bt_coeff_vec;
  {
#ifdef SEQUANT_HAS_EIGEN
    using namespace Eigen;
    // Dimension of permutation matrix is n_particles!
    const auto n = boost::numeric_cast<Eigen::Index>(factorial(n_particles));

    // Permutation matrix
    Eigen::Matrix<double, Dynamic, Dynamic> M(n, n);
    {
      M.setZero();
      size_t n_row = 0;
      svector<int> v(n_particles), v1(n_particles);
      std::iota(v.begin(), v.end(), 0);
      std::iota(v1.begin(), v1.end(), 0);
      do {
        container::svector<double> permutation_vector;
        do {
          permutation_vector.push_back(
              std::pow(-2, sequant::count_cycles(v1, v)));
        } while (std::next_permutation(v.begin(), v.end()));
        Eigen::VectorXd pv_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
            permutation_vector.data(), permutation_vector.size());
        M.row(n_row) = pv_eig;
        ++n_row;
      } while (std::next_permutation(v1.begin(), v1.end()));
      M *= std::pow(-1, n_particles);
    }

    // Normalization constant
    double scalar;
    {
      auto nonZero = [&threshold](const double& d) {
        using std::abs;
        return abs(d) > threshold;
      };

      // Solve system of equations
      SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
      container::svector<double> eig_vals(eig_solver.eigenvalues().size());
      VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
          eig_solver.eigenvalues();

      double non0count =
          std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
      scalar = eig_vals.size() / non0count;
    }

    // Find Pseudo Inverse, get 1st row only
    MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
    container::svector<double> bt_coeff_dvec;
    bt_coeff_dvec.resize(pinv.rows());
    VectorXd::Map(&bt_coeff_dvec[0], bt_coeff_dvec.size()) =
        pinv.row(0) * scalar;
    bt_coeff_vec.reserve(bt_coeff_dvec.size());
    ranges::for_each(bt_coeff_dvec, [&bt_coeff_vec, threshold](double c) {
      bt_coeff_vec.emplace_back(to_rational(c, threshold));
    });

//    std::cout << "n_particles = " << n_particles << "\n bt_coeff_vec = ";
//    std::copy(bt_coeff_vec.begin(), bt_coeff_vec.end(),
//              std::ostream_iterator<rational>(std::cout, " "));
//    std::cout << "\n";
#else
    // hardwire coefficients for n_particles = 1, 2, 3
    switch (n_particles) {
      case 1:
        bt_coeff_vec = {ratio(1, 2)};
        break;
      case 2:
        bt_coeff_vec = {ratio(1, 3), ratio(1, 6)};
        break;
      case 3:
        bt_coeff_vec = {ratio(17, 120), ratio(-1, 120), ratio(-1, 120),
                        ratio(-7, 120), ratio(-7, 120), ratio(-1, 120)};
        break;
      default:
        throw std::runtime_error(
            "biorthogonal_transform requires Eigen library for n_particles > "
            "3.");
    }
#endif
  }

  // Transformation maps
  container::svector<container::map<Index, Index>> bt_maps;
  {
    container::svector<Index> idx_list(ext_index_groups.size());

    for (std::size_t i = 0; i != ext_index_groups.size(); ++i) {
      idx_list[i] = *ext_index_groups[i].begin();
    }

    const container::svector<Index> const_idx_list = idx_list;

    do {
      container::map<Index, Index> map;
      auto const_list_ptr = const_idx_list.begin();
      for (auto& i : idx_list) {
        map.emplace(*const_list_ptr, i);
        const_list_ptr++;
      }
      bt_maps.push_back(map);
    } while (std::next_permutation(idx_list.begin(), idx_list.end()));
  }

  // If this assertion fails, change the threshold parameter
  assert(bt_coeff_vec.size() == bt_maps.size());

  // Checks if the replacement map is a canonical sequence
  auto is_canonical = [](const container::map<Index, Index>& idx_map) {
    bool canonical = true;
    for (auto&& pair : idx_map)
      if (pair.first != pair.second) return false;
    return canonical;
  };

  // Scale transformed expressions and append
  Sum bt_expr{};
  auto coeff_it = bt_coeff_vec.begin();
  for (auto&& map : bt_maps) {
    const auto v = *coeff_it;
    if (is_canonical(map))
      bt_expr.append(ex<Constant>(v) * expr->clone());
    else
      bt_expr.append(ex<Constant>(v) *
                     sequant::transform_expr(expr->clone(), map));
    coeff_it++;
  }
  ExprPtr result = std::make_shared<Sum>(bt_expr);
  return result;
}

}  // namespace sequant
