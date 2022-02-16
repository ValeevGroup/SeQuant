#include "spin.hpp"

#include <SeQuant/core/tensor.hpp>
#include <unordered_map>

namespace sequant {

ExprPtr transform_expression(const ExprPtr& expr,
                             const std::map<Index, Index>& index_replacements,
                             double scaling_factor) {
  if (expr->is<Constant>()) {
    return ex<Constant>(scaling_factor) * expr;
  }

  auto transform_tensor = [&index_replacements](const Tensor& tensor) {
    auto result = std::make_shared<Tensor>(tensor);
    result->transform_indices(index_replacements);
    ranges::for_each(result->const_braket(),
                     [](const Index& idx) { idx.reset_tag(); });
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
    for (auto&& term : *expr) {
      if (term->is<Constant>())
        result->append(ex<Constant>(scaling_factor) * term);
      else if (term->is<Tensor>()) {
        auto transformed_tensor =
            ex<Constant>(scaling_factor) * transform_tensor(term->as<Tensor>());
        result->append(transformed_tensor);
      } else if (term->is<Product>()) {
        auto transformed_product = transform_product(term->as<Product>());
        result->append(transformed_product);
      }
    }
    return result;
  } else
    return nullptr;
}

ExprPtr swap_bra_ket(const ExprPtr& expr) {
  if (expr->is<Constant>()) return expr;

  // Lambda for tensor
  auto tensor_swap = [](const Tensor& tensor) {
    auto result =
        Tensor(tensor.label(), tensor.ket(), tensor.bra(), tensor.symmetry(),
               tensor.braket_symmetry(), tensor.particle_symmetry());
    return ex<Tensor>(result);
  };

  // Lambda for product
  auto product_swap = [&tensor_swap](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) result->append(tensor_swap(term->as<Tensor>()));
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
      if (term->is<Product>())
        result->append(product_swap(term->as<Product>()));
      else if (term->is<Tensor>())
        result->append(tensor_swap(term->as<Tensor>()));
      else if (term->is<Constant>())
        result->append(term);
    }
    return result;
  } else
    throw("Unknown expression type.");
}

ExprPtr append_spin(ExprPtr& expr,
                    const std::map<Index, Index>& index_replacements) {
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
      if (summand->is<Tensor>()) {
        spin_expr->append(add_spin_to_tensor(summand->as<Tensor>()));
      } else if (summand->is<Product>()) {
        spin_expr->append(add_spin_to_product(summand->as<Product>()));
      } else {
        spin_expr->append(summand);
      }
    }
    return spin_expr;
  } else
    return nullptr;
}

ExprPtr remove_spin(ExprPtr& expr) {

  auto remove_spin_from_tensor = [](const Tensor& tensor) {
    container::svector<Index> bra(tensor.bra().begin(),
                                  tensor.bra().end());
    container::svector<Index> ket(tensor.ket().begin(),
                                  tensor.ket().end());
    {

      for (auto&& idx : ranges::views::concat(bra, ket)) {
        auto space = IndexSpace::instance(idx.space().type(),
                                          IndexSpace::nullqns);
        auto subscript_label = idx.label().substr(idx.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(),
                                        subscript_label.end());
        idx = Index::make_label_index(space, subscript_label_ws);
      }
    }
    Tensor result(tensor.label(), bra, ket, tensor.symmetry(),
                      tensor.braket_symmetry());
    return std::make_shared<Tensor>(result);
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
      if (summand->is<Product>()) {
        result->append(remove_spin_from_product(summand->as<Product>()));
      } else if (summand->is<Tensor>()) {
        result->append(remove_spin_from_tensor(summand->as<Tensor>()));
      } else {
        result->append(summand);
      }
    }
    return result;
  } else
    return expr;
}

bool is_tensor_spin_symm(const Tensor& tensor) {
  assert(tensor.bra_rank() == tensor.ket_rank());
  for(auto i = 0; i != tensor.rank(); ++i){
    if(tensor.bra()[i].space().qns() != tensor.ket()[i].space().qns())
       return false;
  }
  return true;
}

bool same_spin_tensor(const Tensor& tensor){
  auto braket = tensor.braket();
  auto spin_element = braket[0].space().qns();
  return std::all_of(braket.begin(), braket.end(),
                     [&spin_element](const auto& idx){
    return idx.space().qns() == spin_element;
  });
}

bool can_expand(const Tensor& tensor) {
  assert(tensor.bra_rank() == tensor.ket_rank() && "can_expand(Tensor) failed.");
  if (tensor.bra_rank() != tensor.ket_rank())
    return false;

  // indices with non-qns are not allowed
  for(auto& idx : tensor.const_braket()){
    assert(idx.space().qns() != IndexSpace::nullqns);
  }

  // count alpha indices in bra
  int a_bra = std::count_if(tensor.bra().begin(), tensor.bra().end(),
                            [](const Index& idx){
                              return idx.space().qns() == IndexSpace::alpha;
                            });

  // count alpha indices in ket
  int a_ket = std::count_if(tensor.ket().begin(), tensor.ket().end(),
                            [](const Index& idx){
                              return idx.space().qns() == IndexSpace::alpha;
                            });

  return a_bra == a_ket;
}

ExprPtr expand_antisymm(const Tensor& tensor, bool skip_spinsymm) {
  assert(tensor.bra_rank() == tensor.ket_rank());
  if (tensor.bra_rank() == 1) {
    Tensor new_tensor(tensor.label(), tensor.bra(), tensor.ket(),
                      Symmetry::nonsymm, tensor.braket_symmetry(),
                      tensor.particle_symmetry());
    return std::make_shared<Tensor>(new_tensor);
  }

  // If all indices have the same spin label,
  // return the antisymm tensor
  if(skip_spinsymm) {
    if (same_spin_tensor(tensor))
      return std::make_shared<Tensor>(tensor);
  }

  assert(tensor.bra_rank() > 1);

  auto get_phase = [](const Tensor& t) {
    container::svector<Index> bra(t.bra().begin(), t.bra().end());
    container::svector<Index> ket(t.ket().begin(), t.ket().end());
    IndexSwapper::thread_instance().reset();
    bubble_sort(std::begin(bra), std::end(bra), std::less<Index>{});
    bubble_sort(std::begin(ket), std::end(ket), std::less<Index>{});
    return IndexSwapper::thread_instance().even_num_of_swaps() ? 1 : -1;
  };

  // Generate a sum of asymmetric tensors if the input tensor is antisymmetric
  // and greater than one body otherwise, return the tensor
  if (tensor.symmetry() == Symmetry::antisymm) {
    const auto prefactor = get_phase(tensor);
    container::set<Index> bra_list;
    for (auto& bra_idx : tensor.bra()) bra_list.insert(bra_idx);
    const auto const_bra_list = bra_list;

    container::set<Index> ket_list;
    for (auto&& ket_idx : tensor.ket()) ket_list.insert(ket_idx);

    Sum expr_sum{};
    do {
      auto bra_list2 = bra_list;
      auto new_tensor =
          Tensor(tensor.label(), bra_list, ket_list, Symmetry::nonsymm);

      if (is_tensor_spin_symm(new_tensor)) {
        const auto phase = get_phase(new_tensor);
        auto new_tensor_ptr = ex<Tensor>(new_tensor);
        Product new_tensor_product{};
        new_tensor_product.append(phase, new_tensor_ptr);
        new_tensor_product.scale(prefactor);
        auto new_tensor_product_ptr = ex<Product>(new_tensor_product);
        expr_sum.append(new_tensor_product_ptr);
      }
    } while (std::next_permutation(bra_list.begin(), bra_list.end()));

    return std::make_shared<Sum>(expr_sum);
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
      if (term->is<Tensor>()) temp.append(expand_antisymm(term->as<Tensor>(), skip_spinsymm));
    }
    ExprPtr result = std::make_shared<Product>(temp);
    rapid_simplify(result);
    return result;
  };

  if (expr->is<Product>())
    return expand_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    Sum sum_result{};
    for (auto&& term : *expr) {
      if (term->is<Product>())
        sum_result.append(expand_product(term->as<Product>()));
      else if (term->is<Tensor>())
        sum_result.append(expand_antisymm(term->as<Tensor>(), skip_spinsymm));
      else if (term->is<Constant>())
        sum_result.append(term);
      else
        sum_result.append(nullptr);
    }
    return std::make_shared<Sum>(sum_result);
  } else
    return nullptr;
}

bool has_A_label(const ExprPtr& expr) {
  if (expr->is<Constant>()) return false;

  if (expr->is<Tensor>())
    return expr->as<Tensor>().label() == L"A";
  else if (expr->is<Product>())
    return (expr->as<Product>().factor(0))->as<Tensor>().label() == L"A";
  else if (expr->is<Sum>()) {
    bool result = false;
    for (auto&& term : *expr) {
      if (term->is<Product>())
        result = (term->as<Product>().factor(0))->as<Tensor>().label() == L"A";
      if (term->is<Tensor>())
        result = expr->as<Tensor>().label() == L"A";
      if (result) return result;
    }
    return result;
  } else
    throw("control reached end of check_A_label function.");  // TODO: Catch
                                                              // exception
}

bool has_tensor_label(const ExprPtr& expr, std::wstring label) {
  if (expr->is<Constant>()) return false;

  auto check_product = [&label](const Product& product) {
    bool result = false;
    for (auto& term : product) {
      if (term->is<Tensor>()) {
        if (term->as<Tensor>().label() == label)
          return true;
      }
    }
    return result;
  };

  if (expr->is<Tensor>()) {
    return expr->as<Tensor>().label() == label;
  } else if (expr->is<Product>()) {
    return check_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    bool result = false;
    for (auto&& term : *expr) {
      if (term->is<Product>()) {
        result = check_product(term->as<Product>());
      } else if (term->is<Tensor>()){
        result = term->as<Tensor>().label() == label;
      }
      if (result)
        return true;
    }
    return result;
  } else
    return false;
}

bool has_operator_label(const ExprPtr& expr, std::wstring label) {
  bool result = false;

  if (expr->is<Constant>())
    return result;
  else if (expr->is<Tensor>())
    return expr->as<Tensor>().label() == label;
  else if (expr->is<Product>()) {
    if (expr->as<Product>().factor(0)->is<Tensor>())
      return expr->as<Product>().factor(0)->as<Tensor>().label() == label;
  } else if (expr->is<Sum>()) {
    for (auto&& term : *expr) {
      if (term->is<Product>()) {
        if ((term->as<Product>().factor(0))->is<Tensor>()) {
          result =
              (term->as<Product>().factor(0))->as<Tensor>().label() == label;
        }
      } else if (term->is<Tensor>()) {
        result = term->as<Tensor>().label() == label;
      }
      if (result) return true;
    }
  }
  return result;
}

std::vector<std::map<Index, Index>> A_replacement_map(const Tensor& A) {
  assert(A.label() == L"A");
  assert(A.bra_rank() > 1);
  assert(A.bra().size() == A.ket().size());
  container::svector<int> bra_int_list(A.bra().size());
  std::iota(std::begin(bra_int_list), std::end(bra_int_list), 0);
  auto ket_int_list = bra_int_list;
  std::vector<std::map<Index, Index>> result;
  container::svector<Index> A_braket(A.const_braket().begin(),
                                     A.const_braket().end());

  do {
    do {
      std::map<Index, Index> replacement_map;
      auto A_braket_ptr = A_braket.begin();
      for (auto&& idx : bra_int_list) {
        replacement_map.insert({*A_braket_ptr, A.bra()[idx]});
        A_braket_ptr++;
      }
      for (auto&& idx : ket_int_list) {
        replacement_map.insert({*A_braket_ptr, A.ket()[idx]});
        A_braket_ptr++;
      }
      result.push_back(replacement_map);
    } while (std::next_permutation(bra_int_list.begin(), bra_int_list.end()));
  } while (std::next_permutation(ket_int_list.begin(), ket_int_list.end()));
  return result;
}

ExprPtr remove_tensor_from_product(const Product& product, std::wstring label) {
  auto new_product = std::make_shared<Product>();
  new_product->scale(product.scalar());
  for (auto&& term : product) {
    if (term->is<Tensor>()) {
      auto tensor = term->as<Tensor>();
      if (tensor.label() != label) new_product->append(1, ex<Tensor>(tensor));
    }
  }
  return new_product;
}

ExprPtr remove_tensor_from_expr(const ExprPtr& expr, std::wstring label){

  if (expr->is<Sum>()){
    Sum result{};
    for(auto& term : *expr){
      result.append(remove_tensor_from_expr(term, label));
    }
    return ex<Sum>(result);
  } else if (expr->is<Product>())
    return remove_tensor_from_product(expr->as<Product>(), label);
  else if (expr->is<Tensor>())
    return expr->as<Tensor>().label() == label ? ex<Constant>(1) : expr;
  else if (expr->is<Constant>())
    return expr;
  else
    return nullptr;
}

ExprPtr expand_A_operator(const Product& product) {
  bool has_A_operator = false;

  // Check A and build replacement map
  std::vector<std::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto A = term->as<Tensor>();
      if ((A.label() == L"A") && (A.bra().size() > 1)) {
        has_A_operator = true;
        map_list = A_replacement_map(A);
        break;
      } else if ((A.label() == L"A") && (A.bra().size() == 1)) {
        return remove_tensor_from_product(product, L"A");
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

      IndexSwapper::thread_instance().reset();
      bubble_sort(std::begin(transformed_list), std::end(transformed_list),
                  std::less<Index>{});
      phase = IndexSwapper::thread_instance().even_num_of_swaps() ? 1 : -1;
    }

    Product new_product{};
    new_product.scale(product.scalar());
    auto temp_product = remove_tensor_from_product(product, L"A");
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
  return new_result;
}

ExprPtr expr_symmetrize(const Product& product) {
  auto result = std::make_shared<Sum>();

  // Assumes canonical sequence of tensors in the product
  if (product.factor(0)->as<Tensor>().label() != L"A")
    return std::make_shared<Product>(product);

  // CHECK: A is present and >1 particle
  // GENERATE S tensor
  auto A_tensor = product.factor(0)->as<Tensor>();
  assert(A_tensor.label() == L"A");

  auto A_is_nconserving = A_tensor.bra_rank() == A_tensor.ket_rank();
  if (A_is_nconserving)
    if (A_tensor.bra_rank() == 1)
      return remove_tensor_from_product(product, L"A");
  assert(A_tensor.rank() > 1);

  auto S = Tensor{};
  if (A_is_nconserving) {
    S = Tensor(L"S", A_tensor.bra(), A_tensor.ket(), Symmetry::nonsymm);
  } else {  // A is not nconserving
    size_t n = (A_tensor.bra_rank() < A_tensor.ket_rank())
                   ? A_tensor.bra_rank()
                   : A_tensor.ket_rank();
    container::svector<Index> bra_list(A_tensor.bra().begin(),
                                       A_tensor.bra().end());
    container::svector<Index> ket_list(A_tensor.ket().begin(),
                                       A_tensor.ket().end());
    S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  }

  // Generate replacement maps from a list(could be a bra or a ket)
  // Uses a permuted list of int to generate permutations
  // TODO factor out for reuse
  auto maps_from_list = [](const container::svector<Index, 4>& list) {
    container::svector<int> int_list(list.size());
    std::iota(int_list.begin(), int_list.end(), 0);
    std::vector<std::map<Index, Index>> result;
    do {
      std::map<Index, Index> map;
      auto list_ptr = list.begin();
      for (auto&& i : int_list) {
        map.insert({*list_ptr, list[i]});
        list_ptr++;
      }
      result.push_back(map);
    } while (std::next_permutation(int_list.begin(), int_list.end()));
    assert(result.size() == std::tgamma(list.size() + 1));
    return result;
  };

  // Get phase relative to the canonical order
  // TODO factor out for reuse
  auto get_phase = [](const std::map<Index, Index>& map) {
    bool even;
    container::svector<Index> idx_list;
    for (const auto& [key, val] : map) idx_list.push_back(val);
    IndexSwapper::thread_instance().reset();
    bubble_sort(std::begin(idx_list), std::end(idx_list), std::less<Index>{});
    return IndexSwapper::thread_instance().even_num_of_swaps() ? 1 : -1;
  };

  std::vector<std::map<Index, Index>> maps;
  // CASE 1: n_bra = n_ket on all tensors
  if (A_is_nconserving) {
    maps = maps_from_list(A_tensor.bra());
  } else {
    assert(A_tensor.bra_rank() != A_tensor.ket_rank());
    if (A_tensor.bra_rank() > A_tensor.ket_rank())
      maps = maps_from_list(A_tensor.bra());
    else
      maps = maps_from_list(A_tensor.ket());
  }
  assert(!maps.empty());
  for (auto&& map : maps) {
    Product new_product{};
    new_product.scale(product.scalar());
    new_product.append(get_phase(map), ex<Tensor>(S));
    auto temp_product = remove_tensor_from_product(product, L"A");
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

ExprPtr expr_symmetrize(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  if (expr->is<Product>())
    return expr_symmetrize(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      if (summand->is<Product>())
        result->append(expr_symmetrize(summand->as<Product>()));
      else
        result->append(summand);
    }
    return result;
  } else
    throw("Unknown arg Type for expr_symmetrize.");  // TODO: Catch exception
}

ExprPtr expand_A_operator(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  if (expr->is<Product>())
    return expand_A_operator(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      if (summand->is<Product>())
        result->append(expand_A_operator(summand->as<Product>()));
      else
        result->append(summand);
    }
    return result;
  } else
    throw("Unknown arg Type for expand_A_operator.");
}

std::vector<std::map<Index, Index>> P_replacement_map(const Tensor& P,
                                                      bool keep_canonical,
                                                      bool pair_wise) {
  assert(P.label() == L"P");

  if(pair_wise){
    // Return pair-wise replacements (this is the preferred method)
    // P_ij -> {{i,j},{j,i}}
    // P_ijkl -> {{i,j},{j,i},{k,l},{l,k}}
    // P_ij^ab -> {{i,j},{j,i},{a,b},{b,a}}
    assert(P.bra_rank() % 2 == 0 && P.ket_rank() % 2 == 0);
    std::map<Index, Index> idx_rep;
    for(auto i = 0; i != P.const_braket().size(); i+=2 ){
      idx_rep.insert({P.const_braket().at(i), P.const_braket().at(i + 1)});
      idx_rep.insert({P.const_braket().at(i+1), P.const_braket().at(i)});
    }

    assert(idx_rep.size() == (P.bra_rank() + P.ket_rank()));
    return std::vector<std::map<Index, Index>>{idx_rep};
  }

  size_t size;
  if (P.bra_rank() == 0 || P.ket_rank() == 0)
    size = std::max(P.bra_rank(), P.ket_rank());
  else
    size = std::min(P.bra_rank(), P.ket_rank());

  container::svector<int> int_list(size);
  std::iota(std::begin(int_list), std::end(int_list), 0);
  std::vector<std::map<Index, Index>> result;
  container::svector<Index> P_braket(P.const_braket().begin(),
                                     P.const_braket().end());

  do {
    auto P_braket_ptr = P_braket.begin();
    std::map<Index, Index> replacement_map;
    for (auto&& i : int_list) {
      if(i < P.bra_rank()) {
        replacement_map.insert({*P_braket_ptr, P.bra()[i]});
        P_braket_ptr++;
      }
    }
    for (auto&& i : int_list) {
      if(i < P.ket_rank()) {
        replacement_map.insert({*P_braket_ptr, P.ket()[i]});
        P_braket_ptr++;
      }
    }
    result.push_back(replacement_map);
  } while (std::next_permutation(int_list.begin(), int_list.end()));

  if(!keep_canonical){
    result.erase(result.begin());
  }
  return result;
}

ExprPtr expand_P_operator(const Product& product,
                          bool keep_canonical, bool pair_wise) {
  bool has_P_operator = false;

  // Check P and build a replacement map
  // Assuming a product can have multiple P operators
  std::vector<std::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto P = term->as<Tensor>();
      if ((P.label() == L"P") && (P.bra_rank() > 1 || (P.ket_rank() > 1))) {
        has_P_operator = true;
        auto map = P_replacement_map(P, keep_canonical, pair_wise);
        map_list.insert(map_list.end(), map.begin(), map.end());
      } else if ((P.label() == L"P") && (P.bra_rank() == 1 && (P.ket_rank() == 1))) {
        return remove_tensor_from_product(product, L"P");
      }
    }
  }

  if (!has_P_operator) return std::make_shared<Product>(product);

  auto result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    Product new_product{};
    new_product.scale(product.scalar());
    auto temp_product = remove_tensor_from_product(product, L"P");
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product.append(1, ex<Tensor>(new_tensor));
      }
    }
    result->append(ex<Product>(new_product));
  }  // map_list
  return result;
}

ExprPtr expand_P_operator(const ExprPtr& expr, bool keep_canonical, bool pair_wise) {
  if (expr->is<Constant>() || expr->is<Tensor>())
    return expr;
  else if (expr->is<Product>())
    return expand_P_operator(expr->as<Product>(), keep_canonical, pair_wise);
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& summand : *expr) {
        result->append(expand_P_operator(summand, keep_canonical, pair_wise));
    }
    return result;
  } else
    throw("Unknown arg Type for expand_P_operator.");
}

std::vector<std::map<Index, Index>> S_replacement_maps(const Tensor& S) {
  assert(S.label() == L"S");
  assert(S.bra_rank() > 1);
  assert(S.bra().size() == S.ket().size());
  container::svector<int> int_list(S.bra().size());
  std::iota(std::begin(int_list), std::end(int_list), 0);

  std::vector<std::map<Index, Index>> maps;
  do {
    std::map<Index, Index> map;
    auto S_bra_ptr = S.bra().begin();
    auto S_ket_ptr = S.ket().begin();
    for (auto&& i : int_list) {
      map.insert({*S_bra_ptr, S.bra()[i]});
      ++S_bra_ptr;
      map.insert({*S_ket_ptr, S.ket()[i]});
      ++S_ket_ptr;
    }
    maps.push_back(map);
  } while (std::next_permutation(int_list.begin(), int_list.end()));

  return maps;
}

ExprPtr expand_S_operator(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Tensor>()) return expr;

  auto result = std::make_shared<Sum>();

  // Check if S operator is present
  if (!has_tensor_label(expr, L"S")) return expr;

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };
  expr->visit(reset_idx_tags);

  // Lambda for applying S on products
  auto expand_S_product = [](const Product& product) {
    // check if S is present
    if (!has_tensor_label(ex<Product>(product), L"S"))
      return ex<Product>(product);

    std::vector<std::map<Index, Index>> maps;
    if (product.factor(0)->as<Tensor>().label() == L"S")
      maps = S_replacement_maps(product.factor(0)->as<Tensor>());
    assert(!maps.empty());
    Sum sum{};
    for (auto&& map : maps) {
      Product new_product{};
      new_product.scale(product.scalar());
      auto temp_product =
          remove_tensor_from_product(product, L"S")->as<Product>();
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
      } else if (term->is<Tensor>()) {
        result->append(term);
      } else if (term->is<Constant>()) {
        result->append(term);
      }
    }
  }

  result->visit(reset_idx_tags);
  return result;
}

int count_cycles(const container::svector<int, 6>& vec1,
                 const container::svector<int, 6>& vec2) {
  assert(vec1.size() == vec2.size());
  int n_cycles = 0;
  auto dummy_idx = 99;
  container::svector<int, 6> v = vec1;
  container::svector<int, 6> v1 = vec2;
  for (auto it = v.begin(); it != v.end(); ++it) {
    if (*it != dummy_idx) {
      n_cycles++;
      auto idx = std::distance(v.begin(), it);
      auto it0 = it;
      auto it1 = std::find(v1.begin(), v1.end(), *it0);
      auto idx1 = std::distance(v1.begin(), it1);
      do {
        it0 = std::find(v.begin(), v.end(), v[idx1]);
        it1 = std::find(v1.begin(), v1.end(), *it0);
        idx1 = std::distance(v1.begin(), it1);
        *it0 = dummy_idx;
      } while (idx1 != idx);
    }
  }
  return n_cycles;
}

ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::vector<container::vector<Index>> ext_index_groups) {
  // NOT supported for Proto indices
  auto check_proto_index = [](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(expr->as<Tensor>().const_braket(), [](const Index& idx) {
        assert(!idx.has_proto_indices() &&
               "Proto index not supported in spintrace call.");
      });
    }
  };
  expression->visit(check_proto_index);

  // Symmetrize the expr and keep S operator
  auto symm_and_expand = [](const ExprPtr& expr) {
    auto temp = expr;
    if (has_A_label(temp)) temp = expr_symmetrize(temp);
    temp = expand_antisymm(temp);
    rapid_simplify(temp);
    return temp;
  };
  auto expr = symm_and_expand(expression);

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };
  expr->visit(reset_idx_tags);  // This call is REQUIRED
  expand(expr);                 // This call is REQUIRED
  rapid_simplify(expr);

  // Returns the number of cycles
  auto count_cycles = [](container::svector<Index>& v,
                         container::svector<Index>& v1) {
    assert(v.size() == v1.size());
    size_t n_cycles = 0;
    auto dummy_idx = Index(L"p_50");
    for (auto it = v.begin(); it != v.end(); ++it) {
      if (*it != dummy_idx) {
        // TODO: Throw exception if bra and ket indices don't match
        n_cycles++;
        auto idx = std::distance(v.begin(), it);
        auto it0 = it;
        auto it1 = std::find(v1.begin(), v1.end(), *it0);
        auto idx1 = std::distance(v1.begin(), it1);
        do {
          it0 = std::find(v.begin(), v.end(), v[idx1]);
          it1 = std::find(v1.begin(), v1.end(), *it0);
          idx1 = std::distance(v1.begin(), it1);
          *it0 = dummy_idx;
        } while (idx1 != idx);
      }
    }
    return n_cycles;
  };

  // Lambda for a product
  auto trace_product = [&ext_index_groups,
                        &count_cycles](const Product& product) {
    // TODO: Check symmetry of tensors

    // Remove S if present in a product
    Product temp_product{};
    temp_product.scale(product.scalar());
    if (product.factor(0)->as<Tensor>().label() == L"S") {
      for (auto&& term : product.factors()) {
        if (term->is<Tensor>() && term->as<Tensor>().label() != L"S")
          temp_product.append(term);
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

    auto product_kets = get_ket_indices(temp_product);
    auto product_bras = get_bra_indices(temp_product);

    // Substitute indices from external index list
    if ((*ext_index_groups.begin()).size() == 2)
      ranges::for_each(ext_index_groups,
                       [&product_bras, &product_kets](
                           const container::vector<Index>& idx_pair) {
                         assert(idx_pair.size() == 2);
                         if (idx_pair.size() == 2) {
                           auto it = idx_pair.begin();
                           auto first = *it;
                           it++;
                           auto second = *it;
                           std::replace(product_bras.begin(),
                                        product_bras.end(), first, second);
                           std::replace(product_kets.begin(),
                                        product_kets.end(), first, second);
                         }
                       });

    auto n_cycles = count_cycles(product_kets, product_bras);

    auto result = std::make_shared<Product>(product);
    result->scale(std::pow(2, n_cycles));
    return result;
  };

  if (expr->is<Constant>())
    return expr;
  else if (expr->is<Tensor>())
    return trace_product(
        (ex<Constant>(1.) * expr)->as<Product>());  // expand_all(expr);
  else if (expr->is<Product>())
    return trace_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      if (summand->is<Product>()) {
        result->append(trace_product(summand->as<Product>()));
      } else if (summand->is<Tensor>()) {
        result->append(
            trace_product((ex<Constant>(1.) * summand)->as<Product>()));
      } else  // summand->is<Constant>()
        result->append(summand);
    }
    return result;
  } else
    return nullptr;
}

container::vector<container::vector<Index>> external_indices(const ExprPtr& expr){
  // Generate external index list from Antisymmetrizer
  Tensor A{};
  for(const auto& prod : *expr){
    if(prod->is<Product>()){
      auto tensor = prod->as<Product>().factor(0)->as<Tensor>();
      if(tensor.label() == L"A"){
        A = tensor;
        break;
      }
    }
  }
  assert(A.bra_rank() != 0 && "Could not generate external index groups due to "
         "absence of Anti-symmetrizer (A) operator in expression.");
  assert(A.bra_rank() == A.ket_rank());
  container::vector<container::vector<Index>> ext_index_groups(A.rank());
  for(auto i = 0; i != A.rank(); ++i){
    ext_index_groups[i] = container::vector<Index>{A.ket()[i], A.bra()[i]};
  }
  return ext_index_groups;
}

ExprPtr closed_shell_CC_spintrace(const ExprPtr& expr){
  return closed_shell_spintrace(expr, external_indices(expr));
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

Tensor swap_spin(const Tensor& t){

  auto is_null_qns = [](const Index& i){
    return i.space().qns() == IndexSpace::nullqns;
  };

  // Return tensor if there are no spin labels
  if(std::all_of(t.const_braket().begin(),
                 t.const_braket().end(), is_null_qns)){
    return t;
  }

  // Return new index where the spin-label is flipped
  auto spin_flipped_idx = [] (const Index& idx) {
    assert(idx.space().qns() != IndexSpace::nullqns);

    auto idx_type = idx.space().type();
    auto qns = idx.space().qns() == IndexSpace::alpha ?
              IndexSpace::beta: IndexSpace::alpha;
    auto label_int =
        idx.label().substr(idx.label().find(L'_') + 1);
    std::wstring label_int_ws(label_int.begin(), label_int.end());

    return Index::make_label_index(IndexSpace(idx_type, qns), label_int_ws);
  };

  std::vector<Index> bra(t.rank()), ket(t.rank());

  for(auto i = 0; i < t.rank(); ++i){
    bra.at(i) = spin_flipped_idx(t.bra().at(i));
    ket.at(i) = spin_flipped_idx(t.ket().at(i));
  }

  return {t.label(), bra, ket,
          t.symmetry(), t.braket_symmetry(),
          t.particle_symmetry()};
}

ExprPtr swap_spin(const ExprPtr& expr){
  if(expr->is<Constant>())
    return expr;

  auto swap_tensor = [] (const Tensor& t) {
    return ex<Tensor>(swap_spin(t));
  };

  auto swap_product = [&swap_tensor] (const Product& p){
    Product result{};
    result.scale(p.scalar());
    for(auto& t : p){
      assert(t->is<Tensor>());
      result.append(swap_tensor(t->as<Tensor>()));
    }
    return ex<Product>(result);
  };

  if(expr->is<Tensor>())
    return swap_tensor(expr->as<Tensor>());
  else if (expr->is<Product>())
    return swap_product(expr->as<Product>());
  else if (expr->is<Sum>()){
    Sum result;
    for(auto& term : *expr){
      result.append(swap_spin(term));
    }
    return ex<Sum>(result);
  } else
    return nullptr;
}

ExprPtr merge_operators(const Tensor& O1, const Tensor& O2){
  assert(O1.label() == O2.label());
  assert(O1.symmetry() == O2.symmetry());
  auto bra = ranges::views::concat(O1.bra(), O2.bra());
  auto ket = ranges::views::concat(O1.ket(), O2.ket());
  return ex<Tensor>(Tensor(O1.label(), bra, ket, O1.symmetry()));
}

std::vector<ExprPtr> open_shell_P_op_vector(const Tensor &A){
  assert(A.label() == L"A");

  // N+1 spin-cases for corresponding residual
  std::vector<ExprPtr> result_vector(A.bra_rank() + 1);

  // List of indices
  const auto rank = A.bra_rank();
  std::vector<int> idx(rank);
  auto bra = A.bra();
  auto ket = A.ket();
  std::iota(idx.begin(), idx.end(), 0);

  // Anti-symmetrizer is preserved for all identical spin cases,
  // So return a constant
  result_vector.at(0) = ex<Constant>(1); // all alpha
  result_vector.at(rank) = ex<Constant>(1);  // all beta

  // This loop generates all the remaining spin cases
  for (auto i = 1; i < rank; ++i) {
    std::vector<int> alpha_spin(idx.begin(), idx.end() - i);
    std::vector<int> beta_spin(idx.end() - i, idx.end());

    std::vector<Tensor> P_bra_list, P_ket_list;
    for (auto& j : alpha_spin) {
      for (auto& k : beta_spin) {
        if (!alpha_spin.empty() && !beta_spin.empty()) {
          P_bra_list.emplace_back(Tensor(L"P", {bra.at(j), bra.at(k)}, {}));
          P_ket_list.emplace_back(Tensor(L"P", {}, {ket.at(j), ket.at(k)}));
        }
      }
    }

    // The P4 terms
    if (alpha_spin.size() > 1 && beta_spin.size() > 1) {
      for(auto a = 0; a != alpha_spin.size() - 1; ++a){
        auto i1 = alpha_spin[a];
        for(auto b = a+1; b != alpha_spin.size(); ++b){
          auto i2 = alpha_spin[b];
        for(auto c = 0; c != beta_spin.size() - 1; ++c){
          auto i3 = beta_spin[c];
          for(auto d = c+1; d != beta_spin.size(); ++d) {
            auto i4 = beta_spin[d];
            P_bra_list.emplace_back(
                Tensor(L"P", {bra.at(i1), bra.at(i3),
                              bra.at(i2), bra.at(i4)}, {}));
            P_ket_list.emplace_back(
                Tensor(L"P", {}, {ket.at(i1), ket.at(i3),
                        ket.at(i2), ket.at(i4)}));
          }
        }
        }
      }
    }

    Sum bra_permutations{};
    bra_permutations.append(ex<Constant>(1.));
    Sum ket_permutations{};
    ket_permutations.append(ex<Constant>(1.));

    for(auto& p : P_bra_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      bra_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    for(auto& p : P_ket_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      ket_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    ExprPtr spin_case_result = ex<Sum>(bra_permutations) * ex<Sum>(ket_permutations);
    simplify(spin_case_result);

    // Merge P operators if it encounters alpha_spin product of operators
    for(auto& term : *spin_case_result){
      if(term->is<Product>()) {
        auto P = term->as<Product>();
        if(P.factors().size() == 2){
          auto P1 = P.factor(0)->as<Tensor>();
          auto P2 = P.factor(1)->as<Tensor>();
          term = merge_operators(P1, P2);
        }
      }
    }
    result_vector.at(i) = spin_case_result;
  }
  return result_vector;
}

ExprPtr open_shell_P_operator(const int r, const int n){
  std::vector<int> i_alpha(r), i_beta(n-r);
  assert(!i_alpha.empty() && !i_beta.empty());
  std::iota(i_alpha.begin(), i_alpha.end(), 0);
  std::iota(i_beta.begin(), i_beta.end(), r);
  assert( r < n && n < 7);
  assert( 2*r != n && "Same number of alpha and beta spins not supported at the moment. "
         "Please completely expand A-opertor for this case.");

  // Index list of act_occ, act_uocc
  std::vector<Index> bra{Index(L"i_1"), Index(L"i_2"), Index(L"i_3"),
                         Index(L"i_4"), Index(L"i_5"), Index(L"i_6")};
  std::vector<Index> ket{Index(L"a_1"), Index(L"a_2"), Index(L"a_3"),
                         Index(L"a_4"), Index(L"a_5"), Index(L"a_6")};

  // Dummy index vectors of size 2
  std::vector<Index> bra_dummy{Index(L"i_98"), Index(L"i_99")};
  std::vector<Index> ket_dummy{Index(L"a_98"), Index(L"a_99")};

  // Generate permutation pairs
  std::vector<std::vector<Index>> bras{}, kets{};
  for(auto& i : i_alpha){
    for(auto& j : i_beta){
      bras.insert(bras.end(),std::vector{bra[i], bra[j]});
      kets.insert(kets.end(),std::vector{ket[i], ket[j]});
    }
  }

  std::vector<std::vector<Index>> brakets{};
  for(auto& bra_pair : bras){
    for(auto& ket_pair : kets){
      brakets.insert(brakets.end(),
                       std::vector{bra_pair[0],bra_pair[1],
                                      ket_pair[0],ket_pair[1]});
    }
  }

  Sum result{};
  result.append(ex<Constant>(1));

  for(auto& bra_pair : bras){
    result.append(ex<Constant>(-1) *
                  ex<Tensor>(Tensor(L"P",bra_pair,ket_dummy)));
  }

  for(auto& ket_pair : kets){
    result.append(ex<Constant>(-1) *
                  ex<Tensor>(Tensor(L"P",bra_dummy,ket_pair)));
  }

  for(auto& bk : brakets){
    result.append(ex<Tensor>(Tensor(L"P",{bk[0], bk[1]},{bk[2], bk[3]})));
  }

  // std::wcout << to_latex(result) << "\n";
  return ex<Sum>(result);
}

std::vector<ExprPtr> open_shell_spintrace(const ExprPtr& expr,
                                          const std::vector<std::vector<Index>> ext_index_groups,
                                          const int single_spin_case){

  if(expr->is<Constant>()){
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

  using IndexGroup = container::vector<Index>;
  std::vector<IndexGroup> int_index_groups;
  for (auto&& i : int_idxlist) {
    int_index_groups.emplace_back(IndexGroup(1, i));
  }

  assert(grand_idxlist.size() == int_idxlist.size() + ext_idxlist.size());

  // Add spin label to index
  auto add_spin_label = [] (const Index& idx, const long int& spin_bit){
    auto idx_n = idx.label().substr(idx.label().find(L'_') + 1);
    std::wstring idx_n_ws(idx_n.begin(), idx_n.end());

    auto idx_type = IndexSpace::instance(idx.label()).type();
    auto space = spin_bit == 0 ?
                 IndexSpace::instance(idx_type, IndexSpace::alpha):
                 IndexSpace::instance(idx_type, IndexSpace::beta);

    return Index::make_label_index(space, idx_n_ws);
  };

  // Generate index replacement maps
  auto spin_cases = [&add_spin_label] (const std::vector<IndexGroup>& idx_group) {
    int ncases = std::pow(2, idx_group.size());
    std::vector<std::map<Index, Index>> all_replacements(ncases);

    for (uint64_t i = 0; i != ncases; ++i) {
      std::map<Index, Index> idx_rep;
      for(size_t idxg = 0; idxg != idx_group.size(); ++idxg){
        auto spin_bit = (i << (64 - idxg - 1)) >> 63;
        assert((spin_bit == 0) || (spin_bit == 1));
        for(auto& idx : idx_group[idxg]){
          auto spin_idx = add_spin_label(idx, spin_bit);
          idx_rep.insert({idx, spin_idx});
        }
      }
      all_replacements[i] = idx_rep;
    }
    return all_replacements;
  };

  // External index replacement maps
  auto ext_spin_cases = [&add_spin_label] (const std::vector<IndexGroup>& idx_group){
    std::vector<std::map<Index, Index>> all_replacements;

    std::vector<int> spins(idx_group.size(), 0);
    for(size_t i = 0; i <= idx_group.size(); ++i){
      std::vector<int> spins(idx_group.size(), 0);
      std::fill(spins.end()-i,spins.end(), 1);

      std::map<Index, Index> idx_rep;
      for(size_t j = 0; j != idx_group.size(); ++j){
        for(auto &idx : idx_group[j]) {
          auto spin_idx = add_spin_label(idx, spins[j]);
          idx_rep.insert({idx, spin_idx});
        }
      }
      all_replacements.push_back(idx_rep);
    }
    return all_replacements;
  };

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };

  // Internal and external index replacements are independent
  auto i_rep = spin_cases(int_index_groups);
  auto e_rep = ext_spin_cases(ext_index_groups);

  // For a single spin case, keep only the relevant spin case
  // PS: all alpha indexing start at 0
  if(single_spin_case){
    auto external_replacement_map = e_rep.at(single_spin_case);
    e_rep.clear();
    e_rep.push_back(external_replacement_map);
  }

  // Expand 'A' operator and 'antisymm' tensors
  auto expanded_expr = expand_A_operator(expr);
  expand(expanded_expr);
  rapid_simplify(expanded_expr);
  expanded_expr->visit(reset_idx_tags);

  std::vector<ExprPtr> result{};

  // return true if a product is spin-symmetric
  auto spin_symm_product = [] (const Product& product) {

    std::vector<Index> cBra, cKet; // concat Bra and concat Ket
    for(auto& term : product){
      if(term->is<Tensor>()){
        auto tnsr = term->as<Tensor>();
        cBra.insert(cBra.end(), tnsr.bra().begin(), tnsr.bra().end());
        cKet.insert(cKet.end(), tnsr.ket().begin(), tnsr.ket().end());
      }
    }
    assert(cKet.size() == cBra.size());

    auto i_ket = cKet.begin();
    for(auto& b : cBra){
      if (b.space().qns() != i_ket->space().qns())
        return false;
      ++i_ket;
    }
    return true;
  };

  //
  // SPIN-TRACING algorithm begins here
  //

  // Loop over external index replacement maps
  for(auto& e : e_rep){

    // Add spin labels to external indices
    auto spin_expr = append_spin(expanded_expr, e);
    spin_expr->visit(reset_idx_tags);
    Sum e_result{};

    // Loop over internal index replacement maps
    for(auto& i : i_rep){

      // Add spin labels to internal indices, expand antisymmetric tensors
      auto spin_expr_i = append_spin(spin_expr, i);
      spin_expr_i = expand_antisymm(spin_expr_i, true);
      expand(spin_expr_i);
      spin_expr_i->visit(reset_idx_tags);
      Sum i_result{};

      if(spin_expr_i->is<Tensor>()){
        e_result.append(spin_expr_i);
      } else if(spin_expr_i->is<Product>()){
        if (spin_symm_product(spin_expr_i->as<Product>()))
          e_result.append(spin_expr_i);
      } else if(spin_expr_i->is<Sum>()){
        for(auto& pr : *spin_expr_i){
          if (pr->is<Product>()){
            if (spin_symm_product(pr->as<Product>()))
              i_result.append(pr);
          } else if (pr->is<Tensor>()){
            if (is_tensor_spin_symm(pr->as<Tensor>()))
              i_result.append(pr);
          } else if (pr->is<Constant>()){
            i_result.append(pr);
          } else
            throw("Unknown ExprPtr type.");
        }
        e_result.append(std::make_shared<Sum>(i_result));
      }

    } // loop over internal indices
    result.push_back(std::make_shared<Sum>(e_result));
  } // loop over external indices

  if(single_spin_case){
    assert(result.size() == 1 && "Spin-specific case must return one expression.");
  }

  // Canonicalize and simplify all expressions
  for(auto& expression: result){
    expression->visit(reset_idx_tags);
    canonicalize(expression);
    rapid_simplify(expression);
  }
  return result;
}

std::vector<ExprPtr> open_shell_CC_spintrace(const ExprPtr& expr){
  return open_shell_spintrace(expr, external_indices(expr));
}

ExprPtr spintrace(
    ExprPtr expression,
    container::vector<container::vector<Index>> ext_index_groups) {
  // SPIN TRACE DOES NOT SUPPORT PROTO INDICES YET.
  auto check_proto_index = [](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(expr->as<Tensor>().const_braket(), [](const Index& idx) {
        assert(!idx.has_proto_indices() &&
               "Proto index not supported in spintrace function.");
      });
    }
  };
  expression->visit(check_proto_index);

  if (expression->is<Constant>()) {
    return expression;
  }

  auto spin_trace_tensor = [](const Tensor& tensor) {
    if (can_expand(tensor)) {
      return expand_antisymm(tensor);
    } else
      return ex<Constant>(0);
  };

  auto spin_trace_product = [&spin_trace_tensor](const Product& product) {
    Product spin_product{};
    spin_product.scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        if (can_expand(term->as<Tensor>())) {
          spin_product.append(1, spin_trace_tensor(term->as<Tensor>()));
        } else
          break;
      }
    }
    if (product.size() != spin_product.size()) {
      spin_product.scale(0);
    }
    ExprPtr result = std::make_shared<Product>(spin_product);
    expand(result);
    rapid_simplify(result);
    return result;
  };

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };

  // Most important lambda of this function
  auto trace_product = [&ext_index_groups, &spin_trace_tensor,
                        &spin_trace_product](const Product& expression) {
    auto result = std::make_shared<Sum>();
    ExprPtr expr = std::make_shared<Product>(expression);

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

    // EFV: generate the grand list of index groups by concatenating list of
    // external index EFV: groups with the groups of internal indices (each
    // internal index = 1 group)
    using IndexGroup = container::vector<Index>;
    container::vector<IndexGroup> index_groups;

    for (auto&& i : int_idxlist) {
      index_groups.emplace_back(IndexGroup(1, i));
    }

    index_groups.insert(index_groups.end(), ext_index_groups.begin(),
                        ext_index_groups.end());

    // EFV: for each spincase (loop over integer from 0 to 2^n-1, n=#of index
    // groups)

    const uint64_t nspincases = std::pow(2, index_groups.size());

    for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
         ++spincase_bitstr) {
      // EFV:  assign spin to each index group => make a replacement list
      std::map<Index, Index> index_replacements;

      uint64_t index_group_count = 0;
      for (auto&& index_group : index_groups) {
        auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
        assert((spin_bit == 0) || (spin_bit == 1));

        for (auto&& index : index_group) {
          IndexSpace space;
          space = spin_bit == 0
                      ? IndexSpace::instance(
                            IndexSpace::instance(index.label()).type(),
                            IndexSpace::alpha)
                      : IndexSpace::instance(
                            IndexSpace::instance(index.label()).type(),
                            IndexSpace::beta);

          // TODO: Check if valid for index with no subscripts
          auto subscript_label =
              index.label().substr(index.label().find(L'_') + 1);
          std::wstring subscript_label_ws(subscript_label.begin(),
                                          subscript_label.end());
          Index spin_index = Index::make_label_index(space, subscript_label_ws);
          index_replacements.insert({index, spin_index});
        }
        ++index_group_count;
      }

      auto spin_expr = append_spin(expr, index_replacements);
      // This call to rapid_simplify is required for Tensor case
      rapid_simplify(spin_expr);

      if (spin_expr->is<Tensor>()) {
        auto temp = spin_trace_tensor(spin_expr->as<Tensor>());
        auto spin_removed = remove_spin(temp);
        result->append(spin_removed);
      } else if (spin_expr->is<Product>()) {
        auto temp = spin_trace_product(spin_expr->as<Product>());
        if (temp->size() != 0) {
          result->append(remove_spin(temp));
        }
      } else if (spin_expr->is<Sum>()) {
        for (auto&& summand : *spin_expr) {
          Sum temp{};
          if (summand->is<Tensor>())
            temp.append(spin_trace_tensor(summand->as<Tensor>()));
          else if (summand->is<Product>())
            temp.append(spin_trace_product(summand->as<Product>()));
          else {
            temp.append(summand);
          }
          ExprPtr SumPtr = std::make_shared<Sum>(temp);
          auto spin_removed = remove_spin(SumPtr);
          result->append(spin_removed);
        }
      } else {
        result->append(expr);
      }
    }  // Permutation FOR loop
    return result;
  };

  // Expand antisymmetrizer operator (A) if present in the expression
  if (has_tensor_label(expression, L"A"))
    expression = expand_A_operator(expression);

  if (expression->is<Tensor>()) expression = ex<Constant>(1) * expression;

  if (expression->is<Product>()) {
    return trace_product(expression->as<Product>());
  } else if ((expression->is<Sum>())) {
    auto result = std::make_shared<Sum>();
    for (auto&& term : *expression) {
      if (term->is<Product>())
        result->append(trace_product(term->as<Product>()));
      else if (term->is<Tensor>()) {
        auto term_as_product = ex<Constant>(1) * term;
        result->append(trace_product(term_as_product->as<Product>()));
      } else
        result->append(term);
    }
    result->visit(reset_idx_tags);
    return result;
  } else
    return nullptr;
}  // ExprPtr spintrace

ExprPtr factorize_S_operator(
    const ExprPtr& expression,
    const std::initializer_list<IndexList> ext_index_groups,
    const bool fast_method) {
  // Canonicalize the expression
  ExprPtr expr = expression;
  // canonicalize(expr);

  // If expression has S operator, do nothing and exit
  if (has_operator_label(expr, L"S")) return expr;

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
    S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  }

  // For any order CC residual equation:
  // Generate a list of permutation indices
  // Erase the canonical entry
  auto replacement_maps = S_replacement_maps(S);
  replacement_maps.erase(replacement_maps.begin());

  // Lambda function for index replacement in tensor
  auto transform_tensor = [](const Tensor& tensor,
                             const std::map<Index, Index>& replacement_map) {
    auto result = std::make_shared<Tensor>(tensor);
    result->transform_indices(replacement_map);
    ranges::for_each(result->const_braket(),
                     [](const Index& idx) { idx.reset_tag(); });
    return result;
  };

  Sum result_sum{};
  ///////////////////////////////////////////////
  ///            Fast approach                ///
  ///////////////////////////////////////////////
  // This method is stores hash values of terms in a container
  // for faster run times

  if (fast_method) {
    // summands_hash_list sorted container of hash values of canonicalized
    // summands summands_hash_map unsorted map of (hash_val, summand) pairs
    // container::set<size_t> summands_hash_list;
    container::vector<size_t> summands_hash_list;
    std::unordered_map<size_t, ExprPtr> summands_hash_map;
    for (auto it = expr->begin(); it != expr->end(); ++it) {
      (*it)->canonicalize();
      auto hash = (*it)->hash_value();
      summands_hash_list.push_back(hash);
      summands_hash_map.insert({hash, *it});
    }
    assert(summands_hash_list.size() == expr->size());
    assert(summands_hash_map.size() == expr->size());

    // Symmetrize every summand, assign its hash value to hash1
    // Check if hash1 exist in summands_hash_list
    // if(hash1 present in summands_hash_list) remove hash0, hash1
    // else continue
    int n_symm_terms = 0;
    auto symm_factor = std::tgamma(S.bra_rank() + 1);
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
          ex<Constant>(1.0 / symm_factor) * ex<Tensor>(S) * new_product;

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
                  transform_tensor(t->as<Tensor>(), replacement_map));
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
      auto n_hash_found = ranges::count_if(hash1_list, hash1_found);

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
            if (summand->hash_value() == hash1) summand = ex<Constant>(0.0);
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

    int n_symm_terms = 0;

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
                  transform_tensor(t->as<Tensor>(), replacement_map));
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
        int n_matches = 0;
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
  expand(result);
  canonicalize(result);
  rapid_simplify(result);
  return result;
}

}  // namespace sequant
