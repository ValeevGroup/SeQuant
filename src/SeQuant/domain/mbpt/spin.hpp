//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include "SeQuant/core/tensor.hpp"

namespace sequant {

/// @brief Adds spins to indices in an expression using a map
/// @param expr an expression pointer
/// @param index_replacements a map of pairs containing the index and its
/// replacement
/// @return expr the expression with substituted indices
ExprPtr append_spin(ExprPtr& expr, std::map<Index, Index>& index_replacements) {
  auto add_spin_to_tensor = [&](const Tensor& tensor) {
    auto spin_tensor = std::make_shared<Tensor>(tensor);
    auto pass_mutated = spin_tensor->transform_indices(index_replacements);
    return spin_tensor;
  };

  auto add_spin_to_product = [&](const Product& product) {
    auto spin_product = std::make_shared<Product>();
    spin_product->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        spin_product->append(1, add_spin_to_tensor(term->as<Tensor>()));
      } else
        abort();
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
    return expr;
}

/// @brief Removes spin label from all indices in an expression
/// @param expr an expression pointer with spin labels
/// @return expr an expression pointer without spin labels
ExprPtr remove_spin(ExprPtr& expr) {
  auto remove_spin_from_tensor = [&](const Tensor& tensor) {
    std::vector<Index> bra;  // An unordered list is required
    std::vector<Index> ket;  // An unordered list is required
    {
      for (auto&& idx : tensor.bra()) bra.emplace_back(idx);
      for (auto&& idx : tensor.ket()) ket.emplace_back(idx);
      auto braket_list = ranges::views::concat(bra, ket);

      for (auto&& idx : braket_list) {
        auto space = IndexSpace::instance(
            IndexSpace::instance(idx.label()).type(), IndexSpace::nullqns);
        auto subscript_label = idx.label().substr(idx.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(),
                                        subscript_label.end());
        idx = Index::make_label_index(space, subscript_label_ws);
      }
    }
    auto sft = Tensor(tensor.label(), bra, ket, tensor.symmetry(),
                      tensor.braket_symmetry());
    auto sf_tensor = std::make_shared<Tensor>(sft);
    return sf_tensor;
  };

  auto remove_spin_from_product = [&](const Product& product) {
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

/// @brief Checks the spin symmetry of pairs of indices corresponding to a
/// particle in tensor notation
/// @param tensor a tensor with indices containing spin labels
/// @return true if spin symmetry matches for all pairs of indices
inline bool is_tensor_spin_symm(const Tensor& tensor) {
  bool result = false;
  assert(tensor.bra().size() == tensor.ket().size());
  // For each index check if QNS match.
  // return FALSE if one of the pairs does not match.
  auto iter_ket = tensor.ket().begin();
  for (auto&& index : tensor.bra()) {
    if (IndexSpace::instance(index.label()).qns() ==
        IndexSpace::instance(iter_ket->label()).qns()) {
      result = true;
    } else {
      result = false;
      return result;
    }
    iter_ket++;
  }
  return result;
}

/// @brief Check if the number of alpha spins in bra and ket are equal
/// beta spins will match if total number of indices is the same
/// @param any tensor
/// @return true if number of alpha spins match in bra and ket
inline bool can_expand(const Tensor& tensor) {
  assert(tensor.bra_rank() == tensor.ket_rank());
  auto result = false;
  auto alpha_in_bra = 0;
  auto alpha_in_ket = 0;
  ranges::for_each(tensor.bra(), [&](Index i) {
    if (IndexSpace::instance(i.label()).qns() == IndexSpace::alpha)
      ++alpha_in_bra;
  });
  ranges::for_each(tensor.ket(), [&](Index i) {
    if (IndexSpace::instance(i.label()).qns() == IndexSpace::alpha)
      ++alpha_in_ket;
  });
  if (alpha_in_bra == alpha_in_ket) result = true;
  return result;
}

/// @brief expand an antisymmetric tensor
/// @param tensor a tensor from a product
/// @return an expression pointer containing the sum of expanded terms if
/// antisymmetric or
/// @return an expression pointer containing the tensor otherwise
ExprPtr expand_antisymm(const Tensor& tensor) {
  assert(tensor.bra().size() == tensor.ket().size());

  // Generate a sum of asymmetric tensors if the input tensor is antisymmetric
  // AND more than one body otherwise, return the tensor
  if ((tensor.symmetry() == Symmetry::antisymm) && (tensor.bra().size() > 1)) {
    container::set<Index> bra_list;
    for (auto&& bra_idx : tensor.bra()) bra_list.insert(bra_idx);
    const auto const_bra_list = bra_list;

    container::set<Index> ket_list;
    for (auto&& ket_idx : tensor.ket()) ket_list.insert(ket_idx);

    Sum expr_sum{};
    auto p_count = 0;
    do {
      // Generate tensor with new labels
      auto new_tensor =
          Tensor(tensor.label(), bra_list, ket_list, Symmetry::nonsymm);

      if (is_tensor_spin_symm(new_tensor)) {
        auto new_tensor_ptr = ex<Tensor>(new_tensor);
        Product new_tensor_product{};
        IndexSwapper::thread_instance().reset();
        using std::begin;
        using std::end;
        auto bra_list2 = bra_list;
        bubble_sort(begin(bra_list2), end(bra_list2), std::less<Index>{});
        bubble_sort(begin(ket_list), end(ket_list), std::less<Index>{});
        bool even = IndexSwapper::thread_instance().even_num_of_swaps();
        new_tensor_product.append((even ? 1 : -1), new_tensor_ptr);
        auto new_tensor_product_ptr = ex<Product>(new_tensor_product);
        // std::wcout << __LINE__ << "L "<< (even ? 1 : -1) << to_latex(new_tensor_product_ptr) << "\n";
        expr_sum.append(new_tensor_product_ptr);
      }
      p_count++;
    } while (std::next_permutation(bra_list.begin(), bra_list.end()));

    ExprPtr result = std::make_shared<Sum>(expr_sum);
    // std::wcout << __LINE__ << "L "<< to_latex(result) << "\n\n";
    return result;
  } else {
    ExprPtr result = std::make_shared<Tensor>(tensor);
    return result;
  }
}

/// @brief expands all antisymmetric tensors in a product
/// @param expr an expression pointer to expand
/// @return an expression pointer with expanded tensors
inline ExprPtr expand_antisymm(const ExprPtr& expr) {
  Sum expanded_sum{};
  for (auto&& summand : *expr) {
    const auto scalar_factor = summand->as<Product>().scalar().real();
    Product expanded_product{};
    for (auto&& expr_product : *summand) {
      auto tensor = expr_product->as<sequant::Tensor>();
      auto expanded_ptr = expand_antisymm(tensor);
      expanded_product.append(1, expanded_ptr);
    }
    expanded_product.scale(scalar_factor);
    ExprPtr expanded_product_ptr = std::make_shared<Product>(expanded_product);
    expanded_sum.append(expanded_product_ptr);
  }
  ExprPtr result = std::make_shared<Sum>(expanded_sum);
  return result;
}

/// @brief Check if the Antisymmetrizer operator is present in given expression
/// @param expr input expression
/// @return true if this function finds an A operator
bool check_A_operator(const ExprPtr& expr){

  auto check_tensor = [&] (const Tensor& tensor) {
    if(tensor.label() == L"A")
      return true;
    else
      return false;
  };

  auto check_product = [&] (const Product& product) {
    bool result = false;
    for(auto&& term: product){
      if(term->is<Tensor>()){
        if(check_tensor(term->as<Tensor>())){
          result = true;
          break;
        }
      }
    }
    return result;
  };

  if(expr->is<Constant>()){
    return false;
  } else if (expr->is<Tensor>()){
    return check_tensor(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return check_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    bool result = false;
    for(auto&& term : *expr){
      if(term->is<Product>())
        result = check_product(term->as<Product>());
      if(result)
        return result;
    }
  } else
    throw("Unknown arg type for check_Antisymm_operator.");
}

/// @brief Generates a vector of replacement maps
/// @param A The antisymmetrizer with replacement indices
/// @return
std::vector<std::map<Index, Index>> A_replacement_map(const Tensor& A){

  // Check A, get bra list, get ket list
  container::set<Index> A_bra;
  container::set<Index> A_ket;
  container::svector<Index> A_braket;
  if(A.label() == L"A"){
    for(auto& idx: A.bra()) A_bra.insert(idx);
    for(auto& idx: A.ket()) A_ket.insert(idx);
    for(auto& idx: A.const_braket()) A_braket.push_back(idx);
  } else
    throw("A_replacement_map needs the antisymmetrizer tensor");

  std::vector<std::map<Index, Index>> result;
  do{
    do{
      auto A_bra_copy = A_bra;
      auto A_ket_copy = A_ket;
      std::map<Index, Index> replacement_map;
      auto A_braket_ptr = A_braket.begin();
      for(auto&& idx: A_bra){
        replacement_map.emplace(std::make_pair(*A_braket_ptr, idx));
        A_braket_ptr++;
      }
      for(auto&& idx: A_ket){
        replacement_map.emplace(std::make_pair(*A_braket_ptr, idx));
        A_braket_ptr++;
      }
      result.push_back(replacement_map);
    } while (std::next_permutation(A_bra.begin(), A_bra.end()));
  } while (std::next_permutation(A_ket.begin(), A_ket.end()));

  return result;
}


/// @brief Expand expression with Antisymmetrization operator
/// @param A product term
/// @return expression pointer with A operator applied
ExprPtr expand_A_operator(const Product& product){

  bool has_A_operator = false;

  // Check A and build replacement map
  std::vector<std::map<Index, Index>> map_list;
  for(auto& term : product){
    if(term->is<Tensor>()){
      auto A = term->as<Tensor>();
      if(A.label() == L"A"){
        has_A_operator = true;
        map_list = A_replacement_map(A);
        break;
      }
    }
  }

  if(!has_A_operator) return std::make_shared<Product>(product);

  // substitute and expand
  Product large_product{};
  auto product_of_expanded_terms = std::make_shared<Product>();
  product_of_expanded_terms->scale(product.scalar());
  for(auto& term : product){
  Sum new_sum{};
    if((term->is<Tensor>()) && (term->as<Tensor>().label() != L"A")){
      auto tensor = term->as<Tensor>();
      for(auto&& map : map_list){
        auto new_tensor = tensor;
        new_tensor.transform_indices(map);
        container::svector<Index> new_tensor_braket;
        ranges::for_each(new_tensor.const_braket(), [&] (const Index& idx){new_tensor_braket.push_back(idx);});
        IndexSwapper::thread_instance().reset();
        bubble_sort(std::begin(new_tensor_braket), std::end(new_tensor_braket), std::less<Index>{});
        bool even = IndexSwapper::thread_instance().even_num_of_swaps();
        auto new_tensor_ptr = std::make_shared<Tensor>(new_tensor);
        Product new_product{};
        new_product.append((even ? 1 : -1), new_tensor_ptr);
//         std::wcout << to_latex(new_product) << "\n";
        auto new_product_ptr = std::make_shared<Product>(new_product);
        new_sum.append(new_product_ptr);
      }
//      std::cout << "\nSum: ";
//      std::wcout << to_latex(new_sum) << "\n";
      auto new_sum_ptr = std::make_shared<Sum>(new_sum);
      product_of_expanded_terms->append(1, new_sum_ptr);
    }
  }
  ExprPtr result = product_of_expanded_terms;
  expand(result);
  rapid_simplify(result);
  return result;
}

/// @brief expands all A operators
/// @param expr A product or sum
/// @return expression pointer with A removed.
ExprPtr expand_A_operator(const ExprPtr& expr){
  if(expr->is<Constant>() || expr->is<Tensor>())
    throw("Unknown arg for expand_A_operator");

    if (expr->is<Product>())
      return expand_A_operator(expr->as<Product>());
    else if(expr->is<Sum>()){
      auto result = std::make_shared<Sum>();
      for(auto&& summand: *expr){
        if(summand->is<Product>())
          result->append(expand_A_operator(summand->as<Product>()));
        else
          result->append(summand);
      }
      return result;
    } else
      throw("Unknown arg Type for expand_A_operator.");
}

/// @brief Spin traces a given expression pointer
/// @detailed Given an expression, this function extracts all indices and adds a
/// spin attribute to all the indices in the expression. A map is generated with
/// all possible spin permutations and substituted in the expression. Only the
/// non-zero terms are kept, expanded, the spin labels removed and a sum of all
/// non-zero expressions is returned.
/// @param expr input expression
/// @param ext_index_groups groups of external indices
/// @return the expression with spin integrated out
ExprPtr spintrace(ExprPtr expression,
                  std::initializer_list<IndexList> ext_index_groups = {{}}) {
  const auto tstart = std::chrono::high_resolution_clock::now();

  // TODO: Fix this. Products have to be identified
  // SPIN TRACE DOES NOT SUPPORT PROTO INDICES YET.
  auto check_proto_index = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(
          expr->as<Tensor>().const_braket(),
          [&](const Index& idx) { assert(!idx.has_proto_indices()); });
    }
  };
  expression->visit(check_proto_index);

  if (expression->is<Constant>()) return expression;

  if (expression->is<Tensor>()){

    return expand_antisymm(expression->as<Tensor>());
  }


  auto spin_trace_tensor = [&](const Tensor& tensor) {
    auto spin_tensor = std::make_shared<Tensor>(tensor);
    if (can_expand(tensor)) {
      auto result = expand_antisymm(spin_tensor->as<Tensor>());
      return result;
    } else
      abort();
  };

  auto spin_trace_product = [&](const Product& product) {
    Product spin_product{};
    spin_product.scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        if (can_expand(term->as<Tensor>())) {
          spin_product.append(1, spin_trace_tensor(term->as<Tensor>()));
        }
      } else
        abort();
    }
    if (product.size() != spin_product.size()) spin_product.scale(0);
    ExprPtr result = std::make_shared<Product>(spin_product);
    expand(result);
    rapid_simplify(result);
    return result;
  };

  auto trace_product = [&](const Product& expression) {
    auto result = std::make_shared<Sum>();
    ExprPtr expr = std::make_shared<Product>(expression);

    container::set<Index, Index::LabelCompare> grand_idxlist;
    auto collect_indices = [&](const ExprPtr& expr) {
      if (expr->is<Tensor>()) {
        ranges::for_each(expr->as<Tensor>().const_braket(),
                         [&](const Index& idx) { grand_idxlist.insert(idx); });
      }
    };
    expr->visit(collect_indices);

    container::set<Index> ext_idxlist;
    for (auto&& idxgrp : ext_index_groups) {
      for (auto&& idx : idxgrp) {
        auto result = ext_idxlist.insert(idx);
        assert(result.second);
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
    using IndexGroup = container::svector<Index>;
    container::svector<IndexGroup> index_groups;

    for (auto&& i : int_idxlist) {
      index_groups.emplace_back(IndexGroup(1, i));
    }

    index_groups.insert(index_groups.end(), ext_index_groups.begin(),
                        ext_index_groups.end());

    // EFV: for each spincase (loop over integer from 0 to 2^n-1, n=#of index
    // groups)

    const uint64_t nspincases = std::pow(2, index_groups.size() - 1);

    for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
         ++spincase_bitstr) {
      // EFV:  assign spin to each index group => make a replacement list
      std::map<Index, Index> index_replacements;

      int64_t index_group_count = 0;
      for (auto&& index_group : index_groups) {
        auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
        assert((spin_bit == 0) || (spin_bit == 1));

        for (auto&& index : index_group) {
          IndexSpace space;
          if (spin_bit == 0) {
            space = IndexSpace::instance(
                IndexSpace::instance(index.label()).type(), IndexSpace::alpha);
          } else {
            space = IndexSpace::instance(
                IndexSpace::instance(index.label()).type(), IndexSpace::beta);
          }

          // TODO: Check if valid for index with no subscripts
          auto subscript_label =
              index.label().substr(index.label().find(L'_') + 1);
          std::wstring subscript_label_ws(subscript_label.begin(),
                                          subscript_label.end());

          Index spin_index = Index::make_label_index(space, subscript_label_ws);

          index_replacements.emplace(std::make_pair(index, spin_index));
        }
        ++index_group_count;
      }

      auto all_terms = std::make_shared<Sum>();
      auto spin_expr = append_spin(expr, index_replacements);
      rapid_simplify(spin_expr);

      if (spin_expr->is<Tensor>()) {
        auto temp = spin_trace_tensor(spin_expr->as<Tensor>());
        auto spin_removed = remove_spin(temp);
        result->append(spin_removed);
      } else if (spin_expr->is<Product>()) {
        auto temp = spin_trace_product(spin_expr->as<Product>());
        if (!(temp->size() == 0)) {
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
          expand(SumPtr);
          rapid_simplify(SumPtr);
          auto spin_removed = remove_spin(SumPtr);
          result->append(spin_removed);
        }
      } else {
        result->append(expr);
      }
    }  // Permutation FOR loop

    return result;
  };

  // TODO: Check if 'A' exist in expression
  // IF true: expand A;
  // ELSE continue with spintrace.

  if (expression->is<Product>()) {
    return trace_product(expression->as<Product>());
  } else if ((expression->is<Sum>())) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expression) {
      if (summand->is<Product>())
        result->append(trace_product(summand->as<Product>()));
    }
    return result;
  } else
    return nullptr;
}  // ExprPtr spintrace

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
