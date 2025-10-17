//
// Created by Eduard Valeyev on 9/22/24.
//

#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/antisymmetrizer.hpp>

namespace sequant {

antisymm_element::antisymm_element(ExprPtr ex_) {
  current_product = ex_;
  SEQUANT_ASSERT(ex_->is<Product>());
  Constant::scalar_type starting_constant = 1;
  unsigned long begining_index = 0;
  for (auto it = begin(*ex_); it != end(*ex_); it++) {
    if (it->get()->is<Tensor>()) {
      auto factor = it->get()->as<Tensor>();
      index_group.push_back(
          {begining_index, begining_index + factor.bra_rank()});
      begining_index += factor.bra_rank();
      SEQUANT_ASSERT(factor.bra_rank() == factor.ket_rank());
      for (int i = 0; i < factor.rank(); i++) {
        sorted_bra_indices.push_back(factor.bra()[i]);
        sorted_ket_indices.push_back(factor.ket()[i]);
      }
    } else if (it->get()->is<Constant>()) {
      starting_constant *= it->get()->as<Constant>().value();
    }

    else if (it->get()->is<FNOperator>()) {
      auto factor = it->get()->as<FNOperator>();
      index_group.push_back(
          {begining_index, begining_index + factor.nannihilators()});
      begining_index += factor.nannihilators();
      SEQUANT_ASSERT(factor.ncreators() == factor.nannihilators());
      for (int i = 0; i < factor.rank(); i++) {
        sorted_bra_indices.push_back(factor.annihilators()[i].index());
        sorted_ket_indices.push_back(factor.creators()[i].index());
      }
    } else {
      throw " unknown type of product, factor is not tensor, constant, or NormalOperator with FermiDirac statistics";
    }
  }

  unique_bras_list = gen_antisymm_unique(sorted_bra_indices);
  unique_kets_list = gen_antisymm_unique(sorted_ket_indices);

  auto new_sum = ex<Constant>(0);

  auto summand_exists =
      [&new_sum](ExprPtr ex) {  // check whether a summand has already been
        // generated to screen out same terms.
        bool value = false;
        for (auto summand = new_sum->begin_subexpr();
             summand != new_sum->end_subexpr(); summand++) {
          value =
              ex.get()->as<Product>() ==
              summand->get()->as<Product>();  // ensure that this equality is
          // mathematical and not hash based.
          if (value == true) {
            return value;
          }
        }
        return value;
      };

  for (int i = 0; i < unique_bras_list.size(); i++) {
    for (int j = 0; j < unique_kets_list.size(); j++) {  // product level

      auto new_product =
          ex<Constant>(unique_kets_list[i].first * unique_bras_list[j].first *
                       starting_constant);

      int index_label_pos = 0;
      for (auto it = begin(*ex_); it != end(*ex_); it++) {  // factor level
        if (it->get()->is<Constant>()) {
        }  // constant already captured in the first loop.
        else if (it->get()->is<Tensor>()) {
          auto old_tensor = it->get()->as<Tensor>();
          auto label = old_tensor.label();
          std::vector<Index> new_bras;
          std::vector<Index> new_kets;
          for (auto k = 0; k < old_tensor.rank(); k++) {  // index level
            new_bras.push_back(unique_bras_list[i].second[index_label_pos]);
            new_kets.push_back(unique_kets_list[j].second[index_label_pos]);
            index_label_pos++;
          }
          auto new_tensor = ex<Tensor>(label, bra(std::move(new_bras)),
                                       ket(std::move(new_kets)));
          new_product = new_tensor * new_product;
          new_product->canonicalize();
        }

        else if (it->get()->is<FNOperator>()) {
          auto old_Nop = it->get()->as<FNOperator>();
          std::vector<Index> new_anni;
          std::vector<Index> new_crea;
          for (auto k = 0; k < old_Nop.rank(); k++) {
            new_anni.push_back(unique_bras_list[i].second[index_label_pos]);
            new_crea.push_back(unique_kets_list[j].second[index_label_pos]);
            index_label_pos++;
          }
          auto new_Nop = ex<FNOperator>(cre(new_crea), ann(new_anni));
          new_product = new_product * new_Nop;
          // std::wcout << "product:  " << to_latex(new_product) << std::endl;
          new_product->canonicalize();
        }

        else {
          throw " unknown type of product, factor is not tensor, constant, or NormalOperator with FermiDirac statistics";
        }
      }
      if (!summand_exists(
              new_product)) {  // since products are canonicalized, repeats
        // can be found.
        new_sum = new_sum + new_product;
      }
    }
  }

  result = new_sum;
}

template <typename T>
std::vector<std::pair<int, std::vector<T>>>
antisymm_element::gen_antisymm_unique(std::vector<T> ordered_indices) {
  std::vector<std::pair<int, std::vector<T>>> result;
  // next_permutation needs an ordered list. works for integers quite well.
  // Easiest to map T to an integer corresponding to its position in a vector
  // and then go back.
  std::vector<int> ordered_numbers;
  for (int i = 0; i < ordered_indices.size(); i++) {
    ordered_numbers.push_back(i);
  }

  // return {next_exists, nswaps for this permutation} N.B.
  // nswaps_relative_to_input != nswaps_relative_to_оriginal
  auto swapcounting_tracking_next_permutation =
      [](auto first, auto last) -> std::pair<bool, int> {
    int nswaps = 0;
    using BidirIt = decltype(first);
    if (first == last) return std::make_pair(false, 0);
    BidirIt i = last;
    if (first == --i) return std::make_pair(false, 0);
    while (true) {
      BidirIt i1, i2;
      i1 = i;
      if (*--i < *i1) {
        i2 = last;
        while (!(*i < *--i2))
          ;
        std::iter_swap(i, i2);
        ++nswaps;
        std::reverse(i1, last);
        nswaps += (std::distance(i1, last)) / 2;  // logic from
        // https://en.cppreference.com/w/cpp/algorithm/reverse
        return std::make_pair(true, nswaps);
      }
      if (i == first) {
        std::reverse(first, last);
        nswaps += (std::distance(first, last)) / 2;
        return std::make_pair(false, nswaps);
      }
    }
  };

  std::vector<int> numbers = ordered_numbers;
  result.push_back({1, ordered_indices});

  int total_swaps = 0;  // even # swaps produces positive and odd # of swaps
  // produce negative
  [[maybe_unused]] int counter = 0;

  bool do_next_perm = true;

  while (do_next_perm) {
    auto do_swaps =
        swapcounting_tracking_next_permutation(begin(numbers), end(numbers));
    do_next_perm = do_swaps.first;
    total_swaps += do_swaps.second;
    if (!do_next_perm) {
      break;
    }
    auto is_canonical_sign =
        [this, &total_swaps](const auto& indices) -> std::pair<bool, bool> {
      for (auto& group :
           this->index_group) {  // There is only one sorted possibility in a
        // set (tensor) considering that no index
        // label should be the same.
        if (!is_sorted(indices.begin() + group.first,
                       indices.begin() + group.second))
          return {false,
                  total_swaps % 2 ==
                      0};  // total swaps divisible by 2 means true = positive.
      }
      return {true, total_swaps % 2 == 0};
    };
    // sieve out non-canonical terms
    auto [is_canonical, sign] = is_canonical_sign(numbers);
    if (is_canonical) {
      std::vector<T> return_vec(ordered_indices.size());
      for (int i = 0; i < numbers.size(); i++) {
        return_vec[i] = ordered_indices[numbers[i]];
      }
      if (sign) {
        result.push_back({1, return_vec});
      } else {
        result.push_back({-1, return_vec});
      }
    }
    counter += 1;
  }
  return result;
}

antisymmetrize::antisymmetrize(ExprPtr s) {
  if (s->is<Sum>()) {
    for (auto&& product :
         s->as<Sum>().summands()) {  // for each element in the sum
      if (product->is<Product>()) {
        antisymm_element ele(
            product);  // calculate the sum of all the valid permutations for
        // each term. each object here should only live until
        // this loop ends.
        result = result + ele.result;  // append that to the final list;
      } else {
        result = result + product;
      }
    }
  } else if (s->is<Product>()) {
    antisymm_element answer(ex<Product>(s->as<Product>()));
    result = answer.result;
  } else {
    result = s;
  }
}

namespace antisymm {

int num_closed_loops(std::vector<Index> init_upper,
                     std::vector<Index> init_lower,
                     std::vector<Index> new_upper,
                     std::vector<Index> new_lower) {
  int result = 0;
  auto equal_indices = [](Index i1, Index i2) {
    return i1.label() == i2.label();
  };
  auto where_same_ele = [&](Index i1, std::vector<Index> ref_list) {
    [[maybe_unused]] int hit_counter = 0;
    int where = -1;
    bool in_list = false;
    for (int i = 0; i < ref_list.size(); i++) {
      if (equal_indices(i1, ref_list[i])) {
        hit_counter += 1;
        where = i;
        in_list = true;
      }
    }
    SEQUANT_ASSERT(hit_counter < 2);
    std::pair<int, bool> result(where, in_list);
    return result;
  };

  std::vector<Index> in_loop;  // lower indices already in a loop.
  for (int i = 0; i < new_upper.size(); i++) {
    SEQUANT_ASSERT(new_upper.size() == new_lower.size());
    if (!where_same_ele(new_lower[i], in_loop).second) {
      auto starting_point = where_same_ele(new_upper[i], init_upper).first;
      auto chain_end = init_lower[starting_point];

      bool closed = equal_indices(chain_end, new_lower[i]);
      auto lower_index = new_lower[i];
      while (!closed) {
        auto where_exist = where_same_ele(
            init_upper[where_same_ele(lower_index, init_lower).first],
            new_upper);  // if the  initial upper index is the same particle as
        // the new lower index in question, find it in the new
        // upper list.
        if (!where_exist
                 .second) {  // if the upper particle index originally connected
          // to the lower index in question is not part of a
          // contraction, there is no loop with these.
          break;
        } else {
          lower_index = new_lower[where_exist.first];  // the lower index below
          // the found upper index
          in_loop.push_back(
              lower_index);  // this lower index is part of one loop so it
          // cannot be part of another.
          closed = equal_indices(chain_end,
                                 lower_index);  // is the lower index the same
          // as the end of the chain?
        }
      }
      if (closed) {
        result += 1;
      }
      in_loop.push_back(
          new_lower[i]);  // put the starting lower index in the list
    }
  }
  SEQUANT_ASSERT(in_loop.size() <= init_lower.size());
  SEQUANT_ASSERT(result <= init_lower.size());
  return result;
}

ExprPtr max_similarity(const std::vector<Index>& original_upper,
                       const std::vector<Index>& original_lower,
                       ExprPtr expression) {
  // index pairing is originally understood as a position in the original
  // vectors, but for this case, a map may do better.
  container::map<std::wstring, std::wstring> original_map;
  for (int i = 0; i < original_upper.size(); i++) {
    original_map.emplace(original_upper[i].to_latex(),
                         original_lower[i].to_latex());
  }
  for (auto&& product : expression->as<Sum>().summands()) {
    for (auto&& factor : product->as<Product>().factors()) {
      int og_pairs = 0;
      int new_pairs = 0;
      if (factor->is<Tensor>()) {
        std::vector<Index> current_upper;
        std::vector<Index> current_lower;
        if (factor->as<Tensor>().rank() == 2) {
          for (int i = 0; i < 2; i++) {
            assert(
                original_map.find(factor->as<Tensor>().ket()[i].to_latex()) !=
                original_map.end());
            if (original_map.find(factor->as<Tensor>().ket()[i].to_latex())
                    ->second == factor->as<Tensor>().bra()[i].to_latex()) {
              og_pairs += 1;
            }
            current_upper.push_back(factor->as<Tensor>().ket()[i]);
            current_lower.push_back(factor->as<Tensor>().bra()[i]);
          }
          std::iter_swap(current_lower.begin(), current_lower.begin() + 1);
          for (int i = 0; i < 2; i++) {
            SEQUANT_ASSERT(original_map.find(current_upper[i].to_latex()) !=
                           original_map.end());
            if (original_map.find(current_upper[i].to_latex())->second ==
                current_lower[i].to_latex()) {
              new_pairs += 1;
            }
          }
        }
        if (new_pairs > og_pairs) {
          factor = ex<Constant>(-1) * ex<Tensor>(factor->as<Tensor>().label(),
                                                 bra(std::move(current_lower)),
                                                 ket(std::move(current_upper)));
        }
      } else if (factor->is<FNOperator>()) {
        std::vector<Index> current_upper;
        std::vector<Index> current_lower;
        if (factor->as<FNOperator>().rank() == 2) {
          for (int i = 0; i < 2; i++) {
            SEQUANT_ASSERT(original_map.find(factor->as<FNOperator>()
                                                 .creators()[i]
                                                 .index()
                                                 .to_latex()) !=
                           original_map.end());
            if (original_map
                    .find(factor->as<FNOperator>()
                              .creators()[i]
                              .index()
                              .to_latex())
                    ->second ==
                factor->as<FNOperator>().annihilators()[i].index().to_latex()) {
              og_pairs += 1;
            }
            current_upper.push_back(
                factor->as<FNOperator>().creators()[i].index());
            current_lower.push_back(
                factor->as<FNOperator>().annihilators()[i].index());
          }
          std::iter_swap(current_lower.begin(), current_lower.begin() + 1);
          for (int i = 0; i < 2; i++) {
            SEQUANT_ASSERT(original_map.find(current_upper[i].to_latex()) !=
                           original_map.end());
            if (original_map.find(current_upper[i].to_latex())->second ==
                current_lower[i].to_latex()) {
              new_pairs += 1;
            }
          }
        }
        if (new_pairs > og_pairs) {
          factor =
              ex<Constant>(-1) * ex<FNOperator>(cre(std::move(current_upper)),
                                                ann(std::move(current_lower)));
        }
      }
    }
  }
  simplify(expression);
  return expression;
}

ExprPtr spin_sum(std::vector<Index> original_upper,
                 std::vector<Index> original_lower, ExprPtr expression,
                 bool singlet_state) {
  if (singlet_state) {
    // std::wcout << "before spin sum! :" << std::endl <<
    // to_latex_align(expression) << std::endl;
    auto init_upper = original_upper;
    auto init_lower = original_lower;
    max_similarity(original_upper, original_lower, expression);
    // may need to add separate loop if the result is a single product or
    // Operator/Tensor
    auto return_val = ex<Constant>(0);
    for (auto&& product : expression->as<Sum>().summands()) {
      auto prefactor = ex<Constant>(
          rational{1, 8});  // each term in the 3 body decomp, should have
      // (1 /2^n) prefactor where n is rank of product
      // so product rank always 3 for 3 body decomp
      std::vector<Index> new_upper;
      std::vector<Index> new_lower;
      for (auto&& factor : product->as<Product>().factors()) {
        if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"γ") {
          // prefactor = ex<Constant>(-0.5) *
          // ex<Constant>(factor->as<Tensor>().rank()) * prefactor;
          for (int i = 0; i < factor->as<Tensor>().rank(); i++) {
            new_upper.push_back(factor->as<Tensor>().ket()[i]);
            new_lower.push_back(factor->as<Tensor>().bra()[i]);
          }
          factor = ex<Tensor>(L"Γ", factor->as<Tensor>().bra(),
                              factor->as<Tensor>().ket());
        } else if (factor->is<FNOperator>()) {
          // prefactor = ex<Constant>(-0.5) *
          // ex<Constant>(factor->as<Tensor>().rank()) * prefactor;
          for (int i = 0; i < factor->as<FNOperator>().rank(); i++) {
            new_upper.push_back(factor->as<FNOperator>().creators()[i].index());
            new_lower.push_back(
                factor->as<FNOperator>().annihilators()[i].index());
          }
        }
      }
      // std::wcout << to_latex_align(product) << std::endl;
      int nloops =
          num_closed_loops(init_upper, init_lower, new_upper, new_lower);
      // std::wcout  << "nloops: " << nloops << std::endl;
      if (nloops == 0) {
      } else {
        prefactor = ex<Constant>(pow2(nloops)) * prefactor;
      }
      return_val = (product * prefactor) + return_val;
    }
    return_val->canonicalize();
    return return_val;
  } else {
    throw " non-singlet states not yet supported";
  }
}

}  // namespace antisymm

}  // namespace sequant
