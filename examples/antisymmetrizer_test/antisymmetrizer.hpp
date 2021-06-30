//
// created by Conner Masteran 06/1/2021
//

#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <iostream>
#include <vector>
#include <algorithm>


namespace sequant{


//@brief generates all unique permutations of a product where products differing only by internal tensor antisymmetry are non-unique.
// i.e. a^{i_1 i_2}_{a_1 a_2} = - a^{i_2 i_1}_{a_1 a_2}. RHS is non-unique in this context.
class antisymm_element{
  using IndexGroup = std::pair<size_t, size_t>;

 private:

  std::vector<IndexGroup> index_group;// where each tensor begins and ends. assumes particle conserving ops. needed to keep track of canonical ordering
  std::vector<std::pair<int, std::vector<Index>>> unique_bras_list; // list of unique bra orderings with associated integer for the sign
  std::vector<std::pair<int, std::vector<Index>>> unique_kets_list; // list of unique ket orderings with associated integer for the sign
  ExprPtr current_product; // used to keep track of the original expression recieved by the constructor


  //generates all possible permutations while observing the canonical ordering of each tensor/operator
  //function kept general to work with other data types, algorithm does not require sequant::Index objects
  // @param a list of T objects/indices in the canonical ordering
  // @return a list of all permutations each with an associated sign
  template <typename T>
  std::vector<std::pair<int, std::vector<T>>> gen_antisymm_unique(std::vector<T> ordered_indices) {
    std::vector<std::pair<int, std::vector<T>>> result;
    //next_permutation needs an ordered list. works for integers quite well.
    // Easiest to map T to an integer corresponding to its position in a vector and then go back.
    std::vector<int> ordered_numbers;
    for (int i = 0; i < ordered_indices.size(); i++){
      ordered_numbers.push_back(i);
    }

    // return {next_exists, nswaps for this permutation} N.B. nswaps_relative_to_input != nswaps_relative_to_оriginal
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
          nswaps += (std::distance(i1, last)) / 2; // logic from https://en.cppreference.com/w/cpp/algorithm/reverse
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
    result.push_back({1,ordered_indices});

    int total_swaps = 0; // even # swaps produces positive and odd # of swaps produce negative
    int counter = 0;

    bool do_next_perm = true;

    while (do_next_perm) {
      auto do_swaps = swapcounting_tracking_next_permutation(begin(numbers),end(numbers));
      do_next_perm = do_swaps.first;
      total_swaps += do_swaps.second;
      if (!do_next_perm){break;}
      auto is_canonical_sign = [this, &total_swaps](const auto &indices) -> std::pair<bool,bool> {
        for (auto &group : this->index_group) {  // There is only one sorted possibility in a set (tensor) considering that no index label should be the same.
          if (!is_sorted(indices.begin() + group.first,
                         indices.begin() + group.second))
            return {false, total_swaps % 2 == 0}; // total swaps divisible by 2 means true = positive.
        }
        return {true, total_swaps % 2 == 0};
      };
      // sieve out non-canonical terms
      auto [is_canonical, sign] = is_canonical_sign(numbers);
      if (is_canonical){
        std::vector<T> return_vec(ordered_indices.size());
        for (int i = 0; i < numbers.size(); i++){
          return_vec[i] = ordered_indices[numbers[i]];
        }
        if(sign){result.push_back({1,return_vec});}
        else{result.push_back({-1,return_vec});}

      }
      counter += 1;
    }
    return result;
  }

 public:

  //takes a sequant::ExprPtr and generates all antisymmetric unique permutations of that expression.
  //requires that ex_ is a product expression at this point
  // @param ex_ as product
  // populates a result ExprPtr that the user can grab. result is in general a Sum.
  antisymm_element(ExprPtr ex_){
    current_product = ex_;
    assert(ex_->is<Product>());
    auto starting_constant = 1.0;
    unsigned long begining_index = 0;
    for(auto it = begin(*ex_); it!= end(*ex_); it++) {
      if(it->get()->is<Tensor>()){
        auto factor = it->get()->as<Tensor>();
        index_group.push_back({begining_index, begining_index + factor.bra_rank()});
        begining_index += factor.bra_rank();
        assert(factor.bra_rank() == factor.ket_rank());
        for (int i = 0; i < factor.rank(); i++){
          sorted_bra_indices.push_back(factor.bra()[i]);
          sorted_ket_indices.push_back(factor.ket()[i]);
        }
      }
      else if( it->get()->is<Constant>()){
        starting_constant *= it->get()->as<Constant>().value().real(); // not sure what to do for imaginary part if needed
      }

      else if(it->get()->is<NormalOperator<Statistics::FermiDirac>>()){
        auto factor = it->get()->as<NormalOperator<Statistics::FermiDirac>>();
        index_group.push_back({begining_index, begining_index + factor.nannihilators()});
        begining_index += factor.nannihilators();
        assert(factor.ncreators() == factor.nannihilators());
        for (int i = 0; i < factor.rank(); i++){
          sorted_bra_indices.push_back(factor.creators()[i].index());
          sorted_ket_indices.push_back(factor.annihilators()[i].index());
        }
      }
      else{
        throw " unknown type of product, factor is not tensor, constant, or NormalOperator with FermiDirac statistics";
      }
    }

    unique_bras_list = gen_antisymm_unique(sorted_bra_indices);
    unique_kets_list = gen_antisymm_unique(sorted_ket_indices);

    auto new_sum = ex<Constant>(0.0);

    auto summand_exists = [&new_sum](ExprPtr ex){ // check whether a summand has already been generated to screen out same terms.
        bool value = false;
        for (auto summand = new_sum->begin_subexpr(); summand != new_sum->end_subexpr(); summand++){
          value = ex.get()->as<Product>() == summand->get()->as<Product>(); // ensure that this equality is mathematical and not hash based.
          if (value == true){
            return value;
          }
        }
        return value;
    };

    for (int i  = 0; i < unique_bras_list.size();i++){
      for (int j = 0; j < unique_kets_list.size(); j++){ // product level


        auto new_product = ex<Constant>(unique_kets_list[i].first * unique_bras_list[j].first * starting_constant);

        int index_label_pos = 0;
        for (auto it = begin(*ex_); it!= end(*ex_); it++){ // factor level
          if (it->get()->is<Constant>()){} // constant already captured in the first loop.
          else if (it->get()->is<Tensor>()){
            auto old_tensor = it->get()->as<Tensor>();
            auto label = old_tensor.label();
            std::vector<Index> new_bras;
            std::vector<Index> new_kets;
            for (auto k = 0; k < old_tensor.rank(); k++){ // index level
              new_bras.push_back(unique_bras_list[i].second[index_label_pos]);
              new_kets.push_back(unique_kets_list[j].second[index_label_pos]);
              index_label_pos++;
            }
            auto new_tensor = ex<Tensor>(label, new_bras, new_kets);
            new_product = new_product * new_tensor;
            new_product->canonicalize();
          }
          else if (it->get()->is<NormalOperator<Statistics::FermiDirac>>()){
            auto old_Nop = it->get()->as<NormalOperator<Statistics::FermiDirac>>();
            std::vector<Index> new_cre;
            std::vector<Index> new_ann;
            for (auto k = 0; k < old_Nop.rank(); k++){
              new_cre.push_back(unique_bras_list[i].second[index_label_pos]);
              new_ann.push_back(unique_kets_list[j].second[index_label_pos]);
              index_label_pos++;
            }
            auto new_Nop = ex<NormalOperator<Statistics::FermiDirac>>(new_cre, new_ann);
            new_product = new_product * new_Nop;
            new_product->canonicalize();
          }
          else{
            throw " unknown type of product, factor is not tensor, constant, or NormalOperator with FermiDirac statistics";
          }
        }
        if (!summand_exists(new_product)) { // since products are canonicalized, repeats can be found.
          new_sum = new_sum + new_product;
        }
      }
    }

    result = new_sum;
  }

  std::vector<Index> sorted_bra_indices; // The original order of the upper indices on a given term
  std::vector<Index> sorted_ket_indices; // the original order of the lower indices on a given term
  ExprPtr result;

};

//@brief simple class to call antisymm_element on only products
class antisymmetrize {
 public:
  ExprPtr result = ex<Constant>(0);
  antisymmetrize(ExprPtr s) {
    if (s->is<Sum>()) {
      for (auto it = s->as<Sum>().begin_subexpr(); it != s->as<Sum>().end_subexpr();it++) {// for each element in the sum
        antisymm_element ele(ex<Product>(it->get()->begin_subexpr(), it->get()->end_subexpr()));  // calculate the sum of all the valid permutations for each term. each object here should only live until this loop ends.
        result = result + ele.result;// append that to the final list;
      }
    }
    else if(s->is<Product>()){
      antisymm_element answer(ex<Product>(s->as<Product>()));
      result = answer.result;
    }
    else{result = s;}
  }

  //function which counts number of closed loops from the starting order of upper and lower indices and the contracted indices.
  //since the ordering of the new contracted indices is arbitrary, the algorithm searched for the upper index which would connect to the lower index
  // checks if contracted lower index closes the loop, if not, continue searching until the corresponding upper index is not present, or the loop closes.
  // keep track of which indices are used so that loops are not double counted
  //@param1 initial order of upper indices before contraction, @param 2 initial order of lower indices before contraction, @param 3 set of contracted upper indices, @param 4 set of lower contracted indices.
  // @out the number of loops.
  //TODO Test this function extensively and add more asserts
  int num_closed_loops(std::vector<Index> init_upper, std::vector<Index> init_lower, std::vector<Index> new_upper, std::vector<Index> new_lower){
    int result = 0;
    auto equal_indices = [](Index i1, Index i2){
      return i1.label() == i2.label();
    };
    auto where_same_ele = [&](Index i1, std::vector<Index> ref_list){
      int hit_counter = 0;
      int where;
      bool in_list = false;
      for (int i =0; i < ref_list.size(); i++){
        if (equal_indices(i1, ref_list[i])){
          hit_counter += 1;
          where = i;
          in_list = true;
        }
      }
      assert(hit_counter < 2);
      std::pair<int,bool> result(where, in_list);
      return result;
    };

    std::vector<Index> in_loop; // lower indices already in a loop.
    for (int i =0; i < init_upper.size(); i++){
      if(where_same_ele(new_lower[i],in_loop).second) {
        auto chain_end =
                      init_lower[where_same_ele(new_upper[i],
                                                init_upper).first];

        bool closed = equal_indices(chain_end, new_lower[i]);
        auto lower_index = new_lower[i];
        while (!closed) {
          auto where_exist = where_same_ele(
                              init_upper[where_same_ele(lower_index,
                                        init_lower).first],
                              new_upper);  // if the  initial upper index is the same particle as the new lower index in question, find it in the new upper list.
          if (!where_exist.second) {  // if the upper particle index originally connected to the lower index in question is not part of a contraction, there is no loop with these.
            break;
          } else {
            lower_index =
                new_lower[where_exist.first];  // the lower index below the found upper index
            in_loop.push_back(lower_index); // this lower index is part of one loop so it cannot be part of another.
            closed = equal_indices(chain_end, lower_index);  // is the lower index the same as the end of the chain?
          }
        }
        if (closed) {
          result += 1;
        }
        in_loop.push_back(new_lower[i]); // put the starting lower index in the list
      }
    }
    assert(in_loop.size() <= init_lower.size());
    assert(result <= init_lower.size());
    return result;
  }

  // not a general spin-summing procedure, implementation for a known singlet state for the prefactor rules to apply.
  // TODO: use a generalized spin summation for non-singlet states
  ExprPtr spin_sum(antisymm_element &ele_result, bool singlet_state = true){
    if (singlet_state) {
      auto init_upper = ele_result.sorted_bra_indices;
      auto init_lower = ele_result.sorted_ket_indices;
      //may need to add separate loop if the result is a single product or Operator/Tensor
      auto return_val = ex<Constant>(0);
      for(auto&& product : ele_result.result->as<Sum>().summands()){
          auto prefactor = ex<Constant>(1.0);
          std::vector<Index> new_upper;
          std::vector<Index> new_lower;
          for(auto&& factor : product->as<Product>().factors()){
            if(factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\gamma") {
              prefactor = ex<Constant>(-0.5) *
                          ex<Constant>(factor->as<Tensor>().rank()) * prefactor;

              for (int i = 0; i < factor->as<Tensor>().rank(); i++) {
                new_upper.push_back(factor->as<Tensor>().bra()[i]);
                new_lower.push_back(factor->as<Tensor>().ket()[i]);
              }
              factor->as<Tensor>().label() = L"\\Gamma";
            }
          }
          int nloops = num_closed_loops(init_upper, init_lower, new_upper, new_lower);
          if (nloops == 0){
          }
          else{
            prefactor = ex<Constant>(nloops * 2) * prefactor;
          }
          return_val = (product * prefactor) + return_val;
      }
      return_val->canonicalize();
      return return_val;
    }
    else{throw " non-singlet states not yet supported";}
  }
};
}
