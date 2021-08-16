//
// Created by Conner Masteran on 7/29/21.
//
#include <SeQuant/core/wick.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/expr.hpp>

using namespace sequant;
//various simplifications of tensor-operator expressions based on properties of particular tensors. Also functionality for restricting our operators and densities to the orbital basis set (obs).
namespace simplification{
//in various transformation methods it seems as if the constants are removed or treated separatly from the main transformed hamiltonian expression.
ExprPtr remove_const(const ExprPtr ex_){
  auto new_expression = ex<Constant>(0);
  if (ex_->is<Sum>()){
    for (auto&& product : ex_->as<Sum>().summands()){
      bool has_fnop = false;
      for (auto&& factor : product->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          has_fnop = true;
        }
      }
      if (has_fnop){ new_expression = new_expression + product;}
    }
  }
  return new_expression;
}
// based on the expressions in the paper, the tensor in the A operator never exists with lower indicies unoccupied and upper indicies occoupied
ExprPtr revert_tens_adj(ExprPtr ex_, std::wstring label){
  if (ex_->is<Sum>()){
    for (auto&& product : ex_->as<Sum>().summands()){
      for (auto&& factor : product->as<Product>().factors()){
        if(factor->is<Tensor>() && factor->as<Tensor>().label() == label){
          factor = adjoint(factor);
        }
      }
    }
  }
  else if (ex_->is<Product>()){
    for(auto&& factor : ex_->as<Product>()){
      if(factor->is<Tensor>() && factor->as<Tensor>().label() == label){
        factor = adjoint(factor);
      }
    }
  }
  else if(ex_->is<Product>() && ex_->as<Tensor>().label() == label){
    ex_ = adjoint(ex_);
  }
  return ex_;
}
auto new_idx(int pos, Index idx){
  auto new_number = std::to_wstring(pos + 6);
  auto new_string = IndexSpace::base_key(idx.space()) + L'_' + new_number;
  auto new_label = Index{new_string};
  return new_label;
}

//
Product replace_idx(ExprPtr ex_, Index og, Index newer){
  assert(ex_->is<Product>());
  auto constant = ex_->as<Product>().scalar();
  auto new_product = ex<Constant>(1);
  for(auto&& factor : ex_->as<Product>().factors()){
    if(factor->is<Tensor>()){
      std::vector<Index> new_bras;
      for(auto&& bra : factor->as<Tensor>().bra()){
        if(bra.label() == og.label()){
          new_bras.push_back(newer);
        }
        else{new_bras.push_back(bra);}
      }
      std::vector<Index> new_kets;
      for(auto&& ket : factor->as<Tensor>().ket()){
        if(ket.label() == og.label()){
          new_kets.push_back(newer);
        }
        else{new_kets.push_back(ket);}
      }
      auto new_tensor = ex<Tensor>(factor->as<Tensor>().label(), new_bras, new_kets);
      new_product = new_tensor * new_product;
    }
    if(factor->is<FNOperator>()){
      std::vector<Index> new_cres;
      for (auto&& cre : factor->as<FNOperator>().creators()){
        if(cre.index().label() == og.label()){
          new_cres.push_back(newer);
        }
        else{new_cres.push_back(cre.index());}
      }
      std::vector<Index> new_anns;
      for (auto&& ann : factor->as<FNOperator>().annihilators()){
        if(ann.index().label() == og.label()){
          new_anns.push_back(newer);
        }
        else{new_anns.push_back(ann.index());}
      }
      if (factor->as<FNOperator>().ncreators() == 1){
        auto o1 = make_overlap({L"p_7"},new_anns[0]);
        auto o3 = make_overlap(new_cres[0],{L"p_9"});
        new_product = new_product * o1 * o3;
      }
      else if(factor->as<FNOperator>().ncreators() == 2){
        auto o1 = make_overlap({L"p_7"},new_anns[0]);
        auto o2 = make_overlap({L"p_8"},new_anns[1]);
        auto o3 = make_overlap(new_cres[0],{L"p_9"});
        auto o4 = make_overlap(new_cres[1],{L"p_10"});
        new_product = new_product * o1 * o2 * o3 * o4;
      }
      else{throw "does not handle size > 2";}
      //auto new_op = ex<FNOperator>(new_cres,new_anns);
      // new_product = new_product * new_op;
    }
  }
  auto result = (ex<Constant>(constant) * new_product);
  return result->as<Product>();
}

ExprPtr op_to_tens(ExprPtr ex_){
  assert(ex_->is<FNOperator>());
  std::vector<Index> bra_indices;
  std::vector<Index> ket_indices;
  for(auto&& cre : ex_->as<FNOperator>().creators()){
    bra_indices.push_back(cre.index());
  }
  for(auto&& ann : ex_->as<FNOperator>().annihilators()){
    ket_indices.push_back(ann.index());
  }
  auto label = get_default_context().spbasis() == SPBasis::spinfree ? L"E" : L"a";
  auto result = ex<Tensor>(label, bra_indices, ket_indices);
  return result;
}

ExprPtr tens_to_op(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  auto result = ex<FNOperator>(ex_->as<Tensor>().ket(),ex_->as<Tensor>().bra());
  return result;
}

ExprPtr screen_F_tensors(ExprPtr ex_) {// F tensors must contain contain indices in the bra with space > all. this includes complete, completeunoccupied, and inactiveunoccupied.
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().label() == L"F");
  bool good = false;
  for (auto&& bra : ex_->as<Tensor>().bra()){
    if(bra.space().type() == IndexSpace::complete || bra.space().type() == IndexSpace::complete_unoccupied || bra.space().type() == IndexSpace::inactive_unoccupied){
      good = true;
    }
  }

  for (auto&& ket : ex_->as<Tensor>().ket()){
    if(ket.space().type() == IndexSpace::complete || ket.space().type() == IndexSpace::complete_unoccupied || ket.space().type() == IndexSpace::inactive_unoccupied){
      good = true;
    }
  }
  if(good){
    return ex_;
  }
  else{
    return ex<Constant>(0);
  }
}

auto screen_densities(ExprPtr ex_){// densities probably should be non-zero if each index has a chance to be occupied, in other words, screen out densities containing unoccupied labels.
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().label() == L"\\Gamma");
  bool good = true;
  for (auto&& bra : ex_->as<Tensor>().bra()){
    if (bra.space().type() == IndexSpace::unoccupied || bra.space().type() == IndexSpace::complete_unoccupied){
      good = false;
    }
  }
  for (auto&& ket : ex_->as<Tensor>().ket()){
    if (ket.space().type() == IndexSpace::unoccupied || ket.space().type() == IndexSpace::complete_unoccupied){
      good = false;
    }
  }
  if(good){
    return ex_;
  }
  else{return ex<Constant>(0);}
};

//based on Brillouin's Theory, the fock matrix should be diagonal in the obs. However, it seems that this is not generally true when denisity is not involved and when dealing with two indices outside of the obs.
auto treat_fock(ExprPtr ex_){
  auto new_ex_ = ex<Constant>(0);
  for (auto&& product : ex_->as<Sum>().summands()){
    auto new_product = ex<Constant>(1);
    auto scalar = product->as<Product>().scalar();
    bool contains_density = false;
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        contains_density = true;
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"f" && contains_density && (factor->as<Tensor>().bra()[0].space().type() != IndexSpace::complete_unoccupied || factor->as<Tensor>().ket()[0].space().type() != IndexSpace::complete_unoccupied)){
        auto overlap = make_overlap(factor->as<Tensor>().bra()[0],factor->as<Tensor>().ket()[0]);
        new_product = overlap * factor * new_product;
      }
      else new_product = new_product * factor;
    }
    new_product = ex<Constant>(scalar) * new_product;
    non_canon_simplify(new_product);
    new_ex_ = new_ex_ + new_product;
  }
  non_canon_simplify(new_ex_);
  FWickTheorem wick{new_ex_};
  wick.reduce(new_ex_);
  non_canon_simplify(new_ex_);
  return new_ex_;
}
//TODO generalize for arbitrary number of density matrices.
//TODO generalize for spin-orbital basis
//simplification to deal with hamiltonian based expressions. involving one body h and two body g tensors.
// not rigorous for more than 2 body operators or more than 2 density matrices whose rank must be <= 2.
// unfortunately, simplify(result) and  wick.reduce(result) will recanonicalize the indices.
// enforces the following obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}}
// should allow analysis of multiple expressions who have the same normal order operator prefactor.
ExprPtr hamiltonian_based(ExprPtr exprs){
  exprs = remove_const(exprs); // remove constant terms from the expression
  simplify(exprs);

  auto overlap_expr = ex<Constant>(0); //enforce an overlap each E with elements from
  for (auto&& product : exprs->as<Sum>().summands()){// may be able to make_overlaps manually and apply them to the products. simplify may know what to do with it.
    auto new_product = ex<Constant>(1);
    auto scalar = product->as<Product>().scalar();
    bool first_density = false; // some products have multiple densities
    bool third_density = false;
    for (int it = product->as<Product>().factors().size() - 1; it >= 0; it--){//loop backwards through the products.
      auto factor = product->as<Product>().factor(it);

      if (it == product->as<Product>().factors().size() - 1 && factor->is<FNOperator>() &&factor->as<FNOperator>().rank() == 2){
        std::wstring label_1;
        std::wstring label_2;
        std::wstring label_3;
        std::wstring label_4;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_2 = factor->as<FNOperator>().annihilators()[1].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        label_4 = factor->as<FNOperator>().creators()[1].index().label();
        auto o1 = make_overlap(Index{L"p_1"},Index{label_1});
        auto o2 = make_overlap(Index{L"p_2"},Index{label_2});
        auto o3 = make_overlap(Index{L"p_3"},Index{label_3});
        auto o4 = make_overlap(Index{L"p_4"},Index{label_4});
        new_product = o1 * o2 * o3 * o4 * new_product * op_to_tens(factor);
      }
      else if (it == product->as<Product>().factors().size() - 1 && factor->is<FNOperator>() && factor->as<FNOperator>().rank() == 1){
        std::wstring label_1;
        std::wstring label_3;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        auto o1 = make_overlap(Index{L"p_1"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_3"},Index{label_3});
        new_product = o1 * o3 * new_product * op_to_tens(factor);
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1 && !first_density){
        first_density = true;
        std::wstring label_1;
        std::wstring label_3;
        label_1 = factor->as<Tensor>().ket()[0].label();
        label_3 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index{L"p_10"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_11"},Index{label_3});
        new_product =  o1 * o3 * factor * new_product;
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1 && first_density){
        first_density = false;
        std::wstring label_2;
        std::wstring label_4;
        label_2 = factor->as<Tensor>().ket()[0].label();
        label_4 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index{L"p_9"},Index{label_2});
        auto o3 = make_overlap(Index{L"p_12"},Index{label_4});
        new_product =  o1 * o3 * factor * new_product;
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1 && third_density){
        first_density = false;
        std::wstring label_2;
        std::wstring label_4;
        label_2 = factor->as<Tensor>().ket()[0].label();
        label_4 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index{L"p_13"},Index{label_2});
        auto o3 = make_overlap(Index{L"p_14"},Index{label_4});
        new_product =  o1 * o3 * factor * new_product;
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 2){
        std::wstring label_1;
        std::wstring label_3;
        std::wstring label_2;
        std::wstring label_4;
        label_1 = factor->as<Tensor>().ket()[0].label();
        label_3 = factor->as<Tensor>().bra()[0].label();
        label_2 = factor->as<Tensor>().ket()[1].label();
        label_4 = factor->as<Tensor>().bra()[1].label();
        auto o1 = make_overlap(Index{L"p_5"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_6"},Index{label_3});
        auto o2 = make_overlap(Index{L"p_7"},Index{label_2});
        auto o4 = make_overlap(Index{L"p_8"},Index{label_4});
        new_product =  o1 * o3 * o2 * o4 * factor * new_product;
      }
      else{new_product = new_product * factor;}
    }
    new_product = ex<Constant>(scalar) * new_product;
    overlap_expr = overlap_expr + new_product;
  }

  FWickTheorem wick{overlap_expr};
  wick.reduce(overlap_expr);
  non_canon_simplify(overlap_expr);

  for (auto&& product : overlap_expr->as<Sum>().summands()){
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F"){
        factor = screen_F_tensors(factor);
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        factor = screen_densities(factor);
      }
    }
  }
  non_canon_simplify(overlap_expr);

  auto result = ex<Constant>(0);
  for (auto&& product : overlap_expr->as<Sum>().summands()){
    auto new_product = ex<Constant>(product->as<Product>().scalar());
   // std::vector<std::pair<Index,Index>> og_new;
    for (auto&& factor : product->as<Product>().factors()){ //loop through the factors once to learn about the indices.
      if(factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" || factor->as<Tensor>().label() == L"a")){
        factor = tens_to_op(factor);
        /*int cre_iter = 1;
        for (auto&& cre : factor->as<FNOperator>().creators()){
          og_new.push_back({cre.index(),new_idx(cre_iter, cre.index())});
          cre_iter +=1;

        }
        int ann_iter = 3;
        for(auto&& ann : factor->as<FNOperator>().annihilators()){
          og_new.push_back({ann.index(),new_idx(ann_iter, ann.index())});
          ann_iter +=1;
        }*/
        if(factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().ncreators() == 1) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            new_product = new_product * o1 * o3;
          } else if (factor->as<FNOperator>().ncreators() == 2) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o2 = make_overlap(
                {L"p_8"}, factor->as<FNOperator>().annihilators()[1].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            auto o4 = make_overlap(
                factor->as<FNOperator>().creators()[1].index(), {L"p_10"});
            new_product = new_product * o1 * o2 * o3 * o4;
          }
        }
      }
      else{new_product = factor * new_product;}
    }
    result = new_product + result;
    /*for (auto&& pair : og_new){// loop through again to update the indices
      product->as<Product>() = replace_idx(product, pair.first, pair.second);
    }*/
    //og_new.resize(0);
  }

  non_canon_simplify(result);

  return result;
  //std::wcout << "[H,R]12" << to_latex_align(overlap_expr, 20, 3) << std::endl;

}
//TODO generalize for arbitrary number of density matrices.
//TODO generalize for spin-orbital basis
//simplification to deal with fock based expressions. involving one body fock operator.
// not rigorous for more than 2 body operators or more than 2 density matrices whose rank must be <= 2.
// unfortunately, simplify(result) and wick.reduce(result) will recanonicalize the indices.
// enforces the following obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}}
// should allow analysis of multiple expressions who have the same normal order operator prefactor.

ExprPtr fock_based (ExprPtr exprs){
  exprs = remove_const(exprs);
  auto double_overlap_expr = ex<Constant>(0); //enforce an overlap each E with elements from
  for (auto&& product : exprs->as<Sum>().summands()){// may be able to make_overlaps manually and apply them to the products. simplify may know what to do with it.
    auto new_product = ex<Constant>(1);
    auto scalar = product->as<Product>().scalar();
    bool first_density = false; // some products have multiple densities
    for (int it = product->as<Product>().factors().size() - 1; it >= 0; it--){//loop backwards through the products.
      auto factor = product->as<Product>().factor(it);

      if (it == product->as<Product>().factors().size() - 1 && factor->as<FNOperator>().rank() == 2){
        std::wstring label_1;
        std::wstring label_2;
        std::wstring label_3;
        std::wstring label_4;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_2 = factor->as<FNOperator>().annihilators()[1].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        label_4 = factor->as<FNOperator>().creators()[1].index().label();
        auto o1 = make_overlap(Index{L"p_1"},Index{label_1});
        auto o2 = make_overlap(Index{L"p_2"},Index{label_2});
        auto o3 = make_overlap(Index{L"p_3"},Index{label_3});
        auto o4 = make_overlap(Index{L"p_4"},Index{label_4});
        new_product = o1 * o2 * o3 * o4 * new_product * op_to_tens(factor);
      }
      else if (it == product->as<Product>().factors().size() - 1 && factor->as<FNOperator>().rank() == 1){
        std::wstring label_1;
        std::wstring label_3;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        auto o1 = make_overlap(Index{L"p_1"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_3"},Index{label_3});
        new_product = o1 * o3 * new_product * op_to_tens(factor);
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1 && !first_density){
        first_density = true;
        std::wstring label_1;
        std::wstring label_3;
        label_1 = factor->as<Tensor>().ket()[0].label();
        label_3 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index{L"p_5"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_6"},Index{label_3});
        new_product =  o1 * o3 * factor * new_product;
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1 && first_density){
        std::wstring label_2;
        std::wstring label_4;
        label_2 = factor->as<Tensor>().ket()[0].label();
        label_4 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index{L"p_7"},Index{label_2});
        auto o3 = make_overlap(Index{L"p_8"},Index{label_4});
        new_product =  o1 * o3 * factor * new_product;
      }
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 2){
        std::wstring label_1;
        std::wstring label_3;
        std::wstring label_2;
        std::wstring label_4;
        label_1 = factor->as<Tensor>().ket()[0].label();
        label_3 = factor->as<Tensor>().bra()[0].label();
        label_2 = factor->as<Tensor>().ket()[1].label();
        label_4 = factor->as<Tensor>().bra()[1].label();
        auto o1 = make_overlap(Index{L"p_5"},Index{label_1});
        auto o3 = make_overlap(Index{L"p_6"},Index{label_3});
        auto o2 = make_overlap(Index{L"p_7"},Index{label_2});
        auto o4 = make_overlap(Index{L"p_8"},Index{label_4});
        new_product =  o1 * o3 * o2 * o4 * factor * new_product;
      }
      else{new_product = new_product * factor;}
    }
    new_product = ex<Constant>(scalar) * new_product;
    double_overlap_expr = double_overlap_expr + new_product;
  }

  FWickTheorem wick3{double_overlap_expr};
  wick3.reduce(double_overlap_expr);
  non_canon_simplify(double_overlap_expr);

  //std::wcout << "new: " << std::endl << to_latex_align(double_overlap_expr,20,3) << std::endl;
  for (auto&& product : double_overlap_expr->as<Sum>().summands()){
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F"){
        factor = screen_F_tensors(factor);
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        factor = screen_densities(factor);
      }
    }
  }
  non_canon_simplify(double_overlap_expr);
  double_overlap_expr = treat_fock(double_overlap_expr);
  for (auto&& product : double_overlap_expr->as<Sum>().summands()){
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F"){
        factor = screen_F_tensors(factor);
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        factor = screen_densities(factor);
      }
    }
  }
  non_canon_simplify(double_overlap_expr);

  auto result = ex<Constant>(0);
  for (auto&& product : double_overlap_expr->as<Sum>().summands()){
    auto new_product = ex<Constant>(product->as<Product>().scalar());
    // std::vector<std::pair<Index,Index>> og_new;
    for (auto&& factor : product->as<Product>().factors()){ //loop through the factors once to learn about the indices.
      if(factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" || factor->as<Tensor>().label() == L"a")){
        factor = tens_to_op(factor);
        /*int cre_iter = 1;
        for (auto&& cre : factor->as<FNOperator>().creators()){
          og_new.push_back({cre.index(),new_idx(cre_iter, cre.index())});
          cre_iter +=1;

        }
        int ann_iter = 3;
        for(auto&& ann : factor->as<FNOperator>().annihilators()){
          og_new.push_back({ann.index(),new_idx(ann_iter, ann.index())});
          ann_iter +=1;
        }*/
        if(factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().ncreators() == 1) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            new_product = new_product * o1 * o3;
          } else if (factor->as<FNOperator>().ncreators() == 2) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o2 = make_overlap(
                {L"p_8"}, factor->as<FNOperator>().annihilators()[1].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            auto o4 = make_overlap(
                factor->as<FNOperator>().creators()[1].index(), {L"p_10"});
            new_product = new_product * o1 * o2 * o3 * o4;
          }
        }
      }
      else{new_product = factor * new_product;}
    }
    result = new_product + result;
    /*for (auto&& pair : og_new){// loop through again to update the indices
      product->as<Product>() = replace_idx(product, pair.first, pair.second);
    }*/
    //og_new.resize(0);
  }

  non_canon_simplify(result);

  return result;
}
}
#ifndef SEQUANT_SIMPLIFICATIONS_H
#define SEQUANT_SIMPLIFICATIONS_H

#endif  // SEQUANT_SIMPLIFICATIONS_H
