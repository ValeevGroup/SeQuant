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
template<typename vec_type>
std::pair<bool,int> in_list(Index idx, vec_type ref_list){
  bool inlist = false;
  int where_inlist = 0;
  for (int i = 0; i < ref_list.size(); i++){
    if (idx.label() == ref_list[i].label()){
      inlist = true;
      where_inlist = i;
    }
  }
  std::pair<bool, int> result{inlist,where_inlist};
  return result;
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

ExprPtr overlap_with_obs(ExprPtr ex_){
  auto overlap_expr = ex<Constant>(0); //enforce an overlap each E with elements from
  for (auto&& product : ex_->as<Sum>().summands()){// may be able to make_overlaps manually and apply them to the products. simplify may know what to do with it.
    auto new_product = ex<Constant>(product->as<Product>().scalar());
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
      else if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma" && factor->as<Tensor>().rank() == 1){
        std::wstring label_2;
        std::wstring label_4;
        label_2 = factor->as<Tensor>().ket()[0].label();
        label_4 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)),Index{label_2});
        auto o3 = make_overlap(Index{label_4},Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)));
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
        auto o1 = make_overlap(Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)),Index{label_1});
        auto o3 = make_overlap(Index{label_3}, Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)));
        auto o2 = make_overlap(Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)),Index{label_2});
        auto o4 = make_overlap(Index{label_4},Index::make_tmp_index(IndexSpace::instance(
            IndexSpace::all)));
        new_product =  o1 * o3 * o2 * o4 * factor * new_product;
      }
      else{new_product = new_product * factor;}
    }
    new_product = new_product;
    overlap_expr = overlap_expr + new_product;
  }

  FWickTheorem wick{overlap_expr};
  wick.reduce(overlap_expr);
  simplify(overlap_expr);
  return overlap_expr;
}

using IDX_list = std::initializer_list<Index>;
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
  simplify(new_expression);
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

ExprPtr tens_to_op(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  auto result = ex<FNOperator>(ex_->as<Tensor>().ket(),ex_->as<Tensor>().bra());
  return result;
}
// F tensors must contain contain indices in the bra with space > all. this includes complete, completeunoccupied, and inactiveunoccupied.
// and if one of the particle indices is connected to the obs virtual space, then the other must be from the CABS set. i.e. if G^{a \beta}_{ij} -> G^{a a'}_{ij}
ExprPtr screen_F_tensors(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().label() == L"F");
  auto overlap = ex<Constant>(1);
  bool good = false;
  bool bra_good = false;
  for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++) {
    auto bra = ex_->as<Tensor>().bra()[i];
    if (bra.space().type() == IndexSpace::complete || bra.space().type() == IndexSpace::complete_unoccupied) {
      good = true;
      bra_good = true;
    }
    else if(bra.space().type() == IndexSpace::complete ||
             bra.space().type() == IndexSpace::complete_unoccupied ||
             bra.space().type() == IndexSpace::other_unoccupied){
      good = true;
    }
  }
    for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++){
      auto bra = ex_->as<Tensor>().bra()[i];
    if ((bra.space().type() == IndexSpace::unoccupied || bra.space().type() ==
        IndexSpace::all) && bra_good && (bra.space().type() != IndexSpace::complete_unoccupied)) {  // if one of the upper indices is explicitly outside of CABs, create an overlap with the other index and the CABs space.
      if (i == 0) {
        overlap = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                   IndexSpace::other_unoccupied)),
                               ex_->as<Tensor>().bra()[1]);
      }
      if (i == 1) {
        overlap = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                   IndexSpace::other_unoccupied)),
                               ex_->as<Tensor>().bra()[0]);
      }
    }
  }

  bool ket_good = false;
  for (int j = 0; j <ex_->as<Tensor>().ket().size(); j++) {
    auto ket = ex_->as<Tensor>().ket()[j];
    if (ket.space().type() == IndexSpace::complete ||
        ket.space().type() == IndexSpace::complete_unoccupied) {
      good = true;
      ket_good = true;
    }
    else if (ket.space().type() == IndexSpace::complete ||
        ket.space().type() == IndexSpace::complete_unoccupied ||
        ket.space().type() == IndexSpace::other_unoccupied) {
      good = true;
    }
  }
    for (int j = 0; j <ex_->as<Tensor>().ket().size(); j++){
      auto ket = ex_->as<Tensor>().ket()[j];
    if ((ket.space().type() == IndexSpace::unoccupied || ket.space().type() == IndexSpace::all) && ket_good && (ket.space().type() != IndexSpace::complete_unoccupied)) {
      if (j == 0) {
        overlap = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                   IndexSpace::other_unoccupied)),
                               ex_->as<Tensor>().ket()[1]);
      }
      if (j == 1) {
        overlap = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                   IndexSpace::other_unoccupied)),
                               ex_->as<Tensor>().ket()[0]);
      }
    }
  }
  if(good){
    return ex_ * overlap;
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

//based on Brillouin's Theory, the fock matrix should be block diagonal.
// generalized says that complete unoccupied might be non-zero with obs unocc, but zero with occ and complete unocc.
// lets assue normal Brillouin's Theory.
auto treat_fock(ExprPtr ex_){
  auto new_ex_ = ex<Constant>(0);
  for (auto&& product : ex_->as<Sum>().summands()){
    auto new_product = ex<Constant>(product->as<Product>().scalar());
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"f"){
        auto space = intersection(factor->as<Tensor>().bra()[0].space(), factor->as<Tensor>().ket()[0].space());
        if(space.type().none()){
          new_product = ex<Constant>(0) * new_product;
        }
        else {
          auto bra_index =
              Index::make_tmp_index(IndexSpace::instance(space.type()));
          auto ket_index =
              Index::make_tmp_index(IndexSpace::instance(space.type()));

          auto bra_overlap =
              make_overlap(factor->as<Tensor>().bra()[0], bra_index);
          auto ket_overlap =
              make_overlap(factor->as<Tensor>().ket()[0], ket_index);
          new_product = bra_overlap * ket_overlap * factor * new_product;
        }
      }
      else new_product = new_product * factor;
    }
    //std::wcout << "problematic product: " << to_latex_align(new_product) << std::endl;
    //simplify(new_product);
    new_ex_ = new_ex_ + new_product;
  }
  //simplify(new_ex_);
   FWickTheorem wick{new_ex_};
  wick.reduce(new_ex_);
  non_canon_simplify(new_ex_);

  return new_ex_;
}
//for a pair of tensors returns relavant information used to specially treat certain products.
// return tuple<number_connected indices, connected space, bool used to decide where the indices go, t1 external indices, t2 external indices>
std::tuple<int, IndexSpace::Type,bool,std::vector<Index>,std::vector<Index>> ncon_spa_updom_t1ext_t2ext(Tensor T1, Tensor T2,bool print_ = false){
  auto space = IndexSpace::occupied; // just a default that will not be used.
  std::vector<Index> T1_is;
  std::vector<Index> T2_is;
  // order will matter here for internal indices.
  std::vector<Index> external_T1;
  std::vector<Index> external_T2;
  //
  std::vector<Index> unordered_internal_T1;
  std::vector<Index> unordered_internal_T2;

  std::vector<Index> connected_indices;
  //not only do I need to know the order which the indices will appear, but I also need to know upper and lower dominance. for each tensor.
  bool upper_dominant; // 0 means T1 is upper dominant, 1 means T2 is upper dominant.
  int nconnects = 0;
  for (auto&& bra : T1.bra()){
    T1_is.push_back(bra);
  }
  for (auto&& ket : T1.ket()){
    T1_is.push_back(ket);
  }
  for (auto&& bra : T2.bra()){
    T2_is.push_back(bra);
  }
  for (auto&& ket : T2.ket()){
    T2_is.push_back(ket);
  }
  for (int i1 = 0; i1 < T1_is.size(); i1++){
    for ( int i2 = 0; i2 <T2_is.size(); i2++){
      if(T1_is[i1].label() == T2_is[i2].label()){
        nconnects += 1;
        space = T1_is[i1].space().type();
        connected_indices.push_back(T1_is[i1]);
      }
    }
  }
  if ( nconnects == 0){
    external_T1 = T1_is;
    external_T2 = T2_is;
    return {nconnects, space,upper_dominant, external_T1, external_T2};
  }
auto not_inlist = [&](std::vector<Index> ref, std::vector<Index> pool){
    for (int i = 0; i< ref.size(); i++){
      for (int j = 0; j < pool.size(); j++){
        if(ref[i].label() == pool[j].label()){
          pool.erase(pool.begin() + j);
        }
      }
    }
    return pool;
  };

auto T1_bra = in_list(connected_indices[0],T1.bra());
auto T2_bra = in_list(connected_indices[0],T2.bra());
auto T1_ket = in_list(connected_indices[0],T1.ket());
auto T2_ket = in_list(connected_indices[0],T2.ket());

if(T1_bra.first){
  upper_dominant = 1;
  assert(!T2_bra.first);
  external_T1.push_back(T1.ket()[0]);
  external_T1.push_back(T1.ket()[1]);
  external_T2.push_back(T2.bra()[0]);
  external_T2.push_back(T2.bra()[1]);
  if(nconnects == 1){
    for (int i =0; i < T1.bra().size(); i++){
      if(i != T1_bra.second){
        external_T1.push_back(T1.bra()[i]);
      }
    }
    for (int i =0; i < T2.ket().size(); i++){
      if(i != T2_ket.second){
        external_T2.push_back(T2.ket()[i]);
      }
    }
  }

}
else{
  assert(T1_ket.first);
  assert(T2_bra.first);
  upper_dominant = 0;
  external_T1.push_back(T1.bra()[0]);
  external_T1.push_back(T1.bra()[1]);
  external_T2.push_back(T2.ket()[0]);
  external_T2.push_back(T2.ket()[1]);
  if(nconnects == 1){
    for (int i =0; i < T1.ket().size(); i++){
      if(i != T1_ket.second){
        external_T1.push_back(T1.ket()[i]);
      }
    }
    for (int i =0; i < T2.bra().size(); i++){
      if(i != T2_bra.second){
        external_T2.push_back(T2.bra()[i]);
      }
    }
  }

}
unordered_internal_T2 = not_inlist(connected_indices,T2_is);
assert(unordered_internal_T2.size() < 5);

  //identify if specialized B case: internal indices on G's don't match size
  if(T1.label() == L"F" && T2.label() == L"F" && nconnects == 1 && space == IndexSpace::complete_unoccupied){
    if(unordered_internal_T2.size() == 3){
      external_T1.resize(2);
      external_T2.resize(2);
      if(unordered_internal_T2[0].label() == T2.bra()[0].label() && unordered_internal_T2[1].label() == T2.bra()[1]){
        external_T1[0] = T1.ket()[0];
        external_T1[1] = T1.ket()[1];
        external_T2[0] = T2.bra()[0];
        external_T2[1] = T2.bra()[1];
        upper_dominant = 1;
      }
      else{
        external_T1[0] = T1.bra()[0];
        external_T1[1] = T1.bra()[1];
        external_T2[0] = T2.ket()[0];
        external_T2[1] = T2.ket()[1];
        upper_dominant = 0;
      }
    }
    else if(unordered_internal_T2.size() == 2){
      external_T1.resize(2);
      external_T2.resize(2);
      external_T2[0] = unordered_internal_T2[0];
      external_T2[1] = unordered_internal_T2[1];
      if(unordered_internal_T2[0].label() == T2.bra()[0].label() && unordered_internal_T2[1].label() == T2.bra()[1]){
        external_T1[0] = T1.ket()[0];
        external_T1[1] = T1.ket()[1];
        upper_dominant = 1;
      }
      else{
        external_T1[0] = T1.bra()[0];
        external_T1[1] = T1.bra()[1];
          upper_dominant = 0;}
    }
    else{
      throw "unknown contraction scheme";}
  }
  assert(nconnects <= 2);

  if(print_){
    std::cout << "connections: " << nconnects << std::endl << "upper dominant: " << upper_dominant << std::endl;
    for (int i = 0; i < external_T1.size(); i++){
      std::wcout << " T1: " << external_T1[i].label() << std::endl;
      std::wcout << " T2: " << external_T2[i].label() << std::endl;
    }
    std::wcout << " T1: " << to_latex_align(T1.clone()) << std::endl;
    std::wcout << " T2: " << to_latex_align(T2.clone()) << std::endl;
  }
  return {nconnects, space,upper_dominant, external_T1, external_T2};
}

ExprPtr densities_to_occ(const ExprPtr& ex_){//densities become diagonal two body matricies and have added enforcement to occupied basis set.
  auto result = ex<Constant>(0);
  for (auto&& product : ex_->as<Sum>().summands()){
    auto new_product = ex<Constant>(product->as<Product>().scalar().real());
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->as<Tensor>().label() == L"\\Gamma" || factor->as<Tensor>().label() == L"\\gamma"){
        for (size_t i = 0; i < factor->as<Tensor>().bra().size(); i++) {
          new_product = new_product *
                        make_overlap(factor->as<Tensor>().bra()[i],
                                     Index::make_tmp_index(IndexSpace::instance(
                                         IndexSpace::occupied))) *
                        make_overlap(factor->as<Tensor>().ket()[i],
                                     Index::make_tmp_index(IndexSpace::instance(
                                         IndexSpace::occupied)));
        }
      }
      new_product = factor * new_product;

    }
    result = result + new_product;

  }

  FWickTheorem wick {result};
  wick.reduce(result);
  simplify (result);
  return result;
}
ExprPtr biproduct_intermediate(ExprPtr T1,ExprPtr T2){
  assert (T1->is<Tensor>());
  assert (T2->is<Tensor>());
  auto result = ex<Constant>(1);
  auto [nconnects,space,upper_dominant, external_T1, external_T2] = ncon_spa_updom_t1ext_t2ext(T1->as<Tensor>(),T2->as<Tensor>());
  if (T1->as<Tensor>().label() == L"g" || T2->as<Tensor>().label() == L"g"){
    if (nconnects == 2 && space == IndexSpace::complete_unoccupied){
      //V^pq_ij
      if(upper_dominant == 1) {
        auto V_pqij = ex<Tensor>(L"V", IDX_list{external_T2[0], external_T2[1]},
                                 IDX_list{external_T1[0], external_T1[1]});
        return V_pqij;
      }
      if(upper_dominant == 0) {
        auto V_pqij = ex<Tensor>(L"V", IDX_list{external_T1[0], external_T1[1]},
                                 IDX_list{external_T2[0], external_T2[1]});
        return V_pqij;
      }
    }
    else{
      result = T1 * T2;
    }
  }
  else{
    if (nconnects == 2 && space == IndexSpace::complete_unoccupied){
      //X^kl_ij
      if(upper_dominant == 1) {
        auto X_klij = ex<Tensor>(L"X", IDX_list{external_T2[0], external_T2[1]},
                                 IDX_list{external_T1[0], external_T1[1]});
        result = X_klij;
      }
      if(upper_dominant == 0) {
        auto X_klij = ex<Tensor>(L"X", IDX_list{external_T1[0], external_T1[1]},
                                 IDX_list{external_T2[0], external_T2[1]});
        result = X_klij;
      }
    }
    else if (nconnects == 1 && space == IndexSpace::complete_unoccupied){
      //B^kl_ij
      if(upper_dominant == 1) {
        auto B_klij = ex<Tensor>(L"B", IDX_list{external_T2[0], external_T2[1]},
                                 IDX_list{external_T1[0], external_T1[1]});
        result = B_klij;
      }
      if(upper_dominant == 0) {
        auto B_klij = ex<Tensor>(L"B", IDX_list{external_T1[0], external_T1[1]},
                                 IDX_list{external_T2[0], external_T2[1]});
        result = B_klij;
      }
    }
    //recognition of some of the S tensor can be turned on, although it should not need to be handled specially.
    /*else if(nconnects == 1 && space == IndexSpace::other_unoccupied){
      //Sk^akl_ijb
      if(upper_dominant == 1) {
        auto S_klbija = ex<Tensor>(
            L"S", IDX_list{external_T1[0], external_T1[1], external_T1[2]},
            IDX_list{external_T2[0], external_T2[1], external_T2[2]});
        result = S_klbija;
      }
      if(upper_dominant == 0) {
        auto S_klbija = ex<Tensor>(
            L"S", IDX_list{external_T2[0], external_T2[1], external_T2[2]},
            IDX_list{external_T1[0], external_T1[1], external_T1[2]});
        result = S_klbija;
      }
    }*/
    else if (nconnects == 0){
      //return original expression (no simplifications to be made)
      result = T1 * T2;
    }
    else{
      result = T1 * T2;
    }
  }
  return result;
}
Product find_f12_interms(ExprPtr ex_){
  assert(ex_->is<Product>());
  int counter = 0;
  std::vector<ExprPtr> T1_T2;
  for (auto&& factors : ex_->as<Product>().factors()){
    if(factors->as<Tensor>().label() == L"F" || factors->as<Tensor>().label() == L"g"){
      T1_T2.push_back(factors);
      counter +=1;
    }
  }
  for (auto&& factors : ex_->as<Product>().factors()){//have to loop through again unfourtunately to remove the factors that were combined.
    if(factors->as<Tensor>().label() == L"F" || factors->as<Tensor>().label() == L"g"){
      if(counter == 2){
        factors = ex<Constant>(1);
      }
    }
  }
  assert(T1_T2.size() <= 2);
  if (T1_T2.size() == 2){
    auto result = biproduct_intermediate(T1_T2[0], T1_T2[1]);
    if(result->is<Tensor>() && result->as<Tensor>().label() == L"B"){
      for (auto&& factors : ex_->as<Product>().factors()){//have to loop through again unfourtunately to remove the factors that were combined.
        if(factors->is<Tensor>() && factors->as<Tensor>().label() == L"f"){
          factors = ex<Constant>(1);
        }
      }
    }
    bool intermediate = false;
    if(result->is<Tensor>() && (result->as<Tensor>().label() == L"B" || result->as<Tensor>().label() == L"X" || result->as<Tensor>().label() == L"V")) {
      auto bra_result = ex<Constant>(1.);
      auto ket_result = ex<Constant>(1.);
      for (auto&& bra : result->as<Tensor>().bra()) {
        if (bra.space() == IndexSpace::instance(IndexSpace::occupied)) {
          bra_result = make_overlap(bra, Index::make_tmp_index(IndexSpace::instance(
                                         IndexSpace::active_occupied))) * bra_result;
        }
      }
      for (auto&& ket : result->as<Tensor>().ket()) {
        if (ket.space() == IndexSpace::instance(IndexSpace::occupied)) {
          ket_result = make_overlap(ket, Index::make_tmp_index(IndexSpace::instance(
              IndexSpace::active_occupied))) * ket_result;
        }
      }
      result = bra_result * ket_result * result;
    }

    result = result * ex_;
    FWickTheorem wick {result};
    wick.reduce(result);
    simplify(result);
    return result->as<Product>();
  }
  else{return ex_->as<Product>();}
}

//in hamiltonian based transformations, it is important to retain the original form of the hamiltonian operator. that is h^p_q E^q_p + 1/2 g^{pq}_{rs} E^{rs}_{pq}.
//to achieve this form, the tensor part of the expression must contain overlaps in place of the normal ordered operators.
//here we chose a canonical form for E^{p_7}_{p_9} and E^{p_7 p_8}_{p_9 p_10}
// this also simultaneously partitions the result into one and two body terms.
std::pair<ExprPtr,ExprPtr> fnop_to_overlap(ExprPtr exprs){
  auto one_body_result = ex<Constant>(0);
  auto two_body_result = ex<Constant>(0);
  for (auto&& product : exprs->as<Sum>().summands()){
    auto one_body_product = ex<Constant>(product->as<Product>().scalar());
    auto two_body_product = ex<Constant>(product->as<Product>().scalar());
    for (auto&& factor : product->as<Product>().factors()){
      if(factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" || factor->as<Tensor>().label() == L"a")){
        factor = tens_to_op(factor);
        if(factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().ncreators() == 1) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            one_body_product = one_body_product * o1 * o3;
            two_body_product = two_body_product * ex<Constant>(0);
          } else if (factor->as<FNOperator>().ncreators() == 2) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().annihilators()[0].index());
            auto o2 = make_overlap(
                {L"p_8"}, factor->as<FNOperator>().annihilators()[1].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().creators()[0].index(), {L"p_9"});
            auto o4 = make_overlap(
                factor->as<FNOperator>().creators()[1].index(), {L"p_10"});
            two_body_product = two_body_product * o1 * o2 * o3 * o4;
            one_body_product = one_body_product * ex<Constant>(0);
          }
        }
      }
      else{one_body_product = factor * one_body_product;
        two_body_product = factor * two_body_product;}
    }
    one_body_result = one_body_product + one_body_result;
    two_body_result = two_body_product + two_body_result;
  }
  simplify(one_body_result);
  simplify(two_body_result);
  return {one_body_result, two_body_result};
}
ExprPtr screen_F12_and_density(ExprPtr exprs){
  for (auto&& product : exprs->as<Sum>().summands()){
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F"){
        factor = screen_F_tensors(factor);
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        factor = screen_densities(factor);
      }
    }
  }
  FWickTheorem wick_f{exprs};
  wick_f.reduce(exprs);
  simplify(exprs);
  return exprs;
}

//TODO generalize for spin-orbital basis
//simplification to deal with hamiltonian based expressions. involving one body h and two body g tensors.
// not rigorous for more than 2 body operators or more than 2 density matrices whose rank must be <= 2.
// unfortunately, simplify(result) and  wick.reduce(result) will recanonicalize the indices.
// enforces the following obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}}
// should allow analysis of multiple expressions who have the same normal order operator prefactor.
std::pair<ExprPtr,ExprPtr> hamiltonian_based(ExprPtr exprs){
  exprs = remove_const(exprs);
  exprs = overlap_with_obs(exprs);
  exprs = screen_F12_and_density(exprs);
  exprs = densities_to_occ(exprs);
  for (auto&& product : exprs->as<Sum>().summands()){
    product->as<Product>() = simplification::find_f12_interms(product);
  }
  return fnop_to_overlap(exprs);
}

//TODO generalize for spin-orbital basis
//simplification to deal with fock based expressions. involving one body fock operator.
// not rigorous for more than 2 body operators or more than 2 density matrices whose rank must be <= 2.
// unfortunately, simplify(result) and wick.reduce(result) will recanonicalize the indices.
// enforces the following obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}}
// should allow analysis of multiple expressions who have the same normal order operator prefactor.
std::pair<ExprPtr,ExprPtr> fock_based (ExprPtr exprs){
  exprs = remove_const(exprs);
  exprs = overlap_with_obs(exprs);
  auto final_screen = exprs;
  //special case of F tensor screening
  for (auto&& product : final_screen->as<Sum>().summands()){
    bool one_F_screened = false; // unfourtunatley, screening both Fs at the same time gives problems to canonicalize because overlaps can lead to s^{\alpha_1}_{a'_1} s^{\alpha_1}_{a'_2}. implicitly we know that this means a'_2 == a'_1 and either is fine, but sequant does not know this yet.
    for (auto&& factor : product->as<Product>().factors()){
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F" && !one_F_screened){
        factor = screen_F_tensors(factor);
        non_canon_simplify(factor);
        one_F_screened = true;
      }
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
        factor = screen_densities(factor);
        non_canon_simplify(factor);
      }
    }
  }
  FWickTheorem wick_f{final_screen};
  wick_f.reduce(final_screen);
  non_canon_simplify(final_screen);

  //in some cases, there will now be no contributing terms left so return zero to one and two body.
if(final_screen->is<Constant>()){
  return std::pair<ExprPtr,ExprPtr> {final_screen, final_screen};
}
  final_screen = screen_F12_and_density(final_screen);
  final_screen = treat_fock(final_screen);
  // after enforcing block diagonal fock, the F12 and density tensors are screened again since some terms vanish.
  final_screen = screen_F12_and_density(final_screen);

  //enforce that densities are in the occupied space since they are only non-zero in occ
  final_screen = densities_to_occ(final_screen);
  //find the special f12 intermediates that cannot efficiently be solved directly.
  for (auto&& product : final_screen->as<Sum>().summands()){
    product->as<Product>() = simplification::find_f12_interms(product);
  }
  return fnop_to_overlap(final_screen);
  }
}
#ifndef SEQUANT_SIMPLIFICATIONS_H
#define SEQUANT_SIMPLIFICATIONS_H

#endif  // SEQUANT_SIMPLIFICATIONS_H
