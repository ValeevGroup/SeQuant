//
// Created by Conner Masteran on 7/29/21.
//
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

using namespace sequant;
// various simplifications of tensor-operator expressions based on properties of
// particular tensors. Also functionality for restricting our operators and
// densities to the orbital basis set (obs).
namespace simplification {
template <typename vec_type>
std::pair<bool, int> in_list(Index idx, vec_type ref_list) {
  bool inlist = false;
  int where_inlist = 0;
  for (int i = 0; i < ref_list.size(); i++) {
    if (idx.label() == ref_list[i].label()) {
      inlist = true;
      where_inlist = i;
    }
  }
  std::pair<bool, int> result{inlist, where_inlist};
  return result;
}
// convert a sequant::FNOperator to a sequant::tensor object
ExprPtr op_to_tens(ExprPtr ex_) {
  assert(ex_->is<FNOperator>());
  std::vector<Index> bra_indices;
  std::vector<Index> ket_indices;
  for (auto&& ann : ex_->as<FNOperator>().annihilators()) {
    bra_indices.push_back(ann.index());
  }
  for (auto&& cre : ex_->as<FNOperator>().creators()) {
    ket_indices.push_back(cre.index());
  }
  auto label =
      get_default_context().spbasis() == SPBasis::spinfree ? L"E" : L"a";
  auto result = ex<Tensor>(label, bra_indices, ket_indices);
  return result;
}

// all densities and the Hamiltonian operators are confined to a given orbital
// basis in second quantized notation. thus any index on a Normal Ordered
// operator or density must be confined to the obs.
/// TODO this dictates that the resulting hamiltonian will be in a particular
/// basis.
ExprPtr overlap_with_obs(ExprPtr ex_) {
  //std::wcout << to_latex_align(ex_,20,4) << std::endl;
  auto overlap_expr =
      ex<Constant>(0);  // enforce an overlap each E with elements from
  for (auto&& product :
       ex_->as<Sum>().summands()) {  // may be able to make_overlaps manually
                                     // and apply them to the products. simplify
                                     // may know what to do with it.
    auto new_product = ex<Constant>(product->as<Product>().scalar());
    for (int it = product->as<Product>().factors().size() - 1; it >= 0;
         it--) {  // loop backwards through the products.
      auto factor = product->as<Product>().factor(it);

      if (it == product->as<Product>().factors().size() - 1 &&
          factor->is<FNOperator>() && factor->as<FNOperator>().rank() == 3) {
        std::wstring label_1;
        std::wstring label_2;
        std::wstring label_3;
        std::wstring label_4;
        std::wstring label_5;
        std::wstring label_6;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_2 = factor->as<FNOperator>().annihilators()[1].index().label();
        label_3 = factor->as<FNOperator>().annihilators()[2].index().label();
        label_4 = factor->as<FNOperator>().creators()[0].index().label();
        label_5 = factor->as<FNOperator>().creators()[1].index().label();
        label_6 = factor->as<FNOperator>().creators()[2].index().label();
        auto o1 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_1});
        auto o2 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_2});
        auto o3 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_3});
        auto o4 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_4});
        auto o5 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_5});
        auto o6 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_6});
        new_product =
            o1 * o2 * o3 * o4 * o5 * o6 * new_product * op_to_tens(factor);
      } else if (it == product->as<Product>().factors().size() - 1 &&
                 factor->is<FNOperator>() &&
                 factor->as<FNOperator>().rank() == 2) {
        std::wstring label_1;
        std::wstring label_2;
        std::wstring label_3;
        std::wstring label_4;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_2 = factor->as<FNOperator>().annihilators()[1].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        label_4 = factor->as<FNOperator>().creators()[1].index().label();
        auto o1 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_1});
        auto o2 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_2});
        auto o3 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_3});
        auto o4 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_4});
        new_product = o1 * o2 * o3 * o4 * new_product * op_to_tens(factor);
      } else if (it == product->as<Product>().factors().size() - 1 &&
                 factor->is<FNOperator>() &&
                 factor->as<FNOperator>().rank() == 1) {
        std::wstring label_1;
        std::wstring label_3;
        label_1 = factor->as<FNOperator>().annihilators()[0].index().label();
        label_3 = factor->as<FNOperator>().creators()[0].index().label();
        auto o1 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_1});
        auto o3 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_3});
        new_product = o1 * o3 * new_product * op_to_tens(factor);
      } else if (factor->is<Tensor>() &&
                 factor->as<Tensor>().label() == L"\\Gamma" &&
                 factor->as<Tensor>().rank() == 1) {
        std::wstring label_2;
        std::wstring label_4;
        label_2 = factor->as<Tensor>().ket()[0].label();
        label_4 = factor->as<Tensor>().bra()[0].label();
        auto o1 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_2});
        auto o3 = make_overlap(
            Index{label_4},
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)));
        new_product = o1 * o3 * factor * new_product;
      } else if (factor->is<Tensor>() &&
                 factor->as<Tensor>().label() == L"\\Gamma" &&
                 factor->as<Tensor>().rank() == 2) {
        std::wstring label_1;
        std::wstring label_3;
        std::wstring label_2;
        std::wstring label_4;
        label_1 = factor->as<Tensor>().ket()[0].label();
        label_3 = factor->as<Tensor>().bra()[0].label();
        label_2 = factor->as<Tensor>().ket()[1].label();
        label_4 = factor->as<Tensor>().bra()[1].label();
        auto o1 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_1});
        auto o3 = make_overlap(
            Index{label_3},
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)));
        auto o2 = make_overlap(
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)),
            Index{label_2});
        auto o4 = make_overlap(
            Index{label_4},
            Index::make_tmp_index(IndexSpace::instance(IndexSpace::all)));
        new_product = o1 * o3 * o2 * o4 * factor * new_product;
      } else {
        new_product = new_product * factor;
      }
    }
    new_product = new_product;
    overlap_expr = overlap_expr + new_product;
  }

  FWickTheorem wick{overlap_expr};
  // std::wcout << to_latex_align(overlap_expr,20,2) << std::endl;
  wick.reduce(overlap_expr);
  non_canon_simplify(overlap_expr);
  // std::wcout << to_latex_align(overlap_expr,20,2) << std::endl;
  return overlap_expr;
}

using IDX_list = std::initializer_list<Index>;
// in various transformation methods it seems as if the constants are removed or
// treated separatly from the main transformed hamiltonian expression.
ExprPtr remove_const(const ExprPtr ex_) {
  auto new_expression = ex<Constant>(0);
  if (ex_->is<Sum>()) {
    for (auto&& product : ex_->as<Sum>().summands()) {
      bool has_fnop = false;
      for (auto&& factor : product->as<Product>().factors()) {
        if (factor->is<FNOperator>()) {
          has_fnop = true;
        }
      }
      if (has_fnop) {
        new_expression = new_expression + product;
      }
    }
  }
  non_canon_simplify(new_expression);
  return new_expression;
}

// params ex_ : a product to replace indices on.
//  og: original index in the product to be replaced
//  newer: the new index which replaces the original index.
Product replace_idx(ExprPtr ex_, Index og, Index newer) {
  assert(ex_->is<Product>());
  auto constant = ex_->as<Product>().scalar();
  auto new_product = ex<Constant>(1);
  for (auto&& factor : ex_->as<Product>().factors()) {
    if (factor->is<Tensor>()) {
      std::vector<Index> new_bras;
      for (auto&& bra : factor->as<Tensor>().bra()) {
        if (bra.label() == og.label()) {
          new_bras.push_back(newer);
        } else {
          new_bras.push_back(bra);
        }
      }
      std::vector<Index> new_kets;
      for (auto&& ket : factor->as<Tensor>().ket()) {
        if (ket.label() == og.label()) {
          new_kets.push_back(newer);
        } else {
          new_kets.push_back(ket);
        }
      }
      auto new_tensor =
          ex<Tensor>(factor->as<Tensor>().label(), new_bras, new_kets);
      new_product = new_tensor * new_product;
    }
    if (factor->is<FNOperator>()) {
      std::vector<Index> new_cres;
      for (auto&& cre : factor->as<FNOperator>().creators()) {
        if (cre.index().label() == og.label()) {
          new_cres.push_back(newer);
        } else {
          new_cres.push_back(cre.index());
        }
      }
      std::vector<Index> new_anns;
      for (auto&& ann : factor->as<FNOperator>().annihilators()) {
        if (ann.index().label() == og.label()) {
          new_anns.push_back(newer);
        } else {
          new_anns.push_back(ann.index());
        }
      }
      if (factor->as<FNOperator>().ncreators() == 1) {
        auto o1 = make_overlap({L"p_7"}, new_anns[0]);
        auto o3 = make_overlap(new_cres[0], {L"p_9"});
        new_product = new_product * o1 * o3;
      } else if (factor->as<FNOperator>().ncreators() == 2) {
        auto o1 = make_overlap({L"p_7"}, new_anns[0]);
        auto o2 = make_overlap({L"p_8"}, new_anns[1]);
        auto o3 = make_overlap(new_cres[0], {L"p_9"});
        auto o4 = make_overlap(new_cres[1], {L"p_10"});
        new_product = new_product * o1 * o2 * o3 * o4;
      } else {
        throw "does not handle size > 2";
      }
      // auto new_op = ex<FNOperator>(new_cres,new_anns);
      //  new_product = new_product * new_op;
    }
  }
  auto result = (ex<Constant>(constant) * new_product);
  return result->as<Product>();
}

// convert a sequant::Tensor to a sequant::FNOperator
ExprPtr tens_to_op(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  auto result =
      ex<FNOperator>(ex_->as<Tensor>().ket(), ex_->as<Tensor>().bra());
  return result;
}
// F tensors must contain indices in the bra with space > all. this
// includes complete, completeunoccupied, and inactiveunoccupied. and if one of
// the particle indices is connected to the obs virtual space, then the other
// must be from the CABS set. i.e. if G^{a \beta}_{ij} -> G^{a a'}_{ij}
ExprPtr screen_F_tensors(ExprPtr ex_, int ansatz = 2) {
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().label() == L"F");
  auto overlap = ex<Constant>(1);
  bool good = false;
  bool bra_good = false;
  if (ansatz == 2) {
    for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++) {
      auto bra = ex_->as<Tensor>().bra()[i];
      if (bra.space().type() == IndexSpace::complete ||
          bra.space().type() == IndexSpace::complete_unoccupied) {
        good = true;
        bra_good = true;
      } else if (bra.space().type() == IndexSpace::complete ||
                 bra.space().type() == IndexSpace::complete_unoccupied ||
                 bra.space().type() == IndexSpace::other_unoccupied) {
        good = true;
      }
    }

    for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++) {
      auto bra = ex_->as<Tensor>().bra()[i];
      if ((bra.space().type() == IndexSpace::unoccupied ||
           bra.space().type() == IndexSpace::all) &&
          bra_good &&
          (bra.space().type() !=
           IndexSpace::complete_unoccupied)) {  // if one of the upper indices
                                                // is explicitly outside of
                                                // CABs, create an overlap with
                                                // the other index and the CABs
                                                // space.
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
    for (int j = 0; j < ex_->as<Tensor>().ket().size(); j++) {
      auto ket = ex_->as<Tensor>().ket()[j];
      if (ket.space().type() == IndexSpace::complete ||
          ket.space().type() == IndexSpace::complete_unoccupied) {
        good = true;
        ket_good = true;
      } else if (ket.space().type() == IndexSpace::complete ||
                 ket.space().type() == IndexSpace::complete_unoccupied ||
                 ket.space().type() == IndexSpace::other_unoccupied) {
        good = true;
      }
    }
    for (int j = 0; j < ex_->as<Tensor>().ket().size(); j++) {
      auto ket = ex_->as<Tensor>().ket()[j];
      if ((ket.space().type() == IndexSpace::unoccupied ||
           ket.space().type() == IndexSpace::all) &&
          ket_good && (ket.space().type() != IndexSpace::complete_unoccupied)) {
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
    if (good) {
      return ex_ * overlap;
    } else {
      return ex<Constant>(0);
    }
  } else if (ansatz == 1) {
    bool non_zero = false;
    bool bra_proj_space =
        false;  // perhaps a better way would be to create a child class of
                // tensor for G tensor which can keep track of geminal
                // generating and projector at construction.
    for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++) {
      auto bra = ex_->as<Tensor>().bra()[i];
      if (bra.space().type() == IndexSpace::complete ||
          bra.space().type() == IndexSpace::complete_unoccupied ||
          bra.space().type() == IndexSpace::other_unoccupied) {
        bra_proj_space = true;
        non_zero = true;
      }
    }
    if (bra_proj_space) {
      for (int i = 0; i < ex_->as<Tensor>().bra().size(); i++) {
        auto bra = ex_->as<Tensor>().bra()[i];
        auto overlap1 = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                         IndexSpace::other_unoccupied)),
                                     bra);
        ex_ = overlap1 * ex_;
      }
    }
    bool ket_proj_space =
        false;  // perhaps a better way would be to create a child class of
                // tensor for G tensor which can keep track of geminal
                // generating and projector at construction.
    for (int i = 0; i < ex_->as<Tensor>().ket().size(); i++) {
      auto ket = ex_->as<Tensor>().ket()[i];
      if (ket.space().type() == IndexSpace::complete ||
          ket.space().type() == IndexSpace::complete_unoccupied ||
          ket.space().type() == IndexSpace::other_unoccupied) {
        ket_proj_space = true;
        non_zero = true;
      }
    }
    if (ket_proj_space) {
      for (int i = 0; i < ex_->as<Tensor>().ket().size(); i++) {
        auto ket = ex_->as<Tensor>().ket()[i];
        auto overlap1 = make_overlap(Index::make_tmp_index(IndexSpace::instance(
                                         IndexSpace::other_unoccupied)),
                                     ket);
        ex_ = overlap1 * ex_;
      }
    }
    if (non_zero) {
      return ex_;
    } else {
      return ex<Constant>(0.0);
    }
  }
}

ExprPtr screen_density(
    ExprPtr ex_) {  // densities probably should be non-zero if each index has a
                    // chance to be occupied, in other words, screen out
                    // densities containing unoccupied labels.
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().label() == L"\\Gamma" ||
         ex_->as<Tensor>().label() == L"\\gamma");
  bool occ_space = true;
  for (auto&& bra : ex_->as<Tensor>().bra()) {
    if (bra.space().type() == IndexSpace::unoccupied ||
        bra.space().type() == IndexSpace::complete_unoccupied) {
      occ_space = false;
    }
  }
  for (auto&& ket : ex_->as<Tensor>().ket()) {
    if (ket.space().type() == IndexSpace::unoccupied ||
        ket.space().type() == IndexSpace::complete_unoccupied) {
      occ_space = false;
    }
  }
  if (occ_space) {
    return ex_;
  } else {
    return ex<Constant>(0);
  }
};
ExprPtr screen_densities(ExprPtr ex_) {
  if (ex_->is<Sum>()) {
    for (auto&& product : ex_->as<Sum>().summands()) {
      for (auto&& factor : product->as<Product>()) {
        if(factor->is<Tensor>()) {
          if (factor->as<Tensor>().label() == L"\\Gamma" ||
              factor->as<Tensor>().label() == L"\\gamma") {
            factor = screen_density(factor);
          }
        }
        else{std::wcout <<"problematic factor" << to_latex_align(factor);}
      }
    }
    return ex_;
  } else if (ex_->is<Product>()) {
    for (auto&& factor : ex_->as<Product>()) {
      if (factor->as<Tensor>().label() == L"\\Gamma" ||
          factor->as<Tensor>().label() == L"\\gamma") {
        factor = screen_density(factor);
      }
    }
    return ex_;
  } else if (ex_->is<Tensor>()) {
    return screen_density(ex_);
  } else {
    throw "unsupported operation for screening densities";
  }
}

// based on Brillouin's Theory, the fock matrix should be block diagonal.
//  generalized says that complete unoccupied might be non-zero with obs unocc,
//  but zero with occ and complete unocc. lets assume normal Brillouin's Theory.
auto treat_fock(ExprPtr ex_) {
  auto new_ex_ = ex<Constant>(0);
  for (auto&& product : ex_->as<Sum>().summands()) {
    double real = product->as<Product>().scalar().real();
    auto new_product = ex<Constant>(real);
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"f") {
        // TODO do not assume EBC
        auto space = intersection(factor->as<Tensor>().bra()[0].space(),
                                  factor->as<Tensor>().ket()[0].space());
        if (space.type().none()) {
          new_product = ex<Constant>(0) * new_product;
        } else {
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
      } else
        new_product = new_product * factor;
    }
    new_ex_ = new_ex_ + new_product;
  }
  FWickTheorem wick{new_ex_};
  wick.reduce(new_ex_);
  simplify(new_ex_);

  return new_ex_;
}
// to Identify the relavant F12 intermediates, the number of connections,the
// connected space, and the resulting ket() and bra() of the intermediate tensor
// are needed.
std::tuple<int, IndexSpace::Type, std::vector<Index>, std::vector<Index>, bool>
ncon_spa_extket_extbra(Tensor T1, Tensor T2, bool print_ = false) {
  // connected space. in each example in f12, the connected space is the same
  // between two tensors.
  auto space = IndexSpace::occupied;  // just a default used for construction.
  // depreciated should be a braket function somewhere in Tensor.
  std::vector<Index> T1_is;
  std::vector<Index> T2_is;
  // ordered list of ket and bra indices which construct the resulting
  // intermediate.
  std::vector<Index> external_ket;
  std::vector<Index> external_bra;
  // do the external ket indices correspond to T1?
  //  only need for intermediates and only works for V or X.
  bool T1_ket;
  // unique list of connected indices. list is searched often to see if a given
  // index is connected.
  std::vector<Index> connected_indices;

  int nconnects = 0;
  for (auto&& bra : T1.bra()) {
    T1_is.push_back(bra);
  }
  for (auto&& ket : T1.ket()) {
    T1_is.push_back(ket);
  }
  for (auto&& bra : T2.bra()) {
    T2_is.push_back(bra);
  }
  for (auto&& ket : T2.ket()) {
    T2_is.push_back(ket);
  }
  for (int i1 = 0; i1 < T1_is.size(); i1++) {
    for (int i2 = 0; i2 < T2_is.size(); i2++) {
      if (T1_is[i1].label() == T2_is[i2].label()) {
        nconnects += 1;
        space = T1_is[i1].space().type();
        connected_indices.push_back(T1_is[i1]);
      }
    }
  }
  if (nconnects == 0) {
    external_ket = T1_is;
    external_bra = T2_is;
    std::tuple zero{nconnects, space, external_ket, external_bra, T1_ket};
    return zero;
  }
  // which indices in the T1  bra and ket are connected? what corresponding ket
  // or bra is it connected to in T2?
  //  place the shared unconnected particle indices in the external bra/ket
  //  list.
  for (int i = 0; i < T1.bra().size(); i++) {
    // is the bra T1 index a connected index?
    if (in_list(T1.bra()[i], connected_indices).first) {
      T1_ket = true;
      for (int j = 0; j < T2.ket().size(); j++) {
        if (T2.ket()[j].label() == T1.bra()[i].label()) {
          external_ket.push_back(T1.ket()[i]);
          external_bra.push_back(T2.bra()[j]);
        }
      }
    }
    // is the ket T1 index a connected index?
    else if (in_list(T1.ket()[i], connected_indices).first) {
      T1_ket = false;
      for (int j = 0; j < T2.ket().size(); j++) {
        if (T2.bra()[j].label() == T1.ket()[i].label()) {
          external_ket.push_back(T2.ket()[i]);
          external_bra.push_back(T1.bra()[j]);
        }
      }
    }
    // neither are connected.
    else {
      external_ket.push_back(T1.ket()[i]);
      external_bra.push_back(T1.bra()[i]);
    }
  }
  // loop through T2
  for (int i = 0; i < T2.ket().size(); i++) {
    // if the ket index is connected, do nothing because the external index is
    // already accounted for
    if (in_list(T2.ket()[i], connected_indices).first ||
        in_list(T2.ket()[i], external_ket).first) {
    }
    // if the bra index is connected, do nothing because the external index is
    // already accounted for
    else if (in_list(T2.bra()[i], connected_indices).first ||
             in_list(T2.bra()[i], external_bra).first) {
    }
    // if niether the bra or the ket are connected or made the external lists by
    // now, add them.
    else {
      external_ket.push_back(T2.ket()[i]);
      external_bra.push_back(T2.bra()[i]);
    }
  }

  // identify if specialized B case: 1 connect between complete unoccupied.
  // TODO This logic is not entirely general! exact, general identification of B
  // will require full product screening. and perhaps screening the sum!
  //  the above logic works in the fock approximation of the double commutator.
  //  consider also that a product can only identify a single term in the B, in
  //  fact when the second particle on each tensor is connected, does not show
  //  up
  // since the left and right terms are equivalent for column symmetric G ops.
  if (T1.label() == L"F" && T2.label() == L"F" && nconnects == 1 &&
      space == IndexSpace::complete_unoccupied) {
    external_ket.resize(0);
    external_bra.resize(0);
    bool bra_connected = false;
    bool ket_connected = false;
    for (int i = 0; i < T1.bra().size(); i++) {
      if (in_list(T1.bra()[i], connected_indices).first) {
        bra_connected = true;
      }
    }
    for (int j = 0; j < T1.ket().size(); j++) {
      if (in_list(T1.ket()[j], connected_indices).first) {
        ket_connected = true;
      }
    }
    assert(ket_connected || bra_connected);
    if (ket_connected) {
      external_ket.push_back(T2.ket()[0]);
      external_ket.push_back(T2.ket()[1]);
      external_bra.push_back(T1.bra()[0]);
      external_bra.push_back(T1.bra()[1]);

    } else {
      external_ket.push_back(T1.ket()[0]);
      external_ket.push_back(T1.ket()[1]);
      external_bra.push_back(T2.bra()[0]);
      external_bra.push_back(T2.bra()[1]);
    }
  }
  assert(nconnects <= 2);
  return {nconnects, space, external_ket, external_bra, T1_ket};
}
// densities are enforced to map obs -> occ since the obs includes frozen core
// orbitals.
ExprPtr densities_to_occ(const ExprPtr& ex_) {
  auto result = ex<Constant>(0);
  for (auto&& product : ex_->as<Sum>().summands()) {
    auto new_product = ex<Constant>(product->as<Product>().scalar());
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() &&
          (factor->as<Tensor>().label() == L"\\Gamma" ||
           factor->as<Tensor>().label() == L"\\gamma")) {
        for (size_t i = 0; i < factor->as<Tensor>().bra().size(); i++) {
          if (factor->as<Tensor>().bra()[i].space() == IndexSpace::all) {
            new_product =
                new_product *
                make_overlap(factor->as<Tensor>().bra()[i],
                             Index::make_tmp_index(
                                 IndexSpace::instance(IndexSpace::occupied)));
          }
          if (factor->as<Tensor>().ket()[i].space() == IndexSpace::all) {
            new_product =
                new_product *
                make_overlap(factor->as<Tensor>().ket()[i],
                             Index::make_tmp_index(
                                 IndexSpace::instance(IndexSpace::occupied)));
          }
        }
      }
      new_product = factor * new_product;
    }
    result = result + new_product;
  }

  FWickTheorem wick{result};
  wick.reduce(result);
  simplify(result);
  return result;
}

// constructs a biproduct intermediate tensor from a given two tensors in an
// expression.
ExprPtr biproduct_intermediate(ExprPtr T1, ExprPtr T2) {
  assert(T1->is<Tensor>());
  assert(T2->is<Tensor>());
  auto result = ex<Constant>(1);
  auto [nconnects, space, external_ket, external_bra, T1_ket] =
      ncon_spa_extket_extbra(T1->as<Tensor>(), T2->as<Tensor>());
  if (T1->as<Tensor>().label() == L"g" || T2->as<Tensor>().label() == L"g") {
    // V^pq_ij
    // intermediate decomposition handled by SeQuant so space labels can be
    // properly handled
    if (nconnects == 2 && space == IndexSpace::complete_unoccupied) {
      if (T1_ket) {
        auto GR_ijpq =
            ex<Tensor>(L"GR", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{external_ket[0], external_ket[1]});
        auto F_ijrs =
            ex<Tensor>(L"F", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{L"p_11", L"p_12"});
        auto g_rspq = ex<Tensor>(L"g", IDX_list{L"p_11", L"p_12"},
                                 IDX_list{external_ket[0], external_ket[1]});
        auto F_ijmc =
            ex<Tensor>(L"F", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{L"m_6", L"α'_4"});
        auto g_mcpq = ex<Tensor>(L"g", IDX_list{L"m_6", L"α'_4"},
                                 IDX_list{external_ket[0], external_ket[1]});
        auto F_jicm =
            ex<Tensor>(L"F", IDX_list{external_bra[1], external_bra[0]},
                       IDX_list{L"α'_4", L"m_6"});
        auto g_cmqp = ex<Tensor>(L"g", IDX_list{L"α'_4", L"m_6"},
                                 IDX_list{external_ket[1], external_ket[0]});

        auto V = GR_ijpq - F_ijrs * g_rspq - F_ijmc * g_mcpq - F_jicm * g_cmqp;
        simplify(V);
        return V;
      } else {
        auto GR_pqij =
            ex<Tensor>(L"GR", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{external_ket[0], external_ket[1]});
        auto F_rsij = ex<Tensor>(L"F", IDX_list{L"p_11", L"p_12"},
                                 IDX_list{external_ket[0], external_ket[1]});
        auto g_pqrs =
            ex<Tensor>(L"g", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{L"p_11", L"p_12"});
        auto F_mcij = ex<Tensor>(L"F", IDX_list{L"m_6", L"α'_4"},
                                 IDX_list{external_ket[0], external_ket[1]});
        auto g_pqmc =
            ex<Tensor>(L"g", IDX_list{external_bra[0], external_bra[1]},
                       IDX_list{L"m_6", L"α'_4"});
        auto F_cmji = ex<Tensor>(L"F", IDX_list{L"α'_4", L"m_6"},
                                 IDX_list{external_ket[1], external_ket[0]});
        auto g_qpcm =
            ex<Tensor>(L"g", IDX_list{external_bra[1], external_bra[0]},
                       IDX_list{L"α'_4", L"m_6"});

        auto V = GR_pqij - F_rsij * g_pqrs - F_mcij * g_pqmc - F_cmji * g_qpcm;
        simplify(V);
        return V;
      }
    } else {
      result = T1 * T2;
    }
  } else {
    if (nconnects == 2 && space == IndexSpace::complete_unoccupied) {
      // X^kl_ij
      auto X_klij = ex<Tensor>(L"X", IDX_list{external_bra[0], external_bra[1]},
                               IDX_list{external_ket[0], external_ket[1]});
      result = X_klij;

    } else if (nconnects == 1 && space == IndexSpace::complete_unoccupied) {
      // B^kl_ij
      auto B_klij = ex<Tensor>(L"B", IDX_list{external_bra[0], external_bra[1]},
                               IDX_list{external_ket[0], external_ket[1]});
      result = B_klij;
    } else if (nconnects == 0) {
      // return original expression (no simplifications to be made)
      result = T1 * T2;
    } else {
      result = T1 * T2;
    }
  }
  return result;
}
// identify F12 intermediates
// intermediates we generate contain either 2 F or g tensors.
// those expressions are biproduct intermediate for further screening.
// special case of B intermediate is handled by an additional check for the fock
// operator f
ExprPtr find_F12_interms(ExprPtr ex_) {
  assert(ex_->is<Product>());
  int counter = 0;
  std::vector<ExprPtr> T1_T2;
  for (auto&& factors : ex_->as<Product>().factors()) {
    assert(!factors->is<FNOperator>());
    if (factors->as<Tensor>().label() == L"F" ||
        factors->as<Tensor>().label() == L"g") {
      T1_T2.push_back(factors);
      counter += 1;
    }
  }
  for (auto&& factors :
       ex_->as<Product>()
           .factors()) {  // have to loop through again unfourtunately to remove
                          // the factors that were combined.
    if (factors->as<Tensor>().label() == L"F" ||
        factors->as<Tensor>().label() == L"g") {
      if (counter == 2) {
        factors = ex<Constant>(1);
      }
    }
  }
  assert(T1_T2.size() <= 2);
  if (T1_T2.size() == 2) {
    assert(counter == 2);
    auto result = biproduct_intermediate(T1_T2[0], T1_T2[1]);
    if (result->is<Tensor>() && result->as<Tensor>().label() == L"B") {
      for (auto&& factors :
           ex_->as<Product>()
               .factors()) {  // have to find fock matrix and remove. factor 1/2
                              // because a product only finds 1/2 of the B
                              // tensor, a sum of two products.
        if (factors->is<Tensor>() && factors->as<Tensor>().label() == L"f") {
          factors = ex<Constant>(1.);
        }
      }
    }

    result = result * ex_;
    non_canon_simplify(result);
    return result;
  }
  return ex_;
}

// in hamiltonian based transformations, it is important to retain the original
// form of the hamiltonian operator. that is h^p_q E^q_p + 1/2 g^{pq}_{rs}
// E^{rs}_{pq}. to achieve this form, the tensor part of the expression must
// contain overlaps in place of the normal ordered operators. here we chose a
// canonical form for E^{p_7}_{p_9} and E^{p_7 p_8}_{p_9 p_10} as the external indicies
//  this also simultaneously partitions the result into one and two body terms.
std::pair<ExprPtr, ExprPtr> fnop_to_overlap(ExprPtr exprs) {
  auto one_body_result = ex<Constant>(0);
  auto two_body_result = ex<Constant>(0);
  for (auto&& product : exprs->as<Sum>().summands()) {
    auto one_body_product = ex<Constant>(product->as<Product>().scalar());
    auto two_body_product = ex<Constant>(product->as<Product>().scalar());
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" ||
                                   factor->as<Tensor>().label() == L"a")) {
        factor = tens_to_op(factor);
        if (factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().ncreators() == 1) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().creators()[0].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().annihilators()[0].index(), {L"p_9"});
            one_body_product = one_body_product * o1 * o3;
            two_body_product = two_body_product * ex<Constant>(0);
          } else if (factor->as<FNOperator>().ncreators() == 2) {
            auto o1 = make_overlap(
                {L"p_7"}, factor->as<FNOperator>().creators()[0].index());
            auto o2 = make_overlap(
                {L"p_8"}, factor->as<FNOperator>().creators()[1].index());
            auto o3 = make_overlap(
                factor->as<FNOperator>().annihilators()[0].index(), {L"p_9"});
            auto o4 = make_overlap(
                factor->as<FNOperator>().annihilators()[1].index(), {L"p_10"});
            two_body_product = two_body_product * o1 * o2 * o3 * o4;
            one_body_product = one_body_product * ex<Constant>(0);
          }
        }
      } else {
        one_body_product = factor * one_body_product;
        two_body_product = factor * two_body_product;
      }
    }
    one_body_result = one_body_product + one_body_result;
    two_body_result = two_body_product + two_body_result;
  }
  non_canon_simplify(one_body_result);
  non_canon_simplify(two_body_result);
  return {one_body_result, two_body_result};
}

std::pair<bool, ExprPtr> contains_tens(ExprPtr ex_, std::wstring label) {
  if (!ex_->is<Product>()) {
    ex_ = ex<Constant>(1.) * ex_;
  }
  int it = 0;
  for (auto&& factor : ex_->as<Product>().factors()) {
    if (factor->is<Tensor>() && factor->as<Tensor>().label() == label) {
      return {true, factor};
    }
    it++;
  }
  return {false, ex_};
}

// TODO this should be a generalized procedure since the screening process is
// different for each number of F tensors.
//  I suppose generally, this should be a product level screening, which first
//  finds the number of F tensors and then picks the correct screening method.
// re-implimentation as a recursive function which gets called every time a
// delta is found, simplifies/reduces the product and returns. products are
// const and two deltas acting on the same index makes this difficult. logically
// the product needs to update within its own loop, but it cannot.
// Alternatively, two delta's to the same index need to occur in the same
// product, but that breaks things. work around. make a copy of product which
// can be modified? break out of product loop?
ExprPtr screen_F12_proj(ExprPtr exprs, int ansatz = 2) {
  if (exprs->is<Sum>()) {
    auto return_sum = ex<Constant>(0);
    for (auto&& product : exprs->as<Sum>().summands()) {
      auto new_product = ex<Constant>(product->as<Product>().scalar().real());
      for (auto&& factor : product->as<Product>().factors()) {
        auto temp_factor = ex<Constant>(1.);
        if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F") {
          temp_factor = screen_F_tensors(
              factor,
              ansatz);  // screen F tensors should just provide the delta.
        } else {
          temp_factor = factor;
        }
        auto product_clone = product->clone();
        if (contains_tens(temp_factor, L"s").first) {
          product_clone =
              product_clone * contains_tens(temp_factor, L"s").second;
          FWickTheorem wick_f{product_clone};
          wick_f.reduce(product_clone);
          simplify(product_clone);
          product_clone = screen_F12_proj(product_clone, ansatz);
          return_sum = product_clone + return_sum;
          new_product = ex<Constant>(0.);
          break;
        }
        new_product = new_product * temp_factor;
        non_canon_simplify(new_product);
      }
      return_sum = new_product + return_sum;
    }
    non_canon_simplify(return_sum);
    return return_sum;
  } else if (exprs->is<Product>()) {
    auto new_product = ex<Constant>(exprs->as<Product>().scalar());
    for (auto&& factor : exprs->as<Product>().factors()) {
      auto temp_factor = ex<Constant>(1.);
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F") {
        temp_factor = screen_F_tensors(
            factor, ansatz);  // screen F tensors should just provide the delta.
      } else {
        temp_factor = factor;
      }
      auto product_clone = exprs->clone();
      if (contains_tens(temp_factor, L"s").first) {
        product_clone = product_clone * contains_tens(temp_factor, L"s").second;
        FWickTheorem wick_f{product_clone};
        wick_f.reduce(product_clone);
        simplify(product_clone);
        product_clone = screen_F12_proj(product_clone, ansatz);
        new_product = product_clone;
        break;
      }
      new_product = new_product * temp_factor;
    }
    return new_product;
  } else
    return exprs;
}

ExprPtr FNOPs_to_tens(ExprPtr ex_) {
  if (ex_->is<Sum>()) {
    auto new_sum = ex<Constant>(0);
    for (auto&& product : ex_->as<Sum>().summands()) {
      auto new_product = ex<Constant>(product->as<Product>().scalar());
      for (auto factor : product->as<Product>().factors()) {
        auto new_factor = ex<Constant>(0);
        if (factor->is<FNOperator>()) {
          new_factor = op_to_tens(factor) + new_factor;
          assert(!new_factor->is<FNOperator>());
        } else {
          new_factor = factor + new_factor;
        }
        new_product = new_product * new_factor;
      }
      new_sum = new_product + new_sum;
    }
    non_canon_simplify(new_sum);
    return new_sum;
  } else if (ex_->is<Product>()) {
    for (auto&& factor : ex_->as<Product>().factors()) {
      if (factor->is<FNOperator>()) {
        factor = op_to_tens(factor);
      }
    }
  } else if (ex_->is<FNOperator>()) {
    ex_ = op_to_tens(ex_);
  } else {
    return ex_;
  }
  return ex_;
}
ExprPtr tens_to_FNOps(ExprPtr ex_) {
  if (ex_->is<Sum>()) {
    auto new_sum = ex<Constant>(0);
    for (auto&& product : ex_->as<Sum>().summands()) {
      auto new_product = ex<Constant>(product->as<Product>().scalar());
      for (auto factor : product->as<Product>().factors()) {
        auto new_factor = ex<Constant>(0);
        if (factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" ||
                                     factor->as<Tensor>().label() == L"a")) {
          new_factor = tens_to_op(factor);
        } else {
          new_factor = factor;
        }
        new_product = new_factor * new_product;
      }
      new_sum = new_product + new_sum;
    }
    non_canon_simplify(new_sum);
    return new_sum;
  } else if (ex_->is<Product>()) {
    for (auto&& factor : ex_->as<Product>().factors()) {
      if (factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" ||
                                   factor->as<Tensor>().label() == L"a")) {
        factor = tens_to_op(factor);
      }
    }
  } else if (ex_->is<Tensor>() && (ex_->as<Tensor>().label() == L"E" ||
                                   ex_->as<Tensor>().label() == L"a")) {
    ex_ = tens_to_op(ex_);
  } else {
    return ex_;
  }
  return ex_;
}

// split F12 operator into its 2 components seen in eq 11. of Chem. Phys. 136,
// 084107 (2012).
//  neccessary in some cases where particles get excited from different spaces.
ExprPtr split_F12(ExprPtr exprs) {
  assert(exprs->is<Tensor>());
  assert(exprs->as<Tensor>().label() == L"F");
  auto result = ex<Constant>(0);
  if ((exprs->as<Tensor>().const_braket()[2].space() ==
           sequant::IndexSpace::complete_unoccupied ||
       exprs->as<Tensor>().const_braket()[2].space() ==
           sequant::IndexSpace::other_unoccupied) ||
      exprs->as<Tensor>().const_braket()[3].space() ==
          sequant::IndexSpace::complete_unoccupied ||
      exprs->as<Tensor>().const_braket()[3].space() ==
          sequant::IndexSpace::other_unoccupied) {
    auto T1 =
        ex<Constant>(3. / 8) *
        ex<Tensor>(L"F",
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[0],
                                      exprs->as<Tensor>().const_braket()[1]},
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[2],
                                      exprs->as<Tensor>().const_braket()[3]});
    auto T2 =
        ex<Constant>(1. / 8) *
        ex<Tensor>(L"F",
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[1],
                                      exprs->as<Tensor>().const_braket()[0]},
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[2],
                                      exprs->as<Tensor>().const_braket()[3]});
    result = T1 + T2;
    return result;
  } else {  // otherwise the geminal generating space must be in the upper
            // indices. so include exchange for those.
    assert((exprs->as<Tensor>().const_braket()[0].space() ==
                sequant::IndexSpace::complete_unoccupied ||
            exprs->as<Tensor>().const_braket()[0].space() ==
                sequant::IndexSpace::other_unoccupied) ||
           (exprs->as<Tensor>().const_braket()[1].space() ==
                sequant::IndexSpace::complete_unoccupied ||
            exprs->as<Tensor>().const_braket()[1].space() ==
                sequant::IndexSpace::other_unoccupied));
    auto T1 =
        ex<Constant>(3. / 8) *
        ex<Tensor>(L"F",
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[0],
                                      exprs->as<Tensor>().const_braket()[1]},
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[2],
                                      exprs->as<Tensor>().const_braket()[3]});
    auto T2 =
        ex<Constant>(1. / 8) *
        ex<Tensor>(L"F",
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[0],
                                      exprs->as<Tensor>().const_braket()[1]},
                   std::vector<Index>{exprs->as<Tensor>().const_braket()[3],
                                      exprs->as<Tensor>().const_braket()[2]});
    result = T1 + T2;
    return result;
  }
}

ExprPtr partition_F12(ExprPtr exprs) {
  if (!exprs->is<Sum>()) {
    return exprs;
  }

  for (auto&& product : exprs->as<Sum>().summands()) {
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"F") {
        factor = split_F12(factor);
      }
    }
  }
  simplify(exprs);
  return (exprs);
}

// TODO generalize for spin-orbital basis
// simplification to deal with hamiltonian based expressions. involving one body
// h and two body g tensors.
//  not rigorous for more than 2 body operators or more than 2 density matrices
//  whose rank must be <= 2. unfortunately, simplify(result) and
//  wick.reduce(result) will recanonicalize the indices. enforces the following
//  obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}} should
//  allow analysis of multiple expressions who have the same normal order
//  operator prefactor.
std::pair<ExprPtr, ExprPtr> hamiltonian_based_projector_2(ExprPtr exprs) {
  exprs = FNOPs_to_tens(exprs);
  non_canon_simplify(exprs);
  exprs = screen_densities(exprs);
  non_canon_simplify(exprs);
  exprs = screen_F12_proj(exprs, 2);
  //simplify(exprs);
  //exprs = partition_F12(exprs);
  non_canon_simplify(exprs);
  auto exprs_intmed = ex<Constant>(0.0);
  for (auto&& product : exprs->as<Sum>().summands()) {
    auto new_product = simplification::find_F12_interms(product);
    exprs_intmed = new_product + exprs_intmed;
  }
  non_canon_simplify(exprs_intmed);
  return fnop_to_overlap(exprs_intmed);
}

// here G can only have projection to the alpha and Beta space otherwise
// projector constructs it to be be zero.
std::pair<ExprPtr, ExprPtr> hamiltonian_based_projector_1(ExprPtr exprs) {
  exprs = FNOPs_to_tens(exprs);
  simplify(exprs);
  exprs = partition_F12(exprs);
  simplify(exprs);
  exprs = screen_F12_proj(exprs, 1);
  simplify(exprs);
  auto exprs_intmed = ex<Constant>(0.0);
  for (auto&& product : exprs->as<Sum>().summands()) {
    auto new_product = simplification::find_F12_interms(product);
    exprs_intmed = new_product + exprs_intmed;
  }
  non_canon_simplify(exprs_intmed);
  return fnop_to_overlap(exprs_intmed);
}
// G can only project to alpha and Beta space. still need to use fock based
// expression.
std::pair<ExprPtr, ExprPtr> fock_based_projector_1(ExprPtr exprs) {
  exprs = FNOPs_to_tens(exprs);
  simplify(exprs);
  if (exprs->is<Constant>()) {
    return std::pair<ExprPtr, ExprPtr>{exprs, exprs};
  }
  exprs = partition_F12(exprs);
  auto final_screen = exprs;
  simplify(final_screen);
  // in some cases, there will now be no contributing terms left so return zero
  // to one and two body.
  if (final_screen->is<Constant>()) {
    return std::pair<ExprPtr, ExprPtr>{final_screen, final_screen};
  }
  simplify(final_screen);
  // find the special f12 intermediates that cannot efficiently be solved
  // directly. This seems to work already for the general case!
  auto last_screen = ex<Constant>(0.0);
  for (auto&& product : final_screen->as<Sum>().summands()) {
    auto new_product = simplification::find_F12_interms(product);
    last_screen = last_screen + new_product;
  }
  simplify(last_screen);
  return fnop_to_overlap(last_screen);
}
// TODO generalize for spin-orbital basis
// simplification to deal with fock based expressions. involving one body fock
// operator.
//  not rigorous for more than 2 body operators or more than 2 density matrices
//  whose rank must be <= 2. unfortunately, simplify(result) and
//  wick.reduce(result) will recanonicalize the indices. enforces the following
//  obs convention. E^{p_7}_{p_9} and E^{{p_7}{p_8}}_{{p_9}{p_{10}}} should
//  allow analysis of multiple expressions who have the same normal order
//  operator prefactor.
std::pair<ExprPtr, ExprPtr> fock_based_projector_2(ExprPtr exprs) {
  if (exprs->is<Constant>()) {
    return std::pair<ExprPtr, ExprPtr>{exprs, exprs};
  }
  non_canon_simplify(exprs);
  exprs = FNOPs_to_tens(exprs);
  non_canon_simplify(exprs);
  exprs = screen_densities(exprs);
  non_canon_simplify(exprs);
  //std::wcout << "pre partition expression: " << to_latex_align(exprs,30,3) << std::endl;
  //exprs = partition_F12(exprs);
  auto final_screen = exprs;
  non_canon_simplify(final_screen);
  // in some cases, there will now be no contributing terms left so return zero
  // to one and two body.
  if (final_screen->is<Constant>()) {
    return std::pair<ExprPtr, ExprPtr>{final_screen, final_screen};
  }
  non_canon_simplify(final_screen);
  // find the special f12 intermediates that cannot efficiently be solved
  // directly. This seems to work already for the general case!
  auto last_screen = ex<Constant>(0.0);
  for (auto&& product : final_screen->as<Sum>().summands()) {
    auto new_product = simplification::find_F12_interms(product);
    last_screen = last_screen + new_product;
  }
  non_canon_simplify(last_screen);
  return fnop_to_overlap(last_screen);
}
}  // namespace simplification
#ifndef SEQUANT_SIMPLIFICATIONS_H
#define SEQUANT_SIMPLIFICATIONS_H

#endif  // SEQUANT_SIMPLIFICATIONS_H
