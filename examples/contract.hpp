#ifndef SEQUANT_CONTRACT_HPP
#define SEQUANT_CONTRACT_HPP

#include <vector>
#include <btas/btas.h>

#include <SeQuant/core/timer.hpp>

#include "contractable_tensor.hpp"

const char ADD{ 1};
const char SUB{-1};
struct Perm {
  std::vector<size_t> perm;
  char even_perm{ADD};
};

using STensor = sequant::Tensor;
using CTensor = sequant::contractable::Tensor;
using BTensor = btas::Tensor<double>;

namespace sequant{
contractable::Tensor contract(const contractable::Tensor & t1,
    const contractable::Tensor & t2, double scal=1.0) {

  using ords_vec    = std::vector<size_t>;
  using index_vec   = std::vector<sequant::Index>;

  using ind_ord     = std::pair<sequant::Index, size_t>;
  using ind_ord_map = std::map<ind_ord::first_type, ind_ord::second_type>;

  ind_ord_map index_to_ordinal; // all indices to ordinal map
  ind_ord_map nc_indices;       // non-contracting indices to ordinal map
  size_t count = 0;

  auto update_index_to_ordinal = [&](const index_vec& vec){

    auto map_size = index_to_ordinal.size();

    for (const auto& i: vec) {
      index_to_ordinal.insert(ind_ord(i, ++count));
      auto current_size = index_to_ordinal.size();
      if (map_size < current_size) {
        nc_indices.insert(ind_ord(i, count));
        ++map_size;
      }
      else
        nc_indices.erase(i);
    }
  };

  update_index_to_ordinal(t1.occs());
  update_index_to_ordinal(t1.virts());
  update_index_to_ordinal(t2.occs());
  update_index_to_ordinal(t2.virts());

  auto gen_ords_vec = [&index_to_ordinal](const contractable::Tensor& t){
    ords_vec v;
    for (const auto& i: t.occs())
      v.push_back(index_to_ordinal.at(i));

    for (const auto& i: t.virts())
      v.push_back(index_to_ordinal.at(i));
    return v;
  };

  const ords_vec t1_ords = gen_ords_vec(t1);

  const ords_vec t2_ords = gen_ords_vec(t2);

  ords_vec nc_ords;
  index_vec nc_index_vec;
  for (const auto& i: nc_indices) {
    nc_index_vec.push_back(i.first);
    nc_ords.push_back(i.second);
  }

  // create a contractable::Tensor object to return
  contractable::Tensor result{t1.label(), nc_index_vec};
  auto btensor = BTensor{};
  // call btas::contract
  btas::contract(scal, t1.bt(), t1_ords,
      t2.bt(), t2_ords,
      0.0, btensor, nc_ords);
  result.link_btensor(btensor);
  return result;
} // contract function for contractable::Tensor objects

} // namespace sequant

using btmap = std::map<std::wstring, BTensor*>;
using sequant::contract;

CTensor contract_vec(const std::vector<CTensor>& vct, size_t n, double scal=1.0, size_t b=1);
CTensor eval_equation(const sequant::ExprPtr&, const btmap&);
CTensor eval_product(const sequant::ExprPtr&, const btmap&);
CTensor eval_sum(const sequant::ExprPtr&, const btmap&);

CTensor contract_vec(const std::vector<CTensor>& vct, size_t n, double scal, size_t b) {
  // evaluates left to right
  auto i = n - b;
  if (i == 0) {
    if (scal != 1.0) {
      auto tmp = vct[i];
      tmp.scale_btensor(scal);
      return tmp;
    }
    return vct[i];
  }
  else if (i == 1)
    return contract(vct[i-1], vct[i], scal);
  // else
  return contract(vct[n-b], contract_vec(vct, n, 1.0, b+1), scal);
}

CTensor eval_equation(const sequant::ExprPtr& expr,
                      const btmap& tmap) {

  if (expr->is<sequant::Product>())
    return eval_product(expr, tmap);

  else if (expr->is<sequant::Sum>())
    return eval_sum(expr, tmap);

  else if (expr->is<sequant::Tensor>()) {
    auto ct = CTensor(expr->as<sequant::Tensor>());
    ct.link_btensor(*(tmap.find(ct.label()+L"_"+ct.translate())->second));
    return ct;
  }

  else
    throw "Only know how to handle sum or product!";
}

CTensor eval_product(const sequant::ExprPtr& expr,
                     const btmap& tmap) {

  auto p = expr->as<sequant::Product>();

  // collect the tensors to be contracted
  std::vector<CTensor> fvec;
  for (const auto& f: p.factors()) {
    if (!f->is<STensor>())
      fvec.push_back(eval_equation(f, tmap)); // call equation evaluater if factor is not a tensor
    else {
      auto t = f->as<STensor>();
      if (t.label() != L"A") // "A": antisymmetrizing tensors are ommitted at this point
        fvec.push_back(CTensor{t});
    }
  }
  // linking the suitable btas tensors to the CTensors
  for (auto& f: fvec) {
    f.link_btensor(*(tmap.find(f.label() + L"_" + f.translate())->second));
  }
  return contract_vec(fvec, fvec.size(), std::real(p.scalar()));
}

CTensor eval_sum(const sequant::ExprPtr& expr,
                 const btmap& tmap) {

  auto smands = expr->as<sequant::Sum>().summands();
  // collect the tensors to be summed
  std::vector<CTensor> svec;
  //
  for (auto i=0; i < smands.size(); ++i) {
    if (! smands.at(i)->is<STensor>()) 
      // call equation evaluater if summand is not a tensor
      svec.push_back(eval_equation(smands[i], tmap));
    else {
      auto t = smands[i]->as<STensor>();
      if (t.label() != L"A") { // 'A' tensors are omitted
        svec.push_back(CTensor{t});
        svec[i].link_btensor(*(tmap.find(svec[i].label() + L"_" + svec[i].translate())->second));
      }
    }
  }
  auto result = svec[0];
  for (auto i = 1; i < svec.size(); ++i) {
    result = result + svec[i];
  }
  return result;
}

std::vector<Perm> perm_calc(std::vector<size_t> to_perm,
                                         size_t  size,
                                         size_t cswap = 0, // count swaps
                                         size_t begin = 0) {
  if (begin+1 == size)
    return std::vector<Perm>{Perm{to_perm, (cswap%2==0)? ADD : SUB}}; // even permutations added. odds subtracted

  auto result = std::vector<Perm>{};
  for (auto i = begin; i < size; ++i) {
    std::swap(to_perm[begin], to_perm[i]);
    auto more_result = perm_calc(to_perm, size, (begin==i)? cswap: cswap+1, begin+1);
    for (auto p: more_result) {
      result.push_back(p);
    }
    std::swap(to_perm[begin], to_perm[i]);
  }
  return result;
}

BTensor antisymmetrize(const BTensor& bt) {

  size_t rank = bt.rank();

  if (rank == 2)
    return bt;
  else if (rank%2 != 0)
    throw "Can't handle odd ranked tensors!";

  auto result = BTensor(bt.range());
  result.fill(0);

  std::vector<size_t> to_perm;
  for (auto i = 0; i < (size_t)rank/2; ++i)
    to_perm.push_back(i);
  auto vp = perm_calc(to_perm, (size_t)rank/2);

  // antisymmetrize
  for (const auto& p: vp) {
    for (const auto& q: vp) {
      // permutation of the bra
      auto perm_vec = p.perm;
      // permutation of the ket
      // ket indices = rank/2 + bra indices
      for (auto qq: q.perm)
        perm_vec.push_back((size_t)rank/2  + qq);

      // permute and add
      auto perm_t = BTensor(btas::permute(bt, perm_vec));

      if (p.even_perm*q.even_perm == ADD)
        result += perm_t;
      else
        result -= perm_t;
    } // for q: vp
  } // for p: vp

  return result;
}

#endif /* SEQUANT_CONTRACT_HPP */
