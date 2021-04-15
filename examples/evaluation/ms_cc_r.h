//
// Created by Nakul Teke on 9/12/20.
//

#ifndef SEQUANT_MS_CC_R_H
#define SEQUANT_MS_CC_R_H

#include "../../examples/contract/scf/hartree-fock.h"
#include "../sequant_setup.hpp"

#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <TiledArray/initialize.h>
#include <tiledarray.h>
#include <libint2.hpp>

#include <memory>

auto expanded_T3(container::svector<Index> bra, container::svector<Index> ket){
  assert(bra.size() == 3);
  assert(bra.size() == ket.size());

  auto t3_1 = Tensor(L"t", bra, ket);
  container::svector<Index> ket2 = {ket[1], ket[0], ket[2]};
  auto t3_2 = Tensor(L"t", bra, ket2);
  container::svector<Index> ket3 = {ket[2], ket[1], ket[0]};
  auto t3_3 = Tensor(L"t", bra, ket3);
  auto result = ex<Constant>(2.0) * ex<Tensor>(t3_1) - ex<Tensor>(t3_2) - ex<Tensor>(t3_3);
  return result;
}

// CCSD R1

ExprPtr r1(const int n){
  // 36
  auto F_em = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"i_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});

  // 34
  auto F_im = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});

  // 37'
  auto W_eimn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});

  // 37"
  auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});

  auto tau_mief = ex<Tensor>(L"t", WstrList{L"i_5", L"i_1"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  auto tau_imef = ex<Tensor>(L"t", WstrList{L"i_1", L"i_5"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_6"});
  auto tau_p = ex<Constant>(2.0) * tau_mief->clone() - tau_imef->clone();
//  std::wcout << "tau_mief: " << to_latex(tau_mief) << "\n";
//  std::wcout << "tau_imef: " << to_latex(tau_imef) << "\n";
//  std::wcout << "tau_p: " << to_latex(tau_p) << "\n";

  auto r1 = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) +
      ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"a_1"}) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) -
      F_im->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_5"}) +
      F_em->clone() * (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_1"}, WstrList{L"a_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_1", L"i_5"}, WstrList{L"a_5", L"a_1"}, Symmetry::nonsymm)) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * tau_p->clone() -
      (ex<Constant>(2.0) * W_eimn->clone() - W_iemn->clone()) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_5", L"a_1"}, Symmetry::nonsymm);

  if(n==3){
    Index a(L"a_1"),
        b(L"a_2"),
        e(L"a_5"),
        f(L"a_6"),
        i(L"i_1"),
        j(L"i_2"),
        m(L"i_5"),
        n(L"i_6");

    container::svector<Index> efa = {e, f, a};
    container::svector<Index> eaf = {e, a, f};
    container::svector<Index> mni = {m, n, i};

    auto ccsdt_r1 = ex<Constant>(1./2) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
        (ex<Constant>(2.0) * expanded_T3(efa, mni) - expanded_T3(eaf, mni));
    r1 = r1->clone() + ccsdt_r1;
  }
  // std::wcout << "CC R1: " << to_latex(r1) << "\n";
  return r1;
}

// CCSD R2
ExprPtr r2(const int n){
  //﻿http://dx.doi.org/10.1063/1.4907278

  // 39
  auto W_ijam = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
  // std::wcout << "39 W_ijam: " << to_latex(W_ijam) << "\n\n";

  // 36
  auto F_em = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"i_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "36 F_em: " << to_latex(F_em) << "\n\n";

  // 35
  auto F_ea = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"a_1"}) -
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) -
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_5"}, WstrList{L"a_6", L"a_1"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_6"});
  // std::wcout << "35 F_ea: " << to_latex(F_ea) << "\n\n";

  // 34
  auto F_im = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "34 F_im: " << to_latex(F_im) << "\n\n";

  // 37
  auto W_ejmn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
  // std::wcout << "37 W_ejmn: " << to_latex(W_ejmn) << "\n\n";

  // 37'
  auto W_eimn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "W_eimn " << to_latex(W_eimn) << "\n\n";

  // 37"
//    auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) +
//        ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "W_iemn " << to_latex(W_iemn) << "\n\n";

  // 43 // Antisymmetrized
  auto W_eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
      W_eimn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) +
      ex<Constant>(0.25) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
      ex<Constant>(0.25) *  ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);
  // std::wcout << "W_eima " << to_latex(W_eima) << "\n\n";

  // 44
  auto W_iema = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
      W_iemn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
      ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) -
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"});
  // std::wcout << "W_iema " << to_latex(W_iema) << "\n\n";

  // 49
  auto W_ijmn_temp = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
          (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));
  auto W_ijmn_temp_c = W_ijmn_temp->clone();
  expand(W_ijmn_temp_c);
  // std::wcout << "W_ijmn_temp_c: " << to_latex(W_ijmn_temp_c) << "\n\n";

  std::map<Index, Index> P_ij_mn;
  {
    Index i(L"i_1");
    Index j(L"i_2");
    Index m(L"i_5");
    Index n(L"i_6");
    P_ij_mn.emplace(std::make_pair(i, j));
    P_ij_mn.emplace(std::make_pair(j, i));
    P_ij_mn.emplace(std::make_pair(m, n));
    P_ij_mn.emplace(std::make_pair(n, m));
  }
  // std::wcout << "W_ijmn_temp transformed: " << to_latex(transform_expression(W_ijmn_temp_c->clone(), P_ij_mn)) << "\n\n";

  auto W_ijmn = W_ijmn_temp_c->clone() + transform_expression(W_ijmn_temp_c->clone(), P_ij_mn);
  // std::wcout << "W_ijmn " << to_latex(W_ijmn) << "\n\n";

  //
  auto temp_ab_ = W_iema->clone() * ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm);
  auto temp_ab_c = temp_ab_->clone();
  // std::wcout << "temp_ab_c " << to_latex(temp_ab_c) << "\n\n";
  expand(temp_ab_c);

  std::map<Index, Index> P_ab;
  {
    Index a(L"a_1");
    Index b(L"a_2");
    P_ab.emplace(std::make_pair(a, b));
    P_ab.emplace(std::make_pair(b, a));
  }
  auto temp_ab = ex<Constant>(0.5) * temp_ab_c->clone() + transform_expression(temp_ab_c->clone(), P_ab);
  // std::wcout << "temp_ab " << to_latex(temp_ab) << "\n\n";

  // CCSD Z2
  auto temp1 = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) -
      W_ijam->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}) +
      F_ea->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2"}, WstrList{L"a_5", L"a_2"}) -
      F_im->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_2"}, WstrList{L"a_1", L"a_2"}) +
      ex<Constant>(0.5) * (ex<Constant>(2.0) * W_eima->clone() - W_iema->clone()) *
          (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_2", L"a_5"}, Symmetry::nonsymm)) -
      temp_ab->clone() +
      ex<Constant>(0.5) * W_ijmn->clone() * (ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"}, Symmetry::nonsymm)) +
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) *
          (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}, Symmetry::nonsymm));

  auto temp1_c = temp1->clone();
  expand(temp1_c);

  auto ccsd_z2 = temp1_c->clone();
  if(n == 3){
    ExprPtr ccsdt_r2;
    {
      auto W_efam = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) -
          ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"});
      // std::wcout << "W_efam " << to_latex(W_efam) << "\n\n";

      Index a(L"a_1"),
          b(L"a_2"),
          e(L"a_5"),
          f(L"a_6"),
          i(L"i_1"),
          j(L"i_2"),
          m(L"i_5"),
          n(L"i_6");

      container::svector<Index> mij = {m, i, j};
      container::svector<Index> min = {m, i, n};
      container::svector<Index> eab = {e, a, b};
      container::svector<Index> feb = {f, e, b};

      ccsdt_r2 = ex<Constant>(0.5) * F_em->clone() * expanded_T3(mij, eab) +
          W_efam->clone() * expanded_T3(mij, feb) -
          W_ejmn->clone() * expanded_T3(min, eab);
      // std::wcout << "ccsdt_r2: " << ccsdt_r2->size() << " " << to_latex(ccsdt_r2) << "\n\n";
    }
    ccsd_z2 = ccsd_z2->clone() + ccsdt_r2;
  }

  ccsd_z2 = ex<Tensor>(L"S", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}) * ccsd_z2->clone(); // + transform_expression(temp1, P_ab_ij);
  expand(ccsd_z2);
  // std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";

  ccsd_z2 = swap_bra_ket(ccsd_z2);
  // std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";

  canonicalize(ccsd_z2);
  rapid_simplify(ccsd_z2);
  // std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";

  simplify(ccsd_z2);
  return ccsd_z2;
}

// CCSDT R3
ExprPtr r3(const int r){

  // 46'
  auto W___ejam = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) -
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);
  // std::wcout << "46\' W___ejam " << to_latex(W___ejam) << "\n\n";

  // 45'
  auto W___ejmb = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"a_2"}, Symmetry::nonsymm) +
      ex<Constant>(0.5) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_2"}, WstrList{L"a_6", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_2"}, WstrList{L"a_2", L"a_6"}, Symmetry::nonsymm)) -
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2", L"i_6"}, WstrList{L"a_6", L"a_2"}, Symmetry::nonsymm);
  // std::wcout << "45\' W___ejmb " << to_latex(W___ejmb) << "\n\n";

  // 37
  auto W_ejmn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
  // std::wcout << "37 W_ejmn " << to_latex(W_ejmn) << "\n\n";

  // 41
  ExprPtr W_ejab;
  {
    std::map<Index, Index> P_ab;
    {
      Index a(L"a_1");
      Index b(L"a_2");
      P_ab.emplace(std::make_pair(a, b));
      P_ab.emplace(std::make_pair(b, a));
    }

    auto g_tau = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) *
        (ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_6", L"a_2"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}));
    auto P_g_tau = ex<Constant>(0.5) * g_tau->clone() + transform_expression(g_tau->clone(), P_ab);
    // std::wcout << "P_g_tau " << to_latex(P_g_tau) << "\n\n";

    Index a(L"a_1"),
        b(L"a_2"),
        e(L"a_5"),
        f(L"a_6"),
        i(L"i_1"),
        j(L"i_2"),
        m(L"i_5"),
        n(L"i_6");

    container::svector<Index> nmj = {n, m, j};
    container::svector<Index> fab = {f, a, b};

    // 41
    W_ejab = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) -
        W___ejam->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}) -
        W___ejmb->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) +
        W_ejmn->clone() * (ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"})) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) +
        ex<Constant>(0.5) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_6", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_6", L"a_2"}, Symmetry::nonsymm) -
                ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"})) -
        P_g_tau->clone() -
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * expanded_T3(nmj, fab);
    // std::wcout << "41 W_ejab " << to_latex(W_ejab) << "\n\n";
  }

  // 36
  auto F_em = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"i_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "36 F_em " << to_latex(F_em) << "\n\n";


  // 40
  ExprPtr W_ijam;
  {
    // 39
    auto W__ijam = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
        ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));

    // 49* W__ijnm
    ExprPtr W_ijnm;
    {
      auto W_ijnm_temp = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) +
          ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
          ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) *
              (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));
      auto W_ijnm_temp_c = W_ijnm_temp->clone();
      expand(W_ijnm_temp_c);
      // std::wcout << "W_ijnm_temp_c: " << to_latex(W_ijnm_temp_c) << "\n";

      std::map<Index, Index> P_ij_mn;
      {
        Index i(L"i_1");
        Index j(L"i_2");
        Index m(L"i_5");
        Index n(L"i_6");
        P_ij_mn.emplace(std::make_pair(i, j));
        P_ij_mn.emplace(std::make_pair(j, i));
        P_ij_mn.emplace(std::make_pair(m, n));
        P_ij_mn.emplace(std::make_pair(n, m));
      }
      W_ijnm = W_ijnm_temp_c->clone() + transform_expression(W_ijnm_temp_c->clone(), P_ij_mn);
    }


    Index a(L"a_1"), b(L"a_2"), e(L"a_5"), f(L"a_6"),
        i(L"i_1"), j(L"i_2"), m(L"i_5"), n(L"i_6");

    ExprPtr P_g_t2;
    {
      auto g_t2 = ex<Tensor>(L"g", WstrList{L"i_2", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_5", L"a_1"}, Symmetry::nonsymm);
      std::map<Index, Index> P_ij;
      P_ij.emplace(std::make_pair(i, j));
      P_ij.emplace(std::make_pair(j, i));
      P_g_t2 = ex<Constant>(1./2) * g_t2->clone() + transform_expression(g_t2->clone(), P_ij);
    }

    container::svector<Index> nij = {n, i, j};
    container::svector<Index> fae = {f, a, e};

    W_ijam = W__ijam->clone() +
             ex<Constant>(1./2) *
                 (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
                 (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_5"}, Symmetry::nonsymm)) -
             P_g_t2->clone() -
             W_ijnm->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}, Symmetry::nonsymm) +
             F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_5"}, Symmetry::nonsymm) +
             ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * expanded_T3(nij, fae);
  }

  // 35
  auto F_ea = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"a_1"}) -
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) -
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_5"}, WstrList{L"a_6", L"a_1"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_6"});
  // std::wcout << "35 F_ea " << to_latex(F_ea) << "\n\n";

  // 34
  auto F_im = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "35 F_im " << to_latex(F_im) << "\n\n";

  // 37'
  auto W_eimn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "37\' W_eimn " << to_latex(W_eimn) << "\n\n";

  // 37"
  auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "37\" W_iemn " << to_latex(W_iemn) << "\n\n";

  ExprPtr W_eima; // 47
  {
    // 43
    auto W__eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
        W_eimn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) +
        ex<Constant>(0.25) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
        ex<Constant>(0.25) *  ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);

    // 45
    auto W___eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) +
        ex<Constant>(0.5) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);


    W_eima = W__eima + ex<Constant>(0.5) * W___eima - ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm);
  }

  ExprPtr W_iema; // 48
  {
    // 44
    auto W__iema = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
        W_iemn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
        ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) -
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"});

    // 46
    auto W___iema = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm); // <======

    W_iema = W__iema + ex<Constant>(0.5) * W___iema - ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm);
  }

  Index a(L"a_1"),
      b(L"a_2"),
      c(L"a_3"),
      d(L"a_4"),
      e(L"a_5"),
      f(L"a_6"),

      i(L"i_1"),
      j(L"i_2"),
      k(L"i_3"),
      l(L"i_4"),
      m(L"i_5"),
      n(L"i_6");

  container::svector<Index> mjk = {m, j, k};
  container::svector<Index> ebc = {e, b, c};

  ExprPtr P_W_T3;
  {
    std::map<Index, Index> P_ab;
    P_ab.emplace(std::make_pair(a, b));
    P_ab.emplace(std::make_pair(b, a));

    auto term = W_iema->clone() * ex<Tensor>(L"t", WstrList{ L"i_2", L"i_5", L"i_3"}, WstrList{L"a_5", L"a_2", L"a_3"});
    expand(term);
    P_W_T3 = ex<Constant>(0.5) * term->clone() + transform_expression(term->clone(), P_ab);
  }

  ExprPtr W_ijmn;
  {
    auto W_ijmn_temp = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
            (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));
    auto W_ijmn_temp_c = W_ijmn_temp->clone();
    expand(W_ijmn_temp_c); // This has to be expanded for correct expression transformation

    std::map<Index, Index> P_ij_mn;
      P_ij_mn.emplace(std::make_pair(i, j));
      P_ij_mn.emplace(std::make_pair(j, i));
      P_ij_mn.emplace(std::make_pair(m, n));
      P_ij_mn.emplace(std::make_pair(n, m));

    W_ijmn = W_ijmn_temp_c->clone() + transform_expression(W_ijmn_temp_c->clone(), P_ij_mn);
  }

  ExprPtr W_efab;
  {
    std::map<Index, Index> P_ef_ab;
    P_ef_ab.emplace(std::make_pair(e, f));
    P_ef_ab.emplace(std::make_pair(f, e));
    P_ef_ab.emplace(std::make_pair(a, b));
    P_ef_ab.emplace(std::make_pair(b, a));

    auto terms = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) -
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{ L"i_5"}, WstrList{L"a_2"}) +
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
        (ex<Tensor>(L"t", WstrList{ L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"}));
    expand(terms);
    W_efab = terms->clone() + transform_expression(terms->clone(), P_ef_ab);
  }

  auto ccsdt_r3 =
      W_ejab->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_3"}, WstrList{L"a_5", L"a_3"}) -
      W_ijam->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_3"}, WstrList{L"a_2", L"a_3"}) +
      ex<Constant>(0.5) * F_ea->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2", L"i_3"}, WstrList{L"a_5", L"a_2", L"a_3"}) -
      ex<Constant>(0.5) * F_im->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}) +
      ex<Constant>(0.25) * (ex<Constant>(2.0) * W_eima->clone() - W_iema->clone()) * expanded_T3(mjk, ebc) -
      P_W_T3 +
      ex<Constant>(0.5) * W_ijmn * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_6", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}) +
      ex<Constant>(0.5) * W_efab * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2", L"i_3"}, WstrList{L"a_5", L"a_6", L"a_3"});

  expand(ccsdt_r3);
  ccsdt_r3 = swap_bra_ket(ccsdt_r3);
  ccsdt_r3 = ex<Tensor>(L"S", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}) * ccsdt_r3->clone();
  simplify(ccsdt_r3);
  return ccsdt_r3;
}


#endif  // SEQUANT_MS_CC_R_H
