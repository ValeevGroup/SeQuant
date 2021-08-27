//
// Created by Nakul Teke on 03/06/20.
//
// Taken from "cc_btas.cpp"

#ifndef SEQUANT_HAS_BTAS
#error "SEQUANT_HAS_BTAS should be defined when building cc_btas"
#endif

#include <iostream>
#include <memory>

#include <btas/btas.h>
#include <btas/tensorview.h>

#include "../sequant_setup.hpp"

#include "interpret/interpreted_tensor.hpp"
#include "scf/hartree-fock.h"
#include "interpret/contract.hpp"
#include <SeQuant/domain/evaluate/eval_expr.hpp>

using BTensor = btas::Tensor<double>;
using BTensorPtr = std::shared_ptr<BTensor>;

int main(int argc, char* argv[]) {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::microseconds;
  /* using std::uniform_real_distribution; */
  /* using std::mt19937; */
  //
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  try {
    // read geometry from a file; by default read from h2o.xyz, else take
    // filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename);

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i) nelectron += atoms[i].atomic_number;
    auto ndocc = nelectron / 2;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij * xij + yij * yij + zij * zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    std::cout << "** Nuclear repulsion energy = " << enuc << std::endl;

    auto shells = make_sto3g_basis(atoms);
    size_t nao = 0;
    for (auto s = 0; s < shells.size(); ++s) nao += shells[s].size();
    size_t nvirt = nao - ndocc;

    std::cout << "nao  : " << nao << "\n";
    std::cout << "ndocc: " << ndocc << "\n";
    std::cout << "nvirt: " << nvirt << "\n";

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    Matrix S = compute_1body_ints(shells, Operator::overlap);
    Matrix T = compute_1body_ints(shells, Operator::kinetic);
    Matrix V = compute_1body_ints(shells, Operator::nuclear, atoms);
    Matrix H = T + V;
    T.resize(0, 0);
    V.resize(0, 0);

    const auto use_hcore_guess = false;
    // use core Hamiltonian eigenstates to guess density?
    // set to true to match the result of versions 0, 1, and 2 of the code
    // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
    // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    Matrix D, C;
    Eigen::VectorXd eps;    // eigenvalues
    if (use_hcore_guess) {  // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      // eps = gen_eig_solver.eigenvalues();
      C = gen_eig_solver.eigenvectors();
      std::cout << "\n\tInitial C Matrix:\n";
      std::cout << C << std::endl;

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    } else {  // SOAD as the guess density, assumes STO-nG basis
      D = compute_soad(atoms);
    }

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rmsd = 0.0;
    auto ediff = 0.0;
    auto ehf = 0.0;
    Matrix Fock_matrix;  // capture fock matrix for ccsd calculations
    printf(
        "\n Iter        E(elec)              E(tot)               "
        "Delta(E)             RMS(D)         Time(s)\n");
    do {
      const auto tstart = high_resolution_clock::now();
      ++iter;
      auto ehf_last = ehf;
      auto D_last = D;
      auto F = H;
      F += compute_2body_fock(shells, D);
      Fock_matrix = F;

      // solve F C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
      C = gen_eig_solver.eigenvectors();
      eps = gen_eig_solver.eigenvalues();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();

      // compute HF energy
      ehf = 0.0;
      for (auto i = 0; i < nao; i++)
        for (auto j = 0; j < nao; j++) ehf += D(i, j) * (H(i, j) + F(i, j));

      // compute difference with last iteration
      ediff = ehf - ehf_last;
      rmsd = (D - D_last).norm();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf,
             ehf + enuc, ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    libint2::finalize();  // done with libint

    auto hf_final = ehf + enuc;
    std::cout << std::endl;
    printf("** Hartree-Fock energy = %20.12f a.u.\n", hf_final);

    // compute eri on MO basis in physicist's notation
    const auto ints_mo_spatial = compute_mo_ints(C, shells);
    size_t nocc = ndocc;  // in place of ndocc use nocc

    auto emp2 = 0.0;
    for (int a = ndocc; a < nao; a++) {
      for (int b = ndocc; b < nao; b++) {
        for (int i = 0; i < ndocc; i++) {
          for (int j = 0; j < ndocc; j++) {
            auto temp = ints_mo_spatial(i, j, a, b) *
                        (2.0 * ints_mo_spatial(i, j, a, b) -
                         ints_mo_spatial(j, i, a, b));
            emp2 += temp / (eps(i) + eps(j) - eps(a) - eps(b));
          }
        }
      }
    }

    printf("** MP2 energy = %20.12f a.u.\n** Total energy = %20.12f a.u.\n",
           emp2, hf_final + emp2);

#define CCSD_BTAS 0
#define CCSDT 0
    std::cout << "\n"
              << "***************************\n"
              << "Calculating ECCSD with BTAS\n"
              << "***************************\n"
              << std::endl;

    // D_ov and D_oovv
    auto D_ov = std::make_shared<BTensor>(nocc, nvirt);
    auto D_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);

#if CCSDT
    auto D_ooovvv =
        std::make_shared<BTensor>(nocc, nocc, nocc, nvirt, nvirt, nvirt);
#endif

    // Fock matrix tensors
    auto F_oo = std::make_shared<BTensor>(nocc, nocc);
    auto F_ov = std::make_shared<BTensor>(nocc, nvirt);
    auto F_vv = std::make_shared<BTensor>(nvirt, nvirt);

    // Integral tensors
    auto G_oooo = std::make_shared<BTensor>(nocc, nocc, nocc, nocc);
    auto G_vvvv = std::make_shared<BTensor>(nvirt, nvirt, nvirt, nvirt);
    auto G_ovvv = std::make_shared<BTensor>(nocc, nvirt, nvirt, nvirt);
    auto G_ooov = std::make_shared<BTensor>(nocc, nocc, nocc, nvirt);
    auto G_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
    auto G_vvoo = std::make_shared<BTensor>(nvirt, nvirt, nocc, nocc);
    auto G_ovov = std::make_shared<BTensor>(nocc, nvirt, nocc, nvirt);

    //  F_oo, F_ov, F_vv, D_ov, D_oovv, G_oovv, G_vvoo, G_ovov
    for (auto i = 0; i < nocc; ++i) {
      (*F_oo)(i, i) = Fock_matrix(i, i);
      for (auto a = 0; a < nvirt; ++a) {
        (*F_ov)(i, a) = Fock_matrix(i, a + nocc);
        (*D_ov)(i, a) = Fock_matrix(i, i) - Fock_matrix(a + nocc, a + nocc);
        for (auto j = 0; j < nocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            (*F_vv)(a, b) = Fock_matrix(a + nocc, b + nocc);
            (*D_oovv)(i, j, a, b) = (*D_ov)(i, a) + Fock_matrix(j, j) -
                                    Fock_matrix(b + nocc, b + nocc);
            (*G_oovv)(i, j, a, b) = ints_mo_spatial(i, j, a + nocc, b + nocc);
            (*G_vvoo)(a, b, i, j) = ints_mo_spatial(a + nocc, b + nocc, i, j);
            (*G_ovov)(i, a, j, b) = ints_mo_spatial(i, a + nocc, j, b + nocc);
          }
        }
      }
    }

    // G_ovvv
    for (auto i = 0; i < nocc; ++i)
      for (auto a = 0; a < nvirt; ++a)
        for (auto b = 0; b < nvirt; ++b)
          for (auto c = 0; c < nvirt; ++c)
            (*G_ovvv)(i, a, b, c) =
                ints_mo_spatial(i, a + nocc, b + nocc, c + nocc);

    // G_ooov
    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto a = 0; a < nvirt; ++a)
            (*G_ooov)(i, j, k, a) = ints_mo_spatial(i, j, k, a + nocc);

    // G_oooo
    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto l = 0; l < nocc; ++l)
            (*G_oooo)(i, j, k, l) = ints_mo_spatial(i, j, k, l);

    // G_vvvv
    for (auto a = 0; a < nvirt; ++a)
      for (auto b = 0; b < nvirt; ++b)
        for (auto c = 0; c < nvirt; ++c)
          for (auto d = 0; d < nvirt; ++d)
            (*G_vvvv)(a, b, c, d) =
                ints_mo_spatial(a + nocc, b + nocc, c + nocc, d + nocc);

#if CCSDT
    // D_ooovvv
    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto a = 0; a < nvirt; ++a)
            for (auto b = 0; b < nvirt; ++b)
              for (auto c = 0; c < nvirt; ++c)
                (*D_ooovvv)(i, j, k, a, b, c) =
                    (*D_oovv)(i, j, a, b) + (*D_ov)(k, c);
#endif

    // singles and doubles for CCSD calculation
    auto t_ov = std::make_shared<BTensor>(nocc, nvirt);
    auto t_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
    (*t_ov).fill(0.0);
    (*t_oovv).fill(0.0);

    // and triples if CCSDT required
#if CCSDT
    auto t_ooovvv =
        std::make_shared<BTensor>(nocc, nocc, nocc, nvirt, nvirt, nvirt);
    (*t_ooovvv).fill(0.0);
#endif

    // the map that is required while evaluating sequant expressions
    std::map<std::wstring, BTensorPtr> btensor_map;
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_oo", F_oo));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_ov", F_ov));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_vv", F_vv));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_oooo", G_oooo));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_vvvv", G_vvvv));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ovvv", G_ovvv));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ooov", G_ooov));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_oovv", G_oovv));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_vvoo", G_vvoo));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ovov", G_ovov));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"t_ov", t_ov));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"t_oovv", t_oovv));
#if CCSDT
    btensor_map.insert(
        std::pair<std::wstring, BTensorPtr>(L"t_ooovvv", t_ooovvv));
#endif

    // global sequant setup...
    std::locale::global(std::locale("en_US.UTF-8"));
    std::wcout.precision(std::numeric_limits<double>::max_digits10);
    std::wcerr.precision(std::numeric_limits<double>::max_digits10);
    std::wcout.sync_with_stdio(true);
    std::wcerr.sync_with_stdio(true);
    std::wcout.imbue(std::locale("en_US.UTF-8"));
    std::wcerr.imbue(std::locale("en_US.UTF-8"));
    std::wcout.sync_with_stdio(true);
    std::wcerr.sync_with_stdio(true);
    sequant::detail::OpIdRegistrar op_id_registrar;

    sequant::mbpt::set_default_convention();

    TensorCanonicalizer::register_instance(
        std::make_shared<DefaultTensorCanonicalizer>());

    //
    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto eccsd = 0.0;
    Logger::get_instance().wick_stats = false;
    auto ccsd_r = cceqvec{2, 2}(true, true, true, true, true);

    std::initializer_list<IndexList> external_indices = {{}};


    // SPIN TRACE THE RESIDUAL
    std::vector<ExprPtr> cc_r(ccsd_r.size());
    for (size_t i = 1; i < ccsd_r.size(); ++i){
      // std::wcout << "R" << i << ":\n" << to_latex(ccsd_r[i]) << std::endl;
      if(i == 1)
        external_indices = {{L"i_1", L"a_1"}};
      else if(i == 2)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};

      auto spin_removed_term = spintrace(ccsd_r[i], external_indices);
      expand(spin_removed_term);
      rapid_simplify(spin_removed_term);
      auto result = std::make_shared<Sum>();
      auto j_iter = 0;
      for (auto&& term: *spin_removed_term){
        ++j_iter;
        if(j_iter > 1)
          canonicalize(term);
        result->append(term);
      }
      rapid_simplify(spin_removed_term);
      cc_r[i] = result;
    }
    expand(cc_r[2]);
    rapid_simplify(cc_r[2]);
    canonicalize(cc_r[2]);
    // std::wcout << "R2 traced:\n" << to_latex(cc_r[2]) << "\n";

#define SIMPLIFIED_R2 1
#if SIMPLIFIED_R2
    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};

    // 1/3 R + 1/6 R' for simpler equations
    auto temp_expr = transform_expression(cc_r[2], idxmap);

    // std::wcout << to_latex(cc_r[2]) << "\n\n";
    auto simpler_R2 =
        ex<Constant>(1.0 / 3.0) * cc_r[2] + ex<Constant>(1.0 / 6.0) * temp_expr;

    expand(simpler_R2);
    rapid_simplify(simpler_R2);
    canonicalize(simpler_R2);
    rapid_simplify(simpler_R2);
    // std::wcout << "CCSD R2 traced:\n" << to_latex(simpler_R2) << "\n";

    // Trick to fool canonicalizer and reduce the number of terms
    auto P = ex<Tensor>(L"P", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}, Symmetry::nonsymm);
    auto P_R2 = P * simpler_R2;

    expand(P_R2);
    rapid_simplify(P_R2);
    canonicalize(P_R2);
    rapid_simplify(P_R2);
    std::wcout << "P2 times CCSD R2 traced:\n" << to_latex(expand_P_operator(P_R2)) << "\n";
#endif
    using sequant::factorize::factorize_expr;
    using sequant::interpret::antisymmetrize;
    using sequant::interpret::eval_equation;

#if 1 // Factorization
    // factorize CCSD equations
    auto index_size_map = std::make_shared<ispace_map>(ispace_map{});
    index_size_map->insert(
        ispace_pair(sequant::IndexSpace::active_occupied, nocc));
    index_size_map->insert(
        ispace_pair(sequant::IndexSpace::active_unoccupied, nvirt));

    for (auto i = 1; i < cc_r.size(); ++i){
      cc_r[i] = factorize_expr(cc_r.at(i), index_size_map, true);
      // std::wcout << to_latex(cc_r[i]) << std::endl;
    }
#endif

    auto cc_r1 = std::make_shared<sequant::evaluate::EvalTensor>(cc_r[1]);
    cc_r1->fill_btas_indices();

    // auto cc_r2 = std::make_shared<sequant::evaluate::EvalTensor>(cc_r[2]);
    auto cc_r2 = std::make_shared<sequant::evaluate::EvalTensor>(P_R2);
    cc_r2->fill_btas_indices();

    sequant::container::map<sequant::evaluate::hash_type, size_t>
        hash_counts_r1, hash_counts_r2;
    sequant::evaluate::fill_hash_counts(cc_r1, hash_counts_r1);
    sequant::evaluate::fill_hash_counts(cc_r2, hash_counts_r2);

    auto generate_hash_map_entry = [](std::wstring& label, std::wstring& spaces,
                                      std::shared_ptr<BTensor>& tnsr_ptr) {
      auto space_vector = container::svector<IndexSpace::Type>{};
      for (auto sp : spaces) {
        if (sp == 'o')
          space_vector.push_back(IndexSpace::active_occupied);
        else if (sp == 'v')
          space_vector.push_back(IndexSpace::active_unoccupied);
        else
          throw std::logic_error(
              "Use 'o' for active_occupied and 'v' for active_unoccupied "
              "IndexSpace!");
      }
      auto hval = sequant::evaluate::DataTensorSpecs(label, space_vector)
                      .get_hash_value();
      return sequant::evaluate::hash_to_dtensor_map<BTensor>::value_type(
          hval, tnsr_ptr);
    };

    std::cout << __LINE__ << "L " << std::endl;
    auto hash_to_ptr = sequant::evaluate::hash_to_dtensor_map<BTensor>();
    for (auto& item : btensor_map) {
      std::wstring label = std::wstring{item.first[0]};
      std::wstring spaces =
          std::wstring{item.first.substr(2, std::wstring::npos)};
      auto entry = generate_hash_map_entry(label, spaces, item.second);
      hash_to_ptr.insert(entry);
    }

    std::cout << __LINE__ << "L " << std::endl;
    auto context_r1 =
        sequant::evaluate::EvalContext<BTensor>(hash_to_ptr, hash_counts_r1);
    auto context_r2 =
        sequant::evaluate::EvalContext<BTensor>(hash_to_ptr, hash_counts_r2);

    auto start = high_resolution_clock::now();
    do {
      ++iter;
      // old method
      // auto R1 = eval_equation(cc_r[1], btensor_map).tensor();
      // auto R2 = eval_equation(cc_r[2], btensor_map).tensor();
      // auto R3 = eval_equation(cc_r[3], btensor_map).tensor();

      // new method
      auto R1 = sequant::evaluate::eval_evtensor(cc_r1, context_r1);
      auto R2 = sequant::evaluate::eval_evtensor(cc_r2, context_r2);

      // R2 = antisymmetrize(R2);
      // R3 = antisymmetrize(R3);

      std::cout << "using BTAS,    iter " << iter << std::endl;
      /* std::cout << "norm(R1) = " << std::sqrt(btas::dot(R1, R1)) <<
       * std::endl; */
      /* std::cout << "norm(R2) = " << std::sqrt(btas::dot(R2, R2)) <<
       * std::endl; */

      // save previous norms
      auto norm_last = std::sqrt(btas::dot(*t_oovv, *t_oovv));

      //////////////
      // Updating amplitudes
      // update t_ov
      for (auto i = 0; i < nocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) {
          (*t_ov)(i, a) += R1(i, a) / (*D_ov)(i, a);
        }
      }

      //
      // update t_oovv
      for (auto i = 0; i < nocc; ++i) {
        for (auto j = 0; j < nocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              (*t_oovv)(i, j, a, b) += R2(i, j, a, b) / (*D_oovv)(i, j, a, b);
            }
          }
        }
      }

#if CCSDT
      //
      // update t_ooovvv
       for (auto i = 0; i < nocc; ++i)
         for (auto j = 0; j < nocc; ++j)
           for (auto k = 0; k < nocc; ++k)
             for (auto a = 0; a < nvirt; ++a)
               for (auto b = 0; b < nvirt; ++b)
                 for (auto c = 0; c < nvirt; ++c)
                   t_ooovvv(i,j,k,a,b,c)
                     += R3(i,j,k,a,b,c)/D_ooovvv(i,j,k,a,b,c);
#endif

      std::cout << "norm(t_ov)    " << std::sqrt(btas::dot(*t_ov, *t_ov))
                << std::endl;
      std::cout << "norm(t_oovv)  " << std::sqrt(btas::dot(*t_oovv, *t_oovv))
                << std::endl;

      auto eccsd_last = eccsd;

      // calculating energy
      BTensor temp_tensor;
      enum { i, j, a, b };
      btas::contract(1.0, *G_oovv, {i, j, a, b}, *t_ov, {i, a}, 0.0,
                     temp_tensor, {j, b});
      //
      eccsd = 0.5 * btas::dot(temp_tensor, *t_ov) +
              0.25 * btas::dot(*G_oovv, *t_oovv) + btas::dot(*F_ov, *t_ov);
      printf("E(CC) is: %20.12f\n\n", eccsd);

      normdiff = norm_last - sqrt(btas::dot(*t_oovv, *t_oovv));
      ediff = eccsd_last - eccsd;

    } while ((fabs(normdiff) > conv || fabs(ediff) > conv) && (iter < maxiter));

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "\nOut of loop after " << iter << " iterations.\n"
              << "\nTime: " << duration.count() << " microseconds."
              << std::endl;

  }  // end of try block; if any exceptions occurred, report them and exit
     // cleanly
  catch (const char* ex) {
    std::cerr << "caught exception: " << ex << std::endl;
    return 1;
  } catch (std::string& ex) {
    std::cerr << "caught exception: " << ex << std::endl;
    return 1;
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "caught unknown exception\n";
    return 1;
  }
  return 0;
}
