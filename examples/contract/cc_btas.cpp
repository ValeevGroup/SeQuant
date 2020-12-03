#include "../sequant_setup.hpp"

#include "scf/hartree-fock.h"
#include "interpret/interpreted_tensor.hpp"
#include "interpret/contract.hpp"

#include <SeQuant/domain/evaluate/eval_expr.hpp>

#ifndef SEQUANT_HAS_BTAS
# error "SEQUANT_HAS_BTAS should be defined when building cc_btas"
#endif
#include <btas/btas.h>
#include <btas/tensorview.h>

#include <memory>
#include <iostream>

using BTensor = btas::Tensor<double>;
using BTensorPtr = std::shared_ptr<BTensor>;

int main(int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::wcout;
  using std::endl;
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::microseconds;
  /* using std::uniform_real_distribution; */
  /* using std::mt19937; */
  using std::vector;
  //
  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  try {
    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename);

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    auto ndocc = nelectron / 2;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij*xij + yij*yij + zij*zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

    auto shells = make_sto3g_basis(atoms);
    size_t nao = 0;
    for (auto s=0; s<shells.size(); ++s)
      nao += shells[s].size();
    size_t nvirt = nao -ndocc;

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute overlap integrals
    auto S = compute_1body_ints(shells, Operator::overlap);
    cout << "\n\tOverlap Integrals:\n";
    cout << S << endl;

    // compute kinetic-energy integrals
    auto T = compute_1body_ints(shells, Operator::kinetic);
    cout << "\n\tKinetic-Energy Integrals:\n";
    cout << T << endl;

    // compute nuclear-attraction integrals
    Matrix V = compute_1body_ints(shells, Operator::nuclear, atoms);
    cout << "\n\tNuclear Attraction Integrals:\n";
    cout << V << endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    cout << "\n\tCore Hamiltonian:\n";
    cout << H << endl;

    // T and V no longer needed, free up the memory
    T.resize(0,0);
    V.resize(0,0);

    /*** =========================== ***/
    /*** build initial-guess density ***/
    /*** =========================== ***/

    const auto use_hcore_guess = false;  // use core Hamiltonian eigenstates to guess density?
    // set to true to match the result of versions 0, 1, and 2 of the code
    // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
    // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    Matrix D;
    Matrix C;
    Eigen::VectorXd C_v; // eigenvalues
    if (use_hcore_guess) { // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();
      cout << "\n\tInitial C Matrix:\n";
      cout << C << endl;

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    }
    else {  // SOAD as the guess density, assumes STO-nG basis
      D = compute_soad(atoms);
    }

    cout << "\n\tInitial Density Matrix:\n";
    cout << D << endl;

    /*** =========================== ***/
    /*** main iterative loop         ***/
    /*** =========================== ***/

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rmsd = 0.0;
    auto ediff = 0.0;
    auto ehf = 0.0;
    Matrix fock_mat; // capture fock matrix for ccsd calculations
    do {
      const auto tstart = high_resolution_clock::now();
      ++iter;

      // Save a copy of the energy and the density
      auto ehf_last = ehf;
      auto D_last = D;

      // build a new Fock matrix
      auto F = H;
      //F += compute_2body_fock_simple(shells, D);
      F += compute_2body_fock(shells, D);
      fock_mat = F;

      if (iter == 1) {
        cout << "\n\tFock Matrix:\n";
        cout << F << endl;
      }

      // solve F C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
      auto eps = gen_eig_solver.eigenvalues();
      C = gen_eig_solver.eigenvectors();
      C_v = gen_eig_solver.eigenvalues();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();

      // compute HF energy
      ehf = 0.0;
      for (auto i = 0; i < nao; i++)
        for (auto j = 0; j < nao; j++)
          ehf += D(i,j) * (H(i,j) + F(i,j));

      // compute difference with last iteration
      ediff = ehf - ehf_last;
      rmsd = (D - D_last).norm();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        printf(
            "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n");
      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf, ehf + enuc,
          ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    libint2::finalize(); // done with libint

    auto hf_final = ehf + enuc;
    cout << endl;
    printf("** Hartree-Fock energy = %20.12f\n", hf_final);

    // compute eri on MO basis in physicist's notation
    const auto ints_mo_spatial = compute_mo_ints(C, shells);

    // convert ERI to Spin Basis
    //
    // from here on number of occupied and virtual orbitals
    // is doubled as we are now working on spin basis
    nao         *= 2;
    nvirt       *= 2;
    size_t nocc  = 2*ndocc; // in place of ndocc use nocc
    //
    const auto ints_mo_spin = [&]{
      BTensor result(nao, nao, nao, nao);
      for (auto r = 0; r < nao; ++r) {
        for (auto s = 0; s < nao; ++s) {
          for (auto p = 0; p < nao; ++p) {
            for (auto q = 0; q < nao; ++q) {
              //
              // finding the spatial orbital position
              // corresponding to a spin orbital
              //
              // spin    spatial
              // ---------------
              // 0,1        -> 0
              // 2,3        -> 1
              // 4,5        -> 2
              // 2*n, 2*n+1 -> n
              // x          -> floor(x/2)
              //
              size_t p_i = floor(p/2);
              size_t q_i = floor(q/2);
              size_t r_i = floor(r/2);
              size_t s_i = floor(s/2);

              auto col_int = 0.0;
              auto exc_int = 0.0;

              // orbs: 0 1 2 3 4 ...
              // spin: α β α β α ...
              // 
              //  <01|23>: coulomb integral

              if ((r%2 == p%2) && (s%2 == q%2))
                col_int  += ints_mo_spatial(r_i, s_i, p_i, q_i);

              if ((r%2 == q%2) && (s%2 == p%2))
                exc_int  += ints_mo_spatial(r_i, s_i, q_i, p_i);

              result(r,s,p,q) = col_int - exc_int;
            }
          }
        }
      }
      return result;
    }();


    // convert Fock matrix to Spin Basis
    const auto fock_spin = [&]{
      BTensor result(nao, nao);
      for (auto i = 0; i < nao; ++i) {
        for (auto j = 0; j < nao; ++j)
          result(i,j) = (i == j)? C_v(floor(i/2)): 0;
      }
      return result;
    }();

    cout << "\n"
         << "***************************\n"
         << "Calculating EMP2 using BTAS\n"
         << "***************************\n" << endl;

    auto emp2 = 0.0;

    for (auto r = 0; r < nocc; ++r) {
      for (auto s = 0; s < nocc; ++s) {
        for (auto p = nocc; p < nao; ++p) {
          for (auto q = nocc; q < nao; ++q) {

            auto calc = ints_mo_spin(r,s,p,q);

            calc *= calc;

            emp2 += calc/(fock_spin(r,r)
                + fock_spin(s,s)
                - fock_spin(p,p)
                - fock_spin(q,q)); } } } }
    emp2 /= 4.0;

    cout << "E(MP2): "
         << emp2
         << "\nFinal energy is : " << hf_final + emp2 << endl;
    /******************************************************/

    cout << "\n"
         << "*****************************\n"
         << "Calculating ECCSD using BTAS\n"
         << "*****************************\n" << endl;

    // D_ov and D_oovv
    auto D_ov = std::make_shared<BTensor>(nocc, nvirt);
    auto D_oovv  = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
    auto D_ooovvv = std::make_shared<BTensor>(nocc, nocc, nocc, nvirt, nvirt, nvirt);

    // Fock matrix tensors
    auto  Fock_oo = std::make_shared<BTensor>(nocc, nocc);
    auto  Fock_ov = std::make_shared<BTensor>(nocc, nvirt);
    auto  Fock_vv = std::make_shared<BTensor>(nvirt, nvirt);

    // Integral tensors
    auto  G_oooo = std::make_shared<BTensor>(nocc, nocc, nocc, nocc);
    auto  G_vvvv = std::make_shared<BTensor>(nvirt, nvirt, nvirt, nvirt);
    auto  G_ovvv = std::make_shared<BTensor>(nocc, nvirt, nvirt, nvirt);
    auto  G_ooov = std::make_shared<BTensor>(nocc, nocc, nocc, nvirt);
    auto  G_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
    auto  G_ovov = std::make_shared<BTensor>(nocc, nvirt, nocc, nvirt);

    for (auto i = 0; i < nocc; ++i) {
      (*Fock_oo)(i,i) = fock_spin(i,i);
      for (auto a = 0; a < nvirt; ++a) {
        (*Fock_ov)(i,a) = fock_spin(i,a+nocc);
        (*D_ov)(i,a) = fock_spin(i,i) - fock_spin(a+nocc,a+nocc);
        for (auto j = 0; j < nocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            (*Fock_vv)(a,b) = fock_spin(a+nocc,b+nocc);
            (*D_oovv)(i,j,a,b) = (*D_ov)(i,a) + fock_spin(j,j) - fock_spin(b+nocc,b+nocc);
            (*G_oovv)(i,j,a,b) = ints_mo_spin(i,j,a+nocc,b+nocc);
            (*G_ovov)(i,a,j,b) = ints_mo_spin(i,a+nocc,j,b+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto a = 0; a < nvirt; ++a) {
        for (auto b = 0; b < nvirt; ++b) {
          for (auto c = 0; c < nvirt; ++c) {
            (*G_ovvv)(i,a,b,c) = ints_mo_spin(i,a+nocc,b+nocc,c+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto j = 0; j < nocc; ++j) {
        for (auto k = 0; k < nocc; ++k) {
          for (auto a = 0; a < nvirt; ++a) {
            (*G_ooov)(i,j,k,a) = ints_mo_spin(i,j,k,a+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto j = 0; j < nocc; ++j) {
        for (auto k = 0; k < nocc; ++k) {
          for (auto l = 0; l < nocc; ++l) {
            (*G_oooo)(i,j,k,l) = ints_mo_spin(i,j,k,l);
          }
        }
      }
    }
    for (auto a = 0; a < nvirt; ++a) {
      for (auto b = 0; b < nvirt; ++b) {
        for (auto c = 0; c < nvirt; ++c) {
          for (auto d = 0; d < nvirt; ++d) {
            (*G_vvvv)(a,b,c,d) = ints_mo_spin(a+nocc,b+nocc,c+nocc,d+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto j = 0; j < nocc; ++j) {
        for (auto k = 0; k < nocc; ++k) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              for (auto c = 0; c < nvirt; ++c) {
                (*D_ooovvv)(i,j,k,a,b,c) = (*D_oovv)(i,j,a,b) + (*D_ov)(k,c);
              }
            }
          }
        }
      }
    }

    //
    // singles and doubles for CCSD calculation
    auto t_ov   = std::make_shared<BTensor>(nocc, nvirt);
    auto t_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
    (*t_ov).fill(0.0);
    (*t_oovv).fill(0.0);
    // and triples if CCSDT required
    auto t_ooovvv = std::make_shared<BTensor>(nocc, nocc, nocc, nvirt, nvirt, nvirt);
    (*t_ooovvv).fill(0.0);

    // printing norms
    // cout << "\n"
    //      << "norm(Fock_oo) " << std::sqrt(btas::dot(Fock_oo,Fock_oo)) <<endl;
    // cout << "norm(Fock_ov) " << std::sqrt(btas::dot(Fock_ov,Fock_ov)) <<endl;
    // cout << "norm(Fock_vv) " << std::sqrt(btas::dot(Fock_vv,Fock_vv)) <<endl;
    // cout << "norm(G_oooo)  " << std::sqrt(btas::dot(G_oooo,G_oooo))   <<endl;
    // cout << "norm(G_vvvv)  " << std::sqrt(btas::dot(G_vvvv,G_vvvv))   <<endl;
    // cout << "norm(G_ovvv)  " << std::sqrt(btas::dot(G_ovvv,G_ovvv))   <<endl;
    // cout << "norm(G_ooov)  " << std::sqrt(btas::dot(G_ooov,G_ooov))   <<endl;
    // cout << "norm(G_oovv)  " << std::sqrt(btas::dot(G_oovv,G_oovv))   <<endl;
    // cout << "norm(G_ovov)  " << std::sqrt(btas::dot(G_ovov,G_ovov))   <<endl;
    // cout << "norm(t_ov)    " << std::sqrt(btas::dot(t_ov,t_ov))       <<endl;
    // cout << "norm(t_oovv)  " << std::sqrt(btas::dot(t_oovv,t_oovv))   <<endl;
    //

    // the map that is required while evaluating sequant expressions
    std::map<std::wstring, BTensorPtr> btensor_map;
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_oo",   Fock_oo));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_ov",   Fock_ov));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"f_vv",   Fock_vv));
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_oooo", G_oooo)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_vvvv", G_vvvv)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ovvv", G_ovvv)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ooov", G_ooov)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_oovv", G_oovv)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"g_ovov", G_ovov)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"t_ov",   t_ov)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"t_oovv", t_oovv)); 
    btensor_map.insert(std::pair<std::wstring, BTensorPtr>(L"t_ooovvv", t_ooovvv)); 

    // global sequant setup...
    std::setlocale(LC_ALL, "en_US.UTF-8");
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
    auto cc_r = cceqvec{ 2, 2 }(true, true, true, true, true);

    using sequant::interpret::antisymmetrize;
    using sequant::interpret::eval_equation;
    using sequant::factorize::factorize_expr;

    // factorize CCSD equations
    auto index_size_map = std::make_shared<ispace_map>(ispace_map{});
    index_size_map->insert(ispace_pair(sequant::IndexSpace::active_occupied, nocc));
    index_size_map->insert(ispace_pair(sequant::IndexSpace::active_unoccupied, nvirt));

    for (auto i=1; i<cc_r.size(); ++i)
      cc_r[i] = factorize_expr(cc_r.at(i), index_size_map, true);

    auto cc_r1 = std::make_shared<sequant::evaluate::EvalTensor>(cc_r[1]);
    cc_r1->fill_btas_indices();
    auto cc_r2 = std::make_shared<sequant::evaluate::EvalTensor>(cc_r[2]);
    cc_r2->fill_btas_indices();

    sequant::container::map<sequant::evaluate::hash_type, size_t>  hash_counts_r1, hash_counts_r2;
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
      return sequant::evaluate::hash_to_dtensor_map<BTensor>::value_type(hval, tnsr_ptr);
    };

    auto hash_to_ptr = sequant::evaluate::hash_to_dtensor_map<BTensor>();
    for (auto& item: btensor_map){
      std::wstring label = std::wstring{item.first[0]};
      std::wstring spaces = std::wstring{item.first.substr(2, std::wstring::npos)};
      auto entry = generate_hash_map_entry(label, spaces, item.second);
      hash_to_ptr.insert(entry);
    }

    auto context_r1 = sequant::evaluate::EvalContext<BTensor>(hash_to_ptr, hash_counts_r1);
    auto context_r2 = sequant::evaluate::EvalContext<BTensor>(hash_to_ptr, hash_counts_r2);

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

      R2 = antisymmetrize(R2);
      // R3 = antisymmetrize(R3);

      cout << "using BTAS,    iter " << iter << endl;
      /* cout << "norm(R1) = " << std::sqrt(btas::dot(R1, R1)) << endl; */
      /* cout << "norm(R2) = " << std::sqrt(btas::dot(R2, R2)) << endl; */

      // save previous norms
      auto norm_last  = std::sqrt(btas::dot(*t_oovv, *t_oovv));

      //////////////
      // Updating amplitudes
      // update t_ov
      for (auto i = 0; i < nocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) {
          (*t_ov)(i, a) += R1(i, a)/(*D_ov)(i,a); } }

      //
      // update t_oovv
      for (auto i = 0; i < nocc; ++i) {
        for (auto j = 0; j < nocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              (*t_oovv)(i, j, a, b) += R2(i, j, a, b)/(*D_oovv)(i,j,a,b); } } } }

      //
      // update t_ooovvv
      // for (auto i = 0; i < nocc; ++i) {
      //   for (auto j = 0; j < nocc; ++j) {
      //     for (auto k = 0; k < nocc; ++k) {
      //       for (auto a = 0; a < nvirt; ++a) {
      //         for (auto b = 0; b < nvirt; ++b) {
      //           for (auto c = 0; c < nvirt; ++c) {
      //             t_ooovvv(i,j,k,a,b,c)
      //               += R3(i,j,k,a,b,c)/D_ooovvv(i,j,k,a,b,c); } } } } } }

      cout << "norm(t_ov)    " << std::sqrt(btas::dot(*t_ov,*t_ov))       <<endl;
      cout << "norm(t_oovv)  " << std::sqrt(btas::dot(*t_oovv,*t_oovv))   <<endl;

      auto eccsd_last = eccsd;

      // calculating energy
      BTensor temp_tensor;
      enum {i, j, a, b};
      btas::contract(1.0, *G_oovv, {i,j,a,b}, *t_ov, {i,a}, 0.0, temp_tensor, {j,b});
      //
      eccsd  = 0.5*btas::dot(temp_tensor, *t_ov)
        + 0.25*btas::dot(*G_oovv, *t_oovv)
        + btas::dot(*Fock_ov, *t_ov);
      printf("E(CC) is: %20.12f\n\n", eccsd);

      normdiff = norm_last - sqrt(btas::dot(*t_oovv, *t_oovv));
      ediff    = eccsd_last - eccsd;

    } while ((fabs(normdiff) > conv || fabs(ediff) > conv) && (iter < maxiter));

    auto stop  = high_resolution_clock::now();
    auto  duration = duration_cast<microseconds>(stop - start);
    cout << "\nOut of loop after " << iter   << " iterations.\n"
         << "\nTime: " << duration.count() << " microseconds." << endl;
  } // end of try block; if any exceptions occurred, report them and exit cleanly
  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }
  return 0;
}
