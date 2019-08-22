#include "contract_test.hpp"

using STensor = sequant::Tensor;
using CTensor = sequant::contractable::Tensor;
using BTensor = btas::Tensor<double>;

int main(int argc, char *argv[])
{
  // global setup...
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

  using std::cout;
  using std::cerr;
  using std::wcout;
  using std::endl;
  using std::uniform_real_distribution;
  using std::mt19937;
  using std::vector;
  //
  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  //
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::microseconds;
  auto start = high_resolution_clock::now();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  /* wcout << "Time: " << duration.count() << " microseconds." << endl; */
  /* auto dist = uniform_real_distribution<double>{-1.0, 1.0}; */

  //
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

    const auto use_hcore_guess = true;  // use core Hamiltonian eigenstates to guess density?
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
      const auto tstart = std::chrono::high_resolution_clock::now();
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
    printf("** Hartree-Fock energy = %20.12f\n", hf_final);

    BTensor ao_ints_tensor = compute_ao_ints(shells);

    BTensor mo_ints_tensor_chem = compute_mo_ints(C, ao_ints_tensor);

    // transpose mo_ints_tensor_chem into physicist's notation
    BTensor mo_ints_tensor(nao, nao, nao, nao);
    for (auto i = 0; i < nao; ++i) {
      for (auto j = 0; j < nao; ++j) {
        for (auto k = 0; k < nao; ++k) {
          for (auto l = 0; l < nao; ++l) {
            mo_ints_tensor(i,j,k,l) = mo_ints_tensor_chem(i,k,j,l);
          }
        }
      }
    }

    // transform Fock matrix to MO basis: spin
    BTensor fock_spin(2*nao, 2*nao);
    for (auto l = 0; l < 2*nao; ++l) {
      fock_spin(l,l) = C_v(floor(l/2));
    }

    // transform spatial MO ints tensor to spin basis
    BTensor ints_spin(2*nao, 2*nao, 2*nao, 2*nao);
    for (auto r = 0; r < 2*nao; ++r) {
      for (auto s = 0; s < 2*nao; ++s) {
        for (auto p = 0; p < 2*nao; ++p) {
          for (auto q = 0; q < 2*nao; ++q) {
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
            auto p_i = floor(p/2);
            auto q_i = floor(q/2);
            auto r_i = floor(r/2);
            auto s_i = floor(s/2);
            //
            auto col_int = 0.0;
            auto exc_int = 0.0;
            //
            // orbs: 0 1 2 3 4 ...
            // spin: α β α β α ...
            // 
            //  <01|23>: coulomb integral
            if ((r%2 == p%2) && (s%2 == q%2)){
              col_int  += mo_ints_tensor(r_i, s_i, p_i, q_i);
            }
            //  <01|32>: exchange integral
            if ((r%2 == q%2) && (s%2 == p%2)){
              exc_int  += mo_ints_tensor(r_i, s_i, q_i, p_i);
            }
            ints_spin(r,s,p,q) = col_int - exc_int;
          }
        }
      }
    }

    // from here on, the no. of occupied and unoccupied orbs is doubled
    // for we are using spin orbitals
    ndocc = 2*ndocc;
    auto nvirt = 2*nao - ndocc;

    // D_ov and D_oovv
    BTensor D_ov(ndocc, nvirt), D_oovv(ndocc, ndocc, nvirt, nvirt), D_ooovvv(ndocc, ndocc, ndocc, nvirt, nvirt, nvirt);

    // Fock matrix tensors
    BTensor Fock_oo(ndocc, ndocc), Fock_ov(ndocc, nvirt), Fock_vv(nvirt, nvirt);

    // Integral tensors
    BTensor G_oooo(ndocc, ndocc, ndocc, ndocc), G_vvvv(nvirt, nvirt, nvirt, nvirt);
    BTensor G_ovvv(ndocc, nvirt, nvirt, nvirt), G_ooov(ndocc, ndocc, ndocc, nvirt);
    BTensor G_oovv(ndocc, ndocc, nvirt, nvirt), G_ovov(ndocc, nvirt, ndocc, nvirt);
    for (auto i = 0; i < ndocc; ++i) {
      Fock_oo(i,i) = fock_spin(i,i);
      for (auto a = 0; a < nvirt; ++a) {
        Fock_ov(i,a) = fock_spin(i,a+ndocc);
        D_ov(i,a) = fock_spin(i,i) - fock_spin(a+ndocc,a+ndocc);
        for (auto j = 0; j < ndocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            Fock_vv(a,b) = fock_spin(a+ndocc,b+ndocc);
            D_oovv(i,j,a,b) = D_ov(i,a) + fock_spin(j,j) - fock_spin(b+ndocc,b+ndocc);
            G_oovv(i,j,a,b) = ints_spin(i,j,a+ndocc,b+ndocc);
            G_ovov(i,a,j,b) = ints_spin(i,a+ndocc,j,b+ndocc);
          }
        }
      }
    }
    for (auto i = 0; i < ndocc; ++i) {
      for (auto a = 0; a < nvirt; ++a) {
        for (auto b = 0; b < nvirt; ++b) {
          for (auto c = 0; c < nvirt; ++c) {
            G_ovvv(i,a,b,c) = ints_spin(i,a+ndocc,b+ndocc,c+ndocc);
          }
        }
      }
    }
    for (auto i = 0; i < ndocc; ++i) {
      for (auto j = 0; j < ndocc; ++j) {
        for (auto k = 0; k < ndocc; ++k) {
          for (auto a = 0; a < nvirt; ++a) {
            G_ooov(i,j,k,a) = ints_spin(i,j,k,a+ndocc);
          }
        }
      }
    }
    for (auto i = 0; i < ndocc; ++i) {
      for (auto j = 0; j < ndocc; ++j) {
        for (auto k = 0; k < ndocc; ++k) {
          for (auto l = 0; l < ndocc; ++l) {
            G_oooo(i,j,k,l) = ints_spin(i,j,k,l);
          }
        }
      }
    }
    for (auto a = 0; a < nvirt; ++a) {
      for (auto b = 0; b < nvirt; ++b) {
        for (auto c = 0; c < nvirt; ++c) {
          for (auto d = 0; d < nvirt; ++d) {
            G_vvvv(a,b,c,d) = ints_spin(a+ndocc,b+ndocc,c+ndocc,d+ndocc);
          }
        }
      }
    }
    for (auto i = 0; i < ndocc; ++i) {
      for (auto j = 0; j < ndocc; ++j) {
        for (auto k = 0; k < ndocc; ++k) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              for (auto c = 0; c < nvirt; ++c) {
                D_ooovvv(i,j,k,a,b,c) = D_oovv(i,j,a,b) + D_ov(k,c);
              }
            }
          }
        }
      }
    }

    //
    // singles and doubles for CCSD calculation
    auto t_ov   = BTensor(ndocc, nvirt);
    auto t_oovv = BTensor(ndocc, ndocc, nvirt, nvirt);
    t_ov.fill(0.0);
    t_oovv.fill(0.0);
    // and triples
    auto t_ooovvv = BTensor(ndocc, ndocc, ndocc, nvirt, nvirt, nvirt);
    t_ooovvv.fill(0.0);

    // a map that is required while evaluating sequant expressions
    std::map<std::wstring, BTensor*> btensor_map;
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"f_oo",   &Fock_oo));
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"f_ov",   &Fock_ov));
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"f_vv",   &Fock_vv));
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_oooo", &G_oooo)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_vvvv", &G_vvvv)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_ovvv", &G_ovvv)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_ooov", &G_ooov)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_oovv", &G_oovv)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"g_ovov", &G_ovov)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"t_ov",   &t_ov)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"t_oovv", &t_oovv)); 
    btensor_map.insert(std::pair<std::wstring, BTensor*>(L"t_ooovvv", &t_ooovvv)); 

    //
    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto eccsd = 0.0;
    Logger::get_instance().wick_stats = false;
    auto ccs_r = cceqvec{ 3, 3 }(true, true, true, true);

    start = high_resolution_clock::now();
    do {
      ++iter;
      auto r1 = eval_equation(ccs_r[1], btensor_map).bt();
      auto r2 = eval_equation(ccs_r[2], btensor_map).bt();
      r2 = antisymmetrize(r2);
      cout << endl;
      cout << "norm(r1) = " << std::sqrt(btas::dot(r1, r1)) << endl;
      cout << "norm(r2) = " << std::sqrt(btas::dot(r2, r2)) << endl;
      auto r3 = eval_equation(ccs_r[3], btensor_map).bt();
      r3 = antisymmetrize(r3);
      /* cout << "norm(r3) = " << std::sqrt(btas::dot(r3, r3)) << endl; */

      // save previous norms
      auto norm_last  = std::sqrt(btas::dot(t_oovv, t_oovv));

      //////////////
      // Updating amplitudes
      // update t_ov
      for (auto i = 0; i < ndocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) {
          t_ov(i, a) += r1(i, a)/D_ov(i,a);
        }
      }
      //
      // update t_oovv
      for (auto i = 0; i < ndocc; ++i) {
        for (auto j = 0; j < ndocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              t_oovv(i, j, a, b) += r2(i, j, a, b)/D_oovv(i,j,a,b);
            }
          }
        }
      }
      //
      /* // update t_ooovvv */
      for (auto i = 0; i < ndocc; ++i) {
        for (auto j = 0; j < ndocc; ++j) {
          for (auto k = 0; k < ndocc; ++k) {
            for (auto a = 0; a < nvirt; ++a) {
              for (auto b = 0; b < nvirt; ++b) {
                for (auto c = 0; c < nvirt; ++c) {
                  t_ooovvv(i,j,k,a,b,c) += r3(i,j,k,a,b,c)/D_ooovvv(i,j,k,a,b,c);
                }
              }
            }
          }
        }
      }

      auto eccsd_last = eccsd;

      // calculating energy
      BTensor temp_tensor;
      enum {i, j, a, b};
      btas::contract(1.0, G_oovv, {i,j,a,b}, t_ov, {i,a}, 0.0, temp_tensor, {j,b});
      //
      eccsd  = 0.5*btas::dot(temp_tensor, t_ov)
        + 0.25*btas::dot(G_oovv, t_oovv)
        + btas::dot(Fock_ov, t_ov);
      printf("E(CCSDT) is: %20.12f\n", eccsd);
      //
      normdiff = norm_last - sqrt(btas::dot(t_oovv, t_oovv));
      ediff    = eccsd_last - eccsd;
    } while ((fabs(normdiff) > conv || fabs(ediff) > conv) && (iter < maxiter));
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Out of loop after " << iter << " iterations." << endl;
    cout << endl;
    cout << "Time: " << duration.count() << " microseconds" << endl;
    /* printf("E(CCSD) is: %20.12f\n", eccsd); */
  }
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
