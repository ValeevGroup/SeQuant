#include <iostream>
#ifndef SEQUANT_HAS_TILEDARRAY
# error "SEQUANT_HAS_TILEDARRAY should be defined when building cc_tiledarray"
#endif

#include "sequant_setup.hpp"

#include "../contract/scf/hartree-fock.h"
#include "../contract/interpret/interpreted_tensor.hpp"
#include "../contract/interpret/contract.hpp"

int main(int argc, char *argv[])
{

  using std::cout;
  using std::cerr;
  using std::endl;

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
    size_t nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const size_t ndocc = nelectron / 2;

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
    size_t nvirt = nao - ndocc;

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

    // use core Hamiltonian eigenstates to guess density?
    // set to true to match the result of versions 0, 1, and 2 of the code
    // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
    // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    const auto use_hcore_guess = false; 

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
        cout <<
        "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n";
      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf, ehf + enuc,
             ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    libint2::finalize(); // done with libint
     
    auto hf_final = ehf + enuc;
    printf("** Hartree-Fock energy = %20.12f\n", hf_final);

    /*** =========================== ***/
    /***          MP2 energy         ***/
    /*** =========================== ***/
    //

    /* cout << "Coeff matrix is: " << endl; */
    /* cout << C << endl << endl; */

    // Initialize MADWorld
    TA::World& world = TA::initialize(argc, argv);

    cout << endl;
    cout << "*********************************" << endl;
    cout << "Calculating EMP2 using TiledArray" << endl;
    cout << "*********************************" << endl;

    auto mo_ints_tensor     = compute_mo_ints(C, shells, world);
    auto tile_ints_spatial  = mo_ints_tensor.find({0,0,0,0}).get();

    // from here on number of occupied and virtual orbitals
    // is doubled as we are now working on spin basis
    nao         *= 2;
    nvirt       *= 2;
    size_t nocc  = 2*ndocc; // in place of ndocc use nocc

    auto spinify_ints = [&](const TA::Range& range) {
      TA::Tensor<double> tile(range);
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
                col_int  += tile_ints_spatial(r_i, s_i, p_i, q_i);

              if ((r%2 == q%2) && (s%2 == p%2))
                exc_int  += tile_ints_spatial(r_i, s_i, q_i, p_i);

              tile(r,s,p,q) = col_int - exc_int;
            }
          }
        }
      }
      return tile;
    };

    TA::TArrayD mo_ints_tensor_spin(world,
                                    TA::TiledRange{{0, nao}, {0, nao},
                                    {0, nao}, {0, nao}});

    auto tile_spin_ints = mo_ints_tensor_spin.world().taskq.add(spinify_ints,
                mo_ints_tensor_spin.trange().make_tile_range(0));
    *(mo_ints_tensor_spin.begin()) = tile_spin_ints;

    auto spinify_fock = [&](const TA::Range& range){
      TA::Tensor<double> tile(range);
      for (auto i = 0; i < nao; ++i){
        for (auto j = 0; j < nao; ++j)
          tile(i,j) = (i == j)? C_v(floor(i/2)): 0;
      }
      return tile;
    };

    TA::TArrayD fock_spin(world, TA::TiledRange{{0, nao}, {0, nao}});
    auto tile_fock = fock_spin.world().taskq.add(spinify_fock,
                           fock_spin.trange().make_tile_range(0));
    *(fock_spin.begin()) = tile_fock;

    // See Step 4 at
    // http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4
    //

    auto tile_ints_spin = mo_ints_tensor_spin.find({0,0,0,0}).get();
    auto tile_fock_spin = fock_spin.find({0,0}).get();
    auto emp2 = 0.0;

    for (auto r = 0; r < nocc; ++r) {
        for (auto s = 0; s < nocc; ++s) {
            for (auto p = nocc; p < nao; ++p) {
                for (auto q = nocc; q < nao; ++q) {

                  auto calc = tile_ints_spin(r,s,p,q);

                  calc *= calc;

                  emp2 += calc/(tile_fock_spin(r,r)
                      + tile_fock_spin(s,s)
                      - tile_fock_spin(p,p)
                      - tile_fock_spin(q,q));
                }
            }
        }
    }
    emp2 /= 4.0;
    cout << "\nThe MP2 correlation energy: " << emp2 << endl;
    cout << "Final energy is : " << hf_final + emp2 << endl;

    TA::TiledRange   tr_oo{{0, nocc},  {0, nocc}};
    TA::TiledRange   tr_ov{{0, nocc},  {0, nvirt}};
    TA::TiledRange   tr_vv{{0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_oooo{{0, nocc},  {0, nocc},  {0, nocc},  {0, nocc}};
    TA::TiledRange tr_ooov{{0, nocc},  {0, nocc},  {0, nocc},  {0, nvirt}};
    TA::TiledRange tr_oovv{{0, nocc},  {0, nocc},  {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_ovov{{0, nocc},  {0, nvirt}, {0, nocc},  {0, nvirt}};
    TA::TiledRange tr_ovvv{{0, nocc},  {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vvvv{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    cout << "\nCout 1" << endl;
    //
    auto D_ov    = TA::TArrayD(world, tr_ov);
    auto D_oovv  = TA::TArrayD(world, tr_oovv);
    auto Fock_oo = TA::TArrayD(world, tr_oo);
    auto Fock_ov = TA::TArrayD(world, tr_ov);
    auto Fock_vv = TA::TArrayD(world, tr_vv);
    auto G_oooo  = TA::TArrayD(world, tr_oooo);
    auto G_ooov  = TA::TArrayD(world, tr_ooov);
    auto G_oovv  = TA::TArrayD(world, tr_oovv);
    auto G_ovov  = TA::TArrayD(world, tr_ovov);
    auto G_ovvv  = TA::TArrayD(world, tr_ovvv);
    auto G_vvvv  = TA::TArrayD(world, tr_vvvv);
    cout << "\nCout 2" << endl;

    auto tile_D_ov    = TA::TArrayD::value_type(   D_ov.trange().make_tile_range(0));
    auto tile_D_oovv  = TA::TArrayD::value_type( D_oovv.trange().make_tile_range(0));
    auto tile_Fock_oo = TA::TArrayD::value_type(Fock_oo.trange().make_tile_range(0));
    auto tile_Fock_ov = TA::TArrayD::value_type(Fock_ov.trange().make_tile_range(0));
    auto tile_Fock_vv = TA::TArrayD::value_type(Fock_vv.trange().make_tile_range(0));
    auto tile_G_oooo  = TA::TArrayD::value_type( G_oooo.trange().make_tile_range(0));
    auto tile_G_ooov  = TA::TArrayD::value_type( G_ooov.trange().make_tile_range(0));
    auto tile_G_oovv  = TA::TArrayD::value_type( G_oovv.trange().make_tile_range(0));
    auto tile_G_ovov  = TA::TArrayD::value_type( G_ovov.trange().make_tile_range(0));
    auto tile_G_ovvv  = TA::TArrayD::value_type( G_ovvv.trange().make_tile_range(0));
    auto tile_G_vvvv  = TA::TArrayD::value_type( G_vvvv.trange().make_tile_range(0));
    cout << "\nCout 3" << endl;

    for (auto i = 0; i < nocc; ++i) {
      tile_Fock_oo(i,i) = tile_fock_spin(i,i);
      for (auto a = 0; a < nvirt; ++a) {
        tile_Fock_ov(i,a) = tile_fock_spin(i,a+nocc);
        tile_D_ov(i,a) = tile_fock_spin(i,i) - tile_fock_spin(a+nocc,a+nocc);
        for (auto j = 0; j < nocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            tile_Fock_vv(a,b) = tile_fock_spin(a+nocc,b+nocc);
            tile_D_oovv(i,j,a,b) = tile_D_ov(i,a) + tile_fock_spin(j,j)
                                                - tile_fock_spin(b+nocc,b+nocc);
            tile_G_oovv(i,j,a,b) = tile_ints_spin(i,j,a+nocc,b+nocc);
            tile_G_ovov(i,a,j,b) = tile_ints_spin(i,a+nocc,j,b+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto a = 0; a < nvirt; ++a) {
        for (auto b = 0; b < nvirt; ++b) {
          for (auto c = 0; c < nvirt; ++c) {
            tile_G_ovvv(i,a,b,c) = tile_ints_spin(i,a+nocc,b+nocc,c+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto j = 0; j < nocc; ++j) {
        for (auto k = 0; k < nocc; ++k) {
          for (auto a = 0; a < nvirt; ++a) {
            tile_G_ooov(i,j,k,a) = tile_ints_spin(i,j,k,a+nocc);
          }
        }
      }
    }
    for (auto i = 0; i < nocc; ++i) {
      for (auto j = 0; j < nocc; ++j) {
        for (auto k = 0; k < nocc; ++k) {
          for (auto l = 0; l < nocc; ++l) {
            tile_G_oooo(i,j,k,l) = tile_ints_spin(i,j,k,l);
          }
        }
      }
    }
    for (auto a = 0; a < nvirt; ++a) {
      for (auto b = 0; b < nvirt; ++b) {
        for (auto c = 0; c < nvirt; ++c) {
          for (auto d = 0; d < nvirt; ++d) {
            tile_G_vvvv(a,b,c,d) = tile_ints_spin(a+nocc,b+nocc,c+nocc,d+nocc);
          }
        }
      }
    }

    cout << "\nCout 4" << endl;

    *D_ov.begin()    = tile_D_ov   ;
    *D_oovv.begin()  = tile_D_oovv ;
    *Fock_oo.begin() = tile_Fock_oo;
    *Fock_ov.begin() = tile_Fock_ov;
    *Fock_vv.begin() = tile_Fock_vv;
    *G_oooo.begin()  = tile_G_oooo ;
    *G_ooov.begin()  = tile_G_ooov ;
    *G_oovv.begin()  = tile_G_oovv ;
    *G_ovov.begin()  = tile_G_ovov ;
    *G_ovvv.begin()  = tile_G_ovvv ;
    *G_vvvv.begin()  = tile_G_vvvv ;

    cout << "\nCout 4" << endl;

    // singles and doubles for CCSD calculations
    auto t_ov   = TA::TArrayD(world, tr_ov);
    auto t_oovv = TA::TArrayD(world, tr_oovv);
    t_ov.fill(0.0);
    t_oovv.fill(0.0);

    cout << "\nCout 5" << endl;

    // a map that is required while evaluating sequant expressions
    using pair_tensor_map = std::pair<std::wstring, TA::TArrayD const *>;
    auto tensor_map = std::map<pair_tensor_map::first_type,
                               pair_tensor_map::second_type>{};
    tensor_map.insert(pair_tensor_map(L"f_oo",   &Fock_oo));
    tensor_map.insert(pair_tensor_map(L"f_ov",   &Fock_ov));
    tensor_map.insert(pair_tensor_map(L"f_vv",   &Fock_vv));
    tensor_map.insert(pair_tensor_map(L"g_oooo", &G_oooo)); 
    tensor_map.insert(pair_tensor_map(L"g_vvvv", &G_vvvv)); 
    tensor_map.insert(pair_tensor_map(L"g_ovvv", &G_ovvv)); 
    tensor_map.insert(pair_tensor_map(L"g_ooov", &G_ooov)); 
    tensor_map.insert(pair_tensor_map(L"g_oovv", &G_oovv)); 
    tensor_map.insert(pair_tensor_map(L"g_ovov", &G_ovov)); 
    tensor_map.insert(pair_tensor_map(L"t_ov",   &t_ov)); 
    tensor_map.insert(pair_tensor_map(L"t_oovv", &t_oovv)); 

    cout << "\nCout 6" << endl;

    //
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
    Logger::get_instance().wick_stats = false;

    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto eccsd = 0.0;
    Logger::get_instance().wick_stats = false;
    auto ccs_r = cceqvec{ 2, 2 }(true, true, true, true);

    cout << "\nCout 7" << endl;

    using sequant::interpret::antisymmetrize;
    using sequant::interpret::eval_equation;

    auto start = std::chrono::high_resolution_clock::now();

    cout << "\nCout 7" << endl;

    do {
      ++iter;
      auto r1 = eval_equation(ccs_r[1], tensor_map);
      cout << "\nCout 8" << endl;
      auto r2 = eval_equation(ccs_r[2], tensor_map);
      cout << "\nCout 9" << endl;
      cout << "norm(r1) = " << r1.tensor()("i,j").norm() << endl;
      cout << "norm(r2) = " << r2.tensor()("i,j,k,l").norm() << endl;
      /* antisymmetrize(r2, r2.tensor().trange().rank()); */

      cout << endl;
    } while(false);

    TA::finalize();
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
