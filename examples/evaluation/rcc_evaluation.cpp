//
// Created by Nakul Teke on 3/18/20.
//
#include "../sequant_setup.hpp"
#include "../../examples/contract/scf/hartree-fock.h"

#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <TiledArray/initialize.h>
#include <tiledarray.h>
#include <libint2.hpp>

#include <memory>

#define CCSDT_eval 0

int main(int argc, char* argv[]) {
  using std::cerr;
  using std::cout;
  using std::endl;

  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  try {
    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take
    // filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename);

    // count the number of electrons
    size_t nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i) nelectron += atoms[i].atomic_number;
    const size_t ndocc = nelectron / 2;

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
    cout << "Nuclear repulsion energy = " << enuc << " a.u." << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

    auto shells = make_sto3g_basis(atoms);
    size_t nao = 0;
    for (auto s = 0; s < shells.size(); ++s) nao += shells[s].size();
    size_t nvirt = nao - ndocc;
    std::cout << "nao  : " << nao << "\t";
    std::cout << "ndocc: " << ndocc << "\t";
    std::cout << "nvirt: " << nvirt << "\n";

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    auto S = compute_1body_ints(shells, Operator::overlap);
    auto T = compute_1body_ints(shells, Operator::kinetic);
    Matrix V = compute_1body_ints(shells, Operator::nuclear, atoms);
    Matrix H = T + V;
    T.resize(0, 0);
    V.resize(0, 0);

    /*** =========================== ***/
    /*** build initial-guess density ***/
    /*** =========================== ***/

    // use core Hamiltonian eigenstates to guess density?
    // set to true to match the result of versions 0, 1, and 2 of the code
    // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
    // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    const auto use_hcore_guess = false;

    Matrix D, C;
    Eigen::VectorXd eps;    // eigenvalues
    if (use_hcore_guess) {  // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    } else {  // SOAD as the guess density, assumes STO-nG basis
      D = compute_soad(atoms);
    }

    /*** =========================== ***/
    /***  SCF iterative loop         ***/
    /*** =========================== ***/

    const auto maxiter = 200;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rmsd = 0.0;
    auto ediff = 0.0;
    auto ehf = 0.0;
    Matrix Fock_matrix;  // capture fock matrix for ccsd calculations
    cout << "\n Iter        E(elec)              E(tot)               "
            "Delta(E)             RMS(D)         Time(s)\n";
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Save a copy of the energy and the density
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
    cout << endl;
    printf("** Hartree-Fock energy = %20.12f a.u.\n", hf_final);

    /*** =========================== ***/
    /***          MP2 energy         ***/
    /*** =========================== ***/

    // Initialize MADWorld
    TA::World& world = TA::initialize(argc, argv);

    auto mo_ints_tensor = compute_mo_ints(C, shells, world);
    auto tile_ints_spatial = mo_ints_tensor.find({0, 0, 0, 0}).get();

    auto TA_fock = [&](const TA::Range& range) {
      TA::Tensor<double> tile(range);
      for (auto i = 0; i < nao; ++i) {
        for (auto j = 0; j < nao; ++j)
          tile(i, j) = (i == j) ? eps(i) : 0.0;
      }
      return tile;
    };

    TA::TArrayD fock(world, TA::TiledRange{{0, nao}, {0, nao}});
    auto tile_fock_temp = fock.world().taskq.add(
        TA_fock, fock.trange().make_tile_range(0));
    *(fock.begin()) = tile_fock_temp;

    cout << "\n"
         << "*********************************\n"
         << "Calculating EMP2 using TiledArray\n"
         << "*********************************\n"
         << endl;

    auto emp2 = 0.0;

    auto tile_ints = mo_ints_tensor.find({0, 0, 0, 0}).get();
    auto tile_fock = fock.find({0, 0}).get();

    for (auto r = 0; r < ndocc; ++r) {
      for (auto s = 0; s < ndocc; ++s) {
        for (auto p = ndocc; p < nao; ++p) {
          for (auto q = ndocc; q < nao; ++q) {
            auto calc = tile_ints(r, s, p, q) *
                   ( 2.0 * tile_ints(r, s, p, q) - tile_ints(r, s, q, p));
            emp2 += calc / (tile_fock(r,r) + tile_fock(s,s) - tile_fock(q,q) - tile_fock(p,p));
          }
        }
      }
    }

    cout << "E(MP2): " << emp2 << "\nFinal energy is : " << hf_final + emp2 << " a.u."
         << endl;

    cout << "\n"
         << "***********************************\n"
         << "Calculating ECCSD using TiledArray\n"
         << "***********************************\n"
         << endl;

    cout << "Initializing tensors..." << endl;

    TA::TiledRange tr_oo{{0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_ov{{0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_vo{{0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_vv{{0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_oooo{{0, ndocc}, {0, ndocc}, {0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_ooov{{0, ndocc}, {0, ndocc}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_oovo{{0, ndocc}, {0, ndocc}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_oovv{{0, ndocc}, {0, ndocc}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_ovov{{0, ndocc}, {0, nvirt}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_ovoo{{0, ndocc}, {0, nvirt}, {0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_ovvo{{0, ndocc}, {0, nvirt}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_vovo{{0, nvirt}, {0, ndocc}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_voov{{0, nvirt}, {0, ndocc}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_ovvv{{0, ndocc}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vvvv{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vovv{{0, nvirt}, {0, ndocc}, {0, nvirt}, {0, nvirt}};

    TA::TiledRange tr_vvoo{{0, nvirt}, {0, nvirt}, {0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_vooo{{0, nvirt}, {0, ndocc}, {0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_vvvo{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_vvov{{0, nvirt}, {0, nvirt}, {0, ndocc}, {0, nvirt}};

#if CCSDT_eval
     TA::TiledRange tr_ooovvv{{0, ndocc},  {0, ndocc},  {0, ndocc},
                              {0, nvirt}, {0, nvirt}, {0, nvirt}};

#endif
    auto D_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    auto D_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
#if CCSDT_eval
    auto D_ooovvv = std::make_shared<TA::TArrayD>(world, tr_ooovvv);
#endif
    auto Fock_oo = std::make_shared<TA::TArrayD>(world, tr_oo);
    auto Fock_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    auto Fock_vo = std::make_shared<TA::TArrayD>(world, tr_vo);
    auto Fock_vv = std::make_shared<TA::TArrayD>(world, tr_vv);
    auto G_oooo = std::make_shared<TA::TArrayD>(world, tr_oooo);
    auto G_ooov = std::make_shared<TA::TArrayD>(world, tr_ooov);
    auto G_oovo = std::make_shared<TA::TArrayD>(world, tr_oovo);
    auto G_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    auto G_ovov = std::make_shared<TA::TArrayD>(world, tr_ovov);
    auto G_ovoo = std::make_shared<TA::TArrayD>(world, tr_ovoo);
    auto G_ovvo = std::make_shared<TA::TArrayD>(world, tr_ovvo);
    auto G_vovo = std::make_shared<TA::TArrayD>(world, tr_vovo);
    auto G_voov = std::make_shared<TA::TArrayD>(world, tr_voov);
    auto G_ovvv = std::make_shared<TA::TArrayD>(world, tr_ovvv);
    auto G_vvvv = std::make_shared<TA::TArrayD>(world, tr_vvvv);
    auto G_vovv = std::make_shared<TA::TArrayD>(world, tr_vovv);

    auto G_vvoo = std::make_shared<TA::TArrayD>(world, tr_vvoo);
    auto G_vooo = std::make_shared<TA::TArrayD>(world, tr_vooo);
    auto G_vvvo = std::make_shared<TA::TArrayD>(world, tr_vvvo);
    auto G_vvov = std::make_shared<TA::TArrayD>(world, tr_vvov);

    (*D_ov).fill(0.0);
    (*D_oovv).fill(0.0);
#if CCSDT_eval
    (*D_ooovvv).fill(0.0);
#endif
    (*Fock_oo).fill(0.0);
    (*Fock_ov).fill(0.0);
    (*Fock_vo).fill(0.0);
    (*Fock_vv).fill(0.0);
    (*G_oooo).fill(0.0);
    (*G_ooov).fill(0.0);
    (*G_oovo).fill(0.0);
    (*G_oovv).fill(0.0);
    (*G_ovov).fill(0.0);
    (*G_ovvo).fill(0.0);
    (*G_ovoo).fill(0.0);
    (*G_vovo).fill(0.0);
    (*G_voov).fill(0.0);
    (*G_ovvv).fill(0.0);
    (*G_vvvv).fill(0.0);
    (*G_vovv).fill(0.0);

    (*G_vvoo).fill(0.0);
    (*G_vooo).fill(0.0);
    (*G_vvvo).fill(0.0);
    (*G_vvov).fill(0.0);

    auto tile_D_ov = (*D_ov).find({0, 0}).get();
    auto tile_D_oovv = (*D_oovv).find({0, 0, 0, 0}).get();
#if CCSDT_eval
    auto tile_D_ooovvv = (*D_ooovvv).find({0, 0, 0, 0, 0, 0}).get();
#endif
    auto tile_Fock_oo = (*Fock_oo).find({0, 0}).get();
    auto tile_Fock_ov = (*Fock_ov).find({0, 0}).get();
    auto tile_Fock_vo = (*Fock_vo).find({0, 0}).get();
    auto tile_Fock_vv = (*Fock_vv).find({0, 0}).get();
    auto tile_G_oooo = (*G_oooo).find({0, 0, 0, 0}).get();
    auto tile_G_ooov = (*G_ooov).find({0, 0, 0, 0}).get();
    auto tile_G_oovo = (*G_oovo).find({0, 0, 0, 0}).get();
    auto tile_G_oovv = (*G_oovv).find({0, 0, 0, 0}).get();
    auto tile_G_ovov = (*G_ovov).find({0, 0, 0, 0}).get();
    auto tile_G_ovoo = (*G_ovoo).find({0, 0, 0, 0}).get();
    auto tile_G_ovvo = (*G_ovvo).find({0, 0, 0, 0}).get();
    auto tile_G_vovo = (*G_vovo).find({0, 0, 0, 0}).get();
    auto tile_G_voov = (*G_voov).find({0, 0, 0, 0}).get();
    auto tile_G_ovvv = (*G_ovvv).find({0, 0, 0, 0}).get();
    auto tile_G_vvvv = (*G_vvvv).find({0, 0, 0, 0}).get();
    auto tile_G_vovv = (*G_vovv).find({0, 0, 0, 0}).get();

    auto tile_G_vvoo = (*G_vvoo).find({0, 0, 0, 0}).get();
    auto tile_G_vooo = (*G_vooo).find({0, 0, 0, 0}).get();
    auto tile_G_vvvo = (*G_vvvo).find({0, 0, 0, 0}).get();
    auto tile_G_vvov = (*G_vvov).find({0, 0, 0, 0}).get();

    for (auto i = 0; i < ndocc; ++i) {
      tile_Fock_oo(i, i) = tile_fock(i, i);
      for (auto a = 0; a < nvirt; ++a) {
        tile_Fock_ov(i, a) = tile_fock(i, a + ndocc);
        tile_Fock_vo(a, i) = tile_fock(a + ndocc, i);

        tile_D_ov(i, a) =
            tile_fock(i, i) - tile_fock(a + ndocc, a + ndocc);
        for (auto j = 0; j < ndocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            tile_Fock_vv(a, b) = tile_fock(a + ndocc, b + ndocc);
            tile_D_oovv(i, j, a, b) = tile_D_ov(i, a) + tile_fock(j, j) -
                tile_fock(b + ndocc, b + ndocc);
            tile_G_oovv(i, j, a, b) = tile_ints_spatial(i, j, a + ndocc, b + ndocc);
            tile_G_ovov(i, a, j, b) = tile_ints_spatial(i, a + ndocc, j, b + ndocc);
            tile_G_ovvo(i, a, b, j) = tile_ints_spatial(i, a + ndocc, b + ndocc, j);
            tile_G_vovo(a, i, b, j) = tile_ints_spatial(a + ndocc, i, b + ndocc, j);
            tile_G_voov(a, i, j, b) = tile_ints_spatial(a + ndocc, i, j, b + ndocc);

            tile_G_vvoo(a, b, i, j) = tile_ints_spatial(a + ndocc,b + ndocc, i, j);
          }
        }
      }
    }

    for (auto i = 0; i < ndocc; ++i)
      for (auto a = 0; a < nvirt; ++a)
        for (auto b = 0; b < nvirt; ++b)
          for (auto c = 0; c < nvirt; ++c){
            tile_G_ovvv(i, a, b, c) =
                tile_ints_spatial(i, a + ndocc, b + ndocc, c + ndocc);
            tile_G_vovv(a, i, b, c) =
                tile_ints_spatial(a + ndocc, i, b + ndocc, c + ndocc);

            tile_G_vvvo(a, b, c, i) = tile_ints_spatial(a + ndocc, b + ndocc, c + ndocc, i);
            tile_G_vvov(a, b, i, c) = tile_ints_spatial(a + ndocc, b + ndocc, i, c + ndocc);
          }

    for (auto i = 0; i < ndocc; ++i)
      for (auto j = 0; j < ndocc; ++j)
        for (auto k = 0; k < ndocc; ++k)
          for (auto a = 0; a < nvirt; ++a){
            tile_G_ooov(i, j, k, a) = tile_ints_spatial(i, j, k, a + ndocc);
            tile_G_oovo(i, j, a, k) = tile_ints_spatial(i, j, a + ndocc, k);

            tile_G_ovoo(i, a, j, k) = tile_ints_spatial(i, a + ndocc, j, k);
            tile_G_vooo(a, i, j, k) = tile_ints_spatial(a + ndocc, i, j, k);
          }

    for (auto i = 0; i < ndocc; ++i)
      for (auto j = 0; j < ndocc; ++j)
        for (auto k = 0; k < ndocc; ++k)
          for (auto l = 0; l < ndocc; ++l)
            tile_G_oooo(i, j, k, l) = tile_ints_spatial(i, j, k, l);

    for (auto a = 0; a < nvirt; ++a)
      for (auto b = 0; b < nvirt; ++b)
        for (auto c = 0; c < nvirt; ++c)
          for (auto d = 0; d < nvirt; ++d)
            tile_G_vvvv(a, b, c, d) =
                tile_ints_spatial(a + ndocc, b + ndocc, c + ndocc, d + ndocc);

#if CCSDT_eval
    for (auto i = 0; i < ndocc; ++i)
      for (auto j = 0; j < ndocc; ++j)
        for (auto k = 0; k < ndocc; ++k)
          for (auto a = 0; a < nvirt; ++a)
            for (auto b = 0; b < nvirt; ++b)
              for (auto c = 0; c < nvirt; ++c)
                tile_D_ooovvv(i, j, k, a, b, c) =
                    tile_D_oovv(i, j, a, b) + tile_D_ov(k, c);
#endif

    //
    // amplitudes for coupled-cluster calculations
    auto t_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    (*t_ov).fill(0.);
    auto t_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    (*t_oovv).fill(0.);
    auto t_vvoo = std::make_shared<TA::TArrayD>(world, tr_vvoo);
    (*t_vvoo).fill(0.);

#if CCSDT_eval
    auto t_ooovvv = std::make_shared<TA::TArrayD>(world, tr_ooovvv);
    (*t_ooovvv).fill(0.0);
#endif

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

#if CCSDT_eval // CCSDT
    auto cc_r = cceqvec{3, 3}(true, true, true, true);
#else // CCSD
    auto cc_r = cceqvec{2, 2}(true, true, true, true);
#endif

    // SPIN TRACE THE RESIDUAL
    std::vector<ExprPtr> cc_st_r(cc_r.size());
    for (size_t i = 1; i < cc_r.size(); ++i){
      const auto tstart = std::chrono::high_resolution_clock::now();
      std::initializer_list<IndexList> external_indices = {{}};
      if(i == 1)
        external_indices = {{L"i_1", L"a_1"}};
      else if(i == 2)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
      else if(i == 3)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};

      cc_st_r[i] = spintrace(cc_r[i], external_indices);
      canonicalize(cc_st_r[i]);
      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;
      printf("CC R%lu size: %lu time: %5.3f sec.\n", i, cc_st_r[i]->size(), time_elapsed.count());
    }

    // Use Biorthogonal transformation for simpler CCSD residual equations
    bool biorthogonal_transformation = true;
    if(biorthogonal_transformation){
      // CCSD R1
      auto cc_r1_biorthogonal = ex<Constant>(0.5) * cc_st_r[1];
      expand(cc_r1_biorthogonal);
      rapid_simplify(cc_r1_biorthogonal);
      cc_st_r[1] = cc_r1_biorthogonal;

      // CCSD R2: 1/3 R2 + 1/6 R2' for simpler equations
      std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                       {Index{L"i_2"}, Index{L"i_1"}}};

      auto temp_expr = transform_expression(cc_st_r[2], idxmap);
      auto biorthogonal_R2 =
              ex<Constant>(1.0 / 3.0) * cc_st_r[2] +
              ex<Constant>(1.0 / 6.0) * temp_expr;

      auto P_cc_r2 = ex<Constant>(0.5) *
          ex<Tensor>(L"P", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}, Symmetry::nonsymm) *
          biorthogonal_R2;
      expand(P_cc_r2);
      canonicalize(P_cc_r2);

      // TODO: Do not expand P operator; evaluate R2 and permute instead
      auto cc_r2_biorthogonal = expand_P_operator(P_cc_r2);
      rapid_simplify(cc_r2_biorthogonal);
      cc_st_r[2] = cc_r2_biorthogonal;
    } else {
      expand(cc_st_r[1]);
      rapid_simplify(cc_st_r[1]);
      expand(cc_st_r[2]);
      rapid_simplify(cc_st_r[2]);
    }

#if CCSDT_eval
    if(biorthogonal_transformation){
      // Biorthogonal transformation of CC R3 equation (http://arxiv.org/abs/1805.00565 Eq. 31)

      using IdxMap = std::map<Index, Index>;
      IdxMap m1 = {{Index{L"i_1"}, Index{L"i_1"}},
                   {Index{L"i_2"}, Index{L"i_3"}},
                   {Index{L"i_3"}, Index{L"i_2"}}}; // ijk -> ikj
      IdxMap m2 = {{Index{L"i_1"}, Index{L"i_2"}},
                   {Index{L"i_2"}, Index{L"i_1"}},
                   {Index{L"i_3"}, Index{L"i_3"}}}; // ijk -> jik
      IdxMap m3 = {{Index{L"i_1"}, Index{L"i_2"}},
                   {Index{L"i_2"}, Index{L"i_3"}},
                   {Index{L"i_3"}, Index{L"i_1"}}}; // ijk -> jki
      IdxMap m4 = {{Index{L"i_1"}, Index{L"i_3"}},
                   {Index{L"i_2"}, Index{L"i_1"}},
                   {Index{L"i_3"}, Index{L"i_2"}}}; // ijk -> kij
      IdxMap m5 = {{Index{L"i_1"}, Index{L"i_3"}},
                   {Index{L"i_2"}, Index{L"i_2"}},
                   {Index{L"i_3"}, Index{L"i_1"}}}; // ijk -> kji

      auto p1 = transform_expression(cc_st_r[3], m1,-1.0/120.0);
      auto p2 = transform_expression(cc_st_r[3], m2,-1.0/120.0);
      auto p3 = transform_expression(cc_st_r[3], m3,-7.0/120.0);
      auto p4 = transform_expression(cc_st_r[3], m4,-7.0/120.0);
      auto p5 = transform_expression(cc_st_r[3], m5,-1.0/120.0);

      auto biorthogonal_R3 =
          ex<Constant>(17/120) * cc_st_r[3] + p1 + p2 + p3 + p4 + p5;

      auto P_cc_r3 = ex<Constant>(1./6.) *
          ex<Tensor>(L"P", WstrList{L"a_1", L"a_2", L"a_3"}, WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::nonsymm) *
          biorthogonal_R3;
      expand(P_cc_r3);
      canonicalize(P_cc_r3);

      // TODO: Do not expand P operator; evaluate R3 and permute instead
      auto cc_r3_biorthogonal = expand_P_operator(P_cc_r3);
      rapid_simplify(cc_r3_biorthogonal);
      cc_st_r[3] = cc_r3_biorthogonal;
    } else {
      expand(cc_st_r[3]);
      rapid_simplify(cc_st_r[3]);
    }
#endif

    std::vector<std::shared_ptr<TA::TArrayD>> data_tensors = {Fock_oo, Fock_ov, Fock_vo, Fock_vv,
              G_oooo, G_ooov, G_oovo, G_oovv, G_ovov, G_ovvo, G_vovo, G_voov, G_ovvv, G_vovv, G_vvvv,
              G_vvoo, G_vooo, G_vvvo, G_vvov, G_ovoo,
              t_ov, t_oovv, t_vvoo};

    std::vector<sequant::ExprPtr> seq_tensors = {
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1"}, {L"i_1"})), // f_oo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1",}, {L"a_1"})), // f_ov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"a_1"}, {L"i_1",})), // f_vo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"a_1",}, {L"a_2",})), // f_vv

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"i_4"})), //oooo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"a_1"})), //ooov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"a_1",L"i_3"})), //oovo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"a_1",L"a_2"})), //oovv
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"i_2",L"a_2"})), //ovov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"a_2",L"i_2"})), //ovvo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"i_1"}, {L"a_2",L"i_2"})), //vovo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"i_1"}, {L"i_2",L"a_2"})), //voov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"a_2",L"a_3"})), //ovvv
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1", L"i_1"}, {L"a_2",L"a_3"})), //vovv
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"a_3",L"a_4"})), //vvvv

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"i_1",L"i_2"})), //vvoo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"i_1"}, {L"i_2",L"i_3"})), //vooo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"a_3",L"i_1"})), //vvvo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"i_1",L"a_3"})), //vvov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"i_2",L"i_3"})), //ovoo

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",}, {L"a_1"})), //ov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",L"i_2"}, {L"a_1",L"a_2"})), //oovv
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"a_1",L"a_2"}, {L"i_1",L"i_2"})) //vvoo
    };

#if CCSDT_eval
    data_tensors.push_back(t_ooovvv);
    seq_tensors.push_back(std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1", L"i_2", L"i_3"}, {L"a_1", L"a_2", L"a_3"})));
#endif

    using sequant::evaluate::EvalTree;
    using ContextMapType = container::map<sequant::evaluate::HashType, 
          std::shared_ptr<TA::TArrayD>>;
    // whether to consider T_{i j}^{a b} and T_{a b}^{i j} kind of tensors
    // equivalent or not.
    bool swap_braket_labels = false;

    ContextMapType context;

    assert(data_tensors.size() == seq_tensors.size());
    for (auto i = 0; i < seq_tensors.size(); ++i)
        context.insert(ContextMapType::value_type(
                EvalTree(seq_tensors.at(i), swap_braket_labels).hash_value(),
                data_tensors.at(i)));


    auto r1_tree = EvalTree(cc_st_r[1], swap_braket_labels);
    auto r2_tree = EvalTree(cc_st_r[2], swap_braket_labels);
#if CCSDT_eval
    auto r3_tree = EvalTree(cc_st_r[3], swap_braket_labels);
#endif

    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto ecc = 0.0;
    cout << "Using TiledArray..." << endl;
    cout << "Iter   norm(t_ov)    norm(t_oovv)     ΔE(CC)          E(CC)       time(s)" << endl;
    cout << "============================================================================" << endl;
    auto CC_start = std::chrono::high_resolution_clock::now();
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;
      auto R1 = r1_tree.evaluate(context);
      auto R2_ = r2_tree.evaluate(context);

      TA::TArrayD R2;
      R2("i,j,a,b") = R2_("a,b,i,j");

#if CCSDT_eval
      auto R3 = r3_tree.evaluate(context);
#endif

      auto tile_R1       = R1.find({0,0}).get();
      auto tile_t_ov     = (*t_ov).find({0,0}).get();

      auto tile_R2       = R2.find({0,0,0,0}).get();
      auto tile_t_oovv   = (*t_oovv).find({0,0,0,0}).get();
      auto tile_t_vvoo   = (*t_vvoo).find({0,0,0,0}).get();

#if CCSDT_eval
      auto tile_R3       = R3.find({0,0,0,0,0,0}).get();
      auto tile_t_ooovvv   = (*t_ooovvv).find({0,0,0,0,0,0}).get();
#endif

      // save previous norm
      auto norm_last = std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));

      // Updating amplitudes

      // update t_ov
      for (auto i = 0; i < ndocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) { // TODO: Derive equation for regular T1
          tile_t_ov(i,a) += tile_R1(i,a)/tile_D_ov(i,a); } }

      // update t_oovv
      for (auto i = 0; i < ndocc; ++i) {
        for (auto j = 0; j < ndocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) { // TODO: Derive equation for regular T2
              tile_t_oovv(i,j,a,b) += tile_R2(i,j,a,b)/tile_D_oovv(i,j,a,b); } } } }

# if CCSDT_eval
      // update t_ooovvv
      for (auto i = 0; i < ndocc; ++i) {
        for (auto j = 0; j < ndocc; ++j) {
          for (auto k = 0; k < ndocc; ++k) {
            for (auto a = 0; a < nvirt; ++a) {
              for (auto b = 0; b < nvirt; ++b) {
                for (auto c = 0; c < nvirt; ++c) {
                  tile_t_ooovvv(i,j,k,a,b,c) +=
                      tile_R3(i,j,k,a,b,c)/tile_D_ooovvv(i,j,k,a,b,c); } } } } } }
#endif

      auto ecc_last = ecc;

      // Calculate CCSD contribution to correlation energy
      {
        TA::TArrayD temp_tensor1, temp_tensor2;
        temp_tensor1("j,b") = (*G_oovv)("i,j,a,b")*(*t_ov)("i,a");
        temp_tensor2("j,b") = (*G_oovv)("i,j,b,a")*(*t_ov)("i,a");

        // Spin-free energy equation
        ecc = (*Fock_ov)("i,a").dot((*t_ov)("i,a")) +
            temp_tensor1("j,b").dot((*t_ov)("j,b")) +
            (*G_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b"));
        ecc *= 2.0;
        ecc -= temp_tensor2("j,b").dot((*t_ov)("j,b"));
        ecc -= (*G_oovv)("i,j,a,b").dot((*t_oovv)("j,i,a,b"));
      }

      // convergence variables
      auto norm_t1 = std::sqrt((*t_ov)("i,j").dot((*t_ov)("i,j")));
      auto norm_t2 = std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));
      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;
      printf("%2d    %4.8f     %4.8f     %4.8f     %4.12f   %5.5f\n", iter,
              norm_t1, norm_t2, ediff, ecc, time_elapsed.count());
      normdiff = norm_last - norm_t2;
      ediff    = ecc_last - ecc;
    } while((fabs(normdiff) > conv || fabs(ediff) > conv) && (iter < maxiter));

    TA::finalize();

    auto CC_stop = std::chrono::high_resolution_clock::now();
    auto time_elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(CC_stop - CC_start);

    cout << "\nOut of loop after "
        << iter << " iterations.\n"
        << "Time: "
        << time_elapsed.count() << " μs"
        << endl;
    cout << "Total energy = " << enuc + ehf + ecc << " a.u." << endl;

    { // Check water molecule, sto-3g
      double ccsd_correlation = -0.070680451962;
      double ccsdt_correlation = -0.070813170670;
#if CCSDT_eval
      assert(fabs(ccsdt_correlation - ecc) < 1e-10);
#else
      assert(fabs(ccsd_correlation - ecc) < 1e-10);
#endif
    }

  }  // end of try block

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  } catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  } catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  } catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }
  return 0;
}

