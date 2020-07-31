//
// Created by Nakul Teke on 3/18/20.
//
#include "../../examples/contract/scf/hartree-fock.h"
#include "../sequant_setup.hpp"

#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <TiledArray/initialize.h>
#include <tiledarray.h>
#include <libint2.hpp>

#include <memory>

#define CCSDT_eval 0

container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold);
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const std::initializer_list<IndexList> ext_index_groups);

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

    const auto maxiter = 100;
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
         << "Coupled Cluster using TiledArray\n"
         << "***********************************\n"
         << endl;

    cout << "Initializing tensors..." << endl;

    TA::TiledRange tr_oo{{0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_ov{{0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_vv{{0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_oooo{{0, ndocc}, {0, ndocc}, {0, ndocc}, {0, ndocc}};
    TA::TiledRange tr_ooov{{0, ndocc}, {0, ndocc}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_oovo{{0, ndocc}, {0, ndocc}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_oovv{{0, ndocc}, {0, ndocc}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_ovov{{0, ndocc}, {0, nvirt}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_ovvo{{0, ndocc}, {0, nvirt}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_vovo{{0, nvirt}, {0, ndocc}, {0, nvirt}, {0, ndocc}};
    TA::TiledRange tr_voov{{0, nvirt}, {0, ndocc}, {0, ndocc}, {0, nvirt}};
    TA::TiledRange tr_ovvv{{0, ndocc}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vvvv{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vovv{{0, nvirt}, {0, ndocc}, {0, nvirt}, {0, nvirt}};

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
    auto Fock_vv = std::make_shared<TA::TArrayD>(world, tr_vv);
    auto G_oooo = std::make_shared<TA::TArrayD>(world, tr_oooo);
    auto G_ooov = std::make_shared<TA::TArrayD>(world, tr_ooov);
    auto G_oovo = std::make_shared<TA::TArrayD>(world, tr_oovo);
    auto G_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    auto G_ovov = std::make_shared<TA::TArrayD>(world, tr_ovov);
    auto G_ovvo = std::make_shared<TA::TArrayD>(world, tr_ovvo);
    auto G_vovo = std::make_shared<TA::TArrayD>(world, tr_vovo);
    auto G_voov = std::make_shared<TA::TArrayD>(world, tr_voov);
    auto G_ovvv = std::make_shared<TA::TArrayD>(world, tr_ovvv);
    auto G_vvvv = std::make_shared<TA::TArrayD>(world, tr_vvvv);
    auto G_vovv = std::make_shared<TA::TArrayD>(world, tr_vovv);

    (*D_ov).fill(0.0);
    (*D_oovv).fill(0.0);
#if CCSDT_eval
    (*D_ooovvv).fill(0.0);
#endif
    (*Fock_oo).fill(0.0);
    (*Fock_ov).fill(0.0);
    (*Fock_vv).fill(0.0);
    (*G_oooo).fill(0.0);
    (*G_ooov).fill(0.0);
    (*G_oovo).fill(0.0);
    (*G_oovv).fill(0.0);
    (*G_ovov).fill(0.0);
    (*G_ovvo).fill(0.0);
    (*G_vovo).fill(0.0);
    (*G_voov).fill(0.0);
    (*G_ovvv).fill(0.0);
    (*G_vvvv).fill(0.0);
    (*G_vovv).fill(0.0);

    auto tile_D_ov = (*D_ov).find({0, 0}).get();
    auto tile_D_oovv = (*D_oovv).find({0, 0, 0, 0}).get();
#if CCSDT_eval
    auto tile_D_ooovvv = (*D_ooovvv).find({0, 0, 0, 0, 0, 0}).get();
#endif
    auto tile_Fock_oo = (*Fock_oo).find({0, 0}).get();
    auto tile_Fock_ov = (*Fock_ov).find({0, 0}).get();
    auto tile_Fock_vv = (*Fock_vv).find({0, 0}).get();
    auto tile_G_oooo = (*G_oooo).find({0, 0, 0, 0}).get();
    auto tile_G_ooov = (*G_ooov).find({0, 0, 0, 0}).get();
    auto tile_G_oovo = (*G_oovo).find({0, 0, 0, 0}).get();
    auto tile_G_oovv = (*G_oovv).find({0, 0, 0, 0}).get();
    auto tile_G_ovov = (*G_ovov).find({0, 0, 0, 0}).get();
    auto tile_G_ovvo = (*G_ovvo).find({0, 0, 0, 0}).get();
    auto tile_G_vovo = (*G_vovo).find({0, 0, 0, 0}).get();
    auto tile_G_voov = (*G_voov).find({0, 0, 0, 0}).get();
    auto tile_G_ovvv = (*G_ovvv).find({0, 0, 0, 0}).get();
    auto tile_G_vvvv = (*G_vvvv).find({0, 0, 0, 0}).get();
    auto tile_G_vovv = (*G_vovv).find({0, 0, 0, 0}).get();

    for (auto i = 0; i < ndocc; ++i) {
      tile_Fock_oo(i, i) = tile_fock(i, i);
      for (auto a = 0; a < nvirt; ++a) {
        tile_Fock_ov(i, a) = tile_fock(i, a + ndocc);
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
          }

    for (auto i = 0; i < ndocc; ++i)
      for (auto j = 0; j < ndocc; ++j)
        for (auto k = 0; k < ndocc; ++k)
          for (auto a = 0; a < nvirt; ++a){
            tile_G_ooov(i, j, k, a) = tile_ints_spatial(i, j, k, a + ndocc);
            tile_G_oovo(i, j, a, k) = tile_ints_spatial(i, j, a + ndocc, k);
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

    cout << "\n"
         << "***********************************\n";
#if CCSDT_eval // CCSDT
    cout << "            CCSDT\n";
    auto cc_r = cceqvec{3, 3}(true, true, true, true);
    // auto cc_r = cceqvec{4, 4}(true, true, true, true);
#else // CCSD
    cout << "            CCSD\n";
    auto cc_r = cceqvec{2, 2}(true, true, true, true);
    // cout << "            CCSDTQ\n";
    // auto cc_r = cceqvec{4, 4}(true, true, true, true);
#endif
    cout << "***********************************\n" << endl;

//    bool expand_S = false;
//    bool bt = false;
//    bool factorize_S = false;

    bool expand_S = true;
    bool bt = true;
    bool factorize_S = true;

    // std::wcout << to_latex(cc_r[1]) << std::endl;
    // std::wcout << to_latex(cc_r[2]) << std::endl;

    // SPIN TRACE THE CC RESIDUAL EQUATIONS
    std::vector<ExprPtr> cc_st_r(cc_r.size());
    for (int i = 1; i < cc_r.size(); ++i){
      // std::wcout << "CC R" << i << ":\n" << to_latex_align(cc_r[i], 20, 5) << "\n";
      std::wcout << "CC R" << i << ":\n" << to_latex(cc_r[i]) << "\n";
      const auto tstart = std::chrono::high_resolution_clock::now();
      std::initializer_list<IndexList> external_indices = {{}};
      if(i == 1)
        external_indices = {{L"i_1", L"a_1"}};
      else if(i == 2)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
      else if(i == 3)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
      else if(i == 4)
        external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}, {L"i_4", L"a_4"}};

      cc_st_r[i] = spintrace(cc_r[i], external_indices);
      // cc_st_r[i] = closed_shell_spintrace(cc_r[i], external_indices);
      canonicalize(cc_st_r[i]);
      if(i<3){
        // std::wcout <<  __LINE__ << to_latex_align(cc_st_r[i], 20, 5) << "\n";
        std::wcout <<  __LINE__ << "\n" << to_latex(cc_st_r[i]) << "\n";
      }
      printf("R%d Spin-orbit: %lu terms;\nSPINTRACED: With S operator: %lu;", i, cc_r[i]->size(), cc_st_r[i]->size());
      if(expand_S){
        cc_st_r[i] = expand_S_operator(cc_st_r[i]);
        rapid_simplify(cc_st_r[i]);
        canonicalize(cc_st_r[i]);
        printf("S expanded: %lu\n", cc_st_r[i]->size());
        // std::wcout << __LINE__ << to_latex(cc_st_r[i]) << "\n" << std::endl;

        if(bt){

          // Checks if the replacement map is a canonical sequence
          auto is_canonical = [&] (const std::map<Index, Index>& idx_map){
            bool canonical = true;
            for(auto&& pair: idx_map) if(pair.first != pair.second) return false;
            return canonical;
          };

          // Get coefficients and replacement maps
          auto btc = biorthogonal_tran_coeff(external_indices.size(), 1.e-12);

          // TODO: external_indices are getting modified before.
          if(i == 1)
            external_indices = {{L"i_1", L"a_1"}};
          else if(i == 2)
            external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
          else if(i == 3)
            external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
          else if(i == 4)
            external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}, {L"i_4", L"a_4"}};

          auto idx_map = biorthogonal_tran_idx_map(external_indices);
          assert(btc.size() == idx_map.size());

          // Append scale and append all transformations
          Sum bt_expr{};
          auto btc_ptr = btc.begin();
          for(auto&& map: idx_map){
            ExprPtr transformed_expr{};
            if(is_canonical(map))
              transformed_expr = ex<Constant>(*btc_ptr) * cc_st_r[i];
            else
              transformed_expr = ex<Constant>(*btc_ptr) * transform_expression(cc_st_r[i], map);
            btc_ptr++;
            bt_expr.append(transformed_expr);
          }
          ExprPtr bt_expr_p = std::make_shared<Sum>(bt_expr);
//          if(i == 3)
//           bt_expr_p = ex<Tensor>(L"S", WstrList{L"i_1", L"i_2", L"i_3"},
//                                  WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) * bt_expr_p;
//
          expand(bt_expr_p);
          canonicalize(bt_expr_p);
          rapid_simplify(bt_expr_p);
          cc_st_r[i] = bt_expr_p;
//          std::wcout << __LINE__ << "\n" << to_latex_align(cc_st_r[i], 20, 5) << "\n";
          printf("Biorthogonal transform: %lu\n", cc_st_r[i]->size());

          if(factorize_S && i > 1){
            if(i == 1)
              external_indices = {{L"i_1", L"a_1"}};
            else if(i == 2)
              external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
            else if(i == 3)
              external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
            else if(i == 4)
              external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}, {L"i_4", L"a_4"}};

            cc_st_r[i] = factorize_S_operator(cc_st_r[i], external_indices, true);
            printf("Factorize S: %lu\n", cc_st_r[i]->size());
          }
        } else {
          // TODO: Scale residual equations
        }
      }
      // canonicalize(cc_st_r[i]);
      // rapid_simplify(cc_st_r[i]);
      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;
      printf("CC R%d size: %lu time: %5.3f sec.\n\n", i, cc_st_r[i]->size(), time_elapsed.count());
      // std::wcout << "R" << i << ":\n" << to_latex_align(cc_st_r[i], 20, 5) << "\n";
    }

    return 0;
#if 0 // Check if S is factorized correctly (get back 490 terms of CCSDT R3)
    {
      auto result = std::make_shared<Sum>();
      int terms_with_S = 0;
      int counter = 0;
      std::cout << "\ncc_st_r[3] size: " << cc_st_r[3]->size() << "\n";
      for (auto&& summand : *cc_st_r[3]) {
        ++counter;
        auto term = summand->clone();
        if (has_tensor_label(term, L"S")) {
          std::wcout << counter << ": " << to_latex(term) << "\n";
          ++terms_with_S;

          ExprPtr S_expanded = expand_S_operator(term);
          int s1 = S_expanded->size();
          std::wcout << to_latex(S_expanded) << "\n";

          auto temp2 = factorize_S_operator(
              S_expanded,
              {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}},
              true);
          int s2 = temp2->as<Sum>().size();
          std::wcout << to_latex(temp2) << "\n\n";

          result->append(temp2);
        } else {
          result->append(term);
        }
//        break;
      }
      std::cout << "\nTerms_with_S: " << terms_with_S;
      std::cout << "\nResult->size(): " << result->size();
    }
    return 0;
#endif

#if 0 // Check terms from CCSDT R3, where S is missing
    {
      int i_term = 0;
      for(auto&& summand: *cc_st_r[3]){
        ++i_term;
        if(summand->is<Product>()){
          if(summand->as<Product>().factor(0)->as<Tensor>().label() != L"S")
            std::wcout << i_term << ": " << to_latex(summand) << "\n";
        }
      }
    }
    return 0;
#endif

#if 0
    {
      auto idx_map = biorthogonal_tran_idx_map({{L"i_1", L"a_1"}, {L"i_2",L"a_2"}, {L"i_3", L"a_3"}});

      auto is_canonical = [&] (const std::map<Index, Index>& idx_map){
        bool canonical = true;
        for(auto&& pair: idx_map) if(pair.first != pair.second) return false;
        return canonical;
      };

      // Get coefficients and replacement maps
      auto btc = biorthogonal_tran_coeff(3, 1.e-12);

      int i_term = 0;
      // Check spintrace of each term
      for (auto&& summand : *cc_r[3]) {
        ++i_term;
        auto term = summand->clone();
        std::wcout << "Term: " << to_latex(term) << "\n";

        auto st_term = closed_shell_spintrace(
            term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
        std::wcout << "Spin traced: " << st_term->size() << " "
                   << to_latex(st_term) << "\n";

        auto S_expanded =
            expand_S_operator(st_term);  // , {{L"i_1", L"a_1"}, {L"i_2",L"a_2"}, {L"i_3", L"a_3"}}, true);
        std::wcout << "S_expanded: " << S_expanded->size() << " "
                   << to_latex(S_expanded) << "\n";

//        Sum bt_expr{};
//        auto btc_ptr = btc.begin();
//        for(auto&& map: idx_map){
//          ExprPtr transformed_expr{};
//          if(is_canonical(map))
//            transformed_expr = ex<Constant>(*btc_ptr) * S_expanded;
//          else
//            transformed_expr = ex<Constant>(*btc_ptr) * transform_expression(S_expanded, map);
//          btc_ptr++;
//          bt_expr.append(transformed_expr);
//        }
//        ExprPtr bt_expr_p = std::make_shared<Sum>(bt_expr);
//        std::wcout << "bt_expr_p: " << bt_expr_p->size() << " "
//                   << to_latex(bt_expr_p) << "\n";

        std::cout << "\n";

        if(i_term > 5) break;
      }
    }
#endif
    // return 0;

#if !CCSDT_eval
    std::wcout << "CCSD R1: " << to_latex(cc_st_r[1]) << "\n";
    std::wcout << "CCSD R2: " << to_latex(cc_st_r[2]) << "\n";
#endif

//    cc_st_r[1] = ex<Constant>(0.5) * cc_st_r[1];
//    expand(cc_st_r[1]);
//    canonicalize(cc_st_r[1]);
//    rapid_simplify(cc_st_r[1]);
//
//    cc_st_r[2] = ex<Constant>(0.25) * cc_st_r[2];
//    expand(cc_st_r[2]);
//    canonicalize(cc_st_r[2]);
//    rapid_simplify(cc_st_r[2]);

#if CCSDT_eval
    // cc_st_r[3] = ex<Constant>(5.0) * cc_st_r[3];
    // expand(cc_st_r[3]);
    // rapid_simplify(cc_st_r[3]);
#endif

    std::vector<std::shared_ptr<TA::TArrayD>> data_tensors = {Fock_oo, Fock_ov, Fock_vv,
        G_oooo, G_ooov, G_oovo,  G_oovv,  G_ovov,  G_ovvo, G_vovo, G_voov, G_ovvv, G_vovv, G_vvvv,
        t_ov, t_oovv};

    std::vector<sequant::ExprPtr> seq_tensors = {
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1"}, {L"i_1"})), // f_oo
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1",}, {L"a_1"})), // f_ov
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
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"i_1"}, {L"a_2",L"a_3"})), //vovv
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"a_3",L"a_4"})), //vvvv

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1"}, {L"a_1"})), //ov
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",L"i_2"}, {L"a_1",L"a_2"})) //oovv
    };

#if CCSDT_eval
    data_tensors.push_back(t_ooovvv);
    seq_tensors.push_back(std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1", L"i_2", L"i_3"}, {L"a_1", L"a_2", L"a_3"}))); //ooovvv
#endif

    using evaluate::HashType;
    using evaluate::EvalTree;
    using ContextMapType =
    sequant::container::map<HashType, std::shared_ptr<TA::TArrayD>>;

    ContextMapType context_map;

    assert(data_tensors.size() == seq_tensors.size());
    for (auto i = 0; i < seq_tensors.size(); ++i) {
      auto hash_val = EvalTree(seq_tensors.at(i)).hash_value();
      context_map.insert(ContextMapType::value_type(hash_val, data_tensors.at(i)));
    }

    printf("R1 size: %lu\n",cc_st_r[1]->size());
    std::wcout << to_latex(cc_st_r[1]) << "\n";

    printf("R2 size: %lu\n",cc_st_r[2]->size());
    std::wcout << to_latex(cc_st_r[2]) << "\n";

    // printf("R3 size: %lu\n",cc_st_r[3]->size());
    // std::wcout << to_latex(cc_st_r[3]) << "\n";

    bool swap_braket_labels = true;
    auto r1_tree = EvalTree(cc_st_r[1], swap_braket_labels);
    auto r2_tree = EvalTree(cc_st_r[2], swap_braket_labels);
#if CCSDT_eval
    canonicalize(cc_st_r[3]);
    rapid_simplify(cc_st_r[3]);
    printf("R3 size: %lu\n",cc_st_r[3]->size());
    // std::wcout << "CCSDT R3:\n" << to_latex_align(cc_st_r[3], 20) << "\n";
/*
    {
      auto r3_no_S = std::make_shared<Sum>();
      auto r3_with_S = std::make_shared<Sum>();
      for(auto&& summand: *cc_st_r[3]){
        if(summand->as<Product>().factor(0)->as<Tensor>().label() != L"S")
          r3_no_S->append(summand);
        else{
          r3_with_S->append(summand);
        }
      }
      std::wcout << "CCSDT R3 terms that don't have S operator: Size: "
                 << r3_no_S->size() << "\n"
                 << to_latex_align(r3_no_S, 20, 5) << std::endl;

      std::cout << r3_no_S->size() << "\n";
      auto new_result = factorize_S_operator(r3_no_S,  {{L"i_1", L"a_1"}, {L"i_3", L"a_3"}}, true);
      std::cout << "ia, kc pairs: " << new_result->size() << "\n";

      cc_st_r[3] = new_result + r3_with_S;
      expand(cc_st_r[3]);
      canonicalize(cc_st_r[3]);
      rapid_simplify(cc_st_r[3]);
      printf("R3 size: %lu\n",cc_st_r[3]->size());
      std::wcout << "\n" << to_latex_align(cc_st_r[3], 20, 5) << "\n";

      auto temp = expand_S_operator(cc_st_r[3]);
      cc_st_r[3] = temp;
      expand(cc_st_r[3]);
      canonicalize(cc_st_r[3]);
      rapid_simplify(cc_st_r[3]);
      printf("R3 size: %lu\n",cc_st_r[3]->size());
    }
*/
    auto r3_tree = EvalTree(cc_st_r[3], swap_braket_labels);
#endif
    return 0;

    const auto cc_conv = conv * 1e2;
    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto ecc = 0.0;
    bool diis = false;

    int start_diis = 3;
    int diis_size = 3;
    std::deque<TA::TArrayD> e1_vec, e2_vec, e3_vec;
    std::deque<TA::TArrayD> t1_vec, t2_vec, t3_vec;
    TA::TArrayD t_ov_prev;
    TA::TArrayD t_oovv_prev;
    TA::TArrayD t_ooovvv_prev;

    cout << "Iter   norm(t_ov)    norm(t_oovv)     ΔE(CC)          E(CC)       time(s)" << endl;
    cout << "============================================================================" << endl;
    auto CC_start = std::chrono::high_resolution_clock::now();
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      if(diis && iter > 1){

        // Error vector
        TA::TArrayD e1;
        e1("i,a") = (*t_ov)("i,a") - t_ov_prev("i,a");
        TA::TArrayD diis_t1;
        diis_t1 = TiledArray::clone(*t_ov);

        TA::TArrayD e2;
        e2("i,j,a,b") = (*t_oovv)("i,j,a,b") - t_oovv_prev("i,j,a,b");
        TA::TArrayD diis_t2;
        diis_t2 = TiledArray::clone(*t_oovv);

        // Vector of error matrices
        e1_vec.push_back(e1);
        e2_vec.push_back(e2);

        t1_vec.push_back(diis_t1);
        t2_vec.push_back(diis_t2);

        int B_size = std::min(iter, diis_size) - 1;

        if(e1_vec.size() > B_size){
          e1_vec.pop_front();
          e2_vec.pop_front();
          t1_vec.pop_front();
          t2_vec.pop_front();
        }

#if CCSDT_eval
        TA::TArrayD e3;
        e3("i,j,k,a,b,c") = (*t_ooovvv)("i,j,k,a,b,c") - t_ooovvv_prev("i,j,k,a,b,c");
        TA::TArrayD diis_t3;
        diis_t3 = TiledArray::clone(*t_ooovvv);

        e3_vec.push_back(e3);
        t3_vec.push_back(diis_t3);

        if(e1_vec.size() > B_size){
          e3_vec.pop_front();
          t3_vec.pop_front();
        }
#endif

        if(iter > start_diis){
          // Fill 'B' matrix
          Matrix B_t1(B_size + 1, B_size + 1);
          B_t1.setConstant(-1.0);
          B_t1(B_size, B_size) = 0.0;
          Matrix B_t2 = B_t1;

          for(int i = 0; i < B_size; ++i){
            for(int j = 0; j < B_size; ++j){
              B_t1(i,j) = TA::dot_product(e1_vec[i], e1_vec[j]);
              B_t2(i,j) = TA::dot_product(e2_vec[i], e2_vec[j]);
            }
          }

          // Solve, get coefficients
          Eigen::VectorXd rhs(B_size + 1);
          rhs.setZero();
          rhs(B_size) = -1;

          Eigen::ColPivHouseholderQR<Matrix> solver_t1(B_t1);
          Eigen::VectorXd coeff_t1 = solver_t1.solve(rhs);

          Eigen::ColPivHouseholderQR<Matrix> solver_t2(B_t2);
          Eigen::VectorXd coeff_t2 = solver_t2.solve(rhs);

          // Update T amplitudes
          TA::TArrayD new_t1(world, tr_ov);
          new_t1.fill(0.0);
          TA::TArrayD new_t2(world, tr_oovv);
          new_t2.fill(0.0);
          for(int i = 0; i < B_size; ++i){
            new_t1("i,a") += t1_vec[i]("i,a") * coeff_t1[i];
            new_t2("i,j,a,b") += t2_vec[i]("i,j,a,b") * coeff_t2[i];
          }

          (*t_ov)("i,a") = new_t1("i,a");
          (*t_oovv)("i,j,a,b") = new_t2("i,j,a,b");

#if CCSDT_eval
          Matrix B_t3(B_size + 1, B_size + 1);
          B_t3.setConstant(-1.0);
          B_t3(B_size, B_size) = 0.0;

          for(int i = 0; i < B_size; ++i){
            for(int j = 0; j < B_size; ++j) {
              B_t3(i,j) = TA::dot_product(e3_vec[i], e3_vec[j]);
            }
          }

          Eigen::ColPivHouseholderQR<Matrix> solver_t3(B_t3);
          Eigen::VectorXd coeff_t3 = solver_t3.solve(rhs);

          TA::TArrayD new_t3(world, tr_ooovvv);
          new_t3.fill(0.0);
          for(int i = 0; i < B_size; ++i){
            new_t3("i,j,k,a,b,c") += t3_vec[i]("i,j,k,a,b,c") * coeff_t3[i];
          }
          (*t_ooovvv)("i,j,k,a,b,c") = new_t3("i,j,k,a,b,c");
#endif
        }
      }

      auto R1 = r1_tree.evaluate(context_map);
      auto R2 = r2_tree.evaluate(context_map);
#if CCSDT_eval
      auto R3 = r3_tree.evaluate(context_map);
#endif

      auto tile_R1       = R1.find({0,0}).get();
      auto tile_t_ov     = (*t_ov).find({0,0}).get();

      auto tile_R2       = R2.find({0,0,0,0}).get();
      auto tile_t_oovv   = (*t_oovv).find({0,0,0,0}).get();

#if CCSDT_eval
      auto tile_R3       = R3.find({0,0,0,0,0,0}).get();
      auto tile_t_ooovvv   = (*t_ooovvv).find({0,0,0,0,0,0}).get();
#endif

      if(diis){
        t_ov_prev = (*t_ov).clone();
        t_oovv_prev = (*t_oovv).clone();
        cout << t_ov_prev;
        cout << t_oovv_prev;
      }

#if CCSDT_eval
      if(diis){
        t_ooovvv_prev = (*t_ooovvv).clone();
        cout << t_ooovvv_prev;
      }
      // cout << (*t_ooovvv) << "\n";
#endif

      // save previous norm
      auto norm_last = std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));

      // Update T amplitudes
      // update t_ov
      for (auto i = 0; i < ndocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) {
          tile_t_ov(i,a) += tile_R1(i,a)/tile_D_ov(i,a); } }

      // update t_oovv
      for (auto i = 0; i < ndocc; ++i) {
        for (auto j = 0; j < ndocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
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
        std::chrono::duration_cast<std::chrono::seconds>(CC_stop - CC_start);

    cout << "\nOut of loop after "
         << iter << " iterations.\n"
         << "Time: "
         << time_elapsed.count() << " sec."
         << endl;
    cout << "Total energy = " << enuc + ehf + ecc << " a.u." << endl;

    // Check water molecule, sto-3g
#if CCSDT_eval
      double ccsdt_correlation = -0.070813170670;
      assert(fabs(ccsdt_correlation - ecc) < 1e-10);
#else
      double ccsd_correlation = -0.070680451962;
      assert(fabs(ccsd_correlation - ecc) < 1e-10);
      /* CCSD ref value obtained with mpqc using this input
  {
      "units": "2010CODATA",
      "atoms": {
          "file_name": "h2o.xyz",
          "sort_input": "true",
          "charge": "0",
          "n_cluster": "1",
          "reblock": "4"
      },
      "obs": {
          "name": "STO-3G",
          "atoms": "$:atoms"
      },
      "wfn_world": {
          "atoms": "$:atoms",
          "basis": "$:obs",
          "screen": "schwarz"
      },
      "scf": {
          "type": "RHF",
          "wfn_world": "$:wfn_world"
      },
      "wfn": {
          "type": "CCSD",
          "wfn_world": "$:wfn_world",
          "export_orbitals": "true",
          "atoms": "$:atoms",
          "ref": "$:scf",
          "reduced_abcd_memory": "true",
          "frozen_core": "false",
          "occ_block_size": "2",
          "unocc_block_size": "2"
      },
      "property": {
          "type": "Energy",
          "precision": "1e-10",
          "wfn": "$:wfn",
          "value": {
              "value": "-75.012759831161077"
          }
      }
  }
       */
#endif
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


/// @brief Find coefficients for biorthogonal transformation
/// @detailed Given the number of external indices, this function calculates the permutation matrix,
/// counts it's eigenvalues to find normalization constant and calculates a pseudoinverse
/// to find the coefficients corresponding to the n! permutations.
/// For details, see: http://arxiv.org/abs/1805.00565
/// @param n_particles Number of external index group
/// @param threshold Cut-off for counting number of non-zero eigen values
/// @return A vector<double> of biorthogonal transformation coefficients
container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold){
  using namespace Eigen;

  int n = std::tgamma(n_particles + 1); // <- Dimension of permutation matrix is n_particles!
  // Permutation matrix
  Eigen::MatrixXd M(n,n);
  {
    M.setZero();
    size_t n_row = 0;
    container::svector<int, 6> v(n_particles), v1(n_particles);
    std::iota(v.begin(), v.end(), 0);
    std::iota(v1.begin(), v1.end(), 0);
    do {
      container::vector<double> permutation_vector;
      do {
        auto cycles = count_cycles(v1, v);
        permutation_vector.push_back(std::pow(-2, cycles));
      } while (std::next_permutation(v.begin(), v.end()));
      Eigen::VectorXd pv_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(permutation_vector.data(),
                                                                             permutation_vector.size());
      M.row(n_row) = pv_eig;
      ++n_row;
    } while (std::next_permutation(v1.begin(), v1.end()));
    M *= std::pow(-1, n_particles);
    // std::cout << "permutation_matrix:\n" << M << "\n";
  }

  // Normalization constant
  double scalar;
  {
    // inline bool nonZero(double d) { return abs(d) > threshold ? true : false; }
    auto nonZero = [&] (const double d) { return abs(d) > threshold ? true : false; };

    // Solve system of equations
    SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
    container::vector<double> eig_vals(eig_solver.eigenvalues().size());
    VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
        eig_solver.eigenvalues();

    double non0count = std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
    scalar = eig_vals.size() / non0count;
  }

  // Find Pseudo Inverse, get 1st row only
  MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
  container::vector<double> result(pinv.rows());
  VectorXd::Map(&result[0], result.size()) = pinv.row(0) * scalar;
  return result;
}

/// @brief Biorthogonal transformation map
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const std::initializer_list<IndexList> ext_index_groups = {{}}){

  //Check size of external index group; catch exception otherwise
  if(ext_index_groups.size() == 0) throw( "Cannot compute index map since " && "ext_index_groups.size() == 0");
  assert(ext_index_groups.size() > 0);

  container::vector<Index> idx_list;
  for(auto&& idx_group : ext_index_groups) idx_list.push_back(*idx_group.begin());

  const container::vector<Index> const_idx_list = idx_list;
  // Do permutations and append to map
  std::vector<std::map<Index, Index>> result;
  do{
    std::map<Index, Index> map;
    auto const_list_ptr = const_idx_list.begin();
    for(auto&& i : idx_list){
      map.emplace(std::make_pair(*const_list_ptr, i));
      const_list_ptr++;
    }
    result.push_back(map);
  } while(std::next_permutation(idx_list.begin(), idx_list.end()));

  return result;
}

