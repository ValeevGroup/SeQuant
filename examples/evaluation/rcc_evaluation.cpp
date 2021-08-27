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

// Use manually coded CC equations
#define MS_CC_EQ 0
#if MS_CC_EQ
#include "ms_cc_r.h"
#endif

// CCSDT or CCSD
#define CCSDT_eval 0

container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold);
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const container::vector<container::vector<Index>> ext_index_groups);
ExprPtr symmetrize_expr(ExprPtr& expr, const container::vector<container::vector<Index>> ext_index_groups);

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

    auto TA_fock = [&nao, &eps](const TA::Range& range) {
      TA::Tensor<double> tile(range);
      for (auto i = 0; i < nao; ++i) {
        for (auto j = 0; j < nao; ++j)
          tile(i, j) = (i == j) ? eps(i) : 0.0;
      }
      return tile;
    };

    std::vector<size_t> tile_boundaries = {0, nao};
    std::vector<TA::TiledRange1> ranges(
        2, TA::TiledRange1(tile_boundaries.begin(), tile_boundaries.end()));
    TA::TiledRange trange(ranges.begin(), ranges.end());

    TA::TArrayD fock(world, trange);
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
            emp2 += calc / (tile_fock(r,r) + tile_fock(s,s) - tile_fock(q,q) - tile_fock(p,p)); } } } }

    cout << "E(MP2): " << emp2 << "\nFinal energy is : " << hf_final + emp2 << " a.u." << endl;

    cout << "\n"
         << "***********************************\n"
         << "Coupled Cluster using TiledArray\n"
         << "***********************************\n"
         << endl;

    container::svector<size_t> r_occ = {0, ndocc};
    container::svector<size_t> r_vir = {0, nvirt};
    auto tr1o = TA::TiledRange1(r_occ.begin(), r_occ.end());
    auto tr1v = TA::TiledRange1(r_vir.begin(), r_vir.end());

    using TR1vec = std::vector<TA::TiledRange1>;
    TR1vec r_oo(2, tr1o);
    TR1vec r_vv(2, tr1v);
    TR1vec r_ov(1, tr1o); r_ov.push_back(tr1v);
    TR1vec r_vo(1, tr1v); r_vo.push_back(tr1o);
    TR1vec r_oooo(4, tr1o);
    TR1vec r_ooov(3, tr1o); r_ooov.push_back(tr1v);
    TR1vec r_oovo(3, tr1o); r_oovo.insert(r_oovo.end()-1, tr1v);
    TR1vec r_oovv(2, tr1o); r_oovv.insert(r_oovv.end(), 2, tr1v);
    TR1vec r_ovov(r_ov.begin(), r_ov.end()); r_ovov.insert(r_ovov.end(), r_ov.begin(), r_ov.end());
    TR1vec r_ovvo(r_ov.begin(), r_ov.end()); r_ovvo.insert(r_ovvo.end(), r_vo.begin(), r_vo.end());
    TR1vec r_vovo(r_vo.begin(), r_vo.end()); r_vovo.insert(r_vovo.end(), r_vo.begin(), r_vo.end());
    TR1vec r_voov(r_vo.begin(), r_vo.end()); r_voov.insert(r_voov.end(), r_ov.begin(), r_ov.end());
    TR1vec r_ovvv(3, tr1v); r_ovvv.insert(r_ovvv.begin(), tr1o);
    TR1vec r_vovv(3, tr1v); r_vovv.insert(r_vovv.begin()+1, tr1o);
    TR1vec r_vvvv(4, tr1v);

    using TTR = TA::TiledRange;
    TTR tr_oo(r_oo.begin(), r_oo.end());
    TTR tr_ov(r_ov.begin(), r_ov.end());
    TTR tr_vv(r_vv.begin(), r_vv.end());
    TTR tr_oooo(r_oooo.begin(), r_oooo.end());
    TTR tr_ooov(r_ooov.begin(), r_ooov.end());
    TTR tr_oovo(r_oovo.begin(), r_oovo.end());
    TTR tr_oovv(r_oovv.begin(), r_oovv.end());
    TTR tr_ovov(r_ovov.begin(), r_ovov.end());
    TTR tr_ovvo(r_ovvo.begin(), r_ovvo.end());
    TTR tr_voov(r_voov.begin(), r_voov.end());
    TTR tr_vovo(r_vovo.begin(), r_vovo.end());
    TTR tr_ovvv(r_ovvv.begin(), r_ovvv.end());
    TTR tr_vovv(r_vovv.begin(), r_vovv.end());
    TTR tr_vvvv(r_vvvv.begin(), r_vvvv.end());

#if CCSDT_eval
    TR1vec r_ooovvv(3,tr1o);
    r_ooovvv.insert(r_ooovvv.end(), 3, tr1v);
    TA::TiledRange tr_ooovvv(r_ooovvv.begin(), r_ooovvv.end());
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
    Logger::get_instance().wick_stats = false;

#if !MS_CC_EQ
    cout << "\n"
         << "***********************************\n";
#if CCSDT_eval // CCSDT
    cout << "            CCSDT\n";
    auto cc_r = cceqvec{3, 3}(true, true, true, true, true);
#else
    cout << "            CCSD\n";
    auto cc_r = cceqvec{2, 2}(true, true, true, true, true);
#endif
    cout << "***********************************\n" << endl;

    /// Make external index
    auto ext_idx_list = [] (const int i_max){
      container::vector<container::vector<Index>> ext_idx_list;

      for(size_t i = 1; i <= i_max; ++i) {
        auto label = std::to_wstring(i);
        auto occ_i = Index::make_label_index(IndexSpace::instance(IndexSpace::active_occupied), label);
        auto virt_i = Index::make_label_index(IndexSpace::instance(IndexSpace::active_unoccupied), label);
        container::vector<Index> pair = {occ_i, virt_i};
        ext_idx_list.push_back(pair);
      }
      return ext_idx_list;
    };

    // SPINspin. TRACE THE CC RESIDUAL EQUATIONS
    std::vector<ExprPtr> cc_st_r(cc_r.size());
    for (int i = 1; i < cc_r.size(); ++i){
      const auto tstart = std::chrono::high_resolution_clock::now();

      auto ext_idx_clone = ext_idx_list(i);
      cc_st_r[i] = closed_shell_spintrace(cc_r[i], ext_idx_clone);
      canonicalize(cc_st_r[i]);
      printf("R%d Spin-orbit: %lu terms;\nSPINTRACED: With S operator: %lu;", i, cc_r[i]->size(), cc_st_r[i]->size());

      for(auto&& term: *cc_st_r[i]){
        if(term->is<Product>()) term = remove_tensor_from_product(term->as<Product>(), L"S");
      }

      // Checks if the replacement map is a canonical sequence
      auto is_canonical = [] (const std::map<Index, Index>& idx_map){
        bool canonical = true;
        for(auto&& pair: idx_map) if(pair.first != pair.second) return false;
        return canonical;
      };

      // Get coefficients and replacement maps
      auto btc = biorthogonal_tran_coeff(ext_idx_clone.size(), 1.e-12);
      auto idx_map = biorthogonal_tran_idx_map(ext_idx_clone);
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
      cc_st_r[i] = std::make_shared<Sum>(bt_expr);

      if(i != 1)
        cc_st_r[i] = symmetrize_expr(cc_st_r[i], ext_idx_clone);
      simplify(cc_st_r[i]);

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;
      printf("CC R%d size: %lu time: %5.3f sec.\n\n", i, cc_st_r[i]->size(), time_elapsed.count());
    }

    auto ccsdt_r1 = cc_st_r[1];
    auto ccsdt_r2 = cc_st_r[2];
#if CCSDT_eval
    auto ccsdt_r3 = cc_st_r[3];
#endif
#else
#if !CCSDT_eval
    int n_cc = 2;
#else
    int n_cc = 3;
    auto ccsdt_r3 = r3(n_cc);
#endif
    auto ccsdt_r1 = r1(n_cc);
    auto ccsdt_r2 = r2(n_cc);
#endif
    simplify(ccsdt_r1);
    simplify(ccsdt_r2);
#if CCSDT_eval
    simplify(ccsdt_r3);
#endif

    std::vector<std::shared_ptr<TA::TArrayD>> data_tensors = {Fock_oo, Fock_ov, Fock_vv,
        G_oooo, G_ooov, G_oovo,  G_oovv,  G_ovov,  G_ovvo, G_vovo, G_ovvv, G_vovv, G_vvvv,
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

#if !CCSDT_eval

#if MS_CC_EQ
    auto ccsd_r1 = r1(2);
    auto ccsd_r2 = r2(2);
    simplify(ccsd_r1);
    simplify(ccsd_r2);
#else
    auto ccsd_r1 = cc_st_r[1];
    auto ccsd_r2 = cc_st_r[2];
#endif
    bool swap_braket_labels = true;
    auto r1_tree = EvalTree(ccsd_r1, swap_braket_labels);
    auto r2_tree = EvalTree(ccsd_r2, swap_braket_labels);
#endif

#if CCSDT_eval
#if MS_CC_EQ
    auto ccsdt_r1 = r1(3);
    auto ccsdt_r2 = r2(3);
    auto ccsdt_r3 = r3(3);
#endif
    bool swap_braket_labels = true;
    auto r1_tree = EvalTree(ccsdt_r1, swap_braket_labels);
    auto r2_tree = EvalTree(ccsdt_r2, swap_braket_labels);
    auto r3_tree = EvalTree(ccsdt_r3, swap_braket_labels);
#endif

    const auto cc_conv = conv * 1e2;
    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto ecc = 0.0;

    cout << "Iter   norm(t_ov)    norm(t_oovv)     ΔE(CC)          E(CC)       time(s)" << endl;
    cout << "============================================================================" << endl;
    auto CC_start = std::chrono::high_resolution_clock::now();
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

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

      auto norm_r1 = std::sqrt(R1("i,j").dot(R1("i,j")));
      auto norm_r2 = std::sqrt(R2("i,j,a,b").dot(R2("i,j,a,b")));

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

      // std::cout << "t1: " << tile_t_ov << "\nt2: " << tile_t_oovv << std::endl;
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
//#if CCSDT_eval
//      double ccsdt_correlation = -0.070813170670;
//      assert(fabs(ccsdt_correlation - ecc) < 1e-10);
//#else
//      double ccsd_correlation = -0.070680451962;
//      assert(fabs(ccsd_correlation - ecc) < 1e-10);
//      /* CCSD ref value obtained with mpqc using this input
//  {
//      "units": "2010CODATA",
//      "atoms": {
//          "file_name": "h2o.xyz",
//          "sort_input": "true",
//          "charge": "0",
//          "n_cluster": "1",
//          "reblock": "4"
//      },
//      "obs": {
//          "name": "STO-3G",
//          "atoms": "$:atoms"
//      },
//      "wfn_world": {
//          "atoms": "$:atoms",
//          "basis": "$:obs",
//          "screen": "schwarz"
//      },
//      "scf": {
//          "type": "RHF",
//          "wfn_world": "$:wfn_world"
//      },
//      "wfn": {
//          "type": "CCSD",
//          "wfn_world": "$:wfn_world",
//          "export_orbitals": "true",
//          "atoms": "$:atoms",
//          "ref": "$:scf",
//          "reduced_abcd_memory": "true",
//          "frozen_core": "false",
//          "occ_block_size": "2",
//          "unocc_block_size": "2"
//      },
//      "property": {
//          "type": "Energy",
//          "precision": "1e-10",
//          "wfn": "$:wfn",
//          "value": {
//              "value": "-75.012759831161077"
//          }
//      }
//  }
//       */
//#endif

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

// Generate S operator from external index list
ExprPtr symmetrize_expr(ExprPtr& expr, const container::vector<container::vector<Index>> ext_index_groups = {{}}){

  container::vector<Index> bra_list, ket_list;
  for(auto&& idx_group : ext_index_groups) {
    bra_list.push_back(*idx_group.begin());
    ket_list.push_back(*(idx_group.begin() + 1));
  }

  assert(bra_list.size() == ket_list.size());
  auto S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  return ex<Tensor>(S) * expr;
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
    // std::cout << "permutation_matrix:\n" << M << std::endl;
  }

  // Normalization constant
  double scalar;
  {
    // inline bool nonZero(double d) { return abs(d) > threshold ? true : false; }
    auto nonZero = [&threshold] (const double d) { return abs(d) > threshold ? true : false; };

    // Solve system of equations
    SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
    container::vector<double> eig_vals(eig_solver.eigenvalues().size());
    VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
        eig_solver.eigenvalues();

    double non0count = std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
    scalar = eig_vals.size() / non0count;
  }
  // std::cout << "scalar: " << scalar << std::endl;
  // Find Pseudo Inverse, get 1st row only
  MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
  container::vector<double> result(pinv.rows());
  VectorXd::Map(&result[0], result.size()) = pinv.row(0) * scalar;
  return result;
}

/// @brief Biorthogonal transformation map
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(
    const container::vector<container::vector<Index>> ext_index_groups = {{}}){
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

