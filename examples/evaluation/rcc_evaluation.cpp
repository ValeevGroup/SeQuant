//
// Created by Nakul Teke on 3/18/20.
//
#include "../../examples/contract/scf/hartree-fock.h"
#include "../../src/SeQuant/domain/evaluate/eval_context.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor_builder.hpp"
#include "../sequant_setup.hpp"

#include <TiledArray/initialize.h>
#include <tiledarray.h>
#include <libint2.hpp>

#include <array>
#include <memory>

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
    cout << "Nuclear repulsion energy = " << enuc << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

    auto shells = make_sto3g_basis(atoms);
    size_t nao = 0;
    for (auto s = 0; s < shells.size(); ++s) nao += shells[s].size();
    size_t nvirt = nao - ndocc;
    std::cout << "nao  : " << nao << "\n";
    std::cout << "ndocc: " << ndocc << "\n";
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
    cout << "\n\n Iter        E(elec)              E(tot)               "
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
    //

    /* cout << "Coeff matrix is: " << endl; */
    /* cout << C << endl << endl; */

    // Initialize MADWorld
    TA::World& world = TA::initialize(argc, argv);

    auto mo_ints_tensor = compute_mo_ints(C, shells, world);
    auto tile_ints_spatial = mo_ints_tensor.find({0, 0, 0, 0}).get();

    // from here on number of occupied and virtual orbitals
    // is doubled as we are now working on spin basis
    // nao *= 2;
    // nvirt *= 2;
    size_t nocc = ndocc;  // in place of ndocc use nocc
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
//    Matrix Fock = eps.asDiagonal();
//    cout << "Fock:\n" << Fock << endl;

    // See Step 4 at
    // http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4
    //

    cout << "\n"
         << "*********************************\n"
         << "Calculating EMP2 using TiledArray\n"
         << "*********************************\n"
         << endl;

    auto emp2 = 0.0;

    auto tile_ints = mo_ints_tensor.find({0, 0, 0, 0}).get();
    auto tile_fock = fock.find({0, 0}).get();
    // auto nocc = ndocc;

    for (auto r = 0; r < nocc; ++r) {
      for (auto s = 0; s < nocc; ++s) {
        for (auto p = ndocc; p < nao; ++p) {
          for (auto q = ndocc; q < nao; ++q) {
            auto calc = tile_ints(r, s, p, q) *
                   ( 2.0 * tile_ints(r, s, p, q) - tile_ints(r, s, q, p));
            // emp2 += calc / (eps(r) + eps(s) - eps(q) - eps(p));
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

    TA::TiledRange tr_oo{{0, nocc}, {0, nocc}};
    TA::TiledRange tr_ov{{0, nocc}, {0, nvirt}};
    TA::TiledRange tr_vv{{0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_oooo{{0, nocc}, {0, nocc}, {0, nocc}, {0, nocc}};
    TA::TiledRange tr_ooov{{0, nocc}, {0, nocc}, {0, nocc}, {0, nvirt}};
    TA::TiledRange tr_oovv{{0, nocc}, {0, nocc}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_ovov{{0, nocc}, {0, nvirt}, {0, nocc}, {0, nvirt}};
    TA::TiledRange tr_ovvv{{0, nocc}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    TA::TiledRange tr_vvvv{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, nvirt}};

    TA::TiledRange tr_vovv{{0, nvirt}, {0, nocc}, {0, nvirt}, {0, nvirt}};

#define CCSDT_eval 0
#if CCSDT_eval 
    // TA::TiledRange tr_ooovvv{{0, nocc},  {0, nocc},  {0, nocc},
    //                          {0, nvirt}, {0, nvirt}, {0, nvirt}};
    //
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
    auto G_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    auto G_ovov = std::make_shared<TA::TArrayD>(world, tr_ovov);
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
    (*G_oovv).fill(0.0);
    (*G_ovov).fill(0.0);
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
    auto tile_G_oovv = (*G_oovv).find({0, 0, 0, 0}).get();
    auto tile_G_ovov = (*G_ovov).find({0, 0, 0, 0}).get();
    auto tile_G_ovvv = (*G_ovvv).find({0, 0, 0, 0}).get();
    auto tile_G_vvvv = (*G_vvvv).find({0, 0, 0, 0}).get();
    auto tile_G_vovv = (*G_vovv).find({0, 0, 0, 0}).get();

    for (auto i = 0; i < nocc; ++i) {
      tile_Fock_oo(i, i) = tile_fock(i, i);
      for (auto a = 0; a < nvirt; ++a) {
        tile_Fock_ov(i, a) = tile_fock(i, a + nocc);
        tile_D_ov(i, a) =
            tile_fock(i, i) - tile_fock(a + nocc, a + nocc);
        for (auto j = 0; j < nocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            tile_Fock_vv(a, b) = tile_fock(a + nocc, b + nocc);
            tile_D_oovv(i, j, a, b) = tile_D_ov(i, a) + tile_fock(j, j) -
                tile_fock(b + nocc, b + nocc);
            tile_G_oovv(i, j, a, b) = tile_ints_spatial(i, j, a + nocc, b + nocc);
            tile_G_ovov(i, a, j, b) = tile_ints_spatial(i, a + nocc, j, b + nocc);
          }
        }
      }
    }

    for (auto i = 0; i < nocc; ++i)
      for (auto a = 0; a < nvirt; ++a)
        for (auto b = 0; b < nvirt; ++b)
          for (auto c = 0; c < nvirt; ++c)
            tile_G_ovvv(i, a, b, c) =
                tile_ints_spatial(i, a + nocc, b + nocc, c + nocc);

    for (auto a = 0; a < nvirt; ++a)
      for (auto i = 0; i < nocc; ++i)
        for (auto b = 0; b < nvirt; ++b)
          for (auto c = 0; c < nvirt; ++c)
            tile_G_vovv(a, i, b, c) =
                tile_ints_spatial(a + nocc, i, b + nocc, c + nocc);


    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto a = 0; a < nvirt; ++a)
            tile_G_ooov(i, j, k, a) = tile_ints_spatial(i, j, k, a + nocc);

    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto l = 0; l < nocc; ++l)
            tile_G_oooo(i, j, k, l) = tile_ints_spatial(i, j, k, l);

    for (auto a = 0; a < nvirt; ++a)
      for (auto b = 0; b < nvirt; ++b)
        for (auto c = 0; c < nvirt; ++c)
          for (auto d = 0; d < nvirt; ++d)
            tile_G_vvvv(a, b, c, d) =
                tile_ints_spatial(a + nocc, b + nocc, c + nocc, d + nocc);
#if CCSDT_eval
    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
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

    auto ccsd_r = cceqvec{2, 2}(true, true, true, true);
    std::initializer_list<IndexList> external_indices = {{}};

    // SPIN TRACE THE RESIDUAL
    std::vector<ExprPtr> cc_r(ccsd_r.size());
    for (size_t i = 1; i < ccsd_r.size(); ++i){
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

#define SIMPLIFIED_R2 1
#if SIMPLIFIED_R2
    // 1/3 R + 1/6 R' for simpler equations
    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};

    auto temp_expr = transform_expression(cc_r[2], idxmap);
    auto simpler_R2 =
        ex<Constant>(1.0 / 3.0) * cc_r[2] + ex<Constant>(1.0 / 6.0) * temp_expr;

    expand(simpler_R2);
    rapid_simplify(simpler_R2);
    canonicalize(simpler_R2);
    rapid_simplify(simpler_R2);

    auto P = ex<Tensor>(L"P", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}, Symmetry::nonsymm);
    auto P_R2 = P * simpler_R2;

    expand(P_R2);
    rapid_simplify(P_R2);
    canonicalize(P_R2);
    rapid_simplify(P_R2);
#endif

    std::vector<std::shared_ptr<TA::TArrayD>> data_tensors = {Fock_oo, Fock_ov, Fock_vv, G_oooo,
                                                              G_ooov,  G_oovv,  G_ovov,  G_ovvv, G_vovv,
                                                              G_vvvv,  t_ov,    t_oovv};

    // TODO: DO I NEED TO GET PERMUTED MAPS ?
    std::vector<sequant::ExprPtr> seq_tensors = {
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1"}, {L"i_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1",}, {L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"a_1",}, {L"a_2",})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"i_4"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"a_1",L"a_2"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"i_2",L"a_2"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"a_2",L"a_3"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1", L"i_1"}, {L"a_2",L"a_3"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"a_3",L"a_4"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",}, {L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",L"i_2"}, {L"a_1",L"a_2"}))
    };

    container::map<ExprPtr, std::shared_ptr<TA::TArrayD>> context_builder;

    assert(data_tensors.size() == seq_tensors.size());
    for (auto i = 0; i < seq_tensors.size(); ++i) {
      context_builder.insert(decltype(context_builder)::value_type(seq_tensors.at(i),
                                                                   data_tensors.at(i)));
    }

    auto context = evaluate::EvalContext(context_builder);
    auto builder = sequant::evaluate::EvalTensorBuilder<TA::TArrayD>();

    std::wcout << "CCSD R1:\n" << to_latex(cc_r[1]) << endl;
    std::wcout << "CCSD R2:\n" << to_latex(cc_r[2]) << endl;

    auto r1_tree = builder.build_tree(cc_r[1]);
    auto r2_tree = builder.build_tree(cc_r[2]);

# if 1
    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto ecc = 0.0;
    cout << "using TiledArray" << endl;
    cout << "iter     norm(t_ov)     norm(t_oovv)    ΔE(CC)    E(CC)" << endl;
    cout << "=======================================================" << endl;
    auto tstart = std::chrono::high_resolution_clock::now();
    do {
      ++iter;
      auto R1 = r1_tree->evaluate(context.get_map());
      cout << "R1: " << R1 << endl;
      auto R2 = r2_tree->evaluate(context.get_map());
      cout << "R2: " << R2 << endl;

      auto tile_R1       = R1.find({0,0}).get();
      auto tile_t_ov     = (*t_ov).find({0,0}).get();

      auto tile_R2       = R2.find({0,0,0,0}).get();
      auto tile_t_oovv   = (*t_oovv).find({0,0,0,0}).get();

      // save previous norm
      auto norm_last = std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));

      //////////////
      // Updating amplitudes

      // update t_ov
      for (auto i = 0; i < nocc; ++i) {
        for (auto a = 0; a < nvirt; ++a) {
          tile_t_ov(i,a) += tile_R1(i,a)/tile_D_ov(i,a); } }

      //
      // update t_oovv
      for (auto i = 0; i < nocc; ++i) {
        for (auto j = 0; j < nocc; ++j) {
          for (auto a = 0; a < nvirt; ++a) {
            for (auto b = 0; b < nvirt; ++b) {
              tile_t_oovv(i,j,a,b) += tile_R2(i,j,a,b)/tile_D_oovv(i,j,a,b); } } } }

      cout << "norm(t_ov) "   << std::sqrt((*t_ov)("i,j").dot((*t_ov)("i,j")))             <<endl;
      cout << "norm(t_oovv) " << std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b"))) <<endl;
      cout << iter << "   " <<  std::sqrt((*t_ov)("i,j").dot((*t_ov)("i,j"))) << "     ";
      cout << std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b"))) << "     ";
      auto ecc_last = ecc;

      // calculating energy
      TA::TArrayD temp_tensor;
      temp_tensor("j,b") = (*G_oovv)("i,j,a,b")*(*t_ov)("i,a");

      // TODO: Spin-free CC expression here
//      ecc = 0.5*temp_tensor("i,a").dot((*t_ov)("i,a"))
//          + 0.25*(*G_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b"))
//          + (*Fock_ov)("i,a").dot((*t_ov)("i,a"));
      ecc = 2.0; // * (*Fock_ov)("i,a").dot((*t_ov)("i,a")) +
            

      normdiff = norm_last - std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));
      ediff    = ecc_last - ecc;
      cout << ediff << "    " << ecc << endl;
    } while((fabs(normdiff) > conv || fabs(ediff) > conv) && (iter < maxiter));
    //} while(false);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto time_elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    std::cout << "\nOut of loop after "
              << iter << " iterations.\n"
              << "\nTime: "
              << time_elapsed.count() << " microseconds"
              << std::endl;

    /* auto r1 = r1_tree->evaluate(context.get_map()); */
    /* std::cout << "norm(r1) = " << std::sqrt(r1("i,a").dot(r1("i,a"))) << std::endl; */

    /* auto r2 = r2_tree->evaluate(context.get_map()); */
    /* std::cout << "norm(r2) = " << std::sqrt(r2("i,j,a,b").dot(r2("i,j,a,b"))) << std::endl; */

    /* std::wcout << "Digraph for R1\n------------\n"; */
    /* std::wcout << r1_tree->to_digraph() << std::endl; */
    /* std::wcout << "Digraph for R2\n------------\n"; */
    /* std::wcout << r2_tree->to_digraph() << std::endl; */
#endif
    TA::finalize();

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

