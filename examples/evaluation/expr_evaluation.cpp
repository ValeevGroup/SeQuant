#include "../../examples/contract/scf/hartree-fock.h"
#include "../sequant_setup.hpp"

#include <SeQuant/domain/evaluate/eval_tree.hpp>

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
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

    auto shells = make_sto3g_basis(atoms);
    size_t nao = 0;
    for (auto s = 0; s < shells.size(); ++s) nao += shells[s].size();
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

    Matrix D;
    Matrix C;
    Eigen::VectorXd C_v;    // eigenvalues
    if (use_hcore_guess) {  // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();
      cout << "\n\tInitial C Matrix:\n";
      cout << C << endl;

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    } else {  // SOAD as the guess density, assumes STO-nG basis
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
    Matrix fock_mat;  // capture fock matrix for ccsd calculations
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Save a copy of the energy and the density
      auto ehf_last = ehf;
      auto D_last = D;

      // build a new Fock matrix
      auto F = H;
      // F += compute_2body_fock_simple(shells, D);
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
        for (auto j = 0; j < nao; j++) ehf += D(i, j) * (H(i, j) + F(i, j));

      // compute difference with last iteration
      ediff = ehf - ehf_last;
      rmsd = (D - D_last).norm();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        cout << "\n\n Iter        E(elec)              E(tot)               "
                "Delta(E)             RMS(D)         Time(s)\n";
      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf,
             ehf + enuc, ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    libint2::finalize();  // done with libint

    auto hf_final = ehf + enuc;
    cout << endl;
    printf("** Hartree-Fock energy = %20.12f\n", hf_final);

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
    nao *= 2;
    nvirt *= 2;
    size_t nocc = 2 * ndocc;  // in place of ndocc use nocc

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
              size_t p_i = floor(p / 2);
              size_t q_i = floor(q / 2);
              size_t r_i = floor(r / 2);
              size_t s_i = floor(s / 2);

              auto col_int = 0.0;
              auto exc_int = 0.0;

              // orbs: 0 1 2 3 4 ...
              // spin: α β α β α ...
              //
              //  <01|23>: coulomb integral

              if ((r % 2 == p % 2) && (s % 2 == q % 2))
                col_int += tile_ints_spatial(r_i, s_i, p_i, q_i);

              if ((r % 2 == q % 2) && (s % 2 == p % 2))
                exc_int += tile_ints_spatial(r_i, s_i, q_i, p_i);

              tile(r, s, p, q) = col_int - exc_int;
            }
          }
        }
      }
      return tile;
    };

    std::vector<size_t> tile_boundaries = {0, nao};
    std::vector<TA::TiledRange1> ranges(
        4, TA::TiledRange1(tile_boundaries.begin(), tile_boundaries.end()));
    TA::TiledRange trange(ranges.begin(), ranges.end());

    TA::TArrayD mo_ints_tensor_spin(world, trange);

    auto tile_spin_ints = mo_ints_tensor_spin.world().taskq.add(
        spinify_ints, mo_ints_tensor_spin.trange().make_tile_range(0));
    *(mo_ints_tensor_spin.begin()) = tile_spin_ints;

    auto spinify_fock = [&](const TA::Range& range) {
      TA::Tensor<double> tile(range);
      for (auto i = 0; i < nao; ++i) {
        for (auto j = 0; j < nao; ++j)
          tile(i, j) = (i == j) ? C_v(floor(i / 2)) : 0;
      }
      return tile;
    };

    std::vector<TA::TiledRange1> ranges2(
        2, TA::TiledRange1(tile_boundaries.begin(), tile_boundaries.end()));
    TA::TiledRange trange2(ranges2.begin(), ranges2.end());

    TA::TArrayD fock_spin(world, trange2);
    auto tile_fock = fock_spin.world().taskq.add(
        spinify_fock, fock_spin.trange().make_tile_range(0));
    *(fock_spin.begin()) = tile_fock;

    // See Step 4 at
    // http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4
    //

    cout << "\n"
         << "*********************************\n"
         << "Calculating EMP2 using TiledArray\n"
         << "*********************************\n"
         << endl;

    auto emp2 = 0.0;

    auto tile_ints_spin = mo_ints_tensor_spin.find({0, 0, 0, 0}).get();
    auto tile_fock_spin = fock_spin.find({0, 0}).get();

    for (auto r = 0; r < nocc; ++r) {
      for (auto s = 0; s < nocc; ++s) {
        for (auto p = nocc; p < nao; ++p) {
          for (auto q = nocc; q < nao; ++q) {
            auto calc = tile_ints_spin(r, s, p, q);

            calc *= calc;

            emp2 += calc / (tile_fock_spin(r, r) + tile_fock_spin(s, s) -
                            tile_fock_spin(p, p) - tile_fock_spin(q, q));
          }
        }
      }
    }
    emp2 /= 4.0;

    cout << "E(MP2): " << emp2 << "\nFinal energy is : " << hf_final + emp2
         << endl;
    /******************************************************/

    cout << "\n"
         << "***********************************\n"
         << "Calculating ECCSD using TiledArray\n"
         << "***********************************\n"
         << endl;

    container::svector<size_t> r_occ = {0, nocc};
    container::svector<size_t> r_vir = {0, nvirt};
    auto tr1o = TA::TiledRange1(r_occ.begin(), r_occ.end());
    auto tr1v = TA::TiledRange1(r_vir.begin(), r_vir.end());

    using TR1vec = std::vector<TA::TiledRange1>;
    TR1vec r_oo(2, tr1o);
    TR1vec r_vv(2, tr1v);
    TR1vec r_ov(1, tr1o); r_ov.push_back(tr1v);
    TR1vec r_vo(1,tr1v); r_vo.push_back(tr1o);
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

    TR1vec r_ooovvv(3,tr1o);
    r_ooovvv.insert(r_ooovvv.end(), 3, tr1v);
    TA::TiledRange tr_ooovvv(r_ooovvv.begin(), r_ooovvv.end());

    auto D_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    auto D_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    // auto D_ooovvv = std::make_shared<TA::TArrayD>(world, tr_ooovvv);
    auto Fock_oo = std::make_shared<TA::TArrayD>(world, tr_oo);
    auto Fock_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    auto Fock_vv = std::make_shared<TA::TArrayD>(world, tr_vv);
    auto G_oooo = std::make_shared<TA::TArrayD>(world, tr_oooo);
    auto G_ooov = std::make_shared<TA::TArrayD>(world, tr_ooov);
    auto G_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
    auto G_ovov = std::make_shared<TA::TArrayD>(world, tr_ovov);
    auto G_ovvv = std::make_shared<TA::TArrayD>(world, tr_ovvv);
    auto G_vvvv = std::make_shared<TA::TArrayD>(world, tr_vvvv);

    (*D_ov).fill(0.0);
    (*D_oovv).fill(0.0);
    // (*D_ooovvv).fill(0.0);
    (*Fock_oo).fill(0.0);
    (*Fock_ov).fill(0.0);
    (*Fock_vv).fill(0.0);
    (*G_oooo).fill(0.0);
    (*G_ooov).fill(0.0);
    (*G_oovv).fill(0.0);
    (*G_ovov).fill(0.0);
    (*G_ovvv).fill(0.0);
    (*G_vvvv).fill(0.0);

    auto tile_D_ov = (*D_ov).find({0, 0}).get();
    auto tile_D_oovv = (*D_oovv).find({0, 0, 0, 0}).get();
    // auto tile_D_ooovvv = (*D_ooovvv).find({0, 0, 0, 0, 0, 0}).get();
    auto tile_Fock_oo = (*Fock_oo).find({0, 0}).get();
    auto tile_Fock_ov = (*Fock_ov).find({0, 0}).get();
    auto tile_Fock_vv = (*Fock_vv).find({0, 0}).get();
    auto tile_G_oooo = (*G_oooo).find({0, 0, 0, 0}).get();
    auto tile_G_ooov = (*G_ooov).find({0, 0, 0, 0}).get();
    auto tile_G_oovv = (*G_oovv).find({0, 0, 0, 0}).get();
    auto tile_G_ovov = (*G_ovov).find({0, 0, 0, 0}).get();
    auto tile_G_ovvv = (*G_ovvv).find({0, 0, 0, 0}).get();
    auto tile_G_vvvv = (*G_vvvv).find({0, 0, 0, 0}).get();

    for (auto i = 0; i < nocc; ++i) {
      tile_Fock_oo(i, i) = tile_fock_spin(i, i);
      for (auto a = 0; a < nvirt; ++a) {
        tile_Fock_ov(i, a) = tile_fock_spin(i, a + nocc);
        tile_D_ov(i, a) =
            tile_fock_spin(i, i) - tile_fock_spin(a + nocc, a + nocc);
        for (auto j = 0; j < nocc; ++j) {
          for (auto b = 0; b < nvirt; ++b) {
            tile_Fock_vv(a, b) = tile_fock_spin(a + nocc, b + nocc);
            tile_D_oovv(i, j, a, b) = tile_D_ov(i, a) + tile_fock_spin(j, j) -
                                      tile_fock_spin(b + nocc, b + nocc);
            tile_G_oovv(i, j, a, b) = tile_ints_spin(i, j, a + nocc, b + nocc);
            tile_G_ovov(i, a, j, b) = tile_ints_spin(i, a + nocc, j, b + nocc);
          }
        }
      }
    }

    for (auto i = 0; i < nocc; ++i)
      for (auto a = 0; a < nvirt; ++a)
        for (auto b = 0; b < nvirt; ++b)
          for (auto c = 0; c < nvirt; ++c)
            tile_G_ovvv(i, a, b, c) =
                tile_ints_spin(i, a + nocc, b + nocc, c + nocc);

    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto a = 0; a < nvirt; ++a)
            tile_G_ooov(i, j, k, a) = tile_ints_spin(i, j, k, a + nocc);

    for (auto i = 0; i < nocc; ++i)
      for (auto j = 0; j < nocc; ++j)
        for (auto k = 0; k < nocc; ++k)
          for (auto l = 0; l < nocc; ++l)
            tile_G_oooo(i, j, k, l) = tile_ints_spin(i, j, k, l);

    for (auto a = 0; a < nvirt; ++a)
      for (auto b = 0; b < nvirt; ++b)
        for (auto c = 0; c < nvirt; ++c)
          for (auto d = 0; d < nvirt; ++d)
            tile_G_vvvv(a, b, c, d) =
                tile_ints_spin(a + nocc, b + nocc, c + nocc, d + nocc);

//    for (auto i = 0; i < nocc; ++i)
//      for (auto j = 0; j < nocc; ++j)
//        for (auto k = 0; k < nocc; ++k)
//          for (auto a = 0; a < nvirt; ++a)
//            for (auto b = 0; b < nvirt; ++b)
//              for (auto c = 0; c < nvirt; ++c)
//                tile_D_ooovvv(i, j, k, a, b, c) =
//                    tile_D_oovv(i, j, a, b) + tile_D_ov(k, c);

    // amplitudes for coupled-cluster calculations
    auto t_ov = std::make_shared<TA::TArrayD>(world, tr_ov);
    auto t_oovv = std::make_shared<TA::TArrayD>(world, tr_oovv);
//    auto t_ooovvv = std::make_shared<TA::TArrayD>(world, tr_ooovvv);
    (*t_ov).fill(0.);
    (*t_oovv).fill(0.);
//    (*t_ooovvv).fill(0.0);

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

    auto cc_r = cceqvec{2, 2}(true, true, true, true);

    std::vector<std::shared_ptr<TA::TArrayD>> data_tensors = {Fock_oo, Fock_ov, Fock_vv, G_oooo,
                                        G_ooov,  G_oovv,  G_ovov,  G_ovvv,
                                        G_vvvv,  t_ov,    t_oovv};

    std::vector<sequant::ExprPtr> seq_tensors = {
        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1"}, {L"i_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"i_1"}, {L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"f", {L"a_1"}, {L"a_2"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"i_4"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"i_3",L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"i_2"}, {L"a_1",L"a_2"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"i_2",L"a_2"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"i_1",L"a_1"}, {L"a_2",L"a_3"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"g", {L"a_1",L"a_2"}, {L"a_3",L"a_4"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",}, {L"a_1"})),

        std::make_shared<sequant::Tensor>(sequant::Tensor(L"t", {L"i_1",L"i_2"}, {L"a_1",L"a_2"}))
    };

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

    auto r1_tree = EvalTree(cc_r[1]);
    auto r2_tree = EvalTree(cc_r[2]);

    iter = 0;
    rmsd = 0.0;
    ediff = 0.0;
    auto normdiff = 0.0;
    auto ecc = 0.0;
    auto tstart = std::chrono::high_resolution_clock::now();
    do {
        ++iter;
        cout << "using TiledArray,    iter: " << iter << endl;
        auto R1 = r1_tree.evaluate(context_map);
        auto R2 = r2_tree.evaluate(context_map);

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

      auto ecc_last = ecc;

      // calculating energy
      TA::TArrayD temp_tensor;
      temp_tensor("j,b") = (*G_oovv)("i,j,a,b")*(*t_ov)("i,a");

      ecc = 0.5*temp_tensor("i,a").dot((*t_ov)("i,a"))
              + 0.25*(*G_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b"))
              + (*Fock_ov)("i,a").dot((*t_ov)("i,a"));

      printf("E(CC) is: %20.12f\n\n", ecc);

      normdiff = norm_last - std::sqrt((*t_oovv)("i,j,a,b").dot((*t_oovv)("i,j,a,b")));
      ediff    = ecc_last - ecc;
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
