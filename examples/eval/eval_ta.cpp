#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/eval/eval_ta.hpp>
#include <SeQuant/domain/eval/read_tensor_ta.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/optimize/optimize.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>

#include <tiledarray.h>
#include <iomanip>
#include <iostream>
#include <range/v3/all.hpp>

template <typename Tensor_t>
struct yield_leaf {
  size_t const no, nv;
  Tensor_t const &G, &F, &t_vo, &t_vvoo;
  yield_leaf(size_t nocc, size_t nvirt, Tensor_t const& fock,
             Tensor_t const& eri, Tensor_t const& ampl_vo,
             Tensor_t const& ampl_vvoo)
      : no{nocc},
        nv{nvirt},
        G{eri},
        F{fock},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo}

  {}

  auto range1_limits(sequant::Tensor const& tensor) {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             auto sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? no : nv;
           });
  }

  Tensor_t operator()(sequant::Tensor const& tensor) {
    if (tensor.label() == L"t") {
      auto rank = tensor.rank();
      assert(rank == 1 || rank == 2);
      return rank == 1 ? t_vo : t_vvoo;
    }

    auto r1_limits = range1_limits(tensor);

    auto trange_vec = r1_limits | ranges::views::transform([](auto x) {
                        return TA::TiledRange1{0, x};
                      }) |
                      ranges::to_vector;

    auto iter_limits =
        r1_limits | ranges::views::transform([this](auto x) {
          return x == no ? std::pair{size_t{0}, no} : std::pair{no, no + nv};
        });

    auto tlabel = tensor.label();
    assert(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G : F;

    auto slice =
        TA::TArrayD{big_tensor.world(),
                    TA::TiledRange{trange_vec.begin(), trange_vec.end()}};
    slice.fill(0);
    auto tile_orig = big_tensor.find(0).get();
    auto tile_dest = slice.find(0).get();

    assert(iter_limits.size() == 2 || iter_limits.size() == 4);
    if (iter_limits.size() == 2) {
      for (auto ii = iter_limits[0].first; ii < iter_limits[0].second; ++ii)
        for (auto jj = iter_limits[1].first; jj < iter_limits[1].second; ++jj) {
          tile_dest(ii - iter_limits[0].first,  //
                    jj - iter_limits[1].first) = tile_orig(ii, jj);
        }
    } else {  // 4 iterations
      for (auto ii = iter_limits[0].first; ii < iter_limits[0].second; ++ii)
        for (auto jj = iter_limits[1].first; jj < iter_limits[1].second; ++jj)
          for (auto kk = iter_limits[2].first; kk < iter_limits[2].second; ++kk)
            for (auto ll = iter_limits[3].first; ll < iter_limits[3].second;
                 ++ll) {
              tile_dest(ii - iter_limits[0].first, jj - iter_limits[1].first,
                        kk - iter_limits[2].first, ll - iter_limits[3].first) =
                  tile_orig(ii, jj, kk, ll);
            }
    }

    return slice;
  }
};  // yield_leaf

// clang-format off
/**
 * <executable> (fock.dat eri.dat | eri.dat fock.dat)
 *
 * .dat format:
 *
 * size_t size_t size_t         # rank, nocc, nvirt
 * double                       # data ------
 * ...                          # data       |
 * ...                          # ....       |  no. of double entries = (nocc+nvirt)^rank
 * ...                          # data       |
 * double                       # data ------
 */
// clang-format on
int main(int argc, char** argv) {
  using sequant::eval::compatible_dims;
  using sequant::eval::make_trange;
  using sequant::eval::read_header;
  using sequant::eval::read_tensor_ta;

  using ranges::views::take;
  using sequant::optimize::optimize;
  using sequant::optimize::tail_factor;

  using std::cout;
  using std::endl;

  std::string_view fock_ifname = argc > 1 ? argv[1] : "fock.dat";
  std::string_view eri_ifname = argc > 2 ? argv[2] : "eri.dat";

  assert(compatible_dims(fock_ifname, eri_ifname));

  auto fock_header = read_header(fock_ifname);
  auto eri_header = read_header(eri_ifname);
  if (fock_header.rank > eri_header.rank) {
    std::swap(fock_ifname, eri_ifname);
    std::swap(fock_header, eri_header);
  }

  assert(fock_header.rank == 2 && "Fock tensor should be rank 2");
  assert(eri_header.rank == 4 && "Eri tensor should be rank 4");

  auto const nocc = fock_header.nocc;
  auto const nvirt = fock_header.nvirt;

  auto& world = TA::initialize(argc, argv);

  auto fock = TA::TArrayD{world, make_trange(fock_header.rank, nocc + nvirt)};

  auto eri = TA::TArrayD{world, make_trange(eri_header.rank, nocc + nvirt)};

  read_tensor_ta(fock_ifname, fock);
  read_tensor_ta(eri_ifname, eri);

  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::mbpt::set_default_convention();

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto D_vo = TA::TArrayD{world, TA::TiledRange{TA::TiledRange1{0, nvirt},
                                                TA::TiledRange1{0, nocc}}};
  auto D_vvoo = TA::TArrayD{
      world,
      TA::TiledRange{TA::TiledRange1{0, nvirt}, TA::TiledRange1{0, nvirt},
                     TA::TiledRange1{0, nocc}, TA::TiledRange1{0, nocc}}};

  D_vo.fill(0);
  D_vvoo.fill(0);
  auto tile_D_vo = D_vo.find(0).get();
  auto tile_D_vvoo = D_vvoo.find(0).get();
  [nocc, nvirt](auto const& fock, auto& dvo_tile, auto& dvvoo_tile) {
    auto tile_fock = fock.find(0).get();
    for (auto a = 0; a < nvirt; ++a)
      for (auto i = 0; i < nocc; ++i) {
        dvo_tile(a, i) = tile_fock(i, i) - tile_fock(nocc + a, nocc + a);
        for (auto b = 0; b < nvirt; ++b)
          for (auto j = 0; j < nocc; ++j) {
            dvvoo_tile(a, b, i, j) = dvo_tile(a, i) + tile_fock(j, j) -
                                     tile_fock(nocc + b, nocc + b);
          }
      }
  }(fock, tile_D_vo, tile_D_vvoo);

  auto t_vo = TA::TArrayD{world, TA::TiledRange{TA::TiledRange1{0, nvirt},
                                                TA::TiledRange1{0, nocc}}};
  auto t_vvoo = TA::TArrayD{
      world,
      TA::TiledRange{TA::TiledRange1{0, nvirt}, TA::TiledRange1{0, nvirt},
                     TA::TiledRange1{0, nocc}, TA::TiledRange1{0, nocc}}};
  t_vo.fill(0);
  t_vvoo.fill(0);

  auto yielder = yield_leaf{nocc, nvirt, fock, eri, t_vo, t_vvoo};

  auto const g_vvoo =
      yielder(sequant::utils::parse_expr(L"g_{a1,a2}^{i1,i2}",
                                         sequant::Symmetry::antisymm)
                  ->as<sequant::Tensor>());
  auto const f_vo = yielder(
      sequant::utils::parse_expr(L"f_{a1}^{i1}", sequant::Symmetry::antisymm)
          ->as<sequant::Tensor>());

  auto cc_r = sequant::eqs::cceqvec{2, 2}(false, true, true, true, true);

  // canonicalize expressions while optimizing
  bool canon = false;
  auto nodes = ranges::views::tail(cc_r) |
               ranges::views::transform([canon](auto const& n) {
                 return optimize(tail_factor(n), canon);
               }) |
               ranges::to_vector;

  auto const& node_r1 = nodes[0];
  auto const& node_r2 = nodes[1];

  // true: leaf tensors (other than 't' tensors) will be cached
  // false: only intermediates will be cached
  auto manager = sequant::eval::make_cache_man<TA::TArrayD>(nodes, true);

  const auto maxiter = 100;
  const auto conv = 1e-12;
  size_t iter = 0;
  double ediff = 0.0;
  auto normdiff = 0.0;
  auto ecc = 0.0;

  auto eval_inst_r1 = sequant::eval::eval_instance_ta{node_r1};
  auto eval_inst_r2 = sequant::eval::eval_instance_ta{node_r2};
  auto start = std::chrono::high_resolution_clock::now();
  do {
    ++iter;
    manager.reset_decaying();

    auto r1 = eval_inst_r1.evaluate_asymm(yielder, manager);
    auto r2 = eval_inst_r2.evaluate_asymm(yielder, manager);

    auto tile_r1 = r1.find(0).get();
    auto tile_r2 = r2.find(0).get();
    auto tile_t_vo = t_vo.find(0).get();
    auto tile_t_vvoo = t_vvoo.find(0).get();

    auto norm_last = t_vvoo.find(0).get().norm();

    // updating T1 and T2
    for (auto i = 0; i < nocc; ++i)
      for (auto a = 0; a < nvirt; ++a) {
        tile_t_vo(a, i) += tile_r1(a, i) / tile_D_vo(a, i);
        for (auto j = 0; j < nocc; ++j)
          for (auto b = 0; b < nvirt; ++b) {
            tile_t_vvoo(a, b, i, j) +=
                tile_r2(a, b, i, j) / tile_D_vvoo(a, b, i, j);
          }
      }

    normdiff = norm_last - t_vvoo.find(0).get().norm();

    // calculating ecc
    auto ecc_last = ecc;
    ecc = 0;
    TA::TArrayD temp;
    temp("b,j") = 0.5 * g_vvoo("a,b,i,j") * t_vo("a,i");

    ecc += temp("b,j").dot(t_vo("b,j"));
    ecc += 0.25 * g_vvoo("a,b,i,j").dot(t_vvoo("a,b,i,j"));
    ecc += f_vo("a,i").dot(t_vo("a,i"));

    cout << "E(CC) = "
         << std::setprecision(std::numeric_limits<double>::max_digits10) << ecc
         << endl;

  } while (iter < maxiter &&
           (std::fabs(normdiff) > conv || std::fabs(ediff) > conv));
  manager.reset_all();

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  cout << "\nOut of loop after " << iter << " iterations.\n"
       << "\nTime: " << duration.count() << " microseconds." << endl;

  return 0;
}
