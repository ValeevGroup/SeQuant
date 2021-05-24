#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/eval/eval_btas.hpp>
#include <SeQuant/domain/eval/read_tensor_btas.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/optimize/optimize.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>

#include <btas/btas.h>
#include <btas/tensor_func.h>
#include <fstream>
#include <range/v3/all.hpp>

#include <chrono>
#include <iomanip>
#include <iostream>

template <typename Tensor_t>
class yield_leaf {
 private:
  Tensor_t const &fock, &eri, &t_vo, &t_vvoo;

  size_t const nocc, nvirt;

  auto range1_limits(sequant::Tensor const& tensor) {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             auto sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? nocc : nvirt;
           });
  }

 public:
  yield_leaf(size_t no, size_t nv, Tensor_t const& F, Tensor_t const& G,
             Tensor_t const& ampl_vo, Tensor_t const& ampl_vvoo)
      : nocc{no},
        nvirt{nv},
        fock{F},
        eri{G},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo} {}

  Tensor_t operator()(sequant::Tensor const& texpr) {
    auto rank = texpr.bra_rank() + texpr.ket_rank();

    if (texpr.label() == L"t") {
      assert(rank == 2 || rank == 4 && "only t_vo and t_vvoo supported");
      return rank == 2 ? t_vo : t_vvoo;
    }

    assert((texpr.label() == L"g" || texpr.label() == L"f") &&
           "unsupported tensor label encountered");

    auto&& big_tensor = texpr.label() == L"g" ? eri : fock;

    auto r1_limits = range1_limits(texpr);
    auto iter_limits = r1_limits | ranges::views::transform([this](auto x) {
                         return x == nocc ? std::pair{size_t{0}, nocc}
                                          : std::pair{nocc, nocc + nvirt};
                       });

    auto slice = Tensor_t{btas::Range{r1_limits | ranges::to_vector}};

    if (iter_limits.size() == 2) {
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          slice(i - loop1.first, j - loop2.first) = big_tensor(i, j);

    } else {  // iter_limits.size() == 4 true
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      auto loop3 = iter_limits[2];
      auto loop4 = iter_limits[3];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          for (auto k = loop3.first; k < loop3.second; ++k)
            for (auto l = loop4.first; l < loop4.second; ++l)
              slice(i - loop1.first,    //
                    j - loop2.first,    //
                    k - loop3.first,    //
                    l - loop4.first) =  //
                  big_tensor(i, j, k, l);
    }

    return slice;
  }
};  // yield_leaf

auto const norm = [](btas::Tensor<double> const& tensor) {
  return std::sqrt(btas::dot(tensor, tensor));
};

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
  using sequant::eval::read_header;
  using sequant::eval::read_tensor_btas;
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

  auto const fock = read_tensor_btas(fock_ifname);
  auto const eri = read_tensor_btas(eri_ifname);

  size_t const nocc = fock_header.nocc;
  size_t const nvirt = fock_header.nvirt;

  auto t_vo = btas::Tensor<double>{btas::Range{nvirt, nocc}};
  auto t_vvoo = btas::Tensor<double>{btas::Range{nvirt, nvirt, nocc, nocc}};
  t_vo.fill(0);
  t_vvoo.fill(0);

  auto d_vo = btas::Tensor<double>{btas::Range{nvirt, nocc}};
  auto d_vvoo = btas::Tensor<double>{btas::Range{nvirt, nvirt, nocc, nocc}};
  d_vo.fill(0.);
  d_vvoo.fill(0.);
  for (auto a = 0; a < nvirt; ++a)
    for (auto i = 0; i < nocc; ++i) {
      d_vo(a, i) = fock(i, i) - fock(nocc + a, nocc + a);
      for (auto b = 0; b < nvirt; ++b)
        for (auto j = 0; j < nocc; ++j)
          d_vvoo(a, b, i, j) =
              d_vo(a, i) + fock(j, j) - fock(nocc + b, nocc + b);
    }

  // ============= SeQuant ==================== //

  using sequant::eqs::cceqvec;
  using sequant::optimize::optimize;
  using sequant::optimize::tail_factor;
  using evxpr_node = sequant::utils::binary_node<sequant::utils::eval_expr>;

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::mbpt::set_default_convention();
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto cc_r = cceqvec{2, 2}(false, true, true, true, true);

  // canonicalize expressions while optimizing
  bool canon = true;
  auto nodes = ranges::views::tail(cc_r) |
               ranges::views::transform([canon](auto const& seqxpr) {
                 return optimize(tail_factor(seqxpr), canon);
               }) |
               ranges::to_vector;

  auto const& r1_node = nodes[0];
  auto const& r2_node = nodes[1];

  auto yielder = yield_leaf{nocc, nvirt, fock, eri, t_vo, t_vvoo};

  auto eval_inst_r1 = sequant::eval::eval_instance_btas{r1_node};
  auto eval_inst_r2 = sequant::eval::eval_instance_btas{r2_node};

  // true: leaf tensors (other than 't' tensors) will be cached
  // false: only intermediates will be cached
  auto manager =
      sequant::eval::make_cache_man<btas::Tensor<double>>(nodes, true);

  auto const g_vvoo =
      yielder(sequant::utils::parse_expr(L"g_{a1,a2}^{i1,i2}",
                                         sequant::Symmetry::antisymm)
                  ->as<sequant::Tensor>());
  auto const f_vo = yielder(
      sequant::utils::parse_expr(L"f_{a1}^{i1}", sequant::Symmetry::antisymm)
          ->as<sequant::Tensor>());

  const auto maxiter = 100;
  const auto conv = 1e-12;

  size_t iter = 0;
  auto ediff = 0.0;
  auto normdiff = 0.0;
  auto ecc = 0.0;

  auto start = std::chrono::high_resolution_clock::now();
  do {
    ++iter;
    manager.reset_decaying();
    auto r1 = eval_inst_r1.evaluate_asymm(yielder, manager);
    auto r2 = eval_inst_r2.evaluate_asymm(yielder, manager);

    auto norm_last = std::sqrt(btas::dot(t_vvoo, t_vvoo));

    // updating T1 and T2
    for (auto a = 0; a < nvirt; ++a)
      for (auto i = 0; i < nocc; ++i) {
        t_vo(a, i) += r1(a, i) / d_vo(a, i);
        for (auto b = 0; b < nvirt; ++b)
          for (auto j = 0; j < nocc; ++j)
            t_vvoo(a, b, i, j) += r2(a, b, i, j) / d_vvoo(a, b, i, j);
      }
    normdiff = norm_last - std::sqrt(btas::dot(t_vvoo, t_vvoo));

    // calculating ecc

    auto ecc_last = ecc;
    decltype(r2) temp;

    btas::contract(1.0, g_vvoo, {'a', 'b', 'i', 'j'}, t_vo, {'a', 'i'}, 0.0,
                   temp, {'b', 'j'});
    ecc = 0.5 * btas::dot(temp, t_vo)         //
          + 0.25 * btas::dot(g_vvoo, t_vvoo)  //
          + btas::dot(f_vo, t_vo);

    ediff = ecc_last - ecc;

    cout << "E(CC) = "
         << std::setprecision(std::numeric_limits<double>::max_digits10) << ecc
         << endl;

    manager.reset_decaying();
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
