#include <SeQuant/core/op.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <clocale>
#include <examples/eval/ta/eval_ta.hpp>
#include <examples/eval/ta/read_tensor_ta.hpp>
#include <iomanip>
#include <iostream>
#include <range/v3/view.hpp>

#include <btas/btas.h>
#include <btas/tensor_func.h>

namespace sequant::eval {

template<typename R, typename F>
void cartesian_foreach(const std::vector<R>& rs, F f) {
  using container::svector;
  using It = decltype(std::begin(rs[0]));
  using T = typename R::value_type;
  svector<It> its, ends;
  for (const auto& r : rs) {
    its.push_back(std::begin(r));
    ends.push_back(std::end(r));
  }
  while (its.front() != ends.front()) {
    svector<T> s;
    s.reserve(its.size());
    for (auto& it : its) {
      s.push_back(*it);
    }
    f(s);
    size_t i = its.size();
    while (i > 0) {
      --i;
      ++its[i];
      if (i == 0) break;
      if (its[i] != ends[i]) break;
      its[i] = std::begin(rs[i]);
    }
  }
}

} // namespace

template <typename Tensor_t>
class yield_leaf {
 private:
  size_t const nocc, nvirt;

  Tensor_t const &fock, &eri, &t_vo, &t_vvoo, &t_vvvooo;

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
             Tensor_t const& ampl_vo, Tensor_t const& ampl_vvoo,
             Tensor_t const& ampl_vvvooo)
      : nocc{no},
        nvirt{nv},
        fock{F},
        eri{G},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo},
        t_vvvooo{ampl_vvvooo} {}

  Tensor_t operator()(sequant::Tensor const& texpr) {
    auto const rank = texpr.bra_rank() + texpr.ket_rank();

    if (texpr.label() == L"t") {
      assert(rank == 2 || rank == 4 ||
             rank == 6 && "only t_vo, t_vvoo, t_vvvooo supported");
      return rank == 2 ? t_vo : rank == 4 ? t_vvoo : t_vvvooo;
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
  auto t_vvvooo =
      btas::Tensor<double>{btas::Range{nvirt, nvirt, nvirt, nocc, nocc, nocc}};
  t_vo.fill(0);
  t_vvoo.fill(0);
  t_vvvooo.fill(0);

  auto d_vo = btas::Tensor<double>{btas::Range{nvirt, nocc}};
  auto d_vvoo = btas::Tensor<double>{btas::Range{nvirt, nvirt, nocc, nocc}};
  auto d_vvvooo =
      btas::Tensor<double>{btas::Range{nvirt, nvirt, nvirt, nocc, nocc, nocc}};
  d_vo.fill(0);
  d_vvoo.fill(0);
  d_vvvooo.fill(0);
  for (auto a = 0; a < nvirt; ++a)
    for (auto i = 0; i < nocc; ++i) {
      d_vo(a, i) = fock(i, i) - fock(nocc + a, nocc + a);
      for (auto b = 0; b < nvirt; ++b)
        for (auto j = 0; j < nocc; ++j) {
          d_vvoo(a, b, i, j) =
              d_vo(a, i) + fock(j, j) - fock(nocc + b, nocc + b);
          for (auto c = 0; c < nvirt; ++c)
            for (auto k = 0; k < nocc; ++k)
              d_vvvooo(a, b, c, i, j, k) =
                  d_vvoo(a, b, i, j) + fock(k, k) - fock(nocc + c, nocc + c);
        }
    }

  // ============= SeQuant ==================== //

  using sequant::eqs::cceqvec;
  using sequant::optimize::optimize;
  using sequant::optimize::tail_factor;

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::mbpt::set_default_convention();
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto cc_r = cceqvec{3, 3}(true, true, true, true, true);

  // canonicalize expressions while optimizing
  bool canon = true;
  auto nodes = ranges::views::tail(cc_r) |
               ranges::views::transform([canon](auto const& seqxpr) {
                 // return optimize(tail_factor(seqxpr), canon);
                 return sequant::to_eval_node(tail_factor(seqxpr));
               }) |
               ranges::to_vector;

  auto const& r1_node = nodes[0];
  auto const& r2_node = nodes[1];
  auto const& r3_node = nodes[2];

  auto yielder = yield_leaf{nocc, nvirt, fock, eri, t_vo, t_vvoo, t_vvvooo};

  // true: leaf tensors (other than 't' tensors) will be cached
  // false: only intermediates will be cached
   auto manager =
      sequant::eval::make_cache_manager<btas::Tensor<double>>(nodes, true);

  auto const g_vvoo = yielder(
      sequant::parse_expr_asymm(L"g_{a1,a2}^{i1,i2}")->as<sequant::Tensor>());
  auto const f_vo =
      yielder(sequant::parse_expr_asymm(L"f_{a1}^{i1}")->as<sequant::Tensor>());

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
    auto r1 = sequant::eval::eval_antisymm(r1_node, yielder, manager);
    auto r2 = sequant::eval::eval_antisymm(r2_node, yielder, manager);
    auto r3 = sequant::eval::eval_antisymm(r3_node, yielder, manager);

    auto norm_last = std::sqrt(btas::dot(t_vvoo, t_vvoo));

    // updating T1 and T2
    for (auto a = 0; a < nvirt; ++a)
      for (auto i = 0; i < nocc; ++i) {
        t_vo(a, i) += r1(a, i) / d_vo(a, i);
        for (auto b = 0; b < nvirt; ++b)
          for (auto j = 0; j < nocc; ++j) {
            t_vvoo(a, b, i, j) += r2(a, b, i, j) / d_vvoo(a, b, i, j);
            for (auto c = 0; c < nvirt; ++c)
              for (auto k = 0; k < nocc; ++k)
                t_vvvooo(a, b, c, i, j, k) +=
                    r3(a, b, c, i, j, k) / d_vvvooo(a, b, c, i, j, k);
          }
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
