//
// Created by Bimal Gaudel on 7/14/21.
//

#ifndef SEQUANT_EVAL_DATA_WORLD_TA_HPP
#define SEQUANT_EVAL_DATA_WORLD_TA_HPP

#include "examples/eval/data_info.hpp"
#include "examples/eval/eval_utils.hpp"

#include <tiledarray.h>
#include <range/v3/view.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant::eval {

template <typename Tensor_t>
class DataWorldTA {
 private:

  size_t const nocc;

  size_t const nvirt;

  std::string const fock_file;

  std::string const eri_file;

  Tensor_t G_pqrs;

  Tensor_t F_pq;

  container::vector<Tensor_t> Ts;

  ///
  /// Read tensor data from a file to a TA::DistArray tensor.
  ///
  /// \param fname Input file name
  /// \param tensor Output TA::DistArray tensor.
  template <typename T>
  static void read_tensor_ta(std::string_view fname, T& tensor) {
    // TODO assert tensor single tiled or extend to handle multiply tiled case
    auto ta_tensor = TA::Tensor<typename T::numeric_type>{tensor.trange().make_tile_range(0)};
    read_tensor(fname, ta_tensor);
    *tensor.begin() = ta_tensor;
  }

  ///
  /// Maps the IndexSpace type of an index in the braket of a tensor
  /// to nocc (for IndexSpace::active_occupied)
  /// or nvirt (for IndexSpace::active_unoccupied)
  ///
  /// \param tensor sequant::Tensor
  /// \return View of an iterable with size_t-type elements.
  ///
  auto range1_limits(sequant::Tensor const& tensor) const {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             const auto& sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? nocc : nvirt;
           });
  }

  ///
  /// From a container of pair-likes returns a TA::TiledRange
  ///
  /// \tparam Iterable type of the container
  /// \param pairs Iterable with std::tuple<size_t, size_t>-like elements.
  /// \return TA::TiledRange
  ///
  template <typename Iterable>
  auto make_trange(Iterable const& pairs) const {
    using ranges::views::transform;
    auto tr1s = pairs
                | transform([](auto&& pair){
                    auto const [lo, hi] = pair;
        return TA::TiledRange1{lo, hi};
    }) | ranges::to_vector;

    return TA::TiledRange{tr1s.begin(), tr1s.end()};
  }

  ///
  /// Update an amplitude tensor using a corresponding residual tensor.
  ///
  /// \param T Amplitude tensor with the layout {virt, virt, ..., occ, occ,}
  /// \param R Residual tensor with the layout {virt, virt, ..., occ, occ,}
  /// \note The ranks and the layouts of T and R must match exactly.
  /// \note The rank is even.
  ///
  void update_single_T(Tensor_t& T, Tensor_t const& R) {
    using namespace ranges::views;
    assert(T.trange().rank() == R.trange().rank()
           && "Incompatible ranks of R and T");

    assert(T.trange().rank() % 2 == 0 && "Odd rank not supported");

    size_t const n = T.trange().rank() / 2;

    auto tile_T = T.find(0).get();
    auto tile_R = R.find(0).get();
    auto tile_F = F_pq.find(0).get();
    auto updater = [&tile_T, &tile_R, &tile_F,
                    n, this](auto const& idx_vec){
      double diag = 0.;
      for (auto x: idx_vec | take(n))
        diag -= tile_F(x + nocc, x + nocc);
      for (auto x: idx_vec | drop(n) | take(n))
        diag += tile_F(x, x);

      tile_T(idx_vec) += tile_R(idx_vec) / diag;
    };

    auto virt_occ = concat(
        repeat_n(iota(size_t{0}, nvirt) | ranges::to_vector, n),
        repeat_n(iota(size_t{0}, nocc)  | ranges::to_vector, n))
                    | ranges::to_vector;

    cartesian_foreach(virt_occ, updater);
  }

 public:
  DataWorldTA(TA::World& world,
              size_t excit, DataInfo const& info):  nocc{info.nocc()},
                                       nvirt{info.nvirt()},
                                       fock_file{info.fock_file()},
                                       eri_file{info.eri_file()} {
    using namespace ranges::views;

    //--------------------------------------------------------
    // read ERI and Fock tensors
    //--------------------------------------------------------
    auto const nobs = nocc + nvirt;
    auto const fock_trange = make_trange(zip(repeat(0), repeat(nobs))
                                         | take(DataInfo::fock_rank));
    auto const eri_trange = make_trange(zip(repeat(0), repeat(nobs))
                                        | take(DataInfo::eri_rank));
    F_pq = Tensor_t{world, fock_trange};
    G_pqrs = Tensor_t{world, eri_trange};
    read_tensor_ta(fock_file, F_pq);
    read_tensor_ta(eri_file, G_pqrs);
    //--------

    //--------------------------------------------------------
    // Initialize amplitude tensors
    //--------------------------------------------------------
    Ts.reserve(excit);
    for (auto i = 0; i < excit; ++i) {
      size_t const bk_rank = i + 1;
      auto bra_tr1s = zip(repeat(0), repeat(nvirt)) | take(bk_rank);
      auto ket_tr1s = zip(repeat(0), repeat(nocc))  | take(bk_rank);
      auto const trange = make_trange(concat(bra_tr1s, ket_tr1s));
      Ts.emplace_back(Tensor_t{world, trange});
      Ts[i].fill(0.);
    }
    //--------
  } // ctor

  Tensor_t operator()(sequant::Tensor const& tensor) const {
    using namespace ranges::views;

    if (tensor.label() == L"t") {
      auto rank = tensor.rank();
      assert(rank <= Ts.size());
      return std::cref(Ts[rank - 1]);
    }

    auto const r1_limits = range1_limits(tensor);
    assert(r1_limits.size() == DataInfo::fock_rank
           || r1_limits.size() == DataInfo::eri_rank);

    auto const iter_limits =  r1_limits
                             | transform([this](auto x) {
                                 return x == nocc ? std::pair{size_t{0}, nocc}
                                                  : std::pair{nocc, nocc + nvirt};
                               });

    auto const tlabel = tensor.label();
    assert(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G_pqrs : F_pq;
    auto slice =
        TA::TArrayD{big_tensor.world(),
                    make_trange(zip(repeat(size_t{0}), r1_limits))};
    slice.fill(0);

    auto tile_orig = big_tensor.find(0).get();
    auto tile_dest = slice.find(0).get();

    if (iter_limits.size() == DataInfo::fock_rank) {
      for (auto ii = iter_limits[0].first; ii < iter_limits[0].second; ++ii)
        for (auto jj = iter_limits[1].first; jj < iter_limits[1].second; ++jj) {
          tile_dest(ii - iter_limits[0].first,  //
                    jj - iter_limits[1].first) = tile_orig(ii, jj);
        }
    } else {  // DataInfo::eri_rank
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

  void update_amplitudes(std::vector<Tensor_t> const& rs) {
    assert(rs.size() == Ts.size()
           && "Unequal number of Rs and Ts!");

    for (auto&& [t, r]: ranges::views::zip(Ts, rs))
      update_single_T(t, r);
  }

  Tensor_t const& amplitude(size_t excitation) const {
    return Ts[excitation-1];
  }

}; // class

} // namespace

#endif  // SEQUANT_EVAL_DATA_WORLD_TA_HPP
