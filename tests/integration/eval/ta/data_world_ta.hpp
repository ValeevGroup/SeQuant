//
// Created by Bimal Gaudel on 7/14/21.
//

#ifndef SEQUANT_EVAL_DATA_WORLD_TA_HPP
#define SEQUANT_EVAL_DATA_WORLD_TA_HPP

#include <data_info.hpp>
#include <eval_utils.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/eval/result.hpp>

#include <tiledarray.h>
#include <range/v3/view.hpp>

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

  container::vector<ResultPtr> Ts;

  mutable container::map<size_t, ResultPtr> cache_;

  ///
  /// Read tensor data from a file to a TA::DistArray tensor.
  ///
  /// \param fname Input file name
  /// \param tensor Output TA::DistArray tensor.
  void read_tensor_ta(std::string_view fname, Tensor_t& tensor) {
    // TODO assert tensor single tiled or extend to handle multiply tiled case
    auto ta_tensor = TA::Tensor<typename Tensor_t::scalar_type>{
        tensor.trange().make_tile_range(0)};
    read_tensor(fname, ta_tensor);
    *tensor.begin() = ta_tensor;
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
    auto tr1s = pairs | transform([](auto&& pair) {
                  auto const [lo, hi] = pair;
                  return TA::TiledRange1{lo, hi};
                }) |
                ranges::to_vector;

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
    assert(T.trange().rank() == R.trange().rank() &&
           "Incompatible ranks of R and T");

    assert(T.trange().rank() % 2 == 0 && "Odd rank not supported");

    size_t const n = T.trange().rank() / 2;

    auto tile_T = T.find(0).get();
    auto tile_R = R.find(0).get();
    auto tile_F = F_pq.find(0).get();
    auto updater = [&tile_T, &tile_R, &tile_F, n, this](auto const& idx_vec) {
      double diag = 0.;
      for (auto x : idx_vec | take(n)) diag -= tile_F(x + nocc, x + nocc);
      for (auto x : idx_vec | drop(n) | take(n)) diag += tile_F(x, x);

      tile_T(idx_vec) += tile_R(idx_vec) / diag;
    };

    auto virt_occ =
        concat(repeat_n(iota(size_t{0}, nvirt) | ranges::to_vector, n),
               repeat_n(iota(size_t{0}, nocc) | ranges::to_vector, n)) |
        ranges::to_vector;

    cartesian_foreach(virt_occ, updater);
  }

 public:
  DataWorldTA(TA::World& world, DataInfo const& info, size_t excit)
      : nocc{info.nocc()},
        nvirt{info.nvirt()},
        fock_file{info.fock_file()},
        eri_file{info.eri_file()} {
    using namespace ranges::views;

    //--------------------------------------------------------
    // read ERI and Fock tensors
    //--------------------------------------------------------
    auto const nobs = nocc + nvirt;
    auto const fock_trange =
        make_trange(zip(repeat(0), repeat(nobs)) | take(DataInfo::fock_rank));
    auto const eri_trange =
        make_trange(zip(repeat(0), repeat(nobs)) | take(DataInfo::eri_rank));
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
      auto ket_tr1s = zip(repeat(0), repeat(nocc)) | take(bk_rank);
      auto const trange = make_trange(concat(bra_tr1s, ket_tr1s));
      Ts.emplace_back(
          eval_result<ResultTensorTA<Tensor_t>>(Tensor_t{world, trange}));
      auto& t = Ts[i]->template get<Tensor_t>();
      t.fill(0);
    }
    //--------
  }  // ctor

  Tensor_t operator()(sequant::Tensor const& tensor) const {
    using namespace ranges::views;

    assert(tensor.label() != L"t");

    auto const r1_limits = range1_limits(tensor, nocc, nvirt);
    assert(r1_limits.size() == DataInfo::fock_rank ||
           r1_limits.size() == DataInfo::eri_rank);

    auto const iter_limits = r1_limits | transform([this](auto x) {
                               return x == nocc ? std::pair{size_t{0}, nocc}
                                                : std::pair{nocc, nocc + nvirt};
                             });

    auto const tlabel = tensor.label();
    assert(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G_pqrs : F_pq;
    auto slice = TA::TArrayD{big_tensor.world(),
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

  ResultPtr operator()(sequant::meta::can_evaluate auto const& n) const {
    using numeric_type = typename Tensor_t::numeric_type;
    if (n->result_type() == ResultType::Scalar) {
      assert(n->expr()->template is<Constant>());
      auto d = n->as_constant().template value<numeric_type>();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    assert(n->result_type() == ResultType::Tensor &&
           n->expr()->template is<Tensor>());

    if (auto t = n->as_tensor(); t.label() == L"t") {
      auto rank = t.rank();
      assert(rank <= Ts.size());
      return Ts[rank - 1];
    }

    auto h = hash::value(*n);
    if (auto exists = cache_.find(h); exists != cache_.end())
      return exists->second;
    else {
      auto tnsr =
          eval_result<ResultTensorTA<Tensor_t>>((*this)(n->as_tensor()));
      auto stored = cache_.emplace(h, std::move(tnsr));
      assert(stored.second && "failed to store tensor");
      return stored.first->second;
    }
  }

  void update_amplitudes(std::vector<Tensor_t> const& rs) {
    using ranges::views::transform;
    using ranges::views::zip;

    assert(rs.size() == Ts.size() && "Unequal number of Rs and Ts!");

    for (auto&& [t, r] : zip(Ts | transform([](auto&& res) -> Tensor_t& {
                               return res->template get<Tensor_t>();
                             }),
                             rs))
      update_single_T(t, r);
  }

  Tensor_t const& amplitude(size_t excitation) const {
    return Ts[excitation - 1]->template get<Tensor_t>();
  }

};  // class

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_DATA_WORLD_TA_HPP
