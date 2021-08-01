//
// Created by Bimal Gaudel on 7/31/21.
//

#ifndef SEQUANT_EVAL_DATA_WORLD_BTAS_HPP
#define SEQUANT_EVAL_DATA_WORLD_BTAS_HPP

#include "examples/eval/data_info.hpp"
#include "examples/eval/eval_utils.hpp"

#include <btas/btas.h>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <range/v3/view.hpp>

namespace sequant::eval::btas {

template <typename Tensor_t>
class DataWorldBTAS {
 private:
  size_t const nocc;

  size_t const nvirt;

  std::string const fock_file;

  std::string const eri_file;

  Tensor_t G_pqrs;

  Tensor_t F_pq;

  container::vector<Tensor_t> Ts;

 public:
  DataWorldBTAS(size_t excit, DataInfo const& info)
      : nocc{info.nocc()},
        nvirt{info.nvirt()},
        fock_file{info.fock_file()},
        eri_file{info.eri_file()} {
    using namespace ranges::views;

    //--------------------------------------------------------
    // read ERI and Fock tensors
    //--------------------------------------------------------
    auto const nobs = nocc + nvirt;
    F_pq = Tensor_t{
        ::btas::Range{repeat_n(nobs, DataInfo::fock_rank) | ranges::to_vector}};
    F_pq.fill(0);
    G_pqrs = Tensor_t{
        ::btas::Range{repeat_n(nobs, DataInfo::eri_rank) | ranges::to_vector}};
    G_pqrs.fill(0);
    read_tensor(fock_file, F_pq);
    read_tensor(eri_file, G_pqrs);
    //--------

    //--------------------------------------------------------
    // Initialize amplitude tensors
    //--------------------------------------------------------
    Ts.reserve(excit);
    for (auto i = 0; i < excit; ++i) {
      size_t const bk_rank = i + 1;

      auto range = ::btas::Range{
          concat(repeat_n(nvirt, bk_rank), repeat_n(nocc, bk_rank)) |
          ranges::to_vector};
      Ts.emplace_back(Tensor_t{range});
      Ts[i].fill(0);
    }
    //--------
  }  // ctor

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
    assert(T.rank() == R.rank() && "Incompatible ranks of R and T");

    assert(T.rank() % 2 == 0 && "Odd rank not supported");

    size_t const n = T.rank() / 2;

    auto const& F = F_pq;
    auto updater = [&T, &R, &F, n, this](auto const& idx_vec) {
      double diag = 0.;
      for (auto x : idx_vec | take(n)) diag -= F(x + nocc, x + nocc);
      for (auto x : idx_vec | drop(n) | take(n)) diag += F(x, x);

      T(idx_vec) += R(idx_vec) / diag;
    };

    auto virt_occ =
        concat(repeat_n(iota(size_t{0}, nvirt) | ranges::to_vector, n),
               repeat_n(iota(size_t{0}, nocc) | ranges::to_vector, n)) |
        ranges::to_vector;

    cartesian_foreach(virt_occ, updater);
  }

  Tensor_t operator()(sequant::Tensor const& tensor) const {
    using namespace ranges::views;

    if (tensor.label() == L"t") {
      auto rank = tensor.rank();
      assert(rank <= Ts.size());
      return std::cref(Ts[rank - 1]);
    }

    auto const r1_limits = range1_limits(tensor, nocc, nvirt);
    assert(r1_limits.size() == DataInfo::fock_rank ||
           r1_limits.size() == DataInfo::eri_rank);

    auto const iter_limits =
        r1_limits | ranges::views::transform([this](auto x) {
          return x == nocc ? std::pair{size_t{0}, nocc}
                           : std::pair{nocc, nocc + nvirt};
        });

    auto const tlabel = tensor.label();
    assert(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G_pqrs : F_pq;

    auto slice = Tensor_t{::btas::Range{r1_limits | ranges::to_vector}};
    if (iter_limits.size() == DataInfo::fock_rank) {
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          slice(i - loop1.first, j - loop2.first) = big_tensor(i, j);

    } else {  // DataInfo::eri_rank
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

  void update_amplitudes(std::vector<Tensor_t> const& rs) {
    assert(rs.size() == Ts.size() && "Unequal number of Rs and Ts!");

    for (auto&& [t, r] : ranges::views::zip(Ts, rs)) update_single_T(t, r);
  }

  Tensor_t const& amplitude(size_t excitation) const {
    return Ts[excitation - 1];
  }

};  // class

}  // namespace sequant::eval::btas

#endif  // SEQUANT_DATA_WORLD_BTAS_HPP
