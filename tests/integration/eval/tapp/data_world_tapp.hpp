//
// TAPP backend data world for integration tests.
// Mirrors btas/data_world_btas.hpp
//

#ifndef SEQUANT_EVAL_DATA_WORLD_TAPP_HPP
#define SEQUANT_EVAL_DATA_WORLD_TAPP_HPP

#include <data_info.hpp>
#include <eval_utils.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/tapp/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tapp/ops.hpp>
#include <SeQuant/core/eval/backends/tapp/result.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <range/v3/view.hpp>

namespace sequant::eval::tapp {

template <typename Tensor_t>
class DataWorldTAPP {
 private:
  size_t const nocc;

  size_t const nvirt;

  std::string const fock_file;

  std::string const eri_file;

  Tensor_t G_pqrs;

  Tensor_t F_pq;

  container::vector<ResultPtr> Ts;

  mutable container::map<size_t, ResultPtr> cache_;

 public:
  DataWorldTAPP(DataInfo const& info, size_t excit)
      : nocc{info.nocc()},
        nvirt{info.nvirt()},
        fock_file{info.fock_file()},
        eri_file{info.eri_file()} {
    using namespace ranges::views;

    //--------------------------------------------------------
    // read ERI and Fock tensors
    //--------------------------------------------------------
    auto const nobs = nocc + nvirt;
    F_pq = Tensor_t{repeat_n(static_cast<int64_t>(nobs), DataInfo::fock_rank) |
                    ranges::to<container::svector<int64_t>>};
    F_pq.fill(0);
    G_pqrs = Tensor_t{repeat_n(static_cast<int64_t>(nobs), DataInfo::eri_rank) |
                      ranges::to<container::svector<int64_t>>};
    G_pqrs.fill(0);
    read_tensor(fock_file, F_pq);
    read_tensor(eri_file, G_pqrs);
    //--------

    //--------------------------------------------------------
    // Initialize amplitude tensors
    //--------------------------------------------------------
    Ts.reserve(excit);
    for (size_t i = 0; i < excit; ++i) {
      size_t const bk_rank = i + 1;

      auto extents = concat(repeat_n(static_cast<int64_t>(nvirt), bk_rank),
                            repeat_n(static_cast<int64_t>(nocc), bk_rank)) |
                     ranges::to<container::svector<int64_t>>;
      Ts.emplace_back(
          eval_result<ResultTensorTAPP<Tensor_t>>(Tensor_t{extents}));
      auto& t = Ts[i]->template get<Tensor_t>();
      t.fill(0);
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
    SEQUANT_ASSERT(T.rank() == R.rank() && "Incompatible ranks of R and T");

    SEQUANT_ASSERT(T.rank() % 2 == 0 && "Odd rank not supported");

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

    SEQUANT_ASSERT(tensor.label() != L"t");

    auto r1_limits = range1_limits(tensor, nocc, nvirt);
    SEQUANT_ASSERT(r1_limits.size() == DataInfo::fock_rank ||
                   r1_limits.size() == DataInfo::eri_rank);

    auto const iter_limits =
        r1_limits | ranges::views::transform([this](auto x) {
          return x == nocc ? std::pair{size_t{0}, nocc}
                           : std::pair{nocc, nocc + nvirt};
        });

    auto const tlabel = tensor.label();
    SEQUANT_ASSERT(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G_pqrs : F_pq;

    auto r1_extents = r1_limits | ranges::views::transform([](size_t v) {
                        return static_cast<int64_t>(v);
                      }) |
                      ranges::to<container::svector<int64_t>>;

    auto slice = Tensor_t{r1_extents};
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

  ResultPtr operator()(meta::can_evaluate auto const& n) const {
    using numeric_type = typename Tensor_t::numeric_type;
    if (n->result_type() == ResultType::Scalar) {
      SEQUANT_ASSERT(n->expr()->template is<Constant>());
      auto d = n->as_constant().template value<numeric_type>();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    SEQUANT_ASSERT(n->result_type() == ResultType::Tensor &&
                   n->expr()->template is<Tensor>());

    if (auto t = n->as_tensor(); t.label() == L"t") {
      auto rank = t.rank();
      SEQUANT_ASSERT(rank <= Ts.size());
      return Ts[rank - 1];
    }
    auto h = hash::value(*n);
    if (auto exists = cache_.find(h); exists != cache_.end())
      return exists->second;
    else {
      auto tnsr =
          eval_result<ResultTensorTAPP<Tensor_t>>((*this)(n->as_tensor()));
      auto stored = cache_.emplace(h, std::move(tnsr));
      SEQUANT_ASSERT(stored.second && "failed to store tensor");
      return stored.first->second;
    }
  }

  void update_amplitudes(std::vector<Tensor_t> const& rs) {
    using ranges::views::transform;
    using ranges::views::zip;

    SEQUANT_ASSERT(rs.size() == Ts.size() && "Unequal number of Rs and Ts!");

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

}  // namespace sequant::eval::tapp

#endif  // SEQUANT_EVAL_DATA_WORLD_TAPP_HPP
