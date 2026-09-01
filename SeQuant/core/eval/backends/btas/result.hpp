#ifndef SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP

#ifdef SEQUANT_HAS_BTAS

#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <btas/btas.h>

#include <complex>

#include <range/v3/view/concat.hpp>
#include <range/v3/view/iota.hpp>

namespace sequant {

// implementation details of the BTAS result backend; prefer sequant::detail
// over an unnamed namespace in a header (see CppCoreGuidelines SF.21 / "Use
// unnamed namespaces in headers ... no" guidance)
namespace detail {

///
/// \brief This function implements the symmetrization of btas::Tensor.
///
/// \param arr The tensor to be symmetrized.
///
/// \pre The rank of the tensor must be even.
///
/// \return The symmetrized btas::Tensor.
///
template <typename... Args>
auto column_symmetrize_btas(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();

  if (rank % 2 != 0)
    throw Exception("This function only supports even-ranked tensors");

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto const lannot = perm;

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm)]() {
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, perm);
    result += temp;
  };

  auto const nparticles = rank / 2;
  symmetric_permutation(SymmetricParticleRange{perm.begin(),               //
                                               perm.begin() + nparticles,  //
                                               nparticles},
                        call_back);
  auto const nf = static_cast<double>(rational{1, factorial(nparticles)});
  btas::scal(nf, result);

  return result;
}

///
/// \brief This function implements the antisymmetrization of btas::Tensor.
///
/// \param arr The tensor to be antisymmetrized
///
/// \param bra_rank The rank of the bra indices
///
/// \return The antisymmetrized btas::Tensor.
///
template <typename... Args>
auto particle_antisymmetrize_btas(btas::Tensor<Args...> const& arr,
                                  size_t bra_rank) {
  using ranges::views::concat;
  using ranges::views::iota;
  size_t const rank = arr.rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;
  const auto lannot = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto process_permutations = [&lannot](const btas::Tensor<Args...>& input_arr,
                                        size_t range_rank, perm_t range_perm,
                                        const perm_t& other_perm, bool is_bra) {
    if (range_rank <= 1) return input_arr;
    btas::Tensor<Args...> result{input_arr.range()};

    auto callback = [&](int parity) {
      const auto annot =
          is_bra ? concat(range_perm, other_perm) | ranges::to<perm_t>()
                 : concat(other_perm, range_perm) | ranges::to<perm_t>();

      typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
      btas::Tensor<Args...> temp;
      btas::permute(input_arr, lannot, temp, annot);
      btas::scal(p_, temp);
      result += temp;
    };

    antisymmetric_permutation(ParticleRange{range_perm.begin(), range_rank},
                              callback);
    return result;
  };
  // Process bra permutations first
  const auto ket_annot = ket_rank == 0 ? perm_t{} : ket_perm;
  auto result = process_permutations(arr, bra_rank, bra_perm, ket_annot, true);

  // Process ket permutations if needed
  const auto bra_annot = bra_rank == 0 ? perm_t{} : bra_perm;
  result = process_permutations(result, ket_rank, ket_perm, bra_annot, false);

  auto const nf = static_cast<double>(
      rational{1, factorial(bra_rank) * factorial(ket_rank)});
  btas::scal(nf, result);

  return result;
}

template <typename... Args>
inline void log_btas(Args const&... args) noexcept {
  log_result("[BTAS] ", args...);
}

/// \brief The batch (Hadamard) axes of a binary contraction: labels that occur
///        in the left operand, the right operand AND the result.
///
/// btas::contract classifies every index as free-left, free-right or contracted
/// (an index shared by both operands is *always* summed); it cannot represent a
/// batch axis -- an index carried through the contraction rather than summed,
/// e.g. an auxiliary/quadrature index in Laplace-transform MP2 or tensor
/// hypercontraction. Such axes are detected here and handled by
/// batched_contract() below.
template <typename Annot>
Annot batch_indices(Annot const& la, Annot const& ra, Annot const& ca) {
  auto has = [](Annot const& v, auto const& x) {
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  Annot batch;
  for (auto const& x : ca)
    if (has(la, x) && has(ra, x)) batch.push_back(x);
  return batch;
}

/// \brief Contracts \p L (labelled \p la) with \p R (labelled \p ra) into a
///        result labelled \p ca in which \p batch labels are Hadamard/batch
///        axes (iterated, not summed).
///
/// Precondition: \p batch is non-empty and equals batch_indices(la, ra, ca).
/// The batch axes are permuted to the front of both operands (making each batch
/// slice a contiguous block in BTAS's row-major storage), and the non-batch
/// remainder of each slice is contracted with plain btas::contract (or, for the
/// degenerate slice shapes, a scaling or a dot). The per-slice results are
/// assembled and permuted back into \p ca order.
template <typename T, typename Annot>
T batched_contract(typename T::numeric_type alpha, T const& L, Annot const& la,
                   T const& R, Annot const& ra, Annot const& ca,
                   Annot const& batch) {
  using numeric_type = typename T::numeric_type;
  auto has = [](Annot const& v, auto const& x) {
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  // non-batch remainder of a label list, preserving order
  auto rest = [&](Annot const& v) {
    Annot r;
    for (auto const& x : v)
      if (!has(batch, x)) r.push_back(x);
    return r;
  };
  Annot const rl = rest(la), rr = rest(ra), rc = rest(ca);

  // permute each operand to [batch..., rest...]
  auto front_batch = [&](Annot const& r) {
    Annot p = batch;
    p.insert(p.end(), r.begin(), r.end());
    return p;
  };
  Annot const la_p = front_batch(rl), ra_p = front_batch(rr);
  T Lp, Rp;
  btas::permute(L, la, Lp, la_p);
  btas::permute(R, ra, Rp, ra_p);

  auto const nb = batch.size();
  auto tail_extents = [](T const& t, std::size_t skip) {
    container::svector<std::size_t> e;
    std::size_t sz = 1;
    for (std::size_t i = skip; i != t.rank(); ++i) {
      e.push_back(t.extent(i));
      sz *= t.extent(i);
    }
    return std::pair{std::move(e), sz};
  };
  // both operands were permuted to [batch..., rest...], so their leading nb
  // axes are the batch axes in the same order; a well-formed annotation gives
  // them identical extents. Guard against a malformed one before slicing (P is
  // taken from the left operand below).
  std::size_t P = 1;
  for (std::size_t i = 0; i != nb; ++i) {
    SEQUANT_ASSERT(Lp.extent(i) == Rp.extent(i) &&
                   "batched_contract: operands disagree on a batch extent");
    P *= Lp.extent(i);
  }
  auto const [rlext, rlsz] = tail_extents(Lp, nb);
  auto const [rrext, rrsz] = tail_extents(Rp, nb);

  // extent of each result-remainder label, looked up from the operands
  container::svector<std::size_t> rcext;
  std::size_t rcsz = 1;
  for (auto const& x : rc) {
    std::size_t ext = 0;
    if (auto it = std::find(rl.begin(), rl.end(), x); it != rl.end())
      ext = rlext[static_cast<std::size_t>(it - rl.begin())];
    else {
      auto jt = std::find(rr.begin(), rr.end(), x);
      SEQUANT_ASSERT(jt != rr.end());
      ext = rrext[static_cast<std::size_t>(jt - rr.begin())];
    }
    rcext.push_back(ext);
    rcsz *= ext;
  }

  // assemble result in [batch..., rc...] layout
  container::svector<std::size_t> cfext;
  for (std::size_t i = 0; i != nb; ++i) cfext.push_back(Lp.extent(i));
  cfext.insert(cfext.end(), rcext.begin(), rcext.end());
  T C{btas::Range{cfext}};

  // extract batch slice b of a [batch..., rest...]-ordered operand as a tensor
  // with the rest extents (never called for an empty remainder -- a scalar
  // slice is read directly, avoiding a rank-0 tensor)
  auto slice = [](T const& src, std::size_t off,
                  container::svector<std::size_t> const& ext, std::size_t sz) {
    T s{btas::Range{ext}};
    std::copy(src.data() + off, src.data() + off + sz, s.data());
    return s;
  };

  for (std::size_t b = 0; b != P; ++b) {
    numeric_type const* lblk = Lp.data() + b * rlsz;
    numeric_type const* rblk = Rp.data() + b * rrsz;

    if (rc.empty()) {
      // every non-batch index is contracted -> a scalar per batch element
      numeric_type d;
      if (rl.empty()) {  // both slices are scalars
        d = lblk[0] * rblk[0];
      } else {
        T Lb = slice(Lp, b * rlsz, rlext, rlsz);
        T Rb = slice(Rp, b * rrsz, rrext, rrsz);
        T Rb2;
        btas::permute(Rb, rr, Rb2, rl);  // align for the dot
        d = btas::dot(Lb, Rb2);
      }
      C.data()[b] = alpha * d;
    } else if (rl.empty() || rr.empty()) {
      // one slice is a scalar -> the other, permuted to rc order and scaled
      numeric_type const s = alpha * (rl.empty() ? lblk[0] : rblk[0]);
      T Cb;
      if (rl.empty())
        btas::permute(slice(Rp, b * rrsz, rrext, rrsz), rr, Cb, rc);
      else
        btas::permute(slice(Lp, b * rlsz, rlext, rlsz), rl, Cb, rc);
      btas::scal(s, Cb);
      std::copy(Cb.data(), Cb.data() + rcsz, C.data() + b * rcsz);
    } else {
      T Lb = slice(Lp, b * rlsz, rlext, rlsz);
      T Rb = slice(Rp, b * rrsz, rrext, rrsz);
      T Cb;
      btas::contract(alpha, Lb, rl, Rb, rr, numeric_type{0}, Cb, rc);
      std::copy(Cb.data(), Cb.data() + rcsz, C.data() + b * rcsz);
    }
  }

  // permute [batch..., rc...] back into the requested result order ca
  T result;
  btas::permute(C, front_batch(rc), result, ca);
  return result;
}

}  // namespace detail

///
/// \brief Result for a tensor value of btas::Tensor type.
/// \tparam T btas::Tensor type. Must be a specialization of btas::Tensor.
///
template <typename T>
class ResultTensorBTAS final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename T::numeric_type;

  explicit ResultTensorBTAS(T arr) : Result{std::move(arr)} {}

 private:
  // TODO make it same as that used by EvalExprBTAS class from eval.hpp file
  using annot_t = container::svector<long>;
  using annot_wrap = Annot<annot_t>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultTensorBTAS<T>>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<ResultTensorBTAS<T>>());
    auto const a = annot_wrap{annot};

    detail::log_btas(detail::ords_to_labels(a.lannot), " + ",
                     detail::ords_to_labels(a.rannot), " = ",
                     detail::ords_to_labels(a.this_annot), "\n");

    T lres, rres;
    btas::permute(get<T>(), a.lannot, lres, a.this_annot);
    btas::permute(other.get<T>(), a.rannot, rres, a.this_annot);
    return eval_result<ResultTensorBTAS<T>>(lres + rres);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto const scalar = other.as<ResultScalar<numeric_type>>().value();
      detail::log_btas(detail::ords_to_labels(a.lannot), " * ", scalar, " = ",
                       detail::ords_to_labels(a.this_annot), "\n");

      T result;
      btas::permute(get<T>(), a.lannot, result, a.this_annot);
      btas::scal(scalar, result);
      return eval_result<ResultTensorBTAS<T>>(std::move(result));
    }

    SEQUANT_ASSERT(other.is<ResultTensorBTAS<T>>());

    if (a.this_annot.empty()) {
      T rres;
      btas::permute(other.get<T>(), a.rannot, rres, a.lannot);
      numeric_type const d = btas::dot(get<T>(), rres);
      detail::log_btas(detail::ords_to_labels(a.lannot), " * ",
                       detail::ords_to_labels(a.rannot), " = ", d, "\n");
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    detail::log_btas(detail::ords_to_labels(a.lannot), " * ",
                     detail::ords_to_labels(a.rannot), " = ",
                     detail::ords_to_labels(a.this_annot), "\n");

    // A batch (Hadamard) axis -- a label shared by both operands AND the result
    // -- is carried through the contraction rather than summed (e.g. an
    // auxiliary/quadrature index in Laplace-MP2 or tensor hypercontraction).
    // btas::contract cannot express one (it always sums an index common to both
    // operands), so route such contractions through batched_contract().
    auto const batch = detail::batch_indices(a.lannot, a.rannot, a.this_annot);
    if (!batch.empty())
      return eval_result<ResultTensorBTAS<T>>(detail::batched_contract<T>(
          numeric_type{1}, get<T>(), a.lannot, other.get<T>(), a.rannot,
          a.this_annot, batch));

    T result;
    btas::contract(numeric_type{1},           //
                   get<T>(), a.lannot,        //
                   other.get<T>(), a.rannot,  //
                   numeric_type{0},           //
                   result, a.this_annot);
    return eval_result<ResultTensorBTAS<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<T>();
    btas::scal(numeric_type(factor), pre);
    return eval_result<ResultTensorBTAS<T>>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);

    detail::log_btas(detail::ords_to_labels(pre_annot), " = ",
                     detail::ords_to_labels(post_annot), "\n");

    T result;
    btas::permute(get<T>(), pre_annot, result, post_annot);
    return eval_result<ResultTensorBTAS<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr apply_transform(
      CanonTransform t, std::array<std::any, 2> const& ann) const override {
    // phase * (elementwise conj) * permute in one pass; conj elided for a
    // real numeric_type (equal annots => the permute is a copy)
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    detail::log_btas(detail::ords_to_labels(post_annot), " = apply_transform(",
                     detail::ords_to_labels(pre_annot), ")\n");
    T result;
    btas::permute(get<T>(), pre_annot, result, post_annot);
    if constexpr (meta::is_complex_v<numeric_type>) {
      if (t.conj)
        for (auto& x : result) x = std::conj(x);
    }
    if (t.phase != 1)
      for (auto& x : result) x = x * numeric_type(t.phase);
    return eval_result<ResultTensorBTAS<T>>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    auto& t = get<T>();
    auto const& o = other.get<T>();
    SEQUANT_ASSERT(t.range() == o.range());

    auto const ann = detail::ords_to_annot(
        ranges::views::iota(size_t{0}, static_cast<size_t>(t.rank())));
    detail::log_btas(ann, " += ", ann, "\n");

    t += o;
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultTensorBTAS<T>>(
        detail::column_symmetrize_btas(get<T>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<ResultTensorBTAS<T>>(
        detail::particle_antisymmetrize_btas(get<T>(), bra_rank));
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    const auto& tensor = get<T>();
    // only count data
    return tensor.range().volume() * sizeof(T);
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_BTAS

#endif  // SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP
