//
// Created by Eduard Valeyev on 3/31/18.
//

#ifndef SEQUANT2_WICK_IMPL_HPP
#define SEQUANT2_WICK_IMPL_HPP

namespace sequant2 {

namespace detail {

struct zero_result : public std::exception {};

/// @brief computes index replacement rules

/// If using orthonormal representation, overlaps are Kronecker deltas, hence
/// summations can be reduced by index replacements. Reducing sums over dummy
/// (internal) indices uses 2 rules:
/// - if a Kronecker delta binds 2 internal indices I and J, replace them with a
/// new internal index
///   representing intersection of spaces of I and J, !!remove delta!!
/// - if a Kronecker delta binds an internal index J and an external index I:
///   - if space of J includes space of I, replace J with I, !!remove delta!!
///   - if space of J is a subset of space of I, replace J with a new internal
///     index representing intersection of spaces of I and J, !!keep the delta!!
/// @throw zero_result if @c product is zero
inline std::map<Index, Index> compute_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::vector<Index> &external_indices) {
  expr_range exrng(product);

  IndexFactory idxfac;
  std::map<Index /* src */, Index /* dst */> result;  // src->dst

  auto make_intersection_index = [&idxfac](const IndexSpace &space1,
                                           const IndexSpace &space2) {
    const auto intersection_space = intersection(space1, space2);
    if (intersection_space == IndexSpace::null_instance()) throw zero_result{};
    return idxfac.make(intersection_space);
  };

  /// adds src->dst or src->intersection(dst,current_dst)
  auto add_rule = [&result, make_intersection_index](
      const Index &src, const Index &dst) {
    auto src_it = result.find(src);
    if (src_it == result.end())  // if brand new, add the rule
      result[src] = dst;
    else {  // else modify the destination of the existing rule to the
      // intersection
      const auto &old_dst = src_it->second;
      if (dst.space() != old_dst.space()) {
        src_it->second = make_intersection_index(old_dst.space(), dst.space());
      }
    }
  };

  /// adds src1->dst and src2->dst; if src1->dst1 and/or src2->dst2 already
  /// exist existing rules are updated to map to the intersection of dst1, dst2
  /// and dst
  auto add_rules = [&result, &idxfac, make_intersection_index](
      const Index &src1, const Index &src2, const Index &dst) {
    auto src1_it = result.find(src1);
    auto src2_it = result.find(src2);
    const auto has_src1_rule = src1_it != result.end();
    const auto has_src2_rule = src2_it != result.end();
    if (!has_src1_rule && !has_src2_rule) {  // if brand new, add the rules
      result.insert({src1, dst});
      assert(!result.empty());
      result.insert({src2, dst});
      assert(!result.empty());
    } else if (has_src1_rule && !has_src2_rule) {
      const auto &old_dst1 = src1_it->second;
      if (dst.space() != old_dst1.space()) {
        src1_it->second =
            make_intersection_index(old_dst1.space(), dst.space());
      }
      result[src2] = src1_it->second;
    } else if (!has_src1_rule && has_src2_rule) {
      const auto &old_dst2 = src2_it->second;
      if (dst.space() != old_dst2.space()) {
        src2_it->second =
            make_intersection_index(old_dst2.space(), dst.space());
      }
      result[src1] = src2_it->second;
    } else {  //
      const auto &old_dst1 = src1_it->second;
      const auto &old_dst2 = src2_it->second;
      const auto new_dst_space =
          (dst.space() != old_dst1.space() || dst.space() != old_dst2.space())
          ? intersection(old_dst1.space(), old_dst2.space(), dst.space())
          : dst.space();
      if (new_dst_space == IndexSpace::null_instance()) throw zero_result{};
      Index new_dst;
      if (new_dst_space == old_dst1.space()) {
        new_dst = old_dst1;
        if (new_dst_space == old_dst2.space() && old_dst2 < new_dst) {
          new_dst = old_dst2;
        }
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == old_dst2.space()) {
        new_dst = old_dst2;
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == dst.space()) {
        new_dst = dst;
      } else
        new_dst = idxfac.make(new_dst_space);
      result[src1] = new_dst;
      result[src2] = new_dst;
    }
  };

  /// this makes the list of replacements ... we do not mutate the expressions
  /// to keep the information about which indices are related
  for (auto it = ranges::begin(exrng); it != ranges::end(exrng);
       ++it) {
    const auto &factor = *it;
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      const auto &tensor = static_cast<const Tensor &>(*factor);
      if (tensor.label() == L"S") {
        assert(tensor.bra().size() == 1);
        assert(tensor.ket().size() == 1);
        const auto &bra = tensor.bra().at(0);
        const auto &ket = tensor.ket().at(0);
        assert(bra != ket);

        const auto bra_is_ext = ranges::find(external_indices, bra) !=
            ranges::end(external_indices);
        const auto ket_is_ext = ranges::find(external_indices, ket) !=
            ranges::end(external_indices);

        const auto intersection_space = intersection(bra.space(), ket.space());
        assert(intersection_space != IndexSpace::null_instance());

        if (!bra_is_ext && !ket_is_ext) {  // int + int
          const auto new_dummy = idxfac.make(intersection_space);
          add_rules(bra, ket, new_dummy);
        } else if (bra_is_ext && !ket_is_ext) {  // ext + int
          if (includes(bra.space(), ket.space())) {
            add_rule(ket, bra);
          } else {
            add_rule(ket, idxfac.make(intersection_space));
          }
        } else if (!bra_is_ext && ket_is_ext) {  // int + ext
          if (includes(ket.space(), bra.space())) {
            add_rule(bra, ket);
          } else {
            add_rule(bra, idxfac.make(intersection_space));
          }
        }
      }
    }
  }

  return result;
}

inline void apply_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const std::map<Index, Index> &const_replrules,
    const container::vector<Index> &external_indices) {
  // to be able to use map[]
  auto &replrules = const_cast<std::map<Index, Index> &>(const_replrules);

  expr_range exrng(product);

  for (auto it = ranges::begin(exrng); it != ranges::end(exrng);) {
    const auto &factor = *it;
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      bool erase_it = false;
      auto &tensor = static_cast<Tensor &>(*factor);

      /// this removes the deltas that are no longer needed
      /// this also ensures that a single pass was sufficient to construct
      /// self-consistent replacement rules
      if (tensor.label() == L"S") {
        const auto &bra = tensor.bra().at(0);
        const auto &ket = tensor.ket().at(0);

        const auto bra_is_ext = ranges::find(external_indices, bra) !=
            ranges::end(external_indices);
        const auto ket_is_ext = ranges::find(external_indices, ket) !=
            ranges::end(external_indices);

#ifndef NDEBUG
        const auto intersection_space = intersection(bra.space(), ket.space());
#endif

        if (!bra_is_ext && !ket_is_ext) {  // int + int
#ifndef NDEBUG
          assert(replrules[bra].space() == replrules[ket].space());
#endif
          erase_it = true;
        } else if (bra_is_ext && !ket_is_ext) {  // ext + int
          if (includes(bra.space(), ket.space())) {
#ifndef NDEBUG
            assert(replrules[ket].space() == bra.space());
#endif
            erase_it = true;
          } else {
#ifndef NDEBUG
            assert(replrules[ket].space() == intersection_space);
#endif
          }
        } else if (!bra_is_ext && ket_is_ext) {  // int + ext
          if (includes(ket.space(), bra.space())) {
#ifndef NDEBUG
            assert(replrules[bra].space() == ket.space());
#endif
            erase_it = true;
          } else {
#ifndef NDEBUG
            assert(replrules[bra].space() == intersection_space);
#endif
          }
        }

        if (erase_it) {
          *it = make<Constant>(1);
        }
      }  // Kronecker delta
      else {
        /// replace
        tensor.transform_indices(const_replrules);
      }

    }
    ++it;
  }
}

/// If using orthonormal representation, resolves Kronecker deltas (=overlaps
/// between indices in orthonormal spaces) in summations
/// @throw zero_result if @c expr is zero
inline void reduce_wick_impl(std::shared_ptr<Product> &expr,
                             const container::vector<Index> &external_indices) {
  if (get_default_context().metric() == IndexSpaceMetric::Unit) {
    const auto replacement_rules =
        compute_index_replacement_rules(expr, external_indices);
//    std::wcout << "reduce_wick_impl(expr, external_indices):\n  expr = " << expr->to_latex() << "\n  external_indices = ";
    ranges::for_each(external_indices, [](auto &index) { std::wcout << index.label() << " "; });
//    std::wcout << "\n  replrules = ";
    ranges::for_each(replacement_rules,
                     [](auto &index) { std::wcout << index.first.label() << "->" << index.second.label() << " "; });
    std::wcout.flush();

    apply_index_replacement_rules(expr, replacement_rules, external_indices);
//    std::wcout << "\n  result = " << expr->to_latex() << std::endl;
  } else
    abort();  // programming error?
}

}  // namespace detail

template<Statistics S>
void WickTheorem<S>::reduce(ExprPtr &expr) const {
  // there are 2 possibilities: expr is a single Product, or it's a Sum of
  // Products
  if (expr->type_id() == Expr::get_type_id<Product>()) {
    auto expr_cast = std::static_pointer_cast<Product>(expr);
    try {
      detail::reduce_wick_impl(expr_cast, external_indices_);
      expr = expr_cast;
    } catch (detail::zero_result &) {
      expr = std::make_shared<Constant>(0);
    }
  } else {
    assert(expr->type_id() == Expr::get_type_id<Sum>());
    for (auto &subexpr : *expr) {
      assert(subexpr->type_id() == Expr::get_type_id<Product>());
      auto subexpr_cast = std::static_pointer_cast<Product>(subexpr);
      try {
        detail::reduce_wick_impl(subexpr_cast, external_indices_);
        subexpr = subexpr_cast;
      }
      catch (detail::zero_result &) {
        subexpr = std::make_shared<Constant>(0);
      }
    }
  }
}

}  // namespace sequant2

#endif  // SEQUANT2_WICK_IMPL_HPP
