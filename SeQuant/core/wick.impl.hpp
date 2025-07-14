//
// Created by Eduard Valeyev on 3/31/18.
//

#ifndef SEQUANT_WICK_IMPL_HPP
#define SEQUANT_WICK_IMPL_HPP

// change to 1 to try TNV2
#define USE_TENSOR_NETWORK_V2 0

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#if USE_TENSOR_NETWORK_V2
#include <SeQuant/core/tensor_network_v2.hpp>
#else
#include <SeQuant/core/tensor_network.hpp>
#endif
#include <SeQuant/core/tensor_network/vertex.hpp>

#ifdef SEQUANT_HAS_EXECUTION_HEADER
#include <execution>
#endif

namespace sequant {

namespace detail {

struct zero_result : public std::exception {};

/// @brief computes index replacement rules

/// If using orthonormal representation, overlaps are Kronecker deltas, hence
/// summations can be reduced by index replacements. Reducing sums over dummy
/// (internal) indices uses 2 rules:
/// - if a Kronecker delta binds 2 internal indices I and J, replace them with a
///   new internal index representing intersection of spaces of I and J,
///   !!remove delta!!
/// - if a Kronecker delta binds an internal index J and an external index I:
///   - if space of J includes space of I, replace J with I, !!remove delta!!
///   - if space of J is a subset of space of I, replace J with a new internal
///     index representing intersection of spaces of I and J, !!keep the delta!!
/// @throw zero_result if @c product is zero for any reason, e.g. because
///        it includes an overlap of 2 indices from nonoverlapping spaces
template <Statistics S>
container::map<Index, Index> compute_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::set<Index> &external_indices,
    const std::set<Index, Index::LabelCompare> &all_indices,
    const std::shared_ptr<const IndexSpaceRegistry> &isr =
        get_default_context(S).index_space_registry()) {
  expr_range exrng(product);

  /// this ensures that all temporary indices have unique *labels* (not just
  /// unique *full labels*)
  auto index_validator = [&all_indices](const Index &idx) {
    return all_indices.find(idx) == all_indices.end();
  };
  IndexFactory idxfac(index_validator);
  container::map<Index /* src */, Index /* dst */> result;  // src->dst

  // computes an index in intersection of space1 and space2
  auto make_intersection_index = [&idxfac, &isr](const IndexSpace &space1,
                                                 const IndexSpace &space2) {
    const auto intersection_space = isr->intersection(space1, space2);
    if (!intersection_space) throw zero_result{};
    return idxfac.make(intersection_space);
  };

  // transfers proto indices from idx (if any) to img
  auto proto = [](const Index &img, const Index &idx) {
    if (idx.has_proto_indices()) {
      if (img.has_proto_indices()) {
        assert(img.proto_indices() == idx.proto_indices());
        return img;
      } else
        return Index(img, idx.proto_indices());
    } else {
      assert(!img.has_proto_indices());
      return img;
    }
  };

  // adds src->dst or src->intersection(dst,current_dst)
  auto add_rule = [&result, &proto, &make_intersection_index](
                      const Index &src, const Index &dst) {
    auto src_it = result.find(src);
    if (src_it == result.end()) {  // if brand new, add the rule
      auto insertion_result = result.emplace(src, proto(dst, src));
      assert(insertion_result.second);
    } else {  // else modify the destination of the existing rule to the
      // intersection
      const auto &old_dst = src_it->second;
      assert(old_dst.proto_indices() == src.proto_indices());
      if (dst.space() != old_dst.space()) {
        src_it->second =
            proto(make_intersection_index(old_dst.space(), dst.space()), src);
      }
    }
  };

  // adds src1->dst and src2->dst; if src1->dst1 and/or src2->dst2 already
  // exist the existing rules are updated to map to the intersection of dst1,
  // dst2 and dst
  auto add_rules = [&result, &idxfac, &proto, &make_intersection_index, &isr](
                       const Index &src1, const Index &src2, const Index &dst) {
    // are there replacement rules already for src{1,2}?
    auto src1_it = result.find(src1);
    auto src2_it = result.find(src2);
    const auto has_src1_rule = src1_it != result.end();
    const auto has_src2_rule = src2_it != result.end();

    // which proto-indices should dst1 and dst2 inherit? a source index without
    // proto indices will inherit its source counterpart's indices, unless it
    // already has its own protoindices: <a_ij|p> = <a_ij|a_ij> (hence replace p
    // with a_ij), but <a_ij|p_kl> = <a_ij|a_kl> != <a_ij|a_ij> (hence replace
    // p_kl with a_kl)
    const auto &dst1_proto =
        !src1.has_proto_indices() && src2.has_proto_indices() ? src2 : src1;
    const auto &dst2_proto =
        !src2.has_proto_indices() && src1.has_proto_indices() ? src1 : src2;

    if (!has_src1_rule && !has_src2_rule) {  // if brand new, add the rules
      auto insertion_result1 = result.emplace(src1, proto(dst, dst1_proto));
      assert(insertion_result1.second);
      auto insertion_result2 = result.emplace(src2, proto(dst, dst2_proto));
      assert(insertion_result2.second);
    } else if (has_src1_rule &&
               !has_src2_rule) {  // update the existing rule for src1
      const auto &old_dst1 = src1_it->second;
      assert(old_dst1.proto_indices() == dst1_proto.proto_indices());
      if (dst.space() != old_dst1.space()) {
        src1_it->second = proto(
            make_intersection_index(old_dst1.space(), dst.space()), dst1_proto);
      }
      result.emplace(src2, src1_it->second);
    } else if (!has_src1_rule &&
               has_src2_rule) {  // update the existing rule for src2
      const auto &old_dst2 = src2_it->second;
      assert(old_dst2.proto_indices() == dst2_proto.proto_indices());
      if (dst.space() != old_dst2.space()) {
        src2_it->second = proto(
            make_intersection_index(old_dst2.space(), dst.space()), dst2_proto);
      }
      result.emplace(src1, src2_it->second);
    } else {  // update both of the existing rules
      const auto &old_dst1 = src1_it->second;
      const auto &old_dst2 = src2_it->second;
      const auto new_dst_space =
          (dst.space() != old_dst1.space() || dst.space() != old_dst2.space())
              ? isr->intersection(
                    isr->intersection(old_dst1.space(), old_dst2.space()),
                    dst.space())
              : dst.space();
      if (!new_dst_space) throw zero_result{};
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
      result.emplace(src1, proto(new_dst, dst1_proto));
      result.emplace(src2, proto(new_dst, dst2_proto));
    }
  };

  /// this makes the list of replacements ... we do not mutate the expressions
  /// to keep the information about which indices are related
  for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
    const auto &factor = *it;
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      const auto &tensor = static_cast<const Tensor &>(*factor);
      if (tensor.label() == overlap_label()) {
        assert(tensor.bra().size() == 1);
        assert(tensor.ket().size() == 1);
        const auto &bra = tensor.bra().at(0);
        const auto &ket = tensor.ket().at(0);
        assert(bra != ket);

        const auto bra_is_ext = ranges::find(external_indices, bra) !=
                                ranges::end(external_indices);
        const auto ket_is_ext = ranges::find(external_indices, ket) !=
                                ranges::end(external_indices);

        const auto intersection_space =
            isr->intersection(bra.space(), ket.space());

        // if overlap's indices are from non-overlapping spaces, return zero
        if (!intersection_space) {
          throw zero_result{};
        }

        if (!bra_is_ext && !ket_is_ext) {  // int + int
          const auto new_dummy = idxfac.make(intersection_space);
          add_rules(bra, ket, new_dummy);
        } else if (bra_is_ext && !ket_is_ext) {  // ext + int
          if (includes(ket.space(), bra.space())) {
            add_rule(ket, bra);
          } else {
            add_rule(ket, idxfac.make(intersection_space));
          }
        } else if (!bra_is_ext && ket_is_ext) {  // int + ext
          if (includes(bra.space(), ket.space())) {
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

/// @return true if made any changes
inline bool apply_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::map<Index, Index> &const_replrules,
    const container::set<Index> &external_indices,
    std::set<Index, Index::LabelCompare> &all_indices,
    const std::shared_ptr<const IndexSpaceRegistry> &isr) {
  // to be able to use map[]
  auto &replrules = const_cast<container::map<Index, Index> &>(const_replrules);

  expr_range exrng(product);

  /// this recursively applies replacement rules until result does not
  /// change removes the deltas that are no longer needed
#ifndef NDEBUG
  // assert that tensors_ indices are not tagged since going to tag indices
  {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        auto &tensor = factor->as<Tensor>();
        assert(ranges::none_of(tensor.const_indices(), [](const Index &idx) {
          return idx.tag().has_value();
        }));
      }
    }
  }
#endif
  bool mutated = false;
  bool pass_mutated = false;
  do {
    pass_mutated = false;

    for (auto it = ranges::begin(exrng); it != ranges::end(exrng);) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        bool erase_it = false;
        auto &tensor = factor->as<Tensor>();

        /// replace indices
        pass_mutated &= tensor.transform_indices(const_replrules);

        if (tensor.label() == overlap_label()) {
          const auto &bra = tensor.bra().at(0);
          const auto &ket = tensor.ket().at(0);

          if (bra.proto_indices() == ket.proto_indices()) {
            const auto bra_is_ext = ranges::find(external_indices, bra) !=
                                    ranges::end(external_indices);
            const auto ket_is_ext = ranges::find(external_indices, ket) !=
                                    ranges::end(external_indices);

#ifndef NDEBUG
            const auto intersection_space =
                isr->intersection(bra.space(), ket.space());
#endif

            if (!bra_is_ext && !ket_is_ext) {  // int + int
#ifndef NDEBUG
              if (replrules.find(bra) != replrules.end() &&
                  replrules.find(ket) != replrules.end())
                assert(replrules[bra].space() == replrules[ket].space());
#endif
              erase_it = true;
            } else if (bra_is_ext && !ket_is_ext) {  // ext + int
              if (isr->intersection(ket.space(), bra.space()) !=
                  IndexSpace::null) {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == bra.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == intersection_space);
#endif
              }
            } else if (!bra_is_ext && ket_is_ext) {  // int + ext
              if (isr->intersection(bra.space(), ket.space()) !=
                  IndexSpace::null) {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == ket.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == intersection_space);
#endif
              }
            } else {  // ext + ext
              if (bra == ket) erase_it = true;
            }

            if (erase_it) {
              pass_mutated = true;
              *it = ex<Constant>(1);
            }
          }  // matching proto indices
        }    // Kronecker delta
      }
      ++it;
    }
    mutated |= pass_mutated;
  } while (pass_mutated);  // keep replacing til fixed point

  // assert that tensors_ indices are not tagged since going to tag indices
  {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        factor->as<Tensor>().reset_tags();
      }
    }
  }

  // update all_indices
  std::set<Index, Index::LabelCompare> all_indices_new;
  ranges::for_each(
      all_indices, [&const_replrules, &all_indices_new](const Index &idx) {
        auto dst_it = const_replrules.find(idx);
        [[maybe_unused]] auto insertion_result = all_indices_new.emplace(
            dst_it != const_replrules.end() ? dst_it->second : idx);
      });
  std::swap(all_indices_new, all_indices);

  return mutated;
}

/// If using orthonormal representation, resolves Kronecker deltas (=overlaps
/// between indices in orthonormal spaces) in summations
/// @throw zero_result if @c expr is zero
template <Statistics S>
void reduce_wick_impl(std::shared_ptr<Product> &expr,
                      const container::set<Index> &external_indices,
                      const Context &ctx) {
  if (ctx.metric() == IndexSpaceMetric::Unit) {
    bool pass_mutated = false;
    do {
      pass_mutated = false;

      // extract current indices
      std::set<Index, Index::LabelCompare> all_indices;
      ranges::for_each(*expr, [&all_indices](const auto &factor) {
        if (factor->template is<Tensor>()) {
          ranges::for_each(factor->template as<const Tensor>().indices(),
                           [&all_indices](const Index &idx) {
                             [[maybe_unused]] auto result =
                                 all_indices.insert(idx);
                           });
        }
      });

      const auto replacement_rules = compute_index_replacement_rules<S>(
          expr, external_indices, all_indices, ctx.index_space_registry());

      if (Logger::instance().wick_reduce) {
        std::wcout << "reduce_wick_impl(expr, external_indices):\n  expr = "
                   << expr->to_latex() << "\n  external_indices = ";
        ranges::for_each(external_indices, [](auto &index) {
          std::wcout << index.label() << " ";
        });
        std::wcout << "\n  replrules = ";
        ranges::for_each(replacement_rules, [](auto &index) {
          std::wcout << to_latex(index.first) << "\\to"
                     << to_latex(index.second) << "\\,";
        });
        std::wcout.flush();
      }

      if (!replacement_rules.empty()) {
        auto isr = ctx.index_space_registry();
        pass_mutated = apply_index_replacement_rules(
            expr, replacement_rules, external_indices, all_indices, isr);
      }

      if (Logger::instance().wick_reduce) {
        std::wcout << "\n  result = " << expr->to_latex() << std::endl;
      }

    } while (pass_mutated);  // keep reducing until stop changing
  } else
    abort();  // programming error?
}

template <Statistics S>
struct NullNormalOperatorCanonicalizerDeregister {
  void operator()(void *) {
    const auto nop_labels = NormalOperator<S>::labels();
    TensorCanonicalizer::deregister_instance(nop_labels[0]);
    TensorCanonicalizer::deregister_instance(nop_labels[1]);
  }
};

}  // namespace detail

inline container::set<Index> extract_external_indices(const Expr &expr) {
  if (ranges::any_of(expr, [](auto &e) { return e.template is<Sum>(); }))
    throw std::invalid_argument(
        "extract_external_indices(expr): expr must be expanded (i.e. no "
        "subexpression can be a Sum)");

  container::map<Index, int64_t> idx_counter;
  auto visitor = [&idx_counter](const auto &expr) {
    auto expr_as_abstract_tensor =
        std::dynamic_pointer_cast<AbstractTensor>(expr);
    if (expr_as_abstract_tensor) {
      ranges::for_each(expr_as_abstract_tensor->_braket(),
                       [&idx_counter](const auto &v) {
                         auto it = idx_counter.find(v);
                         if (it == idx_counter.end()) {
                           idx_counter.emplace(v, 1);
                         } else {
                           it->second++;
                         }
                       });
    }
  };
  expr.visit(visitor);

  return idx_counter |
         ranges::views::filter([](const auto &v) { return v.second == 1; }) |
         ranges::views::transform([](const auto &v) { return v.first; }) |
         ranges::to<container::set<Index>>;
}

template <Statistics S>
ExprPtr WickTheorem<S>::compute(const bool count_only,
                                const bool skip_input_canonicalization) {
  // need to avoid recanonicalization of operators produced by WickTheorem
  // by rapid canonicalization to avoid undoing all the good
  // the NormalOperator<S>::normalize did ... use RAII
  // 1. detail::NullNormalOperatorCanonicalizerDeregister<S> will restore state
  // of tensor canonicalizer
  // 2. this is the RAII object whose destruction will restore state of
  // the tensor canonicalizer
  std::unique_ptr<void, detail::NullNormalOperatorCanonicalizerDeregister<S>>
      raii_null_nop_canonicalizer;
  // 3. this makes the RAII object  ... NOT reentrant, only to be called in
  // top-level WickTheorem after initial canonicalization
  auto disable_nop_canonicalization = [&raii_null_nop_canonicalizer]() {
    if (!raii_null_nop_canonicalizer) {
      const auto nop_labels = NormalOperator<S>::labels();
      assert(nop_labels.size() == 2);
      TensorCanonicalizer::try_register_instance(
          std::make_shared<NullTensorCanonicalizer>(), nop_labels[0]);
      TensorCanonicalizer::try_register_instance(
          std::make_shared<NullTensorCanonicalizer>(), nop_labels[1]);
      raii_null_nop_canonicalizer = decltype(raii_null_nop_canonicalizer)(
          (void *)&raii_null_nop_canonicalizer, {});
    }
  };

  // have an Expr as input? Apply recursively ...
  if (expr_input_) {
    /// expand, then apply recursively to products
    if (Logger::instance().wick_harness)
      std::wcout << "WickTheorem<S>::compute: input (before expand) = "
                 << to_latex_align(expr_input_) << std::endl;
    expand(expr_input_);
    if (Logger::instance().wick_harness)
      std::wcout << "WickTheorem<S>::compute: input (after expand) = "
                 << to_latex_align(expr_input_) << std::endl;
    // if sum, canonicalize and apply to each summand ...
    if (expr_input_->is<Sum>()) {
      if (!skip_input_canonicalization) {
        // initial full canonicalization
        canonicalize(expr_input_);
        assert(!expr_input_->as<Sum>().empty());
      }

      // NOW disable canonicalization of normal operators
      // N.B. even if skipped initial input canonicalization need to disable
      // subsequent nop canonicalization
      disable_nop_canonicalization();

      // parallelize over summands
      auto result = std::make_shared<Sum>();
      std::mutex result_mtx;  // serializes updates of result
      auto summands = expr_input_->as<Sum>().summands();

      // find external_indices if don't have them
      if (!external_indices_) {
        ranges::find_if(summands, [this](const auto &summand) {
          if (summand.template is<Sum>())  // summands must not be a Sum
            throw std::invalid_argument(
                "WickTheorem<S>::compute(expr): expr is a Sum with one of the "
                "summands also a Sum, WickTheorem can only accept a fully "
                "expanded Sum");
          else if (summand.template is<Product>()) {
            external_indices_ = extract_external_indices(
                *(summand.template as_shared_ptr<Product>()));
            return true;
          } else
            return false;
        });
      }

      if (Logger::instance().wick_harness)
        std::wcout << "WickTheorem<S>::compute: input (after canonicalize) has "
                   << summands.size() << " terms = " << to_latex_align(result)
                   << std::endl;

      auto wick_task = [&result, &result_mtx, this,
                        &count_only](const ExprPtr &input) {
        WickTheorem wt(input->clone(), *this);
        auto task_result = wt.compute(
            count_only, /* definitely skip input canonicalization */ true);
        stats() += wt.stats();
        if (task_result) {
          std::scoped_lock<std::mutex> lock(result_mtx);
          result->append(task_result);
        }
      };
      sequant::for_each(summands, wick_task);

      // if the sum is empty return zero
      // if the sum has 1 summand, return it directly
      ExprPtr result_expr = result;
      if (result->summands().size() == 0) {
        result_expr = ex<Constant>(0);
      }
      if (result->summands().size() == 1)
        result_expr = std::move(result->summands()[0]);

      return result_expr;
    }
    // ... else if a product, find NormalOperatorSequence, if any, and compute
    // ...
    else if (expr_input_->is<Product>()) {
      if (!skip_input_canonicalization) {  // canonicalize, unless told to skip
        auto canon_byproduct = expr_input_->rapid_canonicalize();
        assert(canon_byproduct ==
               nullptr);  // canonicalization of Product always returns nullptr
      }
      // NOW disable canonicalization of normal operators
      // N.B. even if skipped initial input canonicalization need to disable
      // subsequent nop canonicalization
      disable_nop_canonicalization();

      // find external_indices if don't have them
      if (!external_indices_) {
        external_indices_ =
            extract_external_indices(*(expr_input_.as_shared_ptr<Product>()));
      } else {
        assert(
            extract_external_indices(*(expr_input_.as_shared_ptr<Product>())) ==
            *external_indices_);
      }

      // split off NormalOperators into input_
      auto first_nop_it = ranges::find_if(
          *expr_input_,
          [](const ExprPtr &expr) { return expr->is<NormalOperator<S>>(); });
      // if have ops, split into nop sequence and cnumber "prefactor"
      if (first_nop_it != ranges::end(*expr_input_)) {
        // extract into prefactor and op sequence
        ExprPtr prefactor =
            ex<CProduct>(expr_input_->as<Product>().scalar(), ExprPtrList{});
        auto nopseq = std::make_shared<NormalOperatorSequence<S>>();
        for (const auto &factor : *expr_input_) {
          if (factor->template is<NormalOperator<S>>()) {
            nopseq->push_back(factor->template as<NormalOperator<S>>());
          } else {
            assert(factor->is_cnumber());
            *prefactor *= *factor;
          }
        }
        init_input(nopseq);

        // compute and record/analyze topological NormalOperator and Index
        // partitions
        if (use_topology_) {
          if (Logger::instance().wick_topology)
            std::wcout
                << "WickTheorem<S>::compute: input to topology computation = "
                << to_latex(expr_input_) << std::endl;

            // construct graph representation of the tensor product
#if USE_TENSOR_NETWORK_V2
          TensorNetworkV2 tn(expr_input_->as<Product>().factors());
          auto g = tn.create_graph();
          const auto &graph = g.bliss_graph;
          const auto &vlabels = g.vertex_labels;
          const auto &vcolors = g.vertex_colors;
          const auto &vtypes = g.vertex_types;
#else
          TensorNetwork tn(expr_input_->as<Product>().factors());
          auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
              tn.make_bliss_graph(
                  {/* need labels to find normal operators */ .make_labels =
                       true,
                   .make_texlabels = false});
#endif
          const auto n = vtypes.size();
          assert(vcolors.size() == n);
          assert(vlabels.size() == n);
          const auto &tn_edges = tn.edges();
          const auto &tn_tensors = tn.tensors();

          if (Logger::instance().wick_topology) {
            std::basic_ostringstream<wchar_t> oss;
#if USE_TENSOR_NETWORK_V2
            graph->write_dot(oss, vlabels);
#else
            graph->write_dot(oss, vlabels);
#endif
            std::wcout
                << "WickTheorem<S>::compute: colored graph produced from TN = "
                << std::endl
                << oss.str() << std::endl;
          }

          // identify vertex indices of NormalOperator objects and Indices
          // 1. list of vertex indices corresponding to NormalOperator objects
          //    on the TN graph and their ordinals in NormalOperatorSequence
          //    N.B. for NormalOperators the vertex indices coincide with
          //    the ordinals
          container::map<size_t, size_t> nop_vidx_ord;
          // 2. list of vertex indices corresponding to Index objects on the TN
          //    graph that appear in NormalOperatorsSequence and
          //    their ordinals therein
          //    N.B. for Index objects the vertex indices do NOT coincide with
          //         the ordinals
          container::map<size_t, size_t> index_vidx_ord;
          {
            const auto &nop_labels = NormalOperator<S>::labels();
            const auto nop_labels_begin = begin(nop_labels);
            const auto nop_labels_end = end(nop_labels);

            using opseq_view_type =
                flattened_rangenest<NormalOperatorSequence<S>>;
            auto opseq_view = opseq_view_type(input_.get());
            const auto opseq_view_begin = ranges::begin(opseq_view);
            const auto opseq_view_end = ranges::end(opseq_view);

            // NormalOperators are not reordered by canonicalization, hence the
            // ordinal can be computed by counting
            std::size_t nop_ord = 0;
            for (size_t v = 0; v != n; ++v) {
              if (vtypes[v] == VertexType::TensorCore &&
                  (std::find(nop_labels_begin, nop_labels_end, vlabels[v]) !=
                   nop_labels_end)) {
                auto insertion_result = nop_vidx_ord.emplace(v, nop_ord++);
                assert(insertion_result.second);
              }
              if (vtypes[v] == VertexType::Index && !input_->empty()) {
                auto &idx = (tn_edges.begin() + v)->idx();
                auto idx_it_in_opseq = ranges::find_if(
                    opseq_view,
                    [&idx](const auto &v) { return v.index() == idx; });
                if (idx_it_in_opseq != opseq_view_end) {
                  const auto ord =
                      ranges::distance(opseq_view_begin, idx_it_in_opseq);
                  auto insertion_result = index_vidx_ord.emplace(v, ord);
                  assert(insertion_result.second);
                }
              }
            }
          }

          // compute and save graph automorphism generators
          std::vector<std::vector<unsigned int>> aut_generators;
          {
            bliss::Stats stats;
            graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

            auto save_aut = [&aut_generators](const unsigned int n,
                                              const unsigned int *aut) {
              aut_generators.emplace_back(aut, aut + n);
            };

            graph->find_automorphisms(
                stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);

            if (Logger::instance().wick_topology) {
              std::basic_ostringstream<wchar_t> oss2;
              bliss::print_auts(aut_generators, oss2, vlabels);
              std::wcout << "WickTheorem<S>::compute: colored graph "
                            "automorphism generators = \n"
                         << oss2.str() << std::endl;
            }
          }

          // Use automorphisms to determine groups of topologically equivalent
          // NormalOperator and Op objects.
          // @param vertices maps vertex indices of the objects to their
          //        ordinals in the sequence of such objects within
          //        the NormalOperatorSequence
          // @param nontrivial_partitions_only if true, only partitions with
          // more than one element, are reported, else even trivial
          // partitions with a single partition will be reported
          // @param vertex_pair_exclude a callable that accepts 2 vertex
          // indices and returns true if the automorphism of this pair
          // of indices is to be ignored
          // @return the \c {vertex_to_partition_idx,npartitions} pair in
          // which \c vertex_to_partition_idx maps vertex indices that are
          // part of nontrivial partitions to their (1-based) partition indices
          auto compute_partitions = [&aut_generators](
                                        const container::map<size_t, size_t>
                                            &vertices,
                                        bool nontrivial_partitions_only,
                                        auto &&vertex_pair_exclude) {
            container::map<size_t, size_t> vertex_to_partition_idx;
            int next_partition_idx = -1;

            // using each automorphism generator
            for (auto &&aut : aut_generators) {
              // skip automorphism generators that involve vertices that are
              // not part of list `vertices`
              // this prevents topology exploitation for spin-free Wick
              // TODO learn how to compute partitions correctly for
              //      spin-free cases
              const auto nv = aut.size();
              bool aut_contains_other_vertices = false;
              for (std::size_t v = 0; v != nv; ++v) {
                const auto v_is_in_aut = v != aut[v];
                if (v_is_in_aut && !vertices.contains(v)) {
                  aut_contains_other_vertices = true;
                  break;
                }
              }
              if (aut_contains_other_vertices) continue;

              // update partitions
              for (auto &&[v1, ord1] : vertices) {
                const auto v2 = aut[v1];
                if (v2 != v1 &&
                    !vertex_pair_exclude(
                        v1, v2)) {  // if the automorphism maps this vertex to
                                    // another ... they both must be in the same
                                    // partition
                  assert(vertices.find(v2) != vertices.end());
                  auto v1_partition_it = vertex_to_partition_idx.find(v1);
                  auto v2_partition_it = vertex_to_partition_idx.find(v2);
                  const bool v1_has_partition =
                      v1_partition_it != vertex_to_partition_idx.end();
                  const bool v2_has_partition =
                      v2_partition_it != vertex_to_partition_idx.end();
                  if (v1_has_partition &&
                      v2_has_partition) {  // both are in partitions? make sure
                                           // they are in the same partition.
                                           // N.B. this may leave gaps in
                                           // partition indices ... no biggie
                    const auto v1_part_idx = v1_partition_it->second;
                    const auto v2_part_idx = v2_partition_it->second;
                    if (v1_part_idx !=
                        v2_part_idx) {  // if they have different partition
                                        // indices, change the larger of the two
                                        // indices to match the lower
                      const auto target_part_idx =
                          std::min(v1_part_idx, v2_part_idx);
                      for (auto &v : vertex_to_partition_idx) {
                        if (v.second == v1_part_idx || v.second == v2_part_idx)
                          v.second = target_part_idx;
                      }
                    }
                  } else if (v1_has_partition) {  // only v1 is in a partition?
                                                  // place v2 in it
                    const auto v1_part_idx = v1_partition_it->second;
                    vertex_to_partition_idx.emplace(v2, v1_part_idx);
                  } else if (v2_has_partition) {  // only v2 is in a partition?
                                                  // place v1 in it
                    const auto v2_part_idx = v2_partition_it->second;
                    vertex_to_partition_idx.emplace(v1, v2_part_idx);
                  } else {  // neither is in a partition? place both in the next
                            // available partition
                    const size_t target_part_idx = ++next_partition_idx;
                    vertex_to_partition_idx.emplace(v1, target_part_idx);
                    vertex_to_partition_idx.emplace(v2, target_part_idx);
                  }
                }
              }
            }
            if (!nontrivial_partitions_only) {
              ranges::for_each(vertices, [&](const auto &vidx_ord) {
                auto &&[vidx, ord] = vidx_ord;
                if (vertex_to_partition_idx.find(vidx) ==
                    vertex_to_partition_idx.end()) {
                  vertex_to_partition_idx.emplace(vidx, ++next_partition_idx);
                }
              });
            }
            const auto npartitions = next_partition_idx;
            return std::make_tuple(vertex_to_partition_idx, npartitions);
          };

          // compute NormalOperator->partition map, convert to partition lists
          // (if any), and register via set_nop_partitions to be used in full
          // contractions
          auto do_not_skip_elements = [](size_t v1, size_t v2) {
            return false;
          };
          auto [nop_vidx2pidx, nop_npartitions] = compute_partitions(
              nop_vidx_ord, /* nontrivial_partitions_only = */ true,
              do_not_skip_elements);

          // converts vertex ordinal to partition key map into a sequence of
          // partitions, each composed of the corresponding ordinals of the
          // vertices in the vertex_list sequence
          // @param vidx2pidx a map from vertex index (in TN) to its
          //        (1-based) partition index
          // @param npartitions the total number of partitions
          // @param vidx_ord ordered sequence of vertex indices, object
          // with vertex index `vidx` will be mapped to ordinal
          // `vidx_ord[vidx]`
          // @return sequence of partitions, sorted by the smallest ordinal
          auto extract_partitions = [](const auto &vidx2pidx,
                                       const auto npartitions,
                                       const auto &vidx_ord) {
            container::svector<container::svector<size_t>> partitions;

            assert(npartitions > -1);
            const size_t max_pidx = npartitions;
            partitions.reserve(max_pidx);

            // iterate over all partition indices ... note that there may be
            // gaps so count the actual partitions
            size_t partition_cnt = 0;
            for (size_t p = 0; p <= max_pidx; ++p) {
              bool p_found = false;
              for (const auto &[vidx, pidx] : vidx2pidx) {
                if (pidx == p) {
                  // !!remember to map the vertex index into the operator
                  // index!!
                  assert(vidx_ord.find(vidx) != vidx_ord.end());
                  const auto ordinal = vidx_ord.find(vidx)->second;
                  if (p_found == false) {  // first time this is found
                    partitions.emplace_back(container::svector<size_t>{
                        static_cast<size_t>(ordinal)});
                  } else
                    partitions[partition_cnt].emplace_back(ordinal);
                  p_found = true;
                }
              }
              if (p_found) ++partition_cnt;
            }

            // sort each partition
            for (auto &partition : partitions) {
              ranges::sort(partition);
            }

            // sort partitions in the order of increasing first element
            ranges::sort(partitions, [](const auto &p1, const auto &p2) {
              return p1.front() < p2.front();
            });

            return partitions;
          };

          if (!nop_vidx2pidx.empty()) {
            container::svector<container::svector<size_t>> nop_partitions;

            nop_partitions = extract_partitions(nop_vidx2pidx, nop_npartitions,
                                                nop_vidx_ord);

            if (Logger::instance().wick_topology) {
              std::wcout
                  << "WickTheorem<S>::compute: topological nop partitions:{\n";
              ranges::for_each(nop_partitions, [](auto &&part) {
                std::wcout << "{";
                ranges::for_each(part,
                                 [](auto &&p) { std::wcout << p << " "; });
                std::wcout << "}";
              });
              std::wcout << "}" << std::endl;
            }

            this->set_nop_partitions(nop_partitions);
          }

          // compute Index->partition map, and convert to partition lists (if
          // any), and check that use_topology_ is compatible with index
          // partitions
          // Index partitions are constructed to *only* include Index
          // objects attached to the bra/ket of any NormalOperator! hence
          // need to use filter in computing partitions
          auto exclude_index_vertex_pair = [&tn_tensors, &tn_edges](size_t v1,
                                                                    size_t v2) {
            // v1 and v2 are vertex indices and also index the edges in the
            // WickGraph
            assert(v1 < tn_edges.size());
            assert(v2 < tn_edges.size());
            const auto &edge1 = *(tn_edges.begin() + v1);
            const auto &edge2 = *(tn_edges.begin() + v2);
            auto connected_to_same_nop =
                [&tn_tensors](const auto &edge1, const auto &edge2) -> bool {
              const auto nt1 =
#if USE_TENSOR_NETWORK_V2
                  edge1.vertex_count();
#else
                  edge1.size();
#endif
              assert(nt1 <= 2);
              const auto nt2 =
#if USE_TENSOR_NETWORK_V2
                  edge2.vertex_count();
#else
                  edge2.size();
#endif
              assert(nt2 <= 2);
              for (auto i1 = 0; i1 != nt1; ++i1) {
                const auto tensor1_ord =
#if USE_TENSOR_NETWORK_V2
                    i1 == 0 ? edge1.first_vertex().getTerminalIndex()
                            : edge1.second_vertex().getTerminalIndex();
#else
                    edge1[i1].tensor_ord;
#endif
                for (auto i2 = 0; i2 != nt2; ++i2) {
                  const auto tensor2_ord =
#if USE_TENSOR_NETWORK_V2
                      i2 == 0 ? edge2.first_vertex().getTerminalIndex()
                              : edge2.second_vertex().getTerminalIndex();
#else
                      edge2[i2].tensor_ord;
#endif
                  if (tensor1_ord == tensor2_ord) {
                    auto tensor_ord = tensor1_ord;
                    const std::shared_ptr<AbstractTensor> &tensor_ptr =
                        tn_tensors.at(tensor_ord);
                    if (std::dynamic_pointer_cast<NormalOperator<S>>(
                            tensor_ptr))
                      return true;
                  }
                }
              }
              return false;
            };
            const bool exclude = !connected_to_same_nop(edge1, edge2);
            return exclude;
          };

          // index_vidx2pidx maps vertex index (see
          // index_vidx_ord) to partition index
          container::map<size_t, size_t> index_vidx2pidx;
          int index_npartitions = -1;
          std::tie(index_vidx2pidx, index_npartitions) = compute_partitions(
              index_vidx_ord, /* nontrivial_partitions_only = */ false,
              /* this is to ensure that each index partition only involves
                 indices attached to bra or to ket of same nop */
              exclude_index_vertex_pair);

          if (!index_vidx2pidx.empty()) {
            container::svector<container::svector<size_t>> index_partitions;

            index_partitions = extract_partitions(
                index_vidx2pidx, index_npartitions, index_vidx_ord);

            if (Logger::instance().wick_topology) {
              std::wcout << "WickTheorem<S>::compute: topological index "
                            "partitions:{\n";
              ranges::for_each(index_vidx2pidx, [&tn_edges](auto &&vidx_pidx) {
                auto &&[vidx, pidx] = vidx_pidx;
                assert(vidx < tn_edges.size());
                auto &idx = (tn_edges.begin() + vidx)->idx();
                std::wcout << "Index " << idx.full_label() << " -> partition "
                           << pidx << "\n";
              });
              std::wcout << "}" << std::endl;
            }

            this->set_op_partitions(index_partitions);
          }
        }

        if (!input_->empty()) {
          if (Logger::instance().wick_contract) {
            std::wcout
                << "WickTheorem<S>::compute: input to compute_nopseq = {\n";
            for (auto &&nop : input_) std::wcout << to_latex(nop) << "\n";
            std::wcout << "}" << std::endl;
          }
          auto result = compute_nopseq(count_only);
          if (result) {  // simplify if obtained nonzero ...
            result = prefactor * result;
            expand(result);
            this->reduce(result);
            rapid_simplify(result);
            canonicalize(result);
            rapid_simplify(
                result);  // rapid_simplify again since canonization may produce
                          // new opportunities (e.g. terms cancel, etc.)
          } else
            result = ex<Constant>(0);
          return result;
        }
      } else {  // product does not include ops
        return expr_input_;
      }
    }  // expr_input_->is<Product>()
    // ... else if NormalOperatorSequence already, compute ...
    else if (expr_input_->is<NormalOperatorSequence<S>>()) {
      abort();  // expr_input_ should no longer be nonnull if constructed with
                // an expression that's a NormalOperatorSequence<S>
      init_input(
          expr_input_.template as_shared_ptr<NormalOperatorSequence<S>>());
      // NB no simplification possible for a bare product w/ full contractions
      // ... partial contractions will need simplification
      return compute_nopseq(count_only);
    } else  // ... else do nothing
      return expr_input_;
  } else  // given a NormalOperatorSequence instead of an expression
    return compute_nopseq(count_only);
  abort();
}

template <Statistics S>
void WickTheorem<S>::reduce(ExprPtr &expr) const {
  if (Logger::instance().wick_reduce) {
    std::wcout << "WickTheorem<S>::reduce: input = "
               << to_latex_align(expr, 20, 1) << std::endl;
  }
  // there are 2 possibilities: expr is a single Product, or it's a Sum of
  // Products
  if (expr->type_id() == Expr::get_type_id<Product>()) {
    auto expr_cast = std::static_pointer_cast<Product>(expr);
    try {
      assert(external_indices_);
      detail::reduce_wick_impl<S>(expr_cast, *external_indices_,
                                  get_default_context(S));
      expr = expr_cast;
    } catch (detail::zero_result &) {
      expr = std::make_shared<Constant>(0);
    }
  } else {
    assert(expr->type_id() == Expr::get_type_id<Sum>());
    for (auto &&subexpr : *expr) {
      assert(subexpr->is<Product>());
      auto subexpr_cast = std::static_pointer_cast<Product>(subexpr);
      try {
        assert(external_indices_);
        detail::reduce_wick_impl<S>(subexpr_cast, *external_indices_,
                                    get_default_context(S));
        subexpr = subexpr_cast;
      } catch (detail::zero_result &) {
        subexpr = std::make_shared<Constant>(0);
      }
    }
  }

  if (Logger::instance().wick_reduce) {
    std::wcout << "WickTheorem<S>::reduce: result = "
               << to_latex_align(expr, 20, 1) << std::endl;
  }
}
template <Statistics S>
WickTheorem<S>::~WickTheorem() {}

}  // namespace sequant

#endif  // SEQUANT_WICK_IMPL_HPP
