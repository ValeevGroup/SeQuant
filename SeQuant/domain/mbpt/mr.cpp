//
// Created by Eduard Valeyev on 2019-02-19.
//

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/domain/mbpt/mr.hpp>

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/wick.hpp>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/reverse_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

namespace sequant {
namespace mbpt {
namespace mr {

qninterval_t ncre(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied || s == IndexSpace::active ||
         s == IndexSpace::active_unoccupied);
  if (s == IndexSpace::active_occupied)
    return qns[0];
  else if (s == IndexSpace::active)
    return qns[2];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[4];
}

qninterval_t ncre(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  if (s == IndexSpace::active_occupied)
    return qns[0];
  else if (s == IndexSpace::active)
    return qns[2];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[4];
}

qninterval_t ncre_occ(qns_t qns) {
  return ncre(qns, IndexSpace::active_occupied);
}

qninterval_t ncre_act(qns_t qns) { return ncre(qns, IndexSpace::active); }

qninterval_t ncre_uocc(qns_t qns) {
  return ncre(qns, IndexSpace::active_unoccupied);
}

qninterval_t ncre(qns_t qns) { return qns[0] + qns[2] + qns[4]; }

qninterval_t nann(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied || s == IndexSpace::active ||
         s == IndexSpace::active_unoccupied);
  if (s == IndexSpace::active_occupied)
    return qns[1];
  else if (s == IndexSpace::active)
    return qns[3];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[5];
}

qninterval_t nann(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  if (s == IndexSpace::active_occupied)
    return qns[1];
  else if (s == IndexSpace::active)
    return qns[3];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[5];
}

qninterval_t nann_occ(qns_t qns) {
  return nann(qns, IndexSpace::active_occupied);
}

qninterval_t nann_act(qns_t qns) { return nann(qns, IndexSpace::active); }

qninterval_t nann_uocc(qns_t qns) {
  return nann(qns, IndexSpace::active_unoccupied);
}

qninterval_t nann(qns_t qns) { return qns[1] + qns[3] + qns[5]; }

qns_t combine(qns_t a, qns_t b) {
  // particle contractions (i.e. above Fermi level; N.B. active are above
  // closed-shell Fermi level)
  const auto ncontr_uocc =
      qninterval_t{0, std::min(ncre(b, IndexSpace::active_unoccupied).upper(),
                               nann(a, IndexSpace::active_unoccupied).upper())};
  const auto ncontr_act =
      qninterval_t{0, std::min(ncre(b, IndexSpace::active).upper(),
                               nann(a, IndexSpace::active).upper())};
  // hole contractions (i.e. below Fermi level)
  const auto ncontr_occ =
      qninterval_t{0, std::min(nann(b, IndexSpace::active_occupied).upper(),
                               ncre(a, IndexSpace::active_occupied).upper())};
  const auto nc_occ =
      nonnegative(ncre(a, IndexSpace::active_occupied) +
                  ncre(b, IndexSpace::active_occupied) - ncontr_occ);
  const auto nc_act = nonnegative(ncre(a, IndexSpace::active) +
                                  ncre(b, IndexSpace::active) - ncontr_act);
  const auto nc_uocc =
      nonnegative(ncre(a, IndexSpace::active_unoccupied) +
                  ncre(b, IndexSpace::active_unoccupied) - ncontr_uocc);
  const auto na_occ =
      nonnegative(nann(a, IndexSpace::active_occupied) +
                  nann(b, IndexSpace::active_occupied) - ncontr_occ);
  const auto na_act = nonnegative(nann(a, IndexSpace::active) +
                                  nann(b, IndexSpace::active) - ncontr_act);
  const auto na_uocc =
      nonnegative(nann(a, IndexSpace::active_unoccupied) +
                  nann(b, IndexSpace::active_unoccupied) - ncontr_uocc);
  return qns_t{nc_occ, na_occ, nc_act, na_act, nc_uocc, na_uocc};
}

}  // namespace mr
}  // namespace mbpt

mbpt::mr::qns_t adjoint(mbpt::mr::qns_t qns) {
  return mbpt::mr::qns_t{nann(qns, IndexSpace::active_occupied),
                         ncre(qns, IndexSpace::active_occupied),
                         nann(qns, IndexSpace::active),
                         ncre(qns, IndexSpace::active),
                         nann(qns, IndexSpace::active_unoccupied),
                         ncre(qns, IndexSpace::active_unoccupied)};
}

namespace mbpt {
namespace mr {

OpMaker::OpMaker(OpType op, std::size_t nbra, std::size_t nket)
    : base_type(op) {
  nket = nket == std::numeric_limits<std::size_t>::max() ? nbra : nket;
  assert(nbra > 0 || nket > 0);

  const auto unocc = IndexSpace::active_maybe_unoccupied;
  const auto occ = IndexSpace::active_maybe_occupied;
  switch (to_class(op)) {
    case OpClass::ex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, unocc);
      ket_spaces_ = decltype(ket_spaces_)(nket, occ);
      break;
    case OpClass::deex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, occ);
      ket_spaces_ = decltype(ket_spaces_)(nket, unocc);
      break;
    case OpClass::gen:
      bra_spaces_ = decltype(bra_spaces_)(nbra, IndexSpace::complete);
      ket_spaces_ = decltype(ket_spaces_)(nket, IndexSpace::complete);
      break;
  }
}

OpMaker::OpMaker(OpType op,
                 const container::svector<IndexSpace::Type>& cre_spaces,
                 const container::svector<IndexSpace::Type>& ann_spaces)
    : base_type(op) {
  bra_spaces_ = cre_spaces;
  ket_spaces_ = ann_spaces;
}

#include <SeQuant/domain/mbpt/mr/op.impl.cpp>

ExprPtr T_act_(std::size_t K) {
  return OpMaker(OpType::t,
                 container::svector<IndexSpace::Type>(K, IndexSpace::active),
                 container::svector<IndexSpace::Type>(K, IndexSpace::active))();
}

ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);
  switch (k) {
    case 1:
      switch (get_default_context().vacuum()) {
        case Vacuum::Physical:
          return OpMaker(OpType::h, 1)();
        case Vacuum::SingleProduct:
          return OpMaker(OpType::f̃, 1)();
        case Vacuum::MultiProduct:
          return OpMaker(OpType::f, 1)();
        default:
          abort();
      }

    case 2:
      return OpMaker(OpType::g, 2)();

    default:
      abort();
  }
}

ExprPtr H(std::size_t k) {
  assert(k > 0 && k <= 2);
  return k == 1 ? H_(1) : H_(1) + H_(2);
}

ExprPtr F() {
  // add \bar{g}^{\kappa x}_{\lambda y} \gamma^y_x with x,y in occ_space_type
  auto make_g_contribution = [](const auto occ_space_type) {
    return mbpt::OpMaker<Statistics::FermiDirac>::make(
        {IndexSpace::complete}, {IndexSpace::complete},
        [=](auto braidxs, auto ketidxs, Symmetry opsymm) {
          auto m1 = Index::make_tmp_index(
              IndexSpace{occ_space_type, IndexSpace::nullqns});
          auto m2 = Index::make_tmp_index(
              IndexSpace{occ_space_type, IndexSpace::nullqns});
          assert(opsymm == Symmetry::antisymm || opsymm == Symmetry::nonsymm);
          if (opsymm == Symmetry::antisymm) {
            braidxs.push_back(m1);
            ketidxs.push_back(m2);
            return ex<Tensor>(to_wstring(mbpt::OpType::g), braidxs, ketidxs,
                              std::vector<Index>{}, Symmetry::antisymm) *
                   ex<Tensor>(to_wstring(mbpt::OpType::RDM), IndexList{m2},
                              IndexList{m1}, IndexList{}, Symmetry::nonsymm);
          } else {  // opsymm == Symmetry::nonsymm
            auto braidx_J = braidxs;
            braidx_J.push_back(m1);
            auto ketidxs_J = ketidxs;
            ketidxs_J.push_back(m2);
            auto braidx_K = braidxs;
            braidx_K.push_back(m1);
            auto ketidxs_K = ketidxs;
            ketidxs_K.emplace(begin(ketidxs_K), m2);
            return (ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_J, ketidxs_J,
                               std::vector<Index>{}, Symmetry::nonsymm) -
                    ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_K, ketidxs_K,
                               std::vector<Index>{}, Symmetry::nonsymm)) *
                   ex<Tensor>(to_wstring(mbpt::OpType::RDM), IndexList{m2},
                              IndexList{m1}, IndexList{}, Symmetry::nonsymm);
          }
        });
  };

  switch (get_default_context().vacuum()) {
    case Vacuum::Physical:
      return OpMaker(OpType::h, 1)() +
             make_g_contribution(IndexSpace::maybe_occupied);  // all occupieds

    case Vacuum::SingleProduct:
      return OpMaker(OpType::f̃, 1)() +
             make_g_contribution(IndexSpace::active);  // actives only

    case Vacuum::MultiProduct:
      return OpMaker(OpType::f, 1)();
    default:
      abort();
  }
}

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  const auto spinorbital =
      get_default_context().spbasis() == SPBasis::spinorbital;
  // convention is to use different label for spin-orbital and spin-free RDM
  const auto rdm_label = spinorbital ? optype2label.at(OpType::RDM) : L"Γ";

  FWickTheorem wick{expr};
  wick.use_topology(use_top)
      .set_nop_connections(nop_connections)
      .full_contractions(false);
  auto result = wick.compute();
  simplify(result);

  // if obtained nontrivial result ...
  if (!result.template is<Constant>()) {
    // need pre-postprocessing unless used extended Wick theorem
    if (get_default_context().vacuum() != Vacuum::MultiProduct) {
      assert(get_default_context().vacuum() != Vacuum::Invalid);

      // replace NormalOperator with RDM in target RDM space:
      // - if Vacuum::Physical: IndexSpace::maybe_occupied
      // - if Vacuum::SingleProduct: IndexSpace::active
      const auto target_rdm_space_type =
          get_default_context().vacuum() == Vacuum::SingleProduct
              ? IndexSpace::active
              : IndexSpace::maybe_occupied;

      // do this in 2 steps (TODO factor out these components?)
      // 1. replace NOPs by RDM
      // 2. project RDM indices onto the target RDM subspace

      // STEP1. replace NOPs by RDM
      auto replace_nop_with_rdm = [&rdm_label, spinorbital](ExprPtr& exptr) {
        auto replace = [&rdm_label, spinorbital](const auto& nop) -> ExprPtr {
          using index_container = container::svector<Index>;
          auto braidxs = nop.annihilators() |
                         ranges::views::transform(
                             [](const auto& op) { return op.index(); }) |
                         ranges::to<index_container>();
          auto ketidxs = nop.creators() |
                         ranges::views::transform(
                             [](const auto& op) { return op.index(); }) |
                         ranges::to<index_container>();
          assert(braidxs.size() ==
                 ketidxs.size());  // need to handle particle # violating case?
          const auto rank = braidxs.size();
          return ex<Tensor>(
              rdm_label, braidxs, ketidxs, index_container{},
              rank > 1 && spinorbital ? Symmetry::antisymm : Symmetry::nonsymm);
        };

        if (exptr.template is<FNOperator>()) {
          exptr = replace(exptr.template as<FNOperator>());
        } else if (exptr.template is<BNOperator>()) {
          exptr = replace(exptr.template as<BNOperator>());
        }
      };
      result->visit(replace_nop_with_rdm, /* atoms_only = */ true);

      // STEP 2: project RDM indices onto the target RDM subspace
      // since RDM indices only make sense within a single TN expand + flatten
      // first, then do the projection individually for each TN
      expand(result);
      // flatten(result);  // TODO where is flatten?
      auto project_rdm_indices_to_target = [&](ExprPtr& exptr) {
        auto impl_for_single_tn = [&](ProductPtr& product_ptr) {
          // enlist all indices and count their instances
          auto for_each_index_in_tn = [](const auto& product_ptr,
                                         const auto& op) {
            ranges::for_each(product_ptr->factors(), [&](auto& factor) {
              auto tensor_ptr =
                  std::dynamic_pointer_cast<AbstractTensor>(factor);
              if (tensor_ptr) {
                ranges::for_each(tensor_ptr->_indices(),
                                 [&](auto& idx) { op(idx, *tensor_ptr); });
              }
            });
          };

          // compute external indices
          container::map<Index, std::size_t> indices_w_counts;
          auto retrieve_indices_with_counts =
              [&indices_w_counts](const auto& idx, auto& /* unused */) {
                auto found_it = indices_w_counts.find(idx);
                if (found_it != indices_w_counts.end()) {
                  found_it->second++;
                } else {
                  indices_w_counts.emplace(idx, 1);
                }
              };
          for_each_index_in_tn(product_ptr, retrieve_indices_with_counts);

          container::set<Index> external_indices =
              indices_w_counts | ranges::views::filter([](auto& idx_cnt) {
                auto& [idx, cnt] = idx_cnt;
                return cnt == 1;
              }) |
              ranges::views::keys | ranges::to<container::set<Index>>;

          // extract RDM-only and all indices
          container::set<Index> rdm_indices;
          std::set<Index, Index::LabelCompare> all_indices;
          auto retrieve_rdm_and_all_indices = [&rdm_indices, &all_indices,
                                               &rdm_label](const auto& idx,
                                                           const auto& tensor) {
            all_indices.insert(idx);
            if (tensor._label() == rdm_label) {
              rdm_indices.insert(idx);
            }
          };
          for_each_index_in_tn(product_ptr, retrieve_rdm_and_all_indices);

          // compute RDM->target replacement rules
          container::map<Index, Index> replacement_rules;
          ranges::for_each(rdm_indices, [&](const Index& idx) {
            const auto target_type =
                idx.space().type().intersection(target_rdm_space_type);
            if (target_type != IndexSpace::nulltype) {
              Index target = Index::make_tmp_index(
                  IndexSpace(target_type, idx.space().qns()));
              replacement_rules.emplace(idx, target);
            }
          });

          if (false) {
            std::wcout << "expr = " << product_ptr->to_latex()
                       << "\n  external_indices = ";
            ranges::for_each(external_indices, [](auto& index) {
              std::wcout << index.label() << " ";
            });
            std::wcout << "\n  replrules = ";
            ranges::for_each(replacement_rules, [](auto& index) {
              std::wcout << to_latex(index.first) << "\\to"
                         << to_latex(index.second) << "\\,";
            });
            std::wcout.flush();
          }

          if (!replacement_rules.empty()) {
            sequant::detail::apply_index_replacement_rules(
                product_ptr, replacement_rules, external_indices, all_indices);
          }
        };

        if (exptr.template is<Product>()) {
          auto product_ptr = exptr.template as_shared_ptr<Product>();
          impl_for_single_tn(product_ptr);
          exptr = product_ptr;
        } else {
          assert(exptr.template is<Sum>());
          auto result = std::make_shared<Sum>();
          for (auto& summand : exptr.template as<Sum>().summands()) {
            assert(summand.template is<Product>());
            auto result_summand = summand.template as<Product>().clone();
            auto product_ptr = result_summand.template as_shared_ptr<Product>();
            impl_for_single_tn(product_ptr);
            result->append(product_ptr);
          }
          exptr = result;
        }
      };
      project_rdm_indices_to_target(result);

      // rename dummy indices that might have been generated by
      // project_rdm_indices_to_target
      // + may combine terms

      // TensorCanonicalizer is given a custom comparer that moves active
      // indices to the front external-vs-internal trait still takes precedence
      auto current_index_comparer =
          TensorCanonicalizer::instance()->index_comparer();
      TensorCanonicalizer::instance()->index_comparer(
          [&](const Index& idx1, const Index& idx2) -> bool {
            const auto idx1_active = idx1.space().type() == IndexSpace::active;
            const auto idx2_active = idx2.space().type() == IndexSpace::active;
            if (idx1_active) {
              if (idx2_active)
                return current_index_comparer(idx1, idx2);
              else
                return true;
            } else {
              if (idx2_active)
                return false;
              else
                return current_index_comparer(idx1, idx2);
            }
          });
      simplify(result);
      TensorCanonicalizer::instance()->index_comparer(
          std::move(current_index_comparer));
    }
  }

  if (Logger::get_instance().wick_stats) {
    std::wcout << "WickTheorem stats: # of contractions attempted = "
               << wick.stats().num_attempted_contractions
               << " # of useful contractions = "
               << wick.stats().num_useful_contractions << std::endl;
  }
  return result;
}

namespace op {

ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);

  switch (k) {
    case 1:
      return ex<op_t>(
          [vacuum = get_default_context().vacuum()]() -> std::wstring_view {
            switch (vacuum) {
              case Vacuum::Physical:
                return optype2label.at(OpType::h);
              case Vacuum::SingleProduct:
                return optype2label.at(OpType::f̃);
              case Vacuum::MultiProduct:
                return optype2label.at(OpType::f);
              default:
                abort();
            }
          },
          [=]() -> ExprPtr { return mbpt::mr::H_(1); },
          [=](qnc_t& qns) {
            qns = combine(qnc_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}},
                          qns);
          });

    case 2:
      return ex<op_t>(
          []() -> std::wstring_view { return optype2label.at(OpType::g); },
          [=]() -> ExprPtr { return mbpt::mr::H_(2); },
          [=](qnc_t& qns) {
            qns = combine(qnc_t{{0, 2}, {0, 2}, {0, 2}, {0, 2}, {0, 2}, {0, 2}},
                          qns);
          });

    default:
      abort();
  }
}

ExprPtr H(std::size_t k) {
  assert(k > 0 && k <= 2);
  return k == 1 ? H_(1) : H_(1) + H_(2);
}

ExprPtr T_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::t); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return mr::T_(K);
      },
      [=](qnc_t& qns) {
        qns = combine(
            qnc_t{
                {0ul, 0ul}, {0ul, K}, {0ul, K}, {0ul, K}, {0ul, K}, {0ul, 0ul}},
            qns);
      });
}

ExprPtr T_act_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::t); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return mr::T_act_(K);
      },
      [=](qnc_t& qns) {
        qns = combine(
            qnc_t{
                {0ul, 0ul}, {0ul, 0ul}, {K, K}, {K, K}, {0ul, 0ul}, {0ul, 0ul}},
            qns);
      });
}

ExprPtr T(std::size_t K) {
  assert(K > 0);

  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + T_(k) : T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::λ); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return mr::Λ_(K);
      },
      [=](qnc_t& qns) {
        qns = combine(
            qnc_t{
                {0ul, K}, {0ul, 0ul}, {0ul, K}, {0ul, K}, {0ul, 0ul}, {0ul, K}},
            qns);
      });
}

ExprPtr Λ(std::size_t K) {
  assert(K > 0);

  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + Λ_(k) : Λ_(k);
  }
  return result;
}

ExprPtr A(std::int64_t K) {
  assert(K != 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::A); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return mr::A(K, K);
      },
      [=](qnc_t& qns) {
        const std::size_t abs_K = std::abs(K);
        if (K < 0)
          qns = combine(qnc_t{{0ul, abs_K},
                              {0ul, 0ul},
                              {0ul, abs_K},
                              {0ul, abs_K},
                              {0ul, 0ul},
                              {0ul, abs_K}},
                        qns);
        else
          qns = combine(qnc_t{{0ul, 0ul},
                              {0ul, abs_K},
                              {0ul, abs_K},
                              {0ul, abs_K},
                              {0ul, abs_K},
                              {0ul, 0ul}},
                        qns);
      });
}

// ExprPtr S(std::int64_t K) {
//   assert(K != 0);
//   return ex<op_t>([]() -> std::wstring_view { return L"S"; },
//                   [=]() -> ExprPtr {
//                     using namespace sequant::mbpt::sr;
//                     return mr::S(K, K);
//                   },
//                   [=](qnc_t& qns) {
//                     const std::size_t abs_K = std::abs(K);
//                     if (K < 0)
//                       qns = combine(qnc_t{{0ul,abs_K}, {0ul,0ul},
//                       {0ul,abs_K}, {0ul,abs_K}, {0ul,0ul}, {0ul,abs_K}},
//                       qns);
//                     else
//                       qns = combine(qnc_t{{0ul,0ul}, {0ul,abs_K},
//                       {0ul,abs_K}, {0ul,abs_K}, {0ul,abs_K}, {0ul,0ul}},
//                       qns);
//                   });
// }

// ExprPtr P(std::int64_t K) {
//   return get_default_context().spbasis() == SPBasis::spinfree ? S(-K) :
//   A(-K);
// }

bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t target_qns,
                    const qns_t source_qns) {
  qns_t qns = source_qns;
  if (op_or_op_product.is<Product>()) {
    const auto& op_product = op_or_op_product.as<Product>();
    for (auto& op_ptr : ranges::views::reverse(op_product.factors())) {
      assert(op_ptr->template is<op_t>());
      const auto& op = op_ptr->template as<op_t>();
      qns = op(qns);
    }
    return qns.overlaps_with(target_qns);
  } else if (op_or_op_product.is<op_t>()) {
    const auto& op = op_or_op_product.as<op_t>();
    qns = op();
    return qns.overlaps_with(target_qns);
  } else
    throw std::invalid_argument(
        "sequant::mbpt::sr::contains_rank(op_or_op_product): op_or_op_product "
        "must be mbpt::sr::op_t or Product thereof");
}

using mbpt::mr::vac_av;

#include <SeQuant/domain/mbpt/vac_av.ipp>

}  // namespace op

}  // namespace mr

// must be defined including op.ipp since it's used there
template <>
bool is_vacuum<mr::qns_t>(mr::qns_t qns) {
  return qns == mr::qns_t{};
}

}  // namespace mbpt

template <Statistics S>
std::wstring to_latex(const mbpt::Operator<mbpt::mr::qns_t, S>& op) {
  using namespace sequant::mbpt;
  using namespace sequant::mbpt::mr;

  auto result = L"{\\hat{" + utf_to_latex(op.label()) + L"}";

  // check if operator has adjoint label, remove if present
  auto base_lbl = sequant::to_wstring(op.label());
  if (base_lbl.back() == adjoint_label) base_lbl.pop_back();

  auto it = label2optype.find(base_lbl);
  OpType optype = OpType::invalid;
  if (it != label2optype.end()) {  // handle special cases
    optype = it->second;
    if (to_class(optype) == OpClass::gen) {
      result += L"}";
      return result;
    }
  }

  // generic operator ... can only handle definite case
  const auto dN = op();
  if (!is_definite(ncre_occ(dN)) || !is_definite(nann_occ(dN)) ||
      !is_definite(ncre_uocc(dN)) || !is_definite(nann_uocc(dN))) {
    throw std::invalid_argument(
        "to_latex(const Operator<qns_t, S>& op): "
        "can only handle  generic operators with definite cre/ann numbers");
  }
  // pure quasiparticle creator/annihilator?
  const auto qprank_cre = nann_occ(dN).lower() + ncre_uocc(dN).lower();
  const auto qprank_ann = ncre_occ(dN).lower() + nann_uocc(dN).lower();
  const auto qppure = qprank_cre == 0 || qprank_ann == 0;
  auto qpaction = to_class(optype);
  if (qppure) {
    if (qprank_cre) {
      // if operator's action implied by the label and actual action agrees, use
      // subscript always
      std::wstring baseline_char = (qpaction != OpClass::deex ? L"_" : L"^");
      if (nann_occ(dN).lower() == ncre_uocc(dN).lower())
        result +=
            baseline_char + L"{" + std::to_wstring(nann_occ(dN).lower()) + L"}";
      else
        result += baseline_char + L"{" + std::to_wstring(nann_occ(dN).lower()) +
                  L"," + std::to_wstring(ncre_uocc(dN).lower()) + L"}";
    } else {
      // if operator's action implied by the label and actual action agrees, use
      // subscript always
      std::wstring baseline_char = (qpaction != OpClass::deex ? L"^" : L"_");
      if (nann_uocc(dN).lower() == ncre_occ(dN).lower()) {
        result +=
            baseline_char + L"{" + std::to_wstring(ncre_occ(dN).lower()) + L"}";
      } else
        result += baseline_char + L"{" + std::to_wstring(ncre_occ(dN).lower()) +
                  L"," + std::to_wstring(nann_uocc(dN).lower()) + L"}";
    }
  } else {  // not pure qp creator/annihilator
    result += L"_{" + std::to_wstring(nann_occ(dN).lower()) + L"," +
              std::to_wstring(ncre_uocc(dN).lower()) + L"}^{" +
              std::to_wstring(ncre_occ(dN).lower()) + L"," +
              std::to_wstring(nann_uocc(dN).lower()) + L"}";
  }
  result += L"}";
  return result;
}

}  // namespace sequant

#include <SeQuant/domain/mbpt/op.ipp>

namespace sequant {
namespace mbpt {
template class Operator<mr::qns_t, Statistics::FermiDirac>;
template class Operator<mr::qns_t, Statistics::BoseEinstein>;
}  // namespace mbpt
}  // namespace sequant
