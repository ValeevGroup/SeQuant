//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "sr.hpp"

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/wick.hpp"
#include "SeQuant/domain/mbpt/formalism.hpp"

namespace sequant {
namespace mbpt {
namespace sr {

inline constexpr size_t fac(std::size_t n) {
  if (n == 1 || n == 0)
    return 1;
  else
    return n * fac(n - 1);
}

make_op::make_op(std::size_t nbra, std::size_t nket, OpType op)
    : nbra_(nbra), nket_(nket), op_(op) {}

ExprPtr make_op::operator()() const {
  const auto unocc =
      get_default_formalism().sum_over_uocc() == SumOverUocc::Complete
          ? IndexSpace::complete_unoccupied
          : IndexSpace::active_unoccupied;
  const auto occ = IndexSpace::active_occupied;
  return (*this)(unocc, occ);
}

ExprPtr make_op::operator()(IndexSpace::Type unocc,
                            IndexSpace::Type occ) const {
  bool antisymm = get_default_formalism().two_body_interaction() ==
                  TwoBodyInteraction::Antisymm;
  bool csv = get_default_formalism().csv_formalism() == CSVFormalism::CSV;

  // not sure what it means to use nonsymmetric operator if nbra != nket
  if (!antisymm) assert(nbra_ == nket_);

  const auto nbra = nbra_;
  const auto nket = nket_;
  OpType op = op_;
  auto make_idx_vector = [](size_t n, IndexSpace::Type spacetype) {
    auto space = IndexSpace::instance(spacetype);
    std::vector<Index> result;
    result.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      result.push_back(Index::make_tmp_index(space));
    }
    return result;
  };
  auto make_depidx_vector = [](size_t n, IndexSpace::Type spacetype,
                               auto&& protoidxs) {
    auto space = IndexSpace::instance(spacetype);
    std::vector<Index> result;
    result.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      result.push_back(Index::make_tmp_index(space, protoidxs, true));
    }
    return result;
  };
  std::vector<Index> braidxs;
  std::vector<Index> ketidxs;
  if (to_class(op) == OpClass::gen) {
    braidxs = make_idx_vector(nbra, IndexSpace::complete);
    ketidxs = make_idx_vector(nket, IndexSpace::complete);
  } else {
    auto make_occidxs = [&make_idx_vector, &occ](size_t n) {
      return make_idx_vector(n, occ);
    };
    auto make_uoccidxs = [csv, &unocc, &make_idx_vector, &make_depidx_vector](
                             size_t n, auto&& occidxs) {
      return csv ? make_depidx_vector(n, unocc, occidxs)
                 : make_idx_vector(n, unocc);
    };
    if (to_class(op) == OpClass::ex) {
      ketidxs = make_occidxs(nket);
      braidxs = make_uoccidxs(nbra, ketidxs);
    } else {
      braidxs = make_occidxs(nbra);
      ketidxs = make_uoccidxs(nket, braidxs);
    }
  }
  const auto mult = antisymm ? fac(nbra) * fac(nket) : fac(nbra);
  const auto opsymm = antisymm ? Symmetry::antisymm : Symmetry::nonsymm;
  return ex<Constant>(1. / mult) *
         ex<Tensor>(to_wstring(op), braidxs, ketidxs, opsymm) *
         ex<FNOperator>(braidxs, ketidxs, get_default_context().vacuum());
}

make_op Op(OpType _Op, std::size_t Nbra, std::size_t Nket) {
  assert(Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nbra > 0 || Nket_ > 0);
  return make_op{Nbra, Nket_, _Op};
}

#include "sr_op.impl.cpp"

ExprPtr H1() {
  return get_default_context().vacuum() == Vacuum::Physical
             ? Op(OpType::h, 1)()
             : Op(OpType::f, 1)();
}

ExprPtr H2() { return Op(OpType::g, 2)(); }

ExprPtr H0mp() {
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  return H1();
}

ExprPtr H1mp() {
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  return H2();
}

ExprPtr F() { return Op(OpType::f, 1)(); }

ExprPtr W() {
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  return H1mp();
}

ExprPtr H() { return H1() + H2(); }

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> op_connections,
               bool use_top) {
  FWickTheorem wick{expr};
  wick.spinfree(false).use_topology(use_top).set_op_connections(op_connections);
  auto result = wick.compute();
  simplify(result);
  if (Logger::get_instance().wick_stats) {
    std::wcout << "WickTheorem stats: # of contractions attempted = "
               << wick.stats().num_attempted_contractions
               << " # of useful contractions = "
               << wick.stats().num_useful_contractions << std::endl;
  }
  return result;
}

namespace op {

ExprPtr H1() {
  return ex<op_t>(
      [vacuum = get_default_context().vacuum()]() -> std::wstring_view {
        return vacuum == Vacuum::Physical ? L"h" : L"f";
      },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return mbpt::sr::H1();
      },
      [=](qns_t& qns) {
        qns += qns_t{{0, 0}, {-1, +1}};
      });
}

ExprPtr H2() {
  return ex<op_t>([]() -> std::wstring_view { return L"g"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return mbpt::sr::H2();
                  },
                  [=](qns_t& qns) {
                    qns += qns_t{{0, 0}, {-2, +2}};
                  });
}

ExprPtr H() { return H1() + H2(); }

ExprPtr T_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"t"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return mbpt::sr::T_(K);
                  },
                  [=](qns_t& qns) {
                    qns += qns_t{size_t{0}, K};
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

ExprPtr A(std::size_t K) {
  assert(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"A"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return mbpt::sr::A(K);
                  },
                  [=](qns_t& qns) {
                    qns += qns_t{size_t{0}, -K};
                  });
}

bool contains_rank(const ExprPtr& op_or_op_product, qns_t::interval_t k) {
  qns_t qns{0, 0};
  if (op_or_op_product.is<Product>()) {
    const auto& op_product = op_or_op_product.as<Product>();
    for (auto& op_ptr : ranges::views::reverse(op_product.factors())) {
      assert(op_ptr->template is<op_t>());
      const auto& op = op_ptr->template as<op_t>();
      qns = op(qns);
    }
    return qns.overlaps(std::array{qns_t::interval_t{0, 0}, k});
  } else if (op_or_op_product.is<op_t>()) {
    const auto& op = op_or_op_product.as<op_t>();
    qns = op(qns);
    return qns.overlaps(std::array{qns_t::interval_t{0, 0}, k});
  } else
    throw std::invalid_argument(
        "sequant::mbpt::sr::contains_rank(op_or_op_product): op_or_op_product "
        "must be mbpt::sr::op_t or Product thereof");
}

bool contains_up_to_rank(const ExprPtr& op_or_op_product,
                         const unsigned long k) {
  return contains_rank(op_or_op_product, qns_t::interval_t{0ul, k});
}

bool contains_rank(const ExprPtr& op_or_op_product, const unsigned long k) {
  return contains_rank(op_or_op_product, qns_t::interval_t{k, k});
}

ExprPtr vac_av(
    ExprPtr expr,
    std::vector<std::pair<std::wstring, std::wstring>> op_connections) {
  auto vac_av_product = [&op_connections](ExprPtr expr) {
    assert(expr.is<Product>());
    // compute connections
    std::vector<std::pair<int, int>> connections;
    {
      std::map<std::wstring, std::vector<int>>
          oplbl2pos;  // maps operator labels to the operator positions in the
                      // product
      int pos = 0;
      for (const auto& factor : expr.as<Product>()) {
        if (factor.is<op_t>()) {
          const auto& op = factor.as<op_t>();
          const std::wstring op_lbl = std::wstring(op.label());
          const auto it = oplbl2pos.find(op_lbl);
          if (it == oplbl2pos.end()) {  // new label
            oplbl2pos.emplace(op_lbl, std::vector<int>{pos});
          } else {
            it->second.emplace_back(pos);
          }
          ++pos;
        } else if (factor.is<FNOperator>() || factor.is<BNOperator>()) {
          ++pos;
        }
      }

      for (const auto& [op1_lbl, op2_lbl] : op_connections) {
        auto it1 = oplbl2pos.find(op1_lbl);
        auto it2 = oplbl2pos.find(op2_lbl);
        if (it1 == oplbl2pos.end() || it2 == oplbl2pos.end())
          continue;  // one of the op labels is not present in the product
        const auto& [dummy1, op1_indices] = *it1;
        const auto& [dummy2, op2_indices] = *it2;
        for (const auto& op1_idx : op1_indices) {
          for (const auto& op2_idx : op2_indices) {
            using std::min;
            using std::max;
            connections.emplace_back(min(op1_idx, op2_idx),
                                     max(op1_idx, op2_idx));
          }
        }
      }
    }

    // lower to tensor form
    auto lower_to_tensor_form = [](ExprPtr& expr) {
      auto op_lowerer = [](ExprPtr& leaf) {
        if (leaf.is<op_t>()) leaf = leaf.as<op_t>().tensor_form();
      };
      expr->visit(op_lowerer, /* atoms only = */ true);
    };
    lower_to_tensor_form(expr);
    expr = simplify(expr);

    // currently topological equivalence of indices within a normal operator is
    // not detected, assumed based on use_topology ... so turn off use of
    // topology if antisymm=false
    const bool use_topology = (get_default_formalism().two_body_interaction() !=
                               TwoBodyInteraction::Antisymm)
                                  ? false
                                  : true;

    // compute VEV
    return mbpt::sr::vac_av(expr, connections, use_topology);
  };

  ExprPtr result;
  if (expr.is<Product>()) {
    return vac_av_product(expr);
  } else if (expr.is<Sum>()) {
    for (const auto& summand : *expr) {
      auto summand_vev = vac_av(summand, op_connections);
      if (result)
        result += summand_vev;
      else
        result = ex<Sum>(ExprPtrList{summand_vev});
    }
    return result;
  } else if (expr.is<op_t>()) {
    return ex<Constant>(
        0.);  // expectation value of a normal-ordered operator is 0
  } else if (expr.is<Constant>()) {
    return expr;  // vacuum is normalized
  }
  throw std::invalid_argument(
      "mpbt::sr::op::vac_av(expr): unknown expression type");
}

}  // namespace op

}  // namespace sr
}  // namespace mbpt
}  // namespace sequant
