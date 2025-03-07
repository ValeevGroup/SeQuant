#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <stack>
#include <utility>
#include <vector>

namespace sequant {

class Tensor;

namespace opt {

ExprPtr tail_factor(ExprPtr const& expr) noexcept {
  if (expr->is<Tensor>())
    return expr->clone();

  else if (expr->is<Product>()) {
    auto scalar = expr->as<Product>().scalar();
    if (scalar == 1 && expr->size() == 2) {
      // product with
      //   -single factor that is a tensor
      //   -scalar is just 1
      //  will not be formed because of this block
      return expr->at(1);
    }
    auto facs = ranges::views::tail(*expr);
    return ex<Product>(Product{scalar, ranges::begin(facs), ranges::end(facs)});
  } else {
    // sum
    auto summands = *expr | ranges::views::transform(
                                [](auto const& x) { return tail_factor(x); });
    return ex<Sum>(Sum{ranges::begin(summands), ranges::end(summands)});
  }
}

void pull_scalar(ExprPtr expr) noexcept {
  if (!expr->is<Product>()) return;
  auto& prod = expr->as<Product>();

  auto scal = prod.scalar();
  for (auto&& x : *expr)
    if (x->is<Product>()) {
      auto& p = x->as<Product>();
      scal *= p.scalar();
      p.scale(1 / p.scalar());
    }

  prod.scale(1 / prod.scalar());
  prod.scale(scal);
}

bool has_only_single_atom(const ExprPtr& term) {
  if (term->is_atom()) {
    return true;
  }

  // Recursively check that all elements in the expression tree have only a
  // single element in them. At this point this means checking for Sum or
  // Product objects that only have a single addend or factor respectively.
  return term->size() == 1 && has_only_single_atom(*term->begin());
}

container::vector<container::vector<size_t>> clusters(Sum const& expr) {
  using ranges::views::tail;
  using ranges::views::transform;
  using hash_t = size_t;
  using pos_t = size_t;
  using stack_t = std::stack<pos_t, container::vector<pos_t>>;

  container::map<hash_t, container::set<pos_t>> positions;
  {
    pos_t pos = 0;
    auto visitor = [&positions, &pos](auto const& n) {
      auto h = hash::value(*n);
      if (auto&& found = positions.find(h); found != positions.end()) {
        found->second.emplace(pos);
      } else {
        positions.emplace(h, decltype(positions)::mapped_type{pos});
      }
    };

    for (auto const& term : expr) {
      auto const node = binarize(term);
      if (has_only_single_atom(term)) {
        visitor(node);
      } else {
        node.visit_internal(visitor);
      }
      ++pos;
    }
  }

  container::map<pos_t, container::vector<pos_t>> connections;
  {
    for (auto const& [_, v] : positions) {
      auto const v0 = ranges::front(v);
      auto const v_ = ranges::views::tail(v) |
                      ranges::to<decltype(connections)::mapped_type>;
      if (auto&& found = connections.find(v0); found != connections.end())
        for (auto p : v_) found->second.push_back(p);
      else
        connections.emplace(v0, v_);
    }
  }
  positions.clear();

  container::vector<container::vector<pos_t>> result;
  {
    container::set<pos_t> visited;
    for (auto k : connections | ranges::views::keys)
      if (!visited.contains(k)) {
        stack_t dfs_stack;
        dfs_stack.push(k);
        container::vector<pos_t> clstr;
        while (!dfs_stack.empty()) {
          auto p = dfs_stack.top();
          dfs_stack.pop();
          if (!visited.contains(p)) {
            clstr.push_back(p);
            visited.emplace(p);
          }
          if (auto&& found = connections.find(p); found != connections.end())
            for (auto p_ : ranges::views::reverse(found->second))
              dfs_stack.push(p_);
        }
        result.emplace_back(std::move(clstr));
      }
  }
  return result;
}

Sum reorder(Sum const& sum) {
  Sum result;

  for (auto const& clstr : clusters(sum))
    for (auto p : clstr) result.append(sum.at(p));

  assert(result.size() == sum.size());
  return result;
}

}  // namespace opt

ExprPtr optimize(ExprPtr const& expr, bool reorder_sum) {
  return opt::optimize(
      expr, [](Index const& ix) { return ix.space().approximate_size(); },
      reorder_sum);
}

ExprPtr density_fit_impl(Tensor const& tnsr, Index const& aux_idx) {
  assert(tnsr.bra_rank() == 2     //
         && tnsr.ket_rank() == 2  //
         && tnsr.aux_rank() == 0);

  auto t1 = ex<Tensor>(L"g", bra({ranges::front(tnsr.bra())}),
                       ket({ranges::front(tnsr.ket())}), aux({aux_idx}));

  auto t2 = ex<Tensor>(L"g", bra({ranges::back(tnsr.bra())}),
                       ket({ranges::back(tnsr.ket())}), aux({aux_idx}));

  return ex<Product>(1, ExprPtrList{t1, t2});
}

ExprPtr density_fit(ExprPtr const& expr, std::wstring const& aux_label) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(*expr | transform([&aux_label](auto&& x) {
      return density_fit(x, aux_label);
    }));

  else if (expr->is<Tensor>()) {
    auto const& g = expr->as<Tensor>();
    if (g.label() == L"g"     //
        && g.bra_rank() == 2  //
        && g.ket_rank() == 2  //
        && ranges::none_of(g.indices(), &Index::has_proto_indices))
      return density_fit_impl(expr->as<Tensor>(), Index(aux_label + L"_1"));
    else
      return expr;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());
    size_t aux_ix = 0;
    for (auto&& f : prod.factors())
      if (f.is<Tensor>() && f.as<Tensor>().label() == L"g") {
        auto const& g = f->as<Tensor>();
        auto g_df = density_fit_impl(
            g, Index(aux_label + L"_" + std::to_wstring(++aux_ix)));
        result.append(1, std::move(g_df), Product::Flatten::Yes);
      } else {
        result.append(1, f, Product::Flatten::No);
      }
    return ex<Product>(std::move(result));
  } else
    return expr;
}

ExprPtr csv_transform_impl(Tensor const& tnsr,
                           std::wstring_view coeff_tensor_label) {
  using ranges::views::transform;

  if (ranges::none_of(tnsr.const_braket(), &Index::has_proto_indices))
    return nullptr;

  ////
  auto drop_protos = [](auto&& ixs) {
    return ixs | transform(&Index::drop_proto_indices);
  };
  ////

  if (tnsr.label() == overlap_label()) {
    assert(tnsr.bra_rank() == 1     //
           && tnsr.ket_rank() == 1  //
           && tnsr.aux_rank() == 0);

    auto&& bra_idx = tnsr.bra().at(0);
    auto&& ket_idx = tnsr.ket().at(0);

    auto dummy_idx = suffix_compare(bra_idx, ket_idx)    //
                         ? bra_idx.drop_proto_indices()  //
                         : ket_idx.drop_proto_indices();

    return ex<Product>(
        1,
        ExprPtrList{ex<Tensor>(coeff_tensor_label,                 //
                               bra({bra_idx}), ket({dummy_idx})),  //
                    ex<Tensor>(coeff_tensor_label,                 //
                               bra({dummy_idx}), ket({ket_idx}))});
  }

  Product result;
  result.append(1, ex<Tensor>(tnsr.label(), bra(drop_protos(tnsr.bra())),
                              ket(drop_protos(tnsr.ket())), tnsr.aux()));

  for (auto&& idx : tnsr.bra())
    if (idx.has_proto_indices())
      result.append(1, ex<Tensor>(coeff_tensor_label, bra({idx}),
                                  ket({idx.drop_proto_indices()}), aux({})));
  for (auto&& idx : tnsr.ket())
    if (idx.has_proto_indices())
      result.append(
          1, ex<Tensor>(coeff_tensor_label, bra({idx.drop_proto_indices()}),
                        ket({idx}), aux({})));

  return ex<Product>(std::move(result));
}

ExprPtr csv_aotransform_impl(Tensor const& tnsr,
                             std::wstring_view coeff_tensor_label) {
  using ranges::views::transform;

  if (ranges::none_of(tnsr.const_braket(), &Index::has_proto_indices))
    return nullptr;

  assert(ranges::none_of(tnsr.aux(), &Index::has_proto_indices));
  assert(get_default_context().index_space_registry());
  assert(get_default_context().index_space_registry()->contains(L"μ"));
  auto aospace = get_default_context().index_space_registry()->retrieve(L"μ");

  Product result;
  container::svector<Index> rbra, rket;

  rbra.reserve(tnsr.bra_rank());
  for (auto&& idx : tnsr.bra()) {
    if (idx.has_proto_indices()) {
      Index aoidx = Index::make_tmp_index(aospace);
      result.append(
          1, ex<Tensor>(coeff_tensor_label, bra({idx}), ket({aoidx}), aux({})));
      rbra.emplace_back(std::move(aoidx));
    } else
      rbra.emplace_back(idx);
  }

  rket.reserve(tnsr.ket_rank());
  for (auto&& idx : tnsr.ket()) {
    if (idx.has_proto_indices()) {
      Index aoidx = Index::make_tmp_index(aospace);
      result.append(
          1, ex<Tensor>(coeff_tensor_label, bra({aoidx}), ket({idx}), aux({})));
      rket.emplace_back(std::move(aoidx));
    } else
      rket.emplace_back(idx);
  }

  auto tnsr_ao =
      ex<Tensor>(tnsr.label(), bra(rbra), ket(rket), aux(), tnsr.symmetry(),
                 tnsr.braket_symmetry(), tnsr.particle_symmetry());

  return ex<Product>(std::move(result));
}

ExprPtr csv_transform(ExprPtr const& expr,
                      std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& csv_tensors) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(*expr                              //
                   | transform([&coeff_tensor_label,  //
                                &csv_tensors](auto&& x) {
                       return csv_transform(x, coeff_tensor_label, csv_tensors);
                     }));
  else if (expr->is<Tensor>()) {
    auto const& tnsr = expr->as<Tensor>();
    if (!ranges::contains(csv_tensors, tnsr.label())) return expr;
    if (ranges::none_of(tnsr.indices(), &Index::has_proto_indices)) return expr;
    return csv_transform_impl(tnsr, coeff_tensor_label);
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());

    for (auto&& f : prod.factors()) {
      auto trans = csv_transform(f, coeff_tensor_label, csv_tensors);
      result.append(1, trans ? trans : f,
                    (f->is<Product>() || f->is<Sum>()) ? Product::Flatten::No
                                                       : Product::Flatten::Yes);
    }

    return ex<Product>(std::move(result));

  } else
    return expr;
}

}  // namespace sequant
