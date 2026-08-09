#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/detail/concepts.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/slotted_index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/overloads.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/swap.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/detail/variant.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/interface.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace sequant::mbpt {

namespace detail {

Index make_index_with_spincase(const Index& idx, mbpt::Spin s) {
  // sanity check: make sure have only one spin label
  SEQUANT_ASSERT(!(idx.label().find(L'↑') != std::wstring::npos &&
                   idx.label().find(L'↓') != std::wstring::npos));

  // to preserve rest of bits first unset spin bit, then set them to the desired
  // state
  auto qns = mbpt::spinannotation_remove(idx.space().qns()).unIon(s);

  IndexSpace space;
  // try looking up space in registry
  const auto label = mbpt::spinannotation_replacе(idx.space().base_key(), s);
  if (auto isr = get_default_context().index_space_registry()) {
    auto* space_ptr = isr->retrieve_ptr(label);
    if (space_ptr && space_ptr->type() == idx.space().type() &&
        space_ptr->qns() == qns) {
      space = *space_ptr;
    }
  }
  // if space not found, construct
  if (!space) {
    space = IndexSpace{label, idx.space().type(), qns,
                       // N.B. assume size does not depend on spin
                       idx.space().approximate_size()};
  }
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_spincase(pidx, s);
  return Index{space, idx.ordinal(), protoindices};
}

// The argument really should be non-const but const semantics are broken
// for the ExprPtr type so we are required to make this const in order
// to be able to use this function everywhere we want to.
void reset_idx_tags(const ExprPtr& expr) {
  expr->visit(
      [](ExprPtr& current) {
        if (current.is<AbstractTensor>()) {
          current.as<AbstractTensor>()._reset_tags();
        }
      },
      true);
}

template <typename Container, typename TraceFunction, typename... Args>
[[nodiscard]] Container wrap_trace(const ResultExpr& expr,
                                   TraceFunction&& tracer, Args&&... args) {
  bool searchForNonEquivalentResults = expr.symmetry() != Symmetry::Nonsymm;
  searchForNonEquivalentResults &=
      expr.bra().size() > 1 || expr.ket().size() > 1;
  const bool brasSameSpace = std::all_of(
      expr.bra().begin(), expr.bra().end(),
      [&](const Index& idx) { return idx.space() == expr.bra()[0].space(); });
  const bool ketsSameSpace = std::all_of(
      expr.ket().begin(), expr.ket().end(),
      [&](const Index& idx) { return idx.space() == expr.ket()[0].space(); });
  searchForNonEquivalentResults &= !brasSameSpace && !ketsSameSpace;

  if (!searchForNonEquivalentResults) {
    ResultExpr traced = expr.clone();
    traced.expression() = tracer(
        traced.expression(),
        traced.index_particle_grouping<container::svector<SlottedIndex>>(),
        std::forward<Args>(args)...);

    traced.set_symmetry(Symmetry::Nonsymm);

    return {std::move(traced)};
  }

  SEQUANT_ASSERT(expr.symmetry() == Symmetry::Antisymm ||
                 expr.symmetry() == Symmetry::Symm);

  // TODO: Do we have to track the sign?
  const bool permuteBra = expr.bra().size() >= expr.ket().size();
  auto permIndices = permuteBra ? expr.bra() : expr.ket();
  const std::size_t unchangedSize =
      permuteBra ? expr.ket().size() : expr.bra().size();

  [[maybe_unused]] auto get_phase = [](auto container) {
    reset_ts_swap_counter<Index>();
    bubble_sort(container.begin(), container.end(), std::less<Index>{});
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  reset_ts_swap_counter<Index>();
  bubble_sort(permIndices.begin(), permIndices.end(), std::less<Index>{});
  const int initialSign = ts_swap_counter_is_even<Index>() ? 1 : -1;
  const auto originalIndices = permIndices;

  container::svector<container::set<std::pair<IndexSpace, IndexSpace>>>
      idxPairings;

  Container resultSet;

  // For next_permutation to work in this context, permIndices must be sorted
  SEQUANT_ASSERT(std::is_sorted(permIndices.begin(), permIndices.end()));

  int sign = initialSign;
  do {
    const int currentSign = sign;
    // std::next_permutation creates one lexicographical permutation after the
    // other, which should imply that the phase should alternate between
    // iterations.
    sign *= -1;
    SEQUANT_ASSERT(currentSign == get_phase(permIndices) * initialSign);

    container::set<std::pair<IndexSpace, IndexSpace>> currentPairing;

    for (std::size_t i = 0; i < unchangedSize; ++i) {
      if (permuteBra) {
        currentPairing.insert(
            std::make_pair(permIndices[i].space(), expr.ket()[i].space()));
      } else {
        currentPairing.insert(
            std::make_pair(expr.bra()[i].space(), permIndices[i].space()));
      }
    }

    for (std::size_t i = unchangedSize; i < permIndices.size(); ++i) {
      currentPairing.insert(
          std::make_pair(permIndices[i].space(), IndexSpace::null));
    }

    if (std::find(idxPairings.begin(), idxPairings.end(), currentPairing) !=
        idxPairings.end()) {
      continue;
    }

    // Found a new index pairing

    ExprPtr expression = expr.expression().clone();

    expression *= ex<Constant>(currentSign);
    expression = simplify(expression);

    ResultExpr result = [&]() {
      SEQUANT_ASSERT(expr.has_label());
      if (permuteBra) {
        return ResultExpr(bra(permIndices), ket(expr.ket()), aux(expr.aux()),
                          expr.symmetry(), expr.braket_symmetry(),
                          expr.column_symmetry(), expr.label(),
                          std::move(expression));
      } else {
        return ResultExpr(bra(expr.bra()), ket(permIndices), aux(expr.aux()),
                          expr.symmetry(), expr.braket_symmetry(),
                          expr.column_symmetry(), expr.label(),
                          std::move(expression));
      }
    }();

    result.expression() = tracer(
        result.expression(),
        result.index_particle_grouping<container::svector<SlottedIndex>>(),
        std::forward<Args>(args)...);

    result.set_symmetry(Symmetry::Nonsymm);

    resultSet.push_back(std::move(result));
  } while (std::next_permutation(permIndices.begin(), permIndices.end()));

  return resultSet;
}

}  // namespace detail

Index make_spinalpha(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::alpha);
};

Index make_spinbeta(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::beta);
};

Index make_spinfree(const Index& idx) {
  return detail::make_index_with_spincase(idx, mbpt::Spin::any);
};

ExprPtr swap_bra_ket(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Lambda for tensor
  auto tensor_swap = [](const Tensor& tensor) {
    return ex<Tensor>(tensor.label(), bra(tensor.ket().value()),
                      ket(tensor.bra().value()), tensor.symmetry(),
                      tensor.braket_symmetry(), tensor.column_symmetry());
  };

  // Lambda for product
  auto product_swap = [&tensor_swap](const Product& product) {
    auto result = std::make_shared<Product>();
    result->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        result->append(1, tensor_swap(term->as<Tensor>()),
                       Product::Flatten::No);
      } else if (term->is<Variable>() || term->is<Constant>()) {
        result->append(1, term);
      } else {
        throw Exception("Invalid Expr type in product_swap: " +
                        term->type_name());
      }
    }
    return result;
  };

  if (expr->is<Tensor>())
    return tensor_swap(expr->as<Tensor>());
  else if (expr->is<Product>())
    return product_swap(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      result->append(swap_bra_ket(term));
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in swap_bra_ket: " + expr->type_name());
  }
}

ExprPtr append_spin(const ExprPtr& expr,
                    const container::map<Index, Index>& index_replacements) {
  auto add_spin_to_tensor = [&index_replacements](const Tensor& tensor) {
    auto spin_tensor = std::make_shared<Tensor>(tensor);
    spin_tensor->transform_indices(index_replacements);
    return spin_tensor;
  };

  auto add_spin_to_product = [&add_spin_to_tensor](const Product& product) {
    auto spin_product = std::make_shared<Product>();
    spin_product->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        spin_product->append(1, add_spin_to_tensor(term->as<Tensor>()));
      } else if (term->is_scalar()) {
        spin_product->append(1, term);
      } else {
        throw Exception(
            "Invalid Expr type in append_spin::add_spin_to_product: " +
            term->type_name());
      }
    }
    return spin_product;
  };

  if (expr->is<Tensor>()) {
    return add_spin_to_tensor(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return add_spin_to_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    auto spin_expr = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      spin_expr->append(append_spin(summand, index_replacements));
    }
    return spin_expr;
  } else if (expr->is<Constant>() || expr->is<Variable>()) {
    return expr;
  }

  throw Exception("Unsupported Expr type in append_spin");
}

ExprPtr remove_spin(const ExprPtr& expr) {
  auto remove_spin_from_tensor = [](const Tensor& tensor) {
    container::svector<Index> b(tensor.bra().begin(), tensor.bra().end());
    container::svector<Index> k(tensor.ket().begin(), tensor.ket().end());
    {
      for (auto&& idx : ranges::views::concat(b, k)) {
        idx = make_spinfree(idx);
      }
    }
    // column symmetry must survive spin removal for triplet doubles
    return ex<Tensor>(tensor.label(), bra(std::move(b)), ket(std::move(k)),
                      tensor.aux(), tensor.symmetry(), tensor.braket_symmetry(),
                      tensor.column_symmetry());
  };

  auto remove_spin_from_product =
      [&remove_spin_from_tensor](const Product& product) {
        auto result = std::make_shared<Product>();
        result->scale(product.scalar());
        for (auto&& term : product) {
          if (term->is<Tensor>()) {
            result->append(1, remove_spin_from_tensor(term->as<Tensor>()));
          } else if (term->is<Constant>() || term->is<Variable>()) {
            result->append(1, term);
          } else {
            throw Exception(
                "Invalid Expr type in remove_spin::remove_spin_from_product: " +
                term->type_name());
          }
        }
        return result;
      };

  if (expr->is<Tensor>()) {
    return remove_spin_from_tensor(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return remove_spin_from_product(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(remove_spin(summand));
    }
    return result;
  } else if (expr->is<Constant>() || expr->is<Variable>()) {
    return expr;
  } else {
    throw Exception("Invalid Expr type in remove_spin: " + expr->type_name());
  }
}

ExprPtr remove_spin_with_relabel(const ExprPtr& expr) {
  auto remove_spin_from_tensor_with_relabel = [](const Tensor& tensor) {
    return remove_spin(ex<Tensor>(tensor));
  };

  auto remove_spin_from_product_with_relabel =
      [](const Product& product) -> ExprPtr {
    // count occurrences in the current product to detect contracted indices.
    container::map<Index, std::size_t> index_occurrences;
    for (const auto& factor : product) {
      if (!factor->is<Tensor>()) continue;
      const auto& tensor = factor->as<Tensor>();
      for (const auto& idx : tensor.const_braket_indices()) {
        ++index_occurrences[idx];
      }
      for (const auto& idx : tensor.aux()) {
        ++index_occurrences[idx];
      }
    }

    // collision detection
    // group distinct spin-labeled indices by their would-be spin-free identity.
    container::map<Index, container::set<Index>> sf_to_originals;
    container::map<IndexSpace, Index::ordinal_type> max_ordinal_per_space;
    for (const auto& [idx, _] : index_occurrences) {
      const auto sf_idx = make_spinfree(idx);
      sf_to_originals[sf_idx].insert(idx);
      if (sf_idx.ordinal()) {
        const auto& sf_space = sf_idx.space();
        const auto ord = sf_idx.ordinal().value();
        auto it = max_ordinal_per_space.find(sf_space);
        if (it == max_ordinal_per_space.end() || ord > it->second) {
          max_ordinal_per_space[sf_space] = ord;
        }
      } else if (!max_ordinal_per_space.contains(sf_idx.space())) {
        max_ordinal_per_space[sf_idx.space()] = 0;
      }
    }

    container::map<Index, Index> relabel_map;
    for (const auto& [sf_idx, originals] : sf_to_originals) {
      if (originals.size() <= 1) continue;

      // keep external index unchanged; relabel only internals.
      auto keep_it = originals.begin();
      for (auto it = originals.begin(); it != originals.end(); ++it) {
        const auto occ_it = index_occurrences.find(*it);
        SEQUANT_ASSERT(occ_it != index_occurrences.end());
        if (occ_it->second <= 1) {
          keep_it = it;
          break;
        }
      }

      for (auto it = originals.begin(); it != originals.end(); ++it) {
        if (it == keep_it) continue;

        const auto occ_it = index_occurrences.find(*it);
        SEQUANT_ASSERT(occ_it != index_occurrences.end());
        if (occ_it->second <= 1) continue;  // do not relabel external indices

        auto& next_ordinal = max_ordinal_per_space[sf_idx.space()];
        ++next_ordinal;
        if (next_ordinal >= Index::min_tmp_index()) {
          throw Exception(
              "remove_spin_with_relabel: exhausted non-reserved "
              "ordinals while resolving spin collisions");
        }

        relabel_map[*it] = Index(it->space(), next_ordinal, it->proto_indices(),
                                 it->symmetric_proto_indices());
      }
    }

    // if (!relabel_map.empty()) {
    //   std::wcout << "relabel_map non-empty for product:\n";
    //   std::wcout << to_latex_align(ex<Product>(product)) << "\n";
    //   for (const auto& [from, to] : relabel_map) {
    //     std::wcout << "  " << from.to_latex() << " -> " << to.to_latex()
    //                << "\n";
    //   }
    // }

    Product relabeled_product = product;
    if (!relabel_map.empty()) {
      for (auto& factor : relabeled_product.factors()) {
        if (factor->is<Tensor>()) {
          auto tensor = factor->as<Tensor>();
          tensor.transform_indices(relabel_map);
          factor = ex<Tensor>(tensor);
        }
      }
    }

    return remove_spin(ex<Product>(std::move(relabeled_product)));
  };

  if (expr->is<Tensor>()) {
    return remove_spin_from_tensor_with_relabel(expr->as<Tensor>());
  } else if (expr->is<Product>()) {
    return remove_spin_from_product_with_relabel(expr->as<Product>());
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(remove_spin_with_relabel(summand));
    }
    return result;
  } else if (expr->is<Constant>() || expr->is<Variable>()) {
    return expr;
  } else {
    throw Exception("Invalid Expr type in remove_spin_with_relabel: " +
                    expr->type_name());
  }
}

bool ms_conserving_columns(const AbstractTensor& tensor) {
  for (const auto& [bra, ket] :
       ranges::zip_view(tensor._bra(), tensor._ket())) {
    if (bra.nonnull() && ket.nonnull()) {
      const auto bra_ms = mbpt::to_spin(bra.space().qns());
      ;
      const auto ket_ms = mbpt::to_spin(ket.space().qns());
      ;
      if (bra_ms != ket_ms) return false;
    }
  }
  return true;
}

bool ms_uniform_tensor(const AbstractTensor& tensor) {
  auto braket = tensor._braket();
  SEQUANT_ASSERT(ranges::empty(braket) == false);
  std::optional<mbpt::Spin> ms;
  return ranges::all_of(braket, [&ms](const auto& idx) {
    if (idx.nonnull()) {
      const auto idx_ms = mbpt::to_spin(idx.space().qns());
      if (ms.has_value())
        return idx_ms == *ms;
      else {
        ms = idx_ms;
        return true;
      }
    } else
      return true;
  });
}

bool can_expand(const AbstractTensor& tensor) {
  [[maybe_unused]] auto all_have_spin =
      ranges::all_of(tensor._braket(), [](const auto& idx) {
        auto idx_spin = mbpt::to_spin(idx.space().qns());
        return idx_spin == mbpt::Spin::alpha || idx_spin == mbpt::Spin::beta;
      });
  SEQUANT_ASSERT(ranges::all_of(tensor._braket(), [](const auto& idx) {
    auto idx_spin = mbpt::to_spin(idx.space().qns());
    return idx_spin == mbpt::Spin::alpha || idx_spin == mbpt::Spin::beta;
  }));

  if (tensor._bra_rank() == tensor._ket_rank()) {
    // particle-conserving: column-wise check
    auto is_alpha = [](const Index& idx) {
      return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha;
    };

    // count alpha indices in bra
    auto a_bra = ranges::count_if(tensor._bra(), is_alpha);

    // count alpha indices in ket
    auto a_ket = ranges::count_if(tensor._ket(), is_alpha);

    return a_bra == a_ket;
  } else {
    // NPC: global count over paired columns only
    auto n_paired = std::min(tensor._bra_rank(), tensor._ket_rank());
    auto is_alpha = [](const Index& idx) {
      return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha;
    };
    auto a_bra = ranges::count_if(tensor._bra() | ranges::views::take(n_paired),
                                  is_alpha);
    auto a_ket = ranges::count_if(tensor._ket() | ranges::views::take(n_paired),
                                  is_alpha);
    return a_bra == a_ket;
  }
}

ExprPtr expand_antisymm(const Tensor& tensor, bool skip_spinsymm) {
  if (tensor.bra_rank() <= 1 && tensor.ket_rank() <= 1) {
    return std::make_shared<Tensor>(Tensor(
        tensor.label(), tensor.bra(), tensor.ket(), tensor.aux(),
        Symmetry::Nonsymm, tensor.braket_symmetry(), tensor.column_symmetry()));
  }

  if (skip_spinsymm && ms_uniform_tensor(tensor)) {
    return std::make_shared<Tensor>(tensor);
  }

  if (tensor.symmetry() != Symmetry::Antisymm) {
    return std::make_shared<Tensor>(tensor);
  }

  auto get_phase = [](const Tensor& t) {
    container::svector<Index> b(t.bra().begin(), t.bra().end());
    container::svector<Index> k(t.ket().begin(), t.ket().end());
    reset_ts_swap_counter<Index>();
    bubble_sort(std::begin(b), std::end(b));
    bubble_sort(std::begin(k), std::end(k));
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  const auto prefactor = get_phase(tensor);
  const bool permute_bra = tensor.bra_rank() >= tensor.ket_rank();

  auto expand = [&]() {
    auto expr_sum = std::make_shared<Sum>();
    container::set<Index> perm_list(
        permute_bra ? tensor.bra().begin() : tensor.ket().begin(),
        permute_bra ? tensor.bra().end() : tensor.ket().end());

    do {
      auto col_symm = tensor.column_symmetry();
      if (tensor.label() == L"R") {
        auto bra_r = tensor.bra_rank();
        auto ket_r = tensor.ket_rank();
        auto diff = (bra_r > ket_r) ? (bra_r - ket_r) : (ket_r - bra_r);
        if (diff >= 2) {
          col_symm = ColumnSymmetry::Nonsymm;
        }
      }

      auto new_tensor =
          permute_bra
              ? Tensor(tensor.label(), bra(perm_list), ket(tensor.ket()),
                       tensor.aux(), Symmetry::Nonsymm,
                       tensor.braket_symmetry(), col_symm)
              : Tensor(tensor.label(), bra(tensor.bra()), ket(perm_list),
                       tensor.aux(), Symmetry::Nonsymm,
                       tensor.braket_symmetry(), col_symm);

      if (!ms_conserving_columns(new_tensor)) continue;
      auto prod = std::make_shared<Product>();
      prod->append(get_phase(new_tensor), ex<Tensor>(new_tensor));
      prod->scale(prefactor);
      expr_sum->append(prod);
    } while (std::next_permutation(perm_list.begin(), perm_list.end()));

    return expr_sum;
  };

  return expand();
}

ExprPtr expand_antisymm(const ExprPtr& expr, bool skip_spinsymm) {
  if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else if (expr->is<Tensor>())
    return expand_antisymm(expr->as<Tensor>(), skip_spinsymm);

  // Product lambda
  auto expand_product = [&skip_spinsymm](const Product& expr) {
    Product temp{};
    temp.scale(expr.scalar());
    for (auto&& term : expr) {
      if (term->is<Tensor>()) {
        temp.append(1, expand_antisymm(term->as<Tensor>(), skip_spinsymm),
                    Product::Flatten::No);
      } else if (term->is_scalar()) {
        temp.append(1, term, Product::Flatten::No);
      } else {
        throw Exception(
            "Invalid Expr type in expand_antisymm::expand_product: " +
            term->type_name());
      }
    }
    ExprPtr result = std::make_shared<Product>(temp);
    rapid_simplify(result);
    return result;
  };

  if (expr->is<Product>())
    return expand_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      result->append(expand_antisymm(term, skip_spinsymm));
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in expand_antisymm: " +
                    expr->type_name());
  }
}

container::svector<container::map<Index, Index>> A_maps(const Tensor& A) {
  SEQUANT_ASSERT(A.label() == reserved::antisymm_label());

  container::svector<std::size_t> bra_indices(A.bra_rank());
  container::svector<std::size_t> ket_indices(A.ket_rank());
  std::iota(bra_indices.begin(), bra_indices.end(), 0);
  std::iota(ket_indices.begin(), ket_indices.end(), 0);

  container::svector<container::map<Index, Index>> result;

  do {
    do {
      container::map<Index, Index> current_replacements;

      for (std::size_t i = 0; i < bra_indices.size(); ++i) {
        current_replacements.emplace(A.bra()[i], A.bra()[bra_indices[i]]);
      }
      for (std::size_t i = 0; i < ket_indices.size(); ++i) {
        current_replacements.emplace(A.ket()[i], A.ket()[ket_indices[i]]);
      }

      result.push_back(std::move(current_replacements));
    } while (std::next_permutation(bra_indices.begin(), bra_indices.end()));
  } while (std::next_permutation(ket_indices.begin(), ket_indices.end()));

  return result;
}

ExprPtr expand_A_op(const ProductPtr& product) {
  bool has_A_operator = false;

  // Check A and build replacement map
  container::svector<container::map<Index, Index>> map_list;
  Tensor A_tensor;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      auto A = term->as<Tensor>();
      if (A.label() == reserved::antisymm_label() && A.bra_rank() <= 1 &&
          A.ket_rank() <= 1) {
        return remove_tensor(product, reserved::antisymm_label());
      } else if ((A.label() == reserved::antisymm_label())) {
        has_A_operator = true;
        A_tensor = A;
        map_list = A_maps(A);
        break;
      }
    }
  }

  if (!has_A_operator) return product;

  const auto nf = rational{
      1, (factorial(A_tensor.bra_rank()) * factorial(A_tensor.ket_rank()))};

  auto new_result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    // Get phase of the transformation
    int phase;
    {
      container::svector<Index> transformed_list;
      for (const auto& [key, val] : map) transformed_list.push_back(val);

      reset_ts_swap_counter<Index>();
      bubble_sort(std::begin(transformed_list), std::end(transformed_list));
      phase = ts_swap_counter_is_even<Index>() ? 1 : -1;
    }

    ProductPtr new_product = std::make_shared<Product>();
    new_product->scale(product->scalar());
    auto temp_product = remove_tensor(product, reserved::antisymm_label());
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product->append(1, ex<Tensor>(new_tensor));
      } else {
        new_product->append(1, term);
      }
    }
    new_product->scale(phase * nf);
    new_result->append(new_product);
  }  // map_list

  detail::reset_idx_tags(new_result);

  return new_result;
}

ExprPtr symmetrize_expr(const ProductPtr& product) {
  auto result = std::make_shared<Sum>();

  // Drops canonical-order assumption; handles arbitrary sequence and variables.
  const auto& factors = product->factors();
  auto it = ranges::find_if(factors, [](const ExprPtr& factor) {
    return factor->is<Tensor>() &&
           factor->as<Tensor>().label() == reserved::antisymm_label();
  });
  if (it == ranges::end(factors)) return product;
  const auto& A_tensor = (*it)->as<Tensor>();
  SEQUANT_ASSERT(A_tensor.label() == reserved::antisymm_label());

  auto A_is_nconserving = A_tensor.bra_rank() == A_tensor.ket_rank();

  if (A_is_nconserving && A_tensor.bra_rank() == 1)
    return remove_tensor(product, reserved::antisymm_label());

  SEQUANT_ASSERT(A_tensor.rank() > 1);

  auto S = Tensor{};
  if (A_is_nconserving) {
    S = Tensor(reserved::symm_label(), A_tensor.bra(), A_tensor.ket(),
               A_tensor.aux(), Symmetry::Nonsymm);
  } else {  // A is N-nonconserving
    auto n = std::min(A_tensor.bra_rank(), A_tensor.ket_rank());
    container::svector<Index> bra_list(A_tensor.bra().begin(),
                                       A_tensor.bra().begin() + n);
    container::svector<Index> ket_list(A_tensor.ket().begin(),
                                       A_tensor.ket().begin() + n);
    S = Tensor(reserved::symm_label(), bra(std::move(bra_list)),
               ket(std::move(ket_list)), A_tensor.aux(), Symmetry::Nonsymm);
  }
  const auto nf = rational{1, factorial(S.ket_rank())};

  // Generate replacement maps from a list of Index type (could be a bra or a
  // ket)
  // Uses a permuted list of int to generate permutations
  // TODO factor out for reuse
  auto maps_from_list = [](const container::svector<Index>& list) {
    container::svector<int> int_list(list.size());
    std::iota(int_list.begin(), int_list.end(), 0);
    container::svector<container::map<Index, Index>> result;
    do {
      container::map<Index, Index> map;
      auto list_ptr = list.begin();
      for (auto&& i : int_list) {
        map.emplace(*list_ptr, list[i]);
        list_ptr++;
      }
      result.push_back(map);
    } while (std::next_permutation(int_list.begin(), int_list.end()));
    SEQUANT_ASSERT(result.size() ==
                   boost::numeric_cast<size_t>(factorial(list.size())));
    return result;
  };

  // Get phase relative to the canonical order
  // TODO factor out for reuse
  auto get_phase = [](const container::map<Index, Index>& map) {
    container::svector<Index> idx_list;
    for (const auto& [key, val] : map) idx_list.push_back(val);
    reset_ts_swap_counter<Index>();
    bubble_sort(std::begin(idx_list), std::end(idx_list));
    return ts_swap_counter_is_even<Index>() ? 1 : -1;
  };

  container::svector<container::map<Index, Index>> maps;
  // CASE 1: n_bra = n_ket on all tensors
  if (A_is_nconserving) {
    maps = maps_from_list(A_tensor.bra());
  } else {
    SEQUANT_ASSERT(A_tensor.bra_rank() != A_tensor.ket_rank());
    maps = A_tensor.bra_rank() > A_tensor.ket_rank()
               ? maps_from_list(A_tensor.bra())
               : maps_from_list(A_tensor.ket());
  }
  SEQUANT_ASSERT(!maps.empty());
  for (auto&& map : maps) {
    Product new_product{};
    new_product.scale(product->scalar() * nf);
    new_product.append(get_phase(map), ex<Tensor>(S));
    auto temp_product = remove_tensor(product, reserved::antisymm_label());
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_product.append(1, ex<Tensor>(new_tensor));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        new_product.append(1, term);
      } else {
        throw Exception("Invalid Expr type in symmetrize_expr: " +
                        term->type_name());
      }
    }
    result->append(ex<Product>(new_product));
  }  // map
  return result;
}

ExprPtr symmetrize_expr(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  if (expr->is<Product>())
    return symmetrize_expr(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(symmetrize_expr(summand));
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in symmetrize_expr: " +
                    expr->type_name());
  }
}

ExprPtr expand_A_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  if (expr->is<Product>())
    return expand_A_op(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      result->append(expand_A_op(summand));
    }
    return result;
  }

  throw Exception("Invalid Expr type in expand_A_op: " + expr->type_name());
}

container::svector<container::map<Index, Index>> P_maps(const Tensor& P) {
  SEQUANT_ASSERT(P.label() == reserved::transposition_label());

  // Return pair-wise replacements
  // P_ij -> {{i,j},{j,i}}
  // P_ijkl \equiv P_ij P_kl -> {{i,j},{j,i},{k,l},{l,k}}
  // P_ij^ab \equiv P_ij P^ab -> {{i,j},{j,i},{a,b},{b,a}}
  SEQUANT_ASSERT(P.bra_rank() % 2 == 0 && P.ket_rank() % 2 == 0);
  container::map<Index, Index> idx_rep;
  auto indices = P.const_braket_indices();
  for (auto it = indices.begin(); it != indices.end(); ranges::advance(it, 2)) {
    auto& idx1 = *it;
    auto it_next = it;
    ++it_next;
    SEQUANT_ASSERT(it_next != indices.end());
    auto& idx2 = *it_next;
    idx_rep.emplace(idx1, idx2);
    idx_rep.emplace(idx2, idx1);
  }

  SEQUANT_ASSERT(idx_rep.size() == (P.bra_net_rank() + P.ket_net_rank()));
  return container::svector<container::map<Index, Index>>{idx_rep};
}

ExprPtr expand_P_op(const ProductPtr& product) {
  bool has_P_operator = false;

  // Check P and build a replacement map
  // Assuming a product can have multiple P operators
  container::svector<container::map<Index, Index>> map_list;
  for (auto& term : product) {
    if (term->is<Tensor>()) {
      const auto& P = term->as<Tensor>();
      if ((P.label() == reserved::transposition_label()) &&
          (P.bra_rank() > 1 || (P.ket_rank() > 1))) {
        has_P_operator = true;
        auto map = P_maps(P);
        map_list.insert(map_list.end(), map.begin(), map.end());
      } else if ((P.label() == reserved::transposition_label()) &&
                 (P.bra_rank() == 1 && (P.ket_rank() == 1))) {
        return remove_tensor(product, reserved::transposition_label());
      }
    }
  }

  if (!has_P_operator) return product;

  auto result = std::make_shared<Sum>();
  for (auto&& map : map_list) {
    ProductPtr new_product = std::make_shared<Product>();
    new_product->scale(product->scalar());
    auto temp_product = remove_tensor(product, reserved::transposition_label());
    for (auto&& term : *temp_product) {
      if (term->is<Tensor>()) {
        auto new_tensor = term->as<Tensor>();
        new_tensor.transform_indices(map);
        new_tensor.reset_tags();
        new_product->append(1, ex<Tensor>(new_tensor));
      } else if (term->is_scalar()) {
        new_product->append(1, term);
      } else {
        throw Exception("Invalid Expr type in expand_P_op: " +
                        term->type_name());
      }
    }
    result->append(new_product);
  }  // map_list

  return result;
}

ExprPtr expand_P_op(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;
  else if (expr->is<Product>())
    return expand_P_op(expr.as_shared_ptr<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto& summand : *expr) {
      result->append(expand_P_op(summand));
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in expand_P_op: " + expr->type_name());
  }
}

container::svector<container::map<Index, Index>> S_replacement_maps(
    const Tensor& S) {
  SEQUANT_ASSERT(S.label() == reserved::symm_label());
  SEQUANT_ASSERT(S.bra_rank() > 1);
  SEQUANT_ASSERT(S.bra().size() == S.ket().size());
  container::svector<int> int_list(S.bra().size());
  std::iota(std::begin(int_list), std::end(int_list), 0);

  container::svector<container::map<Index, Index>> maps;
  do {
    container::map<Index, Index> map;
    auto S_bra_ptr = S.bra().begin();
    auto S_ket_ptr = S.ket().begin();
    for (auto&& i : int_list) {
      map.emplace(*S_bra_ptr, S.bra()[i]);
      ++S_bra_ptr;
      map.emplace(*S_ket_ptr, S.ket()[i]);
      ++S_ket_ptr;
    }
    maps.push_back(map);
  } while (std::next_permutation(int_list.begin(), int_list.end()));

  return maps;
}

ExprPtr S_maps(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>() || expr->is<Tensor>())
    return expr;

  auto result = std::make_shared<Sum>();

  // Check if S operator is present
  if (!has_tensor(expr, reserved::symm_label())) return expr;

  detail::reset_idx_tags(expr);

  // Lambda for applying S on products
  auto expand_S_product = [](const ProductPtr& product) -> ExprPtr {
    // check if S is present
    if (!has_tensor(product, reserved::symm_label())) return product;

    container::svector<container::map<Index, Index>> maps;
    // supports arbitrary sequence and variables
    std::size_t S_tensor_rank;
    for (auto&& factor : product->factors()) {
      if (factor->is<Tensor>() &&
          factor->as<Tensor>().label() == reserved::symm_label()) {
        auto factor_as_tensor = factor.as<Tensor>();
        S_tensor_rank = factor_as_tensor.ket_rank();
        maps = S_replacement_maps(factor_as_tensor);
        break;
      }
    }
    SEQUANT_ASSERT(!maps.empty());
    const auto nf = rational{1, factorial(S_tensor_rank)};

    Sum sum{};
    for (auto&& map : maps) {
      ProductPtr new_product = std::make_shared<Product>();
      new_product->scale(product->scalar() * nf);
      auto temp_product = remove_tensor(product, reserved::symm_label())
                              .as_shared_ptr<Product>();
      for (auto&& term : temp_product) {
        if (term->is<Tensor>()) {
          auto new_tensor = term->as<Tensor>();
          new_tensor.transform_indices(map);
          new_product->append(1, ex<Tensor>(new_tensor));
        } else {
          new_product->append(1, term);
        }
      }
      sum.append(new_product);
    }
    ExprPtr result = std::make_shared<Sum>(sum);
    return result;
  };

  if (expr->is<Product>()) {
    result->append(expand_S_product(expr.as_shared_ptr<Product>()));
  } else if (expr->is<Sum>()) {
    for (auto&& term : *expr) {
      if (term->is<Product>()) {
        result->append(expand_S_product(term.as_shared_ptr<Product>()));
      } else if (term->is<Tensor>() || term->is<Constant>() ||
                 expr->is<Variable>()) {
        result->append(term);
      }
    }
  }

  detail::reset_idx_tags(result);
  return result;
}

template <detail::index_group_range IdxGroups>
ExprPtr closed_shell_spintrace_impl(const ExprPtr& expression,
                                    IdxGroups&& ext_index_groups,
                                    bool full_expansion) {
  // Symmetrize and expression
  // Partially expand the antisymmetrizer and write it in terms of S operator.
  // See symmetrize_expr(expr) function for implementation details. We want an
  // expression with non-symmetric tensors, hence we are partially expanding the
  // antisymmetrizer (A) and fully expanding the anti-symmetric tensors to
  // non-symmetric.
  // full_expansion: it fully expands the antisymmetrizer directly (can be used
  // for v2 eqs, however it is not an optimized way).
  auto partially_or_fully_expand = [&full_expansion](const ExprPtr& expr) {
    auto temp = expr;
    if (has_tensor(temp, reserved::antisymm_label())) {
      if (full_expansion) {
        temp = expand_A_op(temp);
      } else {
        temp = symmetrize_expr(temp);
      }
    }
    temp = expand_antisymm(temp);
    rapid_simplify(temp);
    return temp;
  };
  ExprPtr expr = partially_or_fully_expand(expression);

  // Index tags are cleaned prior to calling the fast canonicalizer
  detail::reset_idx_tags(expr);  // This call is REQUIRED
  expand(expr);                  // This call is REQUIRED
  simplify(expr);  // full simplify to combine terms before count_cycles

  // Lambda for spin-tracing a product term
  // For closed-shell case, a spin-traced result is a product term scaled by
  // 2^{n_cycles}, where n_cycles are counted by the lambda function described
  // above. For every product term, the bra indices on all tensors are merged
  // into a single list, so are the ket indices. External indices are
  // substituted with either one of the index (because the two vectors should be
  // permutations of each other to count cycles). All tensors must be nonsymm.
  auto trace_product = [&ext_index_groups](const Product& product) {
    // Remove S if present in a product
    Product temp_product{};
    temp_product.scale(product.scalar());
    const auto& factors = product.factors();
    auto is_symm_tensor = [](const auto& factor) {
      return factor->template is<Tensor>() &&
             factor->template as<Tensor>().label() == reserved::symm_label();
    };

    if (std::ranges::any_of(factors, is_symm_tensor)) {
      for (auto&& factor : factors) {
        if (!is_symm_tensor(factor)) {
          temp_product.append(1, factor, Product::Flatten::No);
        }
      }
    } else {
      temp_product = product;
    }

    const bool collect_symmetrizer_indices = ext_index_groups.empty();

    auto get_ket_indices = [&](const Product& prod) {
      container::svector<Index> ket_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>() &&
            (collect_symmetrizer_indices ||
             (t->as<Tensor>().label() != reserved::antisymm_label() &&
              t->as<Tensor>().label() != reserved::symm_label()))) {
          const Tensor& tensor = t->as<Tensor>();
          ket_idx.insert(ket_idx.end(), tensor.ket().begin(),
                         tensor.ket().end());
        }
      }
      return ket_idx;
    };
    auto product_kets = get_ket_indices(temp_product);

    auto get_bra_indices = [&](const Product& prod) {
      container::svector<Index> bra_idx;
      for (auto&& t : prod) {
        if (t->is<Tensor>() &&
            (collect_symmetrizer_indices ||
             (t->as<Tensor>().label() != reserved::antisymm_label() &&
              t->as<Tensor>().label() != reserved::symm_label()))) {
          const Tensor& tensor = t->as<Tensor>();
          bra_idx.insert(bra_idx.end(), tensor.bra().begin(),
                         tensor.bra().end());
        }
      }
      return bra_idx;
    };
    auto product_bras = get_bra_indices(temp_product);

    auto substitute_ext_idx = [&product_bras,
                               &product_kets](const auto& idx_pair) {
      SEQUANT_ASSERT(idx_pair.size() == 2);
      const auto& what = get_bra_idx(idx_pair);
      const auto& with = get_ket_idx(idx_pair);
      std::replace(product_bras.begin(), product_bras.end(), what, with);
      std::replace(product_kets.begin(), product_kets.end(), what, with);
    };

    // Substitute indices from external index list
    ranges::for_each(ext_index_groups, substitute_ext_idx);

    auto n_cycles = count_cycles(product_kets, product_bras);

    auto result = std::make_shared<Product>(product);
    result->scale(pow2(n_cycles));
    return result;
  };

  if (expr->is<Constant>() || expr->is<Variable>())
    return expr;
  else if (expr->is<Tensor>())
    return trace_product(
        (ex<Constant>(1) * expr)->as<Product>());  // expand_all(expr);
  else if (expr->is<Product>())
    return trace_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      if (summand->is<Product>()) {
        result->append(trace_product(summand->as<Product>()));
      } else if (summand->is<Tensor>()) {
        result->append(
            trace_product((ex<Constant>(1) * summand)->as<Product>()));
      } else {
        SEQUANT_ASSERT(summand->is<Constant>() || summand->is<Variable>());
        result->append(summand);
      }
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in closed_shell_spintrace: " +
                    expr->type_name());
  }
}

ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<SlottedIndex>>&
        ext_index_groups,
    bool full_expansion) {
  return closed_shell_spintrace_impl(
      expression, as_view_of_index_groups(ext_index_groups), full_expansion);
}

ExprPtr closed_shell_spintrace(const ExprPtr& expression, EmptyInitializerList,
                               bool full_expansion) {
  return closed_shell_spintrace_impl(
      expression, container::svector<container::svector<Index>>{},
      full_expansion);
}

ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool full_expansion) {
  return closed_shell_spintrace_impl(expression, ext_index_groups,
                                     full_expansion);
}

container::svector<ResultExpr> closed_shell_spintrace(const ResultExpr& expr,
                                                      bool full_expansion) {
  using TraceFunction = ExprPtr (*)(
      const ExprPtr&,
      const container::svector<container::svector<SlottedIndex>>&, bool);

  return detail::wrap_trace<container::svector<ResultExpr>>(
      expr, static_cast<TraceFunction>(&closed_shell_spintrace),
      full_expansion);
}

ExprPtr closed_shell_CC_spintrace_v1(ExprPtr const& expr,
                                     ClosedShellCCSpintraceOptions options) {
  SEQUANT_ASSERT(options.method == BiorthogonalizationMethod::V1);
  using ranges::views::transform;

  auto const ext_idxs = external_indices(expr);
  auto st_expr = options.naive_spintrace
                     ? spintrace(expr, ext_idxs)
                     : closed_shell_spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Biorthogonal transformation without factoring out NNS projector
    st_expr = biorthogonal_transform_pre_nnsproject(st_expr, ext_idxs, false);
  }
  simplify(st_expr);

  // if (ext_idxs.size() <= 1){
  //   if (st_expr->is<Sum>())
  //   {
  //     for (auto& term : *st_expr) {
  //       if (term->is<Product>())
  //         term = remove_tensor(term.as_shared_ptr<Product>(),
  //                              reserved::symm_label());
  //     }
  //   }
  // }else if(ext_idxs.size() > 1) {
  //   st_expr = S_maps(st_expr);
  //   st_expr = ex<Constant>(2) * st_expr;
  //   simplify(st_expr);
  // }

  return st_expr;
}

ExprPtr closed_shell_CC_spintrace_v2(ExprPtr const& expr,
                                     ClosedShellCCSpintraceOptions options) {
  SEQUANT_ASSERT(options.method == BiorthogonalizationMethod::V2);
  using ranges::views::transform;
  auto const ext_idxs = external_indices(expr);
  auto st_expr = options.naive_spintrace
                     ? spintrace(expr, ext_idxs)
                     : closed_shell_spintrace(expr, ext_idxs);
  canonicalize(st_expr);

  if (!ext_idxs.empty()) {
    // Biorthogonal transformation with factoring out NNS projector
    st_expr = biorthogonal_transform_pre_nnsproject(st_expr, ext_idxs);
  }
  // if (ext_idxs.size() > 1) {
  //   st_expr = S_maps(st_expr);
  //   simplify(st_expr);
  // }

  simplify(st_expr);
  // std::wcout << "final eqs after symm: "
  //            << sequant::to_latex_align(
  //                   sequant::ex<sequant::Sum>(
  //                       sequant::opt::reorder(st_expr->as<sequant::Sum>())),
  //                   0, 4)
  //            << std::endl;

  return st_expr;
}

ExprPtr closed_shell_CC_spintrace(ExprPtr const& expr,
                                  ClosedShellCCSpintraceOptions options) {
  switch (options.method) {
    case BiorthogonalizationMethod::V1:
      return closed_shell_CC_spintrace_v1(expr, options);
    case BiorthogonalizationMethod::V2:
      return closed_shell_CC_spintrace_v2(expr, options);
    default:
      SEQUANT_ASSERT(false && "unreachable code reached");
      abort();
  }
}

Tensor swap_spin(const Tensor& t) {
  auto is_any_spin = [](const Index& i) {
    return mbpt::to_spin(i.space().qns()) == mbpt::Spin::any;
  };

  // Return tensor if there are no spin labels
  if (std::all_of(t.const_braket_indices().begin(),
                  t.const_braket_indices().end(), is_any_spin)) {
    return t;
  }

  // Return new index where the spin-label is flipped
  auto spin_flipped_idx = [](const Index& idx) {
    SEQUANT_ASSERT(mbpt::to_spin(idx.space().qns()) != mbpt::Spin::any);
    return mbpt::to_spin(idx.space().qns()) == mbpt::Spin::alpha
               ? make_spinbeta(idx)
               : make_spinalpha(idx);
  };

  container::svector<Index> b(t.rank()), k(t.rank());

  for (std::size_t i = 0; i < t.rank(); ++i) {
    b.at(i) = spin_flipped_idx(t.bra().at(i));
    k.at(i) = spin_flipped_idx(t.ket().at(i));
  }

  return {t.label(),    bra(std::move(b)),   ket(std::move(k)),  t.aux(),
          t.symmetry(), t.braket_symmetry(), t.column_symmetry()};
}

ExprPtr swap_spin(const ExprPtr& expr) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  auto swap_tensor = [](const Tensor& t) { return ex<Tensor>(swap_spin(t)); };

  auto swap_product = [&swap_tensor](const Product& p) {
    Product result{};
    result.scale(p.scalar());
    for (auto& t : p) {
      if (t->is<Tensor>()) {
        result.append(1, swap_tensor(t->as<Tensor>()), Product::Flatten::No);
      } else if (t->is<Constant>() || t->is<Variable>()) {
        result.append(1, t, Product::Flatten::No);
      } else {
        throw Exception("Invalid Expr type in swap_spin::swap_product: " +
                        t->type_name());
      }
    }
    return ex<Product>(result);
  };

  if (expr->is<Tensor>())
    return swap_tensor(expr->as<Tensor>());
  else if (expr->is<Product>())
    return swap_product(expr->as<Product>());
  else if (expr->is<Sum>()) {
    Sum result;
    for (auto& term : *expr) {
      result.append(swap_spin(term));
    }
    return ex<Sum>(result);
  } else {
    throw Exception("Invalid Expr type in swap_spin: " + expr->type_name());
  }
}

ExprPtr merge_tensors(const Tensor& O1, const Tensor& O2) {
  SEQUANT_ASSERT(O1.label() == O2.label());
  SEQUANT_ASSERT(O1.symmetry() == O2.symmetry());
  auto b = ranges::views::concat(O1.bra(), O2.bra());
  auto k = ranges::views::concat(O1.ket(), O2.ket());
  auto a = ranges::views::concat(O1.aux(), O2.aux());
  return ex<Tensor>(Tensor(O1.label(), bra(b), ket(k), aux(a), O1.symmetry()));
}

template <detail::index_group_range IdxGroups>
ExprPtr spintrace_impl(
    const ExprPtr& expression, IdxGroups&& ext_index_groups,
    bool spinfree_index_spaces,
    const std::function<bool(uint64_t)>& spincase_predicate = {});

namespace {

std::size_t open_shell_n_spin_cases(std::size_t n_groups, bool npc,
                                    std::size_t n_paired,
                                    bool all_external_spin_assignments) {
  if (all_external_spin_assignments) return std::size_t{1} << n_groups;
  const auto n_unpaired = n_groups - n_paired;
  return npc ? (n_paired + 1) * (n_unpaired + 1) : (n_groups + 1);
}

void open_shell_fill_group_spins(std::size_t sc, std::size_t n_groups, bool npc,
                                 std::size_t n_paired,
                                 bool all_external_spin_assignments,
                                 container::svector<int>& group_spins) {
  group_spins.assign(n_groups, 0);
  if (all_external_spin_assignments) {
    for (std::size_t g = 0; g < n_groups; ++g)
      group_spins[g] = static_cast<int>((sc >> g) & 1u);
    return;
  }
  if (npc) {
    const auto n_paired_beta = sc % (n_paired + 1);
    const auto n_unpaired_beta = sc / (n_paired + 1);
    if (n_paired_beta > 0)
      std::fill(group_spins.begin() + (n_paired - n_paired_beta),
                group_spins.begin() + n_paired, 1);
    if (n_unpaired_beta > 0)
      std::fill(group_spins.begin() + (n_groups - n_unpaired_beta),
                group_spins.end(), 1);
  } else {
    std::fill(group_spins.end() - static_cast<std::size_t>(sc),
              group_spins.end(), 1);
  }
}

}  // namespace

std::vector<ExprPtr> open_shell_A_op(const Tensor& A,
                                     bool all_external_spin_assignments) {
  SEQUANT_ASSERT(A.label() == reserved::antisymm_label());
  const auto n_bra = A.bra_rank();
  const auto n_ket = A.ket_rank();
  const auto n_paired = std::min(n_bra, n_ket);
  const auto n_groups = std::max(n_bra, n_ket);
  const bool npc = (n_bra != n_ket);
  const auto n_spin_cases = open_shell_n_spin_cases(
      n_groups, npc, n_paired, all_external_spin_assignments);

  std::vector<ExprPtr> result(n_spin_cases);

  // first (all alpha) and last (all beta) are always trivial
  result[0] = ex<Constant>(1);
  result[n_spin_cases - 1] = ex<Constant>(1);

  for (std::size_t sc = 1; sc < n_spin_cases - 1; ++sc) {
    // per-group spin assignment: 0 = alpha, 1 = beta
    // groups [0, n_paired) are paired (bra+ket), [n_paired, n_groups) are
    // unpaired
    container::svector<int> group_spin;
    open_shell_fill_group_spins(sc, n_groups, npc, n_paired,
                                all_external_spin_assignments, group_spin);

    // assign spin to bra and ket indices according to their group
    auto assign_spin = [&](auto& indices, std::size_t rank,
                           const auto& originals) {
      for (std::size_t g = 0; g < rank; ++g) {
        indices[g] = group_spin[g] == 0 ? make_spinalpha(originals[g])
                                        : make_spinbeta(originals[g]);
      }
    };

    auto spin_bra = A.bra();
    auto spin_ket = A.ket();
    assign_spin(spin_bra, n_bra, A.bra());
    assign_spin(spin_ket, n_ket, A.ket());

    ranges::for_each(spin_bra, [](const Index& i) { i.reset_tag(); });
    ranges::for_each(spin_ket, [](const Index& i) { i.reset_tag(); });

    result[sc] = ex<Tensor>(Tensor(reserved::antisymm_label(), spin_bra,
                                   spin_ket, A.aux(), Symmetry::Antisymm));
  }

  return result;
}

std::vector<ExprPtr> open_shell_P_op_vector(
    const Tensor& A, bool all_external_spin_assignments) {
  SEQUANT_ASSERT(A.label() == reserved::antisymm_label());
  const auto n_bra = A.bra_rank();
  const auto n_ket = A.ket_rank();
  const auto n_groups = std::max(n_bra, n_ket);
  const auto n_paired = std::min(n_bra, n_ket);
  const bool nonparticle_conserving = (n_bra != n_ket);
  const auto n_spin_cases =
      open_shell_n_spin_cases(n_groups, nonparticle_conserving, n_paired,
                              all_external_spin_assignments);

  std::vector<ExprPtr> result_vector(n_spin_cases);

  for (size_t sc = 0; sc < n_spin_cases; ++sc) {
    container::svector<int> group_spins;
    open_shell_fill_group_spins(sc, n_groups, nonparticle_conserving, n_paired,
                                all_external_spin_assignments, group_spins);

    container::svector<int> alpha_bra_indices, beta_bra_indices;
    for (size_t i = 0; i < n_bra; ++i) {
      size_t g = (i < n_paired) ? i : n_paired + (i - n_paired);
      if (group_spins[g] == 0)
        alpha_bra_indices.push_back(static_cast<int>(i));
      else
        beta_bra_indices.push_back(static_cast<int>(i));
    }

    container::svector<int> alpha_ket_indices, beta_ket_indices;
    for (size_t i = 0; i < n_ket; ++i) {
      size_t g = (i < n_paired) ? i : n_paired + (i - n_paired);
      if (group_spins[g] == 0)
        alpha_ket_indices.push_back(static_cast<int>(i));
      else
        beta_ket_indices.push_back(static_cast<int>(i));
    }

    bool bra_uniform = alpha_bra_indices.empty() || beta_bra_indices.empty();
    bool ket_uniform = alpha_ket_indices.empty() || beta_ket_indices.empty();

    if (bra_uniform && ket_uniform) {
      result_vector[sc] = ex<Constant>(1);
      continue;
    }

    container::svector<Tensor> P_bra_list, P_ket_list;

    for (auto& j : alpha_bra_indices) {
      for (auto& k : beta_bra_indices) {
        if (!alpha_bra_indices.empty() && !beta_bra_indices.empty()) {
          P_bra_list.emplace_back(Tensor(reserved::transposition_label(),
                                         bra{A.bra().at(j), A.bra().at(k)},
                                         ket{}, Symmetry::Symm));
        }
      }
    }
    for (auto& j : alpha_ket_indices) {
      for (auto& k : beta_ket_indices) {
        if (!alpha_ket_indices.empty() && !beta_ket_indices.empty()) {
          P_ket_list.emplace_back(Tensor(reserved::transposition_label(), bra{},
                                         ket{A.ket().at(j), A.ket().at(k)},
                                         Symmetry::Symm));
        }
      }
    }
    // P4 terms
    if (alpha_bra_indices.size() > 1 && beta_bra_indices.size() > 1) {
      for (std::size_t a = 0; a != alpha_bra_indices.size() - 1; ++a) {
        auto i1 = alpha_bra_indices[a];
        for (std::size_t b = a + 1; b != alpha_bra_indices.size(); ++b) {
          auto i2 = alpha_bra_indices[b];
          for (std::size_t c = 0; c != beta_bra_indices.size() - 1; ++c) {
            auto i3 = beta_bra_indices[c];
            for (std::size_t d = c + 1; d != beta_bra_indices.size(); ++d) {
              auto i4 = beta_bra_indices[d];
              P_bra_list.emplace_back(
                  Tensor(reserved::transposition_label(),
                         bra{A.bra().at(i1), A.bra().at(i3), A.bra().at(i2),
                             A.bra().at(i4)},
                         ket{}, Symmetry::Symm));
            }
          }
        }
      }
    }
    if (alpha_ket_indices.size() > 1 && beta_ket_indices.size() > 1) {
      for (std::size_t a = 0; a != alpha_ket_indices.size() - 1; ++a) {
        auto i1 = alpha_ket_indices[a];
        for (std::size_t b = a + 1; b != alpha_ket_indices.size(); ++b) {
          auto i2 = alpha_ket_indices[b];
          for (std::size_t c = 0; c != beta_ket_indices.size() - 1; ++c) {
            auto i3 = beta_ket_indices[c];
            for (std::size_t d = c + 1; d != beta_ket_indices.size(); ++d) {
              auto i4 = beta_ket_indices[d];
              P_ket_list.emplace_back(
                  Tensor(reserved::transposition_label(), bra{},
                         ket{A.ket().at(i1), A.ket().at(i3), A.ket().at(i2),
                             A.ket().at(i4)},
                         Symmetry::Symm));
            }
          }
        }
      }
    }

    Sum bra_permutations{};
    bra_permutations.append(ex<Constant>(1));
    for (auto& p : P_bra_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      bra_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    Sum ket_permutations{};
    ket_permutations.append(ex<Constant>(1));
    for (auto& p : P_ket_list) {
      int prefactor = (p.bra_rank() + p.ket_rank() == 4) ? 1 : -1;
      ket_permutations.append(ex<Constant>(prefactor) * ex<Tensor>(p));
    }

    ExprPtr spin_case_result =
        ex<Sum>(bra_permutations) * ex<Sum>(ket_permutations) /
        (bra_permutations.size() * ket_permutations.size());
    expand(spin_case_result);

    // Merge P operators if it encounters alpha_spin product of operators
    for (auto& term : *spin_case_result) {
      if (term->is<Product>()) {
        const auto& P = term->as<Product>();
        const auto nfactors = P.factors().size();
        SEQUANT_ASSERT(
            nfactors <=
            2);  // constant, single P, or 2 P's (one for alpha, one for beta)
        if (nfactors == 2) {
          const auto scalar = P.scalar();
          auto P1 = P.factor(0)->as<Tensor>();
          auto P2 = P.factor(1)->as<Tensor>();
          term = ex<Constant>(scalar) * merge_tensors(P1, P2);
        }
      }
    }
    result_vector[sc] = spin_case_result;
  }
  return result_vector;
}

template <detail::index_group_range IdxGroups>
std::vector<ExprPtr> open_shell_spintrace_impl(
    const ExprPtr& expr, IdxGroups&& ext_index_groups,
    const std::optional<int>& target_spin_case) {
  if (expr->is<Constant>() || expr->is<Variable>()) {
    return std::vector<ExprPtr>{expr};
  }

  // Grand index list contains both internal and external indices
  container::set<Index, Index::LabelCompare> grand_idxlist =
      get_used_indices<decltype(grand_idxlist)>(expr);

  container::set<Index> ext_idxlist;
  for (const auto& idxgrp : ext_index_groups) {
    for (const Index& current : idxgrp) {
      Index idx = current;
      idx.reset_tag();
      ext_idxlist.insert(std::move(idx));
    }
  }

  container::set<Index> int_idxlist;
  for (auto&& gidx : grand_idxlist) {
    if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
      int_idxlist.insert(gidx);
    }
  }

  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> int_index_groups;
  for (auto&& i : int_idxlist) {
    int_index_groups.emplace_back(IndexGroup(1, i));
  }

  SEQUANT_ASSERT(grand_idxlist.size() ==
                 int_idxlist.size() + ext_idxlist.size());

  // paired groups have size 2 (bra + ket), unpaired have size 1
  const std::size_t n_ext_groups = ext_index_groups.size();
  const std::size_t n_paired = ranges::count_if(
      ext_index_groups, [](const auto& grp) { return grp.size() == 2; });
  const std::size_t n_unpaired = n_ext_groups - n_paired;
  const bool nonparticle_conserving = n_unpaired > 0;

  // make a spin-specific index, orientation is given by spin_bit: 0 =
  // spin-down/beta, 1 = spin-up/alpha
  auto make_spinspecific = [](const Index& idx, const long int& spin_bit) {
    return spin_bit == 0 ? make_spinalpha(idx) : make_spinbeta(idx);
  };

  // Generate index replacement maps
  auto spin_cases = [&make_spinspecific](
                        const container::svector<IndexGroup>& idx_group) {
    const auto ncases = pow2(idx_group.size());
    container::svector<container::map<Index, Index>> all_replacements(ncases);

    for (uint64_t i = 0; i != ncases; ++i) {
      container::map<Index, Index> idx_rep;
      for (size_t idxg = 0; idxg != idx_group.size(); ++idxg) {
        auto spin_bit = (i << (64 - idxg - 1)) >> 63;
        SEQUANT_ASSERT((spin_bit == 0) || (spin_bit == 1));
        for (auto& idx : idx_group[idxg]) {
          auto spin_idx = make_spinspecific(idx, spin_bit);
          idx_rep.emplace(idx, spin_idx);
        }
      }
      all_replacements[i] = idx_rep;
    }
    return all_replacements;
  };

  // External index replacement maps
  auto ext_spin_cases = [&make_spinspecific, nonparticle_conserving, n_paired,
                         n_ext_groups](const auto& idx_groups) {
    container::svector<container::map<Index, Index>> all_replacements;

    const uint64_t ncases =
        nonparticle_conserving
            ? ((n_paired + 1) * (n_ext_groups - n_paired + 1))
            : (n_ext_groups + 1);

    for (uint64_t i = 0; i != ncases; ++i) {
      container::svector<int> spins(idx_groups.size(), 0);

      if (nonparticle_conserving) {
        const std::size_t n_paired_beta = i % (n_paired + 1);
        const std::size_t n_unpaired_beta = i / (n_paired + 1);
        if (n_paired_beta > 0)
          std::fill(spins.begin() + (n_paired - n_paired_beta),
                    spins.begin() + n_paired, 1);
        if (n_unpaired_beta > 0)
          std::fill(spins.begin() + (n_ext_groups - n_unpaired_beta),
                    spins.end(), 1);
      } else {
        std::fill(spins.end() - static_cast<std::size_t>(i), spins.end(), 1);
      }

      container::map<Index, Index> idx_rep;
      for (std::size_t j = 0; j != idx_groups.size(); ++j) {
        for (const Index& idx : idx_groups[j]) {
          auto spin_idx = make_spinspecific(idx, spins[j]);
          idx_rep.emplace(idx, spin_idx);
        }
      }
      all_replacements.push_back(idx_rep);
    }
    return all_replacements;
  };

  // Internal and external index replacements are independent
  auto i_rep = spin_cases(int_index_groups);
  auto e_rep = ext_spin_cases(ext_index_groups);

  // For a single spin case, keep only the relevant spin case
  // PS: all alpha indexing start at 0
  if (target_spin_case) {
    auto external_replacement_map = e_rep.at(*target_spin_case);
    e_rep.clear();
    e_rep.push_back(external_replacement_map);
  }

  // Expand 'A' operator and 'antisymm' tensors
  auto expanded_expr = expand_A_op(expr);
  detail::reset_idx_tags(expanded_expr);
  expand(expanded_expr);
  simplify(expanded_expr);

  std::vector<ExprPtr> result{};

  // return true if a product is spin-symmetric
  // - for purely particle-conserving products (all tensors have bra_rank ==
  //   ket_rank), use the original global bra/ket check.
  // - if any non-particle-conserving tensor is present, fall back to
  //   per-tensor ms_conserving_columns check

  auto spin_symm_product = [](const Product& product) {
    bool has_nonpc_tensor = false;
    for (const auto& term : product) {
      if (!term->is<Tensor>()) continue;
      const auto& tn = term->as<Tensor>();
      if (tn.bra_rank() != tn.ket_rank()) {
        has_nonpc_tensor = true;
        break;
      }
    }

    if (!has_nonpc_tensor) {
      container::svector<Index> cBra, cKet;  // concat Bra and concat Ket
      for (auto& term : product) {
        if (term->is<Tensor>()) {
          const auto& tnsr = term->as<Tensor>();
          cBra.insert(cBra.end(), tnsr.bra().begin(), tnsr.bra().end());
          cKet.insert(cKet.end(), tnsr.ket().begin(), tnsr.ket().end());
        } else if (term->is<Product>() || term->is<Sum>()) {
          throw Exception(
              "Nested Product and Sum not supported in spin_symm_product");
        }
      }
      SEQUANT_ASSERT(cKet.size() == cBra.size());

      auto i_ket = cKet.begin();
      for (auto& b : cBra) {
        if (b.space().qns() != i_ket->space().qns()) return false;
        ++i_ket;
      }
      return true;
    }

    // NPC fallback: check each PC tensor column-wise
    for (const auto& term : product) {
      if (!term->is<Tensor>()) continue;
      const auto& tn = term->as<Tensor>();
      if (tn.bra_rank() != tn.ket_rank()) continue;
      if (!ms_conserving_columns(tn)) return false;
    }
    return true;
  };

  //
  // SPIN-TRACING algorithm begins here
  //

  // Loop over external index replacement maps
  for (auto& e : e_rep) {
    // Add spin labels to external indices
    auto spin_expr = append_spin(expanded_expr, e);
    detail::reset_idx_tags(spin_expr);
    Sum e_result{};

    // Loop over internal index replacement maps
    for (auto& i : i_rep) {
      // Add spin labels to internal indices, expand antisymmetric tensors
      ExprPtr spin_expr_i = append_spin(spin_expr, i);
      spin_expr_i = expand_antisymm(spin_expr_i, true);
      // std::wcout << "after expand_antisymm: " << to_latex_align(spin_expr_i,
      // 0, 4) << '\n';
      expand(spin_expr_i);
      detail::reset_idx_tags(spin_expr_i);
      Sum i_result{};

      if (spin_expr_i->is<Tensor>() || spin_expr_i->is<Constant>() ||
          spin_expr_i->is<Variable>()) {
        e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Product>()) {
        if (spin_symm_product(spin_expr_i->as<Product>()))
          e_result.append(spin_expr_i);
      } else if (spin_expr_i->is<Sum>()) {
        for (auto& pr : *spin_expr_i) {
          if (pr->is<Product>()) {
            if (spin_symm_product(pr->as<Product>())) i_result.append(pr);
          } else if (pr->is<Tensor>()) {
            if (ms_conserving_columns(pr->as<Tensor>())) i_result.append(pr);
          } else if (pr->is<Constant>() || pr->is<Variable>()) {
            i_result.append(pr);
          } else
            throw("Unknown ExprPtr type.");
        }
        e_result.append(std::make_shared<Sum>(i_result));
      }
    }  // loop over internal indices
    // std::wcout << "e_result size after internal loop: " <<
    // e_result.summands().size() << "\n";
    result.push_back(std::make_shared<Sum>(e_result));
  }  // loop over external indices

  if (target_spin_case) {
    SEQUANT_ASSERT(result.size() == 1 &&
                   "Spin-specific case must return one expression.");
  }
  for (auto& expression : result) {
    detail::reset_idx_tags(expression);
    // if (target_spin_case && *target_spin_case == 1) {
    //   std::wcout << "before canon (sc=1): " << to_latex_align(expression)
    //              << "\n";
    // }
    canonicalize(expression);
    // if (target_spin_case && *target_spin_case == 1) {
    //   std::wcout << "after canon (sc=1): " << to_latex_align(expression)
    //              << "\n";
    // }
    rapid_simplify(expression);
  }

  return result;
}

std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<SlottedIndex>>&
        ext_index_groups,
    const std::optional<int>& target_spin_case) {
  return open_shell_spintrace_impl(
      expr, as_view_of_index_groups(ext_index_groups), target_spin_case);
}

std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr, EmptyInitializerList,
    const std::optional<int>& target_spin_case) {
  return open_shell_spintrace_impl(
      expr, container::svector<container::svector<Index>>{}, target_spin_case);
}

std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    const std::optional<int>& target_spin_case) {
  return open_shell_spintrace_impl(expr, ext_index_groups, target_spin_case);
}

std::vector<ExprPtr> open_shell_CC_spintrace(const ExprPtr& expr) {
  SEQUANT_ASSERT(expr->is<Sum>() || expr->is<Product>());
  // Pop the antisymmetrizer A off a copy of the leading term to detect and
  // retrieve it. Energy-like expressions (e.g. CC energy) have none.
  // The input is assumed homogeneous: either every term carries a leading
  // antisymmetrizer (residual) or none does (energy), so inspecting the
  // leading term alone suffices to classify the whole expression.
  ExprPtr leading_term =
      (expr->is<Sum>() ? expr->as<Sum>().summand(0) : expr).clone();
  const auto A_opt = pop_tensor(leading_term, reserved::antisymm_label());

  // Expressions without a leading antisymmetrizer (e.g. CC energy) collapse
  // to a single expression. They are spin-traced directly, but term-by-term.
  if (!A_opt) {
    auto trace_term = [](const ExprPtr& term) -> ExprPtr {
      if (term->is<Constant>() || term->is<Variable>()) return term;
      return open_shell_spintrace(term, {}).at(0);
    };
    const auto terms =
        expr.is<Sum>() ? expr->as<Sum>().summands() : Sum::summands_type{expr};
    auto result = ex<Sum>(terms | ranges::views::transform(trace_term));
    simplify(result);
    return {result};
  }

  const Tensor& A = A_opt.value()->as<Tensor>();
  const size_t n_bra = A.bra_rank();
  const size_t n_ket = A.ket_rank();
  const size_t n_groups = std::max(n_bra, n_ket);
  const auto n_paired = std::min(n_bra, n_ket);
  const auto n_unpaired = n_groups - n_paired;
  const bool nonparticle_conserving = (n_bra != n_ket);
  const size_t n_spin_cases = nonparticle_conserving
                                  ? (n_paired + 1) * (n_unpaired + 1)
                                  : (n_groups + 1);

  auto P_vec = open_shell_P_op_vector(A);
  auto A_vec = open_shell_A_op(A);
  SEQUANT_ASSERT(P_vec.size() == n_spin_cases);
  SEQUANT_ASSERT(A_vec.size() == n_spin_cases);

  const auto ext_groups = external_indices(A);

  std::vector<Sum> concat_terms(n_spin_cases);
  [[maybe_unused]] size_t n_spin_orbital_term = 0;
  for (auto& product_term : expr.is<Sum>()
                                ? std::span{expr.as<Sum>().summands()}
                                : std::span{&expr, 1}) {
    auto term = remove_tensor(product_term.as_shared_ptr<Product>(),
                              reserved::antisymm_label());
    std::vector<ExprPtr> os_st(n_spin_cases);

    for (std::size_t s = 0; s != os_st.size(); ++s) {
      os_st.at(s) = P_vec.at(s) * term;
      expand(os_st.at(s));
      os_st.at(s) = expand_P_op(os_st.at(s));
      os_st.at(s) =
          open_shell_spintrace(os_st.at(s), ext_groups, static_cast<int>(s))
              .at(0);

      bool need_A =
          (std::max(n_bra, n_ket) > 2) && !A_vec.at(s)->is<Constant>();
      if (need_A) {
        os_st.at(s) = A_vec.at(s) * os_st.at(s);
        simplify(os_st.at(s));
        os_st.at(s) = remove_tensor(os_st.at(s), reserved::antisymm_label());
      }
    }

    for (size_t j = 0; j != os_st.size(); ++j) {
      concat_terms.at(j).append(os_st.at(j));
    }
    ++n_spin_orbital_term;
  }

  std::vector<ExprPtr> expr_vec;
  for (auto& spin_case : concat_terms) {
    auto ptr = sequant::ex<Sum>(spin_case);
    expr_vec.push_back(ptr);
  }

  return expr_vec;
}

template <detail::index_group_range IdxGroups>
ExprPtr spintrace_impl(
    const ExprPtr& expression, IdxGroups&& ext_index_groups,
    bool spinfree_index_spaces,
    const std::function<bool(uint64_t)>& spincase_predicate) {
  // Escape immediately if expression is a constant
  if (expression->is<Constant>() || expression->is<Variable>()) {
    return expression;
  }

  if constexpr (assert_enabled()) {
    // Verify that the number of external indices matches the number of indices
    // in ext_index_groups, UNLESS user overrode external definitions in default
    // context
    const auto& copts = get_default_context().canonicalization_options();
    if (!copts.has_value() || !copts->named_indices) {
      auto count_indices = [](const auto& range) {
        auto sizes = range | ranges::views::transform(
                                 [](const auto& list) { return list.size(); });
        return std::accumulate(sizes.begin(), sizes.end(), 0);
      };
      auto determined_externals = external_indices(expression);

      SEQUANT_ASSERT(count_indices(ext_index_groups) ==
                     count_indices(determined_externals));
    }
  }

  // This function must be used for tensors with spin-specific indices only. If
  // the spin-symmetry is conserved: the tensor is expanded; else: zero is
  // returned.
  auto spintrace_tensor = [](const Tensor& tensor) {
    return can_expand(tensor) ? expand_antisymm(tensor) : ex<Constant>(0);
  };

  // This function is used to spin-trace a product terms with spin-specific
  // indices. It checks if all tensors can be expanded and spintraces individual
  // tensors by call to the spin_trace_tensor lambda.
  auto spintrace_product =
      [&spintrace_tensor](const ProductPtr& product) -> ExprPtr {
    ProductPtr spin_product = std::make_shared<Product>();

    // Check if all tensors in this product can be expanded
    // If NOT all tensors can be expanded, return zero
    for (const ExprPtr& expr : product->factors()) {
      if (expr.is<Tensor>()) {
        if (!can_expand(expr.as<Tensor>())) {
          return ex<Constant>(0);
        }
      } else if (expr.is<Sum>() || expr.is<Product>()) {
        throw Exception(
            "Nested sums/products not supported in spin_trace_product");
      }
    }

    spin_product->scale(product->scalar());
    for (const ExprPtr& expr : product->factors()) {
      if (expr.is<Tensor>()) {
        spin_product->append(1, spintrace_tensor(expr.as<Tensor>()));
      } else if (expr.is<Variable>() || expr.is<Constant>()) {
        spin_product->append(1, expr.clone());
      } else {
        // Would need some sort of recursion but it is not clear how that would
        // interact with other code in here yet so prefer to error instead.
        throw Exception(
            "spin_trace_product: Nested products or sums inside of a Product "
            "not supported (yet)");
      }
    }

    ExprPtr result = spin_product;
    expand(result);
    rapid_simplify(result);
    return result;
  };

  // Most important lambda of this function
  auto trace_product = [&ext_index_groups, &spintrace_tensor,
                        &spintrace_product, spinfree_index_spaces,
                        spincase_predicate](const ProductPtr& product) {
    ExprPtr expr = product->clone();
    // List of all indices in the expression
    container::set<Index, Index::LabelCompare> grand_idxlist =
        get_used_indices<decltype(grand_idxlist)>(expr);

    // List of external indices, i.e. indices that are not summed over Einstein
    // style (indices that are not repeated in an expression)
    container::set<Index> ext_idxlist;
    for (const auto& idxgrp : ext_index_groups) {
      for (const Index& current : idxgrp) {
        Index idx = current;
        idx.reset_tag();
        ext_idxlist.insert(std::move(idx));
      }
    }

    // List of internal indices, i.e. indices that are contracted over
    container::set<Index> int_idxlist;
    for (auto&& gidx : grand_idxlist) {
      if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
        int_idxlist.insert(gidx);
      }
    }

    // EFV: generate the grand list of index groups by concatenating list of
    // external index groups with the groups of internal indices (each
    // internal index = 1 group)
    // TODO some internal indices can be a priori placed in the same group, if
    // they refer to the same particle of a spin-free non-antisymmetrized Tensor
    //      so visit all Tensors in the expression and locate such groups of
    //      internal indices before placing the rest into separate groups
    using IndexGroup = container::svector<Index>;
    container::svector<IndexGroup> index_groups;
    // Externals first so external group g lives at bit g of spincase_bitstr;
    // this lets spincase_predicate locate external bits without knowing the
    // per-product internal count.
    for (const auto& group : ext_index_groups) {
      IndexGroup target;
      target.reserve(group.size());
      for (const Index& current : group) {
        target.emplace_back(current);
      }
      index_groups.emplace_back(std::move(target));
    }
    for (auto&& i : int_idxlist) index_groups.emplace_back(IndexGroup(1, i));

    // EFV: for each spincase (loop over integer from 0 to 2^n-1, n=#of index
    // groups)
    SEQUANT_ASSERT(index_groups.size() <= 64);
    const uint64_t nspincases = pow2(index_groups.size());

    auto result = std::make_shared<Sum>();
    for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
         ++spincase_bitstr) {
      if (spincase_predicate && !spincase_predicate(spincase_bitstr)) {
        continue;
      }
      // EFV:  assign spin to each index group => make a replacement list
      container::map<Index, Index> index_replacements;

      uint64_t index_group_count = 0;
      for (auto&& index_group : index_groups) {
        auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
        SEQUANT_ASSERT(spin_bit == 0 || spin_bit == 1);

        for (auto&& index : index_group) {
          index_replacements.emplace(index, spin_bit == 0
                                                ? make_spinalpha(index)
                                                : make_spinbeta(index));
        }
        ++index_group_count;
      }

      // Append spin labels to indices in the expression
      auto spin_expr = append_spin(expr, index_replacements);
      rapid_simplify(spin_expr);  // This call is required for Tensor case
      detail::reset_idx_tags(spin_expr);

      // NB: There are temporaries in the following code to enable
      // printing intermediate expressions.
      if (spin_expr->is<Tensor>()) {
        auto st_expr = spintrace_tensor(spin_expr->as<Tensor>());
        result->append(spinfree_index_spaces ? remove_spin(st_expr) : st_expr);
      } else if (spin_expr->is<Product>()) {
        ExprPtr st_expr = spintrace_product(spin_expr.as_shared_ptr<Product>());
        if (!st_expr->is<Constant>() || st_expr->as<Constant>().value() != 0) {
          result->append(spinfree_index_spaces ? remove_spin(st_expr)
                                               : st_expr);
        }
      } else if (spin_expr->is<Sum>()) {
        for (auto&& summand : *spin_expr) {
          Sum st_expr{};
          if (summand->is<Tensor>())
            st_expr.append(spintrace_tensor(summand->as<Tensor>()));
          else if (summand->is<Product>())
            st_expr.append(spintrace_product(summand.as_shared_ptr<Product>()));
          else {
            st_expr.append(summand);
          }
          auto st_expr_ptr = ex<Sum>(st_expr);
          result->append(spinfree_index_spaces ? remove_spin(st_expr_ptr)
                                               : st_expr_ptr);
        }
      } else {
        result->append(expr);
      }
    }  // Permutation FOR loop
    return result;
  };

  // Expand antisymmetrizer operator (A) if present in the expression
  ExprPtr expr = expression;
  if (has_tensor(expr, reserved::antisymm_label())) expr = expand_A_op(expr);

  if (expr->is<Tensor>()) expr = ex<Constant>(1) * expr;

  ExprPtr result;
  if (expr->is<Product>()) {
    result = trace_product(expr.as_shared_ptr<Product>());
  } else if (expr->is<Sum>()) {
    auto result_sum = std::make_shared<Sum>();
    for (auto&& term : *expr) {
      if (term->is<Product>())
        result_sum->append(trace_product(term.as_shared_ptr<Product>()));
      else if (term->is<Tensor>()) {
        auto term_as_product = ex<Constant>(1) * term;
        result_sum->append(
            trace_product(term_as_product.as_shared_ptr<Product>()));
      } else
        result_sum->append(term);
      result = result_sum;
    }
    return result;
  } else {
    throw Exception("Invalid Expr type in spintrace: " + expr->type_name());
  }

  detail::reset_idx_tags(result);
  return result;
}

ExprPtr spintrace(const ExprPtr& expression,
                  const container::svector<container::svector<SlottedIndex>>&
                      ext_index_groups,
                  bool spinfree_index_spaces) {
  return spintrace_impl(expression, as_view_of_index_groups(ext_index_groups),
                        spinfree_index_spaces);
}

ExprPtr spintrace(const ExprPtr& expression, EmptyInitializerList,
                  bool spinfree_index_spaces) {
  return spintrace_impl(expression,
                        container::svector<container::svector<Index>>{},
                        spinfree_index_spaces);
}

ExprPtr spintrace(
    const ExprPtr& expression,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool spinfree_index_spaces) {
  return spintrace_impl(expression, ext_index_groups, spinfree_index_spaces);
}

container::svector<ResultExpr> spintrace(const ResultExpr& expr,
                                         bool spinfree_index_spaces) {
  using TraceFunction = ExprPtr (*)(
      const ExprPtr&,
      const container::svector<container::svector<SlottedIndex>>&, bool);

  return detail::wrap_trace<container::svector<ResultExpr>>(
      expr, static_cast<TraceFunction>(&spintrace), spinfree_index_spaces);
}

namespace {

ExprPtr triplet_adapt_amplitudes(const ExprPtr& spin_labeled,
                                 bool amp_no_swap = false) {
  auto adapt_product = [amp_no_swap](const Product& p) -> ExprPtr {
    // sign of the M_S = 0 triplet coupling for a spin-labeled line:
    // T_{pq} = a_{pq}(alpha) - a_{pq}(beta)
    auto spin_line_sign = [](const Index& idx) {
      const auto s = mbpt::to_spin(idx.space().qns());
      SEQUANT_ASSERT(s == mbpt::Spin::alpha || s == mbpt::Spin::beta);
      return s == mbpt::Spin::alpha ? 1 : -1;
    };

    auto out = std::make_shared<Product>();
    out->scale(p.scalar());

    for (const auto& f : p) {
      if (!f->is<Tensor>() || f->as<Tensor>().label() != L"R") {
        SEQUANT_ASSERT(
            (!f->is<Tensor>() || f->as<Tensor>().label() != L"L") &&
            "triplet spin adaptation supports right eigenvectors (R) only");
        out->append(1, f, Product::Flatten::No);
        continue;
      }

      const auto& t = f->as<Tensor>();
      SEQUANT_ASSERT(t.bra_rank() == t.ket_rank() &&
                     "triplet spin adaptation supports particle-conserving "
                     "amplitudes only");
      const auto rank = t.bra_rank();

      if (rank == 1) {
        out->scale(spin_line_sign(t.bra().at(0)));
        out->append(1, f, Product::Flatten::No);
      } else if (rank == 2) {
        // antisymmetric spin-orbital tensors were already expanded to
        // Nonsymm representatives by expand_antisymm inside spintrace
        SEQUANT_ASSERT(t.symmetry() == Symmetry::Nonsymm);
        const int s0 = spin_line_sign(t.bra().at(0));
        const int s1 = spin_line_sign(t.bra().at(1));
        const bool same_spin = (s0 == s1);

        // Production (R + R_swap): ColumnSymmetry::Nonsymm keeps R_{a1i1;a2i2}
        // and R_{a2i2;a1i1} distinct. ter_only (amp_no_swap): single R only, so
        // use default ColumnSymmetry::Symm and let canonicalization fix order.
        const auto col_symm =
            amp_no_swap ? ColumnSymmetry::Nonsymm : ColumnSymmetry::Nonsymm;
        Tensor stored(t.label(), bra(t.bra()), ket(t.ket()), t.aux(),
                      Symmetry::Nonsymm, t.braket_symmetry(), col_symm);

        // EFV test: just TE (not TE_swap), emit a single R amplitude (no
        // R_swap),
        if (amp_no_swap) {
          out->scale(s0);
          out->append(1, ex<Tensor>(std::move(stored)), Product::Flatten::No);
          continue;
        }

        // swap columns
        container::svector<Index> swapped_bra{t.bra().at(1), t.bra().at(0)};
        container::svector<Index> swapped_ket{t.ket().at(1), t.ket().at(0)};
        Tensor swapped(t.label(), bra(std::move(swapped_bra)),
                       ket(std::move(swapped_ket)), t.aux(), Symmetry::Nonsymm,
                       t.braket_symmetry(), ColumnSymmetry::Nonsymm);

        auto channel = std::make_shared<Sum>();
        channel->append(ex<Constant>(s0) * ex<Tensor>(std::move(stored)));
        channel->append(ex<Constant>(same_spin ? s0 : -s0) *
                        ex<Tensor>(std::move(swapped)));
        out->append(1, ExprPtr(channel), Product::Flatten::No);
      } else if (rank == 3) {
        SEQUANT_ASSERT(t.symmetry() == Symmetry::Nonsymm);
        if (amp_no_swap)
          throw Exception(
              "triplet_adapt_amplitudes: ter_only/amp_no_swap is a "
              "doubles-only experiment, not implemented for triples");
        // R3 = sum over labels of r * T(c0) E(c1) E(c2) (T rides column 0).
        // In a fixed external spin sector the spin-orbital amplitude collects
        // all 6 column-pair permutations of R, each weighted by the spin sign
        // of the column that lands in R's first (T-carrying) slot — the rank-3
        // generalization of s0*R + s1*R_swap above.
        // I might be able to use only 3 column-pairs to make the eqns more
        // compact. Check it later.
        constexpr std::array<std::array<int, 3>, 6> col_perms{
            {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
        auto channel = std::make_shared<Sum>();
        for (const auto& p : col_perms) {
          container::svector<Index> pbra{t.bra().at(p[0]), t.bra().at(p[1]),
                                         t.bra().at(p[2])};
          container::svector<Index> pket{t.ket().at(p[0]), t.ket().at(p[1]),
                                         t.ket().at(p[2])};
          Tensor permuted(t.label(), bra(std::move(pbra)), ket(std::move(pket)),
                          t.aux(), Symmetry::Nonsymm, t.braket_symmetry(),
                          ColumnSymmetry::Nonsymm);
          channel->append(ex<Constant>(spin_line_sign(t.bra().at(p[0]))) *
                          ex<Tensor>(std::move(permuted)));
        }
        out->append(1, ExprPtr(channel), Product::Flatten::No);
      } else {
        throw Exception(
            "Triplet spin tracing EOM is implemented for singles, doubles and "
            "triples");
      }
    }

    ExprPtr result = out;
    expand(result);
    rapid_simplify(result);
    return result;
  };

  if (spin_labeled->is<Product>())
    return adapt_product(spin_labeled->as<Product>());
  if (spin_labeled->is<Sum>()) {
    auto out = std::make_shared<Sum>();
    for (const auto& t : *spin_labeled)
      out->append(t->is<Product>() ? adapt_product(t->as<Product>()) : t);
    return out;
  }
  return spin_labeled;
}

ExprPtr triplet_doubles_paper_combined_residual(
    const ExprPtr& TE, const container::map<Index, Index>& pair_swap,
    bool te_only = false) {
  // EFV test, te_only, ter_only
  // post-processing for te_only restores the dropped TE_swap. This is what I
  // checked, not what EFV asked me to test. Omega = (3*TE - TE_ps)/16 = TE/4 +
  // (TE_bs + TE_ks)/16 (null-space identity TE + TE_ps + TE_bs + TE_ks = 0).
  if (te_only) return ex<Constant>(ratio(1, 4)) * TE;
  ExprPtr TE_ps = transform_expr(TE, pair_swap);
  ExprPtr tauSymm = TE + TE_ps;
  ExprPtr tauAnti = TE - TE_ps;
  ExprPtr tauSymm_N = ex<Constant>(ratio(1, 8)) * tauSymm;
  ExprPtr tauAnti_N = ex<Constant>(ratio(1, 8)) * tauAnti;
  return ex<Constant>(ratio(1, 2)) * tauSymm_N +
         tauAnti_N;  // paper does not have this 1/2. why? it correct. Faber
                     // paper has the factor.
}

// The 18 TEE ops (6 pairings x 3 T
// positions) have rank 9. The minimal null-space-reduced
// dual of TEE is : (1/80)[ 6 TEE - TEE_ps01 - TEE_ps02 + 2 TEE_ks12 ].
ExprPtr triplet_triples_paper_combined_residual(
    const ExprPtr& TEE, const container::map<Index, Index>& pair_swap_01,
    const container::map<Index, Index>& pair_swap_02,
    const container::map<Index, Index>& ket_swap_12) {
  ExprPtr TEE_ps01 = transform_expr(TEE, pair_swap_01);
  ExprPtr TEE_ps02 = transform_expr(TEE, pair_swap_02);
  ExprPtr TEE_ks12 = transform_expr(TEE, ket_swap_12);
  return ex<Constant>(ratio(1, 160)) * (ex<Constant>(6) * TEE - TEE_ps01 -
                                        TEE_ps02 + ex<Constant>(2) * TEE_ks12);
}

}  // namespace

ExprPtr closed_shell_EOM_triplet_spintrace(
    ExprPtr const& expr, ClosedShellCCSpintraceOptions options) {
  container::svector<container::svector<Index>> ext_groups;
  const auto ext_idxs = external_indices(expr);
  for (const auto& g : ext_idxs) {
    container::svector<Index> grp;
    grp.reserve(g.size());
    for (const auto& s : g) grp.push_back(s.index());
    ext_groups.push_back(std::move(grp));
  }

  const auto n_ext = ext_idxs.size();
  if (n_ext == 0 || n_ext > 3)
    throw Exception(
        "closed_shell_EOM_triplet_spintrace: the explicitly spin-coupled "
        "triplet manifold is implemented for singles, doubles and triples");

  // spintrace_by_sector already called remove_spin_with_relabel on each sector.
  auto sectors = spintrace_by_sector(expr, ext_groups, /*triplet_R=*/true,
                                     options.triplet_amp_no_swap);
  SEQUANT_ASSERT(sectors.size() == (std::size_t{1} << n_ext));

  // V_mu: weight every sector by the sign of the spin of external group 0
  auto V = std::make_shared<Sum>();
  for (std::size_t s = 0; s < sectors.size(); ++s) {
    if (s & 1u)
      V->append(ex<Constant>(-1) * sectors[s].second);
    else
      V->append(sectors[s].second);
  }
  ExprPtr triplet = V;
  canonicalize(triplet);
  simplify(triplet);

  if (n_ext == 1) {
    // metric of the T_{ai} kets is 2*delta -> the rank-1 biorthogonal
    // for now I am calling v1/v2, but later I will just add a normalization
    // factor 1/2
    switch (options.method) {
      case BiorthogonalizationMethod::V1:
        triplet =
            biorthogonal_transform_pre_nnsproject(triplet, ext_idxs, false);
        break;
      case BiorthogonalizationMethod::V2:
        triplet =
            biorthogonal_transform_pre_nnsproject(triplet, ext_idxs, true);
        break;
      default:
        SEQUANT_ASSERT(false && "unreachable");
        abort();
    }
  } else if (n_ext == 2) {
    // paper orthonormal (1)/(2) channels (1/8) folded into one R2 residual
    const auto& g0 = ext_idxs.at(0);
    const auto& g1 = ext_idxs.at(1);
    SEQUANT_ASSERT(g0.size() == 2 && g1.size() == 2);
    container::map<Index, Index> pair_swap;
    const Index b0 = get_bra_idx(g0);
    const Index b1 = get_bra_idx(g1);
    const Index k0 = get_ket_idx(g0);
    const Index k1 = get_ket_idx(g1);
    pair_swap.emplace(b0, b1);
    pair_swap.emplace(b1, b0);
    pair_swap.emplace(k0, k1);
    pair_swap.emplace(k1, k0);

    triplet = triplet_doubles_paper_combined_residual(triplet, pair_swap,
                                                      options.triplet_te_only);
    // the te_only combined residual is a Constant*Sum Product;
    // triplet_maxcoeff_compact only acts on a top-level Sum.
    simplify(triplet);
    if (options.triplet_doubles_compact)
      // Keep the verified orbit representative per hash group so the dropped
      // terms can be rebuilt exactly: paper metric keeps the -3c member of
      // each {c,c,c,-3c} group (triplet_symbolic_reconstruct, or the
      // {1,-1/3,-1/3,-1/3} four-swap reconstruction on the tensor side);
      // te_only keeps the -2c member of each {c,c,-2c} group
      // (TeNnsReconstruction, {1,-1/2,-1/2,0}).
      triplet = triplet_maxcoeff_compact(
          triplet, ext_groups,
          options.triplet_te_only ? TripletOrbitWeightKind::TeNnsReconstruction
                                  : TripletOrbitWeightKind::NnsReconstruction);
  } else {  // n_ext == 3
    if (options.triplet_te_only || options.triplet_amp_no_swap)
      throw Exception(
          "closed_shell_EOM_triplet_spintrace: te_only/ter_only are "
          "doubles-only experiments, not implemented for triples");
    const auto& g0 = ext_idxs.at(0);
    const auto& g1 = ext_idxs.at(1);
    const auto& g2 = ext_idxs.at(2);
    SEQUANT_ASSERT(g0.size() == 2 && g1.size() == 2 && g2.size() == 2);
    auto whole_pair_swap = [](const auto& ga, const auto& gb) {
      container::map<Index, Index> m;
      const Index ba = get_bra_idx(ga);
      const Index bb = get_bra_idx(gb);
      const Index ka = get_ket_idx(ga);
      const Index kb = get_ket_idx(gb);
      m.emplace(ba, bb);
      m.emplace(bb, ba);
      m.emplace(ka, kb);
      m.emplace(kb, ka);
      return m;
    };
    container::map<Index, Index> ket_swap_12;
    {
      const Index k1 = get_ket_idx(g1);
      const Index k2 = get_ket_idx(g2);
      ket_swap_12.emplace(k1, k2);
      ket_swap_12.emplace(k2, k1);
    }
    triplet = triplet_triples_paper_combined_residual(
        triplet, whole_pair_swap(g0, g1), whole_pair_swap(g0, g2), ket_swap_12);
    // the combined residual is a Constant*Sum Product;
    // triplet_maxcoeff_compact only acts on a top-level Sum.
    simplify(triplet);
    if (options.triplet_doubles_compact)
      // Keep one 36-slot-perm-orbit representative per hash group, scaled by
      // the inverse stabilizer order, so the dropped terms can be rebuilt
      // exactly (triplet_symbolic_reconstruct, or the w10[m]/5 36-perm
      // reconstruction of triplet_nns_project on the tensor side).
      triplet = triplet_maxcoeff_compact(triplet, ext_groups);
  }
  simplify(triplet);
  std::wcout << "closed_shell_EOM_triplet_spintrace size: " << triplet->size()
             << " terms\n";
  return triplet;
}

// number-conserving spin trace that keeps external spin sectors separate.
// returns one spin-free expression per external spin string (αα.., αβ.., ..,
// ββ..). summing all returned sectors reproduces generic spintrace() exactly,
// so each entry is a genuine projection of the same result onto one external
// sector.
container::svector<std::pair<std::wstring, ExprPtr>> spintrace_by_sector(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool triplet_R, bool amp_no_swap) {
  ExprPtr work = expr->clone();
  if (has_tensor(work, reserved::antisymm_label())) {
    work = expand_A_op(work);
  }

  const std::size_t n_ext = ext_index_groups.size();
  SEQUANT_ASSERT(n_ext <= 31);

  // spintrace_impl places external groups at the head of its per-product
  // index_groups list, so external group g sits at bit g of spincase_bitstr
  // (independent of the product's internal-index count).
  container::svector<std::pair<std::wstring, ExprPtr>> sectors;
  for (std::uint32_t ext_val = 0; ext_val < (1u << n_ext); ++ext_val) {
    std::function<bool(uint64_t)> predicate = [ext_val, n_ext](uint64_t sc) {
      for (std::size_t g = 0; g < n_ext; ++g) {
        const uint64_t spin_bit = (sc >> g) & 1u;
        const uint64_t want = (ext_val >> g) & 1u;
        if (spin_bit != want) return false;
      }
      return true;
    };

    ExprPtr sector = spintrace_impl(work, ext_index_groups,
                                    /*spinfree_index_spaces=*/false, predicate);
    detail::reset_idx_tags(sector);
    canonicalize(sector);
    simplify(sector);
    if (triplet_R) sector = triplet_adapt_amplitudes(sector, amp_no_swap);
    sector = remove_spin_with_relabel(sector);
    detail::reset_idx_tags(sector);
    canonicalize(sector);
    simplify(sector);

    std::wstring label;
    for (std::size_t g = 0; g < n_ext; ++g)
      label += ((ext_val >> g) & 1u) ? L'β' : L'α';
    sectors.emplace_back(std::move(label), std::move(sector));
  }
  return sectors;
}

std::vector<ExprPtr> open_shell_CC_spintrace_by_sector(const ExprPtr& expr) {
  SEQUANT_ASSERT(expr->is<Sum>() || expr->is<Product>());
  Tensor A = expr.is<Sum>() ? expr->at(0)->at(0)->as<Tensor>()
                            : expr->at(0)->as<Tensor>();
  SEQUANT_ASSERT(A.label() == reserved::antisymm_label());

  container::svector<container::svector<Index>> ext_groups;
  for (const auto& grp : external_indices(expr)) {
    container::svector<Index> indices;
    indices.reserve(grp.size());
    for (const auto& slotted : grp) indices.push_back(slotted.index());
    ext_groups.push_back(std::move(indices));
  }

  auto sectors = spintrace_by_sector(expr, ext_groups);
  std::vector<ExprPtr> result;
  result.reserve(sectors.size());
  for (auto& [label, sector] : sectors) {
    (void)label;
    result.push_back(std::move(sector));
  }
  return result;
}

}  // namespace sequant::mbpt
