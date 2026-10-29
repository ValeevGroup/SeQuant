#include <SeQuant/core/utility/polynomial_division.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/conversion.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant {

namespace {

using IndexMap = container::map<Index, Index>;
using Indices = container::set<Index>;
using scalar_type = Product::scalar_type;

/// Options for canonicalizing individual terms. Named indices are treated as
/// meaningful (as is done when canonicalizing the summands of a Sum) such that
/// equal monomials end up in equal canonical forms.
CanonicalizeOptions term_canonicalization_options() {
  return CanonicalizeOptions::default_options().copy_and_set(
      CanonicalizeOptions::IgnoreNamedIndexLabel::No);
}

/// Options for canonicalizing tensor networks for matching purposes. The
/// structure of the canonical form must not depend on the labels of the named
/// indices such that networks differing only in their index labels can be
/// compared slot by slot.
CanonicalizeOptions matching_canonicalization_options() {
  return CanonicalizeOptions::default_options().copy_and_set(
      CanonicalizeOptions::IgnoreNamedIndexLabel::Yes);
}

/// @returns The number of (non-proto) occurrences of every index in @p expr
container::map<Index, std::int64_t> index_counts(const Expr &expr) {
  container::map<Index, std::int64_t> counts;
  for (const auto &[idx, counter] : get_used_indices_with_counts(expr)) {
    if (counter.nonproto() > 0) {
      counts.emplace(idx, counter.nonproto());
    }
  }
  return counts;
}

/// @returns The indices occurring exactly once in @p term
Indices open_indices(const Expr &term) {
  Indices result;
  for (const auto &[idx, count] : index_counts(term)) {
    if (count == 1) {
      result.insert(idx);
    }
  }
  return result;
}

/// @returns The indices occurring exactly once in any of the terms of @p expr
Indices external_indices(const ExprPtr &expr) {
  if (!expr.is<Sum>()) {
    return open_indices(*expr);
  }
  Indices result;
  for (const ExprPtr &summand : expr.as<Sum>().summands()) {
    for (const Index &idx : open_indices(*summand)) {
      result.insert(idx);
    }
  }
  return result;
}

/// @returns All indices (including proto indices) used in @p expr
Indices used_indices(const Expr &expr) { return get_used_indices(expr); }

/// @returns @p expr as a Product (wrapping it, if necessary)
ExprPtr to_product(const ExprPtr &expr) {
  if (expr.is<Product>()) {
    return expr.clone();
  }
  if (expr.is<Constant>()) {
    return ex<Product>(expr.as<Constant>().value(), ExprPtrList{});
  }
  return ex<Product>(ExprPtrList{expr});
}

/// Canonicalizes the given Product in-place
/// @returns @p term
ExprPtr &canonicalize_term(ExprPtr &term) {
  SEQUANT_ASSERT(term.is<Product>());
  // Merge repeated scalar factors (e.g. x x -> x^2) such that every monomial
  // has a unique representation
  rapid_simplify(term);
  term = to_product(term);
  ExprPtr byproduct = term->canonicalize(term_canonicalization_options());
  if (byproduct) {
    SEQUANT_ASSERT(byproduct.is<Constant>());
    term.as<Product>().scale(byproduct.as<Constant>().value());
  }
  return term;
}

/// @returns The terms of @p expr, each represented as a canonicalized Product
container::svector<ExprPtr> canonical_terms(const ExprPtr &expr) {
  container::svector<ExprPtr> terms;
  auto add = [&terms](const ExprPtr &term) {
    ExprPtr product = to_product(term);
    canonicalize_term(product);
    if (!product.as<Product>().is_zero()) {
      terms.push_back(std::move(product));
    }
  };

  if (expr.is<Sum>()) {
    for (const ExprPtr &summand : expr.as<Sum>().summands()) {
      add(summand);
    }
  } else {
    add(expr);
  }

  return terms;
}

/// @returns Whether @p factor is a scalar factor that participates in
/// (polynomial) division
bool is_scalar_factor(const ExprPtr &factor) {
  return factor.is<Variable>() || factor.is<Power>();
}

/// A scalar factor represented as base^exponent
struct ScalarPower {
  ExprPtr base;
  bool conjugated = false;
  rational exponent = 1;

  explicit ScalarPower(const ExprPtr &factor) {
    if (factor.is<Power>()) {
      const auto &power = factor.as<Power>();
      base = power.base();
      conjugated = power.conjugated();
      exponent = power.exponent();
    } else {
      base = factor;
    }
  }

  bool same_base(const ScalarPower &other) const {
    return conjugated == other.conjugated && *base == *other.base;
  }

  ExprPtr to_expr() const {
    if (exponent == 1 && !conjugated) {
      return base->clone();
    }
    ExprPtr result = ex<Power>(base->clone(), exponent);
    if (conjugated) {
      result.as<Power>().conjugate();
    }
    Power::flatten(result);
    return result;
  }
};

/// Divides the scalar factors @p dividend by @p divisor
/// @returns The remaining scalar factors or std::nullopt if not divisible
std::optional<container::svector<ExprPtr>> divide_scalars(
    const container::svector<ExprPtr> &dividend,
    const container::svector<ExprPtr> &divisor) {
  container::svector<std::optional<ScalarPower>> remaining;
  for (const ExprPtr &factor : dividend) {
    remaining.emplace_back(ScalarPower(factor));
  }

  for (const ExprPtr &factor : divisor) {
    const ScalarPower current(factor);
    auto it = std::ranges::find_if(
        remaining, [&](const std::optional<ScalarPower> &candidate) {
          // N.B.: base^a is considered divisible by base^b if both exponents
          // have the same sign and |a| >= |b|
          return candidate.has_value() && candidate->same_base(current) &&
                 (candidate->exponent < 0) == (current.exponent < 0) &&
                 abs(candidate->exponent) >= abs(current.exponent);
        });
    if (it == remaining.end()) {
      return std::nullopt;
    }

    (*it)->exponent -= current.exponent;
    if ((*it)->exponent == 0) {
      it->reset();
    }
  }

  container::svector<ExprPtr> result;
  for (const std::optional<ScalarPower> &power : remaining) {
    if (power.has_value()) {
      result.push_back(power->to_expr());
    }
  }
  return result;
}

/// Produces fresh indices that are not used anywhere in the expressions
/// involved in a division
class FreshIndices {
 public:
  explicit FreshIndices(Indices used)
      : used_(std::make_shared<Indices>(std::move(used))),
        factory_(std::make_unique<IndexFactory>(
            [used = used_](const Index &idx) { return !used->contains(idx); },
            1)) {}

  Index make(const Index &like) {
    Index result = factory_->make(like);
    used_->insert(result);
    return result;
  }

 private:
  std::shared_ptr<Indices> used_;
  // N.B. IndexFactory is not movable
  std::unique_ptr<IndexFactory> factory_;
};

/// Renames the dummy indices of every term of @p expr such that they don't
/// coincide with any of @p avoid while being as "simple" as possible
ExprPtr tidy_dummies(const ExprPtr &expr, const Indices &avoid) {
  auto tidy_term = [&avoid](const ExprPtr &term) -> ExprPtr {
    const auto counts = index_counts(*term);
    IndexMap map;
    Indices used = used_indices(*term);
    for (const auto &[idx, count] : counts) {
      if (count > 1) {
        // Only rename if necessary
        if (avoid.contains(idx)) {
          map.emplace(idx, Index{});
        }
      }
    }
    if (map.empty()) {
      return term;
    }

    IndexFactory factory(
        [&](const Index &idx) {
          return !avoid.contains(idx) && !used.contains(idx);
        },
        1);
    for (auto &[from, to] : map) {
      to = factory.make(from);
      used.insert(to);
    }

    return transform_expr(*term, map);
  };

  if (expr.is<Sum>()) {
    auto result = std::make_shared<Sum>(ExprPtrList{});
    for (const ExprPtr &summand : expr.as<Sum>().summands()) {
      result->append(tidy_term(summand));
    }
    return result;
  }

  return tidy_term(expr);
}

/// @returns The total degree of the given monomial, i.e. the number of its
/// non-scalar factors plus the sum of the exponents of its scalar factors
rational total_degree(const Product &term) {
  rational result = 0;
  for (const ExprPtr &factor : term.factors()) {
    result += factor.is<Power>() ? factor.as<Power>().exponent() : 1;
  }
  return result;
}

/// A single term of the divisor
struct DivisorTerm {
  /// The term as a Product
  ExprPtr expr;
  scalar_type scalar = 1;
  container::svector<ExprPtr> scalars;
  container::svector<ExprPtr> tensors;

  rational degree() const { return total_degree(expr.as<Product>()); }
};

/// The divisor, prepared for division
struct Divisor {
  container::svector<DivisorTerm> terms;
  Indices externals;
};

/// Brings the (simplified) divisor into the form needed for division.
/// Dummy indices of the divisor are replaced by fresh indices to rule out
/// any clashes with indices of the dividend.
Divisor prepare_divisor(const ExprPtr &divisor, FreshIndices &fresh) {
  Divisor result;

  bool first = true;
  for (ExprPtr &term : canonical_terms(divisor)) {
    const auto &product = term.as<Product>();

    Indices externals = open_indices(product);
    if (first) {
      result.externals = externals;
      first = false;
    } else {
      SEQUANT_ENFORCE(externals == result.externals,
                      "All terms of the divisor must have the same external "
                      "indices");
    }

    IndexMap dummy_map;
    for (const auto &[idx, count] : index_counts(product)) {
      if (count > 1) {
        dummy_map.emplace(idx, fresh.make(idx));
      }
    }

    DivisorTerm div_term;
    div_term.expr = transform_expr(product, dummy_map);
    div_term.scalar = product.scalar();
    for (const ExprPtr &factor : div_term.expr.as<Product>().factors()) {
      if (is_scalar_factor(factor)) {
        div_term.scalars.push_back(factor);
      } else {
        SEQUANT_ENFORCE(factor.is<Tensor>(),
                        "Divisor terms may only consist of tensors, "
                        "variables, constants and powers thereof");
        div_term.tensors.push_back(factor);
      }
    }

    result.terms.push_back(std::move(div_term));
  }

  SEQUANT_ENFORCE(!result.terms.empty(), "Can't divide by zero");

  return result;
}

/// The result of matching the tensors of a divisor term to a subset of the
/// tensors of a dividend term
struct NetworkMatch {
  /// maps external indices of the divisor term to indices of the dividend
  IndexMap sigma;
  /// subset == ratio * transform_expr(divisor_tensors, sigma)
  scalar_type ratio = 1;
};

/// Checks whether @p subset equals (up to a scalar factor and renaming of
/// dummy indices) the network formed by @p divisor_tensors with its external
/// indices relabeled according to @p sigma
/// @returns The ratio subset / relabeled divisor tensors or std::nullopt if
/// they don't match
std::optional<scalar_type> verify_binding(
    const container::svector<ExprPtr> &subset,
    const container::svector<ExprPtr> &divisor_tensors, const IndexMap &sigma) {
  container::set<Index> named;
  for (const auto &[external, image] : sigma) {
    named.insert(image);
  }
  const auto opts = term_canonicalization_options().copy_and_set(
      std::optional<container::set<Index>>(std::move(named)));

  auto canonicalized = [&opts](auto &&tensors) {
    auto product = std::make_shared<Product>(tensors, Product::Flatten::No);
    ExprPtr byproduct = product->canonicalize(opts);
    if (byproduct) {
      product->scale(byproduct.as<Constant>().value());
    }
    return product;
  };

  container::svector<ExprPtr> mapped;
  for (const ExprPtr &tensor : divisor_tensors) {
    mapped.push_back(transform_expr(*tensor, sigma));
  }

  auto lhs = canonicalized(subset);
  auto rhs = canonicalized(mapped);
  if (lhs->is_zero() || rhs->is_zero() ||
      lhs->factors().size() != rhs->factors().size()) {
    return std::nullopt;
  }
  for (std::size_t i = 0; i < lhs->factors().size(); ++i) {
    if (*lhs->factor(i) != *rhs->factor(i)) {
      return std::nullopt;
    }
  }

  return lhs->scalar() / rhs->scalar();
}

/// Attempts to match the tensor network formed by @p subset to the one formed
/// by @p divisor_tensors (up to relabeling of indices)
std::optional<NetworkMatch> match_networks(
    const container::svector<ExprPtr> &subset,
    const container::svector<ExprPtr> &divisor_tensors,
    const Indices &divisor_externals) {
  SEQUANT_ASSERT(subset.size() == divisor_tensors.size());

  if (subset.empty()) {
    return NetworkMatch{.sigma = {}, .ratio = 1};
  }

  const auto opts = matching_canonicalization_options();

  auto canonicalized = [&opts](const container::svector<ExprPtr> &tensors) {
    auto product = std::make_shared<Product>(tensors, Product::Flatten::No);
    ExprPtr byproduct = product->canonicalize(opts);
    if (byproduct) {
      product->scale(byproduct.as<Constant>().value());
    }
    return product;
  };

  auto lhs = canonicalized(subset);
  auto rhs = canonicalized(divisor_tensors);

  if (lhs->is_zero() || rhs->is_zero() ||
      lhs->factors().size() != rhs->factors().size()) {
    return std::nullopt;
  }

  // Derive the index mapping by comparing slot by slot
  IndexMap sigma;
  IndexMap dummies;
  Indices images;
  for (std::size_t i = 0; i < lhs->factors().size(); ++i) {
    const ExprPtr &lhs_factor = lhs->factor(i);
    const ExprPtr &rhs_factor = rhs->factor(i);
    if (!lhs_factor.is<AbstractTensor>() || !rhs_factor.is<AbstractTensor>()) {
      return std::nullopt;
    }

    const auto lhs_slots = slots(lhs_factor.as<AbstractTensor>()) |
                           ranges::to<container::svector<Index>>;
    const auto rhs_slots = slots(rhs_factor.as<AbstractTensor>()) |
                           ranges::to<container::svector<Index>>;
    if (lhs_slots.size() != rhs_slots.size()) {
      return std::nullopt;
    }

    for (std::size_t slot = 0; slot < rhs_slots.size(); ++slot) {
      const Index &lhs_idx = lhs_slots[slot];
      const Index &rhs_idx = rhs_slots[slot];

      if (lhs_idx.nonnull() != rhs_idx.nonnull()) {
        return std::nullopt;
      }
      if (!rhs_idx.nonnull()) {
        continue;
      }
      if (lhs_idx.space() != rhs_idx.space()) {
        return std::nullopt;
      }

      IndexMap &map = divisor_externals.contains(rhs_idx) ? sigma : dummies;
      auto [it, inserted] = map.emplace(rhs_idx, lhs_idx);
      if (inserted) {
        if (!images.insert(lhs_idx).second) {
          // mapping is not injective
          return std::nullopt;
        }
      } else if (it->second != lhs_idx) {
        return std::nullopt;
      }
    }
  }

  // Verify that the mapping indeed maps the divisor onto the subset
  IndexMap full_map = sigma;
  full_map.insert(dummies.begin(), dummies.end());
  ExprPtr mapped = transform_expr(*rhs, full_map);
  const auto &mapped_product = mapped.as<Product>();
  for (std::size_t i = 0; i < lhs->factors().size(); ++i) {
    if (*lhs->factor(i) != *mapped_product.factor(i)) {
      return std::nullopt;
    }
  }

  if (sigma.size() != divisor_externals.size()) {
    return std::nullopt;
  }

  // In case of (permutational) symmetries, there might be multiple equivalent
  // bindings of the divisor's externals of which the canonicalization picks
  // an arbitrary one. Prefer the binding that maps as many of the externals
  // onto themselves as possible (required for IndexMatching::Exact and
  // helpful for grouping in IndexMatching::Relabel).
  Indices images_of_externals;
  for (const auto &[external, image] : sigma) {
    images_of_externals.insert(image);
  }
  IndexMap preferred;
  for (const Index &external : divisor_externals) {
    if (images_of_externals.contains(external)) {
      preferred.emplace(external, external);
    }
  }
  if (preferred.size() != sigma.size()) {
    Indices leftover_images = images_of_externals;
    for (const auto &[external, image] : preferred) {
      leftover_images.erase(image);
    }
    for (const auto &[external, image] : sigma) {
      if (preferred.contains(external)) {
        continue;
      }
      auto it = std::ranges::find_if(leftover_images, [&](const Index &idx) {
        return idx.space() == external.space();
      });
      if (it == leftover_images.end()) {
        break;
      }
      preferred.emplace(external, *it);
      leftover_images.erase(it);
    }
  }
  if (preferred.size() == sigma.size() && preferred != sigma) {
    if (auto ratio = verify_binding(subset, divisor_tensors, preferred)) {
      return NetworkMatch{.sigma = std::move(preferred), .ratio = *ratio};
    }
  }

  return NetworkMatch{.sigma = std::move(sigma),
                      .ratio = lhs->scalar() / rhs->scalar()};
}

/// The result of dividing a single term of the dividend by a single term of
/// the divisor
struct TermMatch {
  /// The quotient term (as a Product)
  ExprPtr quotient;
  /// The binding of the divisor's external indices (only non-identity
  /// entries)
  IndexMap binding;
};

/// Holds the state shared between all term divisions within a division
struct DivisionContext {
  IndexMatching matching;
  /// The external indices of the dividend
  Indices dividend_externals;
  /// Labels to use for those divisor externals that are bound to dummies of a
  /// dividend term but can't use their own label as that one is an external
  /// index of the dividend
  IndexMap alternative_labels;
  FreshIndices fresh;
};

/// Calls @p callback for every subset of positions in @p candidates that
/// consists of tensors with the same labels as @p tensors
void for_each_tensor_subset(
    const container::svector<std::size_t> &candidates,
    const Product::factors_type &factors,
    const container::svector<ExprPtr> &tensors,
    const std::function<void(const container::svector<std::size_t> &)>
        &callback) {
  // Group the divisor's tensors by their "kind" (label and number of slots).
  // For every kind that occurs k times, we have to choose k (distinct)
  // candidates of the same kind. Since the order within such a choice is
  // irrelevant, we directly generate combinations (rather than permutations).
  auto same_kind = [](const Tensor &lhs, const Tensor &rhs) {
    return lhs.label() == rhs.label() && lhs.num_slots() == rhs.num_slots();
  };

  struct Kind {
    /// how many tensors of this kind the divisor contains
    std::size_t count;
    /// positions (in factors) of the candidates of this kind
    container::svector<std::size_t> candidates;
  };
  container::svector<Kind> kinds;
  container::svector<const Tensor *> representatives;
  for (const ExprPtr &tensor : tensors) {
    const auto &current = tensor.as<Tensor>();
    auto it = std::ranges::find_if(representatives, [&](const Tensor *rep) {
      return same_kind(*rep, current);
    });
    if (it != representatives.end()) {
      ++kinds[std::distance(representatives.begin(), it)].count;
      continue;
    }

    Kind kind{.count = 1, .candidates = {}};
    for (std::size_t pos : candidates) {
      if (same_kind(factors[pos].as<Tensor>(), current)) {
        kind.candidates.push_back(pos);
      }
    }
    kinds.push_back(std::move(kind));
    representatives.push_back(&current);
  }

  for (const Kind &kind : kinds) {
    if (kind.candidates.size() < kind.count) {
      return;
    }
  }

  container::svector<std::size_t> current;
  // Chooses the next element for the current kind, starting at candidate
  // position start; chosen is the number of already chosen elements of the
  // current kind
  std::function<void(std::size_t, std::size_t, std::size_t)> recurse =
      [&](std::size_t kind_idx, std::size_t start, std::size_t chosen) {
        if (kind_idx == kinds.size()) {
          auto sorted = current;
          std::ranges::sort(sorted);
          callback(sorted);
          return;
        }

        const Kind &kind = kinds[kind_idx];
        if (chosen == kind.count) {
          recurse(kind_idx + 1, 0, 0);
          return;
        }

        for (std::size_t i = start;
             i + (kind.count - chosen) <= kind.candidates.size(); ++i) {
          current.push_back(kind.candidates[i]);
          recurse(kind_idx, i + 1, chosen + 1);
          current.pop_back();
        }
      };

  recurse(0, 0, 0);
}

/// Divides the dividend term @p term by the divisor term @p divisor_term
/// @returns All possible ways in which the division can be performed (may be
/// empty if @p term is not divisible)
container::svector<TermMatch> divide_term(const ExprPtr &term,
                                          const DivisorTerm &divisor_term,
                                          const Divisor &divisor,
                                          DivisionContext &ctx) {
  const auto &product = term.as<Product>();

  container::svector<std::size_t> tensor_positions;
  container::svector<ExprPtr> scalars;
  for (std::size_t i = 0; i < product.factors().size(); ++i) {
    const ExprPtr &factor = product.factor(i);
    if (is_scalar_factor(factor)) {
      scalars.push_back(factor);
    } else if (factor.is<Tensor>()) {
      tensor_positions.push_back(i);
    }
  }

  auto remaining_scalars = divide_scalars(scalars, divisor_term.scalars);
  if (!remaining_scalars.has_value()) {
    return {};
  }
  if (tensor_positions.size() < divisor_term.tensors.size()) {
    return {};
  }

  const auto term_counts = index_counts(product);

  container::svector<TermMatch> matches;

  for_each_tensor_subset(
      tensor_positions, product.factors(), divisor_term.tensors,
      [&](const container::svector<std::size_t> &positions) {
        container::svector<ExprPtr> subset;
        for (std::size_t pos : positions) {
          subset.push_back(product.factor(pos));
        }

        std::optional<NetworkMatch> network_match =
            match_networks(subset, divisor_term.tensors, divisor.externals);
        if (!network_match.has_value()) {
          return;
        }

        // Determine how the divisor's externals will be bound and which
        // dummies of the term need renaming as a consequence
        IndexMap binding;
        IndexMap renaming;
        Indices homes;
        for (const auto &[external, image] : network_match->sigma) {
          auto count_it = term_counts.find(image);
          SEQUANT_ASSERT(count_it != term_counts.end());
          if (count_it->second == 1) {
            // Bound to an external index of the term
            if (ctx.matching == IndexMatching::Exact && image != external) {
              return;
            }
            binding.emplace(external, image);
            continue;
          }

          // Bound to a dummy index of the term
          Index home = external;
          if (ctx.dividend_externals.contains(external)) {
            if (ctx.matching == IndexMatching::Exact) {
              return;
            }
            auto it = ctx.alternative_labels.find(external);
            if (it == ctx.alternative_labels.end()) {
              it = ctx.alternative_labels
                       .emplace(external, ctx.fresh.make(external))
                       .first;
            }
            home = it->second;
          }

          binding.emplace(external, home);
          homes.insert(home);
          if (image != home) {
            renaming.emplace(image, home);
          }
        }

        // Move dummies of the term that would clash with the home labels out
        // of the way
        for (const Index &home : homes) {
          if (term_counts.contains(home) && !renaming.contains(home)) {
            bool is_own_image = false;
            for (const auto &[external, image] : network_match->sigma) {
              if (image == home && binding.at(external) == home) {
                is_own_image = true;
                break;
              }
            }
            if (!is_own_image) {
              renaming.emplace(home, ctx.fresh.make(home));
            }
          }
        }

        // Assemble quotient
        auto quotient = std::make_shared<Product>();
        quotient->scale(product.scalar() * network_match->ratio /
                        divisor_term.scalar);
        for (const ExprPtr &factor : *remaining_scalars) {
          quotient->append(1, factor, Product::Flatten::No);
        }
        for (std::size_t i = 0; i < product.factors().size(); ++i) {
          const ExprPtr &factor = product.factor(i);
          if (is_scalar_factor(factor) ||
              std::ranges::find(positions, i) != positions.end()) {
            continue;
          }
          quotient->append(1, factor, Product::Flatten::No);
        }

        TermMatch match;
        match.quotient = transform_expr(*quotient, renaming);
        for (const auto &[external, bound] : binding) {
          if (external != bound) {
            match.binding.emplace(external, bound);
          }
        }

        matches.push_back(std::move(match));
      });

  return matches;
}

/// The terms of an expression undergoing division, indexed for quick lookup
/// of proportional terms
class TermTable {
 public:
  enum class State { Pending, Remainder, Dead };

  struct Entry {
    ExprPtr term;
    State state;
    rational degree;
  };

  /// Adds @p term (which must be canonical) to the table, combining it with
  /// an existing proportional term, if any
  void add(ExprPtr term) {
    SEQUANT_ASSERT(term.is<Product>());
    if (term->is_zero()) {
      return;
    }

    if (auto idx = find(term); idx.has_value()) {
      Entry &entry = entries_[*idx];
      // N.B.: must not modify the entry in-place as the term might be shared
      auto combined = std::make_shared<Product>(entry.term.as<Product>());
      combined->add_identical(term);
      entry.term = combined;
      if (combined->is_zero()) {
        entry.state = State::Dead;
      }
      return;
    }

    by_hash_.emplace(term->hash_value(), entries_.size());
    rational term_degree = total_degree(term.as<Product>());
    pending_.emplace(term_degree, entries_.size());
    entries_.push_back(Entry{.term = std::move(term),
                             .state = State::Pending,
                             .degree = std::move(term_degree)});
  }

  /// @returns The position of a live entry proportional to @p term
  std::optional<std::size_t> find(const ExprPtr &term) const {
    auto [begin, end] = by_hash_.equal_range(term->hash_value());
    for (auto it = begin; it != end; ++it) {
      const Entry &entry = entries_[it->second];
      if (entry.state != State::Dead && proportional_to{}(entry.term, term)) {
        return it->second;
      }
    }
    return std::nullopt;
  }

  /// @returns The pending entry of highest degree (if any). As in standard
  /// polynomial division, processing terms in order of decreasing degree
  /// ensures that cancellation terms (which tend to be of lower degree) are
  /// produced before they are processed themselves.
  /// Among entries of the same degree, the one added first is returned.
  std::optional<std::size_t> next() {
    while (!pending_.empty()) {
      const std::size_t idx = pending_.top().second;
      pending_.pop();
      // N.B.: entries are only ever pushed once but they might have changed
      // their state since
      if (entries_[idx].state == State::Pending) {
        return idx;
      }
    }
    return std::nullopt;
  }

  Entry &operator[](std::size_t idx) { return entries_[idx]; }

  /// @returns The sum of all live terms
  ExprPtr sum() const {
    auto result = std::make_shared<Sum>(ExprPtrList{});
    for (const Entry &entry : entries_) {
      if (entry.state != State::Dead) {
        result->append(entry.term);
      }
    }
    return result;
  }

 private:
  std::vector<Entry> entries_;
  std::unordered_multimap<std::size_t, std::size_t> by_hash_;

  /// (degree, position) of entries that have been pending when added
  using QueueItem = std::pair<rational, std::size_t>;
  struct QueueCompare {
    bool operator()(const QueueItem &lhs, const QueueItem &rhs) const {
      // max-heap w.r.t. degree, min-heap w.r.t. position
      if (lhs.first != rhs.first) {
        return lhs.first < rhs.first;
      }
      return lhs.second > rhs.second;
    }
  };
  std::priority_queue<QueueItem, std::vector<QueueItem>, QueueCompare> pending_;
};

/// @returns The number of terms in @p expr
std::size_t num_terms(const ExprPtr &expr) {
  if (expr.is<Sum>()) {
    return expr.as<Sum>().summands().size();
  }
  if (expr.is<Constant>() && expr.as<Constant>().is_zero()) {
    return 0;
  }
  return 1;
}

/// @returns The divisor term relabeled according to @p binding, with its
/// dummy indices renamed such that they don't coincide with any of @p avoid
ExprPtr bound_divisor_term(const DivisorTerm &term, const IndexMap &binding,
                           const Indices &avoid = {}) {
  return tidy_dummies(transform_expr(*term.expr, binding), avoid);
}

/// Thrown if the division doesn't terminate for a given pivot
class NonTermination : public Exception {
 public:
  using Exception::Exception;
};

/// The expansion of divisor * quotient for a candidate quotient
struct QuotientExpansion {
  /// The quotient
  ExprPtr quotient;
  /// All terms of divisor * quotient except the one corresponding to the
  /// divided term
  container::svector<ExprPtr> companions;
};

/// Determines the coefficient of the candidate quotient of @p match such that
/// divisor * quotient reproduces @p term exactly (taking into account that
/// several divisor terms may map onto the same monomial). This is done by
/// expanding divisor * quotient explicitly.
/// @returns The quotient and the remaining terms of divisor * quotient or
/// std::nullopt if the candidate quotient doesn't reproduce @p term
std::optional<QuotientExpansion> expand_quotient(const ExprPtr &term,
                                                 const TermMatch &match,
                                                 const Divisor &divisor) {
  auto unit_quotient = std::make_shared<Product>(match.quotient.as<Product>());
  SEQUANT_ASSERT(!unit_quotient->is_zero());
  unit_quotient->scale(1 / unit_quotient->scalar());
  const ExprPtr unit_quotient_ptr = unit_quotient;
  const Indices quotient_indices = used_indices(*unit_quotient);

  container::svector<ExprPtr> expansion;
  for (const DivisorTerm &divisor_term : divisor.terms) {
    ExprPtr product = ex<Product>(ExprPtrList{
        bound_divisor_term(divisor_term, match.binding, quotient_indices),
        unit_quotient_ptr});
    canonicalize_term(product);
    if (product->is_zero()) {
      continue;
    }

    auto it = std::ranges::find_if(expansion, [&](const ExprPtr &existing) {
      return proportional_to{}(existing, product);
    });
    if (it == expansion.end()) {
      expansion.push_back(std::move(product));
    } else {
      auto combined = std::make_shared<Product>(it->as<Product>());
      combined->add_identical(product);
      *it = combined;
    }
  }

  auto self = std::ranges::find_if(expansion, [&](const ExprPtr &existing) {
    return proportional_to{}(existing, term);
  });
  if (self == expansion.end() || (*self)->is_zero()) {
    return std::nullopt;
  }

  const scalar_type factor =
      term.as<Product>().scalar() / self->as<Product>().scalar();

  QuotientExpansion result;
  unit_quotient->scale(factor);
  result.quotient = unit_quotient;
  for (auto it = expansion.begin(); it != expansion.end(); ++it) {
    if (it == self || (*it)->is_zero()) {
      continue;
    }
    auto companion = std::make_shared<Product>(it->as<Product>());
    companion->scale(factor);
    result.companions.push_back(companion);
  }

  return result;
}

/// Performs the division using the divisor term at position @p pivot as the
/// pivot ("leading") term
DivisionResult divide_with_pivot(const container::svector<ExprPtr> &terms,
                                 const Divisor &divisor, std::size_t pivot,
                                 DivisionMode mode, DivisionContext &ctx) {
  TermTable table;
  for (const ExprPtr &term : terms) {
    table.add(term);
  }

  container::svector<std::pair<IndexMap, container::svector<ExprPtr>>> groups;

  const DivisorTerm &pivot_term = divisor.terms.at(pivot);

  // The number of required iterations grows with the number of terms and
  // their degree (e.g. dividing x^n by x + 1 requires n steps)
  std::size_t degree_bound = 1;
  for (const ExprPtr &term : terms) {
    const double degree = total_degree(term.as<Product>()).convert_to<double>();
    degree_bound =
        std::max(degree_bound, static_cast<std::size_t>(std::ceil(degree)));
  }
  const std::size_t max_iterations =
      64 * (terms.size() + 1) * (divisor.terms.size() + 1) * degree_bound;
  std::size_t iterations = 0;

  while (auto idx = table.next()) {
    if (++iterations > max_iterations) {
      throw NonTermination(
          "Polynomial division did not terminate within the expected number "
          "of iterations");
    }

    const ExprPtr term = table[*idx].term;

    container::svector<TermMatch> matches =
        divide_term(term, pivot_term, divisor, ctx);

    std::optional<IndexMap> binding;
    std::optional<QuotientExpansion> chosen;
    for (const TermMatch &match : matches) {
      std::optional<QuotientExpansion> expansion =
          expand_quotient(term, match, divisor);
      if (!expansion.has_value()) {
        continue;
      }

      if (mode == DivisionMode::Split &&
          !std::ranges::all_of(expansion->companions,
                               [&](const ExprPtr &companion) {
                                 return table.find(companion).has_value();
                               })) {
        continue;
      }

      binding = match.binding;
      chosen = std::move(expansion);
      break;
    }

    if (!chosen.has_value()) {
      table[*idx].state = TermTable::State::Remainder;
      continue;
    }

    table[*idx].state = TermTable::State::Dead;

    auto group_it = std::ranges::find_if(
        groups, [&](const auto &group) { return group.first == *binding; });
    if (group_it == groups.end()) {
      groups.emplace_back(*binding, container::svector<ExprPtr>{});
      group_it = std::prev(groups.end());
    }
    group_it->second.push_back(chosen->quotient);

    for (ExprPtr &companion : chosen->companions) {
      companion.as<Product>().scale(-1);
      table.add(std::move(companion));
    }
  }

  DivisionResult result;
  result.remainder = simplify(ExprPtr(table.sum()));

  for (auto &[binding, quotient_terms] : groups) {
    ExprPtr quotient =
        simplify(ex<Sum>(quotient_terms.begin(), quotient_terms.end()));
    if (quotient.is<Constant>() && quotient.as<Constant>().is_zero()) {
      continue;
    }

    Indices bound_externals;
    for (const Index &external : divisor.externals) {
      auto it = binding.find(external);
      bound_externals.insert(it == binding.end() ? external : it->second);
    }

    // Simplification may have renamed the quotient's dummies such that they
    // coincide with the externals of the (bound) divisor
    quotient = tidy_dummies(quotient, bound_externals);

    auto divisor_expr = std::make_shared<Sum>(ExprPtrList{});
    for (const DivisorTerm &term : divisor.terms) {
      divisor_expr->append(bound_divisor_term(term, binding));
    }

    ExprPtr bound_divisor = divisor_expr;
    if (divisor_expr->summands().size() == 1) {
      bound_divisor = divisor_expr->summand(0);
    }
    bound_divisor = rapid_simplify(bound_divisor);

    result.groups.push_back(
        DivisionGroup{.divisor = tidy_dummies(bound_divisor, bound_externals),
                      .quotient = std::move(quotient),
                      .relabeling = binding});
  }

  return result;
}

/// Assembles Σ_σ factor_σ * Q_σ + R where factor_σ is produced by
/// @p make_factor
ExprPtr assemble(
    const DivisionResult &result,
    const std::function<ExprPtr(const DivisionGroup &)> &make_factor) {
  auto sum = std::make_shared<Sum>(ExprPtrList{});

  for (const DivisionGroup &group : result.groups) {
    ExprPtr factor = make_factor(group);
    ExprPtr quotient = tidy_dummies(group.quotient, external_indices(factor));
    factor = tidy_dummies(factor, used_indices(*quotient));

    sum->append(ex<Product>(ExprPtrList{std::move(factor), std::move(quotient)},
                            Product::Flatten::Once));
  }

  if (!result.exact()) {
    sum->append(result.remainder->clone());
  }

  if (sum->summands().empty()) {
    return ex<Constant>(0);
  }
  if (sum->summands().size() == 1) {
    return sum->summand(0);
  }
  return sum;
}

}  // namespace

bool DivisionResult::exact() const {
  return !remainder ||
         (remainder.is<Constant>() && remainder.as<Constant>().is_zero());
}

ExprPtr DivisionResult::to_expr() const {
  return assemble(
      *this, [](const DivisionGroup &group) { return group.divisor->clone(); });
}

namespace {

/// Implementation of divide()
/// @param reserved Indices that must not be used for any newly introduced index
DivisionResult divide_impl(const ExprPtr &dividend, const ExprPtr &divisor,
                           DivisionMode mode, IndexMatching matching,
                           const Indices &reserved) {
  SEQUANT_ENFORCE(dividend && divisor);

  ExprPtr simplified_dividend = simplify(dividend->clone());
  ExprPtr simplified_divisor = simplify(divisor->clone());

  SEQUANT_ENFORCE(!(simplified_divisor.is<Constant>() &&
                    simplified_divisor.as<Constant>().is_zero()),
                  "Can't divide by zero");

  if (simplified_dividend.is<Constant>() &&
      simplified_dividend.as<Constant>().is_zero()) {
    return DivisionResult{.groups = {}, .remainder = ex<Constant>(0)};
  }

  Indices used = used_indices(*simplified_dividend);
  for (const Index &idx : used_indices(*simplified_divisor)) {
    used.insert(idx);
  }
  used.insert(reserved.begin(), reserved.end());

  DivisionContext ctx{.matching = matching,
                      .dividend_externals = {},
                      .alternative_labels = {},
                      .fresh = FreshIndices(std::move(used))};

  const Divisor prepared_divisor =
      prepare_divisor(simplified_divisor, ctx.fresh);

  const container::svector<ExprPtr> terms =
      canonical_terms(simplified_dividend);
  for (const ExprPtr &term : terms) {
    for (const Index &idx : open_indices(*term)) {
      ctx.dividend_externals.insert(idx);
    }
  }

  // Try pivots in order of decreasing degree
  container::svector<std::size_t> pivots;
  for (std::size_t i = 0; i < prepared_divisor.terms.size(); ++i) {
    pivots.push_back(i);
  }
  std::ranges::stable_sort(pivots, [&](std::size_t lhs, std::size_t rhs) {
    return prepared_divisor.terms[lhs].degree() >
           prepared_divisor.terms[rhs].degree();
  });

  // A pivot that divides another term of the divisor can't be the leading term
  // w.r.t. any monomial order and generally leads to a non-terminating
  // division. Such pivots are only tried as a last resort.
  container::svector<std::size_t> proper_pivots;
  container::svector<std::size_t> improper_pivots;
  {
    const IndexMatching original_matching = ctx.matching;
    ctx.matching = IndexMatching::Relabel;
    for (std::size_t pivot : pivots) {
      bool proper = true;
      for (std::size_t j = 0; j < prepared_divisor.terms.size() && proper;
           ++j) {
        if (j != pivot &&
            !divide_term(prepared_divisor.terms[j].expr,
                         prepared_divisor.terms[pivot], prepared_divisor, ctx)
                 .empty()) {
          proper = false;
        }
      }
      (proper ? proper_pivots : improper_pivots).push_back(pivot);
    }
    ctx.matching = original_matching;
  }

  // The best result is the one with the fewest remainder terms. Among those,
  // prefer the one that has factored the most terms.
  auto quality = [](const DivisionResult &result) {
    std::size_t quotient_terms = 0;
    for (const DivisionGroup &group : result.groups) {
      quotient_terms += num_terms(group.quotient);
    }
    return std::make_pair(
        num_terms(result.remainder),
        std::numeric_limits<std::size_t>::max() - quotient_terms);
  };

  std::optional<DivisionResult> best;
  std::pair<std::size_t, std::size_t> best_quality = {
      std::numeric_limits<std::size_t>::max(),
      std::numeric_limits<std::size_t>::max()};
  std::optional<NonTermination> failure;
  bool any_terminated = false;
  for (const auto *candidates : {&proper_pivots, &improper_pivots}) {
    if (best.has_value()) {
      break;
    }

    for (std::size_t pivot : *candidates) {
      std::optional<DivisionResult> attempt;
      try {
        attempt = divide_with_pivot(terms, prepared_divisor, pivot, mode, ctx);
      } catch (const NonTermination &e) {
        // Not every choice of pivot leads to a terminating division. Other
        // pivots may work.
        failure = e;
        continue;
      }
      any_terminated = true;

      if (mode == DivisionMode::Strict) {
        if (attempt->exact()) {
          return std::move(*attempt);
        }
        continue;
      }

      const auto current_quality = quality(*attempt);
      if (current_quality < best_quality) {
        best_quality = current_quality;
        best = std::move(*attempt);
      }
    }
  }

  if (!any_terminated) {
    SEQUANT_ASSERT(failure.has_value());
    throw *failure;
  }

  if (mode == DivisionMode::Strict) {
    throw Exception("Dividend is not exactly divisible by the divisor");
  }

  SEQUANT_ASSERT(best.has_value());
  return std::move(best.value());
}

}  // namespace

DivisionResult divide(const ExprPtr &dividend, const ExprPtr &divisor,
                      DivisionMode mode, IndexMatching matching) {
  return divide_impl(dividend, divisor, mode, matching, {});
}

ExprPtr factor_out(const ExprPtr &expr, const ExprPtr &factor,
                   DivisionMode mode, IndexMatching matching,
                   FactorStyle style) {
  DivisionResult result = divide(expr, factor, mode, matching);

  switch (style) {
    case FactorStyle::Keep:
      return result.to_expr();
    case FactorStyle::Drop:
      if (result.groups.size() > 1) {
        throw Exception(
            "Can't drop factor as the division produced multiple, differently "
            "indexed, groups");
      }
      if (result.groups.empty()) {
        return ex<Constant>(0);
      }
      return result.groups.front().quotient->clone();
  }

  SEQUANT_UNREACHABLE;
}

ExprPtr factor_out(const ExprPtr &expr, const ExprPtr &factor,
                   const Expr &substitute, DivisionMode mode,
                   IndexMatching matching) {
  // Indices introduced by the division must not clash with indices of the
  // substitute
  DivisionResult result =
      divide_impl(expr, factor, mode, matching, used_indices(substitute));

  return assemble(result, [&substitute](const DivisionGroup &group) {
    // Make sure the substitute's dummies are neither affected by the
    // relabeling nor clash with the relabeled indices
    Indices affected;
    for (const auto &[from, to] : group.relabeling) {
      affected.insert(from);
      affected.insert(to);
    }
    ExprPtr tidied = tidy_dummies(substitute.clone(), affected);
    return transform_expr(*tidied, group.relabeling);
  });
}

ExprPtr factor_out(const ExprPtr &expr, const ResultExpr &factor,
                   DivisionMode mode, IndexMatching matching) {
  ExprPtr substitute;
  if (factor.produces_tensor()) {
    substitute = ex<Tensor>(factor.result_as_tensor());
  } else {
    substitute = ex<Variable>(factor.result_as_variable());
  }

  return factor_out(expr, factor.expression(), *substitute, mode, matching);
}

ResultExpr factor_out(const ResultExpr &expr, const ExprPtr &factor,
                      DivisionMode mode, IndexMatching matching) {
  ResultExpr result = expr.clone();
  result =
      factor_out(expr.expression(), factor, mode, matching, FactorStyle::Keep);
  return result;
}

ResultExpr factor_out(const ResultExpr &expr, const ResultExpr &factor,
                      DivisionMode mode, IndexMatching matching) {
  ResultExpr result = expr.clone();
  result = factor_out(expr.expression(), factor, mode, matching);
  return result;
}

}  // namespace sequant
