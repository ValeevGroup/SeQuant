#include "judge.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <algorithm>
#include <exception>

namespace sequant::deriv::test {

namespace {

/// whether @p x is real by @p data
bool is_real(const AbstractTensor& x, const Data& data) {
  return data.real.contains(std::wstring(strip_adjoint_label(x._label())));
}

/// the first term of @p r that breaks validity, and how; empty if none does.
/// Each term of @p r , read in this order (simplify throws on what (b) and (c)
/// find): (b) no index carries two equal protoindices; (c) an index in more
/// than two non-proto slots is kept, carried by a projector or a protoindex of
/// an index in one slot; (d) its free indices, those a projector carries
/// aside, are @p keep , a protoindex of an index a projector carries being
/// free like one of any other free index; (a) simplify accepts it. Slots are
/// counted by core's get_used_indices_with_counts, the alternatives of a
/// factor that is a Sum apart
std::string validity_violation(const ExprPtr& r,
                               const container::set<Index>& keep) {
  for (const auto& term : summands_of(r)) {
    const auto where = " in " + text_of(term);
    for (const auto& i : get_used_indices(*term)) {
      const auto& p = i.proto_indices();
      if (container::set<Index>(p.begin(), p.end()).size() != p.size())
        return "(b) " + toUtf8(i.full_label()) +
               " carries two equal protoindices" + where;
    }
    // a factor that is a Sum holds alternatives: count each product apart
    for (const auto& product : products_of(term)) {
      const auto counts = get_used_indices_with_counts(*product);
      auto [free, projected] = indices_of_product(product);
      for (const auto& i : projected)
        free.insert(i.proto_indices().begin(), i.proto_indices().end());
      for (const auto& [i, count] : counts) {
        if (count.nonproto() <= 2 || keep.contains(i) || projected.contains(i))
          continue;
        const bool of_one_slot =
            std::ranges::any_of(counts, [&i](const auto& e) {
              return e.second.nonproto() == 1 &&
                     std::ranges::count(e.first.proto_indices(), i) != 0;
            });
        if (!of_one_slot)
          return "(c) " + toUtf8(i.full_label()) +
                 " occupies more than two slots" + where;
      }
      for (const auto& i : keep)
        if (!projected.contains(i) && !free.contains(i))
          return "(d) the kept " + toUtf8(i.full_label()) + " is not free" +
                 where;
      for (const auto& i : free)
        if (!projected.contains(i) && !keep.contains(i))
          return "(d) " + toUtf8(i.full_label()) + " is free and not kept" +
                 where;
    }
    try {
      simplify(term->clone());
    } catch (const std::exception& e) {
      return "(a) simplify throws" + where + ": " + e.what();
    }
  }
  return {};
}

}  // namespace

std::size_t group_order(const ExprPtr& x, const Data& data) {
  if (!x->is<AbstractTensor>()) return 1;
  const auto& t = x->as<AbstractTensor>();
  return rearrangements_of(t, is_real(t, data)).size();
}

std::size_t block_order(const ExprPtr& x, const Data& data) {
  if (!x->is<AbstractTensor>()) return 1;
  const auto& t = x->as<AbstractTensor>();
  const auto arrangements = rearrangements_of(t, is_real(t, data));
  const auto& as_is = arrangements.front().slots;
  return static_cast<std::size_t>(
      std::ranges::count_if(arrangements, [&as_is](const auto& r) {
        for (std::size_t k = 0; k != as_is.size(); ++k)
          if (r.slots[k].basis() != as_is[k].basis()) return false;
        return true;
      }));
}

std::string validity(const ExprPtr& f, const ExprPtr& x, const ExprPtr& r) {
  if (r->is_zero()) return {};
  auto [keep, projected] =
      indices_of_product(products_of(summands_of(f).front()).front());
  for (const auto& i : projected)
    keep.insert(i.proto_indices().begin(), i.proto_indices().end());
  keep.merge(indices_of_leaf(x));
  return validity_violation(r, keep);
}

void check_rule(const ExprPtr& f, const ExprPtr& x, const ExprPtr& expected,
                const Data& data) {
  INFO(text_of(f) << " by " << text_of(x));
  const auto r = sequant::deriv::diff(f, {x});
  INFO("= " << text_of(r));
  if (expected) REQUIRE_THAT(r, EquivalentTo(expected));
  const auto c = compare_tangents(f, r, x, data, {.count = 256, .seed = 1});
  INFO(c.difference);
  REQUIRE(c.difference.empty());
  const auto violation = validity(f, x, r);
  INFO(violation);
  REQUIRE(violation.empty());
}

void check_rule(std::wstring_view f, std::wstring_view x,
                std::wstring_view expected, const Data& data) {
  check_rule(deserialize(f), deserialize(x),
             expected.empty() ? nullptr : deserialize(expected), data);
}

void require_same(const ExprPtr& a, const ExprPtr& b,
                  const container::set<Index>& free, const Data& data) {
  INFO(text_of(a) << "\n  vs\n" << text_of(b));
  const auto c = compare_arrays(a, b, free, data);
  INFO(c.difference);
  REQUIRE(c.difference.empty());
  REQUIRE(c.nonzero != 0);
}

}  // namespace sequant::deriv::test
