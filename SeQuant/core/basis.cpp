#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <unordered_set>

namespace sequant {

IndexBasis::IndexBasis(IndexSpace space, optional_instance basis_instance,
                       std::wstring name, std::optional<std::size_t> extent,
                       IndexBasisMetric metric, std::optional<Field> field)
    : space_(std::move(space)),
      basis_instance_(basis_instance),
      name_(std::move(name)),
      extent_(extent),
      metric_(metric),
      field_(field) {
  SEQUANT_ASSERT(name_.empty() || basis_instance_);
  SEQUANT_ASSERT(!name_.empty() ||
                 (!extent_ && metric_ == IndexBasisMetric::Unit && !field_));
}

std::wstring IndexBasis::instance_suffix() const {
  return basis_instance_ ? L";" + std::to_wstring(*basis_instance_)
                         : std::wstring{};
}

bool same_instance(const IndexBasis& b1, const IndexBasis& b2) noexcept {
  return b1.space() == b2.space() && b1.basis_instance() == b2.basis_instance();
}

bool includes(const IndexBasis& basis, const IndexBasis& subbasis) {
  return includes(basis.space(), subbasis.space()) &&
         (!basis.has_basis_instance() ||
          (basis.basis_instance() == subbasis.basis_instance() &&
           basis.name() == subbasis.name()));
}

bool different_instances(const IndexBasis& b1, const IndexBasis& b2) {
  return b1.has_basis_instance() && b2.has_basis_instance() &&
         (b1.basis_instance() != b2.basis_instance() || b1.name() != b2.name());
}

namespace {

// all attributes of a basis, including those ignored by IndexBasis::operator==
struct FullEqual {
  bool operator()(const IndexBasis& b1, const IndexBasis& b2) const noexcept {
    const auto& s1 = b1.space();
    const auto& s2 = b2.space();
    return s1 == s2 && s1.dimension() == s2.dimension() &&
           s1.field() == s2.field() && b1 == b2 &&
           b1.own_extent() == b2.own_extent() && b1.metric() == b2.metric() &&
           b1.own_field() == b2.own_field();
  }
};

struct FullHash {
  std::size_t operator()(const IndexBasis& b) const {
    std::size_t result = hash_value(b);
    hash::combine(result, b.space().dimension());
    hash::combine(result, static_cast<int>(b.space().field()));
    hash::combine(result, b.own_extent().value_or(0));
    hash::combine(result, static_cast<int>(b.metric()));
    hash::combine(result,
                  b.own_field() ? static_cast<int>(*b.own_field()) : -1);
    return result;
  }
};

struct InternTable {
  std::shared_mutex mutex;
  // node-based, hence references to elements are stable
  std::unordered_set<IndexBasis, FullHash, FullEqual> bases;
};

InternTable& intern_table() {
  // intentionally leaked: Index objects of static storage duration may be
  // destroyed after (or compared with) any static of this translation unit
  static InternTable& table = *new InternTable;
  return table;
}

}  // namespace

const IndexBasis& intern(const IndexBasis& basis) {
  auto& table = intern_table();
  {
    std::shared_lock lock(table.mutex);
    if (auto it = table.bases.find(basis); it != table.bases.end()) return *it;
  }
  std::unique_lock lock(table.mutex);
  return *table.bases.insert(basis).first;
}

const IndexBasis& null_basis() {
  static const IndexBasis* const basis = &intern(IndexBasis{});
  return *basis;
}

}  // namespace sequant
