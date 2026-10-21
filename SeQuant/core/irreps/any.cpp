#include <SeQuant/core/irreps/any.hpp>
#include <SeQuant/core/utility/macros.hpp>

namespace sequant {

// AnyIrrep

AnyIrrep::AnyIrrep(const AnyIrrep &other)
    : m_impl(other.m_impl ? other.m_impl->clone() : nullptr) {}

AnyIrrep::AnyIrrep(AnyIrrep &&other) noexcept = default;

AnyIrrep &AnyIrrep::operator=(const AnyIrrep &other) {
  if (this != &other) {
    m_impl = other.m_impl ? other.m_impl->clone() : nullptr;
  }
  return *this;
}

AnyIrrep &AnyIrrep::operator=(AnyIrrep &&other) noexcept = default;

AnyIrrep::~AnyIrrep() = default;

std::string AnyIrrep::to_string() const {
  SEQUANT_ASSERT(m_impl);
  return m_impl->to_string();
}

bool operator==(const AnyIrrep &lhs, const AnyIrrep &rhs) {
  SEQUANT_ASSERT(lhs.m_impl && rhs.m_impl);
  return lhs.m_impl->equals(*rhs.m_impl);
}

// AnyIrrepSet

AnyIrrepSet::AnyIrrepSet() = default;

AnyIrrepSet::AnyIrrepSet(const AnyIrrepSet &other)
    : m_impl(other.m_impl ? other.m_impl->clone() : nullptr) {}

AnyIrrepSet::AnyIrrepSet(AnyIrrepSet &&other) noexcept = default;

AnyIrrepSet &AnyIrrepSet::operator=(const AnyIrrepSet &other) {
  if (this != &other) {
    m_impl = other.m_impl ? other.m_impl->clone() : nullptr;
  }
  return *this;
}

AnyIrrepSet &AnyIrrepSet::operator=(AnyIrrepSet &&other) noexcept = default;

AnyIrrepSet::~AnyIrrepSet() = default;

bool AnyIrrepSet::is_any() const { return m_impl && m_impl->is_any(); }

bool AnyIrrepSet::empty() const { return !m_impl || m_impl->empty(); }

std::size_t AnyIrrepSet::size() const { return m_impl ? m_impl->size() : 0; }

AnyIrrepSet &AnyIrrepSet::unite(const AnyIrrepSet &other) {
  if (!other.m_impl) {
    return *this;
  }
  if (!m_impl) {
    m_impl = other.m_impl->clone();
    return *this;
  }

  if (!m_impl->unite(*other.m_impl)) {
    throw Exception("Can't unite IRREP sets of different IRREP types: " +
                    to_string() + " and " + other.to_string());
  }

  return *this;
}

std::string AnyIrrepSet::to_string() const {
  return m_impl ? m_impl->to_string() : "{}";
}

void AnyIrrepSet::for_each(
    const std::function<void(const AnyIrrep &)> &f) const {
  if (m_impl) {
    m_impl->for_each(f);
  }
}

bool operator==(const AnyIrrepSet &lhs, const AnyIrrepSet &rhs) {
  if (!lhs.m_impl || !rhs.m_impl) {
    // At least one of the sets is an untyped, empty set
    return lhs.empty() && rhs.empty();
  }

  return lhs.m_impl->equals(*rhs.m_impl);
}

// AnyIrrepContext

AnyIrrepContext::AnyIrrepContext(const AnyIrrepContext &other)
    : m_impl(other.m_impl ? other.m_impl->clone() : nullptr) {}

AnyIrrepContext::AnyIrrepContext(AnyIrrepContext &&other) noexcept = default;

AnyIrrepContext &AnyIrrepContext::operator=(const AnyIrrepContext &other) {
  if (this != &other) {
    m_impl = other.m_impl ? other.m_impl->clone() : nullptr;
  }
  return *this;
}

AnyIrrepContext &AnyIrrepContext::operator=(AnyIrrepContext &&other) noexcept =
    default;

AnyIrrepContext::~AnyIrrepContext() = default;

AnyIrrepSet AnyIrrepContext::irreps_of(const Tensor &tensor) const {
  return m_impl->irreps_of(tensor);
}

AnyIrrepSet AnyIrrepContext::totally_symmetric() const {
  return m_impl->totally_symmetric();
}

AnyIrrepSet AnyIrrepContext::empty_set() const { return m_impl->empty_set(); }

AnyIrrepSet AnyIrrepContext::product(const AnyIrrepSet &lhs,
                                     const AnyIrrepSet &rhs) const {
  return m_impl->product(lhs, rhs);
}

AnyIrrepSet AnyIrrepContext::conjugate(const AnyIrrepSet &set) const {
  return m_impl->conjugate(set);
}

}  // namespace sequant
