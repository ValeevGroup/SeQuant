#include "asy_cost.hpp"
#include <iostream>

namespace sequant {

AsyCost::AsyCostEntry::AsyCostEntry(size_t nocc, size_t nvirt)
    : occ_{nocc}, virt_{nvirt} {}

AsyCost::AsyCostEntry::AsyCostEntry(size_t nocc, size_t nvirt, int count)
    : occ_{nocc}, virt_{nvirt}, count_{count} {}

AsyCost::AsyCostEntry AsyCost::AsyCostEntry::max() {
  return AsyCostEntry{std::numeric_limits<size_t>::max(),
                      std::numeric_limits<size_t>::max(),
                      std::numeric_limits<int>::max()};
}

size_t AsyCost::AsyCostEntry::occ() const { return occ_; }

size_t AsyCost::AsyCostEntry::virt() const { return virt_; }

int AsyCost::AsyCostEntry::count() const { return count_; }

void AsyCost::AsyCostEntry::set_count(int n) const { count_ = n; }

bool AsyCost::AsyCostEntry::operator<(const AsyCost::AsyCostEntry &rhs) const {
  return virt() < rhs.virt() || (virt() == rhs.virt() && occ() < rhs.occ());
}

bool AsyCost::AsyCostEntry::operator==(const AsyCost::AsyCostEntry &rhs) const {
  return occ() == rhs.occ() && virt() == rhs.virt();
}

AsyCost::AsyCost(size_t nocc, size_t nvirt) {
  if (!(nocc == 0 && nvirt == 0)) cost_.emplace(AsyCostEntry{nocc, nvirt});
}

AsyCost::AsyCost(AsyCostEntry c) : AsyCost{0, 0} {
  cost_.emplace(std::move(c));
}

signed long long AsyCost::ops(unsigned short nocc, unsigned short nvirt) const {
  auto total = 0;
  for (auto &&c : cost_) {
    auto temp = 1;
    if (c.occ() > 0) temp *= static_cast<int>(std::pow(nocc, c.occ()));
    if (c.virt() > 0) temp *= static_cast<int>(std::pow(nvirt, c.virt()));
    // 2 * c.count() because matrix operations flops
    total += temp > 1 ? 2 * c.count() * temp : 0;
  }
  return total;
}

AsyCost AsyCost::operator+(const AsyCost &rhs) const {
  auto sum = AsyCost{*this};
  auto &data = sum.cost_;
  for (auto const &c : rhs.cost_) {
    if (auto found = data.find(c); found != data.end()) {
      found->set_count(found->count() + c.count());
      if (found->count() == 0) data.erase(found);
    } else {
      data.emplace(c);
    }
  }
  return sum;
}

AsyCost AsyCost::operator-(const AsyCost &rhs) const {
  auto diff = AsyCost{*this};
  auto &data = diff.cost_;
  for (auto const &c : rhs.cost_) {
    if (auto found = data.find(c); found != data.end()) {
      found->set_count(found->count() - c.count());
      if (found->count() == 0) data.erase(found);
    } else {
      data.emplace(AsyCostEntry{c.occ(), c.virt(), -c.count()});
    }
  }
  return diff;
}

AsyCost &AsyCost::operator+=(const AsyCost &rhs) {
  *this = *this + rhs;
  return *this;
}

AsyCost &AsyCost::operator-=(const AsyCost &rhs) {
  *this = *this - rhs;
  return *this;
}

bool AsyCost::operator<(const AsyCost &rhs) const {
  using ranges::views::reverse;
  using ranges::views::zip;

  for (auto &&[c1, c2] : reverse(zip(cost_, rhs.cost_))) {
    if (c1 < c2)
      return true;
    else if (c1 == c2) {
      if (c1.count() < c2.count())
        return true;
      else if (c1.count() > c2.count())
        return false;
    } else
      return false;
  }

  return cost_.size() < rhs.cost_.size();
}

bool AsyCost::operator==(const AsyCost &rhs) const {
  return cost_.size() == rhs.cost_.size() && !(*this < rhs || rhs < *this);
}

bool AsyCost::operator!=(const AsyCost &rhs) const { return !(*this == rhs); }

bool AsyCost::operator>(const AsyCost &rhs) const {
  return !(*this < rhs || *this == rhs);
}

AsyCost const &AsyCost::max() {
  static const AsyCost max = AsyCost{AsyCostEntry::max()};
  return max;
}

AsyCost const &AsyCost::zero() {
  static const AsyCost zero = AsyCost{0, 0};
  return zero;
}

}  // namespace sequant
