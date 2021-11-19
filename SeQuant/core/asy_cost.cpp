#include "asy_cost.hpp"

namespace sequant {

AsyCostEntry::AsyCostEntry(size_t nocc, size_t nvirt,
                           boost::rational<int> count)
    : occ_{nocc}, virt_{nvirt}, count_{count} {
  if (count_ == 0 || (occ_ == 0 && virt_ == 0)) {
    occ_ = 0;
    virt_ = 0;
    occ_ = 0;
  }
}

AsyCostEntry const &AsyCostEntry::max() {
  static AsyCostEntry const max_cost = AsyCostEntry{
      std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(),
      boost::rational<int>{std::numeric_limits<int>::max(),
                           std::numeric_limits<int>::max()}};
  return max_cost;
}

AsyCostEntry const &AsyCostEntry::zero() {
  static AsyCostEntry const zero_cost = AsyCostEntry{0, 0, 0};
  return zero_cost;
}

size_t AsyCostEntry::occ() const { return occ_; }

size_t AsyCostEntry::virt() const { return virt_; }

boost::rational<int> AsyCostEntry::count() const { return count_; }

void AsyCostEntry::set_count(boost::rational<int> n) const { count_ = n; }

bool AsyCostEntry::operator<(const AsyCostEntry &rhs) const {
  return virt() < rhs.virt() || (virt() == rhs.virt() && occ() < rhs.occ());
}

bool AsyCostEntry::operator==(const AsyCostEntry &rhs) const {
  return occ() == rhs.occ() && virt() == rhs.virt();
}

bool AsyCostEntry::operator!=(const AsyCostEntry &rhs) const {
  return !(*this == rhs);
}

AsyCost::AsyCost(AsyCostEntry c) {
  if (c != AsyCostEntry::zero()) cost_.emplace(c);
}

AsyCost::AsyCost(size_t nocc, size_t nvirt, boost::rational<int> count)
    : AsyCost{AsyCostEntry{nocc, nvirt, count}} {}

boost::rational<long long int> AsyCost::ops(unsigned short nocc,
                                            unsigned short nvirt) const {
  boost::rational<long long int> total = 0;
  for (auto &&c : cost_) {
    auto temp = 1;
    if (c.occ() > 0) temp *= static_cast<int>(std::pow(nocc, c.occ()));
    if (c.virt() > 0) temp *= static_cast<int>(std::pow(nvirt, c.virt()));
    // 2 * c.count() because matrix operations flops
    total += temp > 1 ? 2 * boost::rational_cast<int>(c.count()) * temp : 0;
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
  static const AsyCost zero = AsyCost{AsyCostEntry::zero()};
  return zero;
}

}  // namespace sequant
