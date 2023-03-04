#include "asy_cost.hpp"

namespace sequant {

AsyCost::AsyCostEntry::AsyCostEntry(size_t nocc, size_t nvirt,
                                    boost::rational<int> count)
    : occ_{nocc}, virt_{nvirt}, count_{count} {
  if (count_ == 0 || (occ_ == 0 && virt_ == 0)) {
    occ_ = 0;
    virt_ = 0;
    occ_ = 0;
  }
}

AsyCost::AsyCostEntry const &AsyCost::AsyCostEntry::max() {
  static AsyCostEntry const max_cost = AsyCostEntry{
      std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(),
      boost::rational<int>{std::numeric_limits<int>::max(),
                           std::numeric_limits<int>::max()}};
  return max_cost;
}

AsyCost::AsyCostEntry const &AsyCost::AsyCostEntry::zero() {
  static AsyCostEntry const zero_cost = AsyCostEntry{0, 0, 0};
  return zero_cost;
}

size_t AsyCost::AsyCostEntry::occ() const { return occ_; }

size_t AsyCost::AsyCostEntry::virt() const { return virt_; }

boost::rational<int> AsyCost::AsyCostEntry::count() const { return count_; }

void AsyCost::AsyCostEntry::set_count(boost::rational<int> n) const {
  count_ = n;
}

bool AsyCost::AsyCostEntry::operator<(const AsyCost::AsyCostEntry &rhs) const {
  return virt() < rhs.virt() || (virt() == rhs.virt() && occ() < rhs.occ());
}

bool AsyCost::AsyCostEntry::operator==(const AsyCost::AsyCostEntry &rhs) const {
  return occ() == rhs.occ() && virt() == rhs.virt();
}

bool AsyCost::AsyCostEntry::operator!=(const AsyCost::AsyCostEntry &rhs) const {
  return !(*this == rhs);
}

AsyCost::AsyCost(AsyCostEntry c) {
  if (c != AsyCostEntry::zero()) cost_.emplace(c);
}

AsyCost::AsyCost() : AsyCost{AsyCostEntry::zero()} {}

AsyCost::AsyCost(boost::rational<int> count, size_t nocc, size_t nvirt)
    : AsyCost{AsyCostEntry{nocc, nvirt, count}} {}

AsyCost::AsyCost(size_t nocc, size_t nvirt) : AsyCost{1, nocc, nvirt} {}

double AsyCost::ops(size_t nocc, size_t nvirt) const {
  double total = 0;
  for (auto &&c : cost_) {
    double temp = 1;
    temp *= std::pow(nocc, c.occ());
    temp *= std::pow(nvirt, c.virt());
    total += temp > 1 ? boost::rational_cast<double>(c.count()) * temp : 0;
  }
  return total;
}

AsyCost const &AsyCost::max() {
  static const AsyCost max = AsyCost{AsyCostEntry::max()};
  return max;
}

AsyCost const &AsyCost::zero() {
  static const AsyCost zero = AsyCost{AsyCostEntry::zero()};
  return zero;
}

AsyCost operator+(AsyCost const &lhs, AsyCost const &rhs) {
  auto sum = lhs;
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

AsyCost operator-(AsyCost const &lhs, AsyCost const &rhs) {
  return lhs + (-1 * rhs);
}

AsyCost operator*(AsyCost const &cost, boost::rational<int> scale) {
  auto ac = cost;
  for (auto &c : ac.cost_) c.set_count(c.count() * scale);
  return ac;
}

AsyCost operator*(boost::rational<int> scale, AsyCost const &cost) {
  return cost * scale;
}

AsyCost operator/(AsyCost const &cost, boost::rational<int> scale) {
  return cost * (1 / scale);
}

bool operator<(AsyCost const &lhs, AsyCost const &rhs) {
  using ranges::views::reverse;
  using ranges::views::zip;

  for (auto &&[c1, c2] : reverse(zip(lhs.cost_, rhs.cost_))) {
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

  return lhs.cost_.size() < rhs.cost_.size();
}

bool operator==(AsyCost const &lhs, AsyCost const &rhs) {
  return lhs.cost_.size() == rhs.cost_.size() && !(lhs < rhs || rhs < lhs);
}

bool operator!=(AsyCost const &lhs, AsyCost const &rhs) {
  return !(lhs == rhs);
}

bool operator>(AsyCost const &lhs, AsyCost const &rhs) {
  return !(lhs < rhs || lhs == rhs);
}

}  // namespace sequant
