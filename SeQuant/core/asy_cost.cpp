#include "asy_cost.hpp"
#include "container.hpp"
#include "meta/meta.hpp"
#include "rational.hpp"
#include "wstring.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include <cmath>
#include <limits>

#include <range/v3/iterator.hpp>
#include <range/v3/view.hpp>

namespace sequant {

AsyCost::AsyCostEntry::AsyCostEntry(size_t nocc, size_t nvirt, rational count)
    : occ_{nocc}, virt_{nvirt}, count_{count} {
  if (count_ == 0 || (occ_ == 0 && virt_ == 0)) {
    occ_ = 0;
    virt_ = 0;
    count_ = 0;
  }
}

AsyCost::AsyCostEntry AsyCost::AsyCostEntry::max() {
  return AsyCostEntry{std::numeric_limits<size_t>::max(),
                      std::numeric_limits<size_t>::max(),
                      std::numeric_limits<intmax_t>::max()};
}

AsyCost::AsyCostEntry const &AsyCost::AsyCostEntry::zero() {
  static AsyCostEntry const zero_cost = AsyCostEntry{0, 0, 0};
  return zero_cost;
}

size_t AsyCost::AsyCostEntry::occ() const { return occ_; }

size_t AsyCost::AsyCostEntry::virt() const { return virt_; }

rational AsyCost::AsyCostEntry::count() const { return count_; }

void AsyCost::AsyCostEntry::set_count(rational n) const { count_ = n; }

std::ostream &AsyCost::AsyCostEntry::stream_out_rational(std::ostream &os,
                                                         rational const &r) {
  os << numerator(r);
  if (denominator(r) != rational{1}) {
    os << '/';
    os << denominator(r);
  }
  return os;
}

std::string AsyCost::AsyCostEntry::text() const {
  auto oss = std::ostringstream{};

  if (*this == AsyCostEntry::max()) {
    oss << "max";
  } else if (*this == AsyCostEntry::zero()) {
    oss << "zero";
  } else {
    auto abs_c = abs(count_);
    oss << (count_ < abs_c ? "- " : "");
    if (abs_c == 1) {
      // do nothing
    } else {
      AsyCostEntry::stream_out_rational(oss, abs_c);
      oss << "*";
    }
    oss << (occ_ > 0 ? "O" : "");
    if (occ_ > 1) oss << "^" << occ_;

    oss << (virt_ > 0 ? "V" : "");
    if (virt_ > 1) oss << "^" << virt_;
  }

  return oss.str();
}

std::string AsyCost::AsyCostEntry::to_latex() const {
  auto oss = std::ostringstream{};

  if (*this == AsyCostEntry::max()) {
    oss << "\\texttt{max}";
  } else if (*this == AsyCostEntry::zero()) {
    oss << "\\texttt{zero}";
  } else {
    auto abs_c = abs(count_);
    oss << (count_ < abs_c ? "- " : "");
    bool frac_mode = abs(denominator(count_)) != 1;
    if (!frac_mode && (abs_c != 1)) oss << numerator(count_);
    if (frac_mode) {
      oss << "\\frac{"               //
          << abs(numerator(count_))  //
          << "}{"                    //
          << denominator(count_) << "}";
    }
    oss << (occ_ > 0 ? "O" : "");
    if (occ_ > 1) {
      oss << "^{" << occ_ << "}";
    }
    oss << (virt_ > 0 ? "V" : "");
    if (virt_ > 1) {
      oss << "^{" << virt_ << "}";
    }
  }
  return oss.str();
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

AsyCost::AsyCost(rational count, size_t nocc, size_t nvirt)
    : AsyCost{AsyCostEntry{nocc, nvirt, count}} {}

AsyCost::AsyCost(size_t nocc, size_t nvirt) : AsyCost{1, nocc, nvirt} {}

AsyCost::AsyCost(std::pair<size_t, size_t> const &ov)
    : AsyCost{ov.first, ov.second} {}

double AsyCost::ops(size_t nocc, size_t nvirt) const {
  double total = 0;
  for (auto &&c : cost_) {
    double temp = 1;
    temp *= std::pow(nocc, c.occ());
    temp *= std::pow(nvirt, c.virt());
    total += temp > 1 ? boost::numeric_cast<double>(c.count()) * temp : 0;
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

AsyCost& AsyCost::operator+=(AsyCost const& other) {
  *this = *this + other;
  return *this;
}

AsyCost& AsyCost::operator-=(AsyCost const& other) {
  *this = *this - other;
  return *this;
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

AsyCost operator*(AsyCost const &cost, rational scale) {
  auto ac = cost;
  for (auto &c : ac.cost_) c.set_count(c.count() * scale);
  return ac;
}

AsyCost operator*(rational scale, AsyCost const &cost) { return cost * scale; }

AsyCost operator/(AsyCost const &cost, rational scale) {
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

bool operator<=(AsyCost const& lhs, AsyCost const& rhs) {
  return lhs < rhs || lhs == rhs;
}

bool operator>=(AsyCost const& lhs, AsyCost const& rhs) {
  return lhs > rhs || lhs == rhs;
}

std::wstring AsyCost::to_latex() const {
  auto oss = std::wostringstream{};
  if (cost_.empty())
    oss << 0;
  else {
    //
    // stream out in reverse so that more expensive terms appear first
    auto rev = ranges::views::reverse(cost_);
    oss << sequant::to_wstring(ranges::front(rev).to_latex());
    for (auto &&c : ranges::views::tail(rev))
      oss << L" + " << sequant::to_wstring(c.to_latex());
  }
  return oss.str();
}

std::string AsyCost::text() const {
  if (*this == AsyCost::zero()) return "0";

  std::ostringstream oss{};
  // reversed so that more expensive terms appear first
  auto rev = ranges::views::reverse(cost_);
  oss << ranges::front(rev).text();
  for (auto &&c : ranges::views::tail(rev)) oss << " + " << c.text();

  return oss.str();
}

}  // namespace sequant
