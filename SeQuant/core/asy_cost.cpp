#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include <range/v3/range/operations.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/zip.hpp>

#include <cmath>
#include <limits>

namespace sequant {

AsyCost::AsyCostEntry::AsyCostEntry()
    : exponents_{}, count_{0}, is_max_{false} {}

AsyCost::AsyCostEntry::AsyCostEntry(ExponentMap exponents, rational count)
    : exponents_{std::move(exponents)}, count_{count}, is_max_{false} {
  for (auto it = exponents_.begin(); it != exponents_.end();) {
    if (it->second == 0)
      it = exponents_.erase(it);
    else
      ++it;
  }
  if (count_ == 0 || exponents_.empty()) {
    exponents_.clear();
    count_ = 0;
  }
}

AsyCost::AsyCostEntry AsyCost::AsyCostEntry::max() {
  AsyCostEntry e;
  e.is_max_ = true;
  e.count_ = std::numeric_limits<intmax_t>::max();
  return e;
}

AsyCost::AsyCostEntry const &AsyCost::AsyCostEntry::zero() {
  static AsyCostEntry const zero_cost{};
  return zero_cost;
}

AsyCost::ExponentMap const &AsyCost::AsyCostEntry::exponents() const {
  return exponents_;
}

rational AsyCost::AsyCostEntry::count() const { return count_; }

void AsyCost::AsyCostEntry::set_count(rational n) const { count_ = n; }

bool AsyCost::AsyCostEntry::is_zero() const {
  return !is_max_ && exponents_.empty() && count_ == 0;
}

bool AsyCost::AsyCostEntry::is_max() const { return is_max_; }

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

  if (is_max_) {
    oss << "max";
  } else if (is_zero()) {
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
    // Spaces print in IndexSpace order (by attr, then base_key).
    for (auto const &[space, exp] : exponents_) {
      if (exp == 0) continue;
      oss << toUtf8(space.base_key());
      if (exp > 1) oss << "^" << exp;
    }
  }

  return oss.str();
}

std::string AsyCost::AsyCostEntry::to_latex() const {
  auto oss = std::ostringstream{};

  if (is_max_) {
    oss << "\\texttt{max}";
  } else if (is_zero()) {
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
    for (auto const &[space, exp] : exponents_) {
      if (exp == 0) continue;
      oss << toUtf8(space.base_key());
      if (exp > 1) {
        oss << "^{" << exp << "}";
      }
    }
  }
  return oss.str();
}

bool AsyCost::AsyCostEntry::operator<(const AsyCost::AsyCostEntry &rhs) const {
  // The max() sentinel is the greatest entry.
  if (is_max_ != rhs.is_max_) return !is_max_;

  // Order by exponents, space by space, from the greatest IndexSpace (highest
  // priority) down to the least, using IndexSpace's own ordering. The first
  // space whose exponents differ decides: a larger exponent on a higher-
  // priority space makes the entry larger. Entries with identical exponents
  // on every space are equivalent.
  auto exp_for = [](ExponentMap const &m, IndexSpace const &s) {
    auto it = m.find(s);
    return it == m.end() ? std::size_t{0} : it->second;
  };

  container::set<IndexSpace> spaces;
  for (auto const &kv : exponents_) spaces.insert(kv.first);
  for (auto const &kv : rhs.exponents_) spaces.insert(kv.first);

  for (auto it = spaces.rbegin(); it != spaces.rend(); ++it) {
    auto const l = exp_for(exponents_, *it);
    auto const r = exp_for(rhs.exponents_, *it);
    if (l != r) return l < r;
  }
  return false;
}

bool AsyCost::AsyCostEntry::operator==(const AsyCost::AsyCostEntry &rhs) const {
  return is_max_ == rhs.is_max_ && exponents_ == rhs.exponents_;
}

bool AsyCost::AsyCostEntry::operator!=(const AsyCost::AsyCostEntry &rhs) const {
  return !(*this == rhs);
}

AsyCost::AsyCost(AsyCostEntry c) {
  if (!c.is_zero()) cost_.emplace(std::move(c));
}

AsyCost::AsyCost() : AsyCost{AsyCostEntry::zero()} {}

AsyCost::AsyCost(ExponentMap exponents, rational count)
    : AsyCost{AsyCostEntry{std::move(exponents), count}} {}

double AsyCost::ops(ExtentMap const &extents) const {
  double total = 0;
  for (auto const &c : cost_) {
    if (c.is_max()) return std::numeric_limits<double>::infinity();
    double temp = 1;
    for (auto const &[space, exp] : c.exponents()) {
      auto it = extents.find(space);
      std::size_t const extent = (it != extents.end()) ? it->second : 1;
      temp *= std::pow(static_cast<double>(extent), static_cast<double>(exp));
    }
    total += boost::numeric_cast<double>(c.count()) * temp;
  }
  return total;
}

AsyCost const &AsyCost::max() {
  static const AsyCost m = AsyCost{AsyCostEntry::max()};
  return m;
}

AsyCost const &AsyCost::zero() {
  static const AsyCost z = AsyCost{AsyCostEntry::zero()};
  return z;
}

AsyCost &AsyCost::operator+=(AsyCost const &other) {
  *this = *this + other;
  return *this;
}

AsyCost &AsyCost::operator-=(AsyCost const &other) {
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

bool operator<=(AsyCost const &lhs, AsyCost const &rhs) {
  return lhs < rhs || lhs == rhs;
}

bool operator>=(AsyCost const &lhs, AsyCost const &rhs) {
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
    oss << toUtf16(ranges::front(rev).to_latex());
    for (auto &&c : ranges::views::tail(rev))
      oss << L" + " << toUtf16(c.to_latex());
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
