#ifndef SEQUANT_ASY_COST_HPP
#define SEQUANT_ASY_COST_HPP

#include <SeQuant/core/container.hpp>
#include <boost/rational.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>
#include <sstream>

namespace sequant {

///
/// Represents an symbolic asymptotic cost in terms of active_occupied
/// and the rest orbitals.
/// eg.
///     - AsyCost{2,4} implies scaling of $O^2V^4$. In other words, the cost
///       scales by the second power in the number of active_occupied orbitals
///       and the fourth power in the number of the rest orbitals.
///     - AsyCost{2, 4, boost::rational<int>{1,2}} implies the same scaling as
///       above except the numeric value obtained by substituting $O$ and $V$
///       numbers is then halved.
class AsyCost {
 private:
  class AsyCostEntry {
    size_t occ_; // power of active_occupied
    size_t virt_; // power of the rest orbitals
    mutable boost::rational<int> count_; // count of this asymptotic symbol

   public:
    template <typename Os, typename IntType>
    static Os &stream_out_rational(Os &os, boost::rational<IntType> r) {
      os << r.numerator();
      if (r.denominator() != IntType{1}) {
        os << '/';
        os << r.denominator();
      }
      return os;
    }

    static AsyCostEntry const &max();

    static AsyCostEntry const &zero();

    AsyCostEntry(size_t nocc, size_t nvirt, boost::rational<int> count);

    AsyCostEntry(AsyCostEntry const &) = default;

    AsyCostEntry(AsyCostEntry &&) = default;

    AsyCostEntry &operator=(AsyCostEntry const &) = default;

    AsyCostEntry &operator=(AsyCostEntry &&) = default;

    size_t occ() const;

    size_t virt() const;

    boost::rational<int> count() const;

    void set_count(boost::rational<int> n) const;

    bool operator<(AsyCostEntry const &rhs) const;

    bool operator==(AsyCostEntry const &rhs) const;

    bool operator!=(AsyCostEntry const &rhs) const;

    template <typename String_t>
    String_t text() const {
      auto oss = std::basic_ostringstream<typename String_t::value_type>{};

      if (*this == AsyCostEntry::max()) {
        oss << "max";
      } else if (*this == AsyCostEntry::zero()) {
        oss << "zero";
      } else {
        auto abs_c = boost::abs(count_);
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

    template <typename String_t>
    String_t to_latex() const {
      auto oss = std::basic_ostringstream<typename String_t::value_type>{};

      if (*this == AsyCostEntry::max()) {
        oss << "\\texttt{max}";
      } else if (*this == AsyCostEntry::zero()) {
        oss << "\\texttt{zero}";
      } else {
        auto abs_c = boost::abs(count_);
        oss << (count_ < abs_c ? "- " : "");
        bool frac_mode = abs_c.denominator() != 1;
        if (!frac_mode && (abs_c != 1))
          oss << count_.numerator();
        if (frac_mode) {
          oss << "\\frac{" << std::abs(count_.numerator()) << "}{"
              << count_.denominator() << "}";
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
  };

 private:
  container::set<AsyCostEntry> cost_;

  AsyCost(AsyCostEntry);

 public:
  ///
  /// \return The infinitely scaling cost.
  static AsyCost const &max();

  ///
  /// \return The zero cost.
  static AsyCost const &zero();

  ///
  /// Default construct to zero cost.
  AsyCost();

  ///
  /// \param nocc Asymptotic scaling exponent in the active occupied orbitals.
  /// \param nrest Asymptotic scaling exponent in the rest orbitals.
  /// \param count Rational number of times this cost repeats.
  AsyCost(size_t nocc, size_t nrest, boost::rational<int> count = 1);

  AsyCost(AsyCost const &) = default;

  AsyCost(AsyCost &&) = default;

  AsyCost &operator=(AsyCost const &) = default;

  AsyCost &operator=(AsyCost &&) = default;

  ///
  /// \param nocc Substitute $O$ by nocc.
  /// \param nvirt Substitute $V$ by nvirt.
  /// \return Scaled asymptotic cost.
  [[nodiscard]] boost::rational<long long int> ops(unsigned short nocc,
                                                   unsigned short nvirt) const;

  template <typename String_t>
  String_t to_latex() const {
    auto oss = std::basic_ostringstream<typename String_t::value_type>{};
    // oss << "{";
    if (cost_.empty())
      oss << 0;
    else {
      // stream out in reverse so that more expensive terms appear first
      auto rev = ranges::views::reverse(cost_);
      oss << ranges::front(rev).to_latex<String_t>();
      if (cost_.size() > 1)
        for (auto &&c : ranges::views::tail(rev)) {
          oss << (c.count() > 0 ? " + " : " ") << c.to_latex<String_t>();
        }
    }
    // oss << "}";
    auto str = oss.str();
    return oss.str();
  }

  friend AsyCost operator+(AsyCost const &lhs, AsyCost const &rhs);

  friend AsyCost operator*(AsyCost const &lhs, boost::rational<int> scale);

  friend bool operator<(AsyCost const &lhs, AsyCost const &rhs);

  friend bool operator==(AsyCost const &lhs, AsyCost const &rhs);

  template <typename Os>
  friend Os &operator<<(Os &os, AsyCost const &cost);
};

AsyCost operator+(AsyCost const &lhs, AsyCost const &rhs);

AsyCost operator-(AsyCost const &lhs, AsyCost const &rhs);

AsyCost operator*(AsyCost const &cost, boost::rational<int> scale);

AsyCost operator*(boost::rational<int> scale, AsyCost const &cost);

AsyCost operator/(AsyCost const &cost, boost::rational<int> scale);

bool operator==(AsyCost const &lhs, AsyCost const &rhs);

bool operator!=(AsyCost const &lhs, AsyCost const &rhs);

bool operator<(AsyCost const &lhs, AsyCost const &rhs);

bool operator>(AsyCost const &lhs, AsyCost const &rhs);

template <typename Os>
Os &operator<<(Os &os, AsyCost const &cost) {
  if (cost == AsyCost::zero()) {
    os << 0;
    return os;
  }
  // stream out in reverse so that more expensive terms appear first
  auto rev = ranges::views::reverse(cost.cost_);
  os << ranges::front(rev).text<std::basic_string<typename Os::char_type>>();

  if (cost.cost_.size() > 1)
    for (auto &&c : ranges::views::tail(rev))
      os << (c.count() > 0 ? " + " : " ")
         << c.text<std::basic_string<typename Os::char_type>>();

  return os;
}
}  // namespace sequant

#endif  // SEQUANT_ASY_COST_HPP
