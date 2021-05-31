#ifndef SEQUANT_ASY_COST_HPP
#define SEQUANT_ASY_COST_HPP

#include <SeQuant/core/container.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>
#include <sstream>

namespace sequant {

class AsyCost {
 private:
  class AsyCostEntry {
    size_t occ_;
    size_t virt_;
    mutable int count_ = 1;

   public:
    static AsyCostEntry max();

    AsyCostEntry(size_t nocc, size_t nvirt);

    AsyCostEntry(size_t nocc, size_t nvirt, int count);

    AsyCostEntry(AsyCostEntry const &) = default;

    AsyCostEntry(AsyCostEntry &&) = default;

    AsyCostEntry &operator=(AsyCostEntry const &) = default;

    AsyCostEntry &operator=(AsyCostEntry &&) = default;

    size_t occ() const;

    size_t virt() const;

    int count() const;

    void set_count(int n) const;

    bool operator<(AsyCostEntry const &rhs) const;

    bool operator==(AsyCostEntry const &rhs) const;

    template <typename String_t>
    String_t text() const {
      auto oss = std::basic_ostringstream<typename String_t::value_type>{};
      if (*this == AsyCostEntry::max()) {
        oss << "max";
      } else {
        if (count_ < 0) {
          if (count_ == -1)
            oss << "- ";
          else
            oss << "- " << std::abs(count_) << "*";
        } else if (count_ > 1) {
          oss << count_ << "*";
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
        oss << "max";
      } else {
        if (count_ < 0) oss << "- ";
        oss << "{";
        if (std::abs(count_) != 1) oss << std::abs(count_) << " ";
        oss << (occ_ > 0 ? "O" : "");
        if (occ_ > 1) {
          oss << "^{" << occ_ << "}";
        }
        oss << (virt_ > 0 ? "V" : "");
        if (virt_ > 1) {
          oss << "^{" << virt_ << "}";
        }
        oss << "}";
      }
      return oss.str();
    }
  };

  sequant::container::set<AsyCostEntry> cost_;

  AsyCost(AsyCostEntry);

 public:
  static AsyCost const &max();

  static AsyCost const &zero();

  AsyCost(size_t nocc, size_t nvirt);

  AsyCost(AsyCost const &) = default;

  AsyCost(AsyCost &&) = default;

  AsyCost &operator=(AsyCost const &) = default;

  AsyCost &operator=(AsyCost &&) = default;

  signed long long ops(unsigned short nocc, unsigned short int nvirt) const;

  AsyCost operator+(AsyCost const &rhs) const;

  AsyCost operator-(AsyCost const &rhs) const;

  AsyCost &operator+=(AsyCost const &rhs);

  AsyCost &operator-=(AsyCost const &rhs);

  bool operator==(AsyCost const &rhs) const;

  bool operator!=(AsyCost const &rhs) const;

  bool operator<(AsyCost const &rhs) const;

  bool operator>(AsyCost const &rhs) const;

  template <typename String_t>
  String_t to_latex() const {
    auto oss = std::basic_ostringstream<typename String_t::value_type>{};
    // oss << "{";
    if (cost_.empty())
      oss << 0;
    else {
      oss << ranges::front(cost_).to_latex<String_t>();
      if (cost_.size() > 1)
        for (auto &&c : ranges::views::tail(cost_)) {
          oss << (c.count() > 0 ? " + " : " ") << c.to_latex<String_t>();
        }
    }
    // oss << "}";
    return oss.str();
  }

  template <typename Os>
  friend Os &operator<<(Os &os, AsyCost const &cost);
};

template <typename Os>
Os &operator<<(Os &os, AsyCost const &cost) {
  if (cost.cost_.empty()) {
    os << 0;
    return os;
  }
  os << ranges::front(cost.cost_)
            .text<std::basic_string<typename Os::char_type>>();

  if (cost.cost_.size() > 1)
    for (auto &&c : ranges::views::tail(cost.cost_))
      os << (c.count() > 0 ? " + " : " ")
         << c.text<std::basic_string<typename Os::char_type>>();

  return os;
}
}  // namespace sequant

#endif  // SEQUANT_ASY_COST_HPP
