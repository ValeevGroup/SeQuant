#ifndef SEQUANT_ASY_COST_HPP
#define SEQUANT_ASY_COST_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>

#include <cstddef>
#include <string>

namespace sequant {

///
/// Represents a symbolic asymptotic cost as a polynomial in the sizes of
/// index spaces. A cost is a sum of terms, each of which is a rational
/// multiplier times a product of space sizes raised to integer powers.
/// Spaces are identified by `IndexSpace`; `AsyCost` orders and prints them
/// using `IndexSpace`'s own ordering and `base_key()`, and never consults an
/// `IndexSpaceRegistry`.
///
/// Examples (with `I`, `A` denoting two index spaces):
///   - `AsyCost({{I, 2}, {A, 4}})` represents $I^2 A^4$.
///   - `AsyCost({{I, 2}, {A, 4}}, rational{1, 2})` halves the numeric value
///     above when extents are substituted.
///
class AsyCost {
 public:
  using ExponentMap = container::map<IndexSpace, std::size_t>;
  using ExtentMap = container::map<IndexSpace, std::size_t>;

 private:
  class AsyCostEntry {
    ExponentMap exponents_;   // space -> power; zero exponents are not stored
    mutable rational count_;  // multiplier
    bool is_max_ = false;     // true for the AsyCost::max() sentinel

   public:
    static std::ostream &stream_out_rational(std::ostream &os,
                                             rational const &r);

    static AsyCostEntry max();

    static AsyCostEntry const &zero();

    AsyCostEntry();

    AsyCostEntry(ExponentMap exponents, rational count);

    AsyCostEntry(AsyCostEntry const &) = default;

    AsyCostEntry(AsyCostEntry &&) = default;

    AsyCostEntry &operator=(AsyCostEntry const &) = default;

    AsyCostEntry &operator=(AsyCostEntry &&) = default;

    [[nodiscard]] ExponentMap const &exponents() const;

    [[nodiscard]] rational count() const;

    void set_count(rational n) const;

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] bool is_max() const;

    bool operator<(AsyCostEntry const &rhs) const;

    bool operator==(AsyCostEntry const &rhs) const;

    bool operator!=(AsyCostEntry const &rhs) const;

    [[nodiscard]] std::string text() const;

    [[nodiscard]] std::string to_latex() const;
  };

 private:
  container::set<AsyCostEntry> cost_;

  explicit AsyCost(AsyCostEntry);

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
  /// \param exponents Map from index space to its exponent in this term.
  ///                  Zero exponents may be supplied; they are dropped.
  /// \param count     Rational multiplier; defaults to 1.
  ///
  explicit AsyCost(ExponentMap exponents, rational count = 1);

  AsyCost(AsyCost const &) = default;

  AsyCost(AsyCost &&) = default;

  AsyCost &operator=(AsyCost const &) = default;

  AsyCost &operator=(AsyCost &&) = default;

  ///
  /// Substitute each space in this cost by an extent and evaluate.
  /// \param extents Map from index space to extent (size). Any space appearing
  ///                in this cost but missing from `extents` is treated as
  ///                extent 1 (i.e. contributes a factor of 1).
  /// \return Numerical value of the cost.
  ///
  [[nodiscard]] double ops(ExtentMap const &extents) const;

  [[nodiscard]] std::wstring to_latex() const;

  [[nodiscard]] std::string text() const;

  AsyCost &operator+=(AsyCost const &);

  AsyCost &operator-=(AsyCost const &);

  friend AsyCost operator+(AsyCost const &lhs, AsyCost const &rhs);

  friend AsyCost operator*(AsyCost const &lhs, rational scale);

  friend AsyCost operator*(rational scale, AsyCost const &lhs);

  friend bool operator<(AsyCost const &lhs, AsyCost const &rhs);

  friend bool operator==(AsyCost const &lhs, AsyCost const &rhs);

  friend bool operator<=(AsyCost const &lhs, AsyCost const &rhs);

  friend bool operator>=(AsyCost const &lhs, AsyCost const &rhs);
};

AsyCost operator+(AsyCost const &lhs, AsyCost const &rhs);

AsyCost operator-(AsyCost const &lhs, AsyCost const &rhs);

AsyCost operator*(AsyCost const &cost, rational scale);

AsyCost operator*(rational scale, AsyCost const &cost);

AsyCost operator/(AsyCost const &cost, rational scale);

bool operator==(AsyCost const &lhs, AsyCost const &rhs);

bool operator!=(AsyCost const &lhs, AsyCost const &rhs);

bool operator<(AsyCost const &lhs, AsyCost const &rhs);

bool operator>(AsyCost const &lhs, AsyCost const &rhs);

}  // namespace sequant

#endif  // SEQUANT_ASY_COST_HPP
