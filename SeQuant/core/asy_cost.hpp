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
/// multiplier times a product of space sizes raised to non-negative integer
/// powers. Spaces are identified by `IndexSpace`; `AsyCost` orders and prints
/// them using `IndexSpace`'s own ordering and `base_key()`, and never consults
/// an `IndexSpaceRegistry`.
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
  /// A single term of an AsyCost: a rational prefactor times a
  /// product of index-space sizes raised to per-space exponents, e.g.
  /// `3/2 * O^2 V^4`.
  class AsyCostEntry {
    ExponentMap exponents_;       // space -> power; zero exponents not stored
    mutable rational prefactor_;  // rational multiplier
    bool is_max_ = false;         // true for the AsyCost::max() sentinel

   public:
    /// Write a rational to a stream as `num`, or `num/den` when its
    /// denominator is not 1.
    /// \param os Stream to write to.
    /// \param r Rational to format
    /// \return `os`
    static std::ostream &stream_out_rational(std::ostream &os,
                                             rational const &r);

    /// \return The sentinel entry representing an infinitely scaling cost; it
    ///         compares greater than every non-max entry.
    static AsyCostEntry max();

    /// \return The canonical zero entry (no exponents, zero prefactor).
    static AsyCostEntry const &zero();

    /// Constructors
    AsyCostEntry();

    AsyCostEntry(AsyCostEntry const &) = default;

    AsyCostEntry(AsyCostEntry &&) = default;

    /// \param exponents Map from index space to its exponent. Zero exponents
    ///                  are dropped.
    /// \param prefactor Rational multiplier. If it is zero, or no exponents
    ///                  remain after dropping zeros, the entry collapses to
    ///                  zero.
    AsyCostEntry(ExponentMap exponents, rational prefactor);

    AsyCostEntry &operator=(AsyCostEntry const &) = default;

    AsyCostEntry &operator=(AsyCostEntry &&) = default;

    /// \return The per-space exponents of this term.
    [[nodiscard]] ExponentMap const &exponents() const;

    /// \return The rational multiplier (prefactor) of this term, e.g. `3/2`
    ///         in `3/2 * O^2 V^4`.
    [[nodiscard]] rational prefactor() const;

    /// Set the rational multiplier. `const` because the prefactor is not part
    /// of the term's identity (see \ref operator==).
    void set_prefactor(rational n) const;

    /// \return Whether this is the zero entry.
    [[nodiscard]] bool is_zero() const;

    /// \return Whether this is the \ref max() sentinel.
    [[nodiscard]] bool is_max() const;

    /// Order by the \ref max() sentinel (greatest), then by total polynomial
    /// degree (sum of exponents), then space by space from highest- to
    /// lowest-priority IndexSpace. Independent of the prefactor.
    bool operator<(AsyCostEntry const &rhs) const;

    /// \return Whether two entries share the same monomial (same exponents and
    ///         max-ness). The prefactor is not compared.
    bool operator==(AsyCostEntry const &rhs) const;

    bool operator!=(AsyCostEntry const &rhs) const;

    /// \return A plain-text rendering of this term, e.g. `3/2*O^2V^4`.
    [[nodiscard]] std::string text() const;

    /// \return A LaTeX rendering of this term.
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
  /// \param prefactor Rational multiplier; defaults to 1.
  ///
  explicit AsyCost(ExponentMap exponents, rational prefactor = 1);

  AsyCost(AsyCost const &) = default;

  AsyCost(AsyCost &&) = default;

  AsyCost &operator=(AsyCost const &) = default;

  AsyCost &operator=(AsyCost &&) = default;

  ///
  /// Substitute each space in this cost by an extent and evaluate.
  /// \param extents Map from index space to extent (size). A space appearing
  ///                in this cost but absent from `extents` falls back to its
  ///                `IndexSpace::approximate_size()`. Defaults to empty, in
  ///                which case every space uses its `approximate_size()`.
  /// \return Numerical value of the cost.
  ///
  [[nodiscard]] double ops(ExtentMap const &extents = {}) const;

  /// \return A LaTeX rendering of the whole cost: its terms joined by ` + `,
  ///         most expensive first. The zero cost renders as `0`.
  [[nodiscard]] std::wstring to_latex() const;

  /// \return A plain-text rendering of the whole cost: its terms joined by
  ///         ` + `, most expensive first. The zero cost renders as `0`.
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
