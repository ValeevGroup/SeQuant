#ifndef SEQUANT_ASY_COST_HPP
#define SEQUANT_ASY_COST_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/wstring.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>
#include <sstream>

namespace sequant {

///
/// Represents a symbolic asymptotic cost in terms of active_occupied
/// and the rest orbitals.
/// eg.
///     - AsyCost{2,4} implies scaling of $O^2V^4$. In other words, the cost
///       scales by the second power in the number of active_occupied orbitals
///       and the fourth power in the number of the rest orbitals.
///     - AsyCost{sequant::rational{1,2}, 2, 4} implies the same scaling as
///       above except the numeric value obtained by substituting $O$ and $V$
///       numbers is then halved.
///
class AsyCost {
 private:
  class AsyCostEntry {
    size_t occ_;              // power of active_occupied
    size_t virt_;             // power of the rest orbitals
    mutable rational count_;  // count of this asymptotic symbol

   public:
    static std::ostream &stream_out_rational(std::ostream &os,
                                             rational const &r);

    static AsyCostEntry max();

    static AsyCostEntry const &zero();

    AsyCostEntry(size_t nocc, size_t nvirt, rational count);

    AsyCostEntry(AsyCostEntry const &) = default;

    AsyCostEntry(AsyCostEntry &&) = default;

    AsyCostEntry &operator=(AsyCostEntry const &) = default;

    AsyCostEntry &operator=(AsyCostEntry &&) = default;

    size_t occ() const;

    size_t virt() const;

    rational count() const;

    void set_count(rational n) const;

    bool operator<(AsyCostEntry const &rhs) const;

    bool operator==(AsyCostEntry const &rhs) const;

    bool operator!=(AsyCostEntry const &rhs) const;

    std::string text() const;

    std::string to_latex() const;
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
  /// \param count Rational number of times this cost repeats.
  /// \param nocc Asymptotic scaling exponent in the active occupied orbitals.
  /// \param nvirt Asymptotic scaling exponent in the active unoccupied
  ///              orbitals.
  ///
  AsyCost(rational count, size_t nocc, size_t nvirt);

  ///
  /// \param nocc Asymptotic scaling exponent in the active occupied orbitals.
  /// \param nvirt Asymptotic scaling exponent in the active unoccupied
  ///              orbitals.
  ///
  AsyCost(size_t nocc, size_t nvirt);

  ///
  ///
  /// \param ov A pair of size_ts.
  ///           ov.first is the asymptotic scaling exponent in the active
  ///           occupied orbitals.
  ///           ov.second is that in the active unoccupied orbitals
  ///
  explicit AsyCost(std::pair<size_t, size_t> const &ov);

  AsyCost(AsyCost const &) = default;

  AsyCost(AsyCost &&) = default;

  AsyCost &operator=(AsyCost const &) = default;

  AsyCost &operator=(AsyCost &&) = default;

  ///
  /// \param nocc Substitute $O$ by nocc.
  /// \param nvirt Substitute $V$ by nvirt.
  /// \return Scaled asymptotic cost.
  [[nodiscard]] double ops(size_t nocc, size_t nvirt) const;

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
