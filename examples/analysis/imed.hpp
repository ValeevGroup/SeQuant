
#ifndef SEQUANT_ANALYSIS_IMED_HPP
#define SEQUANT_ANALYSIS_IMED_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/expr.hpp>

namespace sequant {
///
/// Represents an intermediate of sequant expression evaluation.
///
struct Imed {
  ///
  /// \brief The identifier.
  ///
  /// If two objects have same \c key their ExprPtr
  /// representation is equivalent modulo the index labels.
  ///
  size_t key;

  ///
  /// \brief The flops required to compute the intermediate.
  ///
  AsyCost flops;

  ///
  /// \brief The memory required to hold the intermediate.
  ///
  AsyCost memory;

  ///
  /// \brief The position of the terms in a an many-body equation (Sum)
  ///        that contain this intermediate.
  ///
  container::vector<size_t> pos;

  ///
  /// \brief The ExprPtr of the intermediate as seen in different terms.
  ///        @see pos
  ///
  Sum expr;

  explicit Imed(size_t);

  struct CompKey {
    using is_transparent = void;

    [[nodiscard]] bool operator()(Imed const&, Imed const&) const noexcept;

    [[nodiscard]] bool operator()(size_t, Imed const&) const noexcept;

    [[nodiscard]] bool operator()(Imed const&, size_t) const noexcept;
  };

  template <typename Comp = std::less<AsyCost>>
  struct CompFlops {
    using is_transparent = void;

    [[nodiscard]] bool operator()(Imed const& left,
                                  Imed const& right) const noexcept {
      return Comp{}(left.flops, right.flops);
    }

    [[nodiscard]] bool operator()(AsyCost const& cost,
                                  Imed const& right) const noexcept {
      return Comp{}(cost, right.flops);
    }

    [[nodiscard]] bool operator()(Imed const& left,
                                  AsyCost const& cost) const noexcept {
      return Comp{}(left.flops, cost);
    }
  };

  template <typename Comp = std::less<AsyCost>>
  struct CompMemory {
    using is_transparent = void;

    [[nodiscard]] bool operator()(Imed const& left,
                                  Imed const& right) const noexcept {
      return Comp{}(left.memory, right.memory);
    }

    [[nodiscard]] bool operator()(AsyCost const& cost,
                                  Imed const& right) const noexcept {
      return Comp{}(cost, right.memory);
    }

    [[nodiscard]] bool operator()(Imed const& left,
                                  AsyCost const& cost) const noexcept {
      return Comp{}(left.memory, cost);
    }
  };

  template <typename Comp = std::less<size_t>>
  struct CompCount {
    using is_transparent = void;

    [[nodiscard]] bool operator()(Imed const& left,
                                  Imed const& right) const noexcept {
      return Comp{}(left.pos.size(), right.pos.size());
    }

    [[nodiscard]] bool operator()(size_t count,
                                  Imed const& right) const noexcept {
      return Comp{}(count, right.pos.size());
    }

    [[nodiscard]] bool operator()(Imed const& left,
                                  size_t count) const noexcept {
      return Comp{}(left.pos.size(), count);
    }
  };

  ///
  /// \return Returns a predicate that takes an Imed argument and calls Comp
  ///         on the argument's no. of repetitions and \c count:
  ///         Comp{}(Imed.pos.size(), count)
  ///
  template <typename Comp = std::less<size_t>>
  [[nodiscard]] static auto filter_count(size_t count,
                                         Comp const& fun = {}) noexcept {
    return [count, &fun](Imed const& i) -> bool {
      return fun(i.pos.size(), count);
    };
  }

  template <typename Comp = std::less<AsyCost>>
  [[nodiscard]] static auto filter_flops(AsyCost const& cost,
                                         Comp const& fun = {}) noexcept {
    return [&cost, &fun](Imed const& i) -> bool { return fun(i.flops, cost); };
  }

  template <typename Comp = std::less<AsyCost>>
  [[nodiscard]] static auto filter_memory(AsyCost const& cost,
                                          Comp const& fun = {}) noexcept {
    return [&cost, &fun](Imed const& i) -> bool { return fun(i.memory, cost); };
  }
};

[[nodiscard]] bool operator<(Imed const&, Imed const&) noexcept;
}  // namespace sequant

#endif  // SEQUANT_ANALYSIS_IMED_HPP
