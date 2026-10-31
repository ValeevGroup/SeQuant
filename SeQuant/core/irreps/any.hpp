#ifndef SEQUANT_CORE_IRREPS_ANY_HPP
#define SEQUANT_CORE_IRREPS_ANY_HPP

#include <SeQuant/core/irreps/algebra.hpp>
#include <SeQuant/core/irreps/propagation.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace sequant {

class Tensor;
class Variable;

/// Type-erased IRREP (of arbitrary IRREP type). A moved-from object may only
/// be assigned to or destroyed.
class AnyIrrep {
 public:
  template <typename Irrep>
    requires(!std::same_as<std::remove_cvref_t<Irrep>, AnyIrrep>)
  explicit AnyIrrep(Irrep irrep)
      : m_impl(std::make_unique<Model<std::remove_cvref_t<Irrep>>>(
            std::move(irrep))) {}

  AnyIrrep(const AnyIrrep &other);
  AnyIrrep(AnyIrrep &&other) noexcept;
  AnyIrrep &operator=(const AnyIrrep &other);
  AnyIrrep &operator=(AnyIrrep &&other) noexcept;
  ~AnyIrrep();

  /// @returns The string representation of this IRREP as produced by
  /// irrep_to_string
  [[nodiscard]] std::string to_string() const;

  /// @returns A pointer to the wrapped IRREP, if it is of type Irrep.
  /// Otherwise, nullptr.
  template <typename Irrep>
  [[nodiscard]] const Irrep *as() const {
    const auto *model = dynamic_cast<const Model<Irrep> *>(m_impl.get());
    return model ? &model->value : nullptr;
  }

  friend bool operator==(const AnyIrrep &lhs, const AnyIrrep &rhs);

 private:
  struct Concept {
    virtual ~Concept() = default;
    virtual std::unique_ptr<Concept> clone() const = 0;
    virtual std::string to_string() const = 0;
    virtual bool equals(const Concept &other) const = 0;
  };

  template <typename Irrep>
  struct Model final : Concept {
    explicit Model(Irrep value) : value(std::move(value)) {}

    std::unique_ptr<Concept> clone() const override {
      return std::make_unique<Model>(value);
    }
    std::string to_string() const override { return irrep_to_string(value); }
    bool equals(const Concept &other) const override {
      const auto *model = dynamic_cast<const Model *>(&other);
      return model && model->value == value;
    }

    Irrep value;
  };

  std::unique_ptr<Concept> m_impl;
};

/// Type-erased IrrepSet (of arbitrary IRREP type). A default-constructed
/// object represents an empty set whose IRREP type is not yet known. It
/// adopts the IRREP type of the first set it is united with. Empty sets
/// compare equal regardless of their IRREP type.
class AnyIrrepSet {
 public:
  /// Creates an empty set
  AnyIrrepSet();
  template <typename Irrep>
  AnyIrrepSet(IrrepSet<Irrep> set)
      : m_impl(std::make_unique<Model<Irrep>>(std::move(set))) {}

  AnyIrrepSet(const AnyIrrepSet &other);
  AnyIrrepSet(AnyIrrepSet &&other) noexcept;
  AnyIrrepSet &operator=(const AnyIrrepSet &other);
  AnyIrrepSet &operator=(AnyIrrepSet &&other) noexcept;
  ~AnyIrrepSet();

  /// @sa IrrepSet::is_any
  [[nodiscard]] bool is_any() const;
  /// @sa IrrepSet::empty
  [[nodiscard]] bool empty() const;
  /// @sa IrrepSet::size
  [[nodiscard]] std::size_t size() const;

  /// Turns this set into the union of itself and the given set
  /// @throws Exception if the two sets have different IRREP types
  AnyIrrepSet &unite(const AnyIrrepSet &other);

  /// @returns The string representation of this set
  /// @sa sequant::to_string(const IrrepSet<Irrep>&)
  [[nodiscard]] std::string to_string() const;

  /// Calls f(const AnyIrrep&) for every IRREP explicitly contained in this set
  /// (in ascending order)
  void for_each(const std::function<void(const AnyIrrep &)> &f) const;

  /// @returns A pointer to the wrapped set, if it is of type IrrepSet<Irrep>.
  /// Otherwise, nullptr.
  template <typename Irrep>
  [[nodiscard]] const IrrepSet<Irrep> *as() const {
    const auto *model = dynamic_cast<const Model<Irrep> *>(m_impl.get());
    return model ? &model->set : nullptr;
  }

  friend bool operator==(const AnyIrrepSet &lhs, const AnyIrrepSet &rhs);

 private:
  struct Concept {
    virtual ~Concept() = default;
    virtual std::unique_ptr<Concept> clone() const = 0;
    virtual bool is_any() const = 0;
    virtual bool empty() const = 0;
    virtual std::size_t size() const = 0;
    virtual std::string to_string() const = 0;
    virtual bool equals(const Concept &other) const = 0;
    /// @returns Whether the union could be formed (false on type mismatch)
    virtual bool unite(const Concept &other) = 0;
    virtual void for_each(
        const std::function<void(const AnyIrrep &)> &f) const = 0;
  };

  template <typename Irrep>
  struct Model final : Concept {
    explicit Model(IrrepSet<Irrep> set) : set(std::move(set)) {}

    std::unique_ptr<Concept> clone() const override {
      return std::make_unique<Model>(set);
    }
    bool is_any() const override { return set.is_any(); }
    bool empty() const override { return set.empty(); }
    std::size_t size() const override { return set.size(); }
    std::string to_string() const override { return sequant::to_string(set); }
    bool equals(const Concept &other) const override {
      const auto *model = dynamic_cast<const Model *>(&other);
      return model && model->set == set;
    }
    bool unite(const Concept &other) override {
      const auto *model = dynamic_cast<const Model *>(&other);
      if (!model) {
        return false;
      }
      set.unite(model->set);
      return true;
    }
    void for_each(
        const std::function<void(const AnyIrrep &)> &f) const override {
      for (const Irrep &current : set) {
        f(AnyIrrep(current));
      }
    }

    IrrepSet<Irrep> set;
  };

  std::unique_ptr<Concept> m_impl;
};

/// Type-erased IrrepContext (for arbitrary IRREP algebras). It satisfies
/// IrrepPropagationContext and can thus be used with irreps_of. A moved-from
/// object may only be assigned to or destroyed.
class AnyIrrepContext {
 public:
  using set_type = AnyIrrepSet;

  template <IrrepAlgebra A>
  AnyIrrepContext(IrrepContext<A> ctx)
      : m_impl(std::make_unique<Model<A>>(std::move(ctx))) {}

  AnyIrrepContext(const AnyIrrepContext &other);
  AnyIrrepContext(AnyIrrepContext &&other) noexcept;
  AnyIrrepContext &operator=(const AnyIrrepContext &other);
  AnyIrrepContext &operator=(AnyIrrepContext &&other) noexcept;
  ~AnyIrrepContext();

  /// @sa IrrepContext::irreps_of
  [[nodiscard]] AnyIrrepSet irreps_of(const Tensor &tensor) const;
  [[nodiscard]] AnyIrrepSet irreps_of(const Variable &variable) const;
  [[nodiscard]] AnyIrrepSet totally_symmetric() const;
  /// @returns Whether the given set contains nothing but the totally symmetric
  /// IRREP
  /// @throws Exception if the set is not of this context's IRREP type
  [[nodiscard]] bool is_totally_symmetric(const AnyIrrepSet &set) const;
  [[nodiscard]] AnyIrrepSet empty_set() const;
  /// @sa IrrepContext::canonical_irrep
  [[nodiscard]] AnyIrrep canonical_irrep(std::string_view name) const;
  /// @throws Exception if the sets are not of this context's IRREP type
  [[nodiscard]] AnyIrrepSet product(const AnyIrrepSet &lhs,
                                    const AnyIrrepSet &rhs) const;
  /// @throws Exception if the set is not of this context's IRREP type
  [[nodiscard]] AnyIrrepSet conjugate(const AnyIrrepSet &set) const;

  /// @returns A pointer to the wrapped context, if it is of type
  /// IrrepContext<A>. Otherwise, nullptr.
  template <IrrepAlgebra A>
  [[nodiscard]] const IrrepContext<A> *as() const {
    const auto *model = dynamic_cast<const Model<A> *>(m_impl.get());
    return model ? &model->ctx : nullptr;
  }

 private:
  struct Concept {
    virtual ~Concept() = default;
    virtual std::unique_ptr<Concept> clone() const = 0;
    virtual AnyIrrepSet irreps_of(const Tensor &tensor) const = 0;
    virtual AnyIrrepSet irreps_of(const Variable &variable) const = 0;
    virtual AnyIrrepSet totally_symmetric() const = 0;
    virtual bool is_totally_symmetric(const AnyIrrepSet &set) const = 0;
    virtual AnyIrrepSet empty_set() const = 0;
    virtual AnyIrrep canonical_irrep(std::string_view name) const = 0;
    virtual AnyIrrepSet product(const AnyIrrepSet &lhs,
                                const AnyIrrepSet &rhs) const = 0;
    virtual AnyIrrepSet conjugate(const AnyIrrepSet &set) const = 0;
  };

  template <IrrepAlgebra A>
  struct Model final : Concept {
    using Set = typename IrrepContext<A>::set_type;
    using Irrep = typename IrrepContext<A>::irrep_type;

    explicit Model(IrrepContext<A> ctx) : ctx(std::move(ctx)) {}

    std::unique_ptr<Concept> clone() const override {
      return std::make_unique<Model>(ctx);
    }
    AnyIrrepSet irreps_of(const Tensor &tensor) const override {
      return ctx.irreps_of(tensor);
    }
    AnyIrrepSet irreps_of(const Variable &variable) const override {
      return ctx.irreps_of(variable);
    }
    AnyIrrepSet totally_symmetric() const override {
      return ctx.totally_symmetric();
    }
    bool is_totally_symmetric(const AnyIrrepSet &set) const override {
      return ctx.is_totally_symmetric(typed(set));
    }
    AnyIrrepSet empty_set() const override { return ctx.empty_set(); }
    AnyIrrep canonical_irrep(std::string_view name) const override {
      return AnyIrrep(ctx.canonical_irrep(name));
    }
    AnyIrrepSet product(const AnyIrrepSet &lhs,
                        const AnyIrrepSet &rhs) const override {
      return ctx.product(typed(lhs), typed(rhs));
    }
    AnyIrrepSet conjugate(const AnyIrrepSet &set) const override {
      return ctx.conjugate(typed(set));
    }

    /// @returns The typed version of the given set
    static Set typed(const AnyIrrepSet &set) {
      if (const Set *typed_set = set.as<Irrep>()) {
        return *typed_set;
      }
      if (set.empty()) {
        // Empty set of yet-unknown type
        return {};
      }
      throw Exception("IRREP set " + set.to_string() +
                      " is of a type incompatible with the IRREP context");
    }

    IrrepContext<A> ctx;
  };

  std::unique_ptr<Concept> m_impl;
};

static_assert(IrrepPropagationContext<AnyIrrepContext>);

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_ANY_HPP
