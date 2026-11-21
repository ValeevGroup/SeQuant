#ifndef SEQUANT_CORE_EXPORT_STACKCHECKER_HPP
#define SEQUANT_CORE_EXPORT_STACKCHECKER_HPP

#include <SeQuant/core/export/forwarding_generator.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <span>
#include <string_view>
#include <variant>
#include <vector>

namespace sequant {

/// Simulates the memory stack of a backend that only supports stack-like
/// allocations: an object can only be released if it was the most recently
/// allocated object that is still allocated.
class MemoryStack {
 public:
  /// Sets the indices that distinguish otherwise equal tensors (batching)
  void set_batch_indices(std::span<const Index> indices);

  void allocate(const Tensor &tensor);
  void allocate(const Variable &variable);

  /// @param action Name of the action causing the release (for diagnostics)
  /// @throw Exception if the tensor is not the top of the stack
  void release(const Tensor &tensor, std::string_view action);
  /// @param action Name of the action causing the release (for diagnostics)
  /// @throw Exception if the variable is not the top of the stack
  void release(const Variable &variable, std::string_view action);

  /// @throw Exception if there are objects that have not been released
  void require_empty() const;

 private:
  using Object = std::variant<Tensor, Variable>;

  std::vector<Object> m_stack;
  IndexSpecificTensorBlockEqualComparator m_cmp;

  bool equals(const Object &lhs, const Object &rhs) const;
  void release(const Object &object, std::string_view action);
};

/// Wrapper for other Generators that verifies that the received sequence of
/// memory actions (create, load, unload, destroy, persist) obeys a stack-like
/// order: objects must be released in the reverse order of their allocation.
/// All callbacks are forwarded to the wrapped Generator unchanged.
///
/// This is required by backends that only support stack-like memory
/// allocations (e.g. ITF). Violations are reported by throwing an Exception.
template <typename MainGenerator, typename MainContext = MainGenerator::Context>
class StackCheckingGenerator final
    : public ForwardingGenerator<MainGenerator, MainContext> {
 private:
  using Base = ForwardingGenerator<MainGenerator, MainContext>;

 public:
  StackCheckingGenerator(MainGenerator generator)
      : Base(std::move(generator)) {}

  void create(const Tensor &tensor, bool zero_init,
              const MainContext &ctx) override {
    m_stack.allocate(tensor);
    Base::create(tensor, zero_init, ctx);
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const MainContext &ctx) override {
    m_stack.allocate(tensor);
    Base::load(tensor, set_to_zero, ctx);
  }

  void unload(const Tensor &tensor, const MainContext &ctx) override {
    m_stack.release(tensor, "unload");
    Base::unload(tensor, ctx);
  }

  void destroy(const Tensor &tensor, const MainContext &ctx) override {
    m_stack.release(tensor, "destroy");
    Base::destroy(tensor, ctx);
  }

  void persist(const Tensor &tensor, const MainContext &ctx) override {
    m_stack.release(tensor, "persist");
    Base::persist(tensor, ctx);
  }

  void create(const Variable &variable, bool zero_init,
              const MainContext &ctx) override {
    m_stack.allocate(variable);
    Base::create(variable, zero_init, ctx);
  }

  void load(const Variable &variable, bool set_to_zero,
            const MainContext &ctx) override {
    m_stack.allocate(variable);
    Base::load(variable, set_to_zero, ctx);
  }

  void unload(const Variable &variable, const MainContext &ctx) override {
    m_stack.release(variable, "unload");
    Base::unload(variable, ctx);
  }

  void destroy(const Variable &variable, const MainContext &ctx) override {
    m_stack.release(variable, "destroy");
    Base::destroy(variable, ctx);
  }

  void persist(const Variable &variable, const MainContext &ctx) override {
    m_stack.release(variable, "persist");
    Base::persist(variable, ctx);
  }

  void begin_expression(const MainContext &ctx) override {
    if (this->supports_index_batching()) {
      m_stack.set_batch_indices(ctx.batch_indices(ctx.current_expression_id()));
    }
    Base::begin_expression(ctx);
  }

  void end_expression(const MainContext &ctx) override {
    m_stack.require_empty();
    Base::end_expression(ctx);
  }

 private:
  MemoryStack m_stack;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_STACKCHECKER_HPP
