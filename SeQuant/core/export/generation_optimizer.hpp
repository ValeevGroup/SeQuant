#ifndef SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
#define SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <pv/polymorphic_variant.hpp>

#include <algorithm>
#include <cassert>
#include <stack>
#include <variant>

namespace sequant {

namespace {}  // namespace

template <typename Generator, typename Context = Generator::Context>
class GenerationOptimizer : public Generator {
 private:
  using Object = std::variant<Tensor, Variable>;

  enum class OperationType {
    Create,
    Load,
    Zero,
    Unload,
    Destroy,
    Persist,
    Compute,
  };

  enum class MemoryAction {
    Deallocate,
    None,
    Allocate,
  };

  class AbstractOperation {
   public:
    AbstractOperation(Object obj) : m_object(std::move(obj)) {}

    void execute(Generator &generator, const Context &ctx) {
      std::visit([&](const auto &obj) { execute(obj, generator, ctx); },
                 m_object);
    }

    virtual void execute(const Tensor &tensor, Generator &generator,
                         const Context &ctx) = 0;
    virtual void execute(const Variable &variable, Generator &generator,
                         const Context &ctx) = 0;

    virtual OperationType type() const = 0;

    bool cancels(const AbstractOperation &other) const {
      if ((type() == OperationType::Load &&
           other.type() == OperationType::Unload) ||
          (type() == OperationType::Unload &&
           other.type() == OperationType::Load)) {
        return object_equals(other.m_object);
      }

      return false;
    }

    bool pairs_with(const AbstractOperation &operation) const {
      if (!object_equals(operation.m_object)) {
        return false;
      }

      const OperationType other = operation.type();

      switch (type()) {
        case OperationType::Create:
        case OperationType::Load:
          return other == OperationType::Destroy ||
                 other == OperationType::Persist ||
                 other == OperationType::Unload;
        case OperationType::Unload:
        case OperationType::Destroy:
        case OperationType::Persist:
          return other == OperationType::Create || other == OperationType::Load;
        case OperationType::Zero:
        case OperationType::Compute:
          return false;
      }

      // should not be reached
      assert(false);
      return false;
    }

    MemoryAction memory_action() const {
      switch (type()) {
        case OperationType::Create:
        case OperationType::Load:
          return MemoryAction::Allocate;
        case OperationType::Unload:
        case OperationType::Destroy:
        case OperationType::Persist:
          return MemoryAction::Deallocate;
        case OperationType::Zero:
        case OperationType::Compute:
          return MemoryAction::None;
      }

      // should not be reached
      assert(false);
      return MemoryAction::None;
    }

   private:
    Object m_object;

    bool object_equals(const Object &other) const {
      if (m_object.index() != other.index()) {
        return false;
      }

      if (std::holds_alternative<Tensor>(m_object)) {
        const Tensor &lhs = std::get<Tensor>(m_object);
        const Tensor &rhs = std::get<Tensor>(other);

        TensorBlockCompare cmp;
        return !cmp(lhs, rhs) && !cmp(rhs, lhs);
      }

      return std::get<Variable>(m_object) == std::get<Variable>(other);
    }
  };

  class CreateOperation final : public AbstractOperation {
   public:
    CreateOperation(Object obj, bool zero_init)
        : AbstractOperation(std::move(obj)), m_zero_init(zero_init) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::create(tensor, m_zero_init, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::create(variable, m_zero_init, ctx);
    }

    OperationType type() const override { return OperationType::Create; }

   private:
    bool m_zero_init;
  };

  class LoadOperation final : public AbstractOperation {
   public:
    LoadOperation(Object obj, bool set_to_zero)
        : AbstractOperation(std::move(obj)), m_set_to_zero(set_to_zero) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::load(tensor, m_set_to_zero, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::load(variable, m_set_to_zero, ctx);
    }

    OperationType type() const override { return OperationType::Load; }

   private:
    bool m_set_to_zero;
  };

  class ZeroOperation final : public AbstractOperation {
   public:
    ZeroOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::set_to_zero(tensor, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::set_to_zero(variable, ctx);
    }

    OperationType type() const override { return OperationType::Zero; }
  };

  class UnloadOperation final : public AbstractOperation {
   public:
    UnloadOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::unload(tensor, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::unload(variable, ctx);
    }

    OperationType type() const override { return OperationType::Unload; }
  };

  class DestroyOperation final : public AbstractOperation {
   public:
    DestroyOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::destroy(tensor, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::destroy(variable, ctx);
    }

    OperationType type() const override { return OperationType::Destroy; }
  };

  class PersistOperation final : public AbstractOperation {
   public:
    PersistOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::persist(tensor, ctx);
    }

    void execute(const Variable &variable, Generator &generator,
                 const Context &ctx) override {
      generator.Generator::persist(variable, ctx);
    }

    OperationType type() const override { return OperationType::Persist; }
  };

  class ComputeOperation final : public AbstractOperation {
   public:
    ComputeOperation(Object result, ExprPtr expr)
        : AbstractOperation(std::move(result)), m_expr(std::move(expr)) {}

    void execute(const Tensor &result, Generator &generator,
                 const Context &ctx) override {
      assert(m_expr);
      generator.Generator::compute(*m_expr, result, ctx);
    }

    void execute(const Variable &result, Generator &generator,
                 const Context &ctx) override {
      assert(m_expr);
      generator.Generator::compute(*m_expr, result, ctx);
    }

    OperationType type() const override { return OperationType::Compute; }

   private:
    ExprPtr m_expr;
  };

  using Operation =
      pv::polymorphic_variant<AbstractOperation, CreateOperation, LoadOperation,
                              ZeroOperation, UnloadOperation, DestroyOperation,
                              PersistOperation, ComputeOperation>;

 public:
  using Generator::Generator;

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    m_queue.emplace_back(CreateOperation(tensor, zero_init));
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    m_queue.emplace_back(LoadOperation(tensor, set_to_zero));
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    m_queue.emplace_back(ZeroOperation(tensor));
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_queue.emplace_back(UnloadOperation(tensor));
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    m_queue.emplace_back(DestroyOperation(tensor));
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_queue.emplace_back(PersistOperation(tensor));
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    m_queue.emplace_back(CreateOperation(variable, zero_init));
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_queue.emplace_back(LoadOperation(variable, set_to_zero));
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_queue.emplace_back(ZeroOperation(variable));
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_queue.emplace_back(UnloadOperation(variable));
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    m_queue.emplace_back(DestroyOperation(variable));
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_queue.emplace_back(PersistOperation(variable));
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_cache(*this, ctx);
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_cache(*this, ctx);
  }

  void end_expression(const Context &ctx) override {
    // Assumption: Context is the same as for all queued operations
    process_operation_cache(*this, ctx);

    assert(m_queue.empty());
    assert(m_cache.empty());

    Generator::end_expression(ctx);
  }

  void end_export(const Context &ctx) override {
    assert(m_queue.empty());
    assert(m_cache.empty());

    Generator::end_export(ctx);
  }

 private:
  container::svector<Operation> m_queue;
  container::svector<Operation> m_cache;

  template <typename Iterator>
  static Iterator find_unpaired_allocation(Iterator begin, const Iterator end) {
    if (begin == end) {
      return begin;
    }

    std::stack<Operation> deallocations;

    while (begin != end) {
      switch ((*begin)->memory_action()) {
        case MemoryAction::Deallocate:
          deallocations.push(*begin);
          break;
        case MemoryAction::None:
          break;
        case MemoryAction::Allocate: {
          const Operation &op = *begin;
          if (!deallocations.empty() && op->pairs_with(deallocations.top())) {
            deallocations.pop();
          } else {
            return begin;
          }
          break;
        }
      }

      ++begin;
    }

    return end;
  }

  void process_operation_cache(Generator &generator, const Context &ctx) {
    // Remove operations that cancel each other from the queue
    for (std::size_t i = 0; i < m_queue.size(); ++i) {
      const Operation &first = m_queue.at(i);

      for (std::size_t k = i; k < m_queue.size(); ++k) {
        const Operation &second = m_queue.at(k);

        if (first->cancels(second)) {
          m_queue.erase(m_queue.begin() + k);
          m_queue.erase(m_queue.begin() + i);
          // Needed to make the outer loop re-visit the current
          // value of i
          i--;
          break;
        }
      }
    }

    // Move items from the queue into the cache
    for (Operation &op : m_queue) {
      if (op->memory_action() != MemoryAction::Deallocate) {
        m_cache.push_back(std::move(op));
        continue;
      }

      auto op_alloc = std::find_if(
          m_cache.rbegin(), m_cache.rend(),
          [&op](const Operation &other) { return op->pairs_with(other); });

      if (op_alloc == m_cache.rend()) {
        throw std::runtime_error(
            "GenerationOptimizer: Attempt to deallocate object that hasn't "
            "been allocated before!");
      }

      assert((*op_alloc)->memory_action() == MemoryAction::Allocate);

      std::size_t rpos = std::distance(m_cache.rbegin(), op_alloc);
      assert(rpos < m_cache.size());

      // Search for an unpaired allocation in the range (m_cache.end(),
      // op_alloc). That is, between where the allocation for op happens and the
      // end of the cache (starting from the back).
      // In such cases, we need to re-order operations in order to maintain
      // compatibility to stack-based memory models.
      // Example:
      // Load B
      // Load C
      // Drop C
      // Load D
      // Drop B
      // Here, loading of D would be the unpaired allocation we are looking for.
      // This allocation needs to be moved before the allocation in B.
      auto unpaired_alloc =
          find_unpaired_allocation(m_cache.rbegin(), op_alloc);

      while (unpaired_alloc != op_alloc) {
        assert(!(*unpaired_alloc)->pairs_with(op));
        assert(std::find_if(m_cache.begin(), m_cache.end(),
                            [&](const Operation &op) {
                              return op->pairs_with(*unpaired_alloc);
                            }) == m_cache.end());
        assert(rpos > 0);
        assert((*(m_cache.rbegin() + rpos))->pairs_with(op));

        // Move the unpaired allocation to be located before the allocation
        // of op. That is, we do
        // M, A, B, C, X -> X, M, A, B, C
        // where M is the allocation pairing with op and X is the
        // unpaired allocation.
        std::rotate(unpaired_alloc, unpaired_alloc + 1,
                    m_cache.rbegin() + rpos + 1);

        rpos -= 1;
        op_alloc = m_cache.rbegin() + rpos;

        unpaired_alloc = find_unpaired_allocation(m_cache.rbegin(), op_alloc);
      }

      m_cache.push_back(std::move(op));
    }

    m_queue.clear();

    if (!m_cache.empty() && m_cache.front()->pairs_with(m_cache.back())) {
      // We have a full, memory-stack-order-preserving chain of operations
      // -> execute it
      for (Operation &operation : m_cache) {
        operation->execute(generator, ctx);
      }

      m_cache.clear();
    }
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
