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

template <typename MainGenerator, typename MainContext = MainGenerator::Context>
class GenerationOptimizer final : public Generator<MainContext> {
 private:
  using Object = std::variant<Tensor, Variable>;

  enum class OperationType {
    Create,
    Load,
    LoadAndZero,
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

    void execute(MainGenerator &generator, const MainContext &ctx) {
      std::visit([&](const auto &obj) { execute(obj, generator, ctx); },
                 m_object);
    }

    virtual void execute(const Tensor &tensor, MainGenerator &generator,
                         const MainContext &ctx) = 0;
    virtual void execute(const Variable &variable, MainGenerator &generator,
                         const MainContext &ctx) = 0;

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
        case OperationType::LoadAndZero:
          return other == OperationType::Destroy ||
                 other == OperationType::Persist ||
                 other == OperationType::Unload;
        case OperationType::Unload:
        case OperationType::Destroy:
        case OperationType::Persist:
          return other == OperationType::Create ||
                 other == OperationType::Load ||
                 other == OperationType::LoadAndZero;
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
        case OperationType::LoadAndZero:
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

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.create(tensor, m_zero_init, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.create(variable, m_zero_init, ctx);
    }

    OperationType type() const override { return OperationType::Create; }

   private:
    bool m_zero_init;
  };

  class LoadOperation final : public AbstractOperation {
   public:
    LoadOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.load(tensor, false, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.load(variable, false, ctx);
    }

    OperationType type() const override { return OperationType::Load; }
  };

  class LoadAndZeroOperation final : public AbstractOperation {
   public:
    LoadAndZeroOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.load(tensor, true, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.load(variable, true, ctx);
    }

    OperationType type() const override { return OperationType::LoadAndZero; }
  };

  class ZeroOperation final : public AbstractOperation {
   public:
    ZeroOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.set_to_zero(tensor, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.set_to_zero(variable, ctx);
    }

    OperationType type() const override { return OperationType::Zero; }
  };

  class UnloadOperation final : public AbstractOperation {
   public:
    UnloadOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.unload(tensor, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.unload(variable, ctx);
    }

    OperationType type() const override { return OperationType::Unload; }
  };

  class DestroyOperation final : public AbstractOperation {
   public:
    DestroyOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.destroy(tensor, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.destroy(variable, ctx);
    }

    OperationType type() const override { return OperationType::Destroy; }
  };

  class PersistOperation final : public AbstractOperation {
   public:
    PersistOperation(Object obj) : AbstractOperation(std::move(obj)) {}

    void execute(const Tensor &tensor, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.persist(tensor, ctx);
    }

    void execute(const Variable &variable, MainGenerator &generator,
                 const MainContext &ctx) override {
      generator.persist(variable, ctx);
    }

    OperationType type() const override { return OperationType::Persist; }
  };

  class ComputeOperation final : public AbstractOperation {
   public:
    ComputeOperation(Object result, ExprPtr expr)
        : AbstractOperation(std::move(result)), m_expr(std::move(expr)) {}

    void execute(const Tensor &result, MainGenerator &generator,
                 const MainContext &ctx) override {
      assert(m_expr);
      generator.compute(*m_expr, result, ctx);
    }

    void execute(const Variable &result, MainGenerator &generator,
                 const MainContext &ctx) override {
      assert(m_expr);
      generator.compute(*m_expr, result, ctx);
    }

    OperationType type() const override { return OperationType::Compute; }

   private:
    ExprPtr m_expr;
  };

  using Operation =
      pv::polymorphic_variant<AbstractOperation, CreateOperation, LoadOperation,
                              LoadAndZeroOperation, ZeroOperation,
                              UnloadOperation, DestroyOperation,
                              PersistOperation, ComputeOperation>;

 public:
  GenerationOptimizer(MainGenerator generator)
      : m_generator(std::move(generator)) {}

  std::string get_format_name() const override {
    return m_generator.get_format_name() + " (optimized)";
  }

  void create(const Tensor &tensor, bool zero_init,
              const MainContext &ctx) override {
    m_queue.emplace_back(CreateOperation(tensor, zero_init));
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const MainContext &ctx) override {
    if (set_to_zero) {
      m_queue.emplace_back(LoadAndZeroOperation(tensor));
    } else {
      m_queue.emplace_back(LoadOperation(tensor));
    }
  }

  void set_to_zero(const Tensor &tensor, const MainContext &ctx) override {
    m_queue.emplace_back(ZeroOperation(tensor));
  }

  void unload(const Tensor &tensor, const MainContext &ctx) override {
    m_queue.emplace_back(UnloadOperation(tensor));
  }

  void destroy(const Tensor &tensor, const MainContext &ctx) override {
    m_queue.emplace_back(DestroyOperation(tensor));
  }

  void persist(const Tensor &tensor, const MainContext &ctx) override {
    m_queue.emplace_back(PersistOperation(tensor));
  }

  void create(const Variable &variable, bool zero_init,
              const MainContext &ctx) override {
    m_queue.emplace_back(CreateOperation(variable, zero_init));
  }

  void load(const Variable &variable, bool set_to_zero,
            const MainContext &ctx) override {
    if (set_to_zero) {
      m_queue.emplace_back(LoadAndZeroOperation(variable));
    } else {
      m_queue.emplace_back(LoadOperation(variable));
    }
  }

  void set_to_zero(const Variable &variable, const MainContext &ctx) override {
    m_queue.emplace_back(ZeroOperation(variable));
  }

  void unload(const Variable &variable, const MainContext &ctx) override {
    m_queue.emplace_back(UnloadOperation(variable));
  }

  void destroy(const Variable &variable, const MainContext &ctx) override {
    m_queue.emplace_back(DestroyOperation(variable));
  }

  void persist(const Variable &variable, const MainContext &ctx) override {
    m_queue.emplace_back(PersistOperation(variable));
  }

  void compute(const Expr &expression, const Tensor &result,
               const MainContext &ctx) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_cache(ctx);
  }

  void compute(const Expr &expression, const Variable &result,
               const MainContext &ctx) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_cache(ctx);
  }

  void end_expression(const MainContext &ctx) override {
    // Assumption: Context is the same as for all queued operations
    process_operation_cache(ctx);

    assert(m_queue.empty());
    assert(m_cache.empty());

    m_generator.end_expression(ctx);
  }

  void end_export(const MainContext &ctx) override {
    assert(m_queue.empty());
    assert(m_cache.empty());

    m_generator.end_export(ctx);
  }

  /////////////////////////////////////////////////////////
  /////////// Pass-through implementations ////////////////
  /////////////////////////////////////////////////////////

  // clang-format off
  bool supports_named_sections() const override { return m_generator.supports_named_sections(); }
  bool requires_named_sections() const override { return m_generator.requires_named_sections(); }
  DeclarationScope index_declaration_scope() const override { return m_generator.index_declaration_scope(); }
  DeclarationScope variable_declaration_scope() const override { return m_generator.variable_declaration_scope(); }
  DeclarationScope tensor_declaration_scope() const override { return m_generator.tensor_declaration_scope(); }
  std::string represent(const Index &idx, const MainContext &ctx) const override { return m_generator.represent(idx, ctx); }
  std::string represent(const Tensor &tensor, const MainContext &ctx) const override { return m_generator.represent(tensor, ctx); }
  std::string represent(const Variable &variable, const MainContext &ctx) const override { return m_generator.represent(variable, ctx); }
  std::string represent(const Constant &constant, const MainContext &ctx) const override { return m_generator.represent(constant, ctx); }
  void declare(const Index &idx, const MainContext &ctx)  override { m_generator.declare(idx, ctx); }
  void declare(const Variable &variable, UsageSet usage, const MainContext &ctx)  override { m_generator.declare(variable, usage, ctx); }
  void declare(const Tensor &tensor, UsageSet usage, const MainContext &ctx)  override { m_generator.declare(tensor, usage, ctx); }
  void all_indices_declared(std::size_t amount, const MainContext &ctx)  override { m_generator.all_indices_declared(amount, ctx); }
  void all_variables_declared(std::size_t amount, const MainContext &ctx)  override { m_generator.all_variables_declared(amount, ctx); }
  void all_tensors_declared(std::size_t amount, const MainContext &ctx)  override { m_generator.all_tensors_declared(amount, ctx); }
  void begin_declarations(DeclarationScope scope, const MainContext &ctx) override { m_generator.begin_declarations(scope, ctx); }
  void end_declarations(DeclarationScope scope, const MainContext &ctx) override { m_generator.end_declarations(scope, ctx); }
  void insert_comment(const std::string &comment, const MainContext &ctx) override { m_generator.insert_comment(comment, ctx); }
  void begin_named_section(std::string_view name, const MainContext &ctx) override { m_generator.begin_named_section(name, ctx); }
  void end_named_section(std::string_view name, const MainContext &ctx) override { m_generator.end_named_section(name, ctx); }
  void begin_expression(const MainContext &ctx) override { m_generator.begin_expression(ctx); }
  void begin_export(const MainContext &ctx) override { m_generator.begin_export(ctx); }
  std::string get_generated_code() const override { return m_generator.get_generated_code(); }
  // clang-format on

 private:
  MainGenerator m_generator;
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

  std::size_t distance_to_matching_alloc(const Operation &dealloc) const {
    assert(dealloc->memory_action() == MemoryAction::Deallocate);

    for (auto it = find_unpaired_allocation(m_cache.rbegin(), m_cache.rend());
         it != m_cache.rend();
         it = find_unpaired_allocation(it + 1, m_cache.rend())) {
      const Operation &alloc = *it;

      if (alloc->pairs_with(dealloc)) {
        assert(std::distance(m_cache.rbegin(), it) >= 0);
        return std::distance(m_cache.rbegin(), it);
      }
    }

    assert(false);
    throw std::runtime_error("Unmatched deallocation!");
  }

  void process_operation_cache(const MainContext &ctx) {
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

    std::sort(m_queue.begin(), m_queue.end(),
              [&](const Operation &lhs, const Operation &rhs) {
                MemoryAction lhs_action = lhs->memory_action();
                MemoryAction rhs_action = rhs->memory_action();

                if (lhs_action != rhs_action) {
                  if (lhs_action == MemoryAction::None ||
                      rhs_action == MemoryAction::None) {
                    // Don't mess with the order of memory actions and
                    // non-memory actions (e.g. computations)
                    return false;
                  }

                  // Make sure deallocations are listed before allocations
                  return lhs_action == MemoryAction::Deallocate;
                }

                if (lhs_action != MemoryAction::Deallocate) {
                  // We only want to reorder deallocations
                  return false;
                }

                // The deallocation whose matching allocation is closest, is
                // prioritized
                return distance_to_matching_alloc(lhs) <
                       distance_to_matching_alloc(rhs);
              });

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
        operation->execute(m_generator, ctx);
      }

      m_cache.clear();
    }
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
