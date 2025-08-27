#ifndef SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
#define SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <pv/polymorphic_variant.hpp>

#include <algorithm>
#include <cassert>
#include <stack>
#include <utility>
#include <variant>

namespace sequant {

/// Wrapper for other Generators that will perform optimizations on received
/// actions before passing them along to the wrapped generator for actual code
/// generation.
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

    /// @returns Whether this operation cancels the provided operation. That is,
    /// executing them both in sequence will result in a no-op
    bool cancels(const AbstractOperation &other) const {
      if ((type() == OperationType::Load &&
           other.type() == OperationType::Unload) ||
          (type() == OperationType::Unload &&
           other.type() == OperationType::Load)) {
        return object_equals(other.m_object);
      }

      return false;
    }

    /// @returns Whether this operation forms a semantic pair with the provided
    /// one
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

      SEQUANT_UNREACHABLE;
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

      SEQUANT_UNREACHABLE;
    }

    std::string to_string() const {
      std::string str;
      switch (type()) {
        case OperationType::Create:
          str += "Create";
          break;
        case OperationType::Load:
          str += "Load";
          break;
        case OperationType::LoadAndZero:
          str += "Load&Zero";
          break;
        case OperationType::Zero:
          str += "Zero";
          break;
        case OperationType::Unload:
          str += "Unload";
          break;
        case OperationType::Destroy:
          str += "Destroy";
          break;
        case OperationType::Persist:
          str += "Persist";
          break;
        case OperationType::Compute:
          str += "Compute";
          break;
      }

      assert(!str.empty());

      str += " ";

      if (std::holds_alternative<Tensor>(m_object)) {
        const Tensor &tensor = std::get<Tensor>(m_object);
        str += toUtf8(tensor.label());
        str += ":";
        for (const Index &idx : tensor.indices()) {
          str += toUtf8(idx.space().reduce_key(idx.label()));
        }
      } else {
        str += toUtf8(std::get<Variable>(m_object).label());
      }

      return str;
    }

    bool operator==(const AbstractOperation &other) const {
      return type() == other.type() && object_equals(other.m_object);
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

        TensorBlockEqualComparator cmp;
        return cmp(lhs, rhs);
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
              const MainContext &) override {
    m_queue.emplace_back(CreateOperation(tensor, zero_init));
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const MainContext &) override {
    if (set_to_zero) {
      m_queue.emplace_back(LoadAndZeroOperation(tensor));
    } else {
      m_queue.emplace_back(LoadOperation(tensor));
    }
  }

  void set_to_zero(const Tensor &tensor, const MainContext &) override {
    m_queue.emplace_back(ZeroOperation(tensor));
  }

  void unload(const Tensor &tensor, const MainContext &) override {
    m_queue.emplace_back(UnloadOperation(tensor));
  }

  void destroy(const Tensor &tensor, const MainContext &) override {
    m_queue.emplace_back(DestroyOperation(tensor));
  }

  void persist(const Tensor &tensor, const MainContext &) override {
    m_queue.emplace_back(PersistOperation(tensor));
  }

  void create(const Variable &variable, bool zero_init,
              const MainContext &) override {
    m_queue.emplace_back(CreateOperation(variable, zero_init));
  }

  void load(const Variable &variable, bool set_to_zero,
            const MainContext &) override {
    if (set_to_zero) {
      m_queue.emplace_back(LoadAndZeroOperation(variable));
    } else {
      m_queue.emplace_back(LoadOperation(variable));
    }
  }

  void set_to_zero(const Variable &variable, const MainContext &) override {
    m_queue.emplace_back(ZeroOperation(variable));
  }

  void unload(const Variable &variable, const MainContext &) override {
    m_queue.emplace_back(UnloadOperation(variable));
  }

  void destroy(const Variable &variable, const MainContext &) override {
    m_queue.emplace_back(DestroyOperation(variable));
  }

  void persist(const Variable &variable, const MainContext &) override {
    m_queue.emplace_back(PersistOperation(variable));
  }

  void compute(const Expr &expression, const Tensor &result,
               const MainContext &) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_queue();
  }

  void compute(const Expr &expression, const Variable &result,
               const MainContext &) override {
    m_queue.emplace_back(ComputeOperation(result, expression.clone()));
    process_operation_queue();
  }

  void end_expression(const MainContext &ctx) override {
    // Assumption: Context is the same as for all queued operations
    process_operation_queue();
    process_operation_cache(ctx);

    assert(m_queue.empty());
    assert(m_cache.empty());
    assert(m_paired.empty());

    m_generator.end_expression(ctx);
  }

  void end_export(const MainContext &ctx) override {
    assert(m_queue.empty());
    assert(m_cache.empty());
    assert(m_paired.empty());

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
  container::svector<std::pair<std::size_t, std::size_t>> m_paired;

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

  void process_operation_queue() {
    // If two subsequent elements cancel each other, they can be erased without
    // having to worry about potentially violating the stack-like ordering of
    // memory operations (it will be preserved).
    for (std::size_t i = 0; i < m_queue.size() - 1; ++i) {
      if (m_queue[i]->cancels(m_queue.at(i + 1))) {
        m_queue.erase(m_queue.begin() + i + 1);
        m_queue.erase(m_queue.begin() + i);
        // Re-visit the previous value of i in the next iteration
        // (keep in mind that i gets incremented at the end of the loop)
        static_assert(static_cast<decltype(i)>(-1) + 1 == 0);
        i = std::max(i - 2, static_cast<decltype(i)>(-1));
      }
    }

    // If, within the given block, we have two operations that cancel each other
    // but which do not appear in subsequent order, we **might** be able to
    // remove them both. However, this will affect mess up the stack-like
    // ordering of memory operations in case that the cancelled operations are
    // themselves memory operations. Hence, some care needs to be employed when
    // attempting this. At this point, we don't have enough context information
    // available to make this choice as this will require a bidirectional flow
    // of information through the operation chain. This is because we are
    // attempting to remove a pair of operations. See also:
    // Saabas & Uustalu, Electron. Notes Theor. Comput. Sci., 190 (2007)
    // DOI: 10.1016/j.entcs.2007.02.063
    for (std::size_t i = 0; i < m_queue.size() - 1; ++i) {
      const Operation &first = m_queue.at(i);

      for (std::size_t k = i + 1; k < m_queue.size(); ++k) {
        const Operation &second = m_queue.at(k);

        if (first->pairs_with(second)) {
          auto pair = std::make_pair(i + m_cache.size(), k + m_cache.size());
          assert(std::ranges::find(m_paired, pair.first,
                                   &decltype(m_paired)::value_type::first) ==
                 m_paired.end());
          assert(std::ranges::find(m_paired, pair.second,
                                   &decltype(m_paired)::value_type::second) ==
                 m_paired.end());

          m_paired.push_back(std::move(pair));
          break;
        }
      }
    }

    m_cache.insert(m_cache.end(), std::make_move_iterator(m_queue.begin()),
                   std::make_move_iterator(m_queue.end()));
    m_queue.clear();
  }

  void optimize_operation_cache() {
    std::size_t erased = 0;

    // Process paired operations that **might** be removed, provided we can (and
    // want to) reorder other operations to retain a consistent stack.
    // A rigorous optimization would likely require methods as those described
    // in the following publication:
    // Saabas & Uustalu, Electron. Notes Theor. Comput. Sci., 190 (2007)
    // DOI: 10.1016/j.entcs.2007.02.063
    // In particular, this would almost certainly require explicit tracking of
    // the memory stack during the individual operations.
    for (auto [first_idx, second_idx] : m_paired) {
      // Account for previously erased pairs (this simple method only works in
      // combination with the assumptions under which we are currently doing
      // erasures - see below)
      assert(first_idx >= erased);
      assert(second_idx >= erased);
      first_idx -= erased;
      second_idx -= erased;

      assert(first_idx < m_cache.size());
      assert(second_idx < m_cache.size());
      assert(first_idx < second_idx);

      [[maybe_unused]] const Operation &first = m_cache.at(first_idx);
      [[maybe_unused]] const Operation &second = m_cache.at(second_idx);
      assert(first->pairs_with(second));

      if (second_idx - first_idx > 2) {
        // There is more than one intermittent operation between the pair. We
        // skip these cases for now as the required reordering seems too
        // complicated for the time being.
        continue;
      }

      const Operation &intermittent = m_cache.at(first_idx + 1);
      const MemoryAction intermittent_action = intermittent->memory_action();
      if (intermittent_action == MemoryAction::None) {
        // This case requires some more thought -> skip for now
        continue;
      }

      if (intermittent_action == MemoryAction::Allocate) {
        // Should work similar to Deallocate but hasn't been implemented yet
        continue;
      }

      assert(intermittent_action == MemoryAction::Deallocate);
      // If this was an allocation, the stack model would already be violated
      assert(first->memory_action() == MemoryAction::Deallocate);

      auto alloc_it = find_unpaired_allocation(
          m_cache.rbegin() + m_cache.size() - 1 - first_idx + 1,
          m_cache.rend());
      assert(alloc_it != m_cache.rend());
      assert(std::distance(alloc_it, m_cache.rend()) > 0);
      std::size_t first_alloc_idx = std::distance(alloc_it, m_cache.rend()) - 1;
      assert(m_cache.at(first_alloc_idx)->pairs_with(first));

      alloc_it = find_unpaired_allocation(alloc_it + 1, m_cache.rend());
      assert(alloc_it != m_cache.rend());
      assert(std::distance(alloc_it, m_cache.rend()) > 0);
      std::size_t intermittent_alloc_idx =
          std::distance(alloc_it, m_cache.rend()) - 1;
      assert(m_cache.at(intermittent_alloc_idx)->pairs_with(intermittent));

      assert(intermittent_alloc_idx < first_alloc_idx);
      if (first_alloc_idx - intermittent_alloc_idx == 1) {
        // The allocations happen in subsequently -> we can simply swap them
        std::swap(m_cache.at(intermittent_alloc_idx),
                  m_cache.at(first_alloc_idx));
        m_cache.erase(m_cache.begin() + second_idx);
        m_cache.erase(m_cache.begin() + first_idx);
        erased += 2;
        continue;
      }

      // At this point we could do a more sophisticated analysis of whether or
      // not to reorder operations in such a way that would allow for the
      // current pair of operations to be erased. However, this would require
      // sophisticated analysis including a cost function to decide whether or
      // not to move an allocation before a computation that doesn't need the
      // allocated tensor. So there is a tradeoff between increased memory usage
      // of the generated code and reduced I/O redundancies.
    }
  }

  void process_operation_cache(const MainContext &ctx) {
    assert(m_queue.empty());
    assert(m_cache.empty() || m_cache.front()->pairs_with(m_cache.back()));

    optimize_operation_cache();

    for (Operation &op : m_cache) {
      op->execute(m_generator, ctx);
    }

    m_cache.clear();
    m_paired.clear();
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
