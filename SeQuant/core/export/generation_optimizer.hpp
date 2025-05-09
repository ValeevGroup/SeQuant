#ifndef SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
#define SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <pv/polymorphic_variant.hpp>

#include <cassert>
#include <ranges>
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

    bool erases(const AbstractOperation &other) const {
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
          return false;
      }

      // should not be reached
      assert(false);
      return false;
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

  using Operation =
      pv::polymorphic_variant<AbstractOperation, CreateOperation, LoadOperation,
                              ZeroOperation, UnloadOperation, DestroyOperation,
                              PersistOperation>;

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
    // Assumption: Context is the same as for all queued operations
    process_operation_queue(*this, ctx);
    Generator::compute(expression, result, ctx);
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    // Assumption: Context is the same as for all queued operations
    process_operation_queue(*this, ctx);
    Generator::compute(expression, result, ctx);
  }

  void end_expression(const Context &ctx) override {
    // Assumption: Context is the same as for all queued operations
    process_operation_queue(*this, ctx);
    Generator::end_expression(ctx);
  }

  void end_export(const Context &ctx) override {
    assert(m_queue.empty());
    Generator::end_export(ctx);
  }

 private:
  container::svector<Operation> m_queue;
  container::svector<container::svector<Operation>> m_blocks;

  void process_operation_queue(Generator &generator, const Context &ctx) {
    // Remove operations that cancel each other
    for (std::size_t i = 0; i < m_queue.size(); ++i) {
      const Operation &first = m_queue.at(i);

      for (std::size_t k = i; k < m_queue.size(); ++k) {
        const Operation &second = m_queue.at(k);

        if (first->erases(second)) {
          m_queue.erase(m_queue.begin() + k);
          m_queue.erase(m_queue.begin() + i);
          i--;
          break;
        }
      }
    }

    for (Operation &operation : m_queue) {
      for (container::svector<Operation> &block :
           std::ranges::reverse_view(m_blocks)) {
        auto it = std::find_if(
            block.rbegin(), block.rend(),
            [&](const Operation &op) { return op.pairs_with(operation); });

        if (it == block.rend()) {
          continue;
        }

		auto distance = std::distance(block.rbegin(), it);
		if (distance > 0) {
			// We might have a load order conflict
			// Also true if we find the pairing for operation in a block that is not m_blocks.back()
		}

        Operation op = std::move(*it);
        // TODO: Execute op immediately and operation after the computation that
        // has triggered this processing

        block.erase(it);
        break;
      }

      // TODO: Removing an empty block should trigger execution of all
      // computations that depend on this block (only)
      std::remove_if(m_blocks.begin(), m_blocks.end(),
                     [](const container::svector<Operation> &block) {
                       return block.empty();
                     });
    }

    // for (Operation &operation : m_queue) {
    //  operation->execute(generator, ctx);
    //}

    if (!m_queue.empty()) {
      m_blocks.push_back(std::move(m_queue));
      m_queue.clear();
    }
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_GENERATIONOPTIMIZER_HPP
