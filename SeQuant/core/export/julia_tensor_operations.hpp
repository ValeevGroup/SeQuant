#ifndef SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP
#define SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstdlib>
#include <string>
#include <string_view>

namespace sequant {

/// Context class for the @see JuliaTensorOperationsGenerator generator
class JuliaTensorOperationsGeneratorContext : public ExportContext {
 public:
  using TagMap = std::map<IndexSpace, std::string>;
  using DimMap = std::map<IndexSpace, std::string>;

  JuliaTensorOperationsGeneratorContext() = default;
  ~JuliaTensorOperationsGeneratorContext() = default;
  JuliaTensorOperationsGeneratorContext(TagMap index_tags, DimMap index_dims);

  std::string get_tag(const IndexSpace &space) const;

  std::string get_dim(const IndexSpace &space) const;

  std::string get_tags(const Tensor &tensor) const;

  void set_tag(const IndexSpace &space, std::string tag);
  void set_dim(const IndexSpace &space, std::string dim);

 protected:
  TagMap m_index_tags;
  DimMap m_index_dims;
};

/// Code generator for the Julia language, which will make use of the
/// TensorOperations.jl package for tensor contractions
template <typename Context = JuliaTensorOperationsGeneratorContext>
class JuliaTensorOperationsGenerator : public Generator<Context> {
 public:
  JuliaTensorOperationsGenerator() = default;
  ~JuliaTensorOperationsGenerator() = default;

  std::string get_format_name() const override {
    return "Julia (TensorOperations)";
  }

  bool supports_named_sections() const override { return false; }

  bool requires_named_sections() const override { return false; }

  DeclarationScope index_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope variable_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope tensor_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  std::string represent(const Index &idx, const Context &ctx) const override {
    if (idx.has_proto_indices()) {
      throw std::runtime_error("Proto Indices are not (yet) supported!");
    }

    return toUtf8(idx.full_label());
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    const auto &indices = tensor.const_indices();
    std::string representation = tensor_name(tensor, ctx);

    representation += "[ ";

    for (std::size_t i = 0; i < indices.size(); ++i) {
      representation += represent(indices[i], ctx);

      if (i + 1 < indices.size()) {
        representation += ", ";
      }
    }

    representation += " ]";

    return representation;
  }

  std::string represent(const Variable &variable,
                        const Context &ctx) const override {
    return toUtf8(variable.label());
  }

  std::string represent(const Constant &constant,
                        const Context &ctx) const override {
    std::stringstream sstream;

    if (constant.value().imag() != 0) {
      sstream << "(" << constant.value().real();
      if (constant.value().imag() < 0) {
        sstream << "-" << (-1 * constant.value().imag()) << "im)";
      } else {
        sstream << "+" << constant.value().imag() << "im)";
      }
    } else {
      sstream << constant.value().real();
    }

    return sstream.str();
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "In Julia tensors can't be created without being initialized");
    }

    m_generated += zero_initialization(tensor, ctx) + "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    m_generated += tensor_name(tensor, ctx);
    m_generated += " = ";

    if (!set_to_zero) {
      m_generated += "deserialize(\"" + tensor_name(tensor, ctx) + ".jlbin\")";
    } else {
      m_generated += zero_initialization(tensor, ctx);
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    // In Julia setting to zero essentially has same effect as creating a zero
    // tensor
    create(tensor, true, ctx);
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += tensor_name(tensor, ctx) + " = nothing\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    unload(tensor, ctx);
    // Remove potential serialized version of this tensor from disk
    m_generated += "rm( \"" + tensor_name(tensor, ctx) + +".jlbin\"" + ")\n";
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += "return " + tensor_name(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Julia doesn't support declaring a variable without initializing it");
    }

    m_generated += represent(variable, ctx) + " = 0.0\n";
  }
  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = deserialize(\"" +
                   represent(variable, ctx) + ".jlbin\")";
    m_generated += "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = 0.0\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = nothing\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    unload(variable, ctx);
    // Remove any potentially serialized version from disk
    m_generated += "rm(\"" + represent(variable, ctx) + ".jlbin\")\n";
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_generated += "return " + represent(variable, ctx) + "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_generated += "@tensor " + represent(result, ctx) +
                   " += " + to_julia_expr(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += "@tensor " + represent(result, ctx) +
                   " += " + to_julia_expr(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    (void)idx;
    (void)ctx;
  }

  void declare(const Variable &variable, UsageSet usage,
               const Context &ctx) override {
    (void)variable;
    (void)ctx;
  }

  void declare(const Tensor &tensor, UsageSet usage,
               const Context &ctx) override {
    (void)tensor;
    (void)ctx;
  }

  void all_indices_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void all_variables_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void all_tensors_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void begin_declarations(DeclarationScope scope, const Context &ctx) override {
    (void)scope;
    (void)ctx;
  }

  void end_declarations(DeclarationScope scope, const Context &ctx) override {
    (void)scope;
    (void)ctx;
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    m_generated += "# " + comment + "\n";
  }

  void begin_named_section(std::string_view name, const Context &ctx) override {
    std::abort();
  }

  void end_named_section(std::string_view name, const Context &ctx) override {
    std::abort();
  }

  void begin_expression(const Context &ctx) override {
    if (!m_generated.empty() && !m_generated.ends_with("\n\n")) {
      m_generated += "\n";
    }
  }

  void end_expression(const Context &ctx) override {}

  void begin_export(const Context &ctx) override { m_generated.clear(); }

  void end_export(const Context &ctx) override {}

  std::string get_generated_code() const override { return m_generated; }

 protected:
  std::string m_generated;

  std::string to_julia_expr(const Expr &expr, const Context &ctx) const {
    if (expr.is<Tensor>()) {
      return represent(expr.as<Tensor>(), ctx);
    } else if (expr.is<Variable>()) {
      return represent(expr.as<Variable>(), ctx);
    } else if (expr.is<Constant>()) {
      return represent(expr.as<Constant>(), ctx);
    } else if (expr.is<Product>()) {
      const Product &product = expr.as<Product>();
      std::string repr;

      if (!product.scalar().is_identity()) {
        repr += represent(Constant(product.scalar()), ctx) + " * ";
      }

      for (std::size_t i = 0; i < product.size(); ++i) {
        repr += to_julia_expr(*product.factor(i), ctx);

        if (i + 1 < product.size()) {
          repr += " * ";
        }
      }

      return repr;
    } else if (expr.is<Sum>()) {
      const Sum &sum = expr.as<Sum>();
      std::string repr;

      for (std::size_t i = 0; i < sum.size(); ++i) {
        repr += to_julia_expr(*sum.summand(i), ctx);

        if (i + 1 < sum.size()) {
          repr += " + ";
        }
      }

      return repr;
    }

    throw std::runtime_error("Unsupported expression type in to_julia_expr");
  }

  std::string tensor_name(const Tensor &tensor, const Context &ctx) const {
    std::string representation = toUtf8(tensor.label());

    representation += "_";
    representation += ctx.get_tags(tensor);

    return representation;
  }

  std::string zero_initialization(const Tensor &tensor,
                                  const Context &ctx) const {
    std::string code = tensor_name(tensor, ctx);
    code += " = zeros(Float64";

    for (const Index &idx : tensor.const_indices()) {
      code += ", ";
      code += ctx.get_dim(idx.space());
    }

    code += ")";

    return code;
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP
