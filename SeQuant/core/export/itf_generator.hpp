#ifndef SEQUANT_CORE_EXPORT_ITFGENERATOR_HPP
#define SEQUANT_CORE_EXPORT_ITFGENERATOR_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cassert>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>

namespace sequant {

class ItfGeneratorContext : public ExportContext {
 public:
  using TagMap = std::map<IndexSpace, std::string>;
  using NameMap = std::map<IndexSpace, std::string>;

  ItfGeneratorContext() = default;
  ~ItfGeneratorContext() = default;

  std::string index_name(const IndexSpace &space, std::size_t ordinal) const;

  std::string get_name(const IndexSpace &space) const;

  std::string get_tag(const IndexSpace &space) const;

  std::optional<std::string> import_name(const Tensor &tensor) const;
  std::optional<std::string> import_name(const Variable &variable) const;

 private:
  TagMap m_index_tags;
  NameMap m_space_names;
};

template <typename Context = ItfGeneratorContext>
class ItfGenerator : public Generator<Context> {
 public:
  ItfGenerator() = default;
  ~ItfGenerator() = default;

  std::string get_format_name() const override { return "ITF"; }

  bool supports_named_sections() const override { return true; }

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
      throw std::runtime_error("ITF doesn't support proto indices");
    }

    const std::size_t ordinal = idx.ordinal();

    std::string name = ctx.index_name(idx.space(), ordinal);
    assert(!name.empty());

    if (name.size() > 1) {
      name = "{" + name + "}";
    }

    return name;
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = toUtf8(tensor.label());

    const auto &indices = tensor.const_indices();

    if (!indices.empty()) {
      representation += ":";
      for (const Index &idx : indices) {
        representation += ctx.get_tag(idx.space());
      }
    }

    representation += "[";

    for (std::size_t i = 0; i < indices.size(); ++i) {
      representation += represent(indices[i], ctx);
    }

    representation += "]";

    return representation;
  }

  std::string represent(const Variable &variable,
                        const Context &ctx) const override {
    return toUtf8(variable.label()) + "[]";
  }

  std::string represent(const Constant &constant,
                        const Context &ctx) const override {
    std::stringstream sstream;
    if (constant.value().imag() != 0) {
      sstream << "(" << constant.value().real();
      if (constant.value().imag() < 0) {
        sstream << " - j" << (-1 * constant.value().imag());
      } else {
        sstream << " + j" << constant.value().imag();
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
          "Can't create ITF tensor without setting it to zero");
    }

    m_generated += "alloc " + represent(tensor, ctx) + "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (set_to_zero) {
      create(tensor, true, ctx);
      return;
    }

    m_generated += "load " + represent(tensor, ctx) + "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    m_generated += "." + represent(tensor, ctx) + " := 0 * One[]\n";
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += "drop " + represent(tensor, ctx) + "\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    // There exists no real notion of deleting/destroying a given tensor
    unload(tensor, ctx);
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += "store " + represent(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Can't create ITF variable without setting it to zero");
    }

    m_generated += "alloc " + represent(variable, ctx) + "\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += "load " + represent(variable, ctx) + "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated += "." + represent(variable, ctx) + " := 0 * One[]\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += "drop " + represent(variable, ctx) + "\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    // There exists no real notion of deleting/destroying a given variable
    unload(variable, ctx);
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_generated += "store " + represent(variable, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += "." + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_generated += "." + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    m_generated += "index-space: " + represent(idx, ctx) + ", " +
                   ctx.get_name(idx.space()) + ", " + ctx.get_tag(idx.space()) +
                   "\n";
  }

  void declare(const Variable &variable, const Context &ctx) override {
    m_generated += "tensor: " + represent(variable, ctx) + ", ";

    std::optional<std::string> import_name = ctx.import_name(variable);
    if (import_name.has_value()) {
      m_generated += import_name.value();
    } else {
      m_generated += "!Create{type:scalar}";
    }

    m_generated += "\n";
  }

  void declare(const Tensor &tensor, const Context &ctx) override {
    m_generated += "tensor: " + represent(tensor, ctx) + ", ";

    std::optional<std::string> import_name = ctx.import_name(tensor);
    if (import_name.has_value()) {
      m_generated += import_name.value();
    } else {
      // TODO: we can use type:plain for tensors that don't need to be persisted
      m_generated += "!Create{type:disk}";
    }

    m_generated += "\n";
  }

  void all_indices_declared(std::size_t amount, const Context &ctx) override {
    if (amount > 0) {
      m_generated += "\n\n";
    }
  }

  void all_variables_declared(std::size_t amount, const Context &ctx) override {
    if (amount > 0 && !ctx.inside_named_section()) {
      m_generated += "\n";
    }
  }

  void all_tensors_declared(std::size_t amount, const Context &ctx) override {
    if (amount > 0 && !ctx.inside_named_section()) {
      m_generated += "\n";
    }
  }

  void begin_declarations(DeclarationScope scope, const Context &ctx) override {
    if (scope == DeclarationScope::Global) {
      m_generated += "---- decl\n";
    }
  }

  void end_declarations(DeclarationScope scope, const Context &ctx) override {
    if (scope == DeclarationScope::Global) {
      m_generated += "\n\n";
    }
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    m_generated += "// " + comment + "\n";
  }

  void begin_named_section(std::string_view name, const Context &ctx) override {
    m_generated += "---- code(\"" + std::string(name) + "\")\n";
  }

  void end_named_section(std::string_view name, const Context &ctx) override {
    m_generated += "\n\n";
  }

  void begin_expression(const Context &ctx) override {
    if (!m_generated.empty() && !m_generated.ends_with("\n\n") &&
        !m_generated.ends_with(")\n")) {
      m_generated += "\n";
    }
  }

  void end_expression(const Context &ctx) override {}

  void begin_export(const Context &ctx) override { m_generated.clear(); }

  void end_export(const Context &ctx) override { m_generated += "---- end\n"; }

  std::string get_generated_code() const override { return m_generated; }

 private:
  std::string m_generated;

  std::string stringify(const Expr &expr, const Context &ctx) const {
    if (expr.is<Tensor>()) {
      return represent(expr.as<Tensor>(), ctx);
    } else if (expr.is<Variable>()) {
      return represent(expr.as<Variable>(), ctx);
    } else if (expr.is<Constant>()) {
      return represent(expr.as<Constant>(), ctx);
    } else if (expr.is<Product>()) {
      const Product &product = expr.as<Product>();

      if (product.factors().size() > 2) {
        throw std::runtime_error("ITF can only handle binary contractions");
      }

      std::string repr;

      if (!product.scalar().is_identity()) {
        repr += represent(Constant(product.scalar()), ctx) + " * ";
      }

      for (std::size_t i = 0; i < product.size(); ++i) {
        repr += stringify(*product.factor(i), ctx);

        if (i + 1 < product.size()) {
          repr += " ";
        }
      }

      return repr;
    } else if (expr.is<Sum>()) {
      // TODO: Handle on-the-fly antisymmetrization (K[abij] - K[baij])
      throw std::runtime_error("ITF doesn't support explicit summation");
    }

    throw std::runtime_error(
        "Unsupported expression type in ItfGenerator::compute");
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_ITFGENERATOR_HPP
