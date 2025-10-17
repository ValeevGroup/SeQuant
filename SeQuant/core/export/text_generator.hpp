#ifndef SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP
#define SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

#include <sstream>
#include <string>
#include <string_view>

namespace sequant {

/// Context for the TextGenerator
struct TextGeneratorContext : ExportContext {};

/// A dummy generator producing a plain text representation of the code. Mostly
/// intended for having a convenient backend for tests available but it is also
/// very useful for debugging or other cases in which a human-readable version
/// of the code is required.
template <typename Context = TextGeneratorContext>
class TextGenerator : public Generator<Context> {
 public:
  TextGenerator() = default;
  ~TextGenerator() = default;

  std::string get_format_name() const override { return "Plain Text"; }

  bool supports_named_sections() const override { return true; }

  bool requires_named_sections() const override { return false; }

  DeclarationScope index_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope variable_declaration_scope() const override {
    return DeclarationScope::Section;
  }

  DeclarationScope tensor_declaration_scope() const override {
    return DeclarationScope::Section;
  }

  std::string represent(const Index &idx, const Context &) const override {
    return toUtf8(idx.label());
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = toUtf8(tensor.label()) + "[";

    using namespace std::literals;
    representation += tensor.const_indices() |
                      ranges::views::transform([&](const Index &idx) {
                        return represent(idx, ctx);
                      }) |
                      ranges::views::join(", "s) | ranges::to<std::string>();

    representation += "]";

    return representation;
  }

  std::string represent(const Variable &variable,
                        const Context &) const override {
    return toUtf8(variable.label());
  }

  std::string represent(const Constant &constant,
                        const Context &) const override {
    std::stringstream sstream;
    if (constant.value().imag() != 0) {
      sstream << "(" << constant.value().real();
      if (constant.value().imag() < 0) {
        sstream << " - i" << (-1 * constant.value().imag());
      } else {
        sstream << " + i" << constant.value().imag();
      }
    } else {
      sstream << constant.value().real();
    }
    return sstream.str();
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    m_generated += m_indent + "Create " + represent(tensor, ctx);
    if (zero_init) {
      m_generated += " and initialize to zero";
    }
    m_generated += "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    m_generated += m_indent + "Load " + represent(tensor, ctx);

    if (set_to_zero) {
      m_generated += " and set it to zero";
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    m_generated +=
        m_indent + "Setting " + represent(tensor, ctx) + " to zero\n";
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += m_indent + "Unload " + represent(tensor, ctx) + "\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    m_generated += m_indent + "Delete " + represent(tensor, ctx) + "\n";
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += m_indent + "Persist " + represent(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    m_generated += m_indent + "Create " + represent(variable, ctx);
    if (zero_init) {
      m_generated += " and initialize to zero";
    }
    m_generated += "\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += m_indent + "Load " + represent(variable, ctx);
    if (set_to_zero) {
      m_generated += " and set it to zero";
    }
    m_generated += "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated +=
        m_indent + "Setting " + represent(variable, ctx) + " to zero\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += m_indent + "Unload " + represent(variable, ctx) + "\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    m_generated += m_indent + "Delete " + represent(variable, ctx) + "\n";
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_generated += m_indent + "Persist " + represent(variable, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += m_indent + "Compute " + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_generated += m_indent + "Compute " + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    m_generated += m_indent + "Declare index " + represent(idx, ctx) + "\n";
  }

  void declare(const Variable &variable, UsageSet /*usage*/,
               const Context &ctx) override {
    if (ctx.inside_named_section()) {
      if (m_generated.back() != '(') {
        m_generated += ", ";
      }

      m_generated += "variable " + represent(variable, ctx);
    } else {
      m_generated +=
          m_indent + "Declare variable " + represent(variable, ctx) + "\n";
    }
  }

  void declare(const Tensor &tensor, UsageSet /*usage*/,
               const Context &ctx) override {
    if (ctx.inside_named_section()) {
      if (m_generated.back() != '(') {
        m_generated += ", ";
      }

      m_generated += "tensor " + represent(tensor, ctx);
    } else {
      m_generated +=
          m_indent + "Declare tensor " + represent(tensor, ctx) + "\n";
    }
  }

  void all_indices_declared(std::size_t amount, const Context &) override {
    if (amount > 0) {
      m_generated += "\n";
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

  void begin_declarations(DeclarationScope, const Context &) override {}

  void end_declarations(DeclarationScope scope, const Context &ctx) override {
    if (scope == DeclarationScope::Section && ctx.inside_named_section()) {
      m_generated += ")\n";
    }
  }

  void insert_comment(const std::string &comment, const Context &) override {
    m_generated += m_indent + "// " + comment + "\n";
  }

  void begin_named_section(std::string_view name, const Context &) override {
    m_generated += m_indent + "section " + std::string(name) + "(";
    m_indent += "  ";
  }

  void end_named_section(std::string_view /*name*/, const Context &) override {
    SEQUANT_ASSERT(m_indent.size() >= 2);
    m_indent = m_indent.substr(2, std::string::npos);
    m_generated += m_indent + "end section\n";
  }

  void begin_expression(const Context &) override {
    if (!m_generated.empty() && !m_generated.ends_with("\n\n") &&
        !m_generated.ends_with(")\n")) {
      m_generated += "\n";
    }
  }

  void begin_export(const Context &) override { m_generated.clear(); }

  void end_export(const Context &) override {}

  void end_expression(const Context &) override {}

  std::string get_generated_code() const override { return m_generated; }

 private:
  std::string m_generated;
  std::string m_indent;

  std::string stringify(const Expr &expr, const Context &ctx) const {
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
        repr += represent(Constant(product.scalar()), ctx) + " ";
      }

      for (std::size_t i = 0; i < product.size(); ++i) {
        repr += stringify(*product.factor(i), ctx);

        if (i + 1 < product.size()) {
          repr += " ";
        }
      }

      return repr;
    } else if (expr.is<Sum>()) {
      const Sum &sum = expr.as<Sum>();
      std::string repr;

      for (std::size_t i = 0; i < sum.size(); ++i) {
        repr += stringify(*sum.summand(i), ctx);

        if (i + 1 < sum.size()) {
          repr += " + ";
        }
      }

      return repr;
    }

    throw std::runtime_error(
        "Unsupported expression type in TextGenerator::compute");
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP
