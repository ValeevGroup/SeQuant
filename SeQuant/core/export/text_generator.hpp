#ifndef SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP
#define SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <sstream>
#include <string>

namespace sequant {

struct TextGeneratorContext : ExportContext {};

template <typename Context = TextGeneratorContext>
class TextGenerator : public Generator<Context> {
 public:
  TextGenerator() = default;
  ~TextGenerator() = default;

  std::string get_format_name() const override { return "Plain Text"; }

  std::string represent(const Index &idx, const Context &ctx) const override {
    return toUtf8(idx.label());
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = toUtf8(tensor.label()) + "[";
    const auto &indices = tensor.const_braket();

    for (std::size_t i = 0; i < indices.size(); ++i) {
      representation += represent(indices[i], ctx);

      if (i + 1 < indices.size()) {
        representation += ", ";
      }
    }

    representation += "]";

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
    m_generated += "Create " + represent(tensor, ctx);
    if (zero_init) {
      m_generated += " and initialize to zero";
    }
    m_generated += "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    m_generated += "Load " + represent(tensor, ctx);

    if (set_to_zero) {
      m_generated += " and set it to zero";
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    m_generated += "Setting " + represent(tensor, ctx) + " to zero\n";
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += "Unload " + represent(tensor, ctx) + "\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    m_generated += "Delete " + represent(tensor, ctx) + "\n";
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += "Persist " + represent(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    m_generated += "Create " + represent(variable, ctx);
    if (zero_init) {
      m_generated += " and initialize to zero";
    }
    m_generated += "\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += "Load " + represent(variable, ctx);
    if (set_to_zero) {
      m_generated += " and set it to zero";
    }
    m_generated += "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated += "Setting " + represent(variable, ctx) + " to zero\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += "Unload " + represent(variable, ctx) + "\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    m_generated += "Delete " + represent(variable, ctx) + "\n";
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_generated += "Persist " + represent(variable, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += "Compute " + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_generated += "Compute " + represent(result, ctx) +
                   " += " + stringify(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    m_generated += "Declare index " + represent(idx, ctx) + "\n";
  }

  void declare(const Variable &variable, const Context &ctx) override {
    m_generated += "Declare variable " + represent(variable, ctx) + "\n";
  }

  void declare(const Tensor &tensor, const Context &ctx) override {
    m_generated += "Declare tensor " + represent(tensor, ctx) + "\n";
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    m_generated += "// " + comment + "\n";
  }

  void insert_blank_lines(std::size_t count, const Context &ctx) override {
    for (std::size_t i = 0; i < count; ++i) {
      m_generated += "\n";
    }
  }

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
