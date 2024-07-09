#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <sstream>
#include <string>

namespace sequant {

TextGenerator::~TextGenerator() = default;

std::string TextGenerator::get_format_name() const { return "Plain Text"; }

std::string represent_index(const Index &idx) {
  return toUtf8(idx.full_label());
}

std::string TextGenerator::represent(const Tensor &tensor,
                                     const Context &ctx) const {
  std::string representation = toUtf8(tensor.label()) + "[";
  const auto &indices = tensor.const_braket();

  for (std::size_t i = 0; i < indices.size(); ++i) {
    representation += represent_index(indices[i]);

    if (i + 1 < indices.size()) {
      representation += ", ";
    }
  }

  representation += "]";

  return representation;
}

std::string TextGenerator::represent(const Variable &variable,
                                     const Context &ctx) const {
  return toUtf8(variable.label());
}

std::string TextGenerator::represent(const Constant &constant,
                                     const Context &ctx) const {
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

void TextGenerator::create(const Tensor &tensor, bool zero_init,
                           const Context &ctx) {
  m_generated += "Create " + represent(tensor);
  if (zero_init) {
    m_generated += " and initialize to zero";
  }
  m_generated += "\n";
}

void TextGenerator::load(const Tensor &tensor, bool set_to_zero,
                         const Context &ctx) {
  m_generated += "Load " + represent(tensor);

  if (set_to_zero) {
    m_generated += " and set it to zero";
  }

  m_generated += "\n";
}

void TextGenerator::set_to_zero(const Tensor &tensor, const Context &ctx) {
  m_generated += "Setting " + represent(tensor) + " to zero\n";
}

void TextGenerator::unload(const Tensor &tensor, const Context &ctx) {
  m_generated += "Unload " + represent(tensor) + "\n";
}

void TextGenerator::destroy(const Tensor &tensor, const Context &ctx) {
  m_generated += "Delete " + represent(tensor) + "\n";
}

void TextGenerator::persist(const Tensor &tensor, const Context &ctx) {
  m_generated += "Persist " + represent(tensor) + "\n";
}

void TextGenerator::create(const Variable &variable, bool zero_init,
                           const Context &ctx) {
  m_generated += "Create " + represent(variable);
  if (zero_init) {
    m_generated += " and initialize to zero";
  }
  m_generated += "\n";
}

void TextGenerator::load(const Variable &variable, bool set_to_zero,
                         const Context &ctx) {
  m_generated += "Load " + represent(variable);
  if (set_to_zero) {
    m_generated += " and set it to zero";
  }
  m_generated += "\n";
}

void TextGenerator::set_to_zero(const Variable &variable, const Context &ctx) {
  m_generated += "Setting " + represent(variable) + " to zero\n";
}

void TextGenerator::unload(const Variable &variable, const Context &ctx) {
  m_generated += "Unload " + represent(variable) + "\n";
}

void TextGenerator::destroy(const Variable &variable, const Context &ctx) {
  m_generated += "Delete " + represent(variable) + "\n";
}

void TextGenerator::persist(const Variable &variable, const Context &ctx) {
  m_generated += "Persist " + represent(variable) + "\n";
}

void TextGenerator::compute(const Expr &expression, const Tensor &result,
                            const Context &ctx) {
  m_generated +=
      "Compute " + represent(result) + " += " + stringify(expression) + "\n";
}

void TextGenerator::compute(const Expr &expression, const Variable &result,
                            const Context &ctx) {
  m_generated +=
      "Compute " + represent(result) + " += " + stringify(expression) + "\n";
}

void TextGenerator::declare(const Index &idx, const Context &ctx) {
  m_generated += "Declare index " + toUtf8(idx.full_label()) + "\n";
}

void TextGenerator::declare(const Variable &variable, const Context &ctx) {
  m_generated += "Declare variable " + represent(variable) + "\n";
}

void TextGenerator::declare(const Tensor &tensor, const Context &ctx) {
  m_generated += "Declare tensor " + represent(tensor) + "\n";
}

void TextGenerator::insert_comment(const std::string &comment,
                                   const Context &ctx) {
  m_generated += "// " + comment + "\n";
}

void TextGenerator::insert_blank_lines(std::size_t count, const Context &ctx) {
  for (std::size_t i = 0; i < count; ++i) {
    m_generated += "\n";
  }
}

std::string TextGenerator::get_generated_code() const { return m_generated; }

std::string TextGenerator::stringify(const Expr &expr) const {
  if (expr.is<Tensor>()) {
    return represent(expr.as<Tensor>());
  } else if (expr.is<Variable>()) {
    return represent(expr.as<Variable>());
  } else if (expr.is<Constant>()) {
    return represent(expr.as<Constant>());
  } else if (expr.is<Product>()) {
    const Product &product = expr.as<Product>();
    std::string repr;

    if (!product.scalar().is_identity()) {
      repr += represent(Constant(product.scalar())) + " ";
    }

    for (std::size_t i = 0; i < product.size(); ++i) {
      repr += stringify(*product.factor(i));

      if (i + 1 < product.size()) {
        repr += " ";
      }
    }

    return repr;
  } else if (expr.is<Sum>()) {
    const Sum &sum = expr.as<Sum>();
    std::string repr;

    for (std::size_t i = 0; i < sum.size(); ++i) {
      repr += stringify(*sum.summand(i));

      if (i + 1 < sum.size()) {
        repr += " + ";
      }
    }

    return repr;
  }

  throw std::runtime_error(
      "Unsupported expression type in TextGenerator::compute");
}

}  // namespace sequant
