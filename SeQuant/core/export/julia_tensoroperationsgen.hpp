#ifndef SEQUANT_CORE_EXPORT_JULIATENSOROPERATIONSGEN_HPP
#define SEQUANT_CORE_EXPORT_JULIATENSOROPERATIONSGEN_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>

namespace sequant {

class JuliaTensorOperationsGenContext : public ExportContext {
 public:
  JuliaTensorOperationsGenContext() = default;
  ~JuliaTensorOperationsGenContext() = default;
  JuliaTensorOperationsGenContext(std::map<IndexSpace, std::string> index_tags,
                                  std::map<IndexSpace, std::string> index_dims,
                                  bool print_intermediate_comments);

  std::string get_tag(const IndexSpace &indexspace) const;

  std::string get_dim(const IndexSpace &indexspace) const;

  std::string get_tags(const Tensor &tensor) const;

  bool print_intermediate_comments() const;

 protected:
  std::map<IndexSpace, std::string> m_index_tags;
  std::map<IndexSpace, std::string> m_index_dims;
  bool m_print_intermediate_comments;
};

template <typename Context = JuliaTensorOperationsGenContext>
class JuliaTensorOperationsGen : public Generator<Context> {
 public:
  JuliaTensorOperationsGen() = default;
  ~JuliaTensorOperationsGen() = default;

  std::string get_format_name() const override {
    return "Julia TensorOperations";
  }

  std::string represent(const Index &idx, const Context &ctx) const override {
    std::string labelstring = toUtf8(idx.full_label());
    if (idx.has_proto_indices()) {
      throw std::runtime_error(
          "Incompatible (proto) index detected. Proto Indices are not (yet) "
          "supported!");
    }
    return labelstring;
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    const auto &indices = tensor.const_braket();
    std::string representation = represent_tensor_noindex(tensor, ctx);
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
          "Creation without initialization is not supported.");
    }
    m_generated += represent_tensor_noindex(tensor, ctx);
    const auto &indices = tensor.const_braket();
    std::vector<std::string> dimStrings(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      dimStrings[i] = ctx.get_dim(indices[i].space());
    }
    std::string rhs = "zeros(Float64";
    for (std::size_t i = 0; i < indices.size(); ++i) {
      rhs += ", " + dimStrings[i];
    }
    rhs += ")";
    m_generated += " = " + rhs + "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (!set_to_zero) {
      m_generated += represent_tensor_noindex(tensor, ctx) +
                     " = deserialize(\"" +
                     represent_tensor_noindex(tensor, ctx) + ".jlbin\")";
    }

    else {
      m_generated += represent_tensor_noindex(tensor, ctx);
      const auto &indices = tensor.const_braket();
      std::vector<std::string> dimStrings(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        dimStrings[i] = ctx.get_dim(indices[i].space());
      }
      std::string rhs = "zeros(Float64";
      for (std::size_t i = 0; i < indices.size(); ++i) {
        rhs += ", " + dimStrings[i];
      }
      rhs += ")";
      m_generated += " = " + rhs;
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    // In Julia setting to zero essentially has same effect as creating a zero
    // tensor
    create(tensor, true, ctx);
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += represent_tensor_noindex(tensor, ctx) + " = nothing\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    m_generated += represent_tensor_noindex(tensor, ctx) + " = nothing\n";
    m_generated +=
        "rm( \"" + represent_tensor_noindex(tensor, ctx) + +".jlbin\"" + ")\n";
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += "return " + represent_tensor_noindex(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = 0.0";
    m_generated += "\n";
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
    m_generated += represent(variable, ctx) + " = nothing\n";
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
    if (ctx.print_intermediate_comments()) {
      std::string comment_line = "Declare index " + represent(idx, ctx);
      insert_comment(comment_line, ctx);
    }
  }

  void declare(const Variable &variable, const Context &ctx) override {
    if (ctx.print_intermediate_comments()) {
      std::string comment_line = "Declare Variable " + represent(variable, ctx);
      insert_comment(comment_line, ctx);
    }
  }

  void declare(const Tensor &tensor, const Context &ctx) override {
    if (ctx.print_intermediate_comments()) {
      std::string representation = toUtf8(tensor.label());
      const auto &indices = tensor.const_braket();
      std::vector<std::string> dimStrings(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        dimStrings[i] = ctx.get_dim(indices[i].space());
      }

      std::string rhs = "zeros(Float64";
      for (std::size_t i = 0; i < indices.size(); ++i) {
        rhs += ", " + dimStrings[i];
      }
      rhs += ")";
      representation += " = " + rhs;
      m_generated += "#" + representation + "\n";
    }
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    m_generated += "# " + comment + "\n";
  }

  void insert_blank_lines(std::size_t count, const Context &ctx) override {
    for (std::size_t i = 0; i < count; ++i) {
      m_generated += "\n";
    }
  }

  std::string get_generated_code() const override { return m_generated; }

 protected:
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

  std::string represent_tensor_noindex(const Tensor &tensor,
                                       const Context &ctx) const {
    std::string representation = toUtf8(tensor.label());
    std::string tagstring = ctx.get_tags(tensor);
    representation += "_";
    representation += tagstring;
    return representation;
  }

  std::string m_generated;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIATENSOROPERATIONSGEN_HPP
