#ifndef SEQUANT_CORE_EXPORT_JULIAITENSORGEN_HPP
#define SEQUANT_CORE_EXPORT_JULIAITENSORGEN_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/julia_tensoroperationsgen.hpp>
#include <SeQuant/core/export/utils.hpp>
namespace sequant {

class JuliaITensorGenContext : public JuliaTensorOperationsGenContext {
 public:
  using JuliaTensorOperationsGenContext::JuliaTensorOperationsGenContext;
};

template <typename Context = JuliaITensorGenContext>
class JuliaITensorGen : public JuliaTensorOperationsGen<Context> {
 private:
  using Base = JuliaTensorOperationsGen<Context>;

 public:
  JuliaITensorGen() = default;
  ~JuliaITensorGen() = default;

  std::string get_format_name() const override { return "Julia ITensor"; }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = this->represent_tensor_noindex(tensor, ctx);
    return representation;
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Creation without initialization is not supported.");
    }
    this->m_generated += this->represent_tensor_noindex(tensor, ctx);
    const auto &indices = tensor.const_braket();
    std::vector<std::string> dimStrings(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      dimStrings[i] = ctx.get_dim(indices[i].space());
    }
    std::string rhs = "ITensor(zeros(Float64";
    for (std::size_t i = 0; i < indices.size(); ++i) {
      rhs += ", " + dimStrings[i];
    }
    rhs += ")";
    for (std::size_t i = 0; i < indices.size(); ++i) {
      rhs += ", " + Base::represent(indices[i], ctx);
    }
    rhs += ")";
    this->m_generated += " = " + rhs;
    if (ctx.print_intermediate_comments()) {
      this->m_generated += "\n";
    } else {
      this->m_generated += "\n";
    }
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (!set_to_zero) {
      const auto &indices = tensor.const_braket();
      this->m_generated += this->represent_tensor_noindex(tensor, ctx) +
                           " = ITensor(deserialize(\"" +
                           this->represent_tensor_noindex(tensor, ctx) +
                           ".jlbin\"),";
      for (std::size_t i = 0; i < indices.size(); ++i) {
        this->m_generated += Base::represent(indices[i], ctx);
        if (i + 1 < indices.size()) {
          this->m_generated += ", ";
        }
      }
      this->m_generated += ")";
    }

    else {
      this->m_generated += this->represent_tensor_noindex(tensor, ctx);
      const auto &indices = tensor.const_braket();
      std::vector<std::string> dimStrings(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        dimStrings[i] = ctx.get_dim(indices[i].space());
      }
      std::string rhs = "ITensor(zeros(Float64";
      for (std::size_t i = 0; i < indices.size(); ++i) {
        rhs += ", " + dimStrings[i];
      }
      rhs += ")";
      for (std::size_t i = 0; i < indices.size(); ++i) {
        rhs += ", " + Base::represent(indices[i], ctx);
      }

      rhs += ")";

      this->m_generated += " = " + rhs;
    }

    this->m_generated += "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    this->m_generated += "tmpvar = 0.0\n";
    this->m_generated += Base::represent(variable, ctx) + " = ITensor(tmpvar)";
    this->m_generated += "\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    this->m_generated += "tmpvar = deserialize(\"" +
                         Base::represent(variable, ctx) + ".jlbin\")\n";
    this->m_generated += Base::represent(variable, ctx) + " = ITensor(tmpvar)";
    this->m_generated += "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    this->m_generated += JuliaITensorGen<Context>::represent(result, ctx) +
                         " += " + Base::to_julia_expr(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    this->m_generated += Base::represent(result, ctx) +
                         " += " + Base::to_julia_expr(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    std::string currentindex = toUtf8(idx.full_label());
    this->m_generated += currentindex + "= Index(";
    this->m_generated += ctx.get_dim(idx.space());
    this->m_generated += ", \"" + currentindex + "\")\n";
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIAITENSORGEN_HPP
