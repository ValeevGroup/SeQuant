#ifndef SEQUANT_CORE_EXPORT_JULIA_ITENSOR_HPP
#define SEQUANT_CORE_EXPORT_JULIA_ITENSOR_HPP

#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

namespace sequant {

/// Context class for the JuliaITensorGenerator
class JuliaITensorGeneratorContext
    : public JuliaTensorOperationsGeneratorContext {
 public:
  using JuliaTensorOperationsGeneratorContext::
      JuliaTensorOperationsGeneratorContext;
};

/// Generator for producing Julia code using the ITensor framework
/// - https://julialang.org/
/// - https://itensor.github.io/ITensors.jl
template <typename Context = JuliaITensorGeneratorContext>
class JuliaITensorGenerator : public JuliaTensorOperationsGenerator<Context> {
 private:
  using Base = JuliaTensorOperationsGenerator<Context>;

 public:
  JuliaITensorGenerator() = default;
  ~JuliaITensorGenerator() = default;

  std::string get_format_name() const override { return "Julia (ITensor)"; }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    return Base::tensor_name(tensor, ctx);
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "In Julia tensors can't be created without being initialized");
    }

    Base::m_generated += Base::tensor_name(tensor, ctx);
    Base::m_generated += " = ITensor(zeros(Float64";

    for (const Index &idx : tensor.const_indices()) {
      Base::m_generated += ", ";
      Base::m_generated += ctx.get_dim(idx.space());
    }

    Base::m_generated += "), ";
    Base::m_generated += index_set(tensor, ctx);
    Base::m_generated += ")\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (set_to_zero) {
      // In Julia setting a tensor to zero has to be done by overwriting it with
      // a zero tensor
      create(tensor, true, ctx);
      return;
    }

    Base::m_generated += Base::tensor_name(tensor, ctx);
    Base::m_generated += " = ITensor(deserialize(\"";
    Base::m_generated += Base::tensor_name(tensor, ctx);
    Base::m_generated += ".jlbin\"), ";
    Base::m_generated += index_set(tensor, ctx);
    Base::m_generated += ")\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Julia doesn't support declaring a variable without initializing it");
    }

    Base::m_generated += "tmpvar = 0.0\n";
    Base::m_generated +=
        Base::represent(variable, ctx) + " = ITensor(tmpvar)\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    if (set_to_zero) {
      Base::m_generated += "tmpvar = 0\n";
    } else {
      Base::m_generated += "tmpvar = deserialize(\"" +
                           Base::represent(variable, ctx) + ".jlbin\")\n";
    }

    Base::m_generated +=
        Base::represent(variable, ctx) + " = ITensor(tmpvar)\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    Base::m_generated += Base::tensor_name(result, ctx) +
                         " += " + Base::to_julia_expr(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    Base::m_generated += Base::represent(result, ctx) +
                         " += " + Base::to_julia_expr(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    std::string currentindex = toUtf8(idx.full_label());
    Base::m_generated += currentindex;
    Base::m_generated += " = Index(";
    Base::m_generated += ctx.get_dim(idx.space());
    Base::m_generated += ", \"" + currentindex + "\")\n";
  }

  void all_indices_declared(std::size_t amount, const Context &) override {
    if (amount > 0) {
      Base::m_generated += "\n";
    }
  }

 private:
  std::string index_set(const Tensor &tensor, const Context &ctx) const {
    using namespace std::literals;
    return tensor.const_indices() |
           ranges::views::transform(
               [&](const Index &idx) { return Base::represent(idx, ctx); }) |
           ranges::views::join(", "s) | ranges::to<std::string>();
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIA_ITENSOR_HPP
