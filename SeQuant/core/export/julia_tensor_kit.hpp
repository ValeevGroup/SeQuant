#ifndef SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP
#define SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP

#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

namespace sequant {

// Context class for the JuliaTensorKitGenerator
class JuliaTensorKitGeneratorContext
    : public JuliaTensorOperationsGeneratorContext {
 public:
  using JuliaTensorOperationsGeneratorContext::
      JuliaTensorOperationsGeneratorContext;
};

/// Generator for producing Julia code using the ITensor framework
/// - https://julialang.org/
/// - https://jutho.github.io/TensorKit.jl
template <typename Context = JuliaTensorKitGeneratorContext>
class JuliaTensorKitGenerator : public JuliaTensorOperationsGenerator<Context> {
 private:
  using Base = JuliaTensorOperationsGenerator<Context>;

 public:
  JuliaTensorKitGenerator() = default;
  ~JuliaTensorKitGenerator() = default;

  std::string get_format_name() const override { return "Julia (TensorKit)"; }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "In Julia tensors can't be created without being initialized");
    }

    Base::m_generated += Base::tensor_name(tensor, ctx);
    Base::m_generated += " = TensorMap(zeros(Float64";

    for (const Index &idx : tensor.const_indices()) {
      std::string dim = ctx.get_dim(idx.space());

      Base::m_generated += ", ";
      Base::m_generated += dim;
    }

    Base::m_generated += "), ";
    Base::m_generated += domain(tensor, ctx);
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
    Base::m_generated += " = TensorMap(deserialize(\"";
    Base::m_generated += Base::tensor_name(tensor, ctx);
    Base::m_generated += ".jlbin\"), ";
    Base::m_generated += domain(tensor, ctx);
    Base::m_generated += ")\n";
  }

 private:
  std::string domain(const Tensor &tensor, const Context &ctx) const {
    std::string domain;

    const auto &indices = tensor.const_indices();
    const std::size_t braRank = tensor.bra_rank();

    if (braRank == 0 && indices.size() > 0) {
      throw std::runtime_error(
          "It is not (yet) clear how to represent a zero-dimensional domain "
          "for a tensor in TensorKit");
    }

    for (std::size_t i = 0; i < indices.size(); ++i) {
      const Index &idx = indices[i];

      domain += "ℝ^";
      domain += ctx.get_dim(idx.space());

      if (i == braRank - 1) {
        domain += ", ";
      } else if (i + 1 < indices.size()) {
        domain += " ⊗ ";
      }
    }

    return domain;
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP
