#ifndef SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP
#define SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP

#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <range/v3/view/enumerate.hpp>

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
    Base::create_wrapped(tensor, zero_init, ctx, "TensorMap",
                         domain(tensor, ctx));
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    Base::load_wrapped(tensor, set_to_zero, ctx, "TensorMap",
                       domain(tensor, ctx));
  }

 private:
  std::string domain(const Tensor &tensor, const Context &ctx) const {
    std::string domain;

    const std::size_t braRank = tensor.bra_rank();
    const std::size_t num_indices = tensor.num_indices();

    if (braRank == 0 && num_indices > 0) {
      throw Exception(
          "It is not (yet) clear how to represent a zero-dimensional domain "
          "for a tensor in TensorKit");
    }

    for (const auto &[i, idx] : ranges::views::enumerate(tensor.indices())) {
      domain += "ℝ^";
      domain += ctx.get_dim(idx.space());

      if (i == braRank - 1) {
        domain += ", ";
      } else if (i + 1 < num_indices) {
        domain += " ⊗ ";
      }
    }

    return domain;
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIA_TENSOR_KIT_HPP
