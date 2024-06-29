#ifndef SEQUANT_CORE_EXPORT_JULIATENSORKITGEN_HPP
#define SEQUANT_CORE_EXPORT_JULIATENSORKITGEN_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/julia_tensoroperationsgen.hpp>
namespace sequant {

class JuliaTensorKitGenContext : public JuliaTensorOperationsGenContext {
 public:
  using JuliaTensorOperationsGenContext::JuliaTensorOperationsGenContext;
};

template <typename Context = JuliaTensorKitGenContext>
class JuliaTensorKitGen : public JuliaTensorOperationsGen<Context> {
 public:
  JuliaTensorKitGen() = default;
  ~JuliaTensorKitGen() = default;

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
    std::string rhs = "zeros(Float64";
    for (std::size_t i = 0; i < indices.size(); ++i) {
      rhs += ", " + dimStrings[i];
    }
    rhs += ")";
    std::string rhs2 = "TensorMap(" + rhs + ", ";
    for (std::size_t i = 0; i < indices.size() / 2; ++i) {
      if (i == (indices.size() / 2) - 1) {
        rhs2 += "ℝ^" + dimStrings[i] + ",";
      } else {
        rhs2 += "ℝ^" + dimStrings[i] + " ⊗ ";
      }
    }
    for (std::size_t i = indices.size() / 2; i < indices.size(); ++i) {
      if (i == indices.size() - 1) {
        rhs2 += "ℝ^" + dimStrings[i] + ")";
      } else {
        rhs2 += "ℝ^" + dimStrings[i] + " ⊗ ";
      }
    }
    this->m_generated += " = " + rhs2 + "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (!set_to_zero) {
      this->m_generated += this->represent_tensor_noindex(tensor, ctx) +
                           " = TensorMap(deserialize(\"" +
                           this->represent_tensor_noindex(tensor, ctx) +
                           ".jlbin\"),";
      const auto &indices = tensor.const_braket();
      std::vector<std::string> dimStrings(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        dimStrings[i] = ctx.get_dim(indices[i].space());
      }
      for (std::size_t i = 0; i < indices.size() / 2; ++i) {
        if (i == (indices.size() / 2) - 1) {
          this->m_generated += "ℝ^" + dimStrings[i] + ",";
        } else {
          this->m_generated += "ℝ^" + dimStrings[i] + " ⊗ ";
        }
      }
      for (std::size_t i = indices.size() / 2; i < indices.size(); ++i) {
        if (i == indices.size() - 1) {
          this->m_generated += "ℝ^" + dimStrings[i] + ")";
        } else {
          this->m_generated += "ℝ^" + dimStrings[i] + " ⊗ ";
        }
      }
    }

    else {
      this->m_generated += this->represent_tensor_noindex(tensor, ctx);
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
      std::string rhs2 = "TensorMap(" + rhs + ", ";
      for (std::size_t i = 0; i < indices.size() / 2; ++i) {
        if (i == (indices.size() / 2) - 1) {
          rhs2 += "ℝ^" + dimStrings[i] + ",";
        } else {
          rhs2 += "ℝ^" + dimStrings[i] + " ⊗ ";
        }
      }
      for (std::size_t i = indices.size() / 2; i < indices.size(); ++i) {
        if (i == indices.size() - 1) {
          rhs2 += "ℝ^" + dimStrings[i] + ")";
        } else {
          rhs2 += "ℝ^" + dimStrings[i] + " ⊗ ";
        }
      }
      this->m_generated += " = " + rhs2 + "\n";
    }

    this->m_generated += "\n";
  }
  std::string get_format_name() const override { return "Julia TensorKit"; }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIATENSORKITGEN_HPP
