#ifndef SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP
#define SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/reordering_context.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

#include <cstdlib>
#include <ranges>
#include <string>
#include <string_view>

namespace sequant {

/// Context class for the JuliaTensorOperationsGenerator
class JuliaTensorOperationsGeneratorContext : public ReorderingContext {
 public:
  using TagMap = std::map<IndexSpace, std::string>;
  using DimMap = std::map<IndexSpace, std::string>;

  JuliaTensorOperationsGeneratorContext()
      : ReorderingContext(MemoryLayout::ColumnMajor){};
  ~JuliaTensorOperationsGeneratorContext() = default;
  JuliaTensorOperationsGeneratorContext(TagMap index_tags, DimMap index_dims);

  std::string get_tag(const IndexSpace &space) const;

  std::string get_dim(const IndexSpace &space) const;

  std::string get_tags(const Tensor &tensor) const;

  void set_tag(const IndexSpace &space, std::string tag);
  void set_dim(const IndexSpace &space, std::string dim);

 protected:
  TagMap m_index_tags;
  DimMap m_index_dims;
};

/// Generator for producing Julia code using the ITensor framework
/// - https://julialang.org/
/// - https://jutho.github.io/TensorOperations.jl
template <typename Context = JuliaTensorOperationsGeneratorContext>
class JuliaTensorOperationsGenerator : public Generator<Context> {
 public:
  JuliaTensorOperationsGenerator() = default;
  ~JuliaTensorOperationsGenerator() = default;

  std::string get_format_name() const override {
    return "Julia (TensorOperations)";
  }

  bool supports_named_sections() const override { return false; }

  bool requires_named_sections() const override { return false; }

  bool supports_index_batching() const override { return false; }

  DeclarationScope index_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope variable_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope tensor_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  PrunableScalars prunable_scalars() const override {
    return PrunableScalars::All;
  }

  std::string represent(const Index &idx, const Context &) const override {
    if (idx.has_proto_indices()) {
      throw Exception("Proto Indices are not (yet) supported!");
    }

    return toUtf8(idx.label()) + instance_suffix(idx);
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = tensor_name(tensor, ctx);

    representation += "[ ";

    using namespace std::literals;
    representation += tensor.const_indices() |
                      ranges::views::transform([&](const Index &idx) {
                        return represent(idx, ctx);
                      }) |
                      ranges::views::join(", "sv) | ranges::to<std::string>();

    representation += " ]";

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
        sstream << "-" << (-1 * constant.value().imag()) << "im)";
      } else {
        sstream << "+" << constant.value().imag() << "im)";
      }
    } else {
      sstream << constant.value().real();
    }

    return sstream.str();
  }

  std::string represent(const Power &power, const Context &ctx) const override {
    return detail::format_power(power, to_julia_expr(*power.base(), ctx), "^",
                                /*double_slash*/ true, &wrap_conj);
  }

  std::string represent(const TensorPower &power,
                        const Context &ctx) const override {
    SEQUANT_ASSERT(power.base()->is<Tensor>());
    return detail::format_tensor_power(
        power, tensor_name(power.base()->as<Tensor>(), ctx), ".^",
        /*double_slash*/ true);
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw Exception(
          "In Julia tensors can't be created without being initialized");
    }

    m_generated += zero_initialization(tensor, ctx) + "\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    if (!set_to_zero) {
      m_generated += tensor_name(tensor, ctx);
      m_generated += " = ";
      m_generated += "deserialize(\"" + tensor_name(tensor, ctx) + ".jlbin\")";
    } else {
      // zero_initialization already emits "<name> = zeros(...)"
      m_generated += zero_initialization(tensor, ctx);
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    // In Julia setting to zero essentially has same effect as creating a zero
    // tensor
    create(tensor, true, ctx);
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += tensor_name(tensor, ctx) + " = nothing\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    unload(tensor, ctx);
    // Remove potential serialized version of this tensor from disk
    m_generated += "rm( \"" + tensor_name(tensor, ctx) + +".jlbin\"" + ")\n";
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    m_generated += "return " + tensor_name(tensor, ctx) + "\n";
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw Exception(
          "Julia doesn't support declaring a variable without initializing it");
    }

    m_generated += represent(variable, ctx) + " = 0.0\n";
  }
  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = ";

    if (set_to_zero)
      m_generated += "0.0";
    else {
      m_generated += "deserialize(\"" + represent(variable, ctx) + ".jlbin\")";
    }
    m_generated += "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = 0.0\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += represent(variable, ctx) + " = nothing\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    unload(variable, ctx);
    // Remove any potentially serialized version from disk
    m_generated += "rm(\"" + represent(variable, ctx) + ".jlbin\")\n";
  }

  void persist(const Variable &variable, const Context &ctx) override {
    m_generated += "return " + represent(variable, ctx) + "\n";
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    if (const TensorPower *power = find_tensor_power(expression)) {
      m_generated += tensor_name(result, ctx) + " .+= " +
                     broadcast_power(expression, *power, result, ctx) + "\n";
      return;
    }

    m_generated += "@tensor " + represent(result, ctx) +
                   " += " + to_julia_expr(expression, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += "@tensor " + represent(result, ctx) +
                   " += " + to_julia_expr(expression, ctx) + "\n";
  }

  void declare(const Index &, const Context &) override {}

  void declare(const Variable &, UsageSet, const Context &) override {}

  void declare(const Tensor &, UsageSet, const Context &) override {}

  void all_indices_declared(std::size_t, const Context &) override {}

  void all_variables_declared(std::size_t, const Context &) override {}

  void all_tensors_declared(std::size_t, const Context &) override {}

  void begin_declarations(DeclarationScope, const Context &) override {}

  void end_declarations(DeclarationScope, const Context &) override {}

  void insert_comment(const std::string &comment, const Context &) override {
    m_generated += "# " + comment + "\n";
  }

  void begin_named_section(std::string_view, const Context &) override {
    SEQUANT_ABORT("This function should not have been called");
  }

  void end_named_section(std::string_view, const Context &) override {
    SEQUANT_ABORT("This function should not have been called");
  }

  void begin_expression(const Context &) override {
    if (!m_generated.empty() && !m_generated.ends_with("\n\n")) {
      m_generated += "\n";
    }
  }

  void end_expression(const Context &) override {}

  void begin_export(const Context &) override { m_generated.clear(); }

  void end_export(const Context &) override {}

  std::string get_generated_code() const override { return m_generated; }

 protected:
  std::string m_generated;

  std::string to_julia_expr(const Expr &expr, const Context &ctx) const {
    return detail::stringify_expr(
        expr, " * ",
        [this, &ctx](const auto &e) {
          // through the base so that no overload is hidden in a derived class
          return static_cast<const Generator<Context> &>(*this).represent(e,
                                                                          ctx);
        },
        "to_julia_expr");
  }

  static std::string wrap_conj(std::string s) {
    return "conj(" + std::move(s) + ")";
  }

  /// @return the element-wise power in @p expression if it is a TensorPower or
  /// a Product with a TensorPower factor, else nullptr
  static const TensorPower *find_tensor_power(const Expr &expression) {
    if (expression.is<TensorPower>()) return &expression.as<TensorPower>();
    if (expression.is<Product>()) {
      for (const ExprPtr &factor : expression.as<Product>().factors()) {
        if (factor->is<TensorPower>()) return &factor->as<TensorPower>();
      }
    }
    return nullptr;
  }

  /// @return the broadcast form of @p expression, the element-wise @p power
  /// (of an operand laid out as @p result) times the scalars of @p expression,
  /// as `[<scalars> * ]<operand> .^ <exponent>`, the operand permuted to the
  /// index order of @p result if that differs
  std::string broadcast_power(const Expr &expression, const TensorPower &power,
                              const Tensor &result, const Context &ctx) const {
    std::string scalars = scalar_prefix(expression, power, ctx);
    SEQUANT_ASSERT(power.base()->is<Tensor>());
    const Tensor &operand = power.base()->as<Tensor>();
    return scalars + detail::format_tensor_power(
                         power,
                         permuted(tensor_name(operand, ctx), operand,
                                  result.const_indices()),
                         " .^ ", /*double_slash*/ true);
  }

  /// @return the scalar factors of @p expression other than @p nonscalar, as
  /// `<scalars> * `, or empty if there are none
  /// @throw Exception if @p expression has a non-scalar factor other than
  ///        @p nonscalar
  std::string scalar_prefix(const Expr &expression, const Expr &nonscalar,
                            const Context &ctx) const {
    std::string scalars;
    if (expression.is<Product>()) {
      const Product &product = expression.as<Product>();
      if (!product.scalar().is_identity()) {
        scalars += represent(Constant(product.scalar()), ctx) + " * ";
      }
      for (const ExprPtr &factor : product.factors()) {
        if (factor.get() == &nonscalar) continue;
        if (!factor->is_scalar()) {
          throw Exception(
              "Julia TensorOperations: an element-wise power cannot be part of "
              "a contraction");
        }
        scalars += to_julia_expr(*factor, ctx) + " * ";
      }
    }
    return scalars;
  }

  /// @return @p operand_str (rendering @p operand) permuted by `permutedims`
  /// so that its modes are @p target in that order, or as is if they already
  /// are
  template <std::ranges::input_range Indices>
  static std::string permuted(std::string operand_str, const Tensor &operand,
                              Indices &&target) {
    // permutedims(A, p) has A's dimension p[k] as its k-th dimension
    std::string perm;
    bool is_identity = true;
    std::size_t k = 0;
    for (const Index &idx : target) {
      ++k;
      std::size_t pos = 0;
      [[maybe_unused]] bool found = false;
      for (const Index &operand_idx : operand.const_indices()) {
        ++pos;
        if (operand_idx == idx) {
          found = true;
          break;
        }
      }
      SEQUANT_ASSERT(found);
      is_identity = is_identity && pos == k;
      if (!perm.empty()) perm += ", ";
      perm += std::to_string(pos);
    }
    SEQUANT_ASSERT(k == operand.num_indices());
    if (is_identity) return operand_str;
    return "permutedims(" + std::move(operand_str) + ", (" + perm + "))";
  }

  /// emits `<name> = <wrapper>(zeros(Float64, <dims>...), <shape>)`
  void create_wrapped(const Tensor &tensor, bool zero_init, const Context &ctx,
                      std::string_view wrapper, const std::string &shape) {
    if (!zero_init) {
      throw Exception(
          "In Julia tensors can't be created without being initialized");
    }

    m_generated += tensor_name(tensor, ctx);
    m_generated += " = ";
    m_generated += wrapper;
    m_generated += "(zeros(Float64";

    for (const Index &idx : tensor.const_indices()) {
      m_generated += ", ";
      m_generated += dim(idx, ctx);
    }

    m_generated += "), ";
    m_generated += shape;
    m_generated += ")\n";
  }

  /// emits `<name> = <wrapper>(deserialize("<name>.jlbin"), <shape>)`, or
  /// the zero-initialization of create_wrapped if @p set_to_zero
  void load_wrapped(const Tensor &tensor, bool set_to_zero, const Context &ctx,
                    std::string_view wrapper, const std::string &shape) {
    if (set_to_zero) {
      // In Julia setting a tensor to zero has to be done by overwriting it with
      // a zero tensor
      create_wrapped(tensor, true, ctx, wrapper, shape);
      return;
    }

    m_generated += tensor_name(tensor, ctx);
    m_generated += " = ";
    m_generated += wrapper;
    m_generated += "(deserialize(\"";
    m_generated += tensor_name(tensor, ctx);
    m_generated += ".jlbin\"), ";
    m_generated += shape;
    m_generated += ")\n";
  }

  /// @return `_N` for an index in basis instance `N`, else empty
  static std::string instance_suffix(const Index &idx) {
    const auto tag = detail::basis_instance_tag(idx);
    return tag.empty() ? tag : "_" + tag;
  }

  /// @return the dimension variable of @p idx: the space's, suffixed by the
  /// basis instance (if any) since each basis has its own extent
  static std::string dim(const Index &idx, const Context &ctx) {
    return detail::dim_name(ctx.get_dim(idx.space()), idx);
  }

  std::string tensor_name(const Tensor &tensor, const Context &ctx) const {
    std::string representation = toUtf8(tensor.label());

    representation += "_";
    representation += ctx.get_tags(tensor);

    return representation;
  }

  std::string zero_initialization(const Tensor &tensor,
                                  const Context &ctx) const {
    std::string code = tensor_name(tensor, ctx);
    code += " = zeros(Float64";

    for (const Index &idx : tensor.const_indices()) {
      code += ", ";
      code += dim(idx, ctx);
    }

    code += ")";

    return code;
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_JULIA_TENSOR_OPERATIONS_HPP
