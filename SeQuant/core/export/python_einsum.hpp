#ifndef SEQUANT_CORE_EXPORT_PYTHON_EINSUM_HPP
#define SEQUANT_CORE_EXPORT_PYTHON_EINSUM_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace sequant {

/// Backend selection for Python einsum generator
enum class PythonEinsumBackend {
  NumPy,   ///< Generate code using numpy.einsum
  PyTorch  ///< Generate code using torch.einsum
};

/// Context class for the PythonEinsumGenerator
class PythonEinsumGeneratorContext : public ExportContext {
 public:
  using ShapeMap = std::map<IndexSpace, std::string>;

  PythonEinsumGeneratorContext(
      PythonEinsumBackend backend = PythonEinsumBackend::NumPy)
      : m_backend(backend) {}
  ~PythonEinsumGeneratorContext() = default;
  PythonEinsumGeneratorContext(
      ShapeMap index_shapes,
      PythonEinsumBackend backend = PythonEinsumBackend::NumPy)
      : m_index_shapes(std::move(index_shapes)), m_backend(backend) {}

  /// Get the dimension/shape for a given index space
  std::string get_shape(const IndexSpace &space) const {
    auto it = m_index_shapes.find(space);
    if (it != m_index_shapes.end()) {
      return it->second;
    }
    // Default fallback - use a simple naming scheme
    return "dim_" + toUtf8(space.base_key());
  }

  /// Get the shape tuple for a tensor
  std::string get_shape_tuple(const Tensor &tensor) const {
    std::string shape = "(";
    bool first = true;
    for (const Index &idx : tensor.const_indices()) {
      if (!first) shape += ", ";
      shape += get_shape(idx.space());
      first = false;
    }
    shape += ")";
    return shape;
  }

  /// Set the dimension/shape for a given index space
  void set_shape(const IndexSpace &space, std::string shape) {
    m_index_shapes[space] = std::move(shape);
  }

  /// Get the backend being used
  PythonEinsumBackend backend() const { return m_backend; }

  /// Set the backend
  void set_backend(PythonEinsumBackend backend) { m_backend = backend; }

  /// Get the module prefix (np. or torch.)
  std::string module_prefix() const {
    return m_backend == PythonEinsumBackend::NumPy ? "np." : "torch.";
  }

 protected:
  ShapeMap m_index_shapes;
  PythonEinsumBackend m_backend;
};

/// Generator for producing Python code using numpy.einsum or torch.einsum
/// - NumPy: https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
/// - PyTorch: https://pytorch.org/docs/stable/generated/torch.einsum.html
template <typename Context = PythonEinsumGeneratorContext>
class PythonEinsumGenerator : public Generator<Context> {
 public:
  PythonEinsumGenerator() = default;
  ~PythonEinsumGenerator() = default;

  std::string get_format_name() const override { return "Python (einsum)"; }

  bool supports_named_sections() const override { return true; }

  bool requires_named_sections() const override { return false; }

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

  std::string represent(const Index &idx, const Context &ctx) const override {
    if (idx.has_proto_indices()) {
      throw std::runtime_error("Proto Indices are not (yet) supported!");
    }

    // Convert index label to a single character for einsum notation
    // For einsum, we need single characters
    std::string label = toUtf8(idx.full_label());

    // Try to extract a single character representation
    // If the label is longer, use its first character (this may need
    // customization)
    if (label.empty()) {
      throw std::runtime_error("Empty index label");
    }

    // Return first character (lowercase for einsum convention)
    char c = label[0];
    if (std::isupper(static_cast<unsigned char>(c))) {
      c = std::tolower(static_cast<unsigned char>(c));
    }

    (void)ctx;  // Suppress unused parameter warning
    return std::string(1, c);
  }

  std::string represent(const Tensor &tensor,
                        [[maybe_unused]] const Context &ctx) const override {
    // For Python variable names, sanitize the tensor label
    return sanitize_python_name(tensor.label());
  }

  std::string represent(const Variable &variable,
                        const Context &) const override {
    return sanitize_python_name(variable.label());
  }

  std::string represent(const Constant &constant,
                        const Context &) const override {
    std::stringstream sstream;

    if (constant.value().imag() != 0) {
      sstream << "(";
      sstream << constant.value().real();
      if (constant.value().imag() < 0) {
        sstream << constant.value().imag() << "j";
      } else {
        sstream << "+" << constant.value().imag() << "j";
      }
      sstream << ")";
    } else {
      sstream << constant.value().real();
    }

    return sstream.str();
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Python tensors must be initialized when created");
    }

    // Use Fortran order to match Eigen::Tensor's column-major layout
    m_generated += m_indent + represent(tensor, ctx) + " = " +
                   ctx.module_prefix() + "zeros(" +
                   ctx.get_shape_tuple(tensor) + ", order='F')\n";
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
    m_generated += m_indent + represent(tensor, ctx) + " = ";

    if (set_to_zero) {
      // Use Fortran order to match Eigen::Tensor's column-major layout
      m_generated += ctx.module_prefix() + "zeros(" +
                     ctx.get_shape_tuple(tensor) + ", order='F')";
    } else {
      // Load from file (using pickle or similar)
      if (ctx.backend() == PythonEinsumBackend::NumPy) {
        m_generated +=
            ctx.module_prefix() + "load('" + represent(tensor, ctx) + ".npy')";
      } else {
        m_generated += "torch.load('" + represent(tensor, ctx) + ".pt')";
      }
    }

    m_generated += "\n";
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
    if (ctx.backend() == PythonEinsumBackend::NumPy) {
      m_generated += m_indent + represent(tensor, ctx) + ".fill(0)\n";
    } else {
      m_generated += m_indent + represent(tensor, ctx) + ".zero_()\n";
    }
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
    m_generated += m_indent + "del " + represent(tensor, ctx) + "\n";
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
    unload(tensor, ctx);
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
    if (ctx.backend() == PythonEinsumBackend::NumPy) {
      m_generated += m_indent + ctx.module_prefix() + "save('" +
                     represent(tensor, ctx) + ".npy', " +
                     represent(tensor, ctx) + ")\n";
    } else {
      m_generated += m_indent + "torch.save(" + represent(tensor, ctx) + ", '" +
                     represent(tensor, ctx) + ".pt')\n";
    }
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
      throw std::runtime_error(
          "Python variables must be initialized when created");
    }

    m_generated += m_indent + represent(variable, ctx) + " = 0.0\n";
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
    m_generated += m_indent + represent(variable, ctx) + " = ";

    if (set_to_zero) {
      m_generated += "0.0";
    } else {
      // Load from file
      if (ctx.backend() == PythonEinsumBackend::NumPy) {
        m_generated += ctx.module_prefix() + "load('" +
                       represent(variable, ctx) + ".npy')";
      } else {
        m_generated += "torch.load('" + represent(variable, ctx) + ".pt')";
      }
    }

    m_generated += "\n";
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
    m_generated += m_indent + represent(variable, ctx) + " = 0.0\n";
  }

  void unload(const Variable &variable, const Context &ctx) override {
    m_generated += m_indent + "del " + represent(variable, ctx) + "\n";
  }

  void destroy(const Variable &variable, const Context &ctx) override {
    unload(variable, ctx);
  }

  void persist(const Variable &variable, const Context &ctx) override {
    if (ctx.backend() == PythonEinsumBackend::NumPy) {
      m_generated += m_indent + ctx.module_prefix() + "save('" +
                     represent(variable, ctx) + ".npy', " +
                     represent(variable, ctx) + ")\n";
    } else {
      m_generated += m_indent + "torch.save(" + represent(variable, ctx) +
                     ", '" + represent(variable, ctx) + ".pt')\n";
    }
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
    m_generated += m_indent + represent(result, ctx) +
                   " += " + to_einsum_expr(expression, result, ctx) + "\n";
  }

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
    m_generated += m_indent + represent(result, ctx) +
                   " += " + to_python_scalar_expr(expression, ctx) + "\n";
  }

  void declare(const Index &idx, const Context &ctx) override {
    // Indices don't need explicit declaration in Python
    (void)idx;
    (void)ctx;
  }

  void declare(const Variable &variable, UsageSet usage,
               const Context &ctx) override {
    // Variables are dynamically typed in Python, no declaration needed
    (void)variable;
    (void)usage;
    (void)ctx;
  }

  void declare(const Tensor &tensor, UsageSet usage,
               const Context &ctx) override {
    // Tensors are dynamically typed in Python, no declaration needed
    (void)tensor;
    (void)usage;
    (void)ctx;
  }

  void all_indices_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void all_variables_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void all_tensors_declared(std::size_t amount, const Context &ctx) override {
    (void)amount;
    (void)ctx;
  }

  void begin_declarations(DeclarationScope scope, const Context &ctx) override {
    (void)scope;
    (void)ctx;
  }

  void end_declarations(DeclarationScope scope, const Context &ctx) override {
    (void)scope;
    (void)ctx;
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    (void)ctx;
    m_generated += m_indent + "# " + comment + "\n";
  }

  void begin_named_section(std::string_view name, const Context &ctx) override {
    (void)ctx;
    m_generated += m_indent + "def " + std::string(name) + "():\n";
    m_indent += "    ";
  }

  void end_named_section(std::string_view /*name*/,
                         const Context &ctx) override {
    (void)ctx;
    SEQUANT_ASSERT(m_indent.size() >= 4);
    m_indent = m_indent.substr(0, m_indent.size() - 4);
    m_generated += "\n";
  }

  void begin_expression(const Context &ctx) override {
    (void)ctx;
    if (!m_generated.empty() && !m_generated.ends_with("\n\n")) {
      m_generated += "\n";
    }
  }

  void end_expression(const Context &ctx) override { (void)ctx; }

  void begin_export(const Context &ctx) override {
    (void)ctx;
    m_generated.clear();
  }

  void end_export(const Context &ctx) override { (void)ctx; }

  std::string get_generated_code() const override { return m_generated; }

 protected:
  std::string m_generated;
  std::string m_indent;

  /// Sanitize a label to be a valid Python identifier
  std::string sanitize_python_name(std::wstring_view label) const {
    std::string name = toUtf8(label);

    // Replace invalid characters with underscores
    for (char &c : name) {
      if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') {
        c = '_';
      }
    }

    // Ensure it doesn't start with a number
    if (!name.empty() && std::isdigit(static_cast<unsigned char>(name[0]))) {
      name = "_" + name;
    }

    return name;
  }

  /// Get or create a single-character einsum index for a given Index
  std::string get_einsum_index(
      const Index &idx, [[maybe_unused]] const Context &ctx,
      boost::unordered::unordered_map<std::string, std::string> &index_map,
      boost::unordered::unordered_set<std::string> &used_chars) const {
    std::string full_label = toUtf8(idx.full_label());

    // Check if we've already assigned a character to this index
    auto it = index_map.find(full_label);
    if (it != index_map.end()) {
      return it->second;
    }

    // Try to use the first character of the label
    char base_char =
        full_label.empty()
            ? 'i'
            : std::tolower(static_cast<unsigned char>(full_label[0]));
    std::string candidate(1, base_char);

    // Check if this character is already used (O(1) lookup)
    bool found = used_chars.find(candidate) != used_chars.end();

    // If already used, find next available character
    if (found) {
      // Use a sequence: a-z, then A-Z
      constexpr char chars[] =
          "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
      constexpr auto num_chars = sizeof(chars) / sizeof(chars[0]) - 1;

      for (int i = 0; i < num_chars; ++i) {
        candidate = std::string(1, chars[i]);
        // O(1) lookup instead of O(n) iteration
        if (used_chars.find(candidate) == used_chars.end()) {
          found = false;
          break;
        }
      }

      if (found) {
        // Fallback: use numbered indices (not standard einsum, but necessary)
        throw std::runtime_error("Too many unique indices for einsum notation");
      }
    }

    // Store in both maps for fast forward and reverse lookups
    index_map[full_label] = candidate;
    used_chars.insert(candidate);
    return candidate;
  }

  /// Convert a tensor to its einsum subscript notation
  std::string tensor_to_einsum_subscript(
      const Tensor &tensor, const Context &ctx,
      boost::unordered::unordered_map<std::string, std::string> &index_map,
      boost::unordered::unordered_set<std::string> &used_chars) const {
    std::string subscript;
    for (const Index &idx : tensor.const_indices()) {
      subscript += get_einsum_index(idx, ctx, index_map, used_chars);
    }
    return subscript;
  }

  /// Convert an expression to an einsum call
  std::string to_einsum_expr(const Expr &expr, const Tensor &result,
                             const Context &ctx) const {
    // Create local index mapping for this expression
    boost::unordered::unordered_map<std::string, std::string> index_map;
    boost::unordered::unordered_set<std::string> used_chars;

    std::string einsum_spec;
    std::vector<std::string> tensor_names;
    std::string scalar_factor;

    // Extract scalar prefactor and tensors from the expression
    extract_einsum_components(expr, einsum_spec, tensor_names, scalar_factor,
                              ctx, index_map, used_chars);

    // Add output specification
    einsum_spec +=
        "->" + tensor_to_einsum_subscript(result, ctx, index_map, used_chars);

    // Build einsum call
    std::string einsum_call =
        ctx.module_prefix() + "einsum('" + einsum_spec + "'";

    for (const std::string &name : tensor_names) {
      einsum_call += ", " + name;
    }

    if (ctx.backend() == PythonEinsumBackend::NumPy) {
      einsum_call += ", optimize=True";
    }

    einsum_call += ")";

    // Apply scalar factor if present
    if (!scalar_factor.empty() && scalar_factor != "1" &&
        scalar_factor != "1.0") {
      einsum_call = scalar_factor + " * " + einsum_call;
    }

    return einsum_call;
  }

  /// Extract einsum components from an expression
  void extract_einsum_components(
      const Expr &expr, std::string &einsum_spec,
      std::vector<std::string> &tensor_names, std::string &scalar_factor,
      const Context &ctx,
      boost::unordered::unordered_map<std::string, std::string> &index_map,
      boost::unordered::unordered_set<std::string> &used_chars) const {
    if (expr.is<Tensor>()) {
      const Tensor &tensor = expr.as<Tensor>();
      if (!einsum_spec.empty()) einsum_spec += ",";
      einsum_spec +=
          tensor_to_einsum_subscript(tensor, ctx, index_map, used_chars);
      tensor_names.push_back(represent(tensor, ctx));
    } else if (expr.is<Variable>()) {
      // Treat variable as a scalar
      if (scalar_factor.empty()) {
        scalar_factor = represent(expr.as<Variable>(), ctx);
      } else {
        scalar_factor += " * " + represent(expr.as<Variable>(), ctx);
      }
    } else if (expr.is<Constant>()) {
      std::string const_repr = represent(expr.as<Constant>(), ctx);
      if (scalar_factor.empty()) {
        scalar_factor = const_repr;
      } else {
        scalar_factor += " * " + const_repr;
      }
    } else if (expr.is<Product>()) {
      const Product &product = expr.as<Product>();

      // Handle scalar part
      if (!product.scalar().is_identity()) {
        std::string scalar_repr = represent(Constant(product.scalar()), ctx);
        if (scalar_factor.empty()) {
          scalar_factor = scalar_repr;
        } else {
          scalar_factor += " * " + scalar_repr;
        }
      }

      // Handle tensor factors
      for (std::size_t i = 0; i < product.size(); ++i) {
        extract_einsum_components(*product.factor(i), einsum_spec, tensor_names,
                                  scalar_factor, ctx, index_map, used_chars);
      }
    } else if (expr.is<Sum>()) {
      // For sums, we can't use a single einsum call
      // This should be handled at a higher level by generating multiple
      // statements
      throw std::runtime_error(
          "Sum expressions should be handled by generating multiple einsum "
          "calls");
    }
  }

  /// Convert expression to Python scalar expression (for variable results)
  std::string to_python_scalar_expr(const Expr &expr,
                                    const Context &ctx) const {
    if (expr.is<Variable>()) {
      return represent(expr.as<Variable>(), ctx);
    } else if (expr.is<Constant>()) {
      return represent(expr.as<Constant>(), ctx);
    } else if (expr.is<Product>()) {
      // This should involve tensor contractions
      // For a scalar result, we need to contract all indices
      // This is an einsum with no output indices

      // Create local index mapping for this expression
      boost::unordered::unordered_map<std::string, std::string> index_map;
      boost::unordered::unordered_set<std::string> used_chars;

      std::string einsum_spec;
      std::vector<std::string> tensor_names;
      std::string scalar_factor;

      extract_einsum_components(expr, einsum_spec, tensor_names, scalar_factor,
                                ctx, index_map, used_chars);

      // For scalar result, use "->" or "->..."
      einsum_spec += "->";

      std::string einsum_call =
          ctx.module_prefix() + "einsum('" + einsum_spec + "'";

      for (const std::string &name : tensor_names) {
        einsum_call += ", " + name;
      }

      if (ctx.backend() == PythonEinsumBackend::NumPy) {
        einsum_call += ", optimize=True";
      }

      einsum_call += ")";

      if (!scalar_factor.empty() && scalar_factor != "1" &&
          scalar_factor != "1.0") {
        einsum_call = scalar_factor + " * " + einsum_call;
      }

      return einsum_call;
    } else if (expr.is<Sum>()) {
      const Sum &sum = expr.as<Sum>();
      std::string result = "(";

      for (std::size_t i = 0; i < sum.size(); ++i) {
        if (i > 0) result += " + ";
        result += to_python_scalar_expr(*sum.summand(i), ctx);
      }

      result += ")";
      return result;
    }

    throw std::runtime_error(
        "Unsupported expression type for Python scalar expression");
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_PYTHON_EINSUM_HPP
