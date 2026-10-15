#ifndef SEQUANT_CORE_IO_SERIALIZATION_FUNCTION_REGISTRY_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_FUNCTION_REGISTRY_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <typeindex>
#include <utility>
#include <vector>

namespace sequant::io::serialization {

/// A v2 `~name(...)` function call handed to a registered deserializer.
struct FunctionCall {
  std::wstring name;
  /// the text between the parentheses as written (UTF-8); the deserializer
  /// parses it, e.g. with split_arguments(), v2::from_string() and
  /// parse_action_index()
  std::string arguments;
};

/// builds an Expr from a parsed function call
using FunctionDeserializer = std::function<ExprPtr(
    const FunctionCall &, const DeserializationOptions &)>;
/// emits the arguments of the `~name(...)` representation of an Expr, i.e. the
/// text that FunctionCall::arguments holds when deserializing it; may be empty
/// (nullptr) if the function supports deserialization only
using FunctionSerializer =
    std::function<std::wstring(const Expr &, const SerializationOptions &)>;

/// Registry backing the v2 `~name(...)` extension syntax. It maps a function
/// name to a deserialization builder and an Expr subtype to a serialization
/// emitter. Core registers the built-in `fop`/`bop` (Operator<S>) and
/// `fnopseq`/`bnopseq` (NormalOperatorSequence<S>) handlers; downstream code
/// registers handlers for its own Expr subclasses through register_function().
class FunctionRegistry {
 public:
  static FunctionRegistry &instance();

  /// Registers a handler for function \p name that builds/consumes objects of
  /// the Expr subtype \p type. A null \p serializer registers deserialization
  /// support only. Re-registering a name or type replaces the previous handler.
  void register_function(std::wstring name, std::type_index type,
                         FunctionDeserializer deserializer,
                         FunctionSerializer serializer = nullptr);

  /// @return the Expr built from \p call
  /// @throw SerializationError if `call.name` is not registered
  ExprPtr deserialize(const FunctionCall &call,
                      const DeserializationOptions &options) const;

  /// @return the `~name(...)` serialization of \p expr, or std::nullopt if its
  /// dynamic type has no registered serializer
  std::optional<std::wstring> serialize(
      const Expr &expr, const SerializationOptions &options) const;

 private:
  FunctionRegistry();

  struct Entry {
    std::wstring name;
    FunctionDeserializer deserializer;
    FunctionSerializer serializer;
  };

  // name -> entry (deserialization lookup)
  std::vector<Entry> entries_;
  // Expr type -> index into entries_ (serialization lookup)
  std::vector<std::pair<std::type_index, std::size_t>> by_type_;
};

/// Convenience wrapper for FunctionRegistry::instance().register_function().
void register_function(std::wstring name, std::type_index type,
                       FunctionDeserializer deserializer,
                       FunctionSerializer serializer = nullptr);

/// @return the comma-separated arguments in \p arguments (see
/// FunctionCall::arguments) without surrounding whitespace; a comma inside
/// brackets or a quoted label does not separate arguments. Blank \p arguments
/// have none.
std::vector<std::string_view> split_arguments(std::string_view arguments);

/// @return the index and action of an action-tagged index such as `i1+` or
/// `a1⁻`: a creator is marked by `+`/`⁺`, an annihilator by `-`/`⁻`
/// @throw SerializationError if \p text is not an action-tagged index
std::pair<Index, Action> parse_action_index(
    std::string_view text, const DeserializationOptions &options = {});

}  // namespace sequant::io::serialization

#endif  // SEQUANT_CORE_IO_SERIALIZATION_FUNCTION_REGISTRY_HPP
