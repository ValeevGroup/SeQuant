#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/function_registry.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <mutex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::io::serialization {

namespace {

// Builds a raw Operator<S> from an ordered list of action-tagged indices.
template <Statistics S>
ExprPtr build_operator(const FunctionCall &call,
                       const DeserializationOptions &options) {
  container::svector<Op<S>> ops;
  for (std::string_view arg : split_arguments(call.arguments)) {
    auto [index, action] = parse_action_index(arg, options);
    ops.emplace_back(std::move(index), action);
  }
  return ex<Operator<S>>(std::move(ops));
}

// Builds a NormalOperatorSequence<S> from a list of normal-ordered operators.
template <Statistics S>
ExprPtr build_sequence(const FunctionCall &call,
                       const DeserializationOptions &options) {
  NormalOperatorSequence<S> seq;
  for (std::string_view arg : split_arguments(call.arguments)) {
    ExprPtr nop = v2::from_string<ExprPtr>(arg, options);
    if (!nop->is<NormalOperator<S>>())
      throw SerializationError("'" + std::string(arg) +
                               "' is not a normal-ordered operator of the "
                               "sequence's statistics");
    seq.insert(seq.end(), nop->as<NormalOperator<S>>());
  }
  return ex<NormalOperatorSequence<S>>(std::move(seq));
}

template <Statistics S>
std::wstring write_operator(const Expr &expr,
                            const SerializationOptions &options) {
  return join_strings<std::wstring>(
      expr.as<Operator<S>>(), L", ", [&options](const Op<S> &op) {
        return v2::to_string(op.index(), options) +
               (op.action() == Action::Create ? L"+" : L"-");
      });
}

template <Statistics S>
std::wstring write_sequence(const Expr &expr,
                            const SerializationOptions &options) {
  return join_strings<std::wstring>(expr.as<NormalOperatorSequence<S>>(), L", ",
                                    [&options](const NormalOperator<S> &nop) {
                                      return v2::to_string(nop, options);
                                    });
}

constexpr std::string_view whitespace = " \t\r\n";

std::string_view trim(std::string_view text) {
  const auto begin = text.find_first_not_of(whitespace);
  if (begin == std::string_view::npos) return {};
  return text.substr(begin, text.find_last_not_of(whitespace) - begin + 1);
}

}  // namespace

FunctionRegistry::FunctionRegistry() {
  using enum Statistics;
  register_function(L"fop", std::type_index(typeid(Operator<FermiDirac>)),
                    build_operator<FermiDirac>, write_operator<FermiDirac>);
  register_function(L"bop", std::type_index(typeid(Operator<BoseEinstein>)),
                    build_operator<BoseEinstein>, write_operator<BoseEinstein>);
  register_function(L"fnopseq",
                    std::type_index(typeid(NormalOperatorSequence<FermiDirac>)),
                    build_sequence<FermiDirac>, write_sequence<FermiDirac>);
  register_function(
      L"bnopseq", std::type_index(typeid(NormalOperatorSequence<BoseEinstein>)),
      build_sequence<BoseEinstein>, write_sequence<BoseEinstein>);
}

FunctionRegistry &FunctionRegistry::instance() {
  static FunctionRegistry registry;
  return registry;
}

void FunctionRegistry::register_function(std::wstring name,
                                         std::type_index type,
                                         FunctionDeserializer deserializer,
                                         FunctionSerializer serializer) {
  static std::mutex mtx;
  std::lock_guard<std::mutex> lock(mtx);

  // find or create the entry for this name
  std::size_t idx = entries_.size();
  for (std::size_t i = 0; i < entries_.size(); ++i) {
    if (entries_[i].name == name) {
      idx = i;
      break;
    }
  }
  if (idx == entries_.size()) {
    entries_.push_back(
        {std::move(name), std::move(deserializer), std::move(serializer)});
  } else {
    entries_[idx].deserializer = std::move(deserializer);
    if (serializer) entries_[idx].serializer = std::move(serializer);
  }

  // (re)bind the Expr type to that entry
  for (auto &[t, i] : by_type_) {
    if (t == type) {
      i = idx;
      return;
    }
  }
  by_type_.emplace_back(type, idx);
}

ExprPtr FunctionRegistry::deserialize(
    const FunctionCall &call, const DeserializationOptions &options) const {
  for (const Entry &entry : entries_) {
    if (entry.name == call.name) return entry.deserializer(call, options);
  }
  throw SerializationError("Unknown function '" + toUtf8(call.name) +
                           "' in '~" + toUtf8(call.name) + "(...)'");
}

std::optional<std::wstring> FunctionRegistry::serialize(
    const Expr &expr, const SerializationOptions &options) const {
  const std::type_index type(typeid(expr));
  for (const auto &[t, idx] : by_type_) {
    if (t == type) {
      const FunctionSerializer &serializer = entries_[idx].serializer;
      if (!serializer) return std::nullopt;
      return L"~" + entries_[idx].name + L"(" + serializer(expr, options) +
             L")";
    }
  }
  return std::nullopt;
}

void register_function(std::wstring name, std::type_index type,
                       FunctionDeserializer deserializer,
                       FunctionSerializer serializer) {
  FunctionRegistry::instance().register_function(
      std::move(name), type, std::move(deserializer), std::move(serializer));
}

std::vector<std::string_view> split_arguments(std::string_view arguments) {
  std::vector<std::string_view> result;
  if (trim(arguments).empty()) return result;

  int depth = 0;
  std::size_t begin = 0;
  for (std::size_t i = 0; i < arguments.size(); ++i) {
    switch (arguments[i]) {
      case '(':
      case '[':
      case '{':
        ++depth;
        break;
      case ')':
      case ']':
      case '}':
        --depth;
        break;
      case ',':
        if (depth == 0) {
          result.push_back(trim(arguments.substr(begin, i - begin)));
          begin = i + 1;
        }
        break;
    }
  }
  result.push_back(trim(arguments.substr(begin)));
  return result;
}

std::pair<Index, Action> parse_action_index(
    std::string_view text, const DeserializationOptions &options) {
  text = trim(text);
  // the marks, as UTF-8: '+', '⁺' (U+207A), '-', '⁻' (U+207B)
  constexpr std::pair<std::string_view, Action> marks[] = {
      {"+", Action::Create},
      {"\xE2\x81\xBA", Action::Create},
      {"-", Action::Annihilate},
      {"\xE2\x81\xBB", Action::Annihilate}};
  for (const auto &[mark, action] : marks) {
    if (text.size() > mark.size() && text.ends_with(mark))
      return {v2::from_string<Index>(text.substr(0, text.size() - mark.size()),
                                     options),
              action};
  }
  throw SerializationError(
      "'" + std::string(text) +
      "' is not an index marked as a creator (+) or annihilator (-)");
}

}  // namespace sequant::io::serialization
