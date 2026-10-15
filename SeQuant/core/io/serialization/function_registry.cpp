#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/function_registry.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <mutex>

namespace sequant::io::serialization {

namespace {

// Builds a raw Operator<S> from an ordered list of action-tagged indices.
template <Statistics S>
ExprPtr build_operator(const FunctionCall &call) {
  container::svector<Op<S>> ops;
  ops.reserve(call.arguments.size());
  for (const FunctionArgument &arg : call.arguments) {
    if (!arg.is_action_index())
      throw SerializationError(
          "~op(...) expects action-tagged indices (e.g. i1+, a1-) as "
          "arguments");
    ops.emplace_back(arg.action_index->first, arg.action_index->second);
  }
  return ex<Operator<S>>(std::move(ops));
}

ExprPtr build_op(const FunctionCall &call, const DeserializationOptions &) {
  const Statistics stat = call.statistic.value_or(Statistics::FermiDirac);
  return stat == Statistics::FermiDirac
             ? build_operator<Statistics::FermiDirac>(call)
             : build_operator<Statistics::BoseEinstein>(call);
}

template <Statistics S>
ExprPtr build_sequence(const FunctionCall &call) {
  NormalOperatorSequence<S> seq;
  for (const FunctionArgument &arg : call.arguments) {
    if (arg.is_action_index() || !arg.expr ||
        !arg.expr->is<NormalOperator<S>>())
      throw SerializationError(
          "~nopseq(...) expects normal-ordered operators of a single statistic "
          "as arguments");
    seq.insert(seq.end(), arg.expr->as<NormalOperator<S>>());
  }
  return ex<NormalOperatorSequence<S>>(std::move(seq));
}

ExprPtr build_nopseq(const FunctionCall &call, const DeserializationOptions &) {
  // infer statistics from the explicit suffix, else from the first operator
  Statistics stat = Statistics::FermiDirac;
  if (call.statistic) {
    stat = *call.statistic;
  } else if (!call.arguments.empty() && call.arguments.front().expr) {
    stat = call.arguments.front()
                   .expr->is<NormalOperator<Statistics::BoseEinstein>>()
               ? Statistics::BoseEinstein
               : Statistics::FermiDirac;
  }
  return stat == Statistics::FermiDirac
             ? build_sequence<Statistics::FermiDirac>(call)
             : build_sequence<Statistics::BoseEinstein>(call);
}

}  // namespace

FunctionRegistry::FunctionRegistry() {
  register_function(L"op",
                    std::type_index(typeid(Operator<Statistics::FermiDirac>)),
                    build_op);
  register_function(L"op",
                    std::type_index(typeid(Operator<Statistics::BoseEinstein>)),
                    build_op);
  register_function(
      L"nopseq",
      std::type_index(typeid(NormalOperatorSequence<Statistics::FermiDirac>)),
      build_nopseq);
  register_function(
      L"nopseq",
      std::type_index(typeid(NormalOperatorSequence<Statistics::BoseEinstein>)),
      build_nopseq);
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
      return serializer(expr, options);
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

}  // namespace sequant::io::serialization
