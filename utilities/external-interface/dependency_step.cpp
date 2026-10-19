#include "dependency_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::util::extint {

namespace {

/// Appends @p id to @p entry's dependencies unless it is already present.
void add_dependency(ExportTreeData::Entry &entry, std::size_t id) {
  if (std::ranges::find(entry.dependencies, id) == entry.dependencies.end()) {
    entry.dependencies.push_back(id);
  }
}

/// Invokes @p tensor_fn / @p variable_fn for the tensor/variable represented by
/// every leaf of @p tree. Mirrors the leaf dispatch used when collecting tensor
/// usage during export (see track_usage in SeQuant/core/export/export.hpp).
template <typename TensorFn, typename VariableFn>
void for_each_leaf(const ExportNode<> &tree, TensorFn &&tensor_fn,
                   VariableFn &&variable_fn) {
  tree.visit_leaf([&](const auto &node) {
    const Expr &expr = *node->expr();

    if (expr.is<Tensor>()) {
      tensor_fn(expr.as<Tensor>());
    } else if (expr.is<Variable>()) {
      variable_fn(expr.as<Variable>());
    } else if (expr.is<Power>()) {
      const Power &power = expr.as<Power>();
      if (power.base().is<Tensor>()) {
        tensor_fn(power.base().as<Tensor>());
      } else if (power.base().is<Variable>()) {
        variable_fn(power.base().as<Variable>());
      }
    }
  });
}

}  // namespace

SEQUANT_EXTINT_REGISTER_STEP_TYPE(DependencyStep);

std::string DependencyStep::kind() const { return "dependencies"; }

bool DependencyStep::accepts_options() const { return true; }

bool DependencyStep::requires_options() const { return true; }

void DependencyStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  bool has_autodetect = false;
  bool has_deps = false;

  for (const auto &[key, value] : options.items()) {
    if (key == "autodetect") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      autodetect_ = value.get<bool>();
      has_autodetect = true;
    } else if (key == "deps") {
      if (!value.is_object()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be an object");
      }

      for (const auto &[consumer, producers] : value.items()) {
        std::vector<std::string> producer_ids;

        if (producers.is_string()) {
          producer_ids.emplace_back(producers.get<std::string>());
        } else if (producers.is_array()) {
          for (const nlohmann::json &current : producers) {
            if (!current.is_string()) {
              throw Exception("Producer IDs for '" + consumer + "' in " +
                              kind() + "." + key + " must be strings");
            }

            producer_ids.emplace_back(current.get<std::string>());
          }
        } else {
          throw Exception("Value for '" + consumer + "' in " + kind() + "." +
                          key + " must be a string or an array of strings");
        }

        deps_.emplace_back(consumer, std::move(producer_ids));
      }

      has_deps = true;
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  if (!has_autodetect && !has_deps) {
    throw Exception(kind() +
                    " requires at least one of 'autodetect' or 'deps' to be "
                    "set");
  }
}

std::size_t DependencyStep::run(std::string_view step_id,
                                ExecutionContext &exctx,
                                const std::vector<std::string_view> &inputs) {
  // Gather all input data. The reference_wrapper inside each Data points into
  // the context's store, so annotating the referenced ExportTreeData persists.
  std::vector<ExecutionContext::Data<ProcessingData>> data;
  for (std::string_view current_input : inputs) {
    for (ExecutionContext::Data<ProcessingData> current :
         exctx.get_data(current_input)) {
      data.push_back(std::move(current));
    }
  }

  if (autodetect_) {
    // Index all producible results (block-level for tensors) to their tree ids.
    std::map<Tensor, std::vector<std::size_t>, TensorBlockLessThanComparator>
        tensor_producers;
    std::map<Variable, std::vector<std::size_t>> variable_producers;

    for (ExecutionContext::Data<ProcessingData> &current : data) {
      ExportTreeData &tree_data =
          convert_data<ExportTreeData>(current.data.get());

      for (const ExportTreeData::Entry &entry : tree_data.entries) {
        if (entry.tree->is_tensor()) {
          tensor_producers[entry.tree->as_tensor()].push_back(entry.tree->id());
        } else {
          SEQUANT_ASSERT(entry.tree->is_variable());
          variable_producers[entry.tree->as_variable()].push_back(
              entry.tree->id());
        }
      }
    }

    for (ExecutionContext::Data<ProcessingData> &current : data) {
      ExportTreeData &tree_data =
          convert_data<ExportTreeData>(current.data.get());

      for (ExportTreeData::Entry &entry : tree_data.entries) {
        const std::size_t self_id = entry.tree->id();

        for_each_leaf(
            entry.tree,
            [&](const Tensor &tensor) {
              auto it = tensor_producers.find(tensor);
              if (it == tensor_producers.end()) {
                return;
              }

              for (std::size_t producer_id : it->second) {
                if (producer_id != self_id) {
                  add_dependency(entry, producer_id);
                }
              }
            },
            [&](const Variable &variable) {
              auto it = variable_producers.find(variable);
              if (it == variable_producers.end()) {
                return;
              }

              for (std::size_t producer_id : it->second) {
                if (producer_id != self_id) {
                  add_dependency(entry, producer_id);
                }
              }
            });
      }
    }
  }

  for (const auto &[consumer_id, producer_ids] : deps_) {
    std::vector<std::size_t> producer_root_ids;
    for (const std::string &producer_id : producer_ids) {
      for (const ExecutionContext::Data<const ProcessingData> &producer :
           std::as_const(exctx).get_data(producer_id)) {
        const ExportTreeData &tree_data =
            convert_data<ExportTreeData>(producer.data.get());

        for (const ExportTreeData::Entry &entry : tree_data.entries) {
          producer_root_ids.push_back(entry.tree->id());
        }
      }
    }

    for (ExecutionContext::Data<ProcessingData> consumer :
         exctx.get_data(consumer_id)) {
      ExportTreeData &tree_data =
          convert_data<ExportTreeData>(consumer.data.get());

      for (ExportTreeData::Entry &entry : tree_data.entries) {
        for (std::size_t producer_id : producer_root_ids) {
          add_dependency(entry, producer_id);
        }
      }
    }
  }

  // Alias the (in-place annotated) inputs under this step's id so that they can
  // be referenced downstream like any other step's outputs.
  detail::AliasPropagator propagator(exctx, step_id);
  std::size_t counter = 0;
  for (const ExecutionContext::Data<ProcessingData> &current : data) {
    propagator.record(current.associated_ids, current.associated_group_ids,
                      counter, 1, /* alias_source = */ true);
    ++counter;
  }
  propagator.flush();

  return counter;
}

}  // namespace sequant::util::extint
