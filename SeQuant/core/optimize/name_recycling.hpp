#ifndef SEQUANT_CORE_OPTIMIZE_NAME_RECYCLING_HPP
#define SEQUANT_CORE_OPTIMIZE_NAME_RECYCLING_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cwctype>
#include <functional>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::opt {

/// Determines under which conditions an intermediate may take over the name
/// (and thereby the storage) of another intermediate
enum class ShapeCompatibility {
  /// The tensors must have the exact same block structure: the same number of
  /// bra, ket and aux slots, the same sequence of index bases (including
  /// proto-index bases) and the same symmetries
  ExactBlock,
  /// In addition to ExactBlock, tensors without any symmetry and without null
  /// slots are compatible if the index bases within each of their bra, ket and
  /// aux slot groups are the same up to permutation. The recycling tensor's
  /// indices are then rearranged to match the layout of the recycled tensor.
  IndexMultiset,
};

/// Describes how to rearrange the slots of an object when it adopts a recycled
/// name: slot k (counting bra, ket, aux slots consecutively) of the rearranged
/// tensor holds what has been slot slot_map[k] of the original tensor. An empty
/// map means that the slots remain in place.
using SlotMap = container::svector<std::size_t>;

namespace detail {

/// @returns The label of @p object which must be a Tensor or a Variable
inline std::wstring_view result_label(const Expr &object) {
  if (object.is<Tensor>()) {
    return object.as<Tensor>().label();
  }

  SEQUANT_ASSERT(object.is<Variable>());
  return object.as<Variable>().label();
}

/// @returns @p label without trailing digits (or @p label itself if it consists
/// of only digits)
inline std::wstring_view label_basename(std::wstring_view label) {
  std::size_t size = label.size();
  while (size > 0 && std::iswdigit(label.at(size - 1))) {
    --size;
  }

  return size == 0 ? label : label.substr(0, size);
}

/// @returns Whether @p lhs and @p rhs have the same index basis, taking
/// proto-indices into account
inline bool same_slot_type(const Index &lhs, const Index &rhs) {
  if (lhs.basis() != rhs.basis() ||
      lhs.proto_indices().size() != rhs.proto_indices().size()) {
    return false;
  }

  return std::ranges::equal(lhs.proto_indices(), rhs.proto_indices(),
                            std::equal_to<>{}, &Index::basis, &Index::basis);
}

inline bool same_symmetries(const Tensor &lhs, const Tensor &rhs) {
  return lhs.symmetry() == rhs.symmetry() &&
         lhs.braket_symmetry() == rhs.braket_symmetry() &&
         lhs.column_symmetry() == rhs.column_symmetry();
}

inline bool same_slot_groups(const Tensor &lhs, const Tensor &rhs) {
  return lhs.bra_rank() == rhs.bra_rank() && lhs.ket_rank() == rhs.ket_rank() &&
         lhs.aux_rank() == rhs.aux_rank();
}

/// @returns Whether @p lhs and @p rhs have the same block structure (ignoring
/// their labels)
inline bool same_block(const Tensor &lhs, const Tensor &rhs) {
  if (!same_slot_groups(lhs, rhs) || !same_symmetries(lhs, rhs)) {
    return false;
  }

  auto lhs_slots = lhs.const_slots();
  auto rhs_slots = rhs.const_slots();

  for (auto lit = lhs_slots.begin(), rit = rhs_slots.begin();
       lit != lhs_slots.end(); ++lit, ++rit) {
    SEQUANT_ASSERT(rit != rhs_slots.end());
    if (lit->nonnull() != rit->nonnull() ||
        (lit->nonnull() && !same_slot_type(*lit, *rit))) {
      return false;
    }
  }

  return true;
}

/// @returns Whether the slots of @p tensor may be freely rearranged within
/// their bra, ket and aux groups without changing what @p tensor represents
inline bool freely_permutable(const Tensor &tensor) {
  return tensor.symmetry() == Symmetry::Nonsymm &&
         tensor.braket_symmetry() == BraKetSymmetry::Nonsymm &&
         tensor.column_symmetry() == ColumnSymmetry::Nonsymm &&
         tensor.num_slots() == tensor.num_indices();
}

/// @returns The SlotMap that rearranges @p candidate to match the layout of
/// @p layout, or std::nullopt if this is not possible under @p compatibility.
/// Labels are not taken into account.
inline std::optional<SlotMap> slot_map_onto(const Tensor &layout,
                                            const Tensor &candidate,
                                            ShapeCompatibility compatibility) {
  if (same_block(layout, candidate)) {
    return SlotMap{};
  }

  if (compatibility != ShapeCompatibility::IndexMultiset ||
      !same_slot_groups(layout, candidate) || !freely_permutable(layout) ||
      !freely_permutable(candidate)) {
    return std::nullopt;
  }

  container::svector<Index> layout_slots(layout.const_slots().begin(),
                                         layout.const_slots().end());
  container::svector<Index> candidate_slots(candidate.const_slots().begin(),
                                            candidate.const_slots().end());

  SlotMap slot_map(layout_slots.size());
  container::svector<bool> taken(candidate_slots.size(), false);

  const std::size_t group_ends[] = {
      layout.bra_rank(), layout.bra_rank() + layout.ket_rank(),
      layout.bra_rank() + layout.ket_rank() + layout.aux_rank()};

  std::size_t group_begin = 0;
  for (std::size_t group_end : group_ends) {
    for (std::size_t k = group_begin; k < group_end; ++k) {
      bool found = false;
      for (std::size_t j = group_begin; j < group_end; ++j) {
        if (!taken[j] && same_slot_type(layout_slots[k], candidate_slots[j])) {
          taken[j] = true;
          slot_map[k] = j;
          found = true;
          break;
        }
      }

      if (!found) {
        return std::nullopt;
      }
    }

    group_begin = group_end;
  }

  return slot_map;
}

/// @returns A copy of @p tensor that is labeled @p label and whose slots are
/// rearranged according to @p slot_map
inline Tensor relabeled(const Tensor &tensor, std::wstring_view label,
                        const SlotMap &slot_map) {
  if (slot_map.empty()) {
    Tensor result = tensor;
    result.set_label(std::wstring(label));
    return result;
  }

  container::svector<Index> slots(tensor.const_slots().begin(),
                                  tensor.const_slots().end());
  SEQUANT_ASSERT(slot_map.size() == slots.size());

  container::svector<Index> rearranged;
  rearranged.reserve(slots.size());
  for (std::size_t source : slot_map) {
    rearranged.push_back(slots.at(source));
  }

  auto bra_end = rearranged.begin() + tensor.bra_rank();
  auto ket_end = bra_end + tensor.ket_rank();

  return Tensor(std::wstring(label),
                bra(container::svector<Index>(rearranged.begin(), bra_end)),
                ket(container::svector<Index>(bra_end, ket_end)),
                aux(container::svector<Index>(ket_end, rearranged.end())),
                TensorSymmetries{.perm = tensor.symmetry(),
                                 .braket = tensor.braket_symmetry(),
                                 .hermiticity = tensor.hermiticity(),
                                 .column = tensor.column_symmetry()});
}

/// @returns A copy of @p variable that is labeled @p label
inline Variable relabeled(const Variable &variable, std::wstring_view label) {
  Variable result = variable;
  result.set_label(std::wstring(label));
  return result;
}

/// Brings a tensor into the form in which it is identified (e.g. by
/// reordering its slots the way a code generator would)
using TensorNormalizer = std::function<void(Tensor &)>;

/// Orders Tensors and Variables by the identity of the object they refer to,
/// i.e. tensors are compared by label and block, variables by label. Supports
/// heterogeneous lookup via Expr objects.
struct ResultObjectLess {
  using is_transparent = void;

  /// If set, tensors are compared in the form produced by this normalizer
  TensorNormalizer normalize = {};

  bool operator()(const Expr &lhs, const Expr &rhs) const {
    if (normalize && lhs.is<Tensor>() && rhs.is<Tensor>()) {
      Tensor normalized_lhs = lhs.as<Tensor>();
      Tensor normalized_rhs = rhs.as<Tensor>();
      normalize(normalized_lhs);
      normalize(normalized_rhs);
      return compare(normalized_lhs, normalized_rhs);
    }

    return compare(lhs, rhs);
  }

  bool operator()(const ExprPtr &lhs, const ExprPtr &rhs) const {
    return (*this)(*lhs, *rhs);
  }
  bool operator()(const ExprPtr &lhs, const Expr &rhs) const {
    return (*this)(*lhs, rhs);
  }
  bool operator()(const Expr &lhs, const ExprPtr &rhs) const {
    return (*this)(lhs, *rhs);
  }

 private:
  static bool compare(const Expr &lhs, const Expr &rhs) {
    const bool lhs_tensor = lhs.is<Tensor>();
    const bool rhs_tensor = rhs.is<Tensor>();

    if (lhs_tensor != rhs_tensor) {
      return lhs_tensor;
    }

    if (lhs_tensor) {
      // Same notion of tensor identity as used by the export
      return TensorBlockLessThanComparator{}(lhs.as<Tensor>(),
                                             rhs.as<Tensor>());
    }

    // A variable and its complex conjugate refer to the same object
    return result_label(lhs) < result_label(rhs);
  }
};

/// Invokes @p visitor on every Tensor and Variable that is (or is the base of)
/// the terminal expression @p leaf. If @p visitor returns a non-null ExprPtr,
/// it replaces the visited object.
/// @returns The (potentially new) expression representing @p leaf
///
/// An adjoint-marked tensor is visited in its bare (non-adjoint) form and the
/// replacement is adjointed again.
template <typename Visitor>
ExprPtr map_leaf_objects(const ExprPtr &leaf, Visitor &&visitor) {
  if (leaf->is<Tensor>() && is_adjoint_label(leaf->as<Tensor>().label())) {
    Tensor bare = leaf->as<Tensor>();
    bare.adjoint();

    ExprPtr replacement = visitor(ex<Tensor>(std::move(bare)));
    if (!replacement) {
      return leaf;
    }

    replacement->as<Tensor>().adjoint();
    return replacement;
  }

  if (leaf->is<Tensor>() || leaf->is<Variable>()) {
    ExprPtr replacement = visitor(leaf);
    return replacement ? replacement : leaf;
  }

  if (leaf->is<Power>()) {
    const Power &power = leaf->as<Power>();
    if (power.base()->is<Tensor>() || power.base()->is<Variable>()) {
      ExprPtr replacement = visitor(power.base());
      if (replacement) {
        auto result = ex<Power>(std::move(replacement), power.exponent());
        if (power.conjugated()) {
          result->as<Power>().conjugate();
        }
        return result;
      }
    }
  }

  if (leaf->is<TensorPower>()) {
    const TensorPower &power = leaf->as<TensorPower>();
    if (power.base()->is<Tensor>()) {
      ExprPtr replacement = visitor(power.base());
      if (replacement) {
        return pow(std::move(replacement), power.exponent());
      }
    }
  }

  return leaf;
}

using ObjectMap = std::map<ExprPtr, std::size_t, ResultObjectLess>;

/// A value of an intermediate as tracked by recycle_intermediate_names
struct RecyclableValue {
  /// The intermediate as it is referred to in the input
  ExprPtr original;
  /// Position of the first tree defining this value
  std::size_t first_def = 0;
  /// Position of the last tree (potentially) using this value
  std::size_t last_use = 0;
  /// Whether the value is also defined by a non-recyclable tree and hence has
  /// to keep its name
  bool pinned = false;
  /// The intermediate as it is referred to in the output
  ExprPtr renamed = {};
  /// How to rearrange the slots of the original intermediate
  SlotMap slot_map = {};
};

/// The values referenced by the leaves of a tree and the value defined by it
struct TreeBindings {
  ObjectMap uses;
  std::optional<std::size_t> defines = {};
};

}  // namespace detail

/// Keeps track of the names of intermediates and how often they are currently
/// in use. Names that are no longer in use are handed out to new intermediates
/// of compatible shape and fresh names are generated such that they don't
/// collide with any reserved or previously used name.
///
/// Tensors are identified by their label and their block, i.e. tensors with
/// the same label but different blocks refer to different objects. Hence, the
/// same label can be handed out for tensors of different blocks.
class IntermediateNamePool {
 public:
  /// Generates a label from a basename and a counter (starting at 1)
  using LabelGenerator =
      std::function<std::wstring(std::wstring_view basename, std::size_t)>;

  /// A recycled tensor along with the description of how the slots of the
  /// original tensor have been rearranged
  struct RecycledTensor {
    Tensor tensor;
    SlotMap slot_map;
  };

  /// The outcome of lowest_available
  struct Allocation {
    std::wstring label;
    /// Whether the label refers to an object that has been in use before
    bool recycled;
  };

  /// @returns A LabelGenerator that appends the counter to the basename
  static LabelGenerator default_label_generator() {
    return [](std::wstring_view basename, std::size_t counter) {
      return std::wstring(basename) + std::to_wstring(counter);
    };
  }

  /// @param label_gen Generator for labels
  /// @param normalize If set, tensors are identified by the form produced by
  /// this normalizer, i.e. tensors whose normalized forms have the same label
  /// and block refer to the same object
  explicit IntermediateNamePool(
      LabelGenerator label_gen = default_label_generator(),
      detail::TensorNormalizer normalize = {})
      : label_gen_(std::move(label_gen)),
        objects_(detail::ResultObjectLess{.normalize = std::move(normalize)}) {}

  /// Marks @p label as being in use by something that isn't managed by this
  /// pool such that it will never be handed out
  void reserve(std::wstring label) { reserved_.insert(std::move(label)); }

  /// Registers a use of @p object (a Tensor or Variable)
  void occupy(const Expr &object) { ++entry(object).uses; }

  /// Ends a use of @p object (a Tensor or Variable). Once all of its uses have
  /// ended, the object's name is available for reuse.
  void release(const Expr &object) {
    State &state = entry(object);
    SEQUANT_ASSERT(state.uses > 0);
    --state.uses;
  }

  /// @param basename The basename to derive the label from via the pool's
  /// LabelGenerator
  /// @param candidate The Tensor or Variable the label is meant for
  /// @returns A label that has not been reserved and that, together with the
  /// block of @p candidate, has not been known to this pool before. The
  /// respective object is marked as being in use.
  std::wstring fresh_label(std::wstring_view basename, const Expr &candidate) {
    ExprPtr object = candidate.clone();

    for (std::size_t counter = 1;; ++counter) {
      std::wstring label = label_gen_(basename, counter);
      if (reserved_.contains(label)) {
        continue;
      }

      set_label(*object, label);
      if (objects_.contains(*object)) {
        continue;
      }

      occupy(*object);
      return label;
    }
  }

  /// Determines the first label (in the order of the pool's LabelGenerator)
  /// that, together with the block of @p candidate, refers to an object that
  /// is not in use and for which @p acceptable (if provided) returns true.
  /// Contrary to fresh_label and acquire, the respective object is not marked
  /// as being in use (see occupy).
  /// @param basename The basename to derive labels from
  /// @param candidate The Tensor or Variable the label is meant for
  /// @param acceptable Predicate invoked with @p candidate relabeled to the
  /// label under consideration
  template <typename T>
    requires std::same_as<T, Tensor> || std::same_as<T, Variable>
  Allocation lowest_available(
      std::wstring_view basename, const T &candidate,
      const std::function<bool(const T &)> &acceptable = {}) const {
    T object = candidate;

    for (std::size_t counter = 1;; ++counter) {
      std::wstring label = label_gen_(basename, counter);
      if (reserved_.contains(label)) {
        continue;
      }

      object.set_label(label);
      auto it = objects_.find(static_cast<const Expr &>(object));
      if ((it != objects_.end() && it->second.uses > 0) ||
          (acceptable && !acceptable(object))) {
        continue;
      }

      return Allocation{.label = std::move(label),
                        .recycled = it != objects_.end()};
    }
  }

  /// Searches for an object that is not in use, whose block is compatible
  /// with @p candidate under @p compatibility and for which @p acceptable (if
  /// provided) returns true. Of all such objects, the one that has become known
  /// to this pool first is chosen and marked as being in use.
  /// @returns @p candidate relabeled (and rearranged) to the recycled name, or
  /// std::nullopt if no suitable name is available
  std::optional<RecycledTensor> acquire(
      const Tensor &candidate, ShapeCompatibility compatibility,
      const std::function<bool(const Tensor &)> &acceptable = {}) {
    const Tensor *best = nullptr;
    std::size_t best_order = 0;
    SlotMap best_map;

    for (const auto &[object, state] : objects_) {
      if (state.uses > 0 || !object->is<Tensor>() ||
          (best && state.order >= best_order)) {
        continue;
      }

      const Tensor &freed = object->as<Tensor>();
      std::optional<SlotMap> slot_map =
          detail::slot_map_onto(freed, candidate, compatibility);
      if (!slot_map.has_value() || (acceptable && !acceptable(freed))) {
        continue;
      }

      best = &freed;
      best_order = state.order;
      best_map = std::move(slot_map.value());
    }

    if (!best) {
      return std::nullopt;
    }

    occupy(*best);

    return RecycledTensor{
        .tensor = detail::relabeled(candidate, best->label(), best_map),
        .slot_map = std::move(best_map)};
  }

  /// Variable equivalent of the Tensor overload of acquire
  std::optional<Variable> acquire(
      const Variable &candidate,
      const std::function<bool(const Variable &)> &acceptable = {}) {
    const Variable *best = nullptr;
    std::size_t best_order = 0;

    for (const auto &[object, state] : objects_) {
      if (state.uses > 0 || !object->is<Variable>() ||
          (best && state.order >= best_order)) {
        continue;
      }

      const Variable &freed = object->as<Variable>();
      if (acceptable && !acceptable(freed)) {
        continue;
      }

      best = &freed;
      best_order = state.order;
    }

    if (!best) {
      return std::nullopt;
    }

    occupy(*best);

    return detail::relabeled(candidate, best->label());
  }

 private:
  struct State {
    /// Time at which the object became known to the pool
    std::size_t order;
    /// Number of current uses
    std::size_t uses;
  };

  LabelGenerator label_gen_;
  std::set<std::wstring, std::less<>> reserved_;
  std::map<ExprPtr, State, detail::ResultObjectLess> objects_;

  static void set_label(Expr &object, std::wstring label) {
    if (object.is<Tensor>()) {
      object.as<Tensor>().set_label(std::move(label));
    } else {
      object.as<Variable>().set_label(std::move(label));
    }
  }

  State &entry(const Expr &object) {
    SEQUANT_ASSERT(object.is<Tensor>() || object.is<Variable>());
    auto it = objects_.find(object);
    if (it == objects_.end()) {
      it = objects_
               .emplace(object.clone(),
                        State{.order = objects_.size(), .uses = 0})
               .first;
    }
    return it->second;
  }
};

/// Options for recycle_intermediate_names
template <typename TreeNode>
struct NameRecyclingOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;

  /// Whether the result of the tree at the given position is an intermediate
  /// whose name may be changed (and recycled). Required.
  std::function<bool(std::size_t, const TreeNode &)> is_recyclable;
  /// Whether the (recyclable) tree at the given position starts a new value
  /// of its result, i.e. overwrites rather than adds to a previous value of
  /// the same object. If not provided, only the first occurrence of an object
  /// as a result starts a new value and later occurrences add to it.
  std::function<bool(std::size_t, const TreeNode &)> starts_new_value = {};
  /// Shape compatibility required for recycling a name
  ShapeCompatibility shape_compatibility = ShapeCompatibility::ExactBlock;
  /// Additional requirement for the value defined by the second tree to take
  /// over the name of the value defined by the first tree (e.g. matching
  /// IRREPs or batching)
  std::function<bool(const TreeNode &previous, const TreeNode &current)>
      extra_compatible = {};
  /// Ranges [begin, end) of tree positions that are executed repeatedly (e.g.
  /// inside a loop). A value that is used inside such a range but is defined
  /// outside of it, stays alive until the end of the range.
  /// @warning Loop-carried uses (a use inside such a range that precedes a
  /// redefinition of the same intermediate in that range, i.e. that reads the
  /// value from the previous iteration) are not supported: such a use is
  /// associated with the definition preceding it in the given order. This
  /// situation is not detected.
  std::vector<std::pair<std::size_t, std::size_t>> repeated_regions = {};
  /// Generates names for intermediates from the basename (label without
  /// trailing digits) of their original label and a counter
  IntermediateNamePool::LabelGenerator label_gen =
      IntermediateNamePool::default_label_generator();
  /// If set, tensors are identified by the form produced by this normalizer
  /// rather than by their form in the trees. This is needed if tensors are
  /// later rewritten (e.g. by a code generator reordering their slots) such
  /// that tensors of different blocks may end up referring to the same object.
  detail::TensorNormalizer normalize = {};
  /// Called for every tree that had to be rebuilt in order to reflect the new
  /// names with the original and the rebuilt tree (e.g. to transfer
  /// annotations)
  std::function<void(const TreeNode &original, TreeNode &rebuilt)> on_rebuilt =
      {};
};

/// Statistics about the outcome of recycle_intermediate_names
struct NameRecyclingResult {
  /// Number of distinct values of recyclable intermediates
  std::size_t num_values = 0;
  /// Number of distinct names used to represent these values
  std::size_t num_names = 0;
};

/// Renames intermediates in the given sequence of evaluation trees such that
/// names (and thereby storage) of intermediates that are no longer needed are
/// reused for later intermediates of compatible shape.
///
/// The trees are interpreted as a program that is executed in the given order,
/// where every tree computes its result (the root's Tensor or Variable) from
/// its leaves. Every value of a recyclable intermediate (see
/// NameRecyclingOptions::is_recyclable and
/// NameRecyclingOptions::starts_new_value) is given a new name, which is either
/// recycled from a value whose last use precedes the value's first definition
/// or freshly generated. Names of all other objects are left untouched and are
/// never generated as fresh names. Leaves referring to a recyclable
/// intermediate are associated with the most recent value of that object
/// defined by a previous tree.
///
/// @note Uses are associated with definitions purely by their order in
/// @p trees. In particular, loop-carried uses of recyclable intermediates
/// within NameRecyclingOptions::repeated_regions are not supported.
///
/// @note The trees must not be reordered afterwards as that may invalidate the
/// liveness analysis performed by this function. Furthermore, the same name
/// can refer to different values in the result. Anything identifying objects
/// by name (e.g. common subexpression elimination) must therefore not be
/// applied to the result.
///
/// @param trees Random-access range of evaluation trees, each representing a
/// ResultExpr
/// @param expr_to_tree Functor that can convert ResultExpr objects into
/// evaluation trees
/// @param opts The options to use
/// @returns Statistics about the renaming
template <std::ranges::random_access_range VectorLike, typename Transformer>
  requires meta::eval_node<std::ranges::range_value_t<VectorLike>> &&
           std::regular_invocable<Transformer, ResultExpr>
NameRecyclingResult recycle_intermediate_names(
    VectorLike &trees, const Transformer &expr_to_tree,
    const NameRecyclingOptions<std::ranges::range_value_t<VectorLike>> &opts) {
  using TreeNode = std::ranges::range_value_t<VectorLike>;
  SEQUANT_ASSERT(opts.is_recyclable);

  using Value = detail::RecyclableValue;
  using TreeBindings = detail::TreeBindings;

  const std::size_t num_trees = std::ranges::size(trees);
  auto tree_at = [&](std::size_t i) -> TreeNode & {
    return *(std::ranges::begin(trees) + i);
  };

  auto live_until = [&](std::size_t use, std::size_t def) {
    std::size_t last = use;
    for (const auto &[begin, end] : opts.repeated_regions) {
      const bool use_inside = begin <= use && use < end;
      const bool def_inside = begin <= def && def < end;
      if (use_inside && !def_inside) {
        last = std::max(last, end - 1);
      }
    }
    return last;
  };

  const detail::ResultObjectLess identity{.normalize = opts.normalize};

  std::vector<Value> values;
  std::vector<TreeBindings> bindings(
      num_trees, TreeBindings{.uses = detail::ObjectMap(identity)});
  std::set<std::wstring> reserved;
  detail::ObjectMap current_values(identity);

  // 1. Determine values and their live ranges
  for (std::size_t i = 0; i < num_trees; ++i) {
    const TreeNode &tree = tree_at(i);
    TreeBindings &current = bindings[i];

    (void)to_expr(tree, [&](const ExprPtr &leaf) {
      return detail::map_leaf_objects(leaf, [&](const ExprPtr &object) {
        if (auto it = current_values.find(object); it != current_values.end()) {
          Value &value = values.at(it->second);
          value.last_use =
              std::max(value.last_use, live_until(i, value.first_def));
          current.uses.emplace(object, it->second);
        } else {
          reserved.emplace(detail::result_label(*object));
        }
        return ExprPtr{};
      });
    });

    const ExprPtr result = tree->expr();
    if (!result->is<Tensor>() && !result->is<Variable>()) {
      continue;
    }

    auto it = current_values.find(result);
    const bool recyclable = opts.is_recyclable(i, tree);

    if (it == current_values.end() && !recyclable) {
      reserved.emplace(detail::result_label(*result));
      continue;
    }

    const bool new_value =
        it == current_values.end() ||
        (recyclable && opts.starts_new_value && opts.starts_new_value(i, tree));

    if (new_value) {
      current_values.insert_or_assign(result, values.size());
      values.push_back(
          Value{.original = result, .first_def = i, .last_use = i});
      current.defines = values.size() - 1;
    } else {
      Value &value = values.at(it->second);
      value.last_use = std::max(value.last_use, live_until(i, value.first_def));
      value.pinned |= !recyclable;
      current.defines = it->second;
    }
  }

  // 2. Assign names
  IntermediateNamePool pool(opts.label_gen, opts.normalize);
  for (const Value &value : values) {
    if (value.pinned) {
      reserved.emplace(detail::result_label(*value.original));
    }
  }
  for (const std::wstring &label : reserved) {
    pool.reserve(label);
  }

  // Values that hold a name, ordered by the end of their live range
  std::multimap<std::size_t, std::size_t> active;
  // Value currently (or most recently) holding a given name
  detail::ObjectMap holder(identity);

  auto previous_def_compatible = [&](const Expr &freed, const Value &value) {
    if (!opts.extra_compatible) {
      return true;
    }
    auto it = holder.find(freed);
    SEQUANT_ASSERT(it != holder.end());
    return opts.extra_compatible(tree_at(values.at(it->second).first_def),
                                 tree_at(value.first_def));
  };

  NameRecyclingResult stats;
  std::set<ExprPtr, detail::ResultObjectLess> used_names(identity);

  for (std::size_t v = 0; v < values.size(); ++v) {
    Value &value = values[v];
    if (value.pinned) {
      continue;
    }

    while (!active.empty() && active.begin()->first < value.first_def) {
      const Value &done = values.at(active.begin()->second);
      pool.release(*done.renamed);
      active.erase(active.begin());
    }

    const std::wstring basename(
        detail::label_basename(detail::result_label(*value.original)));

    if (value.original->is<Tensor>()) {
      const Tensor &original = value.original->as<Tensor>();
      auto recycled = pool.acquire(
          original, opts.shape_compatibility, [&](const Tensor &freed) {
            return previous_def_compatible(freed, value);
          });

      if (recycled.has_value()) {
        value.renamed = ex<Tensor>(std::move(recycled->tensor));
        value.slot_map = std::move(recycled->slot_map);
      } else {
        value.renamed = ex<Tensor>(detail::relabeled(
            original, pool.fresh_label(basename, original), SlotMap{}));
      }
    } else {
      const Variable &original = value.original->as<Variable>();
      auto recycled = pool.acquire(original, [&](const Variable &freed) {
        return previous_def_compatible(freed, value);
      });

      value.renamed = ex<Variable>(
          recycled.has_value()
              ? std::move(recycled.value())
              : detail::relabeled(original,
                                  pool.fresh_label(basename, original)));
    }

    holder.insert_or_assign(value.renamed, v);
    used_names.insert(value.renamed);
    active.emplace(value.last_use, v);
    ++stats.num_values;
  }

  stats.num_names = used_names.size();

  // 3. Rewrite trees
  auto rename_object = [&](const ExprPtr &object, const Value &value) {
    if (object->is<Tensor>()) {
      return ex<Tensor>(detail::relabeled(object->as<Tensor>(),
                                          detail::result_label(*value.renamed),
                                          value.slot_map));
    }
    return ex<Variable>(detail::relabeled(
        object->as<Variable>(), detail::result_label(*value.renamed)));
  };

  auto needs_renaming = [&](std::size_t value) {
    return !values.at(value).pinned;
  };

  for (std::size_t i = 0; i < num_trees; ++i) {
    const TreeBindings &current = bindings[i];

    const bool rename_result =
        current.defines.has_value() && needs_renaming(*current.defines);
    const bool rename_uses =
        std::ranges::any_of(current.uses, needs_renaming,
                            [](const auto &pair) { return pair.second; });

    if (!rename_result && !rename_uses) {
      continue;
    }

    TreeNode &tree = tree_at(i);

    ExprPtr expression = to_expr(tree, [&](const ExprPtr &leaf) {
      return detail::map_leaf_objects(leaf, [&](const ExprPtr &object) {
        auto it = current.uses.find(object);
        if (it == current.uses.end() || !needs_renaming(it->second)) {
          return ExprPtr{};
        }
        return rename_object(object, values.at(it->second));
      });
    });

    ExprPtr result = tree->expr();
    if (rename_result) {
      result = rename_object(result, values.at(*current.defines));
    }

    TreeNode rebuilt = [&]() {
      if (result->is<Tensor>()) {
        return expr_to_tree(
            ResultExpr(result->as<Tensor>(), std::move(expression)));
      }
      return expr_to_tree(
          ResultExpr(result->as<Variable>(), std::move(expression)));
    }();

    if (opts.on_rebuilt) {
      opts.on_rebuilt(tree, rebuilt);
    }

    tree = std::move(rebuilt);
  }

  return stats;
}

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_NAME_RECYCLING_HPP
