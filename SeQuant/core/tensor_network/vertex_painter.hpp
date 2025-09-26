#ifndef SEQUANT_TENSOR_NETWORK_VERTEX_PAINTER_H
#define SEQUANT_TENSOR_NETWORK_VERTEX_PAINTER_H

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>

#include <utility>
#include <variant>

namespace sequant {

using ProtoBundle =
    std::decay_t<decltype(std::declval<const Index &>().proto_indices())>;

struct BraGroup {
  explicit BraGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
struct KetGroup {
  explicit KetGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
struct AuxGroup {
  explicit AuxGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
/// column in a group of columns produces single braket slot
struct ColumnGroup {
  /// creates a group of particles of size @p size and symmetry @p symmetry
  explicit ColumnGroup(std::size_t id, std::size_t size = 1,
                       Symmetry symmetry = Symmetry::Nonsymm)
      : id(id), size(size), symmetry(symmetry) {}

  std::size_t id;
  std::size_t size;
  Symmetry symmetry;
};

/// Can be used to assign unique colors to a set of objects. The class
/// automatically ensures that there are no accidental color duplications for
/// objects that actually should have different colors (i.e. this is more than a
/// hash function). It is intended to be used to determine the vertex colors in
/// a colored graph representing a tensor network.
class VertexPainterImpl {
 public:
  using Color = tensor_network::VertexColor;
  using NamedIndexSet = tensor_network::NamedIndexSet;
  using VertexData =
      std::variant<const AbstractTensor *, Index, const ProtoBundle *, BraGroup,
                   KetGroup, AuxGroup, ColumnGroup>;
  using ColorMap = container::map<Color, VertexData>;

  /// \param named_indices indices not in this list will use colors based on
  /// their Index::color() \param distinct_named_indices if false, will use same
  /// color for all named indices that have same Index::color(), else will use
  /// distinct color for each named_index
  VertexPainterImpl(const NamedIndexSet &named_indices,
                    bool distinct_named_indices = true);

  const ColorMap &used_colors() const;

  Color operator()(const Index &idx);
  Color operator()(const ProtoBundle &bundle);

  Color operator()(const BraGroup &group);
  Color operator()(const KetGroup &group);
  Color operator()(const AuxGroup &group);
  Color operator()(const ColumnGroup &group);

  /// will apply @p shade to all subsequent colors.
  /// use this to shade colors of all tensor components by the color of its type
  void apply_shade(std::size_t shade);

  /// unlike the deprecated operator() computes a more complete hash
  /// (see VertexPainter::to_hash_value(const AbstractTensor&)) and sets
  /// it as the shade to all subsequent color computation.
  /// Use this before coloring other vertices of this the tensor.
  Color apply_shade(const AbstractTensor &t);

  void reset_shade();

 protected:
  ColorMap used_colors_;
  const NamedIndexSet &named_indices_;
  bool distinct_named_indices_ = true;
  std::optional<std::size_t> salt_;

  /// @return the salt value used for hashing; if not set by apply_shade use
  /// default_salt_
  std::size_t salt() const { return salt_ ? *salt_ : default_salt; }

  // Due to the way we compute the input color, different colors might only
  // differ by a value of 1. This is fine for the algorithmic purpose (after
  // all, colors need only be different - by how much is irrelevant), but
  // sometimes we'll want to use those colors as actual colors to show to a
  // human being. In those cases, having larger differences makes it easier to
  // recognize different colors. Therefore, we hash-combine with an
  // arbitrarily chosen salt with the goal that this will uniformly spread out
  // all input values and therefore increase color differences.
  constexpr static std::size_t default_salt = 0x43d2c59cb15b73f0;

  // combines hashes, injecting salt between first and second hashes
  std::size_t to_hash_value(
      std::initializer_list<std::size_t> hash_values) const;

  /// @brief computes color of a tensor

  /// @note The color should ideally depend on the "type" of the tensor. Unless
  /// we assume that tensor labels are unique (this is too restrictive), the
  /// type is defined by many things, including label, order (i.e. number of
  /// bra, ket, and aux slots), types of each slot, symmetries etc. For tensors
  /// with symmetries in general it is not possible to compute color uniquely
  /// without doing some sort of canonicalization of the slots. To avoid
  /// the need for canonicalization the color will only depend on the
  /// attrributes that are invariant wrt to slot permutations. This means label,
  /// order, and symmetries. More elaborate versions are possible, but have to
  /// stop somewhere.
  /// @sa to_hash_value()
  /// @warning this includes the default_salt
  std::size_t to_hash_value(const AbstractTensor &tensor) const;

  /// produces a (not necessaryly unique) Color from an input size_t-sized color
  Color to_color(std::size_t color) const;

  /// produces a (not necessaryly unique) Color from several input size_t-sized
  /// hash values
  Color to_color(std::initializer_list<std::size_t> hash_values) const;

  template <typename T>
  Color ensure_uniqueness(Color color, const T &val) {
    auto it = used_colors_.find(color);
    while (it != used_colors_.end() && !may_have_same_color(it->second, val)) {
      // Color collision: val was computed to have the same color
      // as another object, but these objects do not compare equal (for
      // the purpose of color assigning).
      // -> Need to modify color until conflict is resolved.
      color++;
      it = used_colors_.find(color);
    }

    if (it == used_colors_.end()) {
      // We have not yet seen this color before -> add it to cache
      if constexpr (std::is_same_v<T, AbstractTensor> ||
                    std::is_same_v<T, ProtoBundle>) {
        used_colors_[color] = &val;
      } else {
        used_colors_[color] = val;
      }
    }

    return color;
  }

  bool may_have_same_color(const VertexData &data,
                           const AbstractTensor &tensor);
  bool may_have_same_color(const VertexData &data, const BraGroup &group);
  bool may_have_same_color(const VertexData &data, const KetGroup &group);
  bool may_have_same_color(const VertexData &data, const AuxGroup &group);
  bool may_have_same_color(const VertexData &data, const ColumnGroup &group);
  bool may_have_same_color(const VertexData &data, const Index &idx);
  bool may_have_same_color(const VertexData &data, const ProtoBundle &bundle);
};

template <typename TN>
class VertexPainter : public VertexPainterImpl {
 public:
  using VertexPainterImpl::VertexPainterImpl;
};

// Template specializations vor TNv1 and TNv2, which still require
// operator()(const AbstractTensor &) (refactoring would require changing tests
// as well)
class TensorNetworkV1;
class TensorNetworkV2;

template <>
class VertexPainter<TensorNetworkV1> : public VertexPainterImpl {
 public:
  using VertexPainterImpl::VertexPainterImpl;

  using VertexPainterImpl::operator();

  VertexPainterImpl::Color operator()(const AbstractTensor &tensor) {
    Color color = to_color(hash::value(label(tensor)));

    return ensure_uniqueness(color, tensor);
  }
};

template <>
class VertexPainter<TensorNetworkV2> : public VertexPainter<TensorNetworkV1> {
 public:
  using VertexPainter<TensorNetworkV1>::VertexPainter;
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_VERTEX_PAINTER_H
