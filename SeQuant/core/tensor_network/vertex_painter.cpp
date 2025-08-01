#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/tensor_network/vertex_painter.hpp>

namespace sequant {

VertexPainter::VertexPainter(const VertexPainter::NamedIndexSet &named_indices,
                             bool distinct_named_indices)
    : used_colors_(),
      named_indices_(named_indices),
      distinct_named_indices_(distinct_named_indices) {}

std::size_t VertexPainter::to_hash_value(const AbstractTensor &tensor) const {
  auto hashes = {hash::value(tensor._label()),
                 hash::value(tensor._bra_rank()),
                 hash::value(tensor._ket_rank()),
                 hash::value(tensor._aux_rank()),
                 hash::value(tensor._symmetry()),
                 hash::value(tensor._particle_symmetry()),
                 hash::value(tensor._braket_symmetry())};

  return to_hash_value(hashes);
}

VertexPainter::Color VertexPainter::operator()(const AbstractTensor &tensor) {
  Color color = to_color(hash::value(label(tensor)));

  return ensure_uniqueness(color, tensor);
}

VertexPainter::Color VertexPainter::operator()(const BraGroup &group) {
  Color color = to_color(group.id + 0xff);

  return ensure_uniqueness(color, group);
}

VertexPainter::Color VertexPainter::operator()(const KetGroup &group) {
  Color color = to_color(group.id + 0xff00);

  return ensure_uniqueness(color, group);
}

VertexPainter::Color VertexPainter::operator()(const AuxGroup &group) {
  Color color = to_color(group.id + 3 * 0xff0000);

  return ensure_uniqueness(color, group);
}

VertexPainter::Color VertexPainter::operator()(const ParticleGroup &group) {
  Color color;
  if (group.size == 1) {  // legacy coloring works for groups of size 1
    color = to_color(group.id);
  } else {
    color = to_color(
        {group.id, group.size, static_cast<std::size_t>(group.symmetry)});
  }

  return ensure_uniqueness(color, group);
}

void VertexPainter::apply_shade(std::size_t shade) { salt_ = shade; }

VertexPainter::Color VertexPainter::apply_shade(const AbstractTensor &t) {
  const auto hash = to_hash_value(t);
  Color color = to_color(hash);
  salt_ = to_hash_value(t);
  return ensure_uniqueness(color, t);
}

void VertexPainter::reset_shade() { salt_.reset(); }

VertexPainter::Color VertexPainter::operator()(const Index &idx) {
  auto it = named_indices_.find(idx);

  std::size_t pre_color;
  if (it == named_indices_.end()) {
    // anonymous index
    pre_color = idx.color();
  } else {
    if (distinct_named_indices_) {
      pre_color = static_cast<decltype(pre_color)>(
          std::distance(named_indices_.begin(), it));
    } else {  // base colors on Index::color(), but shift to keep distinct from
              // unnamed indices
      pre_color = idx.color() + 0xabcd;
    }
  }
  // shift
  pre_color += 0xaa;

  return ensure_uniqueness(to_color(pre_color), idx);
}

VertexPainter::Color VertexPainter::operator()(const ProtoBundle &bundle) {
  Color color = to_color(Index::proto_indices_color(bundle));

  return ensure_uniqueness(color, bundle);
}

VertexPainter::Color VertexPainter::to_color(std::size_t color) const {
  return to_color({color});
}

std::size_t VertexPainter::to_hash_value(
    std::initializer_list<std::size_t> hash_values) const {
  assert(hash_values.size() > 0);
  std::size_t hash = *(hash_values.begin());
  hash::combine(hash, salt());
  for (auto it = hash_values.begin() + 1; it != hash_values.end(); ++it) {
    hash::combine(hash, *it);
  }
  return hash;
}

VertexPainter::Color VertexPainter::to_color(
    std::initializer_list<std::size_t> colors) const {
  std::size_t color = to_hash_value(colors);

  if constexpr (sizeof(Color) >= sizeof(std::size_t)) {
    return color;
  } else {
    // Need to somehow fit the color into a lower precision integer. In the
    // general case, this is necessarily a lossy conversion. We make the
    // assumption that the input color is
    // - a hash, or
    // - computed from some object ID
    // In the first case, we assume that the used hash function has a uniform
    // distribution or if there is a bias, the bias is towards lower numbers.
    // This allows us to simply reuse the lower x bits of the hash as a new hash
    // (where x == CHAR_BIT * sizeof(VertexColor)). In the second case we assume
    // that such values never exceed the possible value range of VertexColor so
    // that again, we can simply take the lower x bits of color and in this case
    // even retain the numeric value representing the color. Handily, this is
    // exactly what happens when we perform a conversion into a narrower type.
    // We only have to make sure that the underlying types are unsigned as
    // otherwise the behavior is undefined.
    static_assert(std::is_unsigned_v<VertexPainter::Color>,
                  "Narrowing conversion are undefined for signed integers");
    static_assert(std::is_unsigned_v<std::size_t>,
                  "Narrowing conversion are undefined for signed integers");
    return static_cast<Color>(color);
  }
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const AbstractTensor &tensor) {
  return std::holds_alternative<const AbstractTensor *>(data) &&
         label(*std::get<const AbstractTensor *>(data)) == label(tensor);
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const BraGroup &group) {
  return std::holds_alternative<BraGroup>(data) &&
         std::get<BraGroup>(data).id == group.id;
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const KetGroup &group) {
  return std::holds_alternative<KetGroup>(data) &&
         std::get<KetGroup>(data).id == group.id;
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const AuxGroup &group) {
  return std::holds_alternative<AuxGroup>(data) &&
         std::get<AuxGroup>(data).id == group.id;
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const ParticleGroup &group) {
  return std::holds_alternative<ParticleGroup>(data) &&
         std::get<ParticleGroup>(data).id == group.id;
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const Index &idx) {
  if (!std::holds_alternative<Index>(data)) {
    return false;
  }

  const Index &lhs = std::get<Index>(data);

  auto it1 = named_indices_.find(lhs);
  auto it2 = named_indices_.find(idx);

  if (distinct_named_indices_ && it1 != it2) {
    // Either one index is named and the other is not or both are named, but
    // are different indices
    return false;
  }

  return lhs.color() == idx.color();
}

bool VertexPainter::may_have_same_color(const VertexData &data,
                                        const ProtoBundle &bundle) {
  return std::holds_alternative<const ProtoBundle *>(data) &&
         Index::proto_indices_color(*std::get<const ProtoBundle *>(data)) ==
             Index::proto_indices_color(bundle);
}

}  // namespace sequant
