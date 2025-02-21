#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/expr_diff.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <bitset>
#include <cassert>
#include <climits>
#include <sstream>
#include <string>

namespace sequant {

// Top-level diff means a diff of the object instance itself without
// regard for any contained subexpressions

template <typename T>
std::string to_string(const Complex<T> &c) {
  std::stringstream stream;

  if (c.imag() == 0) {
    stream << c.real();
  } else if (c.real() == 0) {
    stream << c.imag() << "*i";
  } else if (c.imag() < 0) {
    stream << "(" << c.real() << " - " << (-c.imag()) << "*i)";
  } else {
    stream << "(" << c.real() << " + " << c.imag() << "*i)";
  }

  return stream.str();
}

std::string toplevel_diff(const Constant &lhs, const Constant &rhs) {
  if (lhs == rhs) {
    return {};
  }

  return to_string(lhs.value()) + " vs. " + to_string(rhs.value());
}

std::string toplevel_diff(const Variable &lhs, const Variable &rhs) {
  if (lhs == rhs) {
    return {};
  }

  if (lhs.label() != rhs.label()) {
    return to_string(lhs.label()) + " vs. " + to_string(rhs.label());
  }

  return (lhs.conjugated() ? "conjugated"
                           : "non-conjugated" + std::string(" vs. ")) +
         (rhs.conjugated() ? "conjugated" : "non-conjugated");
}

std::string toplevel_diff(const Index &lhs, const Index &rhs);

template <typename LRange, typename RRange>
std::string diff_indices(const LRange &lhs, const RRange &rhs) {
  auto lhs_size = std::distance(std::begin(lhs), std::end(lhs));
  auto rhs_size = std::distance(std::begin(rhs), std::end(rhs));

  if (lhs_size != rhs_size) {
    return std::to_string(lhs_size) + " indices vs. " +
           std::to_string(rhs_size) + " indices";
  }

  auto lhs_it = std::begin(lhs);
  auto rhs_it = std::begin(rhs);

  std::string diff;

  for (std::size_t i = 0; i < static_cast<std::size_t>(lhs_size); ++i) {
    const Index &lhs_idx = *lhs_it++;
    const Index &rhs_idx = *rhs_it++;

    std::string subdiff = toplevel_diff(lhs_idx, rhs_idx);

    if (subdiff.empty()) {
      continue;
    }

    if (!diff.empty()) {
      diff += ", ";
    }

    diff += "#" + std::to_string(i + 1) + ": " + subdiff;
  }

  return diff;
}

std::string diff_spaces(const IndexSpace &lhs, const IndexSpace &rhs) {
  if (lhs == rhs) {
    return {};
  }

  const auto &lhs_attrs = lhs.attr();
  const auto &rhs_attrs = rhs.attr();

  std::stringstream stream;

  using AttrSet = std::bitset<sizeof(std::uint32_t) * CHAR_BIT>;

  if (lhs_attrs.type() != rhs_attrs.type()) {
    stream << "Types differ: " << AttrSet(lhs.type().to_int32()) << " vs. "
           << AttrSet(rhs.type().to_int32());
  } else if (lhs_attrs.qns() != rhs_attrs.qns()) {
    stream << "QNs differ: " << AttrSet(lhs.qns().to_int32()) << " vs. "
           << AttrSet(rhs.qns().to_int32());
  } else if (lhs.base_key() != rhs.base_key()) {
    stream << "Base key differs: " << to_string(lhs.base_key()) << " vs. "
           << to_string(rhs.base_key());
  } else if (lhs.approximate_size() != rhs.approximate_size()) {
    stream << "Size differs: " << std::to_string(lhs.approximate_size())
           << " vs. " << std::to_string(rhs.approximate_size());
  } else {
    assert(false);
    throw std::runtime_error("Indeterminate space difference");
  }

  assert(!stream.str().empty());
  return stream.str();
}

std::string toplevel_diff(const Index &lhs, const Index &rhs) {
  if (lhs == rhs) {
    return {};
  }

  if (lhs.full_label() != rhs.full_label()) {
    return to_string(lhs.full_label()) + " vs. " + to_string(rhs.full_label());
  }

  if (lhs.space() != rhs.space()) {
    // No string representation of spaces, unfortunately
    return "Spaces differ: " + diff_spaces(lhs.space(), rhs.space());
  }

  if (lhs.has_proto_indices() != rhs.has_proto_indices()) {
    return (lhs.has_proto_indices() ? "with" : "without") +
           std::string(" vs. ") +
           (rhs.has_proto_indices() ? "with" : "without") + " proto-indices";
  }

  if (lhs.proto_indices() != rhs.proto_indices()) {
    return "Proto indices differ: " +
           diff_indices(lhs.proto_indices(), rhs.proto_indices());
  }

  if (lhs.tag() != rhs.tag()) {
    return "Different tags";
  }

  // We have run out of ideas of what to check
  assert(false);
  throw std::runtime_error("Indeterminate index difference");
}

std::string toplevel_diff(const Tensor &lhs, const Tensor &rhs) {
  if (lhs == rhs) {
    return {};
  }

  if (lhs.label() != rhs.label()) {
    return "Names differ: " + to_string(lhs.label()) + " vs. " +
           to_string(rhs.label());
  }

  if (lhs.indices().size() != rhs.indices().size()) {
    return std::to_string(lhs.indices().size()) + " indices vs. " +
           std::to_string(rhs.indices().size()) + " indices";
  }

  if (lhs.symmetry() != rhs.symmetry()) {
    return "Symmetry differs: " + to_string(to_wstring(lhs.symmetry())) +
           " vs. " + to_string(to_wstring(rhs.symmetry()));
  }
  if (lhs.particle_symmetry() != rhs.particle_symmetry()) {
    return "Particle-Symmetry differs: " +
           to_string(to_wstring(lhs.particle_symmetry())) + " vs. " +
           to_string(to_wstring(rhs.particle_symmetry()));
  }
  if (lhs.braket_symmetry() != rhs.braket_symmetry()) {
    return "BraKet-Symmetry differs: " +
           to_string(to_wstring(lhs.braket_symmetry())) + " vs. " +
           to_string(to_wstring(rhs.braket_symmetry()));
  }

  if (lhs.bra() != rhs.bra()) {
    return "Bra indices differ: " + diff_indices(lhs.bra(), rhs.bra());
  }
  if (lhs.ket() != rhs.ket()) {
    return "Ket indices differ: " + diff_indices(lhs.bra(), rhs.bra());
  }
  if (lhs.aux() != rhs.aux()) {
    return "Aux indices differ: " + diff_indices(lhs.ket(), rhs.ket());
  }

  // Really, this shouldn't produce an empty diff as the objects compare as
  // non-equal but we have run out of ideas of what to check
  assert(false);
  throw std::runtime_error("Indeterminate tensor difference");
}

std::string toplevel_diff(const Sum &lhs, const Sum &rhs) {
  // There is no way two Sum objects can be different on the top-level
  return {};
}

std::string toplevel_diff(const Product &lhs, const Product &rhs) {
  if (lhs.scalar() != rhs.scalar()) {
    return "Prefactor differs: " +
           toplevel_diff(Constant(lhs.scalar()), Constant(rhs.scalar()));
  }

  return {};
}

std::string diff(const Expr &lhs, const Expr &rhs) {
  if (lhs == rhs) {
    return {};
  }

  if (lhs.type_id() != rhs.type_id()) {
    return std::string("Types differ: ") + typeid(lhs).name() + " (" +
           std::to_string(lhs.type_id()) + " vs. " + typeid(rhs).name() +
           std::to_string(rhs.type_id());
  }

  auto lhs_begin = std::begin(lhs);
  auto lhs_end = std::end(lhs);
  auto rhs_begin = std::begin(rhs);
  [[maybe_unused]] auto rhs_end = std::end(rhs);

  auto lhs_size = std::distance(lhs_begin, lhs_end);
  auto rhs_size = std::distance(lhs_begin, lhs_end);

  if (lhs_size != rhs_size) {
    return "Sizes differ: " + std::to_string(lhs_size) + " vs. " +
           std::to_string(rhs_size);
  }

  std::string diff_str;

  for (std::size_t i = 0; i < static_cast<std::size_t>(lhs_size); ++i) {
    const Expr &lhs_nested = *(*lhs_begin++);
    const Expr &rhs_nested = *(*rhs_begin++);

    std::string nested_diff = diff(lhs_nested, rhs_nested);

    if (nested_diff.empty()) {
      continue;
    }

    if (diff_str.empty()) {
      diff_str += "Subexpression diff begin:\n";
    }

    diff_str +=
        "Sub-Expr #" + std::to_string(i + 1) + ":\n" + nested_diff + "\n";
  }

  if (!diff_str.empty()) {
    diff_str += "Subexpression diff end";
    return diff_str;
  }

  if (lhs.is<Sum>()) {
    diff_str = toplevel_diff(lhs.as<Sum>(), rhs.as<Sum>());
  } else if (lhs.is<Product>()) {
    diff_str = toplevel_diff(lhs.as<Product>(), rhs.as<Product>());
  } else if (lhs.is<Tensor>()) {
    diff_str = toplevel_diff(lhs.as<Tensor>(), rhs.as<Tensor>());
  } else if (lhs.is<Constant>()) {
    diff_str = toplevel_diff(lhs.as<Constant>(), rhs.as<Constant>());
  } else if (lhs.is<Variable>()) {
    diff_str = toplevel_diff(lhs.as<Variable>(), rhs.as<Variable>());
  } else {
    assert(false);
    throw std::runtime_error("Unsupported expr type");
  }

  return diff_str;
}

}  // namespace sequant
