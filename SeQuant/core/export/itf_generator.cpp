#include <SeQuant/core/export/itf_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cctype>
#include <optional>
#include <span>
#include <string>

namespace sequant {

std::string ItfGeneratorContext::index_name(const IndexSpace &space,
                                            std::size_t ordinal) const {
  assert(ordinal >= m_idx_offset);
  ordinal -= m_idx_offset;

  std::string base_key = toUtf8(space.base_key());

  char limit = [&]() -> char {
    auto it = m_index_label_limits.find(space);
    if (it == m_index_label_limits.end()) {
      // Note that in between capital and lowercase letters there are symbols
      // [,\,],^,_ and `. Hence, this default is not quite safe but should work
      // in most cases anyway
      return 'z';
    }

    return it->second;
  }();

  if (base_key.empty() || base_key.size() > 1 ||
      base_key[0] + ordinal > limit) {
    return base_key + std::to_string(ordinal);
  }

  // Merge base key and ordinal into a single letter. That is
  // i1 -> i
  // i2 -> j
  // i3 -> k
  // etc.
  char target = base_key[0] + ordinal;

  return std::string(1, target);
}

std::string ItfGeneratorContext::get_name(const IndexSpace &space) const {
  auto it = m_space_names.find(space);

  if (it == m_space_names.end()) {
    throw std::runtime_error("No name known for index space '" +
                             toUtf8(space.base_key()) + "'");
  }

  return it->second;
}

std::string ItfGeneratorContext::get_tag(const IndexSpace &space) const {
  auto it = m_tags.find(space);

  if (it == m_tags.end()) {
    throw std::runtime_error("No tag known for index space '" +
                             toUtf8(space.base_key()) + "'");
  }

  return it->second;
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Tensor &tensor) const {
  auto it = m_tensor_imports.find(tensor);

  if (it == m_tensor_imports.end()) {
    return {};
  }

  return it->second;
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Variable &variable) const {
  auto it = m_variable_imports.find(variable);

  if (it == m_variable_imports.end()) {
    return {};
  }

  return it->second;
}

void ItfGeneratorContext::set_name(const IndexSpace &space, std::string name) {
  m_space_names[space] = std::move(name);
}

void ItfGeneratorContext::set_tag(const IndexSpace &space, std::string tag) {
  m_tags[space] = std::move(tag);
}

void ItfGeneratorContext::set_import_name(const Tensor &tensor,
                                          std::string name) {
  m_tensor_imports[tensor] = std::move(name);
}

void ItfGeneratorContext::set_import_name(const Variable &variable,
                                          std::string name) {
  m_variable_imports[variable] = std::move(name);
}

bool ItfGeneratorContext::rewrite(Tensor &tensor) const {
  if (tensor.label() == L"R1" || tensor.label() == L"R2" ||
      tensor.label() == L"T1" || tensor.label() == L"T1s" ||
      tensor.label() == L"T2g") {
    return false;
  }

  bool modified = false;

  if (tensor.label() == L"g" && tensor.const_indices().size() == 4 &&
      tensor.bra_rank() == 2 && tensor.ket_rank() == 2) {
    // These are the 2-electron integrals. They posses intra-particle braket
    // symmetry, that is in g^{pq}_{rs} we can exchange p with r, without also
    // swapping q with s. This symmetry can currently not be represented within
    // SeQuant, so we have to use this hacky approach in order to make use of
    // it.

    container::svector<Index> bra = tensor.bra();
    container::svector<Index> ket = tensor.ket();

    for (std::size_t i = 0; i < 2; ++i) {
      if (needs_swap(bra.at(i).space(), ket.at(i).space())) {
        std::swap(bra[i], ket[i]);
      }
    }

    bool swap_cols = false;
    if (bra[0].space() != bra[1].space()) {
      swap_cols = needs_swap(bra[0].space(), bra[1].space());
    } else if (ket[0].space() != ket[1].space()) {
      swap_cols = needs_swap(ket[0].space(), ket[1].space());
    }

    if (swap_cols) {
      std::swap(bra[0], bra[1]);
      std::swap(ket[0], ket[1]);
    }

    std::wstring label;

    // Look at the index space patterns to figure out whether
    // this is a K or a J integral. If the previously attempted sorting
    // of index spaces can be improved by switching the second and third
    // index, do that and thereby produce a J tensor. Otherwise, we retain
    // the index sequence as-is and thereby produce a K tensor.
    // There are some explicit exceptions for J-tensors though.
    Index *p1_1 = nullptr;
    Index *p1_2 = nullptr;
    Index *p2_1 = nullptr;
    Index *p2_2 = nullptr;
    if (is_exceptional_J(bra, ket) ||
        needs_swap(bra[1].space(), ket[0].space())) {
      label = L"J";
      std::swap(bra[1], ket[0]);

      p1_1 = &bra[0];
      p2_1 = &ket[0];
      p1_2 = &bra[1];
      p2_2 = &ket[1];
    } else {
      label = L"K";

      p1_1 = &bra[0];
      p2_1 = &bra[1];
      p1_2 = &ket[0];
      p2_2 = &ket[1];
    }

    // Go through the symmetries again to try and produce the most canonical
    // index ordering possible without breaking the index-space ordering
    // established up to this point.
    // This is a purely cosmetic change, but it is very useful for testing
    // purposes to have  unique representation of the integrals.
    if (p1_1->space() == p1_2->space() && *p1_2 < *p1_1) {
      std::swap(*p1_1, *p1_2);
    }
    if (p2_1->space() == p2_2->space() && *p2_2 < *p2_1) {
      std::swap(*p2_1, *p2_2);
    }
    if (p1_1->space() == p2_1->space() && p1_2->space() == p2_2->space()) {
      if (*p2_1 < *p1_1 || (*p1_1 == *p2_1 && *p2_2 < *p1_2)) {
        std::swap(*p1_1, *p2_1);
        std::swap(*p1_2, *p2_2);
      }
    }

    tensor = Tensor(std::move(label), sequant::bra(std::move(bra)),
                    sequant::ket(std::move(ket)), aux(), tensor.symmetry(),
                    tensor.braket_symmetry(), tensor.particle_symmetry());

    modified = true;
  }

  modified |= ReorderingContext::rewrite(tensor);

  return modified;
}

bool ItfGeneratorContext::is_exceptional_J(std::span<Index> bra,
                                           std::span<Index> ket) const {
  assert(bra.size() == 2);
  assert(ket.size() == 2);

  // integrals with 3 external (virtual) indices ought to be converted to
  // J-integrals
  // Here, we generalize this to all integrals for which the bra indices and
  // the first ket index are of the same space, the ket indices are of different
  // spaces and the additional ket space compares less than the space of the
  // bras.
  bool bras_are_same = bra[0].space() == bra[1].space();
  bool kets_are_different = ket[0].space() != ket[1].space();
  bool kets_are_ordered = is_ordered(ket[0].space(), ket[1].space());
  bool first_ket_same_as_bra = ket[0].space() == bra[0].space();

  return bras_are_same && kets_are_different && kets_are_ordered &&
         first_ket_same_as_bra;
}

std::size_t ItfGeneratorContext::index_id_offset() const {
  return m_idx_offset;
}

void ItfGeneratorContext::set_index_id_offset(std::size_t offset) {
  m_idx_offset = offset;
}

}  // namespace sequant
