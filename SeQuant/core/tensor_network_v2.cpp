//
// Created by Eduard Valeyev on 2019-02-26.
//

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/tag.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/vertex_painter.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/utility/tuple.hpp>
#include <SeQuant/core/wstring.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

auto edge2index(const TensorNetworkV2::Edge &e) -> const Index & {
  return e.idx();
};

struct FullLabelIndexLocator {
  std::wstring_view label;
  FullLabelIndexLocator(std::wstring_view label) : label(std::move(label)) {}

  bool operator()(const TensorNetworkV2::Edge &edge) const {
    return edge.idx().full_label() == label;
  }

  bool operator()(const Index &idx) const { return idx.full_label() == label; }
};

bool tensors_commute(const AbstractTensor &lhs, const AbstractTensor &rhs) {
  // tensors commute if their colors are different or either one of them
  // is a c-number
  return !(color(lhs) == color(rhs) && !is_cnumber(lhs) && !is_cnumber(rhs));
}

struct TensorBlockCompare {
  bool operator()(const AbstractTensor &lhs, const AbstractTensor &rhs) const {
    if (label(lhs) != label(rhs)) {
      return label(lhs) < label(rhs);
    }

    if (bra_rank(lhs) != bra_rank(rhs)) {
      return bra_rank(lhs) < bra_rank(rhs);
    }
    if (ket_rank(lhs) != ket_rank(rhs)) {
      return ket_rank(lhs) < ket_rank(rhs);
    }
    if (aux_rank(lhs) != aux_rank(rhs)) {
      return aux_rank(lhs) < aux_rank(rhs);
    }

    auto lhs_indices = indices(lhs);
    auto rhs_indices = indices(rhs);

    for (auto lhs_it = lhs_indices.begin(), rhs_it = rhs_indices.begin();
         lhs_it != lhs_indices.end() && rhs_it != rhs_indices.end();
         ++lhs_it, ++rhs_it) {
      if (lhs_it->space() != rhs_it->space()) {
        return lhs_it->space() < rhs_it->space();
      }
    }

    // Tensors are identical
    return false;
  }
};

/// Compares tensors based on their label and orders them according to the order
/// of the given cardinal tensor labels. If two tensors can't be discriminated
/// via their label, they are compared based on regular
/// AbstractTensor::operator< or based on their tensor block (the spaces of
/// their indices) - depending on the configuration. If this doesn't
/// discriminate the tensors, they are considered equal
template <typename CardinalLabels>
struct CanonicalTensorCompare {
  const CardinalLabels &labels;
  bool blocks_only;

  CanonicalTensorCompare(const CardinalLabels &labels, bool blocks_only)
      : labels(labels), blocks_only(blocks_only) {}

  void set_blocks_only(bool blocks_only) { this->blocks_only = blocks_only; }

  bool operator()(const AbstractTensorPtr &lhs_ptr,
                  const AbstractTensorPtr &rhs_ptr) const {
    assert(lhs_ptr);
    assert(rhs_ptr);
    const AbstractTensor &lhs = *lhs_ptr;
    const AbstractTensor &rhs = *rhs_ptr;

    if (!tensors_commute(lhs, rhs)) {
      return false;
    }

    const auto get_label = [](const auto &t) {
      if (label(t).back() == adjoint_label) {
        // grab base label if adjoint label is present
        return label(t).substr(0, label(t).size() - 1);
      }
      return label(t);
    };

    const auto lhs_it = std::find(labels.begin(), labels.end(), get_label(lhs));
    const auto rhs_it = std::find(labels.begin(), labels.end(), get_label(rhs));

    if (lhs_it != rhs_it) {
      // At least one of the tensors is a cardinal one
      // -> Order by the occurrence in the cardinal label list
      return std::distance(labels.begin(), lhs_it) <
             std::distance(labels.begin(), rhs_it);
    }

    // Either both are the same cardinal tensor or none is a cardinal tensor
    if (blocks_only) {
      TensorBlockCompare cmp;
      return cmp(lhs, rhs);
    } else {
      return lhs < rhs;
    }
  }

  template <typename AbstractTensorPtr_T>
  bool operator()(const AbstractTensorPtr_T &lhs_ptr,
                  const AbstractTensorPtr_T &rhs_ptr) const {
    return (*this)(lhs_ptr.first, rhs_ptr.first);
  }
};

TensorNetworkV2::Vertex::Vertex(Origin origin, std::size_t terminal_idx,
                                std::size_t index_slot, Symmetry terminal_symm)
    : origin(origin),
      terminal_idx(terminal_idx),
      index_slot(index_slot),
      terminal_symm(terminal_symm) {}

TensorNetworkV2::Origin TensorNetworkV2::Vertex::getOrigin() const {
  return origin;
}

std::size_t TensorNetworkV2::Vertex::getTerminalIndex() const {
  return terminal_idx;
}

std::size_t TensorNetworkV2::Vertex::getIndexSlot() const { return index_slot; }

Symmetry TensorNetworkV2::Vertex::getTerminalSymmetry() const {
  return terminal_symm;
}

bool TensorNetworkV2::Vertex::operator<(const Vertex &rhs) const {
  if (terminal_idx != rhs.terminal_idx) {
    return terminal_idx < rhs.terminal_idx;
  }

  // Both vertices belong to same tensor and are both non-aux? -> they must have
  // same symmetry
  assert(origin != Origin::Aux || rhs.origin != Origin::Aux ||
         terminal_symm == rhs.terminal_symm);

  if (origin != rhs.origin) {
    return origin < rhs.origin;
  }

  // We only take the index slot into account for non-symmetric tensors
  if (terminal_symm == Symmetry::nonsymm) {
    return index_slot < rhs.index_slot;
  } else {
    return false;
  }
}

bool TensorNetworkV2::Vertex::operator==(const Vertex &rhs) const {
  // Slot position is only taken into account for non_symmetric tensors
  const std::size_t lhs_slot =
      (terminal_symm == Symmetry::nonsymm) * index_slot;
  const std::size_t rhs_slot =
      (rhs.terminal_symm == Symmetry::nonsymm) * rhs.index_slot;

  // sanity check that bra and ket have same symmetry
  assert(origin == Origin::Aux || rhs.origin == Origin::Aux ||
         terminal_idx != rhs.terminal_idx ||
         terminal_symm == rhs.terminal_symm);

  return terminal_idx == rhs.terminal_idx && lhs_slot == rhs_slot &&
         origin == rhs.origin;
}

std::size_t TensorNetworkV2::Graph::vertex_to_index_idx(
    std::size_t vertex) const {
  assert(vertex_types.at(vertex) == VertexType::Index);

  std::size_t index_idx = 0;
  for (std::size_t i = 0; i <= vertex; ++i) {
    if (vertex_types[i] == VertexType::Index) {
      ++index_idx;
    }
  }

  assert(index_idx > 0);

  return index_idx - 1;
}

std::size_t TensorNetworkV2::Graph::vertex_to_tensor_idx(
    std::size_t vertex) const {
  assert(vertex_types.at(vertex) == VertexType::TensorCore);

  std::size_t tensor_idx = 0;
  for (std::size_t i = 0; i <= vertex; ++i) {
    if (vertex_types[i] == VertexType::TensorCore) {
      ++tensor_idx;
    }
  }

  assert(tensor_idx > 0);

  return tensor_idx - 1;
}

template <typename ArrayLike, typename Permutation>
auto permute(const ArrayLike &vector, const Permutation &perm) {
  using std::size;
  auto sz = size(vector);
  std::decay_t<decltype(vector)> pvector(sz);
  for (size_t i = 0; i != sz; ++i) pvector[perm[i]] = vector[i];
  return pvector;
}

template <typename ReplacementMap>
void apply_index_replacements(AbstractTensor &tensor,
                              const ReplacementMap &replacements) {
#ifndef NDEBUG
  // assert that tensors' indices are not tagged since going to tag indices
  assert(ranges::none_of(
      indices(tensor), [](const Index &idx) { return idx.tag().has_value(); }));
#endif

  bool pass_mutated;
  do {
    pass_mutated = transform_indices(tensor, replacements);
  } while (pass_mutated);  // transform till stops changing

  reset_tags(tensor);
}

template <typename ArrayLike, typename ReplacementMap>
void apply_index_replacements(ArrayLike &tensors,
                              const ReplacementMap &replacements) {
  for (auto &tensor : tensors) {
    apply_index_replacements(*tensor, replacements);
  }
}

template <typename Container>
void order_to_indices(Container &container) {
  std::vector<std::size_t> indices;
  indices.resize(container.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::sort(indices.begin(), indices.end(),
            [&container](std::size_t lhs, std::size_t rhs) {
              return container[lhs] < container[rhs];
            });
  // Overwrite container contents with indices
  std::copy(indices.begin(), indices.end(), container.begin());
}

template <bool stable, typename Container, typename Comparator>
void sort_via_indices(Container &container, const Comparator &cmp) {
  std::vector<std::size_t> indices;
  indices.resize(container.size());
  std::iota(indices.begin(), indices.end(), 0);

  if constexpr (stable) {
    std::stable_sort(indices.begin(), indices.end(), cmp);
  } else {
    std::sort(indices.begin(), indices.end(), cmp);
  }

  // Bring elements in container into the order given by indices
  // (the association is container[k] = container[indices[k]])
  // -> implementation from https://stackoverflow.com/a/838789

  for (std::size_t i = 0; i < container.size(); ++i) {
    if (indices[i] == i) {
      // This element is already where it is supposed to be
      continue;
    }

    // Find the offset of the index pointing to i
    // -> since we are going to change the content of the vector at position i,
    // we have to update the index-mapping referencing i to point to the new
    // location of the element that used to be at position i
    std::size_t k;
    for (k = i + 1; k < container.size(); ++k) {
      if (indices[k] == i) {
        break;
      }
    }
    std::swap(container[i], container[indices[i]]);
    std::swap(indices[i], indices[k]);
  }
}

void TensorNetworkV2::canonicalize_graph(const NamedIndexSet &named_indices) {
  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV2::canonicalize_graph: input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    std::wcout << std::endl;
  }

  if (!have_edges_) {
    init_edges();
  }

  const auto is_anonymous_index = [named_indices](const Index &idx) {
    return named_indices.find(idx) == named_indices.end();
  };

  // index factory to generate anonymous indices
  IndexFactory idxfac(is_anonymous_index, 1);

  // make the graph
  Graph graph = create_graph(&named_indices);
  // graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);

  if (Logger::instance().canonicalize_input_graph) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);
  }

  // canonize the graph
  bliss::Stats stats;
  graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
  const unsigned int *canonize_perm =
      graph.bliss_graph->canonical_form(stats, nullptr, nullptr);

  if (Logger::instance().canonicalize_dot) {
    std::wcout << "Canonicalization permutation:\n";
    for (std::size_t i = 0; i < graph.vertex_labels.size(); ++i) {
      std::wcout << i << " -> " << canonize_perm[i] << "\n";
    }
    std::wcout << "Canonicalized graph:\n";
    bliss::Graph *cgraph = graph.bliss_graph->permute(canonize_perm);
    cgraph->write_dot(std::wcout, {}, true);
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    std::wcout << "with our labels:\n";
    cgraph->write_dot(std::wcout, cvlabels);
    delete cgraph;
  }

  container::map<std::size_t, std::size_t> tensor_idx_to_vertex;
  container::map<std::size_t, container::svector<std::size_t, 3>>
      tensor_idx_to_particle_order;
  container::map<std::size_t, std::size_t> index_idx_to_vertex;
  std::size_t tensor_idx = 0;
  std::size_t index_idx = 0;

  for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
    switch (graph.vertex_types[vertex]) {
      case VertexType::Index:
        index_idx_to_vertex[index_idx] = vertex;
        index_idx++;
        break;
      case VertexType::Particle: {
        assert(tensor_idx > 0);
        const std::size_t base_tensor_idx = tensor_idx - 1;
        assert(symmetry(*tensors_.at(base_tensor_idx)) == Symmetry::nonsymm);
        tensor_idx_to_particle_order[base_tensor_idx].push_back(
            canonize_perm[vertex]);
        break;
      }
      case VertexType::TensorCore:
        tensor_idx_to_vertex[tensor_idx] = vertex;
        tensor_idx++;
        break;
      case VertexType::TensorBra:
      case VertexType::TensorKet:
      case VertexType::TensorAux:
      case VertexType::SPBundle:
        break;
    }
  }

  assert(index_idx_to_vertex.size() ==
         edges_.size() + pure_proto_indices_.size());
  assert(tensor_idx_to_vertex.size() == tensors_.size());
  assert(tensor_idx_to_particle_order.size() <= tensors_.size());

  // order_to_indices(index_order);
  for (auto &current : tensor_idx_to_particle_order) {
    order_to_indices(current.second);
  }

  container::map<Index, Index> idxrepl;
  // Sort edges so that their order corresponds to the order of indices in the
  // canonical graph
  // Use this ordering to relabel anonymous indices
  const auto index_sorter = [&index_idx_to_vertex, &canonize_perm](
                                std::size_t lhs_idx, std::size_t rhs_idx) {
    const std::size_t lhs_vertex = index_idx_to_vertex.at(lhs_idx);
    const std::size_t rhs_vertex = index_idx_to_vertex.at(rhs_idx);

    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  sort_via_indices<false>(edges_, index_sorter);

  for (const Edge &current : edges_) {
    const Index &idx = current.idx();

    if (!is_anonymous_index(idx)) {
      continue;
    }

    idxrepl.insert(std::make_pair(idx, idxfac.make(idx)));
  }

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetworkV2::canonicalize_graph: replacing "
                 << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  apply_index_replacements(tensors_, idxrepl);

  // Perform particle-1,2-swaps as indicated by the graph canonization
  for (std::size_t i = 0; i < tensors_.size(); ++i) {
    AbstractTensor &tensor = *tensors_[i];
    const std::size_t num_particles =
        std::min(bra_rank(tensor), ket_rank(tensor));

    auto it = tensor_idx_to_particle_order.find(i);
    if (it == tensor_idx_to_particle_order.end()) {
      assert(num_particles == 0 || symmetry(*tensors_[i]) != Symmetry::nonsymm);
      continue;
    }

    const auto &particle_order = it->second;
    auto bra_indices = tensor._bra();
    auto ket_indices = tensor._ket();

    assert(num_particles == particle_order.size());

    // Swap indices column-wise
    idxrepl.clear();
    for (std::size_t col = 0; col < num_particles; ++col) {
      if (particle_order[col] == col) {
        continue;
      }

      idxrepl.insert(
          std::make_pair(bra_indices[col], bra_indices[particle_order[col]]));
      idxrepl.insert(
          std::make_pair(ket_indices[col], ket_indices[particle_order[col]]));
    }

    if (!idxrepl.empty()) {
      if (Logger::instance().canonicalize) {
        for (const auto &idxpair : idxrepl) {
          std::wcout
              << "TensorNetworkV2::canonicalize_graph: permuting particles in "
              << to_latex(tensor) << " by replacing " << to_latex(idxpair.first)
              << " with " << to_latex(idxpair.second) << std::endl;
        }
      }
      apply_index_replacements(tensor, idxrepl);
    }
  }

  // Bring tensors into canonical order (analogously to how we reordered
  // indices), but ensure to respect commutativity!
  const auto tensor_sorter = [this, &canonize_perm, &tensor_idx_to_vertex](
                                 std::size_t lhs_idx, std::size_t rhs_idx) {
    const AbstractTensor &lhs = *tensors_[lhs_idx];
    const AbstractTensor &rhs = *tensors_[rhs_idx];

    if (!tensors_commute(lhs, rhs)) {
      return false;
    }

    const std::size_t lhs_vertex = tensor_idx_to_vertex.at(lhs_idx);
    const std::size_t rhs_vertex = tensor_idx_to_vertex.at(rhs_idx);

    // Commuting tensors are sorted based on their canonical order which is
    // given by the order of the corresponding vertices in the canonical graph
    // representation
    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  sort_via_indices<true>(tensors_, tensor_sorter);

  // The tensor reordering and index relabelling made the current set of edges
  // invalid
  edges_.clear();
  have_edges_ = false;

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV2::canonicalize_graph: tensors after "
                  "canonicalization\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
  }
}

ExprPtr TensorNetworkV2::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels, bool fast,
    const NamedIndexSet *named_indices_ptr) {
  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV2::canonicalize(" << (fast ? "fast" : "slow")
               << "): input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    std::wcout << "cardinal_tensor_labels = ";
    ranges::for_each(cardinal_tensor_labels,
                     [](auto &&i) { std::wcout << i << L" "; });
    std::wcout << std::endl;
  }

  if (!have_edges_) {
    init_edges();
  }

  // initialize named_indices by default to all external indices
  const auto &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  if (!fast) {
    // The graph-based canonization is required in call cases in which there are
    // indistinguishable tensors present in the expression. Their order and
    // indexing can only be determined via this rigorous canonization.
    canonicalize_graph(named_indices);
  }

  // Ensure each individual tensor is written in the way that its tensor
  // block (== order of index spaces) is canonical
  ExprPtr byproduct = canonicalize_individual_tensor_blocks(named_indices);

  CanonicalTensorCompare<decltype(cardinal_tensor_labels)> tensor_sorter(
      cardinal_tensor_labels, true);

  using ranges::begin;
  using ranges::end;
  using ranges::views::zip;
  auto tensors_with_ordinals = zip(tensors_, tensor_input_ordinals_);
  bubble_sort(begin(tensors_with_ordinals), end(tensors_with_ordinals),
              tensor_sorter);

  init_edges();

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV2::canonicalize(" << (fast ? "fast" : "slow")
               << "): tensors after initial sort\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
  }

  // helpers to filter named ("external" in traditional use case) / anonymous
  // ("internal" in traditional use case)
  auto is_named_index = [&](const Index &idx) {
    return named_indices.find(idx) != named_indices.end();
  };
  auto is_anonymous_index = [&](const Index &idx) {
    return named_indices.find(idx) == named_indices.end();
  };

  // Sort edges based on the order of the tensors they connect
  std::stable_sort(edges_.begin(), edges_.end(),
                   [&is_named_index](const Edge &lhs, const Edge &rhs) {
                     // Sort first by index's character (named < anonymous),
                     // then by Edge (not by Index's full label) ... this
                     // automatically puts named indices first
                     const bool lhs_is_named = is_named_index(lhs.idx());
                     const bool rhs_is_named = is_named_index(rhs.idx());

                     if (lhs_is_named == rhs_is_named) {
                       return lhs < rhs;
                     } else {
                       return lhs_is_named;
                     }
                   });

  // index factory to generate anonymous indices
  // -> start reindexing anonymous indices from 1
  IndexFactory idxfac(is_anonymous_index, 1);

  container::map<Index, Index> idxrepl;

  // Use the new order of edges as the canonical order of indices and relabel
  // accordingly (but only anonymous indices, of course)
  for (std::size_t i = named_indices.size(); i < edges_.size(); ++i) {
    const Index &index = edges_[i].idx();
    assert(is_anonymous_index(index));
    Index replacement = idxfac.make(index);
    idxrepl.emplace(std::make_pair(index, replacement));
  }

  // Done computing canonical index replacement list

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetworkV2::canonicalize(" << (fast ? "fast" : "slow")
                 << "): replacing " << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  apply_index_replacements(tensors_, idxrepl);

  byproduct *= canonicalize_individual_tensors(named_indices);

  // We assume that re-indexing did not change the canonical order of tensors
  assert(std::is_sorted(tensors_.begin(), tensors_.end(), tensor_sorter));
  // However, in order to produce the most aesthetically pleasing result, we now
  // reorder tensors based on the regular AbstractTensor::operator<, which takes
  // the explicit index labelling of tensors into account.
  tensor_sorter.set_blocks_only(false);
  std::stable_sort(tensors_.begin(), tensors_.end(), tensor_sorter);

  have_edges_ = false;

  assert(byproduct->is<Constant>());
  return (byproduct->as<Constant>().value() == 1) ? nullptr : byproduct;
}

TensorNetworkV2::SlotCanonicalizationMetadata
TensorNetworkV2::canonicalize_slots(
    const container::vector<std::wstring> &cardinal_tensor_labels,
    const NamedIndexSet *named_indices_ptr,
    TensorNetworkV2::SlotCanonicalizationMetadata::named_index_compare_t
        named_index_compare) {
  if (!named_index_compare)
    named_index_compare = [](const auto &idxptr_slottype_1,
                             const auto &idxptr_slottype_2) -> bool {
      const auto &[idxptr1, slottype1] = idxptr_slottype_1;
      const auto &[idxptr2, slottype2] = idxptr_slottype_2;
      return idxptr1->space() < idxptr2->space();
    };

  TensorNetworkV2::SlotCanonicalizationMetadata metadata;

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV2::canonicalize_slots(): input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    std::wcout << "cardinal_tensor_labels = ";
    ranges::for_each(cardinal_tensor_labels,
                     [](auto &&i) { std::wcout << i << L" "; });
    std::wcout << std::endl;
  }

  if (!have_edges_) {
    init_edges();
  }

  // initialize named_indices by default to all external indices
  const auto &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;
  metadata.named_indices = named_indices;

  // helper to filter named ("external" in traditional use case) / anonymous
  // ("internal" in traditional use case)
  auto is_named_index = [&](const Index &idx) {
    return named_indices.find(idx) != named_indices.end();
  };

  // make the graph
  // only slots (hence, attr) of named indices define their color, so
  // distinct_named_indices = false
  Graph graph =
      create_graph(&named_indices, /* distinct_named_indices = */ false);
  // graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);

  if (Logger::instance().canonicalize_input_graph) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);
  }

  // canonize the graph
  bliss::Stats stats;
  graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
  const unsigned int *canonize_perm =
      graph.bliss_graph->canonical_form(stats, nullptr, nullptr);

  metadata.graph =
      std::shared_ptr<bliss::Graph>(graph.bliss_graph->permute(canonize_perm));

  if (Logger::instance().canonicalize_dot) {
    std::wcout << "Canonicalization permutation:\n";
    for (std::size_t i = 0; i < graph.vertex_labels.size(); ++i) {
      std::wcout << i << " -> " << canonize_perm[i] << "\n";
    }
    std::wcout << "Canonicalized graph:\n";
    metadata.graph->write_dot(std::wcout, {}, true);
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    std::wcout << "with our labels:\n";
    metadata.graph->write_dot(std::wcout, cvlabels);
  }

  // maps index ordinal to vertex ordinal
  container::map<std::size_t, std::size_t> index_idx_to_vertex;
  std::size_t index_idx = 0;
  for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
    if (graph.vertex_types[vertex] == VertexType::Index) {
      index_idx_to_vertex[index_idx] = vertex;
      index_idx++;
    }
  }
  assert(index_idx_to_vertex.size() ==
         edges_.size() + pure_proto_indices_.size());

  // produce canonical list of named indices
  {
    using ord_cord_it_t =
        std::tuple<size_t, size_t, NamedIndexSet::const_iterator>;
    using cord_set_t = container::set<ord_cord_it_t, detail::tuple_less<1>>;

    auto grand_index_list = ranges::views::concat(
        edges_ | ranges::views::transform(edge2index), pure_proto_indices_);

    // for each named index type (as defined by named_index_compare) maps its
    // ptr in grand_index_list to its ordinal in grand_index_list + canonical
    // ordinal + its iterator in metadata.named_indices
    container::map<std::pair<const Index *, IndexSlotType>, cord_set_t,
                   decltype(named_index_compare)>
        idx2cord(named_index_compare);

    // collect named indices and sort them on the fly
    size_t idx_ord = 0;
    auto grand_index_list_end = grand_index_list.end();
    for (auto git = grand_index_list.begin(); git != grand_index_list_end;
         ++git) {
      const auto &idx = *git;

      if (is_named_index(idx)) {
        const auto vertex_ord = index_idx_to_vertex[idx_ord];
        const auto named_indices_it = metadata.named_indices.find(idx);
        assert(named_indices_it != metadata.named_indices.end());

        // find the entry for this index type
        IndexSlotType slot_type;
        if (idx_ord < edges_.size()) {
          auto edge_it = edges_.begin();
          std::advance(edge_it, idx_ord);
          // for named indices edges are always disconnected
          assert(edge_it->vertex_count() == 1);
          if (edge_it->first_vertex().getOrigin() == Origin::Aux)
            slot_type = IndexSlotType::TensorAux;
          else if (edge_it->first_vertex().getOrigin() == Origin::Bra)
            slot_type = IndexSlotType::TensorBra;
          else {
            assert(edge_it->first_vertex().getOrigin() == Origin::Ket);
            slot_type = IndexSlotType::TensorKet;
          }
        } else
          slot_type = IndexSlotType::SPBundle;
        const auto idxptr_slottype = std::make_pair(&idx, slot_type);
        auto it = idx2cord.find(idxptr_slottype);

        if (it == idx2cord.end()) {
          bool inserted;
          std::tie(it, inserted) = idx2cord.insert(std::make_pair(
              idxptr_slottype, cord_set_t(cord_set_t::key_compare{})));
          assert(inserted);
        }
        it->second.emplace(idx_ord, canonize_perm[vertex_ord],
                           named_indices_it);
      }
      ++idx_ord;
    }

    // save the result
    for (auto &[idxptr_slottype, cord_set] : idx2cord) {
      for (auto &[idx_ord, idx_ord_can, named_indices_it] : cord_set) {
        metadata.named_indices_canonical.emplace_back(named_indices_it);
      }
    }
    metadata.named_index_compare = std::move(named_index_compare);

  }  // named indices resort to canonical order

  return metadata;
}

TensorNetworkV2::Graph TensorNetworkV2::create_graph(
    const NamedIndexSet *named_indices_ptr, bool distinct_named_indices) const {
  assert(have_edges_);

  // initialize named_indices by default to all external indices
  const NamedIndexSet &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  VertexPainter colorizer(named_indices, distinct_named_indices);

  // core, bra, ket, auxiliary and optionally (for non-symmetric tensors) a
  // particle vertex
  constexpr std::size_t num_tensor_components = 5;

  // results
  Graph graph;
  // We know that at the very least all indices and all tensors will yield
  // vertex representations
  std::size_t vertex_count_estimate = edges_.size() +
                                      pure_proto_indices_.size() +
                                      num_tensor_components * tensors_.size();
  graph.vertex_labels.reserve(vertex_count_estimate);
  graph.vertex_colors.reserve(vertex_count_estimate);
  graph.vertex_types.reserve(vertex_count_estimate);

  container::map<ProtoBundle, std::size_t> proto_bundles;

  container::map<std::size_t, std::size_t> tensor_vertices;
  tensor_vertices.reserve(tensors_.size());

  container::vector<std::pair<std::size_t, std::size_t>> edges;
  edges.reserve(edges_.size() + tensors_.size());

  // Add vertices for tensors
  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    assert(tensor_vertices.find(tensor_idx) == tensor_vertices.end());
    assert(tensors_.at(tensor_idx));
    const AbstractTensor &tensor = *tensors_.at(tensor_idx);

    // Tensor core
    graph.vertex_labels.emplace_back(label(tensor));
    graph.vertex_types.emplace_back(VertexType::TensorCore);
    graph.vertex_colors.push_back(colorizer(tensor));

    const std::size_t tensor_vertex = graph.vertex_labels.size() - 1;
    tensor_vertices.insert(std::make_pair(tensor_idx, tensor_vertex));

    // Create vertices to group indices
    const Symmetry tensor_sym = symmetry(tensor);
    if (tensor_sym == Symmetry::nonsymm) {
      // Create separate vertices for every index
      // Additionally, we need particle vertices to group indices that belong to
      // the same particle (are in the same "column" in the usual tensor
      // notation)
      const std::size_t num_particle_vertices =
          std::min(bra_rank(tensor), ket_rank(tensor));
      const bool is_part_symm =
          particle_symmetry(tensor) == ParticleSymmetry::symm;
      // TODO: How to handle BraKetSymmetry::conjugate?
      const bool is_braket_symm =
          braket_symmetry(tensor) == BraKetSymmetry::symm;

      for (std::size_t i = 0; i < num_particle_vertices; ++i) {
        graph.vertex_labels.emplace_back(L"p_" + std::to_wstring(i + 1));
        graph.vertex_types.push_back(VertexType::Particle);
        // Particles are indistinguishable -> always use same ID
        graph.vertex_colors.push_back(colorizer(ParticleGroup{0}));
        edges.push_back(
            std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
      }

      for (std::size_t i = 0; i < bra_rank(tensor); ++i) {
        const bool is_unpaired_idx = i >= num_particle_vertices;
        const bool color_idx = is_unpaired_idx || !is_part_symm;

        graph.vertex_labels.emplace_back(L"bra_" + std::to_wstring(i + 1));
        graph.vertex_types.push_back(VertexType::TensorBra);
        graph.vertex_colors.push_back(colorizer(BraGroup{color_idx ? i : 0}));

        const std::size_t connect_vertex =
            tensor_vertex + (is_unpaired_idx ? 0 : (i + 1));
        edges.push_back(
            std::make_pair(connect_vertex, graph.vertex_labels.size() - 1));
      }

      for (std::size_t i = 0; i < ket_rank(tensor); ++i) {
        const bool is_unpaired_idx = i >= num_particle_vertices;
        const bool color_idx = is_unpaired_idx || !is_part_symm;

        graph.vertex_labels.emplace_back(L"ket_" + std::to_wstring(i + 1));
        graph.vertex_types.push_back(VertexType::TensorKet);
        if (is_braket_symm) {
          // Use BraGroup for kets as well as they are supposed to be
          // indistinguishable
          graph.vertex_colors.push_back(colorizer(BraGroup{color_idx ? i : 0}));
        } else {
          graph.vertex_colors.push_back(colorizer(KetGroup{color_idx ? i : 0}));
        }

        const std::size_t connect_vertex =
            tensor_vertex + (is_unpaired_idx ? 0 : (i + 1));
        edges.push_back(
            std::make_pair(connect_vertex, graph.vertex_labels.size() - 1));
      }
    } else {
      // Shared set of bra/ket vertices for all indices
      std::wstring suffix = tensor_sym == Symmetry::symm ? L"_s" : L"_a";

      graph.vertex_labels.push_back(L"bra" + suffix);
      graph.vertex_types.push_back(VertexType::TensorBra);
      graph.vertex_colors.push_back(colorizer(BraGroup{0}));
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));

      graph.vertex_labels.push_back(L"ket" + suffix);
      graph.vertex_types.push_back(VertexType::TensorKet);
      // TODO: figure out how to handle BraKetSymmetry::conjugate
      if (braket_symmetry(tensor) == BraKetSymmetry::symm) {
        // Use BraGroup for kets as well as they should be indistinguishable
        graph.vertex_colors.push_back(colorizer(BraGroup{0}));
      } else {
        graph.vertex_colors.push_back(colorizer(KetGroup{0}));
      }
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
    }

    // TODO: handle aux indices permutation symmetries once they are supported
    // for now, auxiliary indices are considered to always be asymmetric
    for (std::size_t i = 0; i < aux_rank(tensor); ++i) {
      graph.vertex_labels.emplace_back(L"aux_" + std::to_wstring(i + 1));
      graph.vertex_types.push_back(VertexType::TensorAux);
      graph.vertex_colors.push_back(colorizer(AuxGroup{i}));
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
    }
  }

  // Now add all indices (edges_ + pure_proto_indices_) to the graph
  container::map<Index, std::size_t> index_vertices;

  for (const Edge &current_edge : edges_) {
    const Index &index = current_edge.idx();
    graph.vertex_labels.push_back(std::wstring(index.full_label()));
    graph.vertex_types.push_back(VertexType::Index);
    graph.vertex_colors.push_back(colorizer(index));

    const std::size_t index_vertex = graph.vertex_labels.size() - 1;

    index_vertices[index] = index_vertex;

    // Handle proto indices
    if (index.has_proto_indices()) {
      // For now we assume that all proto indices are symmetric
      assert(index.symmetric_proto_indices());

      std::size_t proto_vertex;
      if (auto it = proto_bundles.find(index.proto_indices());
          it != proto_bundles.end()) {
        proto_vertex = it->second;
      } else {
        using namespace std::literals;
        // Create a new vertex for this bundle of proto indices
        std::wstring spbundle_label =
            L"<" +
            (ranges::views::transform(
                 index.proto_indices(),
                 [](const Index &idx) { return idx.full_label(); }) |
             ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
            L">";

        graph.vertex_labels.push_back(std::move(spbundle_label));
        graph.vertex_types.push_back(VertexType::SPBundle);
        graph.vertex_colors.push_back(colorizer(index.proto_indices()));

        proto_vertex = graph.vertex_labels.size() - 1;
        proto_bundles.insert(
            std::make_pair(index.proto_indices(), proto_vertex));
      }

      edges.push_back(std::make_pair(index_vertex, proto_vertex));
    }

    // Connect index to the tensor(s) it is connected to
    for (std::size_t i = 0; i < current_edge.vertex_count(); ++i) {
      assert(i <= 1);
      const Vertex &vertex =
          i == 0 ? current_edge.first_vertex() : current_edge.second_vertex();

      assert(tensor_vertices.find(vertex.getTerminalIndex()) !=
             tensor_vertices.end());
      const std::size_t tensor_vertex =
          tensor_vertices.find(vertex.getTerminalIndex())->second;

      // Store an edge connecting the index vertex to the corresponding tensor
      // vertex
      const bool tensor_is_nonsymm =
          vertex.getTerminalSymmetry() == Symmetry::nonsymm;
      const AbstractTensor &tensor = *tensors_[vertex.getTerminalIndex()];
      std::size_t offset;
      if (tensor_is_nonsymm) {
        // We have to find the correct vertex to connect this index to (for
        // non-symmetric tensors each index has its dedicated "group" vertex)

        // Move off the tensor core's vertex
        offset = 1;
        // Move past the explicit particle vertices
        offset += std::min(bra_rank(tensor), ket_rank(tensor));

        if (vertex.getOrigin() > Origin::Bra) {
          offset += bra_rank(tensor);
        }

        offset += vertex.getIndexSlot();
      } else {
        static_assert(static_cast<int>(Origin::Bra) == 1);
        static_assert(static_cast<int>(Origin::Ket) == 2);
        static_assert(static_cast<int>(Origin::Aux) == 3);
        offset = static_cast<std::size_t>(vertex.getOrigin());
      }

      if (vertex.getOrigin() > Origin::Ket) {
        offset += ket_rank(tensor);
      }

      const std::size_t tensor_component_vertex = tensor_vertex + offset;

      assert(tensor_component_vertex < graph.vertex_labels.size());
      edges.push_back(std::make_pair(index_vertex, tensor_component_vertex));
    }
  }

  // also create vertices for pure proto indices
  for (const auto &index : pure_proto_indices_) {
    graph.vertex_labels.push_back(std::wstring(index.full_label()));
    graph.vertex_types.push_back(VertexType::Index);
    graph.vertex_colors.push_back(colorizer(index));

    const std::size_t index_vertex = graph.vertex_labels.size() - 1;

    index_vertices[index] = index_vertex;
  }

  // Add edges between proto index bundle vertices and all vertices of the
  // indices contained in that bundle i.e. if the bundle is {i_1,i_2}, the
  // bundle would be connected with vertices for i_1 and i_2
  for (const auto &[bundle, vertex] : proto_bundles) {
    for (const Index &idx : bundle) {
      auto it = index_vertices.find(idx);

      assert(it != index_vertices.end());
      if (it == index_vertices.end()) {
        std::abort();
      }

      edges.push_back(std::make_pair(it->second, vertex));
    }
  }

  assert(graph.vertex_labels.size() == graph.vertex_colors.size());
  assert(graph.vertex_labels.size() == graph.vertex_types.size());

  // Create the actual BLISS graph object
  graph.bliss_graph =
      std::make_unique<bliss::Graph>(graph.vertex_labels.size());

  for (const std::pair<std::size_t, std::size_t> &current_edge : edges) {
    graph.bliss_graph->add_edge(current_edge.first, current_edge.second);
  }

  for (const auto [vertex, color] :
       ranges::views::enumerate(graph.vertex_colors)) {
    graph.bliss_graph->change_color(vertex, color);
  }

  return graph;
}

void TensorNetworkV2::init_edges() {
  edges_.clear();
  ext_indices_.clear();
  pure_proto_indices_.clear();

  auto idx_insert = [this](const Index &idx, Vertex vertex) {
    if (Logger::instance().tensor_network) {
      std::wcout << "TensorNetworkV2::init_edges: idx=" << to_latex(idx)
                 << " attached to tensor " << vertex.getTerminalIndex() << " ("
                 << vertex.getOrigin() << ") at position "
                 << vertex.getIndexSlot()
                 << " (sym: " << to_wstring(vertex.getTerminalSymmetry()) << ")"
                 << std::endl;
    }

    auto it = std::find_if(edges_.begin(), edges_.end(),
                           FullLabelIndexLocator(idx.full_label()));
    if (it == edges_.end()) {
      edges_.emplace_back(std::move(vertex), idx);
    } else {
      it->connect_to(std::move(vertex));
    }
  };

  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    assert(tensors_[tensor_idx]);
    const AbstractTensor &tensor = *tensors_[tensor_idx];
    const Symmetry tensor_symm = symmetry(tensor);

    auto bra_indices = tensor._bra();
    for (std::size_t index_idx = 0; index_idx < bra_indices.size();
         ++index_idx) {
      idx_insert(bra_indices[index_idx],
                 Vertex(Origin::Bra, tensor_idx, index_idx, tensor_symm));
    }

    auto ket_indices = tensor._ket();
    for (std::size_t index_idx = 0; index_idx < ket_indices.size();
         ++index_idx) {
      idx_insert(ket_indices[index_idx],
                 Vertex(Origin::Ket, tensor_idx, index_idx, tensor_symm));
    }

    auto aux_indices = tensor._aux();
    for (std::size_t index_idx = 0; index_idx < aux_indices.size();
         ++index_idx) {
      // Note: for the time being we don't have a way of expressing
      // permutational symmetry of auxiliary indices so we just assume there is
      // no such symmetry
      idx_insert(aux_indices[index_idx],
                 Vertex(Origin::Aux, tensor_idx, index_idx, Symmetry::nonsymm));
    }
  }

  // extract external indices and all protoindices (since some external indices
  // may be pure protoindices)
  NamedIndexSet proto_indices;
  for (const Edge &current : edges_) {
    assert(current.vertex_count() > 0);
    // External index (== Edge only connected to a single vertex in the
    // network)
    if (current.vertex_count() == 1) {
      if (Logger::instance().tensor_network) {
        std::wcout << "idx " << to_latex(current.idx()) << " is external"
                   << std::endl;
      }

      const auto &[it, inserted] = ext_indices_.insert(current.idx());
      // only scenario where idx is already in ext_indices_ if it were a
      // protoindex of a previously inserted ext index ... check to ensure no
      // accidental duplicates
      if (!inserted) {
        assert(proto_indices.contains(current.idx()));
      }
    }

    // add proto indices to the grand list of proto indices
    for (auto &&proto_idx : current.idx().proto_indices()) {
      if (proto_idx.has_proto_indices())
        throw std::runtime_error(
            "TensorNetworkV2 does not support recursive protoindices");  // for
                                                                         // now
                                                                         // no
                                                                         // recursive
                                                                         // proto
                                                                         // indices
      proto_indices.insert(proto_idx);
    }
  }

  // now identify pure protoindices ...
  for (const Edge &current : edges_) {
    auto it = proto_indices.find(current.idx());
    if (it != proto_indices.end()) proto_indices.erase(it);
  }
  pure_proto_indices_ = std::move(proto_indices);
  if (Logger::instance().tensor_network) {
    for (auto &&idx : pure_proto_indices_) {
      std::wcout << "idx " << to_latex(idx) << " is pure protoindex"
                 << std::endl;
    }
  }

  // ... and add pure protoindices to the external indices
  ext_indices_.reserve(ext_indices_.size() + pure_proto_indices_.size());
  ext_indices_.insert(pure_proto_indices_.begin(), pure_proto_indices_.end());

  have_edges_ = true;
}

container::svector<std::pair<long, long>> TensorNetworkV2::factorize() {
  abort();  // not yet implemented
}

ExprPtr TensorNetworkV2::canonicalize_individual_tensor_blocks(
    const NamedIndexSet &named_indices) {
  return do_individual_canonicalization(
      TensorBlockCanonicalizer(named_indices));
}

ExprPtr TensorNetworkV2::canonicalize_individual_tensors(
    const NamedIndexSet &named_indices) {
  return do_individual_canonicalization(
      DefaultTensorCanonicalizer(named_indices));
}

ExprPtr TensorNetworkV2::do_individual_canonicalization(
    const TensorCanonicalizer &canonicalizer) {
  ExprPtr byproduct = ex<Constant>(1);

  for (auto &tensor : tensors_) {
    auto nondefault_canonizer_ptr =
        TensorCanonicalizer::nondefault_instance_ptr(tensor->_label());
    const TensorCanonicalizer &tensor_canonizer =
        nondefault_canonizer_ptr ? *nondefault_canonizer_ptr : canonicalizer;

    auto bp = canonicalizer.apply(*tensor);

    if (bp) {
      byproduct *= bp;
    }
  }

  return byproduct;
}

}  // namespace sequant
