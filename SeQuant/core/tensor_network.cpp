//
// Created by Eduard Valeyev on 2019-02-26.
//

#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
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
#include <SeQuant/core/wstring.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <type_traits>
#include <numeric>
#include <limits>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

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
    if (auxiliary_rank(lhs) != auxiliary_rank(rhs)) {
      return auxiliary_rank(lhs) < auxiliary_rank(rhs);
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
/// via their label, they are compared based on their tensor block (the spaces
/// of their indices). If this doesn't discriminate the tensors, they are
/// considered equal (note: explicit indexing is NOT compared)
template <typename CardinalLabels>
struct CanonicalTensorCompare {
  const CardinalLabels &labels;

  CanonicalTensorCompare(const CardinalLabels &labels) : labels(labels) {}

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
    // -> Discriminate via tensor block comparison
    TensorBlockCompare cmp;
    return cmp(lhs, rhs);
  }
};

TensorNetwork::Vertex::Vertex(Origin origin, std::size_t terminal_idx,
                              std::size_t index_slot, Symmetry terminal_symm)
    : origin(origin),
      terminal_idx(terminal_idx),
      index_slot(index_slot),
      terminal_symm(terminal_symm) {}

TensorNetwork::Origin TensorNetwork::Vertex::getOrigin() const {
  return origin;
}

std::size_t TensorNetwork::Vertex::getTerminalIndex() const {
  return terminal_idx;
}

std::size_t TensorNetwork::Vertex::getIndexSlot() const { return index_slot; }

Symmetry TensorNetwork::Vertex::getTerminalSymmetry() const {
  return terminal_symm;
}

bool TensorNetwork::Vertex::operator<(const Vertex &rhs) const {
  if (terminal_idx != rhs.terminal_idx) {
    return terminal_idx < rhs.terminal_idx;
  }

  // Both vertices belong to same tensor -> they must have same symmetry
  assert(terminal_symm == rhs.terminal_symm);

  // We only take the index slot into account for non-symmetric tensors
  if (terminal_symm == Symmetry::nonsymm && index_slot != rhs.index_slot) {
    return index_slot < rhs.index_slot;
  }

  // Note: The ordering of index groups must be consistent with the canonical
  // (Abstract)Tensor::operator<
  return origin < rhs.origin;
}

bool TensorNetwork::Vertex::operator==(const Vertex &rhs) const {
  // Slot position is only taken into account for non_symmetric tensors
  const std::size_t lhs_slot =
      (terminal_symm == Symmetry::nonsymm) * index_slot;
  const std::size_t rhs_slot =
      (rhs.terminal_symm == Symmetry::nonsymm) * rhs.index_slot;

  assert(terminal_idx != rhs.terminal_idx ||
         terminal_symm == rhs.terminal_symm);

  return terminal_idx == rhs.terminal_idx && lhs_slot == rhs_slot &&
         origin == rhs.origin;
}

std::size_t TensorNetwork::Graph::vertex_to_index_idx(
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

std::size_t TensorNetwork::Graph::vertex_to_tensor_idx(
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

template <typename ArrayLike, typename ReplacementMap>
void apply_index_replacements(ArrayLike &tensors,
                              const ReplacementMap &replacements) {
#ifndef NDEBUG
  // assert that tensors' indices are not tagged since going to tag indices
  {
    for (const auto &tensor : tensors) {
      assert(ranges::none_of(braket(*tensor), [](const Index &idx) {
        return idx.tag().has_value();
      }));
    }
  }
#endif

  bool pass_mutated = false;
  do {
    pass_mutated = false;
    for (auto &tensor : tensors) {
      pass_mutated |= transform_indices(*tensor, replacements);
    }
  } while (pass_mutated);  // transform till stops changing

  // untag transformed indices (if any)
  {
    for (auto &tensor : tensors) {
      reset_tags(*tensor);
    }
  }
}

void TensorNetwork::canonicalize_graph(const named_indices_t &named_indices) {
  if (Logger::get_instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize_graph: input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    std::wcout << std::endl;
  }

  if (!have_edges_) {
    init_edges();
  }
  // - canonize indices
  // - canonize tensors using canonical list of indices
  // Algorithm sketch:
  // - to canonize indices make a graph whose vertices are the indices as
  // well as the tensors and their terminals.
  //   - Indices with protoindices are connected to their protoindices,
  //   either directly or (if protoindices are symmetric) via a protoindex
  //   vertex.
  //   - Indices are colored by their space, which in general encodes also
  //   the space of the protoindices.
  //   - An anti/symmetric n-body tensor has 2 terminals, each connected to
  //   each other + to n index vertices.
  //   - A nonsymmetric n-body tensor has n terminals, each connected to 2
  //   indices and 1 tensor vertex which is connected to all n terminal
  //   indices.
  //   - tensor vertices are colored by the label+rank+symmetry of the
  //   tensor; terminal vertices are colored by the color of its tensor,
  //     with the color of symm/antisymm terminals augmented by the
  //     terminal's type (bra/ket).
  // - canonize the graph

  auto is_anonymous_index = [named_indices](const Index &idx) {
    return named_indices.find(idx) == named_indices.end();
  };

  // index factory to generate anonymous indices
  IndexFactory idxfac(is_anonymous_index, 1);

  // Clear any potential prior replacements
  idxrepl_.clear();

  // make the graph
  Graph graph = create_graph(&named_indices);
  // graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);

  // TODO: Add logging var for this
  if (false) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);
  }

  // canonize the graph
  bliss::Stats stats;
  graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
  const unsigned int *cl =
      graph.bliss_graph->canonical_form(stats, nullptr, nullptr);

  if (Logger::get_instance().canonicalize_dot) {
    std::wcout << "Canonicalization permutation:\n";
    for (std::size_t i = 0; i < graph.vertex_labels.size(); ++i) {
      std::wcout << i << " -> " << cl[i] << "\n";
    }
    std::wcout << "Canonicalized graph:\n";
    bliss::Graph *cgraph = graph.bliss_graph->permute(cl);
    cgraph->write_dot(std::wcout, {}, true);
    auto cvlabels = permute(graph.vertex_labels, cl);
    std::wcout << "with our labels:\n";
    cgraph->write_dot(std::wcout, cvlabels);
    delete cgraph;
  }

  // make anonymous index replacement list
  {
    // for each color make a replacement list for bringing the indices to
    // the canonical order
    container::set<size_t> colors;
    // maps color to the ordinals of the corresponding indices in edges_ and
    // their canonical ordinals collect colors and anonymous indices sorted by
    // colors
    container::multimap<size_t, std::pair<size_t, size_t>> color2idx;

    std::size_t idx_cnt = 0;
    auto edge_it = edges_.begin();
    for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
      if (graph.vertex_types[vertex] != VertexType::Index) {
        continue;
      }

      auto color = graph.vertex_colors[vertex];
      colors.insert(color);

      // We rely on the fact that the order of index vertices corresponds to
      // the order of corresponding edges
      assert(edge_it != edges_.end());
      const Index &idx = edge_it->idx();

      if (is_anonymous_index(idx)) {
        color2idx.emplace(color, std::make_pair(idx_cnt, cl[vertex]));
      }

      ++idx_cnt;
      ++edge_it;
    }

    // for each color sort anonymous indices by canonical order

    // canonically-ordered list of
    // {index ordinal in edges_, canonical ordinal}
    container::svector<std::pair<size_t, size_t>> idx_can;
    for (const std::size_t color : colors) {
      auto beg = color2idx.lower_bound(color);
      auto end = color2idx.upper_bound(color);

      const auto sz = std::distance(beg, end);

      if (sz > 1) {
        idx_can.resize(sz);

        size_t cnt = 0;
        for (auto it = beg; it != end; ++it, ++cnt) {
          idx_can[cnt] = it->second;
        }

        using std::begin;
        using std::end;
        std::sort(begin(idx_can), end(idx_can),
                  [](const std::pair<size_t, size_t> &a,
                     const std::pair<size_t, size_t> &b) {
                    return a.second < b.second;
                  });

        // make a replacement list by generating new indices in canonical
        // order
        for (auto [orig_idx, _] : idx_can) {
          assert(orig_idx < edges_.size());
          auto edge_it = edges_.begin();
          std::advance(edge_it, orig_idx);

          const auto &idx = edge_it->idx();

          idxrepl_.emplace(std::make_pair(idx, idxfac.make(idx)));
        }
      } else if (sz == 1) {
        // no need for resorting of colors with 1 index only, but still need
        // to replace the index
        const auto edge_it = edges_.begin() + beg->second.first;
        const auto &idx = edge_it->idx();
        idxrepl_.emplace(std::make_pair(idx, idxfac.make(idx)));
      }
      // sz == 0 is possible since some colors in colors refer to tensors
    }
  }  // index repl

  // Bring tensors into canonical order, but ensure to respect commutativity!

  std::vector<std::size_t> tensor_indices;
  tensor_indices.resize(tensors_.size());
  std::iota(tensor_indices.begin(), tensor_indices.end(), 0);

  // TODO: initialize this while iterating over vertices for index (colors)
  container::map<std::size_t, std::size_t> tensor_idx_to_vertex;
  for (std::size_t tensor_idx = 0, vertex = 0;
       vertex < graph.vertex_types.size(); ++vertex) {
    if (graph.vertex_types[vertex] != VertexType::TensorCore) {
      continue;
    }
    tensor_idx_to_vertex[tensor_idx] = vertex;
    tensor_idx++;
  }

  const auto tensor_sorter = [this, &cl, &tensor_idx_to_vertex](
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
    return cl[lhs_vertex] < cl[rhs_vertex];
  };

  std::stable_sort(tensor_indices.begin(), tensor_indices.end(), tensor_sorter);

  assert(tensor_indices.size() == tensors_.size());

  decltype(tensors_) tensors_canonized(tensors_.size());
  for (std::size_t i = 0; i < tensors_.size(); ++i) {
    tensors_canonized[i] = tensors_[tensor_indices[i]];
  }

  tensors_ = std::move(tensors_canonized);

  // Apply the canonical reindexing
  apply_index_replacements(tensors_, idxrepl_);

  have_edges_ = false;

  if (Logger::get_instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize_graph: tensors after "
                  "canonicalization\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
  }
}

ExprPtr TensorNetwork::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels, bool fast,
    const named_indices_t *named_indices_ptr) {
  // initialize named_indices by default to all external indices
  const auto &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  if (Logger::get_instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize(" << (fast ? "fast" : "slow")
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

  if (!fast) {
    // The graph-based canonization is required in call cases in which there are
    // indistinguishable tensors present in the expression. Their order and
    // indexing can only be determined via this rigorous canonization.
    canonicalize_graph(named_indices);
  }

  // Ensure each individual tensor is canonical with the current indexing in
  // order to properly be able to identify tensor blocks
  ExprPtr byproduct = canonicalize_individual_tensors(named_indices);

  const CanonicalTensorCompare<decltype(cardinal_tensor_labels)> tensor_sorter(
      cardinal_tensor_labels);

  std::stable_sort(tensors_.begin(), tensors_.end(), tensor_sorter);

  init_edges();

  if (Logger::get_instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize(" << (fast ? "fast" : "slow")
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

  // Sort edges based on the order of the tensors they connect (instead of based
  // on their (full) label)
  container::svector<Edge> resorted_edges(edges_.begin(), edges_.end());
  std::stable_sort(resorted_edges.begin(), resorted_edges.end(),
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

  idxrepl_.clear();

  // Use the new order of edges as the canonical order of indices and relabel
  // accordingly (but only anonymous indices, of course)
  for (std::size_t i = named_indices.size(); i < resorted_edges.size(); ++i) {
    const Index &index = resorted_edges[i].idx();
    assert(is_anonymous_index(index));
    Index replacement = idxfac.make(index);
    idxrepl_.emplace(std::make_pair(index, replacement));
  }

  // Done computing canonical index replacement list

  if (Logger::get_instance().canonicalize) {
    for (const auto &idxpair : idxrepl_) {
      std::wcout << "TensorNetwork::canonicalize(" << (fast ? "fast" : "slow")
                 << "): replacing " << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  apply_index_replacements(tensors_, idxrepl_);

  byproduct *= canonicalize_individual_tensors(named_indices);

  // We assume that re-indexing did not change the canonical order of tensors
  // If it did, then most likely the (Abstract)Tensor::operator< compares index
  // groups in a different order than Vertex::operator<
  assert(std::is_sorted(tensors_.begin(), tensors_.end(), tensor_sorter));

  have_edges_ = false;

  assert(byproduct->is<Constant>());
  return (byproduct->as<Constant>().value() == 1) ? nullptr : byproduct;
}

TensorNetwork::Graph TensorNetwork::create_graph(
    const named_indices_t *named_indices_ptr) const {
  assert(have_edges_);

  // initialize named_indices by default to all external indices
  const named_indices_t &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  // Colors in the range [ 0, 3 * max_rank + named_indices.size() ) are
  // reserved: Up to max_rank colors can be used for bra indices Up to
  // max_rank colors can be used for ket indices Up to max_rank colors can be
  // used for auxiliary indices Every named index is identified by a unique
  // color
  constexpr std::size_t named_idx_color_start = 3 * max_rank;
  const std::size_t max_reserved_color =
      named_idx_color_start + named_indices.size() - 1;

  // core, bra, ket, auxiliary
  constexpr std::size_t num_tensor_components = 4;

  // results
  Graph graph;
  // We know that at the very least all indices and all tensors will yield
  // vertex representations
  std::size_t vertex_count_estimate =
      edges_.size() + num_tensor_components * tensors_.size();
  graph.vertex_labels.reserve(vertex_count_estimate);
  graph.vertex_colors.reserve(vertex_count_estimate);
  graph.vertex_types.reserve(vertex_count_estimate);

  using proto_bundle_t =
      std::decay_t<decltype(std::declval<const Index &>().proto_indices())>;
  container::map<proto_bundle_t, std::size_t> proto_bundles;

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
    std::wstring_view tensor_label = label(tensor);
    graph.vertex_labels.emplace_back(tensor_label);
    graph.vertex_types.emplace_back(VertexType::TensorCore);
    const std::size_t tensor_color =
        hash::value(tensor_label) + max_reserved_color;
    assert(tensor_color > max_reserved_color);
    graph.vertex_colors.push_back(tensor_color);

    const std::size_t tensor_vertex = graph.vertex_labels.size() - 1;
    tensor_vertices.insert(std::make_pair(tensor_idx, tensor_vertex));

    // Create vertices to group indices
    const Symmetry tensor_sym = symmetry(tensor);
    if (tensor_sym == Symmetry::nonsymm) {
      // Create separate vertices for every index
      assert(bra_rank(tensor) <= max_rank);
      for (std::size_t i = 0; i < bra_rank(tensor); ++i) {
        graph.vertex_labels.emplace_back(L"bra_" + std::to_wstring(i + 1));
        graph.vertex_types.push_back(VertexType::TensorBra);
        graph.vertex_colors.push_back(i);
      }
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));

      assert(ket_rank(tensor) <= max_rank);
      for (std::size_t i = 0; i < ket_rank(tensor); ++i) {
        graph.vertex_labels.emplace_back(L"ket_" + std::to_wstring(i + 1));
        graph.vertex_types.push_back(VertexType::TensorKet);
        graph.vertex_colors.push_back(max_rank + i);
      }
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
    } else {
      // Shared set of bra/ket vertices for all indices
      std::wstring suffix = tensor_sym == Symmetry::symm ? L"_s" : L"_a";

      const std::size_t bra_color = 0;
      graph.vertex_labels.push_back(L"bra" + suffix);
      graph.vertex_types.push_back(VertexType::TensorBra);
      graph.vertex_colors.push_back(bra_color);
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));

      // TODO: figure out how to handle BraKetSymmetry::conjugate
      const std::size_t ket_color =
          braket_symmetry(tensor) == BraKetSymmetry::symm
              ? bra_color
              : bra_color + max_rank;
      graph.vertex_labels.push_back(L"ket" + suffix);
      graph.vertex_types.push_back(VertexType::TensorKet);
      graph.vertex_colors.push_back(ket_color);
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
    }

    // TODO: handle aux indices permutation symmetries
    // for now, auxiliary indices are considered to always be asymmetric
    for (std::size_t i = 0; i < auxiliary_rank(tensor); ++i) {
      graph.vertex_labels.emplace_back(L"aux_" + std::to_wstring(i + 1));
      graph.vertex_types.push_back(VertexType::TensorAux);
      graph.vertex_colors.push_back(2 * max_rank + i);
      edges.push_back(
          std::make_pair(tensor_vertex, graph.vertex_labels.size() - 1));
    }

    // TODO: Add vertex for total or braket grouping?
  }

  // Now add all indices (edges) to the graph
  for (const Edge &current_edge : edges_) {
    const Index &index = current_edge.idx();
    graph.vertex_labels.push_back(index.to_latex());
    graph.vertex_types.push_back(VertexType::Index);

    // Assign index color
    std::size_t idx_color;
    auto named_idx_iter = named_indices.find(index);
    if (named_idx_iter == named_indices.end()) {
      // This is an anonymous index
      idx_color = index.color();
      assert(idx_color > max_reserved_color);
    } else {
      idx_color = static_cast<decltype(idx_color)>(
          std::distance(named_indices.begin(), named_idx_iter));
      idx_color += named_idx_color_start;
    }
    graph.vertex_colors.push_back(idx_color);

    const std::size_t index_vertex = graph.vertex_labels.size() - 1;

    // Handle proto indices
    if (index.has_proto_indices()) {
      // For now we assume that all proto indices are symmetric
      assert(index.symmetric_proto_indices());

      std::size_t proto_vertex;
      if (auto it = proto_bundles.find(index.proto_indices());
          it != proto_bundles.end()) {
        proto_vertex = it->second;
      } else {
        // Create a new vertex for this bundle of proto indices
        std::wstring spbundle_label = L"{";
        for (const Index &proto : index.proto_indices()) {
          spbundle_label += proto.to_latex();
        }
        spbundle_label += L"}";
        graph.vertex_labels.push_back(std::move(spbundle_label));
        graph.vertex_types.push_back(VertexType::SPBundle);
        const std::size_t bundle_color =
            Index::proto_indices_color(index.proto_indices()) +
            max_reserved_color;
        assert(bundle_color);
        graph.vertex_colors.push_back(bundle_color);

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
      // TODO: if we re-introduce a braket-like vertex, we have to add 1 to
      // the tensor vertex
      static_assert(static_cast<int>(Origin::Bra) == 1);
      static_assert(static_cast<int>(Origin::Ket) == 2);
      static_assert(static_cast<int>(Origin::Aux) == 3);
      const bool tensor_is_nonsymm =
          vertex.getTerminalSymmetry() == Symmetry::nonsymm;
      // Determine the index of the vertex for the tensor component we want to
      // connect the current index to. For (anti)symmetric tensors there
      // exists only a single bra, ket or aux vertex per tensor and thus the
      // origin of the index is the only factor to account for. Non-symmetric
      // tensors on the other hand have a bra, ket or aux vertex for each
      // individual index and thus the index's position within the bra, ket or
      // aux group of the tensor has to be accounted for in order to selected
      // the correct vertex to connect to. N.B. conversion from bool to int is
      // always: true -> 1, false -> 0
      const std::size_t tensor_component_vertex =
          tensor_vertex + static_cast<int>(vertex.getOrigin()) +
          tensor_is_nonsymm * vertex.getIndexSlot() *
              (num_tensor_components - /* core */ 1);
      assert(tensor_component_vertex < graph.vertex_labels.size());
      edges.push_back(std::make_pair(index_vertex, tensor_component_vertex));
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

  // compress vertex colors to 32 bits, as required by Bliss, by hashing
  for (std::size_t vertex = 0; vertex < graph.vertex_colors.size(); ++vertex) {
    auto color = graph.vertex_colors[vertex];
    static_assert(sizeof(color) == 8);

    color = (~color) + (color << 18);  // color = (color << 18) - color - 1;
    color = color ^ (color >> 31);
    color = color * 21;  // color = (color + (color << 2)) + (color << 4);
    color = color ^ (color >> 11);
    color = color + (color << 6);
    color = color ^ (color >> 22);

    graph.bliss_graph->change_color(vertex, static_cast<int>(color));
  }

  return graph;
}

void TensorNetwork::init_edges() {
  edges_.clear();
  ext_indices_.clear();

  auto idx_insert = [this](const Index &idx, Vertex vertex) {
    if (Logger::get_instance().tensor_network) {
      std::wcout << "TensorNetwork::init_edges: idx=" << to_latex(idx)
                 << " attached to tensor " << vertex.getTerminalIndex() << " ("
                 << vertex.getOrigin() << ") at position "
                 << vertex.getIndexSlot()
                 << " (sym: " << to_wstring(vertex.getTerminalSymmetry()) << ")"
                 << std::endl;
    }

    auto it = edges_.find(idx.full_label());
    if (it == edges_.end()) {
      edges_.emplace(Edge(std::move(vertex), idx));
    } else {
      it->connect_to(std::move(vertex));
    }
  };

  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    assert(tensors_[tensor_idx]);
    const AbstractTensor &tensor = *tensors_[tensor_idx];
    const Symmetry tensor_symm = symmetry(tensor);

    auto bra_indices = bra(tensor);
    for (std::size_t index_idx = 0; index_idx < bra_indices.size();
         ++index_idx) {
      idx_insert(bra_indices[index_idx],
                 Vertex(Origin::Bra, tensor_idx, index_idx, tensor_symm));
    }

    auto ket_indices = ket(tensor);
    for (std::size_t index_idx = 0; index_idx < ket_indices.size();
         ++index_idx) {
      idx_insert(ket_indices[index_idx],
                 Vertex(Origin::Ket, tensor_idx, index_idx, tensor_symm));
    }

    auto aux_indices = auxiliary(tensor);
    for (std::size_t index_idx = 0; index_idx < aux_indices.size();
         ++index_idx) {
      idx_insert(aux_indices[index_idx],
                 Vertex(Origin::Aux, tensor_idx, index_idx, tensor_symm));
    }
  }

  // extract external indices
  for (const Edge &current : edges_) {
    assert(current.vertex_count() > 0);
    if (current.vertex_count() == 1) {
      // External index (== Edge only connected to a single vertex in the
      // network)
      if (Logger::get_instance().tensor_network) {
        std::wcout << "idx " << to_latex(current.idx()) << " is external"
                   << std::endl;
      }

      bool inserted = ext_indices_.insert(current.idx()).second;
      assert(inserted);
    }
  }

  have_edges_ = true;
}

container::svector<std::pair<long, long>> TensorNetwork::factorize() {
  abort();  // not yet implemented
}

ExprPtr TensorNetwork::canonicalize_individual_tensors(
    const named_indices_t &named_indices) {
  ExprPtr byproduct = ex<Constant>(1);

  // override the default canonicalizer
  DefaultTensorCanonicalizer default_tensor_canonizer(named_indices);
  for (auto &tensor : tensors_) {
    auto nondefault_canonizer_ptr =
        TensorCanonicalizer::nondefault_instance_ptr(tensor->_label());
    TensorCanonicalizer *tensor_canonizer = nondefault_canonizer_ptr
                                                ? nondefault_canonizer_ptr.get()
                                                : &default_tensor_canonizer;

    auto bp = tensor_canonizer->apply(*tensor);

    if (bp) {
      byproduct *= bp;
    }
  }

  return byproduct;
}

}  // namespace sequant
