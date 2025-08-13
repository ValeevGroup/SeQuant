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
#include <SeQuant/core/tensor_network/utils.hpp>
#include <SeQuant/core/tensor_network/vertex_painter.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/utility/swap.hpp>
#include <SeQuant/core/utility/tuple.hpp>
#include <SeQuant/core/wstring.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <type_traits>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

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

  const auto is_anonymous_index = [&named_indices](const Index &idx) {
    return named_indices.find(idx) == named_indices.end();
  };

  // index factory to generate anonymous indices
  IndexFactory idxfac(is_anonymous_index, 1);

  // make the graph
  Graph graph = create_graph(
      {.named_indices = &named_indices,
       .make_labels = Logger::instance().canonicalize_input_graph ||
                      Logger::instance().canonicalize_dot,
       .make_texlabels = Logger::instance().canonicalize_input_graph ||
                         Logger::instance().canonicalize_dot});
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
    cgraph->write_dot(std::wcout, {}, {}, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    std::wcout << "with our labels:\n";
    cgraph->write_dot(std::wcout, cvlabels);
    delete cgraph;
  }

  // maps tensor ordinal -> input vertex ordinal
  std::vector<std::size_t> tensor_idx_to_vertex;
  tensor_idx_to_vertex.reserve(tensors_.size());
  std::size_t tensor_idx = 0;
  // for nonsymmetric tensors only: maps tensor ordinal -> canonical order of
  // its columns'/particles' vertex ordinals
  container::map<std::size_t, container::svector<std::size_t, 3>>
      tensor_idx_to_particle_order;
  std::vector<std::size_t> index_idx_to_vertex;
  index_idx_to_vertex.reserve(edges_.size() + pure_proto_indices_.size());

  for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
    switch (graph.vertex_types[vertex]) {
      case VertexType::Index:
        index_idx_to_vertex.emplace_back(index_idx_to_vertex.size()) = vertex;
        break;
      case VertexType::TensorBraKet: {
        assert(tensor_idx > 0);
        const std::size_t base_tensor_idx = tensor_idx - 1;
        assert(symmetry(*tensors_.at(base_tensor_idx)) == Symmetry::nonsymm);
        tensor_idx_to_particle_order[base_tensor_idx].push_back(
            canonize_perm[vertex]);
        break;
      }
      case VertexType::TensorCore:
        assert(tensor_idx_to_vertex.size() == tensor_idx);
        tensor_idx_to_vertex.emplace_back(tensor_idx_to_vertex.size()) = vertex;
        tensor_idx++;
        break;
      case VertexType::TensorBra:
      case VertexType::TensorKet:
      case VertexType::TensorAux:
      case VertexType::TensorBraBundle:
      case VertexType::TensorKetBundle:
      case VertexType::TensorAuxBundle:
      case VertexType::SPBundle:
        break;
    }
  }

  assert(index_idx_to_vertex.size() ==
         edges_.size() + pure_proto_indices_.size());
  assert(tensor_idx_to_vertex.size() == tensors_.size());
  assert(tensor_idx_to_particle_order.size() <= tensors_.size());

  // sort_then_replace_by_ordinals(index_order);
  for (auto &current : tensor_idx_to_particle_order) {
    sort_then_replace_by_ordinals(current.second);
  }

  container::map<Index, Index> idxrepl;
  auto idxrepl_emplace = [&idxrepl](auto &&from, auto &&to) {
    if (from != to) idxrepl.emplace(std::move(from), std::move(to));
  };

  // Sort edges so that their order corresponds to the order of indices in the
  // canonical graph

  // Use this ordering to relabel anonymous indices
  const auto index_sorter = [&index_idx_to_vertex, &canonize_perm](
                                std::size_t lhs_idx, std::size_t rhs_idx) {
    assert(lhs_idx < index_idx_to_vertex.size());
    const std::size_t lhs_vertex = index_idx_to_vertex[lhs_idx];
    assert(rhs_idx < index_idx_to_vertex.size());
    const std::size_t rhs_vertex = index_idx_to_vertex[rhs_idx];

    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  sort_via_ordinals<OrderType::StrictWeak>(edges_, index_sorter);

  for (const Edge &current : edges_) {
    const Index &idx = current.idx();

    const auto is_named = current.vertex_count() == 1;
    if (is_named) continue;

    idxrepl_emplace(idx, idxfac.make(idx));
  }

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetworkV2::canonicalize_graph: replacing "
                 << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  // The tensor reordering and index relabeling will make edges_ invalid
  edges_.clear();
  have_edges_ = false;

  apply_index_replacements(tensors_, idxrepl, true);

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

      idxrepl_emplace(bra_indices[col], bra_indices[particle_order[col]]);
      idxrepl_emplace(ket_indices[col], ket_indices[particle_order[col]]);
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
      apply_index_replacements(tensor, idxrepl, false);
    }
  }

  // Bring tensors into canonical order (analogously to how we reordered
  // indices); elements that do not commute are equivalent, due to this
  // this specifies a non-strict weak order
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

  sort_via_ordinals<OrderType::Weak>(tensors_, tensor_sorter);

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
  if (Logger::instance().canonicalize) {
    std::wcout << "named_indices = ";
    ranges::for_each(named_indices,
                     [](auto &&i) { std::wcout << i.full_label() << L" "; });
  }

  if (!fast) {
    // The graph-based canonization is required in all cases in which there are
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
    if (index != replacement) idxrepl.emplace(index, std::move(replacement));
  }

  // Done computing canonical index replacement list
  // reset edges since renamings will make them obsolete
  edges_.clear();
  have_edges_ = false;

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetworkV2::canonicalize(" << (fast ? "fast" : "slow")
                 << "): replacing " << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  apply_index_replacements(tensors_, idxrepl, true);

  byproduct *= canonicalize_individual_tensors(named_indices);

  // We assume that re-indexing did not change the canonical order of tensors
  assert(std::is_sorted(tensors_.begin(), tensors_.end(), tensor_sorter));
  // However, in order to produce the most aesthetically pleasing result, we now
  // reorder tensors based on the regular AbstractTensor::operator<, which takes
  // the explicit index labelling of tensors into account.
  tensor_sorter.set_blocks_only(false);
  std::stable_sort(tensors_.begin(), tensors_.end(), tensor_sorter);

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
  Graph graph = create_graph(
      {.named_indices = &named_indices,
       .distinct_named_indices = false,
       .make_labels = Logger::instance().canonicalize_input_graph ||
                      Logger::instance().canonicalize_dot,
       .make_texlabels = Logger::instance().canonicalize_input_graph ||
                         Logger::instance().canonicalize_dot,
       .make_idx_to_vertex = true});
  const auto &idx_to_vertex = graph.idx_to_vertex;
  // graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);

  if (Logger::instance().canonicalize_input_graph) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels,
                                 graph.vertex_texlabels);
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
    metadata.graph->write_dot(std::wcout, {}, {}, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    auto cvtexlabels = permute(graph.vertex_texlabels, canonize_perm);
    std::wcout << "with our labels:\n";
    metadata.graph->write_dot(std::wcout, cvlabels, cvtexlabels);
  }

  // produce canonical list of named indices
  {
    using ord_cord_it_t =
        std::tuple<size_t, size_t, NamedIndexSet::const_iterator>;
    using cord_set_t = container::set<ord_cord_it_t, detail::tuple_less<1>>;

    auto grand_index_list = ranges::views::concat(
        edges_ | ranges::views::transform(&Edge::idx), pure_proto_indices_);

    // for each named index type (as defined by named_index_compare) maps its
    // ptr in grand_index_list to its ordinal in grand_index_list + canonical
    // ordinal + its iterator in metadata.named_indices
    container::map<std::pair<const Index *, IndexSlotType>, cord_set_t,
                   decltype(named_index_compare)>
        idx2cord(named_index_compare);

    // collect named indices and sort them on the fly
    auto grand_index_list_end = grand_index_list.end();
    for (auto [idx_ord, idx] : ranges::views::enumerate(grand_index_list)) {
      if (!is_named_index(idx)) {
        continue;
      }

      const auto named_indices_it = metadata.named_indices.find(idx);
      assert(named_indices_it != metadata.named_indices.end());
      const auto vertex_ord = idx_to_vertex.at(*named_indices_it);

      // find the entry for this index type
      IndexSlotType slot_type;
      if (idx_ord < edges_.size()) {
        auto edge_it = edges_.begin();
        std::advance(edge_it, idx_ord);
        // there are 2 possibilities: its index edge is disconnected or
        // connected ... the latter would only occur if this index is named
        // due to also being a protoindex on one of the named indices!
        if (edge_it->vertex_count() == 1) {
          if (edge_it->first_vertex().getOrigin() == Origin::Aux)
            slot_type = IndexSlotType::TensorAux;
          else if (edge_it->first_vertex().getOrigin() == Origin::Bra)
            slot_type = IndexSlotType::TensorBra;
          else {
            assert(edge_it->first_vertex().getOrigin() == Origin::Ket);
            slot_type = IndexSlotType::TensorKet;
          }
        } else {  // if
          assert(edge_it->vertex_count() == 2);
          slot_type = IndexSlotType::SPBundle;
        }
      } else
        slot_type = IndexSlotType::SPBundle;
      const auto idxptr_slottype = std::make_pair(&idx, slot_type);
      auto it = idx2cord.find(idxptr_slottype);

      if (it == idx2cord.end()) {
        bool inserted;
        std::tie(it, inserted) = idx2cord.emplace(
            idxptr_slottype, cord_set_t(cord_set_t::key_compare{}));
        assert(inserted);
      }

      it->second.emplace(idx_ord, canonize_perm[vertex_ord], named_indices_it);
    }

    // save the result
    for (auto &[idxptr_slottype, cord_set] : idx2cord) {
      for (auto &[idx_ord, idx_ord_can, named_indices_it] : cord_set) {
        metadata.named_indices_canonical.emplace_back(named_indices_it);
      }
    }
    metadata.named_index_compare = std::move(named_index_compare);

  }  // named indices resort to canonical order

  // - For each bra/ket bundle canonical order of slots is the lexicographic
  //   order of the canonicalized vertices representing the contained indices.
  // - Reordering indices into this canonical order incurs a phase change if the
  //   index bundle is antisymmetric.
  // - Determine this phase change by determining the parity of index
  //   permutations required to arrive at canonical form
  metadata.phase = 1;
  container::svector<SwapCountable<std::size_t>> vertices;
  for (const AbstractTensor &tensor : tensors_ | ranges::views::indirect) {
    if (symmetry(tensor) != Symmetry::antisymm) {
      // Only antisymmetric tensors (or rather: their indices) can incur a phase
      // change due to index permutation
      continue;
    }

    // Note that the current assumption is that auxiliary indices don't have
    // permutational symmetry, let alone being antisymmetric. Hence, we don't
    // have to include them in the iteration.
    // Note2: have to create dedicated container to hold ranges as an
    // initializer list will only return const entries upon iteration and one
    // can't iterate over const ranges.
    std::vector index_groups = {tensor._bra(), tensor._ket()};
    for (auto &indices : index_groups) {
      using ranges::size;
      std::size_t n_indices = size(indices);

      if (n_indices < 2) {
        // If there are < 2 indices, no two indices could have been swapped
        continue;
      }

      vertices.clear();
      vertices.reserve(n_indices);

      for (const Index &idx : indices) {
        const std::size_t vertex = idx_to_vertex.at(idx);
        vertices.emplace_back(canonize_perm[vertex]);
      }

      reset_ts_swap_counter<std::size_t>();
      bubble_sort(vertices.begin(), vertices.end());
      if (!ts_swap_counter_is_even<std::size_t>()) {
        // Performed an uneven amount of pairwise exchanges -> this incurs a
        // phase change
        metadata.phase *= -1;
      }
    }
  }

  return metadata;
}

TensorNetworkV2::Graph TensorNetworkV2::create_graph(
    const CreateGraphOptions &options) const {
  assert(have_edges_);

  // initialize named_indices by default to all external indices
  const NamedIndexSet &named_indices = options.named_indices == nullptr
                                           ? this->ext_indices()
                                           : *(options.named_indices);

  VertexPainter colorizer(named_indices, options.distinct_named_indices);

  // core, bra, ket, auxiliary and optionally (for non-symmetric tensors) a
  // particle vertex
  constexpr std::size_t num_tensor_components = 5;

  // results
  Graph graph;
  std::size_t nvertex = 0;
  // We know that at the very least all indices and all tensors will yield
  // vertex representations
  std::size_t vertex_count_estimate = edges_.size() +
                                      pure_proto_indices_.size() +
                                      num_tensor_components * tensors_.size();
  if (options.make_labels) graph.vertex_labels.reserve(vertex_count_estimate);
  if (options.make_texlabels)
    graph.vertex_texlabels.reserve(vertex_count_estimate);
  graph.vertex_colors.reserve(vertex_count_estimate);
  graph.vertex_types.reserve(vertex_count_estimate);

  container::svector<std::pair<ProtoBundle, std::size_t>> proto_bundles;

  // Mapping from the i-th tensor in tensors_ to the ID of the corresponding
  // vertex
  static constexpr std::size_t uninitialized_vertex =
      std::numeric_limits<std::size_t>::max();
  container::svector<std::size_t> tensor_vertices;
  tensor_vertices.resize(tensors_.size(), uninitialized_vertex);

  container::vector<std::pair<std::size_t, std::size_t>> edges;
  edges.reserve(edges_.size() + tensors_.size());

  // Add vertices for tensors
  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    assert(tensor_idx < tensor_vertices.size());
    assert(tensor_vertices[tensor_idx] == uninitialized_vertex);
    assert(tensors_.at(tensor_idx));
    const AbstractTensor &tensor = *tensors_.at(tensor_idx);

    // Tensor core
    const auto tlabel = label(tensor);
    ++nvertex;
    if (options.make_labels) graph.vertex_labels.emplace_back(tlabel);
    if (options.make_texlabels)
      graph.vertex_texlabels.emplace_back(L"$" + utf_to_latex(tlabel) + L"$");
    graph.vertex_types.emplace_back(VertexType::TensorCore);
    graph.vertex_colors.push_back(colorizer(tensor));

    const std::size_t tensor_vertex = nvertex - 1;
    tensor_vertices[tensor_idx] = tensor_vertex;

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
        ++nvertex;
        if (options.make_labels)
          graph.vertex_labels.emplace_back(L"p_" + std::to_wstring(i + 1));
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(std::nullopt);
        graph.vertex_types.push_back(VertexType::TensorBraKet);
        // Particles are indistinguishable -> always use same ID
        graph.vertex_colors.push_back(colorizer(ParticleGroup{0}));
        edges.push_back(std::make_pair(tensor_vertex, nvertex - 1));
      }

      for (std::size_t i = 0; i < bra_rank(tensor); ++i) {
        const bool is_unpaired_idx = i >= num_particle_vertices;
        const bool color_idx = is_unpaired_idx || !is_part_symm;

        ++nvertex;
        if (options.make_labels)
          graph.vertex_labels.emplace_back(L"bra_" + std::to_wstring(i + 1));
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(std::nullopt);
        graph.vertex_types.push_back(VertexType::TensorBra);
        graph.vertex_colors.push_back(colorizer(BraGroup{color_idx ? i : 0}));

        const std::size_t connect_vertex =
            tensor_vertex + (is_unpaired_idx ? 0 : (i + 1));
        edges.push_back(std::make_pair(connect_vertex, nvertex - 1));
      }

      for (std::size_t i = 0; i < ket_rank(tensor); ++i) {
        const bool is_unpaired_idx = i >= num_particle_vertices;
        const bool color_idx = is_unpaired_idx || !is_part_symm;

        ++nvertex;
        if (options.make_labels)
          graph.vertex_labels.emplace_back(L"ket_" + std::to_wstring(i + 1));
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(std::nullopt);
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
        edges.push_back(std::make_pair(connect_vertex, nvertex - 1));
      }
    } else {
      // Shared set of bra/ket vertices for all indices
      std::wstring suffix = tensor_sym == Symmetry::symm ? L"_s" : L"_a";

      ++nvertex;
      if (options.make_labels) graph.vertex_labels.push_back(L"bra" + suffix);
      if (options.make_texlabels)
        graph.vertex_texlabels.emplace_back(std::nullopt);
      graph.vertex_types.push_back(VertexType::TensorBra);
      graph.vertex_colors.push_back(colorizer(BraGroup{0}));
      edges.push_back(std::make_pair(tensor_vertex, nvertex - 1));

      ++nvertex;
      if (options.make_labels) graph.vertex_labels.push_back(L"ket" + suffix);
      if (options.make_texlabels)
        graph.vertex_texlabels.emplace_back(std::nullopt);
      graph.vertex_types.push_back(VertexType::TensorKet);
      // TODO: figure out how to handle BraKetSymmetry::conjugate
      if (braket_symmetry(tensor) == BraKetSymmetry::symm) {
        // Use BraGroup for kets as well as they should be indistinguishable
        graph.vertex_colors.push_back(colorizer(BraGroup{0}));
      } else {
        graph.vertex_colors.push_back(colorizer(KetGroup{0}));
      }
      edges.push_back(std::make_pair(tensor_vertex, nvertex - 1));
    }

    // TODO: handle aux indices permutation symmetries once they are supported
    // for now, auxiliary indices are considered to always be asymmetric
    for (std::size_t i = 0; i < aux_rank(tensor); ++i) {
      ++nvertex;
      if (options.make_labels)
        graph.vertex_labels.emplace_back(L"aux_" + std::to_wstring(i + 1));
      if (options.make_texlabels)
        graph.vertex_texlabels.emplace_back(std::nullopt);
      graph.vertex_types.push_back(VertexType::TensorAux);
      graph.vertex_colors.push_back(colorizer(AuxGroup{i}));
      edges.push_back(std::make_pair(tensor_vertex, nvertex - 1));
    }
  }

  // Now add all indices (edges_ + pure_proto_indices_) to the graph
  container::vector<std::size_t> index_vertices;
  index_vertices.resize(edges_.size() + pure_proto_indices_.size(),
                        uninitialized_vertex);

  for (std::size_t i = 0; i < edges_.size(); ++i) {
    const Edge &current_edge = edges_[i];

    const Index &index = current_edge.idx();
    ++nvertex;
    if (options.make_labels)
      graph.vertex_labels.push_back(std::wstring(index.full_label()));
    using namespace std::string_literals;
    if (options.make_texlabels)
      graph.vertex_texlabels.emplace_back(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.push_back(VertexType::Index);
    graph.vertex_colors.push_back(colorizer(index));

    const std::size_t index_vertex = nvertex - 1;

    assert(index_vertices.at(i) == uninitialized_vertex);
    index_vertices[i] = index_vertex;

    // Handle proto indices
    if (index.has_proto_indices()) {
      // For now we assume that all proto indices are symmetric
      assert(index.symmetric_proto_indices());

      std::size_t proto_vertex;
      if (auto it =
              std::ranges::find(proto_bundles, index.proto_indices(),
                                &decltype(proto_bundles)::value_type::first);
          it != proto_bundles.end()) {
        proto_vertex = it->second;
      } else {
        // Create a new vertex for this bundle of proto indices
        ++nvertex;
        if (options.make_labels) {
          using namespace std::literals;
          std::wstring spbundle_label =
              L"<" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.full_label(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L">";
          graph.vertex_labels.push_back(std::move(spbundle_label));
        }
        if (options.make_texlabels) {
          using namespace std::literals;
          std::wstring spbundle_texlabel =
              L"$\\langle" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.to_latex(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L"\\rangle$";
          graph.vertex_texlabels.push_back(std::move(spbundle_texlabel));
        }
        graph.vertex_types.push_back(VertexType::SPBundle);
        graph.vertex_colors.push_back(colorizer(index.proto_indices()));

        proto_vertex = nvertex - 1;
        proto_bundles.emplace_back(index.proto_indices(), proto_vertex);
      }

      edges.push_back(std::make_pair(index_vertex, proto_vertex));
    }

    // Connect index to the tensor(s) it is connected to
    for (std::size_t i = 0; i < current_edge.vertex_count(); ++i) {
      assert(i <= 1);
      const Vertex &vertex =
          i == 0 ? current_edge.first_vertex() : current_edge.second_vertex();

      assert(vertex.getTerminalIndex() < tensor_vertices.size());
      assert(tensor_vertices[vertex.getTerminalIndex()] !=
             uninitialized_vertex);
      const std::size_t tensor_vertex =
          tensor_vertices[vertex.getTerminalIndex()];

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

      assert(tensor_component_vertex < nvertex);
      edges.push_back(std::make_pair(index_vertex, tensor_component_vertex));
    }
  }

  // also create vertices for pure proto indices
  for (const auto &[i, index] : ranges::views::enumerate(pure_proto_indices_)) {
    ++nvertex;
    if (options.make_labels)
      graph.vertex_labels.push_back(std::wstring(index.full_label()));
    using namespace std::string_literals;
    if (options.make_texlabels)
      graph.vertex_texlabels.push_back(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.push_back(VertexType::Index);
    graph.vertex_colors.push_back(colorizer(index));

    const std::size_t index_vertex = nvertex - 1;

    assert(index_vertices.at(i + edges_.size()) == uninitialized_vertex);
    index_vertices[i + edges_.size()] = index_vertex;
  }

  // Add edges between proto index bundle vertices and all vertices of the
  // indices contained in that bundle i.e. if the bundle is {i_1,i_2}, the
  // bundle would be connected with vertices for i_1 and i_2
  for (const auto &[bundle, vertex] : proto_bundles) {
    for (const Index &idx : bundle) {
      std::size_t idx_vertex = uninitialized_vertex;

      auto it = std::ranges::find(edges_, idx, &Edge::idx);
      if (it != edges_.end()) {
        assert(std::distance(edges_.begin(), it) >= 0);
        idx_vertex = index_vertices.at(std::distance(edges_.begin(), it));
      } else {
        auto pure_it = pure_proto_indices_.find(idx);
        assert(pure_it != pure_proto_indices_.end());

        if (pure_it != pure_proto_indices_.end()) {
          assert(std::distance(pure_proto_indices_.begin(),
                               pure_proto_indices_.end()) >= 0);
          idx_vertex = index_vertices.at(
              std::distance(pure_proto_indices_.begin(), pure_it) +
              edges_.size());
        }
      }

      assert(idx_vertex != uninitialized_vertex);
      if (idx_vertex == uninitialized_vertex) {
        std::abort();
      }

      edges.push_back(std::make_pair(idx_vertex, vertex));
    }
  }

  assert(!options.make_labels || nvertex == graph.vertex_labels.size());
  assert(!options.make_texlabels || nvertex == graph.vertex_texlabels.size());
  assert(nvertex == graph.vertex_colors.size());
  assert(nvertex == graph.vertex_types.size());

  // Create the actual BLISS graph object
  graph.bliss_graph = std::make_unique<bliss::Graph>(nvertex);

  for (const std::pair<std::size_t, std::size_t> &current_edge : edges) {
    graph.bliss_graph->add_edge(current_edge.first, current_edge.second);
  }

  for (const auto [vertex, color] :
       ranges::views::enumerate(graph.vertex_colors)) {
    graph.bliss_graph->change_color(vertex, color);
  }

  if (options.make_idx_to_vertex) {
    assert(index_vertices.size() == edges_.size() + pure_proto_indices_.size());
    graph.idx_to_vertex.reserve(index_vertices.size());

    for (std::size_t i = 0; i < edges_.size(); ++i) {
      graph.idx_to_vertex.emplace(
          std::make_pair(edges_.at(i).idx(), index_vertices.at(i)));
    }
    for (const auto &[i, index] :
         ranges::views::enumerate(pure_proto_indices_)) {
      graph.idx_to_vertex.emplace(
          std::make_pair(index, index_vertices.at(i + edges_.size())));
    }
  }

  return graph;
}

void TensorNetworkV2::init_edges() {
  have_edges_ = false;
  edges_.clear();
  ext_indices_.clear();
  pure_proto_indices_.clear();

  auto idx_insert = [this](const Index &idx, Vertex vertex) {
    if (idx.nonnull() == false) return;

    if (Logger::instance().tensor_network) {
      std::wcout << "TensorNetworkV2::init_edges: idx=" << to_latex(idx)
                 << " attached to tensor " << vertex.getTerminalIndex() << " ("
                 << vertex.getOrigin() << ") at position "
                 << vertex.getIndexSlot()
                 << " (sym: " << to_wstring(vertex.getTerminalSymmetry()) << ")"
                 << std::endl;
    }

    auto it = std::ranges::lower_bound(edges_, idx, Index::FullLabelCompare{},
                                       &Edge::idx);
    if (it == edges_.end() || it->idx() != idx) {
      edges_.emplace(it, std::move(vertex), &idx);
    } else {
      it->connect_to(std::move(vertex));
    }
  };

  std::size_t distinct_index_estimate = 0;
  for (const AbstractTensorPtr &current : tensors_) {
    distinct_index_estimate += bra_rank(*current);
    distinct_index_estimate += ket_rank(*current);
    distinct_index_estimate += aux_rank(*current);
  }
  // For a fully contracted tensor network 1/2 of all indices are unique
  // so that can be regarded as a kind of lower bound
  distinct_index_estimate /= 2;
  edges_.reserve(distinct_index_estimate);

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

      const auto &[it, inserted] = ext_indices_.emplace(current.idx());
      // only scenario where idx is already in ext_indices_ if it were a
      // protoindex of a previously inserted ext index ... check to ensure no
      // accidental duplicates
      if (!inserted) {
        assert(proto_indices.contains(current.idx()));
      }
    }

    // add proto indices to the grand list of proto indices
    for (auto &&proto_idx : current.idx().proto_indices()) {
      // for now no recursive proto indices
      if (proto_idx.has_proto_indices())
        throw std::runtime_error(
            "TensorNetworkV2 does not support recursive protoindices");
      proto_indices.emplace(proto_idx);
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

  // some external indices will have protoindices that are NOT among
  // pure_proto_indices_, e.g.
  // i2 in f_i2^{a2^{i1,i2}} t_{a2^{i1,i2}a3^{i1,i2}}^{i2,i1}
  // is not added to ext_indices_ due to being among doubly-connected edges_
  // and thus is not among pure_proto_indices_, but it needs to be
  // and external index due to a3^{i1,i2} being an external index
  NamedIndexSet ext_proto_indices;
  ranges::for_each(ext_indices_, [&](const auto &idx) {
    ranges::for_each(idx.proto_indices(), [&](const auto &pidx) {
      if (!pure_proto_indices_.contains(
              pidx))  // only add indices that not already in
                      // pure_proto_indices_, which will be added to
                      // ext_indices_ below
        ext_proto_indices.emplace(pidx);
    });
  });
  ext_indices_.reserve(ext_indices_.size() + ext_proto_indices.size());
  ext_indices_.insert(ext_proto_indices.begin(), ext_proto_indices.end());

  // ... and add pure protoindices to the external indices
  ext_indices_.reserve(ext_indices_.size() + pure_proto_indices_.size());
  ext_indices_.insert(pure_proto_indices_.begin(), pure_proto_indices_.end());

  have_edges_ = true;
}

container::svector<std::pair<long, long>> TensorNetworkV2::factorize() {
  abort();  // not yet implemented
}

size_t TensorNetworkV2::SlotCanonicalizationMetadata::hash_value() const {
  return graph->get_hash();
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
    [[maybe_unused]] const TensorCanonicalizer &tensor_canonizer =
        nondefault_canonizer_ptr ? *nondefault_canonizer_ptr : canonicalizer;

    auto bp = canonicalizer.apply(*tensor);

    if (bp) {
      byproduct *= bp;
    }
  }

  return byproduct;
}

}  // namespace sequant
