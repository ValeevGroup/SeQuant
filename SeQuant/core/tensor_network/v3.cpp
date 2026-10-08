//
// Created by Eduard Valeyev on 2025-24-07.
//

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/tag.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/utils.hpp>
#include <SeQuant/core/tensor_network/v3.hpp>
#include <SeQuant/core/tensor_network/vertex_painter.hpp>
#include <SeQuant/core/utility/debug.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/utility/swap.hpp>
#include <SeQuant/core/utility/tuple.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip.hpp>

namespace sequant {

namespace {

using Permutation = std::vector<unsigned int>;

/// @return generators of the automorphisms of @p graph that fix each of
/// @p fixed_vertices
std::vector<Permutation> stabilizer_generators(
    const TensorNetworkV3::Graph &graph,
    const std::vector<unsigned int> &fixed_vertices) {
  const auto nvertices = graph.bliss_graph->get_nof_vertices();
  Permutation identity(nvertices);
  std::iota(identity.begin(), identity.end(), 0);
  std::unique_ptr<bliss::Graph> individualized(
      graph.bliss_graph->permute(identity));
  // each fixed vertex gets a color of its own
  container::set<TensorNetworkV3::Graph::VertexColor> used(
      graph.vertex_colors.begin(), graph.vertex_colors.end());
  auto color = std::numeric_limits<TensorNetworkV3::Graph::VertexColor>::max();
  for (const auto vertex : fixed_vertices) {
    while (used.contains(color)) --color;
    individualized->change_color(vertex, color);
    used.insert(color);
  }
  std::vector<Permutation> generators;
  using hook_t = std::function<void(unsigned int, const unsigned int *)>;
  hook_t hook = [&generators](unsigned int n, const unsigned int *aut) {
    generators.emplace_back(aut, aut + n);
  };
  bliss::Stats stats;
  individualized->find_automorphisms(stats, &bliss::aut_hook<hook_t>, &hook);
  return generators;
}

/// Among the canonical labelings of @p graph, which are @p labeling composed
/// with its automorphisms, selects the one that places the named indices in
/// label order: position by position, in canonical order, each position of a
/// named index gets the one with the smallest label that an automorphism
/// fixing the positions before it can bring there. This makes the labeling
/// a function of the expression when named indices are not told apart by label
/// in @p graph.
/// @param labeling a canonical labeling of @p graph (input vertex to
///        canonical position)
/// @param generators generators of the automorphism group of @p graph
/// @param named_indices the named indices of the network
/// @param label_less the order of index labels, that of the tensor
///        canonicalizers so that a lone tensor is ordered alike; must be a
///        strict total order on @p named_indices, else the placement would
///        depend on the input vertex numbering
/// @return the selected labeling
Permutation order_named_indices_by_label(
    const TensorNetworkV3::Graph &graph, Permutation labeling,
    std::vector<Permutation> generators,
    const TensorNetworkV3::NamedIndexSet &named_indices,
    const tensor_index_comparer_t &label_less) {
  const auto nvertices = labeling.size();
  const auto is_named = [&](unsigned int vertex) {
    return graph.vertex_types[vertex] == VertexType::Index &&
           named_indices.contains(graph.vertex_indices[vertex]);
  };
  std::vector<unsigned int> positions;
  for (unsigned int vertex = 0; vertex != nvertices; ++vertex)
    if (is_named(vertex)) positions.push_back(labeling[vertex]);
  if (positions.size() < 2) return labeling;
  std::sort(positions.begin(), positions.end());

  Permutation vertex_at(nvertices);
  const auto invert = [&] {
    for (unsigned int vertex = 0; vertex != nvertices; ++vertex)
      vertex_at[labeling[vertex]] = vertex;
  };
  invert();
  std::vector<unsigned int> placed;
  for (const auto position : positions) {
    const auto current = vertex_at[position];
    // the orbit of current under the automorphisms that fix the placed
    // vertices, each with an automorphism that maps current to it
    container::map<unsigned int, Permutation> transversal;
    {
      Permutation identity(nvertices);
      std::iota(identity.begin(), identity.end(), 0);
      transversal.emplace(current, std::move(identity));
    }
    std::vector<unsigned int> frontier{current};
    while (!frontier.empty()) {
      const auto vertex = frontier.back();
      frontier.pop_back();
      for (const auto &generator : generators) {
        const auto image = generator[vertex];
        if (transversal.contains(image)) continue;
        const auto &to_vertex = transversal.at(vertex);
        Permutation to_image(nvertices);
        for (unsigned int v = 0; v != nvertices; ++v)
          to_image[v] = generator[to_vertex[v]];
        transversal.emplace(image, std::move(to_image));
        frontier.push_back(image);
      }
    }
    // an automorphism maps named indices to named ones
    SEQUANT_ASSERT(ranges::all_of(
        transversal, [&](const auto &entry) { return is_named(entry.first); }));
    const auto best =
        std::min_element(transversal.begin(), transversal.end(),
                         [&](const auto &a, const auto &b) {
                           return label_less(graph.vertex_indices[a.first],
                                             graph.vertex_indices[b.first]);
                         });
    const auto chosen = best->first;
    SEQUANT_ASSERT(ranges::none_of(transversal, [&](const auto &entry) {
      return entry.first != chosen &&
             !label_less(graph.vertex_indices[chosen],
                         graph.vertex_indices[entry.first]);
    }));
    if (chosen != current) {
      // the inverse of best->second maps chosen to current, so composing the
      // labeling with it puts chosen at this position and leaves the placed
      // vertices where they are
      const auto &to_chosen = best->second;
      Permutation inverse(nvertices);
      for (unsigned int v = 0; v != nvertices; ++v) inverse[to_chosen[v]] = v;
      Permutation composed(nvertices);
      for (unsigned int v = 0; v != nvertices; ++v)
        composed[v] = labeling[inverse[v]];
      labeling = std::move(composed);
      invert();
    }
    placed.push_back(chosen);
    // a trivial orbit means every generator fixes chosen already
    if (transversal.size() > 1)
      generators = stabilizer_generators(graph, placed);
  }
  return labeling;
}

}  // namespace

TensorNetworkV3::Vertex::Vertex(Origin origin, std::size_t terminal_idx,
                                std::size_t index_slot, Symmetry terminal_symm)
    : origin(origin),
      terminal_idx(terminal_idx),
      index_slot(index_slot),
      terminal_symm(terminal_symm) {}

SlotType TensorNetworkV3::Vertex::getOrigin() const { return origin; }

std::size_t TensorNetworkV3::Vertex::getTerminalIndex() const {
  return terminal_idx;
}

std::size_t TensorNetworkV3::Vertex::getIndexSlot() const { return index_slot; }

Symmetry TensorNetworkV3::Vertex::getTerminalSymmetry() const {
  return terminal_symm;
}

bool TensorNetworkV3::Vertex::operator<(const Vertex &rhs) const {
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
  if (terminal_symm == Symmetry::Nonsymm) {
    return index_slot < rhs.index_slot;
  } else {
    return false;
  }
}

bool TensorNetworkV3::Vertex::operator==(const Vertex &rhs) const {
  // Slot position is only taken into account for non_symmetric tensors
  const std::size_t lhs_slot =
      (terminal_symm == Symmetry::Nonsymm) * index_slot;
  const std::size_t rhs_slot =
      (rhs.terminal_symm == Symmetry::Nonsymm) * rhs.index_slot;

  // sanity check that bra and ket have same symmetry
  assert(origin == Origin::Aux || rhs.origin == Origin::Aux ||
         terminal_idx != rhs.terminal_idx ||
         terminal_symm == rhs.terminal_symm);

  return terminal_idx == rhs.terminal_idx && lhs_slot == rhs_slot &&
         origin == rhs.origin;
}

std::size_t TensorNetworkV3::Graph::vertex_to_index_idx(
    std::size_t vertex) const {
  SEQUANT_ASSERT(vertex_types.at(vertex) == VertexType::Index);

  std::size_t index_idx = 0;
  for (std::size_t i = 0; i <= vertex; ++i) {
    if (vertex_types[i] == VertexType::Index) {
      ++index_idx;
    }
  }

  SEQUANT_ASSERT(index_idx > 0);

  return index_idx - 1;
}

std::optional<std::size_t> TensorNetworkV3::Graph::vertex_to_tensor_idx(
    std::size_t vertex) const {
  const auto vertex_type = vertex_types[vertex];
  if (vertex_type == VertexType::Index || vertex_type == VertexType::SPBundle)
    return std::nullopt;

  std::size_t tensor_idx = 0;
  for (std::size_t i = 0; i <= vertex; ++i) {
    if (vertex_types[i] == VertexType::TensorCore) {
      ++tensor_idx;
    }
  }

  SEQUANT_ASSERT(tensor_idx > 0);
  return tensor_idx - 1;
}

ExprPtr TensorNetworkV3::canonicalize_graph(const NamedIndexSet &named_indices,
                                            bool ignore_named_index_labels) {
  int parity = 1;

  if (Logger::instance().canonicalize) {
    std::wostringstream oss;
    oss << "TensorNetworkV3::canonicalize_graph: input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      oss << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    oss << std::endl;
    sequant::wprintf(oss.str());
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
       .distinct_named_indices = !ignore_named_index_labels,
       .make_labels = Logger::instance().canonicalize_input_graph ||
                      Logger::instance().canonicalize_dot,
       .make_texlabels = Logger::instance().canonicalize_input_graph ||
                         Logger::instance().canonicalize_dot});

  if (Logger::instance().canonicalize_input_graph) {
    std::wostringstream oss;
    oss << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(oss, {.labels = graph.vertex_labels});
    sequant::wprintf(oss.str());
  }

  // a network with an automorphism of phase -1 equals minus itself, i.e. is
  // zero; phase is a homomorphism from the automorphism group to {+1,-1}, so
  // if a scored generator has phase -1 the term is zero. Generators that are
  // not scored (see Graph::automorphism_phase) can only hide a zero, never
  // invent one.
  bool has_odd_automorphism = false;
  // with labels ignored, the graph does not tell named indices apart, so the
  // canonical labeling places them only up to its automorphisms; they are then
  // placed by label, which needs the automorphism group
  std::vector<Permutation> automorphisms;
  const bool collect_automorphisms = ignore_named_index_labels;
  const unsigned int *bliss_labeling = canonicalize_graph(
      graph, graph.antisymm_bundles.empty() && !collect_automorphisms
                 ? std::function<void(unsigned int, const unsigned int *)>{}
                 : [&](unsigned int n, const unsigned int *aut) {
                     if (collect_automorphisms)
                       automorphisms.emplace_back(aut, aut + n);
                     if (!graph.antisymm_bundles.empty() &&
                         !has_odd_automorphism &&
                         graph.automorphism_phase(aut, &named_indices) == -1)
                       has_odd_automorphism = true;
                   });

  if (has_odd_automorphism) {
    if (Logger::instance().canonicalize)
      sequant::wprintf(
          "TensorNetworkV3::canonicalize_graph: automorphism of phase -1 "
          "found, the network is zero\n");
    return ex<Constant>(0);
  }

  Permutation labeling(bliss_labeling,
                       bliss_labeling + graph.bliss_graph->get_nof_vertices());
  if (collect_automorphisms && !automorphisms.empty())
    labeling = order_named_indices_by_label(
        graph, std::move(labeling), std::move(automorphisms), named_indices,
        get_default_context_snapshot().index_comparer());
  const unsigned int *canonize_perm = labeling.data();

  if (Logger::instance().canonicalize_dot) {
    std::wostringstream oss;
    oss << "Canonicalization permutation:\n";
    for (std::size_t i = 0; i < graph.vertex_labels.size(); ++i) {
      oss << i << " -> " << canonize_perm[i] << "\n";
    }
    oss << "Canonicalized graph:\n";
    bliss::Graph *cgraph = graph.bliss_graph->permute(canonize_perm);
    cgraph->write_dot(oss, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    oss << "with our labels:\n";
    cgraph->write_dot(oss, {.labels = cvlabels});
    delete cgraph;
    sequant::wprintf(oss.str());
  }

  // maps tensor ordinal -> input vertex ordinal
  std::vector<std::size_t> tensor_idx_to_vertex;
  tensor_idx_to_vertex.reserve(tensors_.size());
  std::size_t tensor_count = 0;

  // for symmetric tensors only: maps tensor ordinal -> canonical order of
  // its bra and ket slots
  container::map<
      std::size_t,
      std::array<std::pair</* permutation parity */ std::optional<int>,
                           container::svector<std::size_t, 4>>,
                 /* bra + ket = */ 2>>
      canonical_slot_order;
  // for nonsymmetric column-symmetric tensors only: maps tensor ordinal ->
  // canonical order of its column bundle vertices (one per column).
  // Populated from VertexType::TensorBraKet (column bundle) vertices below,
  // not from individual bra/ket slot vertices.
  container::map<std::size_t, container::svector<std::size_t, 4>>
      canonical_column_bundle_order;
  // for bra-ket symmetric tensors only: maps tensor ordinal -> canonical order
  // of its bra and ket slot bundle vertices
  container::map<std::size_t, std::array<std::size_t, /* bra + ket = */ 2>>
      canonical_bra_ket_bundle_order;

  std::vector<std::size_t> index_idx_to_vertex;
  index_idx_to_vertex.reserve(edges_.size() + pure_proto_indices_.size());
  std::size_t tensor_braket_vertex_ord =
      0;  // counts encountered braket and column bundle vertices,
          // resets to zero when switching to new tensor

  for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
    const auto vertex_type = graph.vertex_types[vertex];
    switch (vertex_type) {
      case VertexType::Index:
        index_idx_to_vertex.emplace_back(index_idx_to_vertex.size()) = vertex;
        break;

      case VertexType::TensorBra:
      case VertexType::TensorKet: {
        SEQUANT_ASSERT(tensor_count > 0);
        const auto bra = vertex_type == VertexType::TensorBra;
        const std::size_t tensor_ord = tensor_count - 1;
        const AbstractTensor &tensor = *tensors_[tensor_ord];
        const auto symm = symmetry(tensor);
        if (symm == Symmetry::Symm || symm == Symmetry::Antisymm) {
          canonical_slot_order[tensor_ord][bra ? 0 : 1].second.emplace_back(
              canonize_perm[vertex]);
        }
        break;
      }

      case VertexType::TensorBraBundle:
      case VertexType::TensorKetBundle: {
        // The braket swap decision must compare canon positions of the
        // *bundle* vertices, not of individual slot vertices. Using slot
        // vertices (as we used to) overwrote on every iteration → the
        // comparison ended up reading canon_perm of the last bra/ket slot,
        // which is not necessarily invariant under the canonical graph's
        // automorphism group (same-color same-degree vertices can be
        // permuted among themselves in valid canonical labelings).
        SEQUANT_ASSERT(tensor_count > 0);
        const auto bra = vertex_type == VertexType::TensorBraBundle;
        const std::size_t tensor_ord = tensor_count - 1;
        const AbstractTensor &tensor = *tensors_[tensor_ord];
        const auto bksymm = braket_symmetry(tensor);
        if (bksymm != BraKetSymmetry::Nonsymm) {
          canonical_bra_ket_bundle_order[tensor_ord][bra ? 0 : 1] =
              canonize_perm[vertex];
        }
        break;
      }

      case VertexType::TensorBraKet: {
        SEQUANT_ASSERT(tensor_count > 0);
        const std::size_t tensor_ord = tensor_count - 1;
        const AbstractTensor &tensor = *tensors_[tensor_ord];
        const auto symm = symmetry(tensor);
        const auto csymm = column_symmetry(tensor);
        if (symm == Symmetry::Nonsymm && csymm == ColumnSymmetry::Symm &&
            /* skip the first one which connects bra and ket bundles */
            tensor_braket_vertex_ord != 0) {
          canonical_column_bundle_order[tensor_ord].emplace_back(
              canonize_perm[vertex]);
        }
        ++tensor_braket_vertex_ord;
        break;
      }

      case VertexType::TensorCore:
        SEQUANT_ASSERT(tensor_idx_to_vertex.size() == tensor_count);
        tensor_idx_to_vertex.emplace_back(tensor_idx_to_vertex.size()) = vertex;
        ++tensor_count;
        tensor_braket_vertex_ord = 0;
        break;

      case VertexType::TensorAux:
      case VertexType::TensorAuxBundle:
      case VertexType::IndexBundle:
        break;
    }
  }

  SEQUANT_ASSERT(index_idx_to_vertex.size() ==
                 edges_.size() + pure_proto_indices_.size());
  SEQUANT_ASSERT(tensor_idx_to_vertex.size() == tensors_.size());
  SEQUANT_ASSERT(canonical_slot_order.size() <= tensors_.size());

  // canonical slot arrays right now contain vertex ordinals, convert to
  // permutations
  for (auto &[ord, braparslots_ketparslots] : canonical_slot_order) {
    auto &[braparslots, ketparslots] = braparslots_ketparslots;
    braparslots.first = sort_then_replace_by_ordinals(braparslots.second);
    ketparslots.first = sort_then_replace_by_ordinals(ketparslots.second);
  }
  for (auto &[ord, bundles] : canonical_column_bundle_order) {
    sort_then_replace_by_ordinals(bundles);
  }

  container::map<Index, Index> idxrepl;
  auto idxrepl_emplace = [&idxrepl](auto &&from, auto &&to) {
    if (from != to) idxrepl.emplace(std::move(from), std::move(to));
  };

  // Sort edges so that their order corresponds to the order of indices in the
  // canonical graph
  // Use this ordering to relabel anonymous indices
  const auto index_less_than = [&index_idx_to_vertex, &canonize_perm](
                                   std::size_t lhs_idx, std::size_t rhs_idx) {
    SEQUANT_ASSERT(lhs_idx < index_idx_to_vertex.size());
    const std::size_t lhs_vertex = index_idx_to_vertex[lhs_idx];
    SEQUANT_ASSERT(rhs_idx < index_idx_to_vertex.size());
    const std::size_t rhs_vertex = index_idx_to_vertex[rhs_idx];

    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  sort_via_ordinals<OrderType::StrictWeak>(edges_, index_less_than);

  for (const Edge &current : edges_) {
    const Index &idx = current.idx();

    const auto is_named = !is_anonymous_index(idx);
    if (is_named) continue;

    idxrepl_emplace(idx, idxfac.make(idx));
  }

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      sequant::wprintf("TensorNetworkV3::canonicalize_graph: replacing ",
                       to_latex(idxpair.first), " with ",
                       to_latex(idxpair.second), "\n");
    }
  }

  // The tensor reordering and index relabeling will make edges_ invalid
  edges_.clear();
  have_edges_ = false;

  apply_index_replacements(tensors_, idxrepl, true);

  // Permute {bra, ket} or column slots of column-symmetric tensors as
  // indicated by graph canonization
  for (std::size_t i = 0; i < tensors_.size(); ++i) {
    AbstractTensor &tensor = *tensors_[i];

    if (column_symmetry(tensor) != ColumnSymmetry::Symm) continue;
    const auto asymm = symmetry(tensor) == Symmetry::Nonsymm;

    if (asymm) {  // asymmetric tensor? order column slots only

      auto it = canonical_column_bundle_order.find(i);
      if (it == canonical_column_bundle_order.end()) continue;

      auto &sorted_ordinals = it->second;

      tensor._permute_columns(
          std::span(sorted_ordinals.data(), sorted_ordinals.size()));
    } else {  // symmetric/antisymmetric bra
      auto it = canonical_slot_order.find(i);
      if (it == canonical_slot_order.end()) continue;

      auto &[braparslots, ketparslots] = it->second;
      auto &[braparity, braslots] = braparslots;
      auto &[ketparity, ketslots] = ketparslots;

      if (Logger::instance().canonicalize) {
        for (auto bk : {Origin::Bra, Origin::Ket}) {
          const auto bra = bk == Origin::Bra;
          auto &sorted_ordinals = bra ? braslots : ketslots;
          if (!ranges::is_sorted(sorted_ordinals)) {
            sequant::wprintf("TensorNetworkV3::canonicalize_graph: permuting ",
                             (bra ? "bra" : "ket"), " slots in ",
                             to_latex(tensor), ":\n");
            auto indices = bra ? tensor._bra() : tensor._ket();
            for (auto i = 0; i != indices.size(); ++i) {
              sequant::wprintf("  ", to_latex(indices[sorted_ordinals[i]]),
                               " -> ", to_latex(indices[i]), "\n");
            }
            sequant::wprintf("\n");
          }
        }
      }

      tensor._permute_bra(std::span(braslots.data(), braslots.size()));
      tensor._permute_ket(std::span(ketslots.data(), ketslots.size()));

      // parity of slot permutations only matters for antisymmetric tensors
      if (symmetry(tensor) == Symmetry::Antisymm) {
        parity *= braparity.value_or(1) * ketparity.value_or(1);
      }
    }

    // lastly permute bra with ket bundles, if needed
    // TODO extend to support conjugate case
    if (braket_symmetry(tensor) != BraKetSymmetry::Symm) continue;

    // swap bra and ket bundles
    if (canonical_bra_ket_bundle_order[i][0] >
        canonical_bra_ket_bundle_order[i][1]) {
      tensor._swap_bra_ket();
    }
  }

  // Less-than relationship for tensors. Tensors that do not commute are
  // equivalent,i .e.g tensors `a` and `b` are equivalent if
  // `!(a<b) && !(b<a)`).
  // Possibility of non-commutativity breaks transitivity (e.g. given tensor of
  // operators `a` and `b` and a tensor of scalars `c` both `a<c` and `c<b`
  // can be, but this does not imply `a<b`.
  const auto tensor_less_than = [this, &canonize_perm, &tensor_idx_to_vertex](
                                    std::size_t lhs_idx, std::size_t rhs_idx) {
    const AbstractTensor &lhs = *tensors_[lhs_idx];
    const AbstractTensor &rhs = *tensors_[rhs_idx];

    if (!tensors_commute(lhs, rhs)) {
      return false;
    }

    const std::size_t lhs_vertex = tensor_idx_to_vertex[lhs_idx];
    const std::size_t rhs_vertex = tensor_idx_to_vertex[rhs_idx];

    // Commuting tensors are sorted based on their canonical order which is
    // given by the order of the corresponding vertices in the canonical graph
    // representation
    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  tensor_input_ordinals_ =
      sort_via_ordinals<OrderType::Weak>(tensors_, tensor_less_than);

  if (Logger::instance().canonicalize) {
    std::wostringstream oss;
    oss << "TensorNetworkV3::canonicalize_graph: tensors after "
           "canonicalization\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      oss << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    sequant::wprintf(oss.str());
  }

  if (parity < 0)
    return ex<Constant>(-1);
  else
    return {};
}

TensorNetworkV3::TensorNetworkV3(TensorNetworkV3 &&) noexcept = default;
TensorNetworkV3 &TensorNetworkV3::operator=(TensorNetworkV3 &&) noexcept =
    default;

TensorNetworkV3::TensorNetworkV3(const TensorNetworkV3 &other) {
  tensors_.reserve(other.tensors_.size());
  for (const auto &t : other.tensors_) {
    tensors_.emplace_back(t->_clone_shared());
  }
  tensor_input_ordinals_ = other.tensor_input_ordinals_;
}

TensorNetworkV3 &TensorNetworkV3::operator=(
    const TensorNetworkV3 &other) noexcept {
  *this = TensorNetworkV3(other);
  return *this;
}

ExprPtr TensorNetworkV3::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels,
    const CanonicalizeOptions &options) {
  if (Logger::instance().canonicalize) {
    std::wostringstream oss;
    oss << "TensorNetworkV3::canonicalize(" << to_wstring(options.method)
        << "): input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      oss << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    oss << "cardinal_tensor_labels = ";
    ranges::for_each(cardinal_tensor_labels,
                     [&oss](auto &&i) { oss << i << L" "; });
    oss << std::endl;
    sequant::wprintf(oss.str());
  }

  if (!have_edges_) {
    init_edges();
  }

  // initialize named_indices by default to all external indices
  using NamedIndexSet = tensor_network::NamedIndexSet;
  std::shared_ptr<NamedIndexSet> named_indices_sptr;
  if (options.named_indices) {
    named_indices_sptr = std::make_shared<NamedIndexSet>(
        options.named_indices->begin(), options.named_indices->end());
  }
  const auto &named_indices =
      !options.named_indices ? this->ext_indices() : *named_indices_sptr;

  if (Logger::instance().canonicalize) {
    std::wostringstream oss;
    oss << "named_indices = ";
    ranges::for_each(named_indices,
                     [&oss](auto &&i) { oss << i.full_label() << L" "; });
    sequant::wprintf(oss.str());
  }

  ExprPtr byproduct;
  if ((options.method & CanonicalizationMethod::Topological) ==
      CanonicalizationMethod::Topological) {
    // The graph-based canonization is required in all cases in which there are
    // indistinguishable tensors present in the expression. Their order and
    // indexing can only be determined via this rigorous canonization.
    byproduct = canonicalize_graph(
        named_indices, static_cast<bool>(options.ignore_named_index_labels));
    if (byproduct && byproduct->as<Constant>().is_zero()) return byproduct;
  }

  if ((options.method & CanonicalizationMethod::Lexicographic) ==
      CanonicalizationMethod::Lexicographic) {
    // Ensure each individual tensor is written in the way that its tensor
    // block (== order of index spaces) is canonical
    byproduct *= canonicalize_individual_tensor_blocks(named_indices);

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
      std::wostringstream oss;
      oss << "TensorNetworkV3::canonicalize(" << to_wstring(options.method)
          << "): tensors after lexicographic sort\n";
      size_t cnt = 0;
      ranges::for_each(tensors_, [&](const auto &t) {
        oss << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
      });
      sequant::wprintf(oss.str());
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
    // accordingly (but only anonymous indices, of course); a named index need
    // not be an index of this network, so the named edges are skipped by
    // membership, not by count
    for (const auto &edge : edges_) {
      const Index &index = edge.idx();
      if (is_named_index(index)) continue;
      Index replacement = idxfac.make(index);
      if (index != replacement) idxrepl.emplace(index, std::move(replacement));
    }

    // Done computing canonical index replacement list
    // reset edges since renamings will make them obsolete
    edges_.clear();
    have_edges_ = false;

    if (Logger::instance().canonicalize) {
      for (const auto &idxpair : idxrepl) {
        sequant::wprintf(
            "TensorNetworkV3::canonicalize(", to_wstring(options.method),
            "): lexicographic rewrite, replacing ", to_latex(idxpair.first),
            " with ", to_latex(idxpair.second), "\n");
      }
    }

    apply_index_replacements(tensors_, idxrepl, true);

    byproduct *= canonicalize_individual_tensors(named_indices);

    // We assume that re-indexing did not change the canonical order of tensors
    SEQUANT_ASSERT(
        std::is_sorted(tensors_.begin(), tensors_.end(), tensor_sorter));
    // However, in order to produce the most aesthetically pleasing result, we
    // now reorder tensors based on the regular AbstractTensor::operator<, which
    // takes the explicit index labelling of tensors into account.
    tensor_sorter.set_blocks_only(false);
    std::stable_sort(tensors_.begin(), tensors_.end(), tensor_sorter);
  }  // lexicographic canonicalization

  if (byproduct) {
    SEQUANT_ASSERT(byproduct->is<Constant>());
    return (byproduct->as<Constant>().value() == 1) ? nullptr : byproduct;
  } else
    return nullptr;
}

TensorNetworkV3::SlotCanonicalizationMetadata
TensorNetworkV3::canonicalize_slots(
    const container::vector<std::wstring> &cardinal_tensor_labels,
    const NamedIndexSet *named_indices_ptr,
    TensorNetworkV3::SlotCanonicalizationMetadata::named_index_compare_t
        named_index_compare,
    const tensor_network::NamedIndexColorMap *named_index_colors) {
  if (!named_index_compare)
    named_index_compare = [](const auto &idxptr_slottype_1,
                             const auto &idxptr_slottype_2) -> bool {
      const auto &[idxptr1, slottype1] = idxptr_slottype_1;
      const auto &[idxptr2, slottype2] = idxptr_slottype_2;
      return idxptr1->space() < idxptr2->space();
    };

  TensorNetworkV3::SlotCanonicalizationMetadata metadata;

  if (Logger::instance().canonicalize) {
    std::wostringstream oss;
    oss << "TensorNetworkV3::canonicalize_slots(): input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      oss << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    oss << "cardinal_tensor_labels = ";
    ranges::for_each(cardinal_tensor_labels,
                     [&oss](auto &&i) { oss << i << L" "; });
    oss << std::endl;
    sequant::wprintf(oss.str());
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
       .named_index_colors = named_index_colors,
       .distinct_named_indices = false,
       .make_labels = Logger::instance().canonicalize_input_graph ||
                      Logger::instance().canonicalize_dot,
       .make_texlabels = Logger::instance().canonicalize_input_graph ||
                         Logger::instance().canonicalize_dot,
       .make_idx_to_vertex = true});
  const auto &idx_to_vertex = graph.idx_to_vertex;

  if (Logger::instance().canonicalize_input_graph) {
    std::wostringstream oss;
    oss << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(oss, {.labels = graph.vertex_labels,
                                       .xlabels = graph.vertex_xlabels,
                                       .texlabels = graph.vertex_texlabels});
    sequant::wprintf(oss.str());
  }

  const unsigned int *canonize_perm = canonicalize_graph(graph);

  metadata.graph =
      std::shared_ptr<bliss::Graph>(graph.bliss_graph->permute(canonize_perm));

  if (Logger::instance().canonicalize_dot) {
    std::wostringstream oss;
    oss << "Canonicalization permutation:\n";
    for (std::size_t i = 0; i < graph.vertex_labels.size(); ++i) {
      oss << i << " -> " << canonize_perm[i] << "\n";
    }
    oss << "Canonicalized graph:\n";
    metadata.graph->write_dot(oss, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    auto cvtexlabels = permute(graph.vertex_texlabels, canonize_perm);
    oss << "with our labels:\n";
    metadata.graph->write_dot(oss,
                              {.labels = cvlabels, .texlabels = cvtexlabels});
    sequant::wprintf(oss.str());
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
    for (auto [idx_ord, idx] : ranges::views::enumerate(grand_index_list)) {
      if (!is_named_index(idx)) {
        continue;
      }

      const auto named_indices_it = metadata.named_indices.find(idx);
      SEQUANT_ASSERT(named_indices_it != metadata.named_indices.end());
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
          const BraKetSymmetry symm = braket_symmetry(
              *tensors_.at(edge_it->vertex(0).getTerminalIndex()));

          if (edge_it->vertex(0).getOrigin() == Origin::Aux) {
            slot_type = IndexSlotType::TensorAux;
          } else if (symm == BraKetSymmetry::Symm ||
                     edge_it->vertex(0).getOrigin() == Origin::Bra) {
            // Note: we must not distinguis bra and ket indices in case braket
            // symmetry is present Technically, this should (to some degree)
            // also apply to BraKetSymmetry::Conjugate but this TN
            // implementation currently doesn't exploit conjugate braket
            // symmetry (as it is not entirely clear how to handle the required
            // complex conjugation)
            slot_type = IndexSlotType::TensorBra;
          } else {
            SEQUANT_ASSERT(edge_it->vertex(0).getOrigin() == Origin::Ket);
            slot_type = IndexSlotType::TensorKet;
          }
        } else if (edge_it->vertex_count() == 2) {
          // an internal contraction edge that is also named (e.g. a proto
          // index on a named index) connects exactly 2 slots
          slot_type = IndexSlotType::IndexBundle;
        } else {
          // a high-order hyperindex: a named index shared among more than 2
          // tensor slots (e.g. an auxiliary/batching index common to many
          // factors, as in Laplace-transform MP2 or tensor hypercontraction).
          // Supported only for auxiliary indices, which carry no vector-space
          // (bra/ket) character; classify by the aux slot it occupies. A
          // non-auxiliary high-order hyperindex has no well-defined slot type,
          // so hard-error (in every build config, not just assertion-enabled
          // ones) rather than silently mis-canonicalize it.
          bool all_aux = true;
          for (std::size_t v = 0; v != edge_it->vertex_count(); ++v)
            if (edge_it->vertex(v).getOrigin() != Origin::Aux) {
              all_aux = false;
              break;
            }
          if (!all_aux)
            throw Exception(
                "TensorNetworkV3::canonicalize_slots: high-order (shared among "
                ">2 tensor slots) non-auxiliary hyperindices are not "
                "supported");
          slot_type = IndexSlotType::TensorAux;
        }
      } else
        slot_type = IndexSlotType::IndexBundle;
      const auto idxptr_slottype = std::make_pair(&idx, slot_type);
      auto it = idx2cord.find(idxptr_slottype);

      if (it == idx2cord.end()) {
        bool inserted;
        std::tie(it, inserted) = idx2cord.emplace(
            idxptr_slottype, cord_set_t(cord_set_t::key_compare{}));
        SEQUANT_ASSERT(inserted);
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
    if (symmetry(tensor) != Symmetry::Antisymm) {
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

      if (bubble_sort_parity(vertices) == -1) {
        // Performed an uneven amount of pairwise exchanges -> this incurs a
        // phase change
        metadata.phase *= -1;
      }
    }
  }

  return metadata;
}

TensorNetworkV3::Graph TensorNetworkV3::create_graph(
    const CreateGraphOptions &options) const {
  if (!have_edges_) const_cast<TensorNetworkV3 *>(this)->init_edges();

  // initialize named_indices by default to all external indices
  const NamedIndexSet &named_indices = options.named_indices == nullptr
                                           ? this->ext_indices()
                                           : *(options.named_indices);

  VertexPainter colorizer(named_indices, options.distinct_named_indices,
                          options.named_index_colors);

  // results
  Graph graph;
  std::size_t nvertex = 0;

  auto make_label = [&nvertex, &options, &graph](std::wstring lbl) {
    graph.vertex_labels.emplace_back(
        options.label_prepend_ordinal
            ? (std::to_wstring(nvertex) + L": " + std::move(lbl))
            : std::move(lbl));
  };
  auto make_xlabel = [&nvertex, &options, &graph]() {
    if (options.xlabel_maker)
      graph.vertex_xlabels.emplace_back(options.xlabel_maker(nvertex));
  };
  auto make_texlabel = [&nvertex, &options, &graph](std::wstring lbl) {
    if (lbl.empty())
      graph.vertex_texlabels.emplace_back(std::nullopt);
    else
      graph.vertex_texlabels.emplace_back(
          options.texlabel_prepend_ordinal
              ? (std::to_wstring(nvertex) + L": " + std::move(lbl))
              : std::move(lbl));
  };

  // predicting exact vertex count is too much work, make a rough estimate only
  // We know that at the very least all indices and all tensors will yield
  // vertex representations; for tensors estimate the average number of verices
  // at 5
  std::size_t vertex_count_estimate =
      edges_.size() + pure_proto_indices_.size() + 5 * tensors_.size();

  if (options.make_labels) graph.vertex_labels.reserve(vertex_count_estimate);
  if (options.make_xlabels) graph.vertex_xlabels.reserve(vertex_count_estimate);
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

  // computes ordinal of the vertex representing index slot of type
  // slot_type which is slot_ordinal'th (empty or nonempty) slot in the slot
  // bundle to obtain ordinal of the slot vertex add this to to tensor_vertex
  // (i.e. ordinal of the tensor core vertex)
  auto index_slot_offset = [](const AbstractTensor &tensor, SlotType slot_type,
                              std::size_t slot_ordinal) {
    const Symmetry tensor_sym = symmetry(tensor);
    const bool is_symm = tensor_sym != Symmetry::Nonsymm;
    std::size_t offset = 0;
    // number of vertices before first index slot vertex varies with symmetry
    if (is_symm) {
      offset += /* {bra,ket} bundle vertex */ 1 +
                /* bra and ket bundle vertices */ 2;
    } else {
      auto nbraket = std::max(bra_rank(tensor), ket_rank(tensor));
      offset += /* {bra,ket} bundle vertex */ 1 +
                /* {bra_i,ket_i} bundle vertices */ nbraket +
                /* bra and ket bundle vertices */ 2;
    }

    // now count slot vertices
    // N.B. empty slots are NOT skipped to avoid having to map nonempty slot
    // ordinal to overall slot ordinal
    if (slot_type == SlotType::Bra)
      offset += slot_ordinal;
    else if (slot_type == SlotType::Ket)
      offset += bra_rank(tensor) + slot_ordinal;
    else
      offset += bra_rank(tensor) + ket_rank(tensor) + slot_ordinal;

    return offset + 1;  // +1 to account for tensor core vertex
  };

  // Add vertices for tensors
  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    SEQUANT_ASSERT(tensor_idx < tensor_vertices.size());
    SEQUANT_ASSERT(tensor_vertices[tensor_idx] == uninitialized_vertex);
    SEQUANT_ASSERT(tensors_.at(tensor_idx));
    const AbstractTensor &tensor = *tensors_[tensor_idx];

    // Tensor core
    const auto tlabel = label(tensor);
    if (options.make_labels) make_label(std::wstring{tlabel});
    if (options.make_xlabels) make_xlabel();
    if (options.make_texlabels)
      make_texlabel(L"$" + io::latex::utf_to_string(tlabel) + L"$");
    graph.vertex_types.emplace_back(VertexType::TensorCore);
    const auto tensor_color =
        colorizer.apply_shade(tensor);  // subsequent vertices will be shaded by
                                        // the color of this tensor
    graph.vertex_colors.emplace_back(tensor_color);

    const std::size_t tensor_vertex = nvertex;
    tensor_vertices[tensor_idx] = tensor_vertex;
    ++nvertex;

    const Symmetry tensor_sym = symmetry(tensor);
    const bool is_symm = tensor_sym != Symmetry::Nonsymm;
    // max (number of bra slots, number of ket slots) slots, i.e. the number of
    // 1-index and 2-index columns
    const std::size_t num_cols = std::max(bra_rank(tensor), ket_rank(tensor));
    // min (number of bra slots, number of ket slots) slots, i.e. the number of
    // 2-index columns
    const std::size_t num_paired_cols =
        std::min(bra_rank(tensor), ket_rank(tensor));
    const bool is_braket_symm = braket_symmetry(tensor) == BraKetSymmetry::Symm;

    // vertices for braket bundles:
    // - antisymmetric/symmetric tensors only need 1 bundle for {bra,ket}
    // - asymmetric tensors also need 1 bundle for each column, i.e. pair of
    // matching slots {bra_i,ket_i} (including pairs where one of the slots is
    // empty/missing)
    const std::size_t num_braket_vertices = !is_symm ? num_cols + 1 : 1;
    const bool is_col_symm = column_symmetry(tensor) == ColumnSymmetry::Symm;

    // make braket slot bundles first
    for (std::size_t i = 0; i < num_braket_vertices; ++i) {
      if (options.make_labels || options.make_texlabels) {
        std::wstring base_label = L"c";
        std::wstring psuffix;
        if (i == 0) {  // {bra,ket} bundle -> "bk{a,s,}"
          switch (tensor_sym) {
            case Symmetry::Symm:
              base_label += L"s";
              break;
            case Symmetry::Antisymm:
              base_label += L"a";
              break;
            case Symmetry::Nonsymm:
              break;
          }
        } else {
          psuffix = L"_" + std::to_wstring(i);
        }
        if (options.make_labels) make_label(base_label + psuffix);
        if (options.make_texlabels)
          make_texlabel(base_label + ((i != 0) ? (L"\\" + psuffix) : L""));
      }
      if (options.make_xlabels) make_xlabel();
      graph.vertex_types.emplace_back(VertexType::TensorBraKet);

      // If tensor is column-symmetric use same color for all braket vertices,
      // else use different colors
      std::size_t color_id;
      if (i == 0) {  // {bra,ket} bundle -> 0
        color_id = 0;
      } else {
        // {bra_i,ket_i} bundle -> column_symmetric ? 1 : i+1
        color_id = is_col_symm ? 1 : i;
      }
      graph.vertex_colors.emplace_back(colorizer(ColumnGroup{color_id}));

      edges.emplace_back(std::make_pair(tensor_vertex, nvertex));
      ++nvertex;
    }

    // create vertices for bra and ket slot bundles of any symmetry
    {
      for (auto s : {Origin::Bra, Origin::Ket}) {
        const bool bra = s == Origin::Bra;
        const auto size = bra ? bra_rank(tensor) : ket_rank(tensor);
        if (options.make_labels || options.make_texlabels) {
          std::wstring label =
              std::wstring(bra ? L"bra" : L"ket") + std::to_wstring(size) +
              ((tensor_sym == Symmetry::Antisymm)
                   ? L"a"
                   : (tensor_sym == Symmetry::Symm ? L"s" : L""));
          if (options.make_labels) make_label(label);
          if (options.make_texlabels) make_texlabel(label);
        }
        if (options.make_xlabels) make_xlabel();
        graph.vertex_types.emplace_back(bra ? VertexType::TensorBraBundle
                                            : VertexType::TensorKetBundle);
        tensor_network::VertexColor color;
        if (is_braket_symm) {  // if have bra<->ket symmetry (not conj!),
                               // use same color for bra and ket
          color = colorizer(BraGroup{size});
        } else {
          color = bra ? colorizer(BraGroup{size}) : colorizer(KetGroup{size});
        }
        graph.vertex_colors.emplace_back(color);

        const auto braket_vertex = tensor_vertex + 1;
        edges.emplace_back(std::make_pair(braket_vertex, nvertex));
        ++nvertex;
      }
    }

    // Graph::automorphism_phase scores the slots of antisymmetric bundles
    const std::size_t antisymm_bra_bundle = tensor_sym == Symmetry::Antisymm
                                                ? graph.antisymm_bundles.size()
                                                : uninitialized_vertex;
    if (tensor_sym == Symmetry::Antisymm)
      graph.antisymm_bundles.resize(graph.antisymm_bundles.size() + 2);

    // - Create vertex for every index slot, regardless of symmetry
    for (auto &slot_type : {SlotType::Bra, SlotType::Ket}) {
      const auto is_bra = slot_type == SlotType::Bra;
      const auto vertex_type =
          is_bra ? VertexType::TensorBra : VertexType::TensorKet;
      const auto nslots = is_bra ? bra_rank(tensor) : ket_rank(tensor);
      auto slots = is_bra ? tensor._bra() : tensor._ket();
      for (std::size_t i = 0; i < nslots; ++i) {
        // N.B. currently AbstractTensor only supports "left"-aligned bra/ket
        // slot sets (i.e. bra[0] is paired with ket[0], etc.), gaps between
        // occupied slots are occupied by null indices) we need to assign
        // different colors to column slots of different types so must track
        // types of column slots:
        // - if tensor is not column symmetric column slots will already be
        // colored uniquely (by column index)
        // - if tensor is symmetric/antisymmetric column slots have same color
        // - if tensor is column symmetric then assign different colors to
        // column slots of different types (paired vs unpaired)

        // N.B. emtpy slots are not skipped!

        const auto is_paired_col = i < num_paired_cols;
        std::size_t color_id = i;
        if (is_symm)
          color_id = 0;
        else if (is_col_symm) {
          if (is_paired_col)
            color_id = 0;
          else
            color_id = 1;
        }

        if (options.make_labels)
          make_label((is_bra ? L"bra_" : L"ket_") + std::to_wstring(i + 1));
        if (options.make_xlabels) make_xlabel();
        if (options.make_texlabels)
          make_texlabel(std::wstring(is_bra ? L"bra" : L"ket") + L"\\_" +
                        std::to_wstring(i + 1));
        graph.vertex_types.emplace_back(vertex_type);
        // see color_id definition for handling of bra, ket, and column
        // bundle symmetries. if symmetric wrt bra<->ket swap use same color
        // for bra and ket bundles, else use distinct colors
        graph.vertex_colors.emplace_back((is_bra || is_braket_symm)
                                             ? colorizer(BraGroup{color_id})
                                             : colorizer(KetGroup{color_id}));

        // connect to bra bundle vertex, regardless of symmetry
        {
          const std::size_t slot_bundle_vertex_offset =
              /* tensor core vertex */ 1 + /* {bra,ket} bundle vertex */ 1 +
              /* {bra_i,ket_i} bundle vertices */
              (!is_symm ? num_cols : 0);
          const std::size_t slot_vertex =
              tensor_vertex + slot_bundle_vertex_offset +
              /* bra or ket bundle vertex */ (is_bra ? 0 : 1);
          edges.emplace_back(std::make_pair(slot_vertex, nvertex));
        }
        // for asymmetric tensors also connect to the {bra_i,ket_i} column
        // bundle vertex
        if (!is_symm) {
          const std::size_t column_bundle_vertex =
              tensor_vertex +
              /* tensor core vertex */ 1 +
              /* {bra,ket} bundle vertex */ 1 +
              /* {bra_i,ket_i} column bundle vertex */ i;
          edges.emplace_back(std::make_pair(column_bundle_vertex, nvertex));
        }
        // make sure logic in index_slot_offset is correct
        assert(nvertex ==
               tensor_vertex + index_slot_offset(tensor, slot_type, i));
        if (antisymm_bra_bundle != uninitialized_vertex && slots[i].nonnull())
          graph.antisymm_bundles[antisymm_bra_bundle + (is_bra ? 0 : 1)]
              .push_back(nvertex);
        ++nvertex;
      }
    }  // bra+ket slots

    // TODO: handle aux indices permutation symmetries once they are supported
    // for now, auxiliary indices are considered to always be asymmetric
    for (std::size_t i = 0; i < aux_rank(tensor); ++i) {
      if (options.make_labels) make_label(L"aux_" + std::to_wstring(i + 1));
      if (options.make_xlabels) make_xlabel();
      if (options.make_texlabels)
        make_texlabel(std::wstring(L"aux") + L"\\_" + std::to_wstring(i + 1));
      graph.vertex_types.emplace_back(VertexType::TensorAux);
      graph.vertex_colors.emplace_back(colorizer(AuxGroup{i}));
      edges.emplace_back(std::make_pair(tensor_vertex, nvertex));
      ++nvertex;
    }

    colorizer.reset_shade();
  }

  // Now add all indices (edges_ + pure_proto_indices_) to the graph
  container::vector<std::size_t> index_vertices;
  index_vertices.resize(edges_.size() + pure_proto_indices_.size(),
                        uninitialized_vertex);

  const bool strict_braket_symmetry =
      get_default_context().assert_strict_braket_symmetry();
  for (std::size_t i = 0; i < edges_.size(); ++i) {
    const Edge &current_edge = edges_[i];

    const Index &index = current_edge.idx();
    if (options.make_labels) make_label(std::wstring(index.full_label()));
    if (options.make_xlabels) make_xlabel();
    using namespace std::string_literals;
    if (options.make_texlabels) make_texlabel(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.emplace_back(VertexType::Index);
    graph.vertex_colors.emplace_back(colorizer(index));

    const std::size_t index_vertex = nvertex;
    ++nvertex;

    SEQUANT_ASSERT(index_vertices.at(i) == uninitialized_vertex);
    index_vertices[i] = index_vertex;

    // Handle proto indices
    if (index.has_proto_indices()) {
      // For now we assume that all proto indices are symmetric
      SEQUANT_ASSERT(index.symmetric_proto_indices());

      std::size_t proto_vertex;
      if (auto it =
              std::ranges::find(proto_bundles, index.proto_indices(),
                                &decltype(proto_bundles)::value_type::first);
          it != proto_bundles.end()) {
        proto_vertex = it->second;
      } else {
        // Create a new vertex for this bundle of proto indices
        if (options.make_labels) {
          using namespace std::literals;
          std::wstring index_bundle_label =
              L"<" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.full_label(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L">";
          make_label(std::move(index_bundle_label));
        }
        if (options.make_xlabels) make_xlabel();
        if (options.make_texlabels) {
          using namespace std::literals;
          std::wstring index_bundle_texlabel =
              L"$\\langle" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.to_latex(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L"\\rangle$";
          make_texlabel(std::move(index_bundle_texlabel));
        }
        graph.vertex_types.emplace_back(VertexType::IndexBundle);
        graph.vertex_colors.emplace_back(colorizer(index.proto_indices()));

        proto_vertex = nvertex;
        proto_bundles.emplace_back(index.proto_indices(), proto_vertex);
        ++nvertex;
      }

      edges.emplace_back(std::make_pair(index_vertex, proto_vertex));
    }

    // strict bra-ket sanity checks
    {
      if (strict_braket_symmetry) {
        // dummy (anonymous) edges to
        // - involve at most 2 bra and/or ket indices (if BraKetSymmetry::Symm)
        // or 1 bra and 1 ket index
        // - can involve any number of aux indices
        if (current_edge.vertex_count() > 1) {
          // ignore if named index
          if (!this->ext_indices_.contains(current_edge.idx())) {
            std::size_t nbra = 0;
            std::size_t nket = 0;
            [[maybe_unused]] std::size_t naux = 0;
            BraKetSymmetry symm = BraKetSymmetry::Nonsymm;
            for (std::size_t v = 0; v < current_edge.vertex_count(); ++v) {
              const Vertex &vertex = current_edge.vertex(v);
              switch (vertex.getOrigin()) {
                case Origin::Bra:
                  ++nbra;
                  break;
                case Origin::Ket:
                  ++nket;
                  break;
                case Origin::Aux:
                  ++naux;
                  break;
                case Origin::Proto:
                  SEQUANT_UNREACHABLE;
              }

              if (symm != BraKetSymmetry::Symm) {
                // We only care if at least one of the vertices has symmetric
                // braket symm
                symm = braket_symmetry(*tensors_[vertex.getTerminalIndex()]);
              }
            }

            // if braket symmetry == BraKetSymmetry::Symm there is no
            // distinction between bra and ket, but still can have at most 2 of
            // them total if braket symmetry != BraKetSymmetry::Symm at most 1
            // bra and 1 ket can connect to aux
            if (symm == BraKetSymmetry::Symm ? (nbra + nket > 2)
                                             : (nbra > 1 || nket > 1)) {
              throw Exception(
                  "TensorNetworkV3: index " +
                  toUtf8(current_edge.idx().full_label()) +
                  " is contracted between two bra slots (or two ket slots) "
                  "of tensors without bra-ket symmetry; a contraction pairs a "
                  "bra slot with a ket slot");
            }
          }
        }
      }
    }

    // Connect index to the tensor(s) it is connected to
    for (std::size_t i = 0; i < current_edge.vertex_count(); ++i) {
      const Vertex &vertex = current_edge.vertex(i);

      SEQUANT_ASSERT(vertex.getTerminalIndex() < tensor_vertices.size());
      SEQUANT_ASSERT(tensor_vertices[vertex.getTerminalIndex()] !=
                     uninitialized_vertex);
      const std::size_t tensor_vertex =
          tensor_vertices[vertex.getTerminalIndex()];

      // Store an edge connecting the index vertex to the corresponding tensor
      // vertex
      const AbstractTensor &tensor = *tensors_[vertex.getTerminalIndex()];
      const std::size_t offset =
          index_slot_offset(tensor, vertex.getOrigin(), vertex.getIndexSlot());
      const std::size_t tensor_component_vertex = tensor_vertex + offset;

      SEQUANT_ASSERT(tensor_component_vertex < nvertex);
      edges.emplace_back(std::make_pair(index_vertex, tensor_component_vertex));
    }
  }

  // also create vertices for pure proto indices
  for (const auto &[i, index] : ranges::views::enumerate(pure_proto_indices_)) {
    if (options.make_labels) make_label(std::wstring(index.full_label()));
    if (options.make_xlabels) make_xlabel();
    using namespace std::string_literals;
    if (options.make_texlabels) make_texlabel(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.emplace_back(VertexType::Index);
    graph.vertex_colors.emplace_back(colorizer(index));

    const std::size_t index_vertex = nvertex;

    SEQUANT_ASSERT(index_vertices.at(i + edges_.size()) ==
                   uninitialized_vertex);
    index_vertices[i + edges_.size()] = index_vertex;
    ++nvertex;
  }

  // Add edges between proto index bundle vertices and all vertices of the
  // indices contained in that bundle i.e. if the bundle is {i_1,i_2}, the
  // bundle would be connected with vertices for i_1 and i_2
  for (const auto &[bundle, vertex] : proto_bundles) {
    for (const Index &idx : bundle) {
      std::size_t idx_vertex = uninitialized_vertex;

      auto it = std::ranges::find(edges_, idx, &Edge::idx);
      if (it != edges_.end()) {
        SEQUANT_ASSERT(std::distance(edges_.begin(), it) >= 0);
        idx_vertex = index_vertices[std::distance(edges_.begin(), it)];
      } else {
        auto pure_it = pure_proto_indices_.find(idx);
        SEQUANT_ASSERT(pure_it != pure_proto_indices_.end());

        if (pure_it != pure_proto_indices_.end()) {
          SEQUANT_ASSERT(std::distance(pure_proto_indices_.begin(),
                                       pure_proto_indices_.end()) >= 0);
          idx_vertex = index_vertices[std::distance(pure_proto_indices_.begin(),
                                                    pure_it) +
                                      edges_.size()];
        }
      }

      SEQUANT_ASSERT(idx_vertex != uninitialized_vertex);
      if (idx_vertex == uninitialized_vertex) {
        SEQUANT_ABORT("Expected all vertices to be initialized at this point");
      }

      edges.emplace_back(std::make_pair(idx_vertex, vertex));
    }
  }

  SEQUANT_ASSERT(!options.make_labels || nvertex == graph.vertex_labels.size());
  SEQUANT_ASSERT(!options.make_texlabels ||
                 nvertex == graph.vertex_texlabels.size());
  SEQUANT_ASSERT(nvertex == graph.vertex_colors.size());
  SEQUANT_ASSERT(nvertex == graph.vertex_types.size());

  // Create the actual BLISS graph object
  graph.bliss_graph =
      std::make_unique<bliss::Graph>(static_cast<unsigned int>(nvertex));

  for (const std::pair<std::size_t, std::size_t> &current_edge : edges) {
    graph.bliss_graph->add_edge(current_edge.first, current_edge.second);
  }

  for (const auto [vertex, color] :
       ranges::views::enumerate(graph.vertex_colors)) {
    graph.bliss_graph->change_color(vertex, color);
  }

  // what Graph::automorphism_phase needs to know about the vertices
  graph.vertex_indices.resize(nvertex);
  graph.vertex_fixed_for_phase.assign(nvertex, false);
  for (std::size_t v = 0; v != nvertex; ++v) {
    switch (graph.vertex_types[v]) {
      case VertexType::TensorCore:
      case VertexType::TensorAux:
      case VertexType::TensorAuxBundle:
      case VertexType::IndexBundle:
        graph.vertex_fixed_for_phase[v] = true;
        break;
      default:
        break;
    }
  }
  for (std::size_t i = 0; i < edges_.size(); ++i) {
    const Index &idx = edges_[i].idx();
    graph.vertex_indices[index_vertices[i]] = idx;
    if (ext_indices_.contains(idx))
      graph.vertex_fixed_for_phase[index_vertices[i]] = true;
  }
  for (const auto &[i, index] : ranges::views::enumerate(pure_proto_indices_))
    graph.vertex_indices[index_vertices[i + edges_.size()]] = index;
  graph.antisymm_slots.assign(nvertex,
                              {uninitialized_vertex, uninitialized_vertex});
  for (std::size_t b = 0; b != graph.antisymm_bundles.size(); ++b)
    for (std::size_t pos = 0; pos != graph.antisymm_bundles[b].size(); ++pos)
      graph.antisymm_slots[graph.antisymm_bundles[b][pos]] = {b, pos};

  if (options.make_idx_to_vertex) {
    SEQUANT_ASSERT(index_vertices.size() ==
                   edges_.size() + pure_proto_indices_.size());
    graph.idx_to_vertex.reserve(index_vertices.size());

    for (std::size_t i = 0; i < edges_.size(); ++i) {
      graph.idx_to_vertex.emplace(
          std::make_pair(edges_[i].idx(), index_vertices[i]));
    }
    for (const auto &[i, index] :
         ranges::views::enumerate(pure_proto_indices_)) {
      graph.idx_to_vertex.emplace(
          std::make_pair(index, index_vertices[i + edges_.size()]));
    }
  }

  return graph;
}

const unsigned int *TensorNetworkV3::canonicalize_graph(
    const TensorNetworkV3::Graph &graph,
    const std::function<void(unsigned int, const unsigned int *)> &aut_hook) {
  using hook_t = std::function<void(unsigned int, const unsigned int *)>;
  bliss::Stats stats;
  graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
  return graph.bliss_graph->canonical_form(
      stats, aut_hook ? &bliss::aut_hook<const hook_t> : nullptr,
      const_cast<hook_t *>(&aut_hook));
}

int TensorNetworkV3::Graph::automorphism_phase(
    const unsigned int *aut,
    const container::set<Index, Index::FullLabelCompare> *named_indices) const {
  const std::size_t nv = vertex_types.size();
  SEQUANT_ASSERT(vertex_indices.size() == nv &&
                 vertex_fixed_for_phase.size() == nv &&
                 antisymm_slots.size() == nv);
  for (std::size_t v = 0; v != nv; ++v) {
    if (aut[v] == v) continue;
    if (vertex_fixed_for_phase[v]) return 0;
    if (named_indices && vertex_types[v] == VertexType::Index &&
        named_indices->contains(vertex_indices[v]))
      return 0;
  }

  // aut fixes every tensor core, hence maps each antisymmetric bundle onto a
  // bundle of the same tensor (itself, or its partner if the tensor is
  // bra<->ket symmetric); the phase is the product of the parities of the
  // induced maps between slot positions
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  int phase = 1;
  container::svector<std::size_t, 4> perm;
  for (std::size_t b = 0; b != antisymm_bundles.size(); ++b) {
    const auto &vertices = antisymm_bundles[b];
    perm.clear();
    [[maybe_unused]] std::size_t image_bundle = npos;
    for (const auto v : vertices) {
      const auto &[bundle, pos] = antisymm_slots[aut[v]];
      SEQUANT_ASSERT(bundle != npos && bundle / 2 == b / 2);
      SEQUANT_ASSERT(image_bundle == npos || image_bundle == bundle);
      SEQUANT_ASSERT(pos < vertices.size());
      image_bundle = bundle;
      perm.push_back(pos);
    }
    phase *= permutation_parity(std::span(perm));
  }
  return phase;
}

void TensorNetworkV3::init_edges() {
  have_edges_ = false;
  edges_.clear();
  ext_indices_.clear();
  pure_proto_indices_.clear();

  auto idx_insert = [this](const Index &idx, Vertex vertex) {
    // skip null indices
    if (!idx) return;
    if (Logger::instance().tensor_network) {
      std::wostringstream oss;
      oss << "TensorNetworkV3::init_edges: idx=" << to_latex(idx)
          << " attached to tensor " << vertex.getTerminalIndex() << " ("
          << vertex.getOrigin() << ") at position " << vertex.getIndexSlot()
          << " (sym: " << sequant::to_wstring(vertex.getTerminalSymmetry())
          << ")" << std::endl;
      sequant::wprintf(oss.str());
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
    distinct_index_estimate += bra_rank(*current);  // assumes no empty slots
    distinct_index_estimate += ket_rank(*current);  // assumes no empty slots
    distinct_index_estimate += aux_rank(*current);
  }
  // For a fully contracted tensor network 1/2 of all indices are unique
  // so that can be regarded as a kind of lower bound
  distinct_index_estimate /= 2;
  edges_.reserve(distinct_index_estimate);

  for (std::size_t tensor_idx = 0; tensor_idx < tensors_.size(); ++tensor_idx) {
    SEQUANT_ASSERT(tensors_[tensor_idx]);
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
                 Vertex(Origin::Aux, tensor_idx, index_idx, Symmetry::Nonsymm));
    }
  }

  // extract external indices and all protoindices (since some external indices
  // may be pure protoindices)
  NamedIndexSet proto_indices;
  for (const Edge &current : edges_) {
    SEQUANT_ASSERT(current.vertex_count() > 0);
    // External index (== Edge only connected to a single vertex in the
    // network)
    if (current.vertex_count() == 1) {
      if (Logger::instance().tensor_network) {
        sequant::wprintf("idx ", to_latex(current.idx()), " is external\n");
      }

      const auto &[it, inserted] = ext_indices_.emplace(current.idx());
      // only scenario where idx is already in ext_indices_ if it were a
      // protoindex of a previously inserted ext index ... check to ensure no
      // accidental duplicates
      if (!inserted) {
        SEQUANT_ASSERT(proto_indices.contains(current.idx()));
      }
    }

    // add proto indices to the grand list of proto indices
    for (auto &&proto_idx : current.idx().proto_indices()) {
      // for now no recursive proto indices
      if (proto_idx.has_proto_indices())
        throw Exception(
            "TensorNetworkV3 does not support recursive protoindices");
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
      sequant::wprintf("idx ", to_latex(idx), " is pure protoindex\n");
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

  if (Logger::instance().tensor_network) {
    sequant::wprintf("TNV3: computed external indices = ");
    ranges::for_each(ext_indices_, [](auto &index) {
      sequant::wprintf(index.full_label(), " ");
    });
    sequant::wprintf("\n");
  }

  have_edges_ = true;
}

container::svector<std::pair<long, long>> TensorNetworkV3::factorize() {
  SEQUANT_ABORT("TensorNetworkV3::factorize is not yet implemented");
}

size_t TensorNetworkV3::SlotCanonicalizationMetadata::hash_value() const {
  return graph->get_hash64();
}

ExprPtr TensorNetworkV3::canonicalize_individual_tensor_blocks(
    const NamedIndexSet &named_indices) {
  return do_individual_canonicalization(
      TensorBlockCanonicalizer(named_indices));
}

ExprPtr TensorNetworkV3::canonicalize_individual_tensors(
    const NamedIndexSet &named_indices) {
  return do_individual_canonicalization(
      DefaultTensorCanonicalizer(named_indices));
}

ExprPtr TensorNetworkV3::do_individual_canonicalization(
    const TensorCanonicalizer &canonicalizer) {
  ExprPtr byproduct = ex<Constant>(1);

  const auto ctx = get_default_context_snapshot();
  for (auto &tensor : tensors_) {
    auto nondefault_canonizer_ptr =
        ctx.nondefault_tensor_canonicalizer_ptr(tensor->_label());
    const TensorCanonicalizer &tensor_canonizer =
        nondefault_canonizer_ptr ? *nondefault_canonizer_ptr : canonicalizer;

    auto bp = tensor_canonizer.apply(*tensor);

    if (bp) {
      byproduct *= bp;
    }
  }

  return byproduct;
}

}  // namespace sequant
