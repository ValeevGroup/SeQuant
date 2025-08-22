//
// Created by Eduard Valeyev on 2025-24-07.
//

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
#include <SeQuant/core/tensor_network_v3.hpp>
#include <SeQuant/core/utility/swap.hpp>
#include <SeQuant/core/utility/tuple.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

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
  if (terminal_symm == Symmetry::nonsymm) {
    return index_slot < rhs.index_slot;
  } else {
    return false;
  }
}

bool TensorNetworkV3::Vertex::operator==(const Vertex &rhs) const {
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

std::size_t TensorNetworkV3::Graph::vertex_to_index_idx(
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

  assert(tensor_idx > 0);
  return tensor_idx - 1;
}

ExprPtr TensorNetworkV3::canonicalize_graph(
    const NamedIndexSet &named_indices) {
  int parity = 1;

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV3::canonicalize_graph: input tensors\n";
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

  if (Logger::instance().canonicalize_input_graph) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout, {.labels = graph.vertex_labels});
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
    cgraph->write_dot(std::wcout, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    std::wcout << "with our labels:\n";
    cgraph->write_dot(std::wcout, {.labels = cvlabels});
    delete cgraph;
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
  // for nonsymmetric particle-symmetric tensors only: maps tensor ordinal ->
  // canonical order of its braket slots
  container::map<std::size_t, container::svector<std::size_t, 4>>
      canonical_braket_slot_order;
  // for bra-ket symmetric tensors only: maps tensor ordinal -> canonical order
  // of its bra and ket slot bundle vertices
  container::map<std::size_t, std::array<std::size_t, /* bra + ket = */ 2>>
      canonical_bra_ket_bundle_order;

  std::vector<std::size_t> index_idx_to_vertex;
  index_idx_to_vertex.reserve(edges_.size() + pure_proto_indices_.size());
  std::size_t tensor_braket_vertex_ord =
      0;  // counts encountered braket bundle vertices, resets to zero when
          // switching to new tensor

  for (std::size_t vertex = 0; vertex < graph.vertex_types.size(); ++vertex) {
    const auto vertex_type = graph.vertex_types[vertex];
    switch (vertex_type) {
      case VertexType::Index:
        index_idx_to_vertex.emplace_back(index_idx_to_vertex.size()) = vertex;
        break;

      case VertexType::TensorBra:
      case VertexType::TensorKet: {
        assert(tensor_count > 0);
        const auto bra = vertex_type == VertexType::TensorBra;
        const std::size_t tensor_ord = tensor_count - 1;
        const AbstractTensor &tensor = *tensors_[tensor_ord];
        const auto symm = symmetry(tensor);
        if (symm == Symmetry::symm || symm == Symmetry::antisymm) {
          canonical_slot_order[tensor_ord][bra ? 0 : 1].second.emplace_back(
              canonize_perm[vertex]);
        }
        const auto bksymm = braket_symmetry(tensor);
        if (bksymm != BraKetSymmetry::nonsymm) {
          canonical_bra_ket_bundle_order[tensor_ord][bra ? 0 : 1] =
              canonize_perm[vertex];
        }
        break;
      }

      case VertexType::TensorBraKet: {
        assert(tensor_count > 0);
        const std::size_t tensor_ord = tensor_count - 1;
        const AbstractTensor &tensor = *tensors_[tensor_ord];
        const auto symm = symmetry(tensor);
        const auto psymm = particle_symmetry(tensor);
        if (symm == Symmetry::nonsymm && psymm == ParticleSymmetry::symm &&
            /* skip the first one which connects bra and ket bundles */
            tensor_braket_vertex_ord != 0) {
          canonical_braket_slot_order[tensor_ord].emplace_back(
              canonize_perm[vertex]);
        }
        ++tensor_braket_vertex_ord;
        break;
      }

      case VertexType::TensorCore:
        assert(tensor_idx_to_vertex.size() == tensor_count);
        tensor_idx_to_vertex.emplace_back(tensor_idx_to_vertex.size()) = vertex;
        ++tensor_count;
        tensor_braket_vertex_ord = 0;
        break;

      case VertexType::TensorAux:
      case VertexType::TensorBraBundle:
      case VertexType::TensorKetBundle:
      case VertexType::TensorAuxBundle:
      case VertexType::IndexBundle:
        break;
    }
  }

  assert(index_idx_to_vertex.size() ==
         edges_.size() + pure_proto_indices_.size());
  assert(tensor_idx_to_vertex.size() == tensors_.size());
  assert(canonical_slot_order.size() <= tensors_.size());

  // canonical slot arrays right now contain vertex ordinals, convert to
  // permutations
  for (auto &[ord, braparslots_ketparslots] : canonical_slot_order) {
    auto &[braparslots, ketparslots] = braparslots_ketparslots;
    braparslots.first = sort_then_replace_by_ordinals(braparslots.second);
    ketparslots.first = sort_then_replace_by_ordinals(ketparslots.second);
  }
  for (auto &[ord, slots] : canonical_braket_slot_order) {
    sort_then_replace_by_ordinals(slots);
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
    assert(lhs_idx < index_idx_to_vertex.size());
    const std::size_t lhs_vertex = index_idx_to_vertex[lhs_idx];
    assert(rhs_idx < index_idx_to_vertex.size());
    const std::size_t rhs_vertex = index_idx_to_vertex[rhs_idx];

    return canonize_perm[lhs_vertex] < canonize_perm[rhs_vertex];
  };

  sort_via_ordinals<OrderType::StrictWeak>(edges_, index_less_than);

  for (const Edge &current : edges_) {
    const Index &idx = current.idx();

    const auto is_named = current.vertex_count() == 1;
    if (is_named) continue;

    idxrepl_emplace(idx, idxfac.make(idx));
  }

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetworkV3::canonicalize_graph: replacing "
                 << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

  // The tensor reordering and index relabeling will make edges_ invalid
  edges_.clear();
  have_edges_ = false;

  apply_index_replacements(tensors_, idxrepl, true);

  // Permute {bra, ket} or braket slots of particle-symmetric tensors as
  // indicated by graph canonization
  for (std::size_t i = 0; i < tensors_.size(); ++i) {
    AbstractTensor &tensor = *tensors_[i];

    if (particle_symmetry(tensor) != ParticleSymmetry::symm) continue;
    const auto asymm = symmetry(tensor) == Symmetry::nonsymm;

    if (asymm) {  // asymmetric tensor? order braket slots only

      auto it = canonical_braket_slot_order.find(i);
      if (it == canonical_braket_slot_order.end()) continue;

      auto &sorted_ordinals = it->second;

      // the logic of _permute_braket is too complicated to capture here
      // if (Logger::instance().canonicalize) {
      //   if (!ranges::is_sorted(sorted_ordinals)) {
      //     for (const auto &idxpair : idxrepl) {
      //       std::wcout << "TensorNetworkV3::canonicalize_graph: permuting "
      //                     "braket slots in "
      //                  << to_latex(tensor) << ":\n";
      //       for (auto i = 0; i != sorted_ordinals.size(); ++i) {
      //         std::wcout << "  {" <<
      //         to_latex(tensor._bra()[sorted_ordinals[i]])
      //                    << "," <<
      //                    to_latex(tensor._ket()[sorted_ordinals[i]])
      //                    << "} -> {" << to_latex(tensor._bra()[i]) << ","
      //                    << to_latex(tensor._ket()[i]) << "}\n";
      //       }
      //       std::wcout << std::endl;
      //     }
      //   }
      // }

      tensor._permute_braket(
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
            for (const auto &idxpair : idxrepl) {
              std::wcout << "TensorNetworkV3::canonicalize_graph: permuting "
                         << (bra ? "bra" : "ket") << " slots in "
                         << to_latex(tensor) << ":\n";
              auto indices = bra ? tensor._bra() : tensor._ket();
              for (auto i = 0; i != indices.size(); ++i) {
                std::wcout << "  " << to_latex(indices[sorted_ordinals[i]])
                           << " -> " << to_latex(indices[i]) << "\n";
              }
              std::wcout << std::endl;
            }
          }
        }
      }

      tensor._permute_bra(std::span(braslots.data(), braslots.size()));
      tensor._permute_ket(std::span(ketslots.data(), ketslots.size()));

      // parity of slot permutations only matters for antisymmetric tensors
      if (symmetry(tensor) == Symmetry::antisymm) {
        parity *= braparity.value_or(1) * ketparity.value_or(1);
      }
    }

    // lastly permute bra with ket bundles, if needed
    // TODO extend to support comjugate case
    if (braket_symmetry(tensor) != BraKetSymmetry::symm) continue;

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
    std::wcout << "TensorNetworkV3::canonicalize_graph: tensors after "
                  "canonicalization\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
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
    tensors_.emplace_back(std::shared_ptr<AbstractTensor>(t->_clone()));
  }
  tensor_input_ordinals_ = other.tensor_input_ordinals_;
}

TensorNetworkV3 &TensorNetworkV3::operator=(
    const TensorNetworkV3 &other) noexcept {
  *this = TensorNetworkV3(other);
  return *this;
}

ExprPtr TensorNetworkV3::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels, bool fast,
    const NamedIndexSet *named_indices_ptr) {
  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV3::canonicalize(" << (fast ? "fast" : "slow")
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

  ExprPtr byproduct;
  if (!fast) {
    // The graph-based canonization is required in all cases in which there are
    // indistinguishable tensors present in the expression. Their order and
    // indexing can only be determined via this rigorous canonization.
    byproduct = canonicalize_graph(named_indices);
  }

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
    std::wcout << "TensorNetworkV3::canonicalize(" << (fast ? "fast" : "slow")
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
      std::wcout << "TensorNetworkV3::canonicalize(" << (fast ? "fast" : "slow")
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

TensorNetworkV3::SlotCanonicalizationMetadata
TensorNetworkV3::canonicalize_slots(
    const container::vector<std::wstring> &cardinal_tensor_labels,
    const NamedIndexSet *named_indices_ptr,
    TensorNetworkV3::SlotCanonicalizationMetadata::named_index_compare_t
        named_index_compare) {
  if (!named_index_compare)
    named_index_compare = [](const auto &idxptr_slottype_1,
                             const auto &idxptr_slottype_2) -> bool {
      const auto &[idxptr1, slottype1] = idxptr_slottype_1;
      const auto &[idxptr2, slottype2] = idxptr_slottype_2;
      return idxptr1->space() < idxptr2->space();
    };

  TensorNetworkV3::SlotCanonicalizationMetadata metadata;

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetworkV3::canonicalize_slots(): input tensors\n";
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
  // graph.bliss_graph->write_dot(std::wcout, {.labels = graph.vertex_labels});

  if (Logger::instance().canonicalize_input_graph) {
    std::wcout << "Input graph for canonicalization:\n";
    graph.bliss_graph->write_dot(std::wcout,
                                 {.labels = graph.vertex_labels,
                                  .xlabels = graph.vertex_xlabels,
                                  .texlabels = graph.vertex_texlabels});
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
    metadata.graph->write_dot(std::wcout, {.display_colors = true});
    auto cvlabels = permute(graph.vertex_labels, canonize_perm);
    auto cvtexlabels = permute(graph.vertex_texlabels, canonize_perm);
    std::wcout << "with our labels:\n";
    metadata.graph->write_dot(std::wcout,
                              {.labels = cvlabels, .texlabels = cvtexlabels});
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
          if (edge_it->vertex(0).getOrigin() == Origin::Aux)
            slot_type = IndexSlotType::TensorAux;
          else if (edge_it->vertex(0).getOrigin() == Origin::Bra)
            slot_type = IndexSlotType::TensorBra;
          else {
            assert(edge_it->vertex(0).getOrigin() == Origin::Ket);
            slot_type = IndexSlotType::TensorKet;
          }
        } else {  // if
          assert(edge_it->vertex_count() == 2);
          slot_type = IndexSlotType::IndexBundle;
        }
      } else
        slot_type = IndexSlotType::IndexBundle;
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

TensorNetworkV3::Graph TensorNetworkV3::create_graph(
    const CreateGraphOptions &options) const {
  assert(have_edges_);

  // initialize named_indices by default to all external indices
  const NamedIndexSet &named_indices = options.named_indices == nullptr
                                           ? this->ext_indices()
                                           : *(options.named_indices);

  VertexPainter colorizer(named_indices, options.distinct_named_indices);

  constexpr std::size_t num_tensor_components = 5;

  // results
  Graph graph;
  std::size_t nvertex = 0;
  // predicting exact vertex count is too much work, make a rough estimate only
  // We know that at the very least all indices and all tensors will yield
  // vertex representations; for tensors estimate the average number of verices
  // at 5

  std::size_t vertex_count_estimate =
      edges_.size() + pure_proto_indices_.size() + 5 * tensors_.size();
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

  // computes ordinal of the vertex representing index slot of type
  // slot_type which is slot_ordinal'th (empty or nonempty) slot in the slot
  // bundle to obtain ordinal of the slot vertex add this to to tensor_vertex
  // (i.e. ordinal of the tensor core vertex)
  auto index_slot_offset = [](const AbstractTensor &tensor, SlotType slot_type,
                              std::size_t slot_ordinal) {
    const Symmetry tensor_sym = symmetry(tensor);
    const bool is_symm = tensor_sym != Symmetry::nonsymm;
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
    assert(tensor_idx < tensor_vertices.size());
    assert(tensor_vertices[tensor_idx] == uninitialized_vertex);
    assert(tensors_.at(tensor_idx));
    const AbstractTensor &tensor = *tensors_[tensor_idx];

    // Tensor core
    const auto tlabel = label(tensor);
    if (options.make_labels) graph.vertex_labels.emplace_back(tlabel);
    if (options.make_texlabels)
      graph.vertex_texlabels.emplace_back(L"$" + utf_to_latex(tlabel) + L"$");
    graph.vertex_types.emplace_back(VertexType::TensorCore);
    const auto tensor_color =
        colorizer.apply_shade(tensor);  // subsequent vertices will be shaded by
                                        // the color of this tensor
    graph.vertex_colors.emplace_back(tensor_color);

    const std::size_t tensor_vertex = nvertex;
    tensor_vertices[tensor_idx] = tensor_vertex;
    ++nvertex;

    const Symmetry tensor_sym = symmetry(tensor);
    const bool is_symm = tensor_sym != Symmetry::nonsymm;
    // max (number of bra slots, number of ket slots) slots, i.e. the number of
    // 1-index and 2-index columns
    const std::size_t num_particles =
        std::max(bra_rank(tensor), ket_rank(tensor));
    // min (number of bra slots, number of ket slots) slots, i.e. the number of
    // 2-index columns
    const std::size_t num_paired_particles =
        std::max(bra_rank(tensor), ket_rank(tensor));
    const bool is_braket_symm = braket_symmetry(tensor) == BraKetSymmetry::symm;

    // vertices for braket bundles:
    // - antisymmetric/symmetric tensors only need 1 bundle for {bra,ket}
    // - asymmetric tensors also need 1 bundle for each pair of slots
    // {bra_i,ket_i} (including pairs where one of the slots is empty/missing)
    const std::size_t num_braket_vertices = !is_symm ? num_particles + 1 : 1;
    const bool is_part_symm =
        particle_symmetry(tensor) == ParticleSymmetry::symm;

    // make braket slot bundles first
    for (std::size_t i = 0; i < num_braket_vertices; ++i) {
      if (options.make_labels || options.make_texlabels) {
        std::wstring base_label = L"bk";
        std::wstring psuffix;
        if (i == 0) {  // {bra,ket} bundle -> "bk{a,s,}"
          switch (tensor_sym) {
            case Symmetry::symm:
              base_label += L"s";
              break;
            case Symmetry::antisymm:
              base_label += L"a";
              break;
            default:
              assert(tensor_sym != Symmetry::invalid);
          }
        } else {
          psuffix = L"_" + to_wstring(i);
        }
        if (options.make_labels)
          graph.vertex_labels.emplace_back(base_label + psuffix);
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(
              base_label + ((i != 0) ? (L"\\" + psuffix) : L""));
      }
      graph.vertex_types.emplace_back(VertexType::TensorBraKet);

      // If tensor is particle-symmetric use same color for all braket vertices,
      // else use different colors
      std::size_t color_id;
      if (i == 0) {  // {bra,ket} bundle -> 0
        color_id = 0;
      } else {
        // {bra_i,ket_i} bundle -> particle_symmetric ? 1 : i+1
        color_id = is_part_symm ? 1 : i;
      }
      graph.vertex_colors.emplace_back(colorizer(ParticleGroup{color_id}));

      edges.emplace_back(std::make_pair(tensor_vertex, nvertex));
      ++nvertex;
    }

    // create vertices for bra and ket slot bundles of any symmetry
    // N.B. TNV1/TNV2 created such vertices for symmetric/anstisymmetric
    // bra/ket also but did not create index slots. Here we create them
    // even for asymmetric bra/ket
    {
      for (auto s : {Origin::Bra, Origin::Ket}) {
        const bool bra = s == Origin::Bra;
        const auto size = bra ? bra_rank(tensor) : ket_rank(tensor);
        if (options.make_labels) {
          std::wstring label =
              std::wstring(bra ? L"bra" : L"ket") + std::to_wstring(size) +
              ((tensor_sym == Symmetry::antisymm)
                   ? L"a"
                   : (tensor_sym == Symmetry::symm ? L"s" : L""));
          graph.vertex_labels.emplace_back(label);
        }
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(std::nullopt);
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

    // - Create vertex for every index slot, regardless of symmetry (V1 and V2
    // only created slots for antisymmetric/symmetric tensors)
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
        // different colors to braket slots of different types so must track
        // types of braket slots:
        // - if tensor is not particle symmetric braket slots will already be
        // colored uniquely (by particle index)
        // - if tensor is symmetric/antisymmetric braket slots have same color
        // - if tensor is particle symmetric then assign different colors to
        // braket slots of different types (paired vs unpaired)

        // N.B. emtpy slots are not skipped!

        const auto is_paired_particle = i < num_paired_particles;
        std::size_t color_id = i;
        if (is_symm)
          color_id = 0;
        else if (is_part_symm) {
          if (is_paired_particle)
            color_id = 0;
          else
            color_id = 1;
        }

        if (options.make_labels)
          graph.vertex_labels.emplace_back((is_bra ? L"bra_" : L"ket_") +
                                           std::to_wstring(i + 1));
        if (options.make_texlabels)
          graph.vertex_texlabels.emplace_back(
              std::wstring(is_bra ? L"bra" : L"ket") + L"\\_" +
              std::to_wstring(i + 1));
        graph.vertex_types.emplace_back(vertex_type);
        // see color_id definition for handling of bra, ket, and and braket
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
              (!is_symm ? num_particles : 0);
          const std::size_t slot_vertex =
              tensor_vertex + slot_bundle_vertex_offset +
              /* bra or ket bundle vertex */ (is_bra ? 0 : 1);
          edges.emplace_back(std::make_pair(slot_vertex, nvertex));
        }
        // for asymmetric tensors also connect to the {bra_i,ket_i} bundle
        // vertex
        if (!is_symm) {
          const std::size_t braket_bundle_vertex =
              tensor_vertex +
              /* tensor core vertex */ 1 +
              /* {bra,ket} bundle vertex */ 1 +
              /* {bra_i,ket_i} bundle vertex */ i;
          edges.emplace_back(std::make_pair(braket_bundle_vertex, nvertex));
        }
        // make sure logic in index_slot_offset is correct
        assert(nvertex ==
               tensor_vertex + index_slot_offset(tensor, slot_type, i));
        ++nvertex;
      }
    }  // bra+ket slots

    // TODO: handle aux indices permutation symmetries once they are supported
    // for now, auxiliary indices are considered to always be asymmetric
    for (std::size_t i = 0; i < aux_rank(tensor); ++i) {
      if (options.make_labels)
        graph.vertex_labels.emplace_back(L"aux_" + std::to_wstring(i + 1));
      if (options.make_texlabels)
        graph.vertex_texlabels.emplace_back(std::wstring(L"aux") + L"\\_" +
                                            std::to_wstring(i + 1));
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

  for (std::size_t i = 0; i < edges_.size(); ++i) {
    const Edge &current_edge = edges_[i];

    const Index &index = current_edge.idx();
    if (options.make_labels)
      graph.vertex_labels.emplace_back(std::wstring(index.full_label()));
    using namespace std::string_literals;
    if (options.make_texlabels)
      graph.vertex_texlabels.emplace_back(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.emplace_back(VertexType::Index);
    graph.vertex_colors.emplace_back(colorizer(index));

    const std::size_t index_vertex = nvertex;
    ++nvertex;

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
        if (options.make_labels) {
          using namespace std::literals;
          std::wstring index_bundle_label =
              L"<" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.full_label(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L">";
          graph.vertex_labels.emplace_back(std::move(index_bundle_label));
        }
        if (options.make_texlabels) {
          using namespace std::literals;
          std::wstring index_bundle_texlabel =
              L"$\\langle" +
              (ranges::views::transform(
                   index.proto_indices(),
                   [](const Index &idx) { return idx.to_latex(); }) |
               ranges::views::join(L","sv) | ranges::to<std::wstring>()) +
              L"\\rangle$";
          graph.vertex_texlabels.emplace_back(std::move(index_bundle_texlabel));
        }
        graph.vertex_types.emplace_back(VertexType::IndexBundle);
        graph.vertex_colors.emplace_back(colorizer(index.proto_indices()));

        proto_vertex = nvertex;
        proto_bundles.emplace_back(index.proto_indices(), proto_vertex);
        ++nvertex;
      }

      edges.emplace_back(std::make_pair(index_vertex, proto_vertex));
    }

    // Connect index to the tensor(s) it is connected to
    for (std::size_t i = 0; i < current_edge.vertex_count(); ++i) {
      const Vertex &vertex = current_edge.vertex(i);
      if (i >= 2)  // hyperedges can only occur between aux indices
        assert(vertex.getOrigin() == Origin::Aux);

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
      const std::size_t offset =
          index_slot_offset(tensor, vertex.getOrigin(), vertex.getIndexSlot());
      const std::size_t tensor_component_vertex = tensor_vertex + offset;

      assert(tensor_component_vertex < nvertex);
      edges.emplace_back(std::make_pair(index_vertex, tensor_component_vertex));
    }
  }

  // also create vertices for pure proto indices
  for (const auto &[i, index] : ranges::views::enumerate(pure_proto_indices_)) {
    ++nvertex;
    if (options.make_labels)
      graph.vertex_labels.emplace_back(index.full_label());
    using namespace std::string_literals;
    if (options.make_texlabels)
      graph.vertex_texlabels.emplace_back(L"$"s + index.to_latex() + L"$");
    graph.vertex_types.emplace_back(VertexType::Index);
    graph.vertex_colors.emplace_back(colorizer(index));

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
        idx_vertex = index_vertices[std::distance(edges_.begin(), it)];
      } else {
        auto pure_it = pure_proto_indices_.find(idx);
        assert(pure_it != pure_proto_indices_.end());

        if (pure_it != pure_proto_indices_.end()) {
          assert(std::distance(pure_proto_indices_.begin(),
                               pure_proto_indices_.end()) >= 0);
          idx_vertex = index_vertices[std::distance(pure_proto_indices_.begin(),
                                                    pure_it) +
                                      edges_.size()];
        }
      }

      assert(idx_vertex != uninitialized_vertex);
      if (idx_vertex == uninitialized_vertex) {
        std::abort();
      }

      edges.emplace_back(std::make_pair(idx_vertex, vertex));
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

void TensorNetworkV3::init_edges() {
  have_edges_ = false;
  edges_.clear();
  ext_indices_.clear();
  pure_proto_indices_.clear();

  auto idx_insert = [this](const Index &idx, Vertex vertex) {
    // skip null indices
    if (!idx) return;
    if (Logger::instance().tensor_network) {
      std::wcout << "TensorNetworkV3::init_edges: idx=" << to_latex(idx)
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
    distinct_index_estimate += bra_rank(*current);  // assumes no empty slots
    distinct_index_estimate += ket_rank(*current);  // assumes no empty slots
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

container::svector<std::pair<long, long>> TensorNetworkV3::factorize() {
  abort();  // not yet implemented
}

size_t TensorNetworkV3::SlotCanonicalizationMetadata::hash_value() const {
  return graph->get_hash();
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
