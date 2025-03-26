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
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/vertex_painter.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/utility/tuple.hpp>
#include <SeQuant/core/wstring.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <type_traits>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

std::optional<std::size_t> TensorNetwork::GraphData::vertex_to_tensor_cluster(
    std::size_t vertex) const {
  const auto vertex_type = vertex_types.at(vertex);
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

ExprPtr TensorNetwork::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels, bool fast,
    const named_indices_t *named_indices_ptr) {
  ExprPtr canon_byproduct = ex<Constant>(1);
  container::svector<Edge> idx_terminals_sorted;  // to avoid memory allocs

  if (Logger::instance().canonicalize) {
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

  // - resort tensors (this cannot be done in Product::canonicalize since that
  // requires generic analysis of commutativity ... here we are dealing with
  // tensors only, so if either tensor in a pair of adjacent
  // tensors is a cnumber, they can be reordered
  using ranges::views::zip;
  using std::begin;
  using std::end;
  auto tensors_with_ordinals = zip(tensors_, tensor_input_ordinals_);
  bubble_sort(begin(tensors_with_ordinals), end(tensors_with_ordinals),
              [&cardinal_tensor_labels](const auto &first_ptr_and_ord,
                                        const auto &second_ptr_and_ord) {
                const auto &first = *(first_ptr_and_ord.first);
                const auto &second = *(second_ptr_and_ord.first);
                // grab base label if adjoint label is present
                auto base_label = [](const auto &t) {
                  if (label(t).back() == adjoint_label) {
                    return label(t).substr(0, label(t).size() - 1);
                  } else {
                    return label(t);
                  }
                };
                // tensors commute if their colors are different or either one
                // of them is a c-number
                if ((color(first) != color(second)) || is_cnumber(first) ||
                    is_cnumber(second)) {
                  const auto cardinal_tensor_labels_end =
                      end(cardinal_tensor_labels);
                  const auto first_cardinal_it =
                      std::find(begin(cardinal_tensor_labels),
                                end(cardinal_tensor_labels), base_label(first));
                  const auto second_cardinal_it = std::find(
                      begin(cardinal_tensor_labels),
                      end(cardinal_tensor_labels), base_label(second));
                  const auto first_is_cardinal =
                      first_cardinal_it != cardinal_tensor_labels_end;
                  const auto second_is_cardinal =
                      second_cardinal_it != cardinal_tensor_labels_end;
                  if (first_is_cardinal && second_is_cardinal) {
                    if (first_cardinal_it == second_cardinal_it)
                      return first < second;
                    else
                      return first_cardinal_it < second_cardinal_it;
                  } else if (first_is_cardinal)
                    return true;
                  else if (second_is_cardinal)
                    return false;
                  else  // neither is cardinal
                    return first < second;
                } else
                  return false;
              });

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize(" << (fast ? "fast" : "slow")
               << "): tensors after initial sort\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
  }

  if (edges_.empty()) init_edges();

  // initialize named_indices by default to all external indices (these HAVE
  // been computed in init_edges)
  const auto &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  // helpers to filter named ("external" in traditional use case) / anonymous
  // ("internal" in traditional use case)
  auto is_named_index = [&](const Index &idx) {
    return named_indices.find(idx) != named_indices.end();
  };
  auto is_anonymous_index = [&](const Index &idx) {
    return named_indices.find(idx) == named_indices.end();
  };
  // more efficient version of is_anonymous_index
  auto is_anonymous_index_ord = [&](const std::size_t &idx_ord) {
    assert(idx_ord < edges_.size() + pure_proto_indices_.size());
    if (idx_ord < edges_.size()) {
      const auto edge_it = edges_.begin() + idx_ord;
      assert(edge_it->size() > 0);
      return edge_it->size() == 2;
    } else  // pure proto indices are named
      return false;
  };
  auto namedness = [&](const Index &idx) {
    return is_named_index(idx) ? 1 : 0;
  };

  // fast and slow canonizations produce index replacements for anonymous
  // indices
  idxrepl_.clear();
  auto &idxrepl = idxrepl_;

  // index factory to generate anonymous indices
  IndexFactory idxfac(is_anonymous_index,
                      1);  // start reindexing anonymous indices from 1

  if (!fast) {
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

    // make the graph
    // named indices are not renamed, so each gets a distinct color
    auto [graph, vlabels, vtexlabels, vcolors, vtypes] = make_bliss_graph(
        {.named_indices = &named_indices,
         .distinct_named_indices = true,
         .make_labels = Logger::instance().canonicalize_dot,
         .make_texlabels = Logger::instance().canonicalize_dot});

    // canonize the graph
    bliss::Stats stats;
    graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
    const unsigned int *cl = graph->canonical_form(stats, nullptr, nullptr);

    if (Logger::instance().canonicalize_dot) {
      auto permute = [](const auto &vector, auto &&perm) {
        using std::size;
        auto sz = size(vector);
        std::decay_t<decltype(vector)> pvector(sz);
        for (size_t i = 0; i != sz; ++i) pvector[perm[i]] = vector[i];
        return pvector;
      };

      graph->write_dot(std::wcout, vlabels, vtexlabels);

      bliss::Graph *cgraph = graph->permute(cl);
      auto cvlabels = permute(vlabels, cl);
      auto cvtexlabels = permute(vtexlabels, cl);
      cgraph->write_dot(std::wcout, cvlabels, cvtexlabels);
      delete cgraph;
    }

    // make anonymous index replacement list
    {
      // for each color make a replacement list for bringing the indices to
      // the canonical order

      // grand list of colors
      container::set<size_t> colors;
      // maps color to the ordinals of the corresponding
      // anonymous indices + their canonical ordinals
      container::multimap<size_t, std::pair<size_t, size_t>> color2idx;
      // collect colors and anonymous indices sorted by colors
      size_t idx_ord = 0;
      for ([[maybe_unused]] auto &&ttpair : edges_) {
        if (is_anonymous_index_ord(idx_ord)) {
          auto color = vcolors[idx_ord];
          if (colors.find(color) == colors.end()) colors.insert(color);
          color2idx.emplace(color, std::make_pair(idx_ord, cl[idx_ord]));
        }
        ++idx_ord;
      }
      // for each color sort anonymous indices by canonical order
      container::svector<std::pair<size_t, size_t>>
          idx_can;  // canonically-ordered list of {index ordinal in edges_,
                    // canonical ordinal}
      for (auto &&color : colors) {
        auto [beg, end] = color2idx.equal_range(color);
        const auto sz = end - beg;
        // sz == 0 should not be possible since colors contains only anonymous
        // Index colors
        assert(sz != 0);

        // anonymous indices are regenerated using factory
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
          std::size_t ord = 0;
          for (auto &&[idx_ord, idx_ord_can] : idx_can) {
            const auto &idx = (edges_.begin() + idx_ord)->idx();
            const auto new_idx = idxfac.make(idx);
            if (idx != new_idx) idxrepl.emplace(idx, std::move(new_idx));
            ++ord;
          }
        } else if (sz == 1) {  // no need for resorting of colors with 1 index
                               // only, but still need to replace the index
          const auto it = edges_.begin() + std::get<0>(beg->second);
          const auto &idx = it->idx();
          const auto new_idx = idxfac.make(idx);
          if (idx != new_idx) idxrepl.emplace(idx, std::move(new_idx));
        }
      }
    }  // index repl

    // reorder *commuting* tensors_ to canonical order (defined by the core
    // indices)
    {
      decltype(tensors_) tensors_canonized(tensors_.size(), nullptr);

      container::set<size_t> colors;
      container::multimap<size_t, std::pair<size_t, size_t>>
          color2idx;  // maps color to the ordinals of the corresponding
      // tensors in tensors_ + their canonical ordinals given by cl
      // collect colors and tensors sorted by colors
      size_t vtx_cnt = 0;
      size_t tensor_cnt = 0;
      for (auto &&type : vtypes) {
        if (type == VertexType::TensorCore) {  // tensor core vertices were
                                               // created in the order of their
                                               // appearance in tensors_
          auto color = vcolors[vtx_cnt];
          if (colors.find(color) == colors.end()) colors.insert(color);
          color2idx.emplace(color, std::make_pair(tensor_cnt, cl[vtx_cnt]));
          ++tensor_cnt;
        }
        ++vtx_cnt;
      }
      // for each color sort tensors by canonical order
      // this assumes that tensors of different colors always commute
      // (reasonable) this only reorders tensors if they are c-numbers!
      container::svector<std::pair<size_t, size_t>>
          ord_can;  // canonically-ordered list of {ordinal in tensors_,
                    // canonical ordinal}
      container::svector<size_t>
          ord_orig;  // originaly-ordered list of {canonical ordinal}
      for (auto &&color : colors) {
        auto beg = color2idx.lower_bound(color);
        auto end = color2idx.upper_bound(color);
        const auto sz = end - beg;
        assert(sz > 0);
        if (sz > 1) {
          // all tensors of same color are c-numbers or all q-numbers ...
          // inspect the first to determine the type
          const bool cnumber = is_cnumber(*(tensors_.at(beg->second.first)));

          ord_can.resize(sz);
          ord_orig.resize(sz);

          size_t cnt = 0;
          for (auto it = beg; it != end; ++it, ++cnt) {
            ord_can[cnt] = it->second;
            ord_orig[cnt] = it->second.first;
            // assert that all tensors of same color are all c-numbers or all
            // q-numbers
            assert(cnumber == is_cnumber(*(tensors_.at(ord_orig[cnt]))));
          }
          using std::begin;
          using std::end;
          if (cnumber)  // only resort if these are cnumbers
            std::sort(begin(ord_can), end(ord_can),
                      [](const std::pair<size_t, size_t> &a,
                         const std::pair<size_t, size_t> &b) {
                        return a.second < b.second;
                      });
          std::sort(begin(ord_orig), end(ord_orig));
          // write (potentially reordered) tensors to tensors_canonized
          for (std::ptrdiff_t t = 0; t != sz; ++t) {
            tensors_canonized.at(ord_orig[t]) = tensors_.at(ord_can[t].first);
          }

        } else {  // sz = 1
          auto tidx = beg->second.first;
          tensors_canonized.at(tidx) = tensors_.at(tidx);
        }
      }  // colors

      // commit the canonically-ordered list of tensors to tensors_
      using std::swap;
      swap(tensors_canonized, tensors_);

    }  // tensors canonizing

  } else {  // fast approach uses heuristic canonization
    // simpler approach that will work perfectly as long as tensors are
    // distinguishable

    // - reindex anonymous indices using ordering of Edge as the
    // canonical definition of the anonymous index list
    {
      // resort edges_ first by index's character (named<anonymous), then by
      // Edge (not by Index's full label) ... this automatically puts named
      // indices first
      idx_terminals_sorted.resize(edges_.size());
      std::partial_sort_copy(
          begin(edges_), end(edges_), begin(idx_terminals_sorted),
          end(idx_terminals_sorted),
          [&namedness](const Edge &edge1, const Edge &edge2) {
            const auto n1 =
                namedness(edge1.idx());  // 1 -> named, 0 -> anonymous
            const auto n2 = namedness(edge2.idx());
            if (n1 == n2)
              return edge1 < edge2;
            else
              return n1 > n2;
          });

      // make index replacement list for anonymous indices only
      const auto num_named_indices = named_indices.size();
      std::for_each(
          begin(idx_terminals_sorted) + num_named_indices,
          end(idx_terminals_sorted),
          [&idxrepl, &idxfac, &is_anonymous_index](const auto &terminals) {
            const auto &idx = terminals.idx();
            assert(is_anonymous_index(
                idx));  // should only encounter anonymous indices here
            idxrepl.emplace(std::make_pair(idx, idxfac.make(idx)));
          });
    }
  }  // canonical index replacement list computed

  if (Logger::instance().canonicalize) {
    for (const auto &idxpair : idxrepl) {
      std::wcout << "TensorNetwork::canonicalize(" << (fast ? "fast" : "slow")
                 << "): replacing " << to_latex(idxpair.first) << " with "
                 << to_latex(idxpair.second) << std::endl;
    }
  }

#ifndef NDEBUG
  // assert that tensors' indices are not tagged since going to tag indices
  {
    for (const auto &tensor : tensors_) {
      assert(ranges::none_of(braket(*tensor), [](const Index &idx) {
        return idx.tag().has_value();
      }));
    }
  }
#endif
  bool pass_mutated = false;
  [[maybe_unused]] bool mutated = false;
  do {
    pass_mutated = false;
    for (auto &tensor : tensors_) {
      pass_mutated |= transform_indices(*tensor, idxrepl);
    }
    mutated |= pass_mutated;
  } while (pass_mutated);  // transform till stops changing

  // untag transformed indices (if any)
  {
    for (auto &tensor : tensors_) {
      reset_tags(*tensor);
    }
  }

  // - re-canonize tensors
  {
    // override the default canonicalizer
    DefaultTensorCanonicalizer default_tensor_canonizer(named_indices);
    for (auto &tensor : tensors_) {
      auto nondefault_canonizer_ptr =
          TensorCanonicalizer::nondefault_instance_ptr(tensor->_label());
      TensorCanonicalizer *tensor_canonizer =
          nondefault_canonizer_ptr ? nondefault_canonizer_ptr.get()
                                   : &default_tensor_canonizer;
      auto bp = tensor_canonizer->apply(*tensor);
      if (bp) *canon_byproduct *= *bp;
    }
  }
  edges_.clear();
  ext_indices_.clear();

  assert(canon_byproduct->is<Constant>());
  return (canon_byproduct->as<Constant>().value() == 1) ? nullptr
                                                        : canon_byproduct;
}

TensorNetwork::GraphData TensorNetwork::make_bliss_graph(
    const TensorNetwork::BlissGraphOptions &options) const {
  auto make_texlabel = [&](const auto &t) {
    using T = std::remove_cvref_t<decltype(t)>;
    std::wstring result;
    if constexpr (std::is_same_v<T, Index>) {
      result = L"$" + to_latex(t) + L"$";
    }
    return result;
  };

  // must call init_edges() prior to calling this
  if (edges_.empty()) {
    init_edges();
  }

  // initialize named_indices by default to all external indices (these HAVE
  // been computed in init_edges)
  const auto &named_indices = options.named_indices == nullptr
                                  ? this->ext_indices()
                                  : *options.named_indices;

  VertexPainter colorizer(named_indices, options.distinct_named_indices);

  // results
  std::shared_ptr<bliss::Graph> graph;
  const auto nidx = edges_.size() + pure_proto_indices_.size();
  std::vector<std::wstring> vertex_labels(
      options.make_labels ? nidx : 0);  // the size will be updated
  std::vector<std::optional<std::wstring>> vertex_texlabels(
      options.make_texlabels ? nidx : 0);  // the size will be updated
  std::vector<GraphData::VertexColor> vertex_color(
      nidx, 0);                               // the size will be updated
  std::vector<VertexType> vertex_type(nidx);  // the size will be updated

  // compute # of vertices
  size_t nv = 0;
  size_t index_cnt = 0;
  size_t spbundle_cnt = 0;
  // first count vertex indices ... the only complication are symmetric
  // protoindex bundles this will keep track of unique symmetric protoindex
  // bundles
  using protoindex_bundle_t =
      std::decay_t<decltype(std::declval<const Index &>().proto_indices())>;
  container::set<protoindex_bundle_t> symmetric_protoindex_bundles;
  const size_t spbundle_vertex_offset =
      nidx;  // where spbundle vertices will start

  // iterates over edges_, then pure_proto_indices_, this is equivalent to
  // iteration over grand_index_list
  ranges::for_each(edges_, [&](const Edge &edge) {
    const Index &idx = edge.idx();
    ++nv;  // each index is a vertex
    if (options.make_labels) {
      vertex_labels.at(index_cnt) = idx.full_label();
    }
    if (options.make_texlabels) {
      vertex_texlabels.at(index_cnt) = make_texlabel(idx);
    }
    vertex_type.at(index_cnt) = VertexType::Index;
    vertex_color.at(index_cnt) = colorizer(idx);

    // each symmetric proto index bundle will have a vertex ...
    // for now only store the unique protoindex bundles in
    // symmetric_protoindex_bundles, then commit their data to
    // vertex_{labels,type,color} later
    if (idx.has_proto_indices()) {
      assert(idx.symmetric_proto_indices());  // only symmetric protoindices are
                                              // supported right now
      if (symmetric_protoindex_bundles.find(idx.proto_indices()) ==
          symmetric_protoindex_bundles
              .end()) {  // new bundle? make a vertex for it
        auto graph = symmetric_protoindex_bundles.insert(idx.proto_indices());
        assert(graph.second);
      }
    }
    index_cnt++;
  });
  ranges::for_each(pure_proto_indices_, [&](const Index &idx) {
    ++nv;  // each index is a vertex
    if (options.make_labels) {
      vertex_labels.at(index_cnt) = idx.full_label();
    }
    if (options.make_texlabels) {
      vertex_texlabels.at(index_cnt) = make_texlabel(idx);
    }
    vertex_type.at(index_cnt) = VertexType::Index;
    vertex_color.at(index_cnt) = colorizer(idx);
    index_cnt++;
  });

  // now commit protoindex bundle metadata
  ranges::for_each(symmetric_protoindex_bundles, [&](const auto &bundle) {
    assert(!bundle.empty());
    ++nv;  // each symmetric protoindex bundle is a vertex
    if (options.make_labels) {
      std::wstring spbundle_label = L"<";
      const auto end = bundle.end();
      auto it = bundle.begin();
      spbundle_label += it->full_label();
      for (++it; it != end; ++it) {
        spbundle_label += L",";
        spbundle_label += it->full_label();
      }
      spbundle_label += L">";
      vertex_labels.emplace_back(std::move(spbundle_label));
    }
    if (options.make_texlabels) {
      std::wstring spbundle_texlabel = L"$\\langle";
      const auto end = bundle.end();
      auto it = bundle.begin();
      spbundle_texlabel += it->to_latex();
      for (++it; it != end; ++it) {
        spbundle_texlabel += L",";
        spbundle_texlabel += it->to_latex();
      }
      spbundle_texlabel += L"\\rangle$";
      vertex_texlabels.emplace_back(std::move(spbundle_texlabel));
    }
    vertex_type.push_back(VertexType::SPBundle);
    vertex_color.push_back(colorizer(bundle));
    spbundle_cnt++;
  });

  // now account for vertex representation of tensors
  size_t tensor_cnt = 0;
  // this will map to tensor index to the first (core) vertex in its
  // representation
  container::svector<size_t> tensor_vertex_offset(tensors_.size());
  ranges::for_each(tensors_, [&](const auto &t) {
    tensor_vertex_offset.at(tensor_cnt) = nv;
    // each tensor has a core vertex (to be colored by its label)
    ++nv;
    std::wstring tlabel;
    if (options.make_labels || options.make_texlabels) tlabel = label(*t);
    if (options.make_labels) {
      vertex_labels.emplace_back(tlabel);
    }
    if (options.make_texlabels) {
      vertex_texlabels.emplace_back(L"$" + utf_to_latex(tlabel) + L"$");
    }
    vertex_type.emplace_back(VertexType::TensorCore);
    vertex_color.push_back(colorizer(*t));

    // symmetric/antisymmetric tensors are represented by 3 more vertices:
    // - bra
    // - ket
    // - braket (connecting bra and ket to the core)
    auto &tref = *t;
    if (symmetry(tref) != Symmetry::nonsymm) {
      nv += 3;
      if (options.make_labels) {
        vertex_labels.emplace_back(
            std::wstring(L"bra") + to_wstring(bra_rank(tref)) +
            ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      }
      if (options.make_texlabels) {
        vertex_texlabels.emplace_back(std::nullopt);
      }
      vertex_type.push_back(VertexType::TensorBra);
      vertex_color.push_back(colorizer(BraGroup{0}));
      if (options.make_labels) {
        vertex_labels.emplace_back(
            std::wstring(L"ket") + to_wstring(ket_rank(tref)) +
            ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      }
      if (options.make_texlabels) {
        vertex_texlabels.emplace_back(std::nullopt);
      }
      vertex_type.push_back(VertexType::TensorKet);
      if (braket_symmetry(tref) == BraKetSymmetry::symm) {
        // Use BraGroup for kets as well as they are supposed to be
        // indistinguishable
        vertex_color.push_back(colorizer(BraGroup{0}));
      } else {
        vertex_color.push_back(colorizer(KetGroup{0}));
      }
      if (options.make_labels) {
        vertex_labels.emplace_back(
            std::wstring(L"bk") +
            ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      }
      if (options.make_texlabels) {
        vertex_texlabels.emplace_back(std::nullopt);
      }
      vertex_type.push_back(VertexType::Particle);
      // Color bk node in same color as tensor core
      vertex_color.push_back(colorizer(tref));
    }
    // nonsymmetric tensors are represented by 3*rank more vertices (with rank =
    // max(bra_rank(),ket_rank())
    else {
      const auto rank = std::max(bra_rank(tref), ket_rank(tref));
      assert(rank <= max_rank);
      for (size_t p = 0; p != rank; ++p) {
        nv += 3;
        std::wstring pstr;
        if (options.make_labels) pstr = to_wstring(p + 1);
        if (options.make_labels) {
          vertex_labels.emplace_back(std::wstring(L"bra") + pstr);
        }
        if (options.make_texlabels) {
          vertex_texlabels.emplace_back(std::nullopt);
        }
        vertex_type.push_back(VertexType::TensorBra);
        const bool distinguishable_particles =
            particle_symmetry(tref) == ParticleSymmetry::nonsymm;
        vertex_color.push_back(
            colorizer(BraGroup{distinguishable_particles ? p : 0}));
        if (options.make_labels) {
          vertex_labels.emplace_back(std::wstring(L"ket") + pstr);
        }
        if (options.make_texlabels) {
          vertex_texlabels.emplace_back(std::nullopt);
        }
        vertex_type.push_back(VertexType::TensorKet);
        if (braket_symmetry(tref) == BraKetSymmetry::symm) {
          // Use BraGroup for kets as well as they are supposed to be
          // indistinguishable
          vertex_color.push_back(
              colorizer(BraGroup{distinguishable_particles ? p : 0}));
        } else {
          vertex_color.push_back(
              colorizer(KetGroup{distinguishable_particles ? p : 0}));
        }
        if (options.make_labels) {
          vertex_labels.emplace_back(std::wstring(L"bk") + pstr);
        }
        if (options.make_texlabels) {
          vertex_texlabels.emplace_back(std::nullopt);
        }
        vertex_type.push_back(VertexType::Particle);
        vertex_color.push_back(colorizer(tref));
      }
    }
    // aux indices currently do not support any symmetry
    assert(aux_rank(tref) <= max_rank);
    for (size_t p = 0; p != aux_rank(tref); ++p) {
      nv += 1;
      if (options.make_labels) {
        auto pstr = to_wstring(p + 1);
        vertex_labels.push_back(std::wstring(L"aux") + pstr);
      }
      if (options.make_texlabels) {
        vertex_texlabels.emplace_back(std::nullopt);
      }
      vertex_type.push_back(VertexType::TensorAux);
      vertex_color.push_back(colorizer(AuxGroup{p}));
    }

    ++tensor_cnt;
  });

  // allocate graph
  graph = std::make_shared<bliss::Graph>(nv);

  // add edges
  // - each index's degree <= 2 + # of protoindex terminals
  index_cnt = 0;
  ranges::for_each(edges_, [&](const Edge &edge) {
    assert(edge.size() > 0);
    [[maybe_unused]] const auto edge_connected = edge.size() == 2;
    for (int t = 0; t != edge.size(); ++t) {
      const auto &terminal = edge[t];
      const auto tensor_ord = terminal.tensor_ord;
      // which vertex after the core vertex in the tensor does this edge attach
      // to?
      int vertex_ord;
      if (terminal.slot_type != TensorIndexSlotType::Aux) {
        // each vector slot group is represented by 3 vertices: bra, ket, bk
        const auto slot_group_ord_offset = 3 * terminal.slot_group_ord;
        vertex_ord =
            slot_group_ord_offset + static_cast<int>(terminal.slot_type);
      } else {
        const auto &tensor = *tensors_[tensor_ord];
        const auto num_vector_slot_groups =
            tensor._symmetry() == Symmetry::nonsymm
                ? std::max(tensor._bra().size(), tensor._ket().size())
                : 1;
        // aux slot groups appear after the vector slot groups
        const auto slot_group_ord_offset =
            num_vector_slot_groups * 3 +
            (terminal.slot_group_ord - num_vector_slot_groups);
        // each aux slot group is represented by 1 vertex
        vertex_ord = slot_group_ord_offset;
      }
      const size_t slot_vertex_index = tensor_vertex_offset[tensor_ord] +
                                       /* core */ 1 + vertex_ord;
      graph->add_edge(index_cnt, slot_vertex_index);
    }
    // if this index has symmetric protoindex bundles
    const auto &idx = edge.idx();
    if (idx.has_proto_indices()) {
      if (idx.symmetric_proto_indices()) {
        assert(symmetric_protoindex_bundles.find(idx.proto_indices()) !=
               symmetric_protoindex_bundles.end());
        const auto spbundle_idx =
            symmetric_protoindex_bundles.find(idx.proto_indices()) -
            symmetric_protoindex_bundles.begin();
        graph->add_edge(index_cnt, spbundle_vertex_offset + spbundle_idx);
      } else {
        abort();  // nonsymmetric proto indices not supported yet
      }
    }
    ++index_cnt;
  });
  // - link up proto indices, if any ... only symmetric protobundles are
  // supported now
  spbundle_cnt = spbundle_vertex_offset;
  ranges::for_each(symmetric_protoindex_bundles, [&graph, this, &spbundle_cnt](
                                                     const auto &bundle) {
    for (auto &&proto_index : bundle) {
      // proto index either connects tensors (i.e. it's in edges_) OR
      // it's among pure_proto_indices_
      auto edges_it = edges_.find(proto_index.full_label());
      if (edges_it != edges_.end()) {
        const auto proto_index_vertex = edges_it - edges_.begin();
        graph->add_edge(spbundle_cnt, proto_index_vertex);
      } else {
        auto ppidx_it = pure_proto_indices_.find(proto_index);
        assert(ppidx_it != pure_proto_indices_.end());
        const auto proto_index_vertex =
            ppidx_it - pure_proto_indices_.begin() + edges_.size();
        graph->add_edge(spbundle_cnt, proto_index_vertex);
      }
    }
    ++spbundle_cnt;
  });
  // - link up tensors
  tensor_cnt = 0;
  ranges::for_each(
      tensors_, [&graph, &tensor_cnt, &tensor_vertex_offset](const auto &t) {
        const auto vertex_offset = tensor_vertex_offset.at(tensor_cnt);
        // for each braket terminal linker
        auto &tref = *t;
        const size_t nbk = symmetry(tref) == Symmetry::nonsymm
                               ? std::max(bra_rank(tref), ket_rank(tref))
                               : 1;
        for (size_t bk = 1; bk <= nbk; ++bk) {
          const int bk_vertex = vertex_offset + 3 * bk;
          graph->add_edge(vertex_offset, bk_vertex);  // core
          graph->add_edge(bk_vertex - 2, bk_vertex);  // bra
          graph->add_edge(bk_vertex - 1, bk_vertex);  // ket
        }
        // for each aux terminal linker
        const size_t naux = aux_rank(tref);
        for (size_t aux = 1; aux <= naux; ++aux) {
          const int aux_vertex = vertex_offset + 3 * nbk + aux;
          graph->add_edge(vertex_offset, aux_vertex);  // core
        }
        ++tensor_cnt;
      });

  for (const auto [vertex, color] : ranges::views::enumerate(vertex_color)) {
    graph->change_color(vertex, color);
  }

  return {.graph = std::move(graph),
          .vertex_labels = std::move(vertex_labels),
          .vertex_texlabels = std::move(vertex_texlabels),
          .vertex_colors = std::move(vertex_color),
          .vertex_types = std::move(vertex_type)};
}

void TensorNetwork::init_edges() const {
  if (have_edges_) return;

  auto idx_insert = [this](const Index &idx, int tensor_idx,
                           TensorIndexSlotType slot_type, int slot_group_ord) {
    if (Logger::instance().tensor_network) {
      std::wcout << "TensorNetwork::init_edges: idx=" << to_latex(idx)
                 << " attached to tensor " << std::abs(tensor_idx) << "'s "
                 << ((tensor_idx > 0) ? "bra" : "ket") << " via slot group "
                 << slot_group_ord << std::endl;
    }
    edges_t &edges = this->edges_;
    auto it = edges.find(idx.full_label());
    if (it == edges.end()) {
      edges.emplace(Edge::Terminal(tensor_idx, slot_type, slot_group_ord),
                    &idx);
    } else {
      it->connect_to(Edge::Terminal(tensor_idx, slot_type, slot_group_ord));
    }
  };

  int t_idx = 0;
  for (auto &&t : tensors_) {
    const auto t_is_nonsymm = symmetry(*t) == Symmetry::nonsymm;
    size_t slot_group_ord = 0;
    for (const Index &idx : t->_bra()) {
      idx_insert(idx, t_idx, TensorIndexSlotType::Bra,
                 t_is_nonsymm ? slot_group_ord++
                              : /* bra indices of symm/antisymm tensors are part
                                   of same slot group */
                     0);
    }
    slot_group_ord = 0;  // bra and ket slots are grouped together
    for (const Index &idx : t->_ket()) {
      idx_insert(idx, t_idx, TensorIndexSlotType::Ket,
                 t_is_nonsymm ? slot_group_ord++
                              : /* ket indices of symm/antisymm tensors are part
                                   of same slot group */
                     0);
    }
    // aux slot group count starts after last braket slot group
    for (const Index &idx : t->_aux()) {
      // aux slots are not symmetric
      idx_insert(idx, t_idx, TensorIndexSlotType::Aux, slot_group_ord++);
    }
    ++t_idx;
  }

  // extract external indices and all protoindices (since some external indices
  // may be pure protoindices)
  named_indices_t proto_indices;
  for (const auto &terminals : edges_) {
    assert(terminals.size() != 0);
    // External index (== Edge only connected to a single vertex in the
    // network)
    if (terminals.size() == 1) {
      if (Logger::instance().tensor_network) {
        std::wcout << "idx " << to_latex(terminals.idx()) << " is external"
                   << std::endl;
      }
      const auto &[it, inserted] = ext_indices_.emplace(terminals.idx());
      // only scenario where idx is already in ext_indices_ if it were a
      // protoindex of a previously inserted ext index ... check to ensure no
      // accidental duplicates
      if (!inserted) {
        assert(proto_indices.contains(terminals.idx()));
      }
    }

    // add proto indices to the grand list of proto indices
    for (auto &&proto_idx : terminals.idx().proto_indices()) {
      if (proto_idx.has_proto_indices())
        throw std::runtime_error(
            "TensorNetwork does not support recursive protoindices");  // for
                                                                       // now no
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

  // some external indices will have protoindices that are NOT among
  // pure_proto_indices_, e.g.
  // i2 in f_i2^{a2^{i1,i2}} t_{a2^{i1,i2}a3^{i1,i2}}^{i2,i1}
  // is not added to ext_indices_ due to being among doubly-connected edges_
  // and thus is not among pure_proto_indices_, but it needs to be
  // and external index due to a3^{i1,i2} being an external index
  named_indices_t ext_proto_indices;
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

container::svector<std::pair<long, long>> TensorNetwork::factorize() {
  abort();  // not yet implemented
}

size_t TensorNetwork::SlotCanonicalizationMetadata::hash_value() const {
  return graph->get_hash();
}

TensorNetwork::SlotCanonicalizationMetadata TensorNetwork::canonicalize_slots(
    const std::vector<std::wstring> &cardinal_tensor_labels,
    const TensorNetwork::named_indices_t *named_indices_ptr,
    TensorNetwork::SlotCanonicalizationMetadata::named_index_compare_t
        named_index_compare) {
  if (!named_index_compare)
    named_index_compare = [](const auto &idxptr_slottype_1,
                             const auto &idxptr_slottype_2) -> bool {
      const auto &[idxptr1, slottype1] = idxptr_slottype_1;
      const auto &[idxptr2, slottype2] = idxptr_slottype_2;
      return idxptr1->space() < idxptr2->space();
    };

  TensorNetwork::SlotCanonicalizationMetadata metadata;

  if (Logger::instance().canonicalize) {
    std::wcout << "TensorNetwork::canonicalize_slots(): input tensors\n";
    size_t cnt = 0;
    ranges::for_each(tensors_, [&](const auto &t) {
      std::wcout << "tensor " << cnt++ << ": " << to_latex(*t) << std::endl;
    });
    std::wcout << "cardinal_tensor_labels = ";
    ranges::for_each(cardinal_tensor_labels,
                     [](auto &&i) { std::wcout << i << L" "; });
    std::wcout << std::endl;
  }

  if (edges_.empty()) init_edges();

  // initialize named_indices by default to all external indices (these HAVE
  // been computed in init_edges)
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
  auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
      make_bliss_graph({.named_indices = &named_indices,
                        .distinct_named_indices = false,
                        .make_labels = Logger::instance().canonicalize_dot,
                        .make_texlabels = Logger::instance().canonicalize_dot});
  if (Logger::instance().canonicalize_dot) {
    std::wcout << "Input graph for canonicalization:\n";
    graph->write_dot(std::wcout, vlabels, vtexlabels);
  }

  // canonize the graph
  bliss::Stats stats;
  graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
  const unsigned int *cl = graph->canonical_form(stats, nullptr, nullptr);

  metadata.graph = std::shared_ptr<bliss::Graph>(graph->permute(cl));

  if (Logger::instance().canonicalize_dot) {
    auto permute = [](const auto &vector, auto &&perm) {
      using std::size;
      auto sz = size(vector);
      std::decay_t<decltype(vector)> pvector(sz);
      for (size_t i = 0; i != sz; ++i) pvector[perm[i]] = vector[i];
      return pvector;
    };

    auto cvlabels = permute(vlabels, cl);
    auto cvtexlabels = permute(vtexlabels, cl);
    std::wcout << "Canonicalized graph:\n";
    metadata.graph->write_dot(std::wcout, cvlabels, cvtexlabels);
  }

  // produce named indices sorted by named_index_compare first, then by
  // canonical order produced by bliss
  {
    using ord_cord_it_t =
        std::tuple<size_t, size_t, named_indices_t::const_iterator>;
    using cord_set_t = container::set<ord_cord_it_t, detail::tuple_less<1>>;

    auto grand_index_list = ranges::views::concat(
        edges_ | ranges::views::transform(edge2index<Edge>),
        pure_proto_indices_);

    // for each named index type (as defined by named_index_compare) maps its
    // ptr in grand_index_list to its ordinal in grand_index_list + canonical
    // ordinal + its iterator in metadata.named_indices
    container::map<std::pair<const Index *, IndexSlotType>, cord_set_t,
                   decltype(named_index_compare)>
        idx2cord(named_index_compare);
    // collect named indices and sort on the fly
    size_t idx_ord = 0;
    auto grand_index_list_end = grand_index_list.end();
    for (auto git = grand_index_list.begin(); git != grand_index_list_end;
         ++git) {
      const auto &idx = *git;
      if (is_named_index(idx)) {
        const auto named_indices_it = metadata.named_indices.find(idx);
        assert(named_indices_it != metadata.named_indices.end());

        // deduce the slot type occupied by this index
        IndexSlotType slot_type;
        if (idx_ord < edges_.size()) {
          auto edge_it = edges_.begin();
          std::advance(edge_it, idx_ord);
          // there are 2 possibilities: its index edge is disconnected or
          // connected ... the latter would only occur if this index is named
          // due to also being a protoindex on one of the named indices!
          if (edge_it->size() == 1) {
            if (edge_it->second().slot_type == TensorIndexSlotType::Aux)
              slot_type = IndexSlotType::TensorAux;
            else if (edge_it->second().slot_type == TensorIndexSlotType::Bra)
              slot_type = IndexSlotType::TensorBra;
            else  // edge_it->second().slot_type == TensorIndexSlotType::Ket
              slot_type = IndexSlotType::TensorKet;
          } else {  // if
            assert(edge_it->size() == 2);
            slot_type = IndexSlotType::SPBundle;
          }
        } else
          slot_type = IndexSlotType::SPBundle;

        // find the entry for this index type
        const auto idxptr_slottype = std::make_pair(&idx, slot_type);
        auto it = idx2cord.find(idxptr_slottype);
        if (it == idx2cord.end()) {
          bool inserted;
          std::tie(it, inserted) = idx2cord.insert(std::make_pair(
              idxptr_slottype, cord_set_t(cord_set_t::key_compare{})));
          assert(inserted);
        }
        it->second.emplace(idx_ord, cl[idx_ord], named_indices_it);
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

  // compute the phase associated with *slot* canonicalization ...
  // - choose an index-independent order of slots: loop over bra then ket then
  // aux slots of each tensor ... record each Index in the order of its
  // appearance when iterating over slots. This is the input ordinal of this
  // Index. NB We can't just use Index::full_label() to look up Index in edges_
  // since this would produce label-dependent ordinals
  // - canonicalization reorders slots according to the topology-based order of
  // indices.
  // - for each symmetric/antisymmetric bra/ket bundle canonical
  // order of slots is in lexicographic order of the canonical
  // ordinals of the edges connected to them, such slot reordering
  // induces a phase for antisymmetric bundles. determine the phase by
  // computing the parity of the canonical ordinal sequence
  {
    metadata.phase = 1;

    // iterate over tensor bra/ket/aux index slots in the canonical input order,
    // for each antisymmetric bra/ket compute phase
    container::map<Index, std::size_t, FullLabelCompare>
        idx_inord;  // Index -> input ordinal; helps with computing the ordinals
                    // during the traversal
    for (auto &_t : tensors_) {
      assert(std::dynamic_pointer_cast<Tensor>(_t));
      auto t = std::static_pointer_cast<Tensor>(_t);

      // returns an iterator to {Index,inord} pair
      auto index_inord_it = [&](const Index &idx) {
        auto it = idx_inord.find(idx.full_label());
        if (it == idx_inord.end()) {
          const auto inord = idx_inord.size();
          bool inserted;
          std::tie(it, inserted) = idx_inord.emplace(idx, inord);
          assert(inserted);
        }
        return it;
      };

      // computes parity of a bra/ket bundle due to the reordering of its slots
      // involved in the TN canonicalization
      auto input_to_canonical_parity = [&](const auto &idx_rng) {
        using ranges::size;
        const auto sz = size(idx_rng);
        if (sz < 2) {  // no phase for 1-index bundles, but still process the
                       // indices to ensure ordinals are correct
          using ranges::begin;
          index_inord_it(*begin(idx_rng));
          return 1;
        }

        auto parity = [&](auto &rng) {
          reset_ts_swap_counter<std::size_t>();
          using ranges::begin;
          using ranges::end;
          bubble_sort(begin(rng), end(rng));
          return ts_swap_counter_is_even<std::size_t>() ? +1 : -1;
        };

        container::vector<SwapCountable<std::size_t>> ordinals_canonical;
        ordinals_canonical.reserve(sz);
        for (const auto &idx : idx_rng) {
          // use canonical ordinal of the index vertex as the canonical index
          // ordinal, to find it need the ordinal of the corresponding vertex in
          // the input graph
          auto input_vertex_it = this->edges_.find(idx.full_label());
          assert(input_vertex_it != this->edges_.end());
          const auto inord_vertex = input_vertex_it - this->edges_.begin();
          const auto canord_vertex = cl[inord_vertex];

          // to which edge did this get mapped?
          //          std::wcout << "vlabels[ord_original=" << inord_vertex
          //                     << "]=" << vlabels[inord_vertex]
          //                     << " -> ord_canonical=" << canord_vertex <<
          //                     std::endl;
          ordinals_canonical.emplace_back(canord_vertex);
        }

        const auto parity_canonical = parity(ordinals_canonical);
        return parity_canonical;
      };

      // canonical order of slots: bra, then ket, then aux. We only care about
      // indices in tensor slots here, so no need to worry about protoindex
      // slots
      if (t->symmetry() == Symmetry::antisymm) {
        // bra first, then ket
        metadata.phase *= input_to_canonical_parity(t->bra());
        metadata.phase *= input_to_canonical_parity(t->ket());
      } else {  // although we don't need to worry about phases, we still need
                // to process all indices so that the input ordinals are correct
        for (auto &&idx : t->bra()) {
          index_inord_it(idx);
        }
        for (auto &&idx : t->ket()) {
          index_inord_it(idx);
        }
      }
      for (auto &&idx : t->aux()) {
        index_inord_it(idx);
      }
    }
  }

  return metadata;
}

}  // namespace sequant
