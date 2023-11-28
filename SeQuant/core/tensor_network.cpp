//
// Created by Eduard Valeyev on 2019-02-26.
//

#include "tensor_network.hpp"
#include "bliss.hpp"
#include "logger.hpp"

namespace sequant {

ExprPtr TensorNetwork::canonicalize(
    const container::vector<std::wstring> &cardinal_tensor_labels, bool fast,
    const named_indices_t *named_indices_ptr) {
  ExprPtr canon_biproduct = ex<Constant>(1);
  container::svector<Edge> idx_terminals_sorted;  // to avoid memory allocs

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

  // - resort tensors (this cannot be done in Product::canonicalize since that
  // requires analysis of commutativity ... here we are dealing with tensors
  // only and are free to reorder)
  using std::begin;
  using std::end;
  bubble_sort(
      begin(tensors_), end(tensors_),
      [&cardinal_tensor_labels](const auto &first_ptr, const auto &second_ptr) {
        const auto &first = *first_ptr;
        const auto &second = *second_ptr;
        // grab base label if adjoint label is present
        auto base_label = [](const auto &t) {
          if (label(t).back() == adjoint_label) {
            return label(t).substr(0, label(t).size() - 1);
          } else {
            return label(t);
          }
        };
        // tensors commute if their colors are different or either one of them
        // is a c-number
        if ((color(first) != color(second)) || is_cnumber(first) ||
            is_cnumber(second)) {
          const auto cardinal_tensor_labels_end = end(cardinal_tensor_labels);
          const auto first_cardinal_it =
              std::find(begin(cardinal_tensor_labels),
                        end(cardinal_tensor_labels), base_label(first));
          const auto second_cardinal_it =
              std::find(begin(cardinal_tensor_labels),
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

  if (Logger::get_instance().canonicalize) {
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

    auto permute = [](const auto &vector, auto &&perm) {
      using std::size;
      auto sz = size(vector);
      std::decay_t<decltype(vector)> pvector(sz);
      for (size_t i = 0; i != sz; ++i) pvector[perm[i]] = vector[i];
      return pvector;
    };

    // make the graph
    auto [graph, vlabels, vcolors, vtypes] = make_bliss_graph(&named_indices);
    //    graph->write_dot(std::wcout, vlabels);

    // canonize the graph
    bliss::Stats stats;
    graph->set_splitting_heuristic(bliss::Graph::shs_fsm);
    const unsigned int *cl = graph->canonical_form(stats, nullptr, nullptr);

    if (Logger::get_instance().canonicalize_dot) {
      bliss::Graph *cgraph = graph->permute(cl);
      auto cvlabels = permute(vlabels, cl);
      cgraph->write_dot(std::wcout, cvlabels);
      delete cgraph;
    }

    // make anonymous index replacement list
    {
      // for each color make a replacement list for bringing the indices to
      // the canonical order
      container::set<size_t> colors;
      container::multimap<size_t, std::pair<size_t, size_t>>
          color2idx;  // maps color to the ordinals of the corresponding
      // indices in edges_ + their canonical ordinals
      // collect colors and anonymous indices sorted by colors
      size_t idx_cnt = 0;
      for (auto &&ttpair : edges_) {
        auto color = vcolors[idx_cnt];
        if (colors.find(color) == colors.end()) colors.insert(color);
        if (is_anonymous_index(ttpair.idx())) {
          color2idx.emplace(color, std::make_pair(idx_cnt, cl[idx_cnt]));
        }
        ++idx_cnt;
      }
      // for each color sort anonymous indices by canonical order
      container::svector<std::pair<size_t, size_t>>
          idx_can;  // canonically-ordered list of {index ordinal in edges_,
                    // canonical ordinal}
      for (auto &&color : colors) {
        auto beg = color2idx.lower_bound(color);
        auto end = color2idx.upper_bound(color);
        const auto sz = end - beg;
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
          for (auto &&p : idx_can) {
            const auto &idx = (edges_.begin() + p.first)->idx();
            idxrepl.emplace(std::make_pair(idx, idxfac.make(idx)));
          }
        } else if (sz == 1) {  // no need for resorting of colors with 1 index
                               // only, but still need to replace the index
          const auto edge_it = edges_.begin() + beg->second.first;
          const auto &idx = edge_it->idx();
          idxrepl.emplace(std::make_pair(idx, idxfac.make(idx)));
        }
        // sz == 0 is possible since some colors in colors refer to tensors
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
          for (size_t t = 0; t != sz; ++t) {
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

  if (Logger::get_instance().canonicalize) {
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
  bool mutated = false;
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
      if (bp) *canon_biproduct *= *bp;
    }
  }
  edges_.clear();
  ext_indices_.clear();

  assert(canon_biproduct->is<Constant>());
  return (canon_biproduct->as<Constant>().value() == 1) ? nullptr
                                                        : canon_biproduct;
}

std::tuple<std::shared_ptr<bliss::Graph>, std::vector<std::wstring>,
           std::vector<std::size_t>,
           std::vector<typename TensorNetwork::VertexType>>
TensorNetwork::make_bliss_graph(
    const named_indices_t *named_indices_ptr) const {
  // must call init_edges() prior to calling this
  if (edges_.empty()) {
    init_edges();
  }

  // initialize named_indices by default to all external indices (these HAVE
  // been computed in init_edges)
  const auto &named_indices =
      named_indices_ptr == nullptr ? this->ext_indices() : *named_indices_ptr;

  // results
  std::shared_ptr<bliss::Graph> graph;
  std::vector<std::wstring> vertex_labels(
      edges_.size());  // the size will be updated
  std::vector<std::size_t> vertex_color(edges_.size(),
                                        0);  // the size will be updated
  std::vector<VertexType> vertex_type(
      edges_.size());  // the size will be updated

  // N.B. Colors [0, 2 max rank + named_indices.size()) are reserved:
  // 0 - the bra vertex (for particle 0, if bra is nonsymm, or for the entire
  // bra, if (anti)symm) 1 - the bra vertex for particle 1, if bra is nonsymm
  // ...
  // max_rank - the ket vertex (for particle 0, if particle-asymmetric, or for
  // the entire ket, if particle-symmetric) max_rank+1 - the ket vertex for
  // particle 1, if particle-asymmetric
  // ...
  // 2 max_rank - first named index
  // 2 max_rank + 1 - second named index
  // ...
  // N.B. For braket-symmetric tensors the ket vertices use the same indices as
  // the bra vertices
  auto nonreserved_color = [&named_indices](size_t color) -> bool {
    return color >= 2 * max_rank + named_indices.size();
  };

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
      edges_.size();  // where spbundle vertices will start
  ranges::for_each(edges_, [&](const Edge &ttpair) {
    const Index &idx = ttpair.idx();
    ++nv;  // each index is a vertex
    vertex_labels.at(index_cnt) = idx.to_latex();
    vertex_type.at(index_cnt) = VertexType::Index;
    // assign color: named indices use reserved colors
    const auto named_index_it = named_indices.find(idx);
    if (named_index_it ==
        named_indices.end()) {  // anonymous index? use Index::color
      const auto idx_color = idx.color();
      assert(nonreserved_color(idx_color));
      vertex_color.at(index_cnt) = idx_color;
    } else {
      const auto named_index_rank = named_index_it - named_indices.begin();
      vertex_color.at(index_cnt) = 2 * max_rank + named_index_rank;
    }
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
  // now commit protoindex bundle metadata
  ranges::for_each(symmetric_protoindex_bundles, [&](const auto &bundle) {
    ++nv;  // each symmetric protoindex bundle is a vertex
    std::wstring spbundle_label = L"{";
    for (auto &&pi : bundle) {
      spbundle_label += pi.to_latex();
    }
    spbundle_label += L"}";
    vertex_labels.push_back(spbundle_label);
    vertex_type.push_back(VertexType::SPBundle);
    const auto idx_proto_indices_color = Index::proto_indices_color(bundle);
    assert(nonreserved_color(idx_proto_indices_color));
    vertex_color.push_back(idx_proto_indices_color);
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
    const auto tlabel = label(*t);
    vertex_labels.emplace_back(tlabel);
    vertex_type.emplace_back(VertexType::TensorCore);
    const auto t_color = hash::value(tlabel);
    static_assert(sizeof(t_color) == sizeof(unsigned long int));
    assert(nonreserved_color(t_color));
    vertex_color.push_back(t_color);
    // symmetric/antisymmetric tensors are represented by 3 more vertices:
    // - bra
    // - ket
    // - braket (connecting bra and ket to the core)
    auto &tref = *t;
    if (symmetry(tref) != Symmetry::nonsymm) {
      nv += 3;
      vertex_labels.push_back(
          std::wstring(L"bra") + to_wstring(bra_rank(tref)) +
          ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      vertex_type.push_back(VertexType::TensorBra);
      vertex_color.push_back(0);
      vertex_labels.push_back(
          std::wstring(L"ket") + to_wstring(ket_rank(tref)) +
          ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      vertex_type.push_back(VertexType::TensorKet);
      vertex_color.push_back(
          braket_symmetry(tref) == BraKetSymmetry::symm ? 0 : max_rank);
      vertex_labels.push_back(
          std::wstring(L"bk") +
          ((symmetry(tref) == Symmetry::antisymm) ? L"a" : L"s"));
      vertex_type.push_back(VertexType::TensorBraKet);
      vertex_color.push_back(t_color);
    }
    // nonsymmetric tensors are represented by 3*rank more vertices (with rank =
    // max(bra_rank(),ket_rank())
    else {
      const auto rank = std::max(bra_rank(tref), ket_rank(tref));
      assert(rank <= max_rank);
      for (size_t p = 0; p != rank; ++p) {
        nv += 3;
        auto pstr = to_wstring(p + 1);
        vertex_labels.push_back(std::wstring(L"bra") + pstr);
        vertex_type.push_back(VertexType::TensorBra);
        const bool t_is_particle_symmetric =
            particle_symmetry(tref) == ParticleSymmetry::nonsymm;
        const auto bra_color = t_is_particle_symmetric ? p : 0;
        vertex_color.push_back(bra_color);
        vertex_labels.push_back(std::wstring(L"ket") + pstr);
        vertex_type.push_back(VertexType::TensorKet);
        vertex_color.push_back(braket_symmetry(tref) == BraKetSymmetry::symm
                                   ? bra_color
                                   : bra_color + max_rank);
        vertex_labels.push_back(std::wstring(L"bk") + pstr);
        vertex_type.push_back(VertexType::TensorBraKet);
        vertex_color.push_back(t_color);
      }
    }
    ++tensor_cnt;
  });

  // allocate graph
  graph = std::make_shared<bliss::Graph>(nv);

  // add edges
  // - each index's degree <= 2 + # of protoindex terminals
  index_cnt = 0;
  ranges::for_each(edges_, [&](const Edge &ttpair) {
    for (int t = 0; t != 2; ++t) {
      const auto terminal_index = t == 0 ? ttpair.first() : ttpair.second();
      const auto terminal_position =
          t == 0 ? ttpair.first_position() : ttpair.second_position();
      if (terminal_index) {
        const auto tidx = std::abs(terminal_index) - 1;
        const auto ttpos = terminal_position;
        const bool bra = terminal_index > 0;
        const size_t braket_vertex_index = tensor_vertex_offset[tidx] +
                                           /* core */ 1 + 3 * ttpos +
                                           (bra ? 0 : 1);
        graph->add_edge(index_cnt, braket_vertex_index);
      }
    }
    // if this index has symmetric protoindex bundles
    if (ttpair.idx().has_proto_indices()) {
      if (ttpair.idx().symmetric_proto_indices()) {
        assert(
            symmetric_protoindex_bundles.find(ttpair.idx().proto_indices()) !=
            symmetric_protoindex_bundles.end());
        const auto spbundle_idx =
            symmetric_protoindex_bundles.find(ttpair.idx().proto_indices()) -
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
      assert(edges_.find(proto_index.full_label()) != edges_.end());
      const auto proto_index_vertex =
          edges_.find(proto_index.full_label()) - edges_.begin();
      graph->add_edge(spbundle_cnt, proto_index_vertex);
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
        ++tensor_cnt;
      });

  // compress vertex colors to 32 bits, as required by Bliss, by hashing
  size_t v_cnt = 0;
  for (auto &&color : vertex_color) {
    auto hash6432shift = [](size_t key) {
      static_assert(sizeof(key) == 8);
      key = (~key) + (key << 18);  // key = (key << 18) - key - 1;
      key = key ^ (key >> 31);
      key = key * 21;  // key = (key + (key << 2)) + (key << 4);
      key = key ^ (key >> 11);
      key = key + (key << 6);
      key = key ^ (key >> 22);
      return static_cast<int>(key);
    };
    graph->change_color(v_cnt, hash6432shift(color));
    ++v_cnt;
  }

  return {graph, vertex_labels, vertex_color, vertex_type};
}

void TensorNetwork::init_edges() const {
  if (have_edges_) return;

  auto idx_insert = [this](const Index &idx, int tensor_idx, int pos) {
    if (Logger::get_instance().tensor_network) {
      std::wcout << "TensorNetwork::init_edges: idx=" << to_latex(idx)
                 << " attached to tensor " << std::abs(tensor_idx) << "'s "
                 << ((tensor_idx > 0) ? "bra" : "ket") << " at position " << pos
                 << std::endl;
    }
    decltype(edges_) &indices = this->edges_;
    auto it = indices.find(idx.full_label());
    if (it == indices.end()) {
      indices.emplace(Edge(tensor_idx, &idx, pos));
    } else {
      const_cast<Edge &>(*it).connect_to(tensor_idx, pos);
    }
  };

  int t_idx = 1;
  for (auto &&t : tensors_) {
    const auto t_is_nonsymm = symmetry(*t) == Symmetry::nonsymm;
    size_t cnt = 0;
    for (const Index &idx : bra(*t)) {
      idx_insert(idx, t_idx, t_is_nonsymm ? cnt : 0);
      ++cnt;
    }
    cnt = 0;
    for (const Index &idx : ket(*t)) {
      idx_insert(idx, -t_idx, t_is_nonsymm ? cnt : 0);
      ++cnt;
    }
    ++t_idx;
  }

  // extract external indices
  for (const auto &terminals : edges_) {
    assert(terminals.size() != 0);
    if (terminals.size() == 1) {  // external?
      if (Logger::get_instance().tensor_network) {
        std::wcout << "idx " << to_latex(terminals.idx()) << " is external"
                   << std::endl;
      }
      auto insertion_result = ext_indices_.emplace(terminals.idx());
      assert(insertion_result.second);
    }
  }

  have_edges_ = true;
}

container::svector<std::pair<long, long>> TensorNetwork::factorize() {
  abort();  // not yet implemented
}

}  // namespace sequant
