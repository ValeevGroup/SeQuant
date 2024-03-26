//
// Created by Eduard Valeyev on 2019-02-26.
//

#include "wick_graph.hpp"
#include "bliss.hpp"
#include "logger.hpp"

namespace sequant {

std::tuple<std::shared_ptr<bliss::Graph>, std::vector<std::wstring>,
           std::vector<std::size_t>,
           std::vector<VertexType>>
WickGraph::make_bliss_graph(
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

void WickGraph::init_edges() const {
  if (have_edges_) return;

  auto idx_insert = [this](const Index &idx, int tensor_idx, int pos) {
    if (Logger::get_instance().tensor_network) {
      std::wcout << "WickGraph::init_edges: idx=" << to_latex(idx)
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

}  // namespace sequant
