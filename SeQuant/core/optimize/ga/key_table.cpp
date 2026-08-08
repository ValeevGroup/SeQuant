#include <SeQuant/core/optimize/ga/key_table.hpp>

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace sequant::opt::ga {

/// key_table.hpp spells `FaceSet` out rather than aliasing
/// `TensorNetwork::NamedIndexSet`, to stay free of tensor_network.hpp; this is
/// what keeps that spelling honest.
static_assert(std::is_same_v<FaceSet, TensorNetwork::NamedIndexSet>);

int TermTable::bit_of(Index const& ix) const {
  auto it = std::find(index_list.begin(), index_list.end(), ix);
  SEQUANT_ASSERT(it != index_list.end());
  return static_cast<int>(it - index_list.begin());
}

double TermTable::face_size(NodeMask S) const {
  IndexMask m = face_mask(S);
  double v = 1.;
  for (; m; m &= m - 1) v *= extent[std::countr_zero(m)];
  return v;
}

FaceSet TermTable::face_set(NodeMask S) const {
  const IndexMask fm = face_mask(S);
  FaceSet f;
  f.reserve(static_cast<std::size_t>(std::popcount(fm)));
  for (IndexMask m = fm; m; m &= m - 1)
    f.emplace(index_list[std::countr_zero(m)]);
  return f;
}

container::svector<Index> TermTable::canon_face_indices(NodeMask S) const {
  auto const& bits = canon_face_bits[S];
  container::svector<Index> out;
  out.reserve(bits.size());
  for (std::uint8_t b : bits) out.push_back(index_list[b]);
  return out;
}

int KeyTable::kind_of(Index const& ix) const {
  auto it = kind_ids.find({ix.space(), ix.proto_indices().size()});
  SEQUANT_ASSERT(it != kind_ids.end());
  return it->second;
}

namespace {

// Index -> bit position, for the setup-only interning of a canonical face.
// Throws rather than returning a sentinel: asserts are compiled out of the
// release build the optimizer runs in, and a bogus slot id would be a silent
// axis permutation of every shared intermediate.
std::uint8_t intern(TermTable const& T, NodeMask S, Index const& ix) {
  for (IndexMask m = T.face_mask(S); m; m &= m - 1) {
    const int b = std::countr_zero(m);
    if (T.index_list[b] == ix) return static_cast<std::uint8_t>(b);
  }
  auto it = std::find(T.index_list.begin(), T.index_list.end(), ix);
  if (it == T.index_list.end())
    throw std::logic_error(
        "ga::build_key_table: a canonical face slot is not an index of its "
        "term");
  return static_cast<std::uint8_t>(it - T.index_list.begin());
}

// The canonicalize sweep runs in two passes, and WHICH pass is parallel is
// load-bearing: `meta_to_id` hands out key ids in FIRST-ENCOUNTER order, and
// those ids order `Fitness::resolve`'s ascending-key walk, hence its
// floating-point summation order, hence possibly which schedule wins. Any
// scheme that lets threads insert into `meta_to_id` renumbers the table.
// So pass 1 (parallel) only canonicalizes subsets into its own pre-sized slots
// -- touching no shared mutable state, on per-task clones, since Expr and
// Index carry unsynchronized mutable caches -- and pass 2 (serial) walks S
// ascending assigning ids one `try_emplace` at a time. That is a
// permutation-free replay of the serial insertion sequence, so every
// `T.key[S]` is bit-identical at any thread count.
//
// The passes interleave in chunks so a whole term's metas are never held at
// once; pass 2 still walks S ascending ACROSS chunk boundaries, so chunking
// cannot move an id. The chunk size is measured, not guessed.
constexpr std::size_t canon_chunk_subsets = 256;
// below this many subsets in a chunk the thread spawn costs more than the work
constexpr std::size_t canon_min_parallel_subsets = 32;

// `opt::detail::SubNetIdMap` with the graph hash precomputed in pass 1, where
// the graph is a task's private property (`get_hash64` is neither memoized nor
// cheap, and in the hasher it would run serially inside pass 2). Key identity
// is untouched: equality is `ConstGraphCmp::cmp == 0` and ids are still
// `map.size() - 1` on insertion in ascending-S order. The hash's consistency
// with the equality is load-bearing: `cmp` compares vertex count, colors,
// degrees and edge lists, and `get_hash64` hashes exactly that data.
struct HashedMeta {
  TensorNetwork::SlotCanonicalizationMetadata meta;
  /// meta.graph->get_hash64(); also normalizes the graph, which `cmp`'s const
  /// overload relies on having happened (it neither dedups nor sorts edges).
  std::size_t hash = 0;
};
struct HashedMetaHash {
  std::size_t operator()(HashedMeta const& m) const noexcept { return m.hash; }
};
struct HashedMetaEqual {
  bool operator()(HashedMeta const& a, HashedMeta const& b) const {
    return bliss::ConstGraphCmp::cmp(*a.meta.graph, *b.meta.graph) == 0;
  }
};
using MetaIdMap =
    container::unordered_map<HashedMeta, std::size_t, HashedMetaHash,
                             HashedMetaEqual>;

int register_kind(KeyTable& kt, Index const& ix) {
  auto [it, inserted] = kt.kind_ids.try_emplace(
      {ix.space(), ix.proto_indices().size()},
      static_cast<int>(kt.kind_ids.size()));
  return it->second;
}

// One term's tables. Faces follow the prototype exactly (Definition 2.4):
// eff(v) = slots of tensor v plus the proto-indices of its composite slots;
// F(S) = indices of union_{v in S} eff(v) that are external or carried by a
// tensor outside S.
TermTable build_term(KeyTable& kt, ExprPtr const& summand, FaceSet const& ext,
                     std::function<std::size_t(Index const&)> const& ixex,
                     std::function<bool(Tensor const&)> const& is_volatile_leaf,
                     MetaIdMap& meta_to_id) {
  TermTable T;
  T.ext = ext;
  if (summand->is<Tensor>()) {
    T.tensors.push_back(summand);
  } else {
    auto const& prod = summand->as<Product>();
    T.scalar = prod.scalar();
    for (auto const& f : prod.factors()) {
      SEQUANT_ASSERT(f->is<Tensor>());
      T.tensors.push_back(f);
    }
  }
  const std::size_t n = T.n();
  if (n >= 8 * sizeof(NodeMask))
    throw std::invalid_argument("ga::build_key_table: a term has " +
                                std::to_string(n) +
                                " factors; NodeMask holds < " +
                                std::to_string(8 * sizeof(NodeMask)));

  // Amplitude-dependent leaves of this term. Recorded as a mask so cluster
  // volatility is one AND (see TermTable::volatile_mask).
  if (is_volatile_leaf)
    for (std::size_t v = 0; v < n; ++v)
      if (is_volatile_leaf(T.tensors[v]->as<Tensor>()))
        T.volatile_mask |= NodeMask{1} << v;

  // per-node effective index sets, and the term's index universe
  auto note_index = [&](Index const& ix) {
    if (ranges::find(T.index_list, ix) == ranges::end(T.index_list)) {
      T.index_list.push_back(ix);
      T.extent.push_back(static_cast<double>(ixex(ix)));
      T.kind.push_back(register_kind(kt, ix));
    }
  };
  container::svector<IndexMask> eff(n, 0);
  for (std::size_t v = 0; v < n; ++v) {
    auto const& t = T.tensors[v]->as<Tensor>();
    for (auto const& ix : t.const_indices()) {
      for (auto const& p : ix.proto_indices()) note_index(p);
      note_index(ix);
      eff[v] |= IndexMask{1} << T.bit_of(ix);
      for (auto const& p : ix.proto_indices())
        eff[v] |= IndexMask{1} << T.bit_of(p);
    }
  }
  SEQUANT_ASSERT(T.index_list.size() <= 8 * sizeof(IndexMask));

  // `index_rank` MUST come from an actual `std::sort` over the real `Index`
  // objects -- never from bit order, label order or anything else that merely
  // tends to agree with it -- because the L2 pass's bijection enumeration
  // order is exactly that sort order. See TermTable::index_rank.
  {
    const std::size_t ni = T.index_list.size();
    container::svector<std::uint8_t> ord(ni);
    for (std::size_t b = 0; b < ni; ++b) ord[b] = static_cast<std::uint8_t>(b);
    std::sort(ord.begin(), ord.end(), [&](std::uint8_t a, std::uint8_t b) {
      return T.index_list[a] < T.index_list[b];
    });
    T.index_rank.assign(ni, 0);
    for (std::size_t r = 0; r < ni; ++r)
      T.index_rank[ord[r]] = static_cast<std::uint8_t>(r);
    // reading `index_list` through the inverse permutation is ascending
    // `Index` order
    for (std::size_t r = 1; r < ni; ++r)
      SEQUANT_ASSERT(T.index_list[ord[r - 1]] < T.index_list[ord[r]] &&
                     T.index_rank[ord[r]] == r);

    T.proto_mask.assign(ni, 0);
    int kmax = -1;
    for (std::size_t b = 0; b < ni; ++b) {
      for (auto const& p : T.index_list[b].proto_indices())
        T.proto_mask[b] |= IndexMask{1} << T.bit_of(p);
      kmax = std::max(kmax, T.kind[b]);
    }
    // padded to the table's final kind count by build_key_table below
    T.kind_mask.assign(static_cast<std::size_t>(kmax + 1), 0);
    for (std::size_t b = 0; b < ni; ++b)
      T.kind_mask[static_cast<std::size_t>(T.kind[b])] |= IndexMask{1} << b;
  }

  IndexMask ext_bits = 0;
  for (auto const& ix : ext) ext_bits |= IndexMask{1} << T.bit_of(ix);
  container::svector<NodeMask> carr(T.index_list.size(), 0);
  for (std::size_t v = 0; v < n; ++v)
    for (IndexMask m = eff[v]; m; m &= m - 1)
      carr[std::countr_zero(m)] |= NodeMask{1} << v;

  // F(S), stored ONLY as the (disjoint) composite / non-composite bit pair;
  // `face_mask` recombines them, `face_set` materializes the `Index` form.
  const std::size_t n_subsets = std::size_t{1} << n;
  T.outer_mask.resize(n_subsets, 0);
  T.inner_mask.resize(n_subsets, 0);
  for (std::size_t S = 1; S < n_subsets; ++S) {
    IndexMask effS = 0;
    for (NodeMask m = static_cast<NodeMask>(S); m; m &= m - 1)
      effS |= eff[std::countr_zero(m)];
    for (IndexMask m = effS; m; m &= m - 1) {
      const int b = std::countr_zero(m);
      const IndexMask bit = IndexMask{1} << b;
      if (!(ext_bits & bit) && (carr[b] & ~static_cast<NodeMask>(S)) == 0)
        continue;  // fully absorbed inside S and not external
      (T.index_list[b].has_proto_indices() ? T.inner_mask
                                           : T.outer_mask)[S] |= bit;
    }
  }

  // Global canonical-subnet key for every nonempty subset, built from clones
  // (canonicalize_slots mutates lazily-cached state). Named indices are colored
  // by kind only, so the key is label-blind array identity. See
  // canon_chunk_subsets for why the id assignment may not be parallelized.
  T.key.resize(n_subsets, no_key);
  T.canon_face_bits.resize(n_subsets);

  // read-only inside the parallel region (configuring it is the caller's
  // business -- same contract as optimize.cpp's parallel branch)
  auto const& cardinal = TensorCanonicalizer::cardinal_tensor_labels();
  const std::size_t chunk =
      std::min<std::size_t>(n_subsets - 1, canon_chunk_subsets);
  container::vector<HashedMeta> metas(chunk);
  const bool parallel = num_threads() > 1;

  for (std::size_t lo = 1; lo < n_subsets; lo += chunk) {
    const std::size_t hi = std::min(n_subsets, lo + chunk);

    // ---- pass 1: canonicalize S in [lo, hi), each task in its own slot ----
    auto canonicalize_one = [&, lo](std::size_t S) {
      container::vector<ExprPtr> ts;
      for (std::size_t v = 0; v < n; ++v)
        if (S & (std::size_t{1} << v))
          ts.emplace_back(T.tensors[v]->as<Tensor>().clone());
      auto tn = TensorNetwork{ts};
      // `fs` must outlive the read-back below: canonicalize_slots copies it
      // and `named_indices_canonical` points into that copy.
      const FaceSet fs = T.face_set(static_cast<NodeMask>(S));
      auto meta = tn.canonicalize_slots(cardinal, &fs);
      auto& cfb = T.canon_face_bits[S];
      cfb.reserve(meta.named_indices_canonical.size());
      // push_back in the metadata's own order -- this loop, and nothing else,
      // is what makes canon_face_bits an axis ORDER rather than an index set
      for (auto const& it : meta.named_indices_canonical)
        cfb.push_back(intern(T, static_cast<NodeMask>(S), *it));
      auto& slot = metas[S - lo];
      slot.hash = meta.graph->get_hash64();  // and normalizes meta.graph
      slot.meta = std::move(meta);
    };
    if (parallel && hi - lo >= canon_min_parallel_subsets) {
      auto subsets = ranges::views::iota(lo, hi);
      sequant::for_each(subsets, canonicalize_one);
    } else {
      for (std::size_t S = lo; S < hi; ++S) canonicalize_one(S);
    }

    // ---- pass 2: ascending S, one id at a time, exactly the serial order ----
    for (std::size_t S = lo; S < hi; ++S) {
      auto [it, inserted] = meta_to_id.try_emplace(std::move(metas[S - lo]), 0);
      if (inserted) it->second = meta_to_id.size() - 1;
      T.key[S] = it->second;
    }
  }
  return T;
}

}  // namespace

KeyTable build_key_table(
    container::svector<TargetInput> const& targets,
    std::function<std::size_t(Index const&)> const& idx_to_extent,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  KeyTable kt;
  kt.volatility_aware = static_cast<bool>(is_volatile_leaf);
  MetaIdMap meta_to_id;
  for (auto const& tgt : targets) {
    SEQUANT_ASSERT(tgt.summands.size() == tgt.ext.size());
    TargetBlock blk;
    // the L2 tree over this target's summands uses the same mask codec as
    // the per-term trees (fitness.cpp builds its root as (1 << n) - 1)
    if (tgt.summands.empty())
      throw std::invalid_argument(
          "ga::build_key_table: target has no summands; the L2 root mask would "
          "be 0, and the fitness/explain walks over it do not terminate");
    if (tgt.summands.size() >= 8 * sizeof(NodeMask))
      throw std::invalid_argument(
          "ga::build_key_table: target has " +
          std::to_string(tgt.summands.size()) + " summands; NodeMask holds < " +
          std::to_string(8 * sizeof(NodeMask)));
    for (std::size_t s = 0; s < tgt.summands.size(); ++s) {
      blk.terms.push_back(kt.terms.size());
      kt.terms.push_back(build_term(kt, tgt.summands[s], tgt.ext[s],
                                    idx_to_extent, is_volatile_leaf,
                                    meta_to_id));
    }
    kt.targets.push_back(std::move(blk));
  }
  kt.n_keys = meta_to_id.size();
  // Kind ids are handed out across the whole table, so pad every per-term kind
  // mask now: `T.kind_mask[k]` must be addressable for any k find_beta asks.
  for (auto& T : kt.terms) T.kind_mask.resize(kt.kind_ids.size(), 0);
  return kt;
}

}  // namespace sequant::opt::ga
