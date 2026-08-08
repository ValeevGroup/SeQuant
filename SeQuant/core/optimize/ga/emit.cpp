#include <SeQuant/core/optimize/ga/emit.hpp>

namespace sequant::opt::ga {

namespace {

using IndexMap = container::map<Index, Index>;

void collect_indices(ExprPtr const& e, container::svector<Index>& out) {
  if (e->is<Tensor>()) {
    for (auto const& ix : e->as<Tensor>().const_indices()) {
      for (auto const& p : ix.proto_indices())
        if (ranges::find(out, p) == ranges::end(out)) out.push_back(p);
      if (ranges::find(out, ix) == ranges::end(out)) out.push_back(ix);
    }
  } else if (e->is<Product>()) {
    for (auto const& f : e->as<Product>().factors()) collect_indices(f, out);
  } else if (e->is<Sum>()) {
    for (auto const& s : e->as<Sum>().summands()) collect_indices(s, out);
  }
}

ExprPtr apply_map(ExprPtr const& e, IndexMap const& m) {
  if (e->is<Tensor>()) {
    ExprPtr t = e->as<Tensor>().clone();
    t->as<Tensor>().transform_indices(m);
    // Index::transform tags replaced indices to guard against double
    // substitution; clear them so downstream canonicalization sees a clean
    // slate
    t->as<Tensor>().reset_tags();
    return t;
  }
  if (e->is<Product>()) {
    auto const& p = e->as<Product>();
    Product out(p.scalar(), ExprPtrList{});
    for (auto const& f : p.factors())
      out.append(1, apply_map(f, m), Product::Flatten::No);
    return ex<Product>(std::move(out));
  }
  if (e->is<Sum>()) {
    Sum out;
    for (auto const& s : e->as<Sum>().summands()) out.append(apply_map(s, m));
    return ex<Sum>(std::move(out));
  }
  return e;
}

// rename e: the seed map m wires e's face indices into the caller's
// names; every other index becomes a fresh tmp index (non-composites
// first, then composites, whose protos may have moved)
template <typename FacePred>
ExprPtr rewire(ExprPtr const& e, IndexMap m, FacePred const& is_face) {
  container::svector<Index> all;
  collect_indices(e, all);
  for (auto const& ix : all)
    if (!ix.has_proto_indices() && !m.contains(ix) && !is_face(ix))
      m.emplace(ix, Index::make_tmp_index(ix.space()));
  for (auto const& ix : all) {
    if (!ix.has_proto_indices() || m.contains(ix) || is_face(ix)) continue;
    Index::index_vector protos;
    for (auto const& p : ix.proto_indices()) {
      auto it = m.find(p);
      protos.push_back(it == m.end() ? p : it->second);
    }
    m.emplace(ix, Index::make_tmp_index(ix.space(), std::move(protos)));
  }
  if (m.empty()) return e;
  return apply_map(e, m);
}

/// The named leaf tensor of intermediate \p label with slot list \p slots
/// (all slots auxiliary, no declared symmetry, so nothing downstream may
/// legally reorder them: slot position IS the array axis).
ExprPtr named_leaf(std::wstring const& label,
                   container::svector<Index> const& slots) {
  return ex<Tensor>(label, bra(), ket(),
                    aux(container::svector<Index>(slots.begin(), slots.end())),
                    Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                    ColumnSymmetry::Nonsymm);
}

struct Emitter {
  Emitter(Fitness const& f, Schedule const& s) : F(f), sch(s) {}

  Fitness const& F;
  Schedule const& sch;

  /// key -> number of render sites (the producer's own site included)
  container::map<std::size_t, std::size_t> use_count;
  /// key -> assigned intermediate label (shared keys only)
  container::map<std::size_t, std::wstring> names;
  /// definitions in emission (dependency) order
  container::svector<ResultExpr> defs;
  std::size_t name_counter = 1;

  KeyTable const& kt() const { return F.table(); }

  // --- pass 1: count how many sites each keyed array is rendered at -------

  void count_cl(int d, NodeMask S) {
    if (std::popcount(S) == 1) return;
    auto const k = kt().terms[d].key[S];
    auto it = sch.pick.find(k);
    Cluster const p = it != sch.pick.end() ? it->second : Cluster{d, S};
    if (it != sch.pick.end() && ++use_count[k] > 1)
      return;  // the producer's body is walked once
    auto const* cc = sch.forest.terms[p.d].children_of(p.S);
    SEQUANT_ASSERT(cc);
    count_cl(p.d, cc[0]);
    count_cl(p.d, cc[1]);
  }

  void count(Summand const& s) {
    Val const& v = *s.val;
    switch (v.kind) {
      case Val::Cl:
        count_cl(v.d, v.S);
        return;
      case Val::Fx:
        count(Summand{1, v.inner, {}});
        count_cl(v.d, v.V);
        return;
      case Val::Sum:
        count(v.s1);
        count(v.s2);
        return;
    }
    SEQUANT_UNREACHABLE_TOKEN;
  }

  // --- pass 2: render, defining shared keys once and referencing by name --

  bool shared(std::size_t key) const {
    auto it = use_count.find(key);
    return it != use_count.end() && it->second >= 2;
  }

  /// Ensure the definition of shared key \p k (producer \p p) has been
  /// emitted; definitions this one depends on are emitted first (the body
  /// render below references them by name).
  /// \note returns (and captures) the label BY VALUE: `names` is a flat_map,
  /// so the nested inserts of the recursive body render invalidate
  /// references into it.
  std::wstring define(std::size_t k, Cluster p) {
    if (auto it = names.find(k); it != names.end()) return it->second;
    std::wstring const label = std::wstring(named_intermediate_prefix) +
                               std::to_wstring(name_counter++);
    names.emplace(k, label);
    auto const* cc = sch.forest.terms[p.d].children_of(p.S);
    SEQUANT_ASSERT(cc);
    ExprPtr body = ex<Product>(
        Product{1, ExprPtrList{render_cl(p.d, cc[0]), render_cl(p.d, cc[1])},
                Product::Flatten::No});
    auto const slots = kt().terms[p.d].canon_face_indices(p.S);
    SEQUANT_ASSERT(!slots.empty() &&
                   "scalar-valued GA intermediates are not supported");
    defs.emplace_back(named_leaf(label, slots)->as<Tensor>(), std::move(body));
    return label;
  }

  ExprPtr render_cl(int d, NodeMask S) {
    if (std::popcount(S) == 1)
      return kt().terms[d].tensors[std::countr_zero(S)];
    auto const k = kt().terms[d].key[S];
    auto it = sch.pick.find(k);
    if (it != sch.pick.end() && shared(k)) {
      // shared array: define once under a canonical-face head, reference by
      // name here; this site's slots are ITS canonical face order, which
      // zips axis-wise with the definition head's (the same canonical
      // correspondence render_cl's producer substitution used to rewire by)
      std::wstring const label = define(k, it->second);
      return named_leaf(label, kt().terms[d].canon_face_indices(S));
    }
    if (it != sch.pick.end() && it->second != Cluster{d, S}) {
      // single-use array built by an isomorphic cluster elsewhere: inline
      // the producer's subtree with the face renamed through the canonical
      // correspondence (fresh internal labels)
      Cluster const p = it->second;
      auto const& sigmas = F.correspondences(p, Cluster{d, S});
      SEQUANT_ASSERT(!sigmas.empty());
      IndexMap m;
      for (auto const& [from, to] : sigmas.front())
        if (!(from == to)) m.emplace(from, to);
      // F(p) is stored as bits since T-A7; materialize the set here (emission
      // only, once per inlined producer) so the membership test is the same
      // ordered-set lookup it always was
      auto const pface = kt().terms[p.d].face_set(p.S);
      return rewire(render_cl(p.d, p.S), std::move(m), [&](Index const& ix) {
        return pface.find(ix) != pface.end();
      });
    }
    auto const* cc = sch.forest.terms[d].children_of(S);
    SEQUANT_ASSERT(cc);
    return ex<Product>(Product{1,
                               ExprPtrList{render_cl(d, cc[0]),
                                           render_cl(d, cc[1])},
                               Product::Flatten::No});
  }

  // rename e (a rendered value of face `face_of_e` with ambient `amb`) so
  // that its face indices match `target_face` under the tag identification
  // `target_amb`, and all its internal indices are fresh.
  //
  // Since T-A6 an `AmbientMap` is keyed by index BIT of its own term, so this
  // is where the real `Index` objects come back: `index_list[bit]` of
  // `amb.d` / `target_amb.d`. Tags are unique within a map, so the `find_if`
  // below is the same lookup it always was.
  ExprPtr align(ExprPtr const& e, AmbientMap const& amb,
                AmbientMap const& target_amb) const {
    auto const& src = kt().terms[amb.d].index_list;
    auto const& dst = kt().terms[target_amb.d].index_list;
    container::svector<Index> amb_face;
    amb_face.reserve(amb.e.size());
    for (auto const& [b, tag] : amb.e) amb_face.push_back(src[b]);
    IndexMap m;
    for (std::size_t i = 0; i < amb.e.size(); ++i) {
      const AmbientTag tag = amb.e[i].second;
      auto it =
          std::find_if(target_amb.e.begin(), target_amb.e.end(),
                       [tag](auto const& kv) { return kv.second == tag; });
      SEQUANT_ASSERT(it != target_amb.e.end());
      Index const& x2 = amb_face[i];
      Index const& y = dst[it->first];
      if (!(x2 == y)) m.emplace(x2, y);
    }
    return rewire(e, std::move(m), [&](Index const& ix) {
      return ranges::find(amb_face, ix) != ranges::end(amb_face);
    });
  }

  static ExprPtr scaled(Constant::scalar_type c, ExprPtr e) {
    return c == Constant::scalar_type(1) ? e : ex<Constant>(c) * e;
  }

  ExprPtr render(Summand const& s) {
    Val const& v = *s.val;
    switch (v.kind) {
      case Val::Cl:
        return render_cl(v.d, v.S);
      case Val::Fx: {
        // (combined residual) * A_{d,V}; the residual's face is F(d, X) in
        // d's own names, so the wiring is that of d's X*V merge
        ExprPtr inner = render(Summand{1, v.inner, {}});
        return ex<Product>(Product{
            1, ExprPtrList{inner, render_cl(v.d, v.V)}, Product::Flatten::No});
      }
      case Val::Sum: {
        ExprPtr e1 = scaled(v.s1.coeff, render(v.s1));
        ExprPtr e2 = scaled(v.s2.coeff, render(v.s2));
        return e1 + align(e2, v.s2.ambient, v.s1.ambient);
      }
    }
    SEQUANT_UNREACHABLE_TOKEN;
  }

  ExprPtr target(std::size_t t) {
    Summand const& root = sch.roots[t];
    ExprPtr e = scaled(root.coeff, render(root));
    // externals: rename to the target's first term's names. That eta map is
    // exactly the leaf ambient of that term, which the Fitness already holds
    // (it used to be rebuilt here, byte for byte, from `ref.ext`).
    const int ref = static_cast<int>(kt().targets[t].terms.front());
    return align(e, root.ambient, F.leaf_ambient(ref));
  }
};

/// One definition ready for substitution: its slot list and its fully
/// inlined body.
struct InlinedDef {
  container::svector<Index> slots;
  ExprPtr body;
};

ExprPtr substitute(ExprPtr const& e,
                   container::map<std::wstring, InlinedDef> const& defmap) {
  if (e->is<Tensor>()) {
    auto const& t = e->as<Tensor>();
    auto it = defmap.find(std::wstring(t.label()));
    if (it == defmap.end()) return e;
    auto const& def = it->second;
    auto const& use_slots = t.aux();
    SEQUANT_ASSERT(use_slots.size() == def.slots.size());
    // positional face map (definition names -> this use's names). Slot pairs
    // first, IDENTITY PAIRS INCLUDED -- they are the authoritative
    // correspondence and must claim their def index so the proto zip below
    // can never rebind it. Protos of composite slots are zipped in only as
    // gap fillers: normally every proto is a slot in its own right, and a
    // composite's STORED proto order may be normalized (sorted), so its
    // positional zip is not trustworthy where a slot pair exists.
    IndexMap m;
    TensorNetwork::NamedIndexSet dface;
    for (std::size_t k = 0; k < def.slots.size(); ++k) {
      auto const& ds = def.slots[k];
      auto const& us = use_slots[k];
      dface.emplace(ds);
      m.try_emplace(ds, us);
    }
    for (std::size_t k = 0; k < def.slots.size(); ++k) {
      auto const& ds = def.slots[k];
      auto const& us = use_slots[k];
      SEQUANT_ASSERT(ds.proto_indices().size() == us.proto_indices().size());
      for (std::size_t j = 0; j < ds.proto_indices().size(); ++j) {
        auto const& dp = ds.proto_indices()[j];
        auto const& up = us.proto_indices()[j];
        dface.emplace(dp);
        m.try_emplace(dp, up);
      }
    }
    return rewire(def.body, std::move(m), [&](Index const& ix) {
      return dface.find(ix) != dface.end();
    });
  }
  if (e->is<Product>()) {
    auto const& p = e->as<Product>();
    Product out(p.scalar(), ExprPtrList{});
    for (auto const& f : p.factors())
      out.append(1, substitute(f, defmap), Product::Flatten::No);
    return ex<Product>(std::move(out));
  }
  if (e->is<Sum>()) {
    Sum out;
    for (auto const& s : e->as<Sum>().summands())
      out.append(substitute(s, defmap));
    return ex<Sum>(std::move(out));
  }
  return e;
}

}  // namespace

Emission emit_named(Fitness const& fitness, Schedule const& schedule) {
  Emitter em{fitness, schedule};
  for (auto const& root : schedule.roots) em.count(root);
  Emission out;
  for (std::size_t t = 0; t < fitness.table().targets.size(); ++t)
    out.targets.push_back(em.target(t));
  out.definitions = std::move(em.defs);
  return out;
}

container::svector<ExprPtr> inline_definitions(Emission const& em) {
  container::map<std::wstring, InlinedDef> defmap;
  for (auto const& def : em.definitions) {
    // definitions are in dependency order: inline earlier ones into this body
    ExprPtr body = substitute(def.expression(), defmap);
    defmap.emplace(def.label(),
                   InlinedDef{container::svector<Index>(def.aux().begin(),
                                                        def.aux().end()),
                              std::move(body)});
  }
  container::svector<ExprPtr> out;
  for (auto const& t : em.targets) out.push_back(substitute(t, defmap));
  return out;
}

container::svector<ExprPtr> emit(Fitness const& fitness,
                                 Schedule const& schedule) {
  return inline_definitions(emit_named(fitness, schedule));
}

}  // namespace sequant::opt::ga
