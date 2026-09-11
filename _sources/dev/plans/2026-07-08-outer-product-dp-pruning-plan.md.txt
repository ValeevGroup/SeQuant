# Outer-product (disconnected-subgraph) DP pruning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prune the single-term contraction DP so it never spends work on outer
products (bipartitions whose two parts share no contractible index), collapsing
the per-term search to the number of connected subnetworks, with the optimal
factorization unchanged.

**Architecture:** Add two pure helpers -- `contractible_adjacency` (pairwise
"share a non-target index" graph as per-tensor neighbor bitmasks) and
`connected_subsets` (per-subset induced-connectivity) -- plus a thin
`outer_product_connectivity` wrapper that returns the connectivity mask (or an
all-connected mask when pruning is disabled or the full network is a genuine
product). `solve_single_term` computes this mask from the `network`/`tidxs` it
already receives and skips disconnected subsets / relaxes only connected-child
bipartitions; `build_subnet_metadata` takes the same mask and skips
canonicalizing disconnected subsets. No `build_context`/CostModel-concept
signature changes.

**Tech Stack:** SeQuant C++20, `SeQuant/core/optimize/`, Catch2 unit tests.

## Global Constraints

- Repo: SeQuant, branch `feature/eval-predicted-peak-trace`. Branch-only: no
  push, no PR, no remote operations, do not switch branches.
- Commit messages: plain, describing the change. **No `Co-Authored-By` trailers
  or any AI/tool attribution.**
- **Do NOT stage, commit, or modify** `SeQuant/core/tensor_network/slot.hpp`
  (unrelated pre-existing local edit) or `tests/unit/test_eval_dryrun.cpp`
  (an uncommitted local spike). Commit only the files each task names.
- ASCII only in all source and docs. No U+2013 (en-dash), no U+00A0
  (non-breaking space) -- a pre-commit hook rejects them.
- Run `clang-format` before every commit:
  `/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i <files>` (skip
  `.tpp` files; none here).
- Build and test in the existing configured tree `SeQuant/build-test`:
  - build: `cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8`
  - run a tag: `./tests/unit/unit_tests-sequant "[tag]"` from `build-test`.
- **Never run the unfiltered unit suite** (it hangs on an unrelated slow DP
  test). Always pass a Catch2 tag/name filter. **Never run two SeQuant test
  binaries concurrently.**
- Adjacency compares **top-level** (bra/ket/aux) indices only, never
  protoindices. "Contractible" index := an index **not** in the target set
  `tidxs`.
- Loss-free requirement: for every objective, the pruned DP must return the
  byte-identical optimized expression as the unpruned DP (the parity gate,
  Task 4). Existing `[optimize]` tests must continue to pass unchanged.

## Test scaffolding (reuse in every new test case)

Mirror the existing `"optimize"` `TEST_CASE` (`test_optimize.cpp:168-201`)
exactly -- do not invent a new context or parser:

- **Context:** inside each `TEST_CASE`, do
  `auto ctx_resetter = set_scoped_default_context(get_default_context().clone());`
  then `auto reg = get_default_context().mutable_index_space_registry();` and set
  `reg->retrieve_ptr(L"i")->approximate_size(10);` (occ),
  `reg->retrieve_ptr(L"a")->approximate_size(100);` (virt),
  `reg->retrieve_ptr(L"x")->approximate_size(4);` (aux).
- **Spaces:** occupied = `i`, virtual = `a`, aux = `x`.
- **Parse:** `auto parse = [](auto const& s){ return deserialize(s, {.def_perm_symm = Symmetry::Antisymm}); };`
  Shorthand is LaTeX-like: `g_{i3,i4}^{a3,a4}` = bra `i3,i4`, ket `a3,a4`. A single
  tensor string (e.g. `L"g_{i1}^{a1}"`) parses to a tensor expression; build a
  `TensorNetwork` from the `Tensor` factors of a parsed product, mirroring
  `single_term_detail.hpp:447-452` (`for (auto&& f : prod.factors()) if (f->is<Tensor>()) v.push_back(f);`
  then `TensorNetwork{v}`).
- **Include** `<cstdlib>` in the test TU for `setenv`/`unsetenv`.
- **Tensor ordering:** `TensorNetwork` may reorder its tensors, so assert
  **order-independent** properties (edge count, whole-set connectivity, all-zero)
  rather than a specific `adj[i]` value tied to input position.

---

## File Structure

- `SeQuant/core/optimize/single_term_detail.hpp` -- add `contractible_adjacency`,
  `connected_subsets`, `outer_product_connectivity` (namespace
  `sequant::opt::detail`); add the `connected` parameter + skip to
  `build_subnet_metadata`.
- `SeQuant/core/optimize/cost_model.hpp` -- `solve_single_term` computes the mask
  and applies the skip/guard; `AdditiveModel::build_context` computes the mask
  and passes it to `build_subnet_metadata`.
- `tests/unit/test_optimize.cpp` -- new `TEST_CASE`s for the helpers, the parity
  gate, the performance smoke, and the multi-component safety net.

---

### Task 1: `contractible_adjacency` helper

**Files:**
- Modify: `SeQuant/core/optimize/single_term_detail.hpp` (add helper near
  `init_results`, ~line 414; add `#include <cstdlib>` and
  `#include <range/v3/algorithm/contains.hpp>` with the other includes if not
  present)
- Test: `tests/unit/test_optimize.cpp`

**Interfaces:**
- Produces:
  `template <typename TIdxs> container::vector<std::size_t> contractible_adjacency(TensorNetwork const& network, TIdxs const& tidxs)`
  in `sequant::opt::detail`. Returns `adj` of size `network.tensors().size()`;
  `adj[i]` has bit `j` set iff tensors `i` and `j` share a top-level index that
  is not in `tidxs`. Never sets a tensor's own bit.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_optimize.cpp`. Use the Test-scaffolding context/parse;
assert order-independent facts (undirected edge count = half the summed
popcount; whole-set connectivity via `connected_subsets`).

```cpp
TEST_CASE("contractible_adjacency", "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto parse = [](auto const& s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Antisymm});
  };
  auto net_of = [](ExprPtr const& prod) {
    container::vector<ExprPtr> v;
    for (auto&& f : prod->as<Product>().factors())
      if (f->is<Tensor>()) v.push_back(f);
    return TensorNetwork{v};
  };
  auto edge_count = [](container::vector<std::size_t> const& adj) {
    std::size_t s = 0;
    for (auto m : adj) s += static_cast<std::size_t>(std::popcount(m));
    return s / 2;  // each undirected edge counted twice
  };

  // Chain: f-g share a1; g-h share a2. Targets {i1,i2}. 2 edges.
  auto chain = net_of(parse(L"f_{i1}^{a1} g_{a1}^{a2} h_{a2}^{i2}"));
  std::vector<Index> tgt2{Index{L"i_1"}, Index{L"i_2"}};
  auto adj_chain = o::contractible_adjacency(chain, tgt2);
  REQUIRE(adj_chain.size() == 3);
  CHECK(edge_count(adj_chain) == 2);

  // Hyperedge: p,q,r all carry summed a5 -> clique (3 edges).
  auto hyper = net_of(parse(L"p_{i1}^{a5} q_{i2}^{a5} r_{i3}^{a5}"));
  std::vector<Index> tgt3{Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"}};
  auto adj_hyper = o::contractible_adjacency(hyper, tgt3);
  CHECK(edge_count(adj_hyper) == 3);

  // Spectator-only: two tensors share only the target index i1 -> no edges.
  auto spec = net_of(parse(L"u_{i1}^{a1} v_{i1}^{a2}"));
  std::vector<Index> tgt_spec{Index{L"i_1"}, Index{L"a_1"}, Index{L"a_2"}};
  auto adj_spec = o::contractible_adjacency(spec, tgt_spec);
  CHECK(edge_count(adj_spec) == 0);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `SeQuant/build-test`, after adding the case): build fails to compile
(`contractible_adjacency` not declared). That is the expected initial failure.

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8
```
Expected: compile error, `no member named 'contractible_adjacency'`.

- [ ] **Step 3: Write minimal implementation**

In `single_term_detail.hpp`, in `namespace sequant::opt { namespace detail {`,
after `init_results` (~line 414), add:

```cpp
/// \brief Per-tensor adjacency bitmask over "contractible" shared indices.
///
/// `adj[i]` has bit `j` set iff tensors `i` and `j` share at least one top-level
/// (bra/ket/aux) index that is NOT a target index (`tidxs`) -- i.e. an index
/// that is summed somewhere in the term. Protoindices are never compared, so two
/// composites `a<i,j>`, `b<i,j>` that share only occupied protos create no edge.
/// A hyperedge (a contractible index on three-plus tensors) makes all its
/// carriers mutually adjacent. A tensor is never adjacent to itself.
template <typename TIdxs>
inline container::vector<std::size_t> contractible_adjacency(
    TensorNetwork const& network, TIdxs const& tidxs) {
  std::size_t const nt = network.tensors().size();
  container::vector<std::size_t> adj(nt, 0);
  // carrier bitmask per contractible (non-target) top-level index
  container::map<Index, std::size_t, Index::FullLabelCompare> carriers;
  std::size_t i = 0;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    for (auto&& ix : ranges::views::concat(tp->bra(), tp->ket(), tp->aux())) {
      if (ranges::contains(tidxs, ix)) continue;  // target/output: not summed
      carriers[ix] |= (std::size_t{1} << i);
    }
    ++i;
  }
  for (auto&& [ix, cm] : carriers) {
    if (std::popcount(cm) < 2) continue;  // appears on < 2 tensors
    std::size_t rest = cm;
    while (rest) {
      std::size_t const b = static_cast<std::size_t>(std::countr_zero(rest));
      adj[b] |= (cm & ~(std::size_t{1} << b));  // neighbors, excluding self
      rest &= rest - 1;
    }
  }
  return adj;
}
```

Add includes near the top of `single_term_detail.hpp` if absent:
`#include <cstdlib>` and `#include <range/v3/algorithm/contains.hpp>`.

- [ ] **Step 4: Run test to verify it passes**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8 \
  && ./tests/unit/unit_tests-sequant "[pruning]"
```
Expected: builds; `contractible_adjacency` case passes.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/core/optimize/single_term_detail.hpp tests/unit/test_optimize.cpp
git add SeQuant/core/optimize/single_term_detail.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: add contractible_adjacency helper for DP pruning"
```

---

### Task 2: `connected_subsets` + `outer_product_connectivity`

**Files:**
- Modify: `SeQuant/core/optimize/single_term_detail.hpp` (add both helpers after
  `contractible_adjacency`)
- Test: `tests/unit/test_optimize.cpp`

**Interfaces:**
- Consumes: `contractible_adjacency` (Task 1).
- Produces:
  - `container::vector<char> connected_subsets(container::vector<std::size_t> const& adjacency, std::size_t nt)`
    -- `connected[n] == 1` iff the subgraph induced on the set bits of `n` is
    connected; size `1 << nt`; empty/singleton subsets are `1`.
  - `template <typename TIdxs> container::vector<char> outer_product_connectivity(TensorNetwork const& network, TIdxs const& tidxs)`
    -- returns `connected_subsets(contractible_adjacency(network, tidxs), nt)`,
    EXCEPT it returns an all-`1` mask (size `1 << nt`) when the env var
    `SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING` is set OR when the full-network entry
    (`mask.back()`) is `0` (a genuine multi-component product). So callers apply
    one uniform `if (!mask[n]) skip` test and get the unpruned behavior in both
    the disabled and product-term cases.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_optimize.cpp`:

```cpp
TEST_CASE("connected_subsets and outer_product_connectivity",
          "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;

  // Path 0-1-2 : adj[0]={1}, adj[1]={0,2}, adj[2]={1}.
  container::vector<std::size_t> path{0b010, 0b101, 0b010};
  auto c = o::connected_subsets(path, 3);
  REQUIRE(c.size() == 8);
  CHECK(c[0b001] == 1);  // singleton
  CHECK(c[0b011] == 1);  // {0,1} share edge
  CHECK(c[0b110] == 1);  // {1,2} share edge
  CHECK(c[0b101] == 0);  // {0,2} NOT adjacent -> disconnected
  CHECK(c[0b111] == 1);  // {0,1,2} connected via 1

  // Two disjoint components {0,1} and {2,3}: full set disconnected.
  container::vector<std::size_t> prod{0b0010, 0b0001, 0b1000, 0b0100};
  auto cp = o::connected_subsets(prod, 4);
  CHECK(cp[0b0011] == 1);   // {0,1}
  CHECK(cp[0b1100] == 1);   // {2,3}
  CHECK(cp[0b1111] == 0);   // full: disconnected product

  // outer_product_connectivity: env-disabled -> all ones.
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto prod = deserialize(L"f_{i1}^{a1} g_{a1}^{i2}",
                          {.def_perm_symm = Symmetry::Antisymm});
  container::vector<ExprPtr> v;
  for (auto&& f : prod->as<Product>().factors())
    if (f->is<Tensor>()) v.push_back(f);
  TensorNetwork tn{v};
  std::vector<Index> tgt{Index{L"i_1"}, Index{L"i_2"}};
  setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
  auto m_off = o::outer_product_connectivity(tn, tgt);
  unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
  for (auto val : m_off) CHECK(val == 1);
}
```

Add `#include <cstdlib>` to the test's includes for `setenv`/`unsetenv`.

- [ ] **Step 2: Run test to verify it fails**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8
```
Expected: compile error, `connected_subsets` / `outer_product_connectivity` not
declared.

- [ ] **Step 3: Write minimal implementation**

In `single_term_detail.hpp`, after `contractible_adjacency`:

```cpp
/// \brief `connected[n] == 1` iff the subgraph induced on the set bits of subset
/// mask `n` is connected under `adjacency` (from \ref contractible_adjacency).
/// Empty and singleton subsets are connected by definition. Size `1 << nt`.
inline container::vector<char> connected_subsets(
    container::vector<std::size_t> const& adjacency, std::size_t nt) {
  container::vector<char> connected(std::size_t{1} << nt, 0);
  for (std::size_t n = 0; n < connected.size(); ++n) {
    if (std::popcount(n) <= 1) {
      connected[n] = 1;
      continue;
    }
    std::size_t const start = n & (~n + 1);  // lowest set bit
    std::size_t reached = start, frontier = start;
    while (frontier) {
      std::size_t next = 0, f = frontier;
      while (f) {
        next |= adjacency[static_cast<std::size_t>(std::countr_zero(f))];
        f &= f - 1;
      }
      next &= n & ~reached;
      reached |= next;
      frontier = next;
    }
    connected[n] = (reached == n) ? 1 : 0;
  }
  return connected;
}

/// \brief Outer-product pruning mask for a term. Returns
/// \ref connected_subsets over \ref contractible_adjacency, EXCEPT it returns an
/// all-connected mask when pruning is disabled (env
/// `SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING`) or when the full network is itself
/// disconnected (a genuine product term -- never a CC residual summand, by the
/// linked-cluster theorem). Callers then apply one uniform `if (!mask[n]) skip`
/// test and get unpruned behavior in both those cases.
template <typename TIdxs>
inline container::vector<char> outer_product_connectivity(
    TensorNetwork const& network, TIdxs const& tidxs) {
  std::size_t const nt = network.tensors().size();
  std::size_t const sz = std::size_t{1} << nt;
  if (std::getenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING"))
    return container::vector<char>(sz, 1);
  auto mask = connected_subsets(contractible_adjacency(network, tidxs), nt);
  if (!mask.empty() && !mask.back())  // full network disconnected: disable
    return container::vector<char>(sz, 1);
  return mask;
}
```

- [ ] **Step 4: Run test to verify it passes**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8 \
  && ./tests/unit/unit_tests-sequant "[pruning]"
```
Expected: both `[pruning]` cases pass.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/core/optimize/single_term_detail.hpp tests/unit/test_optimize.cpp
git add SeQuant/core/optimize/single_term_detail.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: add connected_subsets + outer_product_connectivity"
```

---

### Task 3: Apply the mask in the DP driver and CSE metadata

**Files:**
- Modify: `SeQuant/core/optimize/cost_model.hpp` (`solve_single_term`
  ~lines 26-45; `AdditiveModel::build_context` ~lines 138-143)
- Modify: `SeQuant/core/optimize/single_term_detail.hpp` (`build_subnet_metadata`
  ~lines 435-463)
- Test: existing `[optimize]` suite (regression) + a smoke assertion added to
  the parity case is deferred to Task 4; here the test is the unchanged suite.

**Interfaces:**
- Consumes: `outer_product_connectivity` (Task 2).
- Produces: `build_subnet_metadata` gains a trailing parameter
  `container::vector<char> const& connected`; it skips subsets with
  `!connected[n]`. `solve_single_term` behavior: disconnected subsets are not
  formed and only connected-child bipartitions are relaxed (both governed by the
  uniform mask).

- [ ] **Step 1: Add the failing regression checkpoint**

No new test code; the gate is that the existing `[optimize]` suite still passes
after wiring. First confirm it passes BEFORE the change (baseline):

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -4 \
  && ./tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -5
```
Expected: all `[optimize]` assertions pass (record the counts).

- [ ] **Step 2: Modify `build_subnet_metadata` (add the skip)**

In `single_term_detail.hpp`, change the signature and loop:

```cpp
inline SubnetMetadata build_subnet_metadata(
    TensorNetwork const& network, container::vector<OptRes>& results,
    container::vector<char> const& connected) {   // NEW trailing param
  SubnetMetadata out;
  out.meta_ids.resize(results.size(), std::numeric_limits<size_t>::max());
  container::unordered_map<TensorNetwork::SlotCanonicalizationMetadata, size_t,
                           SubNetHash, SubNetEqual>
      meta_to_id;

  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    if (!connected[n]) continue;   // NEW: outer-product subset, never an intermediate
    auto ts = bits::on_bits_index(n) | bits::sieve(network.tensors());
    // ... unchanged body ...
  }
  out.unique_meta_costs.resize(meta_to_id.size(), 0.0);
  return out;
}
```

- [ ] **Step 3: Modify `AdditiveModel::build_context` (compute + pass the mask)**

In `cost_model.hpp` (~line 138):

```cpp
    if (subnet_cse) {
      auto connected = outer_product_connectivity(network, tidxs);  // NEW
      auto md = build_subnet_metadata(network, ctx.results, connected);  // pass mask
      ctx.meta_ids = std::move(md.meta_ids);
      ctx.unique_meta_costs = std::move(md.unique_meta_costs);
    }
```

- [ ] **Step 4: Modify `solve_single_term` (skip + guard)**

In `cost_model.hpp` (~line 26):

```cpp
template <class Model, typename TIdxs>
container::vector<typename Model::State> solve_single_term(
    Model const& m, TensorNetwork const& network, TIdxs const& tidxs,
    typename Model::Context& ctx) {
  auto const nt = network.tensors().size();
  container::vector<typename Model::State> st(size_t{1} << nt);
  auto const connected = outer_product_connectivity(network, tidxs);  // NEW
  for (size_t n = 1; n < st.size(); ++n) {
    if (std::popcount(n) == 1) {
      st[n] = m.leaf(ctx, n);
      continue;
    }
    if (!connected[n]) continue;  // NEW: never form a disconnected subset
    typename Model::State acc = m.init(ctx, n);
    for (auto&& [lp, rp] : bits::bipartitions(n))
      if (lp != 0 && rp != 0 && connected[lp] && connected[rp])  // NEW guard
        m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
    st[n] = std::move(acc);
    m.finalize(ctx, n, st);
  }
  return st;
}
```

Remove the now-unused `(void)tidxs;` line (tidxs is used now).

- [ ] **Step 5: Build and run the regression + verify parity by env toggle**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8 \
  && ./tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -5
```
Expected: identical pass counts to Step 1 (pruning is loss-free, so every
existing optimize assertion is unchanged).

- [ ] **Step 6: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/core/optimize/single_term_detail.hpp SeQuant/core/optimize/cost_model.hpp
git add SeQuant/core/optimize/single_term_detail.hpp SeQuant/core/optimize/cost_model.hpp
git commit -m "optimize: prune outer-product subsets in the single-term DP"
```

---

### Task 4: Parity gate across all objectives (primary correctness test)

**Files:**
- Test: `tests/unit/test_optimize.cpp`

**Interfaces:**
- Consumes: `optimize(ExprPtr, OptimizeOptions)` (public) and the env switch
  `SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING` (Task 2).

- [ ] **Step 1: Write the parity test**

Add to `tests/unit/test_optimize.cpp`. It optimizes each term under each
objective twice -- pruning on (default) and off (env) -- and asserts the
optimized expressions are identical (`to_latex`). Includes a connected hyperedge
term so the "Hadamard intermediate is kept" path is exercised.

```cpp
TEST_CASE("outer-product pruning parity (pruned == unpruned)",
          "[optimize][pruning]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);

  auto opts_for = [&](ObjectiveFunction obj) {
    OptimizeOptions o;
    o.objective_function = obj;
    o.idx_to_extent = [](Index const& ix) -> std::size_t {
      return ix.nonnull() ? ix.space().approximate_size() : 1;
    };
    o.batch_policy.is_batchable_index = [](Index const&) { return false; };
    o.batch_policy.batch_target_size = [](Index const&) -> std::size_t {
      return 1;
    };
    return o;
  };

  std::vector<ObjectiveFunction> const objs{
      ObjectiveFunction::DenseFLOPs,
      ObjectiveFunction::DenseSize,
      ObjectiveFunction::DenseSpaceTime,
      ObjectiveFunction::DenseTimeSpace,
      ObjectiveFunction::DenseSpaceTimeBatched,
      ObjectiveFunction::DenseTimeSpaceBatched};

  std::vector<std::wstring> const terms{
      L"1/4 g_{i2,i3}^{a2,a3} t_{a2,a3}^{i2,i3} t_{a1}^{i1}",
      L"x_{i1,i2}^{a3,a4} y_{a1,a2}^{i1,i2} z_{a3,a4}^{a1,a2}",
      L"g_{i1,i2}^{a1,a2} t_{a1}^{i1} t_{a2}^{i2} t_{a3}^{i3}",
      // hyperedge-flavored: three tensors sharing a common summed index a5
      L"p_{i1}^{a5} q_{i2}^{a5} r_{i3}^{a5} s_{a5}^{i4}",
  };

  auto run = [&](std::wstring const& term, ObjectiveFunction obj, bool prune) {
    if (prune)
      unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
    else
      setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
    auto expr = deserialize(term, {.def_perm_symm = Symmetry::Antisymm});
    auto out = optimize(expr, opts_for(obj));
    unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
    return to_latex(out);
  };

  for (std::size_t ti = 0; ti < terms.size(); ++ti)
    for (auto obj : objs) {
      auto pruned = run(terms[ti], obj, true);
      auto unpruned = run(terms[ti], obj, false);
      CAPTURE(ti, static_cast<int>(obj));
      CHECK(pruned == unpruned);
    }
}
```

- [ ] **Step 2: Run and verify it passes**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8 \
  && ./tests/unit/unit_tests-sequant "[pruning]"
```
Expected: all parity checks pass. If any objective's opts need a non-trivial
`inner_pow` or roofline to run, set them in `opts_for` (mirror the defaults used
in the existing `"optimize"` test's `single_term_opt` lambda). A CHECK failure
means a genuine factorization divergence -- investigate before proceeding
(do not weaken the assertion).

- [ ] **Step 3: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i tests/unit/test_optimize.cpp
git add tests/unit/test_optimize.cpp
git commit -m "optimize: parity test for outer-product pruning across objectives"
```

---

### Task 5: Performance smoke + multi-component safety net

**Files:**
- Test: `tests/unit/test_optimize.cpp`

**Interfaces:**
- Consumes: `optimize` + the env switch.

- [ ] **Step 1: Write the tests**

Add to `tests/unit/test_optimize.cpp`:

```cpp
TEST_CASE("outer-product pruning: large connected term optimizes quickly",
          "[optimize][pruning]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  // A connected chain: each adjacent pair shares one summed index, so the whole
  // term is one connected component. The pruned DP explores only connected
  // subsets and finishes fast; the unpruned 3^n enumeration is far slower.
  auto expr = deserialize(
      L"g_{a0,a1}^{a2,a3} t_{a2}^{a4} t_{a3}^{a5} t_{a4}^{a6} t_{a5}^{a7} "
      L"t_{a6}^{a8} t_{a7}^{a9} v_{a8,a9}^{a0,a1}",
      {.def_perm_symm = Symmetry::Antisymm});
  OptimizeOptions o;
  o.objective_function = ObjectiveFunction::DenseTimeSpace;
  o.idx_to_extent = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : 1;
  };
  // Pruning ON (default): must complete (the assertion is that it returns).
  auto out = optimize(expr, o);
  CHECK(out);
}

TEST_CASE("outer-product pruning: multi-component product falls back unpruned",
          "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto net_of = [](ExprPtr const& p) {
    container::vector<ExprPtr> v;
    for (auto&& f : p->as<Product>().factors())
      if (f->is<Tensor>()) v.push_back(f);
    return TensorNetwork{v};
  };

  // Two independent contractions sharing NO summed index: {f,g} over a1,
  // {p,q} over a2. The full adjacency graph has two components.
  auto prod = deserialize(L"f_{i1}^{a1} g_{a1}^{i2} p_{i3}^{a2} q_{a2}^{i4}",
                          {.def_perm_symm = Symmetry::Antisymm});
  auto tn = net_of(prod);
  std::vector<Index> tgt{Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"},
                         Index{L"i_4"}};
  auto mask = o::outer_product_connectivity(tn, tgt);
  for (auto v : mask) CHECK(v == 1);  // disconnected full net -> all-connected

  // And optimize() on the product term must match the env-disabled result.
  OptimizeOptions opt;
  opt.objective_function = ObjectiveFunction::DenseFLOPs;
  opt.idx_to_extent = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : 1;
  };
  auto with = to_latex(optimize(prod, opt));
  setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
  auto without = to_latex(optimize(prod, opt));
  unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
  CHECK(with == without);
}
```

- [ ] **Step 2: Run and verify**

```
cd SeQuant/build-test && make unit_tests-sequant 2>&1 | tail -8 \
  && ./tests/unit/unit_tests-sequant "[pruning]"
```
Expected: both cases pass; the large-term case returns promptly (well under the
test's normal budget).

- [ ] **Step 3: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i tests/unit/test_optimize.cpp
git add tests/unit/test_optimize.cpp
git commit -m "optimize: perf smoke + multi-component safety net for DP pruning"
```

---

## Post-implementation validation (final review, not a task)

After all tasks pass, confirm the real-world win on the C60 term that motivated
this (needs the TA/dryrun build, not `build-test`): the perf-first optimize of
the 13-tensor residual term that previously did not finish in 38 minutes now
completes in seconds. Re-run the local `[dryrun-perfcost]` / batchability-audit
spike (uncommitted, in `test_eval_dryrun.cpp`) and record the wall time. This is
the motivating evidence, verified once at the end; it is not a committed test.

## Self-review notes

- **Spec coverage:** adjacency (Task 1) = spec "Adjacency: share a contractible
  (non-target) index"; connectivity + toggle (Task 2) = spec "Precompute" +
  "Product terms"; driver + CSE skip (Task 3) = spec "The pruning rule" +
  "CSE metadata skips disconnected subsets"; parity (Task 4) = spec "Testing (1)"
  and the "Correctness / loss-free" gate; perf + product safety net (Task 5) =
  spec "Testing (3),(4)".
- The plan folds the spec's explicit `full_connected` guard into
  `outer_product_connectivity` (returns an all-connected mask when the full
  network is disconnected), so the driver/CSE apply one uniform `!mask[n]` test.
  Behavior is identical to the spec sketch.
- Multi-component per-component optimization is NOT implemented (spec non-goal;
  CC residual summands are single-component by linked-cluster). Product terms
  take the unpruned path via the mask.
