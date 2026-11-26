#include "evaluator.hpp"

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/tensor_power.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <random>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace sequant::deriv::test {

Expressions summands_of(const ExprPtr& e) {
  if (!e->is<Sum>()) return {e};
  Expressions result;
  for (const auto& s : e->as<Sum>().summands())
    for (const auto& t : summands_of(s)) result.push_back(t);
  return result;
}

Expressions factors_of(const ExprPtr& term, Constant::scalar_type& scalar) {
  if (!term->is<Product>()) return {term};
  scalar *= term->as<Product>().scalar();
  Expressions result;
  for (const auto& f : term->as<Product>().factors())
    for (const auto& g : factors_of(f, scalar)) result.push_back(g);
  return result;
}

bool is_projector(const Expr& e) {
  return e.is<Tensor>() &&
         (e.as<Tensor>().label() == reserved::antisymm_label() ||
          e.as<Tensor>().label() == reserved::symm_label());
}

ProductIndices indices_of_product(const ExprPtr& product) {
  ProductIndices result;
  const auto counts = get_used_indices_with_counts(*product);
  Constant::scalar_type scalar = 1;
  for (const auto& f : factors_of(product, scalar))
    if (is_projector(*f))
      for (const auto& i : f->as<AbstractTensor>()._slots())
        result.projected.insert(i);
  for (const auto& [i, count] : counts)
    if (count.nonproto() != 2) result.free.insert(i);
  container::set<Index> carried;
  for (const auto& i : result.free)
    carried.insert(i.proto_indices().begin(), i.proto_indices().end());
  result.free.merge(carried);
  return result;
}

Expressions products_of(const ExprPtr& term) {
  Constant::scalar_type scalar = 1;
  if (std::ranges::none_of(factors_of(term, scalar),
                           [](const ExprPtr& f) { return f->is<Sum>(); }))
    return {term};
  auto expanded = term->clone();
  expand(expanded);
  return summands_of(expanded);
}

container::set<Index> externals_of(const ExprPtr& f) {
  const auto product = products_of(summands_of(f).front()).front();
  auto [free, projected] = indices_of_product(product);
  free.merge(projected);
  return free;
}

container::set<Index> indices_of_leaf(const ExprPtr& x) {
  return get_used_indices(*x);
}

/// the symbol of a leaf, apart from its indices: a Variable's label and
/// conjugation, a Tensor's label and ranks
std::wstring symbol_of(const Expr& leaf) {
  if (leaf.is<Variable>()) {
    const auto& v = leaf.as<Variable>();
    return std::wstring(v.label()) + (v.conjugated() ? L"^*" : L"");
  }
  const auto& t = leaf.as<AbstractTensor>();
  return std::wstring(t._label()) + L"{" + std::to_wstring(t._bra_rank()) +
         L"," + std::to_wstring(t._ket_rank()) + L"," +
         std::to_wstring(t._aux_rank()) + L"}";
}

container::map<Index, Index> names_for(
    const container::set<Index>& indices,
    const std::function<std::optional<Index>(const Index&)>& name) {
  container::map<Index, Index> result;
  for (const auto& i : indices)
    if (!i.has_proto_indices())
      if (auto n = name(i)) result.emplace(i, *n);
  for (const auto& i : indices) {
    if (!i.has_proto_indices()) continue;
    const auto own = name(i.drop_proto_indices());
    const bool carried = std::ranges::any_of(
        i.proto_indices(), [&](const Index& p) { return result.contains(p); });
    if (!own && !carried) continue;
    container::vector<Index> protos;
    for (const auto& p : i.proto_indices())
      protos.push_back(result.contains(p) ? result.at(p) : p);
    result.emplace(i, Index(own ? *own : i.drop_proto_indices(), protos,
                            i.symmetric_proto_indices()));
  }
  return result;
}

namespace {

/// @p f with every index named like an index or protoindex of @p x renamed
/// to a temporary named unlike any index of x and of f, protoindices
/// following
ExprPtr with_dummies_apart(const ExprPtr& f, const ExprPtr& x) {
  container::set<std::wstring> in_x;
  for (const auto& i : get_used_indices(*x)) in_x.emplace(i.label());
  const auto in_f = get_used_indices(*f);
  auto taken = in_x;
  for (const auto& i : in_f) taken.emplace(i.label());
  const auto unused = unused_label_validator(std::move(taken));
  const auto apart =
      names_for(in_f, [&in_x, &unused](const Index& i) -> std::optional<Index> {
        if (!in_x.contains(std::wstring(i.label()))) return std::nullopt;
        Index result;
        do result = Index::make_tmp_index(i.basis());
        while (!unused(result));
        return result;
      });
  return apart.empty() ? f : transform_expr(f, apart);
}

}  // namespace

namespace {

/// the values an index of @p space runs over: @p occupied per occupied and
/// three per virtual elementary space, its quantum numbers ignored unless
/// @p quantum_numbers
container::svector<int> range_of(const IndexSpace& space, bool quantum_numbers,
                                 int occupied) {
  container::svector<int> result;
  const auto& registry = *get_default_context().index_basis_registry();
  const auto types = space.type().to_int32();
  const auto qns = quantum_numbers ? space.qns().to_int32() : 0;
  for (int t = 0; t != 31; ++t) {
    if (!((types >> t) & 1)) continue;
    const IndexSpace elementary(L"elementary", TypeAttr(1 << t), space.qns());
    const int count = registry.contains_unoccupied(elementary) ? 3 : occupied;
    for (int q = 0; q != 31; ++q) {
      // a space without quantum numbers is one elementary space per type
      if (qns ? !((qns >> q) & 1) : q != 0) continue;
      for (int k = 0; k != count; ++k) result.push_back(128 * t + 4 * q + k);
    }
  }
  return result;
}

}  // namespace

void add_domain(const Index& i, const Values& values, Domains& domains,
                bool quantum_numbers, int occupied) {
  const std::wstring label(i.label());
  if (values.contains(label) ||
      std::ranges::any_of(domains,
                          [&](const auto& d) { return d.first == label; }))
    return;
  domains.emplace_back(label, range_of(i.space(), quantum_numbers, occupied));
}

namespace {

/// what a comparison made of an element it read
enum class Draw {
  /// nonzero in either side, and equal
  informative,
  /// zero in both sides, or no element (a draw that vanishes by symmetry)
  empty,
  /// the sides differ: stop
  stop
};

/// the most draws a sample takes per informative element it asks for
constexpr std::size_t draws_per_element = 8;

/// calls @p f with each assignment of @p domains , in their order, or, where
/// there are more than sampling.count, with assignments drawn at random until
/// sampling.count of them are informative or draws_per_element times as many
/// are drawn; until @p f says stop. A drawn label takes, half the time, the
/// value of an earlier label over its range of another group, else a value no
/// earlier label of its own group over that range has, if one is left: an
/// element where an antisymmetric tensor vanishes is avoided, one where a
/// Kronecker delta between two tensors holds is met. @p groups gives the group
/// of each label
/// @return whether every assignment was read
bool for_each_sampled(
    const Domains& domains, const container::svector<int>& groups,
    const Sampling& sampling,
    const std::function<Draw(const container::svector<int>&)>& f) {
  container::svector<int> numbers(domains.size());
  std::size_t assignments = 1;
  for (const auto& [label, range] : domains)
    assignments = std::min(assignments * range.size(), sampling.count + 1);
  if (sampling.count != 0 && assignments > sampling.count) {
    std::mt19937_64 random(sampling.seed);
    std::size_t informative = 0;
    const auto draws =
        sampling.draws ? sampling.draws : draws_per_element * sampling.count;
    for (std::size_t n = 0; n != draws && informative != sampling.count; ++n) {
      for (std::size_t k = 0; k != domains.size(); ++k) {
        const auto& range = domains[k].second;
        container::svector<int> others, taken;
        for (std::size_t j = 0; j != k; ++j)
          if (domains[j].second == range)
            (groups[j] == groups[k] ? taken : others).push_back(numbers[j]);
        if (!others.empty() && random() % 2) {
          numbers[k] = others[random() % others.size()];
          continue;
        }
        container::svector<int> free;
        for (const auto v : range)
          if (std::ranges::find(taken, v) == taken.end()) free.push_back(v);
        const auto& from = free.empty() ? range : free;
        numbers[k] = from[random() % from.size()];
      }
      switch (f(numbers)) {
        case Draw::stop:
          return false;
        case Draw::informative:
          ++informative;
          break;
        case Draw::empty:
          break;
      }
    }
    return false;
  }
  bool going = true;
  std::function<void(std::size_t)> each = [&](std::size_t k) {
    if (!going) return;
    if (k == domains.size()) {
      going = f(numbers) != Draw::stop;
      return;
    }
    for (const auto v : domains[k].second) {
      numbers[k] = v;
      each(k + 1);
    }
  };
  each(0);
  return true;
}

/// the group of each label of @p domains : the position in @p groups of the
/// set that holds an index of that label or one it is a protoindex of; 0 if
/// none does
container::svector<int> groups_of(
    const Domains& domains,
    const container::svector<container::set<Index>>& groups) {
  container::svector<int> result(domains.size(), 0);
  for (std::size_t k = 0; k != domains.size(); ++k)
    for (std::size_t g = 0; g != groups.size(); ++g)
      for (const auto& i : groups[g]) {
        bool holds = i.label() == domains[k].first;
        for (const auto& p : i.proto_indices())
          holds = holds || p.label() == domains[k].first;
        if (holds) result[k] = static_cast<int>(g);
      }
  return result;
}

using Number = std::complex<double>;

std::uint64_t mixed(std::uint64_t h) {
  h += 0x9e3779b97f4a7c15ULL;
  h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
  h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
  return h ^ (h >> 31);
}

std::uint64_t combined(std::uint64_t a, std::uint64_t b) {
  return mixed(a ^ mixed(b));
}

std::uint64_t hash_of(std::wstring_view text) {
  return std::hash<std::wstring_view>{}(text);
}

/// a number in [-1, 1) per hash, real if @p real
Number number_of(std::uint64_t h, bool real) {
  auto unit = [](std::uint64_t x) {
    return static_cast<double>(x >> 11) * 0x1.0p-52 - 1;
  };
  const auto a = mixed(h), b = mixed(a);
  return {unit(a), real ? 0. : unit(b)};
}

/// a slot of a factor: the label numbers of its index and protoindices
struct Slot {
  int own = 0;
  container::svector<int> protos = {};
  bool symmetric = false;
};

using Assignment = std::vector<int>;
using Addresses = container::svector<std::uint64_t>;

/// the values of the protoindices of @p s , sorted if symmetric
container::svector<int> protos_at(const Slot& s, const Assignment& v) {
  container::svector<int> result;
  for (const auto p : s.protos) result.push_back(v[p]);
  if (s.symmetric) std::ranges::sort(result);
  return result;
}

std::uint64_t address(const Slot& s, const Assignment& v) {
  std::uint64_t result = v[s.own] + 1;
  if (s.protos.empty()) return result;
  for (const auto p : protos_at(s, v)) result = result * 4099 + p + 1;
  return result;
}

/// the slots of a data factor at an assignment, arranged into the one order
/// its symmetry allows, as the data is stored
struct Canonical {
  Addresses b, k, x;
  int sign = 1;
  /// an antisymmetric element with two equal slots
  bool vanishes = false;
  /// read at bra and ket exchanged; its own exchange
  bool swapped = false, diagonal = false;
};

/// the values of the labels a factor reads, packed: 12 bits per label
using Key = std::array<std::uint64_t, 2>;

struct KeyHash {
  std::size_t operator()(const Key& k) const {
    return std::hash<std::uint64_t>{}(k[0] ^ (k[1] * 0x9e3779b97f4a7c15ULL));
  }
};

/// a data factor at an assignment of its labels: its canonical slots and its
/// stored element
struct Stored {
  Canonical canonical;
  std::complex<double> datum;
};

/// a factor of a term, read once and evaluated at each assignment
struct Factor {
  enum class Kind { number, kronecker, overlap, data } kind = Kind::number;
  Number number = {};  // Kind::number
  /// Kind::number: a Variable, its value unconjugated
  bool variable = false;
  std::uint64_t key = 0;
  /// the conjugate spelling: the conjugate of the data at the slots stored,
  /// or of the value of the Variable
  bool conjugated = false;
  bool real = false;
  Symmetry symmetry = Symmetry::Nonsymm;
  BraKetSymmetry braket = BraKetSymmetry::Nonsymm;
  ColumnSymmetry column = ColumnSymmetry::Nonsymm;
  container::svector<Slot> bra = {}, ket = {}, aux = {};
  /// an overlap of unit metric within one basis instance, its protoindices
  /// aside: is_identity_at() decides by their values
  bool unit_overlap = false;
  /// a Power or TensorPower: the element raised, conjugated if flagged
  std::optional<rational> exponent = {};
  bool power_conjugated = false;
  /// the label numbers of its slots and their protoindices, once each
  container::svector<int> labels = {};
  /// what stored() read, per assignment of the labels
  mutable std::unordered_map<Key, Stored, KeyHash> read = {};
};

/// the label numbers of the slots of @p f and of their protoindices, once each
container::svector<int> labels_of(const Factor& f) {
  container::svector<int> result;
  auto add = [&result](int n) {
    if (std::ranges::find(result, n) == result.end()) result.push_back(n);
  };
  for (const auto* slots : {&f.bra, &f.ket, &f.aux})
    for (const auto& s : *slots) {
      add(s.own);
      for (const auto p : s.protos) add(p);
    }
  return result;
}

/// whether the overlap @p f is the identity at @p v : it has unit metric and
/// the values of its protoindices agree, as multisets where they are
/// symmetric, else in order; elsewhere it is data with its declared symmetry
bool is_identity_at(const Factor& f, const Assignment& v) {
  return f.unit_overlap && protos_at(f.bra[0], v) == protos_at(f.ket[0], v);
}

/// @p b and @p k arranged as @p f 's symmetry allows into one order, with the
/// sign of the arrangement; false if an antisymmetric element vanishes
bool arranged(const Factor& f, Addresses& b, Addresses& k, int& sign) {
  auto sorted = [&sign](Addresses& a, bool signed_) {
    for (std::size_t i = 1; i < a.size(); ++i)
      for (std::size_t j = i; j > 0 && a[j] < a[j - 1]; --j) {
        std::swap(a[j], a[j - 1]);
        if (signed_) sign = -sign;
      }
    return !signed_ || std::ranges::adjacent_find(a) == a.end();
  };
  switch (f.symmetry) {
    case Symmetry::Antisymm:
      return sorted(b, true) && sorted(k, true);
    case Symmetry::Symm:
      sorted(b, false);
      sorted(k, false);
      return true;
    case Symmetry::Nonsymm:
      break;
  }
  if (f.column != ColumnSymmetry::Symm) return true;
  const auto n = std::min(b.size(), k.size());
  container::svector<std::pair<std::uint64_t, std::uint64_t>> columns;
  for (std::size_t c = 0; c != n; ++c) columns.emplace_back(b[c], k[c]);
  std::ranges::sort(columns);
  for (std::size_t c = 0; c != n; ++c) std::tie(b[c], k[c]) = columns[c];
  std::sort(b.begin() + n, b.end());
  std::sort(k.begin() + n, k.end());
  return true;
}

Number raised(Number z, const rational& exponent) {
  if (denominator(exponent) != 1)
    return std::pow(z, exponent.convert_to<double>());
  const auto n = numerator(exponent).convert_to<long>();
  Number result = 1;
  for (long k = 0; k != std::abs(n); ++k) result *= z;
  return n < 0 ? 1. / result : result;
}

/// a number and its derivatives by the perturbed element z and by its
/// conjugate z̄, as Wirtinger derivatives: exact for the polynomials of the
/// equations
struct Dual {
  Number value = 0, by = 0, by_conjugate = 0;
  Dual() = default;
  Dual(Number v) : value(v) {}
  Dual(Number v, Number d, Number d_bar)
      : value(v), by(d), by_conjugate(d_bar) {}
};

Dual operator*(const Dual& a, const Dual& b) {
  return {a.value * b.value, a.value * b.by + a.by * b.value,
          a.value * b.by_conjugate + a.by_conjugate * b.value};
}
Dual& operator*=(Dual& a, const Dual& b) { return a = a * b; }
Dual& operator+=(Dual& a, const Dual& b) {
  a.value += b.value;
  a.by += b.by;
  a.by_conjugate += b.by_conjugate;
  return a;
}

/// d conj(w)/dz = conj(dw/dz̄)
Dual conj(const Dual& a) {
  return {std::conj(a.value), std::conj(a.by_conjugate), std::conj(a.by)};
}
Number conj(Number a) { return std::conj(a); }

bool is_zero(Number z) { return z == 0.; }
bool is_zero(const Dual& z) {
  return z.value == 0. && z.by == 0. && z.by_conjugate == 0.;
}

Dual raised(const Dual& z, const rational& exponent) {
  const auto value = raised(z.value, exponent);
  if (z.by == 0. && z.by_conjugate == 0.) return value;
  const auto slope =
      exponent.convert_to<double>() * raised(z.value, exponent - 1);
  return {value, slope * z.by, slope * z.by_conjugate};
}

Canonical canonical_of(const Factor& f, const Assignment& v) {
  Canonical c;
  for (const auto& s : f.bra) c.b.push_back(address(s, v));
  for (const auto& s : f.ket) c.k.push_back(address(s, v));
  for (const auto& s : f.aux) c.x.push_back(address(s, v));
  if (!arranged(f, c.b, c.k, c.sign)) {
    c.vanishes = true;
    return c;
  }
  if (f.braket != BraKetSymmetry::Nonsymm) {
    auto b2 = c.k, k2 = c.b;
    int sign2 = c.sign;
    arranged(f, b2, k2, sign2);
    const auto first = std::tie(c.b, c.k), second = std::tie(b2, k2);
    const auto ranks = std::pair(c.b.size(), c.k.size()),
               ranks2 = std::pair(b2.size(), k2.size());
    if (std::tie(ranks2, second) < std::tie(ranks, first)) {
      std::swap(c.b, b2);
      std::swap(c.k, k2);
      std::swap(c.sign, sign2);
      c.swapped = true;
    } else if (ranks == ranks2 && first == second) {
      // the element is its own exchanged element: real if Hermitian
      c.diagonal = true;
    }
  }
  return c;
}

/// the stored element of @p f at @p c , unconjugated
Number datum(const Factor& f, const Canonical& c) {
  if (c.vanishes) return 0;
  auto h = combined(f.key, c.b.size() * 1024 + c.k.size());
  for (const auto& a : {&c.b, &c.k, &c.x}) {
    h = combined(h, a->size());
    for (const auto e : *a) h = combined(h, e);
  }
  auto z = number_of(h, f.real || f.braket == BraKetSymmetry::Symm);
  if (f.braket == BraKetSymmetry::Conjugate) {
    if (c.diagonal) z = z.real();
    if (c.swapped) z = std::conj(z);
  }
  return z * static_cast<double>(c.sign);
}

/// canonical_of() and datum() of @p f at @p v , computed once per assignment
/// of its labels (into @p scratch where they do not fit a Key)
const Stored& stored(const Factor& f, const Assignment& v, Stored& scratch) {
  Key key{};
  bool packed = f.labels.size() <= 10;
  for (std::size_t j = 0; packed && j != f.labels.size(); ++j) {
    const auto value = v[f.labels[j]];
    packed = value >= 0 && value < 4096;
    key[j / 5] |= static_cast<std::uint64_t>(value) << (12 * (j % 5));
  }
  auto compute = [&f, &v](Stored& s) {
    s.canonical = canonical_of(f, v);
    s.datum = datum(f, s.canonical);
  };
  if (!packed) {
    compute(scratch);
    return scratch;
  }
  const auto [it, added] = f.read.try_emplace(key);
  if (added) compute(it->second);
  return it->second;
}

/// one element of a symbol raised by ε, with its symmetry partners
struct Perturbation {
  std::uint64_t key = 0;
  /// a Variable; else the element at these canonical slots
  bool variable = false;
  Addresses b = {}, k = {}, x = {};
  int sign = 1;
  bool swapped = false;
  /// the perturbed spelling is the conjugate one: its derivative is by z̄
  bool conjugated = false;
  /// a real symbol: z̄ is z
  bool real = false;
};

/// whether the data factor @p f at @p c is the perturbed element of @p p
bool is_perturbed(const Factor& f, const Canonical& c, const Perturbation& p) {
  return !p.variable && !c.vanishes && f.key == p.key && c.b == p.b &&
         c.k == p.k && c.x == p.x;
}

template <typename N>
N value_of(const Factor& f, const Assignment& v, const Perturbation* p) {
  N z{};
  auto data = [&]() -> N {
    Stored scratch;
    const auto& [c, element] = stored(f, v, scratch);
    N result = element;
    if constexpr (std::is_same_v<N, Dual>) {
      if (p && is_perturbed(f, c, *p)) {
        const Number slope = static_cast<double>(c.sign * p->sign);
        // the elements read at the other exchange of a Hermitian tensor are
        // conjugates of the perturbed one
        if (f.braket == BraKetSymmetry::Conjugate && c.swapped != p->swapped)
          result.by_conjugate = slope;
        else
          result.by = slope;
      }
    }
    return f.conjugated ? conj(result) : result;
  };
  switch (f.kind) {
    case Factor::Kind::number:
      z = f.number;
      if (f.variable) {
        if constexpr (std::is_same_v<N, Dual>)
          if (p && p->variable && p->key == f.key) z.by = 1;
        if (f.conjugated) z = conj(z);
      }
      break;
    case Factor::Kind::kronecker:
      z = Number(address(f.bra[0], v) == address(f.ket[0], v) ? 1. : 0.);
      break;
    case Factor::Kind::overlap:
      if (is_identity_at(f, v)) {
        z = Number(v[f.bra[0].own] == v[f.ket[0].own] ? 1. : 0.);
        break;
      }
      z = data();
      break;
    case Factor::Kind::data:
      z = data();
      break;
  }
  if (f.exponent) z = raised(z, *f.exponent);
  return f.power_conjugated ? conj(z) : z;
}

/// a projector of a term: the normalized sum over the permutations of its
/// slots, signed for Â, of columns for Ŝ
struct Projection {
  bool antisymmetric = true;
  container::svector<int> bra = {}, ket = {};
};

/// every permutation of 0, …, n-1 with its sign, by core's
/// next_permutation_parity
container::svector<std::pair<container::svector<std::size_t>, int>>
signed_permutations(std::size_t n) {
  container::svector<std::size_t> order(n);
  for (std::size_t k = 0; k != n; ++k) order[k] = k;
  container::svector<std::pair<container::svector<std::size_t>, int>> result;
  int parity = 0;
  do {
    result.emplace_back(order, parity ? -1 : 1);
  } while (next_permutation_parity(parity, order.begin(), order.end()));
  return result;
}

}  // namespace

container::svector<Rearrangement> rearrangements_of(const AbstractTensor& t,
                                                    bool real) {
  using Slots = container::svector<Index>;
  Slots bra, ket, aux, all;
  for (const auto& i : t._bra()) bra.push_back(i);
  for (const auto& i : t._ket()) ket.push_back(i);
  for (const auto& i : t._aux()) aux.push_back(i);
  for (const auto& i : t._slots()) all.push_back(i);
  container::svector<std::pair<Slots, Slots>> orientations{{bra, ket}};
  const auto braket = t._braket_symmetry();
  // bra and ket of unequal ranks, or of rank 0, exchanged address the same
  // stored element
  if (bra.size() == ket.size() && !bra.empty() &&
      (braket == BraKetSymmetry::Symm ||
       (real && braket == BraKetSymmetry::Conjugate)))
    orientations.emplace_back(ket, bra);
  auto permuted = [](const Slots& slots,
                     const container::svector<std::size_t>& by) {
    Slots result;
    for (const auto k : by) result.push_back(slots[k]);
    return result;
  };
  // the protoindices that are no slots, and their orders that keep the
  // protoindex list of every slot
  Slots protoindices;
  for (const auto& i : all)
    for (const auto& p : i.proto_indices())
      if (std::ranges::find(all, p) == all.end() &&
          std::ranges::find(protoindices, p) == protoindices.end())
        protoindices.push_back(p);
  container::svector<Slots> orders;
  for (const auto& [by, sign] : signed_permutations(protoindices.size())) {
    const auto order = permuted(protoindices, by);
    const bool keeps = std::ranges::all_of(all, [&](const Index& i) {
      Slots before(i.proto_indices().begin(), i.proto_indices().end()), after;
      for (const auto& p : before) {
        const auto k = static_cast<std::size_t>(
            std::ranges::find(protoindices, p) - protoindices.begin());
        after.push_back(k == protoindices.size() ? p : order[k]);
      }
      if (i.symmetric_proto_indices()) {
        std::ranges::sort(before);
        std::ranges::sort(after);
      }
      return before == after;
    });
    if (keeps) orders.push_back(order);
  }
  container::svector<Rearrangement> result;
  auto emit = [&](const Slots& b, const Slots& k, int sign) {
    Slots slots = b;
    slots.insert(slots.end(), k.begin(), k.end());
    slots.insert(slots.end(), aux.begin(), aux.end());
    for (const auto& order : orders)
      result.push_back({.slots = slots, .protoindices = order, .sign = sign});
  };
  for (const auto& [b, k] : orientations) {
    if (t._symmetry() != Symmetry::Nonsymm) {
      const bool anti = t._symmetry() == Symmetry::Antisymm;
      for (const auto& [by_bra, bra_sign] : signed_permutations(b.size()))
        for (const auto& [by_ket, ket_sign] : signed_permutations(k.size()))
          emit(permuted(b, by_bra), permuted(k, by_ket),
               anti ? bra_sign * ket_sign : 1);
    } else if (t._column_symmetry() == ColumnSymmetry::Symm) {
      const auto paired = std::min(b.size(), k.size());
      // the paired columns permuted by p, the unpaired slots by u
      auto apply = [paired](const Slots& slots,
                            const container::svector<std::size_t>& p,
                            const container::svector<std::size_t>& u) {
        Slots result;
        for (const auto j : p) result.push_back(slots[j]);
        for (std::size_t j = paired; j < slots.size(); ++j)
          result.push_back(slots[paired + u[j - paired]]);
        return result;
      };
      const auto unpaired = std::max(b.size(), k.size()) - paired;
      for (const auto& [p, p_sign] : signed_permutations(paired))
        for (const auto& [u, u_sign] : signed_permutations(unpaired))
          emit(apply(b, p, u), apply(k, p, u), 1);
    } else {
      emit(b, k, 1);
    }
  }
  return result;
}

namespace {

/// a term read for evaluation: label numbers below `externals` are the
/// external labels, the k-th summed label has number externals + k
struct Term {
  Number scalar = 1;
  container::svector<Factor> factors = {};
  /// the projectors applied to the term, over external labels
  container::svector<Projection> projections = {};
  std::size_t externals = 0;
  container::svector<container::svector<int>> domains = {};
  /// the factors that can be evaluated before any summed label is assigned,
  /// and once the k-th is
  container::svector<std::size_t> first = {};
  container::svector<container::svector<std::size_t>> after = {};
};

/// reads the terms of an expression with fixed external labels
class Reader {
 public:
  Reader(const container::svector<std::wstring>& externals, const Data& data,
         bool quantum_numbers)
      : data_(data), quantum_numbers_(quantum_numbers) {
    for (const auto& l : externals)
      externals_.emplace(l, static_cast<int>(externals_.size()));
  }

  container::svector<Term> terms(const ExprPtr& e) {
    container::svector<Term> result;
    for (const auto& summand : summands_of(e)) {
      Constant::scalar_type scalar = 1;
      const auto factors = factors_of(summand, scalar);
      if (std::ranges::any_of(factors,
                              [](const ExprPtr& f) { return f->is<Sum>(); })) {
        auto expanded = summand->clone();
        expand(expanded);
        for (auto& t : terms(expanded)) result.push_back(std::move(t));
        continue;
      }
      result.push_back(term(Constant(scalar).value<Number>(), factors));
    }
    return result;
  }

 private:
  const Data& data_;
  bool quantum_numbers_;
  container::map<std::wstring, int> externals_;
  // per term
  container::map<std::wstring, int> summed_;
  container::svector<Index> summed_indices_;

  Term term(Number scalar, const Expressions& factors) {
    summed_.clear();
    summed_indices_.clear();
    Term t{.scalar = scalar, .externals = externals_.size()};
    // summed labels of the δ's and unit overlaps first, so that they cut the
    // sum early
    for (int pass = 0; pass != 2; ++pass)
      for (const auto& f : factors) {
        const bool kronecker =
            f->is<Tensor>() &&
            (f->as<Tensor>().label() == reserved::kronecker_label() ||
             f->as<Tensor>().label() == reserved::overlap_label());
        if (kronecker != (pass == 0)) continue;
        for (const auto& i : get_used_indices(*f)) {
          name(i);
          for (const auto& p : i.proto_indices()) name(p);
        }
      }
    for (const auto& f : factors) {
      if (data_.projectors && is_projector(*f)) {
        t.projections.push_back(projection(f->as<AbstractTensor>()));
        continue;
      }
      t.factors.push_back(factor(f));
    }
    for (const auto& i : summed_indices_)
      t.domains.push_back(
          range_of(i.space(), quantum_numbers_, data_.occupied));
    t.after.resize(t.domains.size());
    for (std::size_t k = 0; k != t.factors.size(); ++k) {
      int last = -1;
      auto see = [&last](const Slot& s) {
        last = std::max(last, s.own);
        for (const auto p : s.protos) last = std::max(last, p);
      };
      const auto& g = t.factors[k];
      for (const auto* slots : {&g.bra, &g.ket, &g.aux})
        for (const auto& s : *slots) see(s);
      if (last < static_cast<int>(t.externals))
        t.first.push_back(k);
      else
        t.after[last - t.externals].push_back(k);
    }
    return t;
  }

  int name(const Index& i) {
    const std::wstring label(i.label());
    if (const auto it = externals_.find(label); it != externals_.end())
      return it->second;
    if (const auto it = summed_.find(label); it != summed_.end())
      return it->second;
    const auto n = static_cast<int>(externals_.size() + summed_.size());
    summed_.emplace(label, n);
    summed_indices_.push_back(i);
    return n;
  }

  Slot slot(const Index& i) {
    Slot result{.own = name(i), .symmetric = i.symmetric_proto_indices()};
    for (const auto& p : i.proto_indices()) result.protos.push_back(name(p));
    return result;
  }

  template <typename Range>
  container::svector<Slot> slots(Range&& range) {
    container::svector<Slot> result;
    for (const auto& i : range) result.push_back(slot(i));
    return result;
  }

  Projection projection(const AbstractTensor& t) {
    Projection result{.antisymmetric =
                          t._label() == reserved::antisymm_label()};
    auto external = [this](const Index& i) {
      const auto n = name(i);
      if (n >= static_cast<int>(externals_.size()))
        throw Exception(
            "the evaluator applies a projector over external "
            "labels only, not over " +
            toUtf8(i.full_label()));
      return n;
    };
    for (const auto& i : t._bra()) result.bra.push_back(external(i));
    for (const auto& i : t._ket()) result.ket.push_back(external(i));
    if (!result.antisymmetric && result.bra.size() != result.ket.size())
      throw Exception("the evaluator applies Ŝ of equal bra and ket ranks");
    return result;
  }

 public:
  Factor factor(const ExprPtr& f) {
    auto result = unlabeled_factor(f);
    result.labels = labels_of(result);
    return result;
  }

 private:
  Factor unlabeled_factor(const ExprPtr& f) {
    if (f->is<Power>()) {
      const auto& p = f->as<Power>();
      auto result = factor(p.base());
      result.exponent = p.exponent();
      result.power_conjugated = p.conjugated();
      return result;
    }
    if (f->is<TensorPower>()) {
      auto result = factor(f->as<TensorPower>().base());
      result.exponent = f->as<TensorPower>().exponent();
      return result;
    }
    if (f->is<Constant>()) return {.number = f->as<Constant>().value<Number>()};
    if (f->is<Variable>()) {
      const auto& x = f->as<Variable>();
      return {.number = number_of(hash_of(x.label()),
                                  data_.real.contains(std::wstring(x.label()))),
              .variable = true,
              .key = hash_of(x.label()),
              .conjugated = x.conjugated(),
              .real = data_.real.contains(std::wstring(x.label()))};
    }
    if (!f->is<AbstractTensor>() || !f->is_cnumber())
      throw Exception("the evaluator has no value of " + toUtf8(f->to_latex()));
    const auto& t = f->as<AbstractTensor>();
    const auto label = t._label();
    Factor result{.kind = Factor::Kind::data,
                  .symmetry = t._symmetry(),
                  .braket = t._braket_symmetry(),
                  .column = t._column_symmetry(),
                  .bra = slots(t._bra()),
                  .ket = slots(t._ket()),
                  .aux = slots(t._aux())};
    const auto base = strip_adjoint_label(label);
    result.key = hash_of(base);
    result.real = data_.real.contains(std::wstring(base));
    if (is_adjoint_label(label)) {
      std::swap(result.bra, result.ket);
      result.conjugated = true;
    }
    const bool kronecker = label == reserved::kronecker_label(),
               overlap = label == reserved::overlap_label();
    if ((kronecker || overlap) &&
        (result.bra.size() != 1 || result.ket.size() != 1 ||
         !result.aux.empty()))
      throw Exception(
          "the evaluator reads a δ or s of one bra and one ket "
          "slot, not " +
          toUtf8(f->to_latex()));
    if (kronecker) result.kind = Factor::Kind::kronecker;
    if (overlap) {
      result.kind = Factor::Kind::overlap;
      // the protoindices are compared by value, in is_identity_at()
      result.unit_overlap = is_kronecker_equivalent(
          t._bra()[0].drop_proto_indices(), t._ket()[0].drop_proto_indices());
    }
    return result;
  }
};

template <typename N>
N summed(const Term& t, std::size_t k, N partial, Assignment& v,
         const Perturbation* p) {
  if (k == t.domains.size()) return partial;
  N total{};
  for (const auto value : t.domains[k]) {
    v[t.externals + k] = value;
    auto z = partial;
    for (const auto f : t.after[k]) {
      z *= value_of<N>(t.factors[f], v, p);
      if (is_zero(z)) break;
    }
    if (!is_zero(z)) total += summed(t, k + 1, z, v, p);
  }
  return total;
}

/// the value of the term @p t at @p v , its external labels assigned, without
/// its projectors
template <typename N>
N unprojected_value_of(const Term& t, Assignment v, const Perturbation* p) {
  v.resize(t.externals + t.domains.size());
  N z = t.scalar;
  for (const auto f : t.first) z *= value_of<N>(t.factors[f], v, p);
  return is_zero(z) ? z : summed(t, 0, z, v, p);
}

/// the value of the term @p t at @p v with its projections from the
/// @p k -th on applied: the normalized sum of the values at the permuted
/// assignments
template <typename N>
N projected_value_of(const Term& t, std::size_t k, const Assignment& v,
                     const Perturbation* p) {
  if (k == t.projections.size()) return unprojected_value_of<N>(t, v, p);
  const auto& proj = t.projections[k];
  N total{};
  double weight = 0;
  auto permuted = v;
  if (proj.antisymmetric) {
    const auto bras = signed_permutations(proj.bra.size());
    const auto kets = signed_permutations(proj.ket.size());
    weight = 1. / static_cast<double>(bras.size() * kets.size());
    for (const auto& [by_bra, bra_sign] : bras)
      for (const auto& [by_ket, ket_sign] : kets) {
        for (std::size_t i = 0; i != proj.bra.size(); ++i)
          permuted[proj.bra[i]] = v[proj.bra[by_bra[i]]];
        for (std::size_t i = 0; i != proj.ket.size(); ++i)
          permuted[proj.ket[i]] = v[proj.ket[by_ket[i]]];
        auto z = projected_value_of<N>(t, k + 1, permuted, p);
        z *= N(static_cast<double>(bra_sign * ket_sign));
        total += z;
      }
  } else {
    const auto columns = signed_permutations(proj.bra.size());
    weight = 1. / static_cast<double>(columns.size());
    for (const auto& [by, sign] : columns) {
      for (std::size_t i = 0; i != proj.bra.size(); ++i) {
        permuted[proj.bra[i]] = v[proj.bra[by[i]]];
        permuted[proj.ket[i]] = v[proj.ket[by[i]]];
      }
      total += projected_value_of<N>(t, k + 1, permuted, p);
    }
  }
  total *= N(weight);
  return total;
}

/// the value of @p terms at @p externals , the values of the external labels
template <typename N = Number>
N value_of(const container::svector<Term>& terms,
           const container::svector<int>& externals,
           const Perturbation* p = nullptr) {
  N result{};
  const Assignment v(externals.begin(), externals.end());
  for (const auto& t : terms) result += projected_value_of<N>(t, 0, v, p);
  return result;
}

/// whether the indices of @p es and @p free differ in their quantum numbers:
/// if not, each space is read without them, which changes no comparison and
/// halves the values of a spin-orbital index
bool quantum_numbers_differ(std::initializer_list<ExprPtr> es,
                            const container::set<Index>& free) {
  container::set<std::int32_t> qns;
  auto add = [&qns](const Index& i) {
    qns.insert(i.space().qns().to_int32());
    for (const auto& p : i.proto_indices())
      qns.insert(p.space().qns().to_int32());
  };
  for (const auto& i : free) add(i);
  for (const auto& e : es)
    e->visit(
        [&add](const ExprPtr& leaf) {
          if (leaf->is<AbstractTensor>())
            for (const auto& i : leaf->as<AbstractTensor>()._slots()) add(i);
        },
        /*atoms_only=*/true);
  return qns.size() > 1;
}

}  // namespace

std::complex<double> element_of(const ExprPtr& e, const Values& values,
                                const Data& data) {
  container::svector<std::wstring> labels;
  container::svector<int> numbers;
  for (const auto& [label, value] : values) {
    labels.push_back(label);
    numbers.push_back(value);
  }
  container::set<Index> free;
  e->visit(
      [&](const ExprPtr& leaf) {
        if (!leaf->is<AbstractTensor>()) return;
        for (const auto& i : leaf->as<AbstractTensor>()._slots())
          if (values.contains(std::wstring(i.label()))) free.insert(i);
      },
      /*atoms_only=*/true);
  const auto qns = quantum_numbers_differ({e}, free);
  return value_of(Reader(labels, data, qns).terms(e), numbers);
}

namespace {

/// calls @p f with each element of @p a and of @p b over the labels of
/// @p free and their protoindices that @p sampling reads, and a description
/// of where it is, until @p f says stop
/// @return whether every element was read
bool for_each_element(
    const ExprPtr& a, const ExprPtr& b, const container::set<Index>& free,
    const Data& data, const Sampling& sampling,
    const std::function<Draw(Number, Number,
                             const std::function<std::string()>&)>& f) {
  const auto qns = quantum_numbers_differ({a, b}, free);
  Domains domains;
  for (const auto& i : free) {
    add_domain(i, {}, domains, qns, data.occupied);
    for (const auto& p : i.proto_indices())
      add_domain(p, {}, domains, qns, data.occupied);
  }
  container::svector<std::wstring> labels;
  for (const auto& [label, range] : domains) labels.push_back(label);
  const auto as = Reader(labels, data, qns).terms(a);
  const auto bs = Reader(labels, data, qns).terms(b);
  return for_each_sampled(
      domains, groups_of(domains, sampling.groups), sampling,
      [&](const container::svector<int>& numbers) {
        const auto where = [&]() {
          std::wstring at;
          for (std::size_t j = 0; j != labels.size(); ++j)
            at += (j ? L", " : L"") + labels[j] + L"=" +
                  std::to_wstring(numbers[j]);
          return toUtf8(at);
        };
        return f(value_of(as, numbers), value_of(bs, numbers), where);
      });
}

bool close(Number x, Number y) {
  return std::abs(x - y) <= 1e-12 * std::max(1., std::abs(x));
}

bool is_nonzero(Number x, Number y) {
  return std::abs(x) > 1e-12 || std::abs(y) > 1e-12;
}

std::string printed(Number x) {
  return "(" + std::to_string(x.real()) + ", " + std::to_string(x.imag()) + ")";
}

}  // namespace

Comparison compare_arrays(const ExprPtr& a, const ExprPtr& b,
                          const container::set<Index>& free, const Data& data,
                          const Sampling& sampling) {
  Comparison result;
  result.every_element = for_each_element(
      a, b, free, data, sampling, [&](Number x, Number y, const auto& where) {
        ++result.elements;
        const bool nonzero = is_nonzero(x, y);
        if (nonzero) ++result.nonzero;
        if (close(x, y)) return nonzero ? Draw::informative : Draw::empty;
        result.difference =
            "at " + where() + ": " + printed(x) + " against " + printed(y);
        return Draw::stop;
      });
  return result;
}

namespace {

/// the perturbation of the element of @p x at @p numbers , the values of
/// @p labels , its labels; empty if that element is no variable: it
/// vanishes (an antisymmetric one with two equal slots), or it is an overlap
/// within one basis, which is δ
std::optional<Perturbation> perturbation_of(
    const ExprPtr& x, const container::svector<std::wstring>& labels,
    const container::svector<int>& numbers, const Data& data) {
  Reader reader(labels, data, /*quantum_numbers=*/false);
  const auto f = reader.factor(x);
  Perturbation p{.key = f.key,
                 .variable = f.variable,
                 .conjugated = f.conjugated,
                 .real = f.real};
  if (f.variable) return p;
  const Assignment v(numbers.begin(), numbers.end());
  if (f.kind == Factor::Kind::overlap && is_identity_at(f, v))
    return std::nullopt;
  const auto c = canonical_of(f, v);
  if (c.vanishes) return std::nullopt;
  p.b = c.b;
  p.k = c.k;
  p.x = c.x;
  p.sign = c.sign;
  p.swapped = c.swapped;
  return p;
}

/// d/dz of @p terms at @p numbers , z the element @p p raises
Number tangent_at(const container::svector<Term>& terms,
                  const container::svector<int>& numbers,
                  const Perturbation& p) {
  // a term without a factor of the raised symbol adds nothing to the
  // derivative
  auto may_raise = [&p](const Factor& f) {
    if (f.key != p.key) return false;
    return p.variable ? f.kind == Factor::Kind::number && f.variable
                      : f.kind == Factor::Kind::data ||
                            f.kind == Factor::Kind::overlap;
  };
  Dual z;
  const Assignment v(numbers.begin(), numbers.end());
  for (const auto& t : terms)
    if (std::ranges::any_of(t.factors, may_raise))
      z += projected_value_of<Dual>(t, 0, v, &p);
  if (p.real) return z.by + z.by_conjugate;
  return p.conjugated ? z.by_conjugate : z.by;
}

/// the labels of @p domains in order
container::svector<std::wstring> labels_of(const Domains& domains) {
  container::svector<std::wstring> result;
  for (const auto& [label, range] : domains) result.push_back(label);
  return result;
}

/// an element of @p x that its symmetry makes of another, with its sign
/// relative to it: the values of the labels of x that read it, if any does
/// (a slot carries its protoindices along, which a label of x need not
/// follow)
struct Image {
  std::optional<container::svector<int>> numbers;
  int sign = 1;
};

/// the distinct arrangements of @p x , @p real or not, that its symmetry
/// makes of the one at @p numbers , values of the labels of @p domains , that
/// one first; none if it is equal to minus itself. An arrangement is an
/// element with the values of the protoindices that are no slots, in order:
/// two orders of those values are two arrangements of one element. Values
/// that read an element must lie in the domains of their labels: a block of
/// other spins is no element of x
container::svector<Image> orbit_of(const ExprPtr& x, const Domains& domains,
                                   const container::svector<int>& numbers,
                                   bool real) {
  const auto labels = labels_of(domains);
  if (!x->is<AbstractTensor>()) return {Image{numbers, 1}};
  const auto& v = x->as<AbstractTensor>();
  auto position = [&labels](const Index& i) {
    return std::ranges::find(labels, std::wstring(i.label())) - labels.begin();
  };
  // a slot's value and those of its protoindices, sorted if symmetric
  auto address = [&](const Index& i, const container::svector<int>& at) {
    container::svector<int> result{at[position(i)]};
    container::svector<int> protos;
    for (const auto& p : i.proto_indices()) protos.push_back(at[position(p)]);
    if (i.symmetric_proto_indices()) std::ranges::sort(protos);
    result.insert(result.end(), protos.begin(), protos.end());
    return result;
  };
  // the addresses of the slots, then the values of the protoindices that
  // are no slots
  using Arrangement = container::svector<container::svector<int>>;
  auto arrangement = [&](const Rearrangement& r,
                         const container::svector<int>& at) {
    Arrangement result;
    for (const auto& i : r.slots) result.push_back(address(i, at));
    container::svector<int> protos;
    for (const auto& p : r.protoindices) protos.push_back(at[position(p)]);
    result.push_back(std::move(protos));
    return result;
  };
  const auto arrangements = rearrangements_of(v, real);
  const auto& as_is = arrangements.front();
  container::map<Arrangement, int> signs;
  container::svector<std::pair<Arrangement, Image>> images;
  for (const auto& r : arrangements) {
    const auto element = arrangement(r, numbers);
    const auto [it, added] = signs.emplace(element, r.sign);
    if (!added) {
      if (it->second != r.sign) return {};
      continue;
    }
    auto image = numbers;
    for (std::size_t k = 0; k != r.slots.size(); ++k)
      image[position(as_is.slots[k])] = numbers[position(r.slots[k])];
    for (std::size_t k = 0; k != r.protoindices.size(); ++k)
      image[position(as_is.protoindices[k])] =
          numbers[position(r.protoindices[k])];
    const auto read = arrangement(as_is, image);
    bool in_domains = true;
    for (std::size_t k = 0; k != image.size(); ++k)
      in_domains =
          in_domains && std::ranges::find(domains[k].second, image[k]) !=
                            domains[k].second.end();
    images.emplace_back(
        element, Image{read == element && in_domains ? std::optional(image)
                                                     : std::nullopt,
                       r.sign});
  }
  container::svector<Image> result;
  for (auto& [element, image] : images) result.push_back(std::move(image));
  return result;
}

}  // namespace

Comparison compare_tangents(const ExprPtr& f, const ExprPtr& r,
                            const ExprPtr& x, const Data& data,
                            const Sampling& sampling) {
  auto applied = data;
  applied.projectors = true;
  const auto g = with_dummies_apart(f, x);
  const auto free = externals_of(f);
  const auto x_indices = indices_of_leaf(x);
  auto all = free;
  for (const auto& i : x_indices) all.insert(i);
  const auto qns = quantum_numbers_differ({g, r, x}, all);
  Domains f_domains, x_domains;
  for (const auto& i : free) {
    add_domain(i, {}, f_domains, qns, data.occupied);
    for (const auto& p : i.proto_indices())
      add_domain(p, {}, f_domains, qns, data.occupied);
  }
  for (const auto& i : x_indices)
    add_domain(i, {}, x_domains, qns, data.occupied);
  const auto f_labels = labels_of(f_domains), x_labels = labels_of(x_domains);
  auto r_labels = f_labels;
  r_labels.insert(r_labels.end(), x_labels.begin(), x_labels.end());
  const auto f_terms = Reader(f_labels, applied, qns).terms(g);
  const auto r_terms = Reader(r_labels, applied, qns).terms(r);

  auto domains = f_domains;
  domains.insert(domains.end(), x_domains.begin(), x_domains.end());
  container::svector<int> groups(f_domains.size(), 0);
  groups.resize(domains.size(), 1);
  Comparison result;
  result.every_element = for_each_sampled(
      domains, groups, sampling, [&](const container::svector<int>& all) {
        const container::svector<int> fv(all.begin(),
                                         all.begin() + f_labels.size()),
            xv(all.begin() + f_labels.size(), all.end());
        const auto p = perturbation_of(x, x_labels, xv, applied);
        if (!p) return Draw::empty;
        const auto tangent = tangent_at(f_terms, fv, *p);
        auto numbers = all;
        const auto at_J = value_of(r_terms, numbers);
        const auto orbit = orbit_of(x, x_domains, xv, p->real);
        // the operations of x: the rearrangements orbit_of() reads
        const auto operations = static_cast<double>(
            x->is<AbstractTensor>()
                ? rearrangements_of(x->as<AbstractTensor>(), p->real).size()
                : 1);
        // the result at J times the number of distinct partners of J equals
        // the number of operations of x times the derivative along J and its
        // partners
        const Number expected = tangent * operations;
        const Number found = static_cast<double>(orbit.size()) * at_J;
        std::string asymmetric;
        for (const auto& [image, sign] : orbit) {
          if (!image || !asymmetric.empty()) continue;
          std::copy(image->begin(), image->end(), numbers.begin() + fv.size());
          const auto at_image = value_of(r_terms, numbers);
          if (!close(static_cast<double>(sign) * at_J, at_image))
            asymmetric =
                "; the result lacks the symmetry of x: " + printed(at_J) +
                " against " + printed(at_image);
        }
        ++result.elements;
        const bool nonzero = is_nonzero(expected, found);
        if (nonzero) ++result.nonzero;
        if (close(expected, found) && asymmetric.empty())
          return nonzero ? Draw::informative : Draw::empty;
        std::wstring at;
        for (std::size_t j = 0; j != r_labels.size(); ++j)
          at +=
              (j ? L", " : L"") + r_labels[j] + L"=" + std::to_wstring(all[j]);
        result.difference = "at " + toUtf8(at) + ": tangent times " +
                            std::to_string(static_cast<int>(operations)) +
                            " operations " + printed(expected) +
                            " against the result times " +
                            std::to_string(orbit.size()) + " partners " +
                            printed(found) + asymmetric;
        return Draw::stop;
      });
  return result;
}

std::string text_of(const ExprPtr& e) {
  try {
    return toUtf8(serialize(e));
  } catch (const Exception&) {
    return toUtf8(e->to_latex());
  }
}

}  // namespace sequant::deriv::test
