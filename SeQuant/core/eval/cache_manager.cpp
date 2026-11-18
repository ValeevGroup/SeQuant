//
// Created by Bimal Gaudel on 9/26/23.
//

#include <SeQuant/core/eval/cache_manager.hpp>

#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/transform.hpp>

namespace sequant::eval {

// the state of the diagnostic meters (see cache_manager.hpp) lives here, so
// that it is process-wide rather than per shared library

std::size_t& AccessClock::counter() noexcept {
  static std::size_t c = 0;
  return c;
}

std::unordered_map<std::size_t, std::size_t>&
AccessClock::last_access_map() noexcept {
  static std::unordered_map<std::size_t, std::size_t> m;
  return m;
}

std::size_t& LookupMeter::calls() noexcept {
  static std::size_t c = 0;
  return c;
}

std::uint64_t& LookupMeter::nanos() noexcept {
  static std::uint64_t n = 0;
  return n;
}

std::function<void()>& PhaseTimer::fence_hook() noexcept {
  static std::function<void()> h;
  return h;
}

std::size_t& PhaseTimer::barrier_count() noexcept {
  static std::size_t c = 0;
  return c;
}

std::uint64_t& PhaseTimer::barrier_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::map<std::string, std::pair<std::uint64_t, std::uint64_t>>&
PhaseTimer::acc() noexcept {
  static auto* const m =
      new std::map<std::string, std::pair<std::uint64_t, std::uint64_t>>();
  return *m;
}

int& EvalImplTimeline::depth() noexcept {
  static thread_local int d = 0;
  return d;
}

std::uint64_t& EvalImplTimeline::body_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::uint64_t& EvalImplTimeline::gap_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::size_t& EvalImplTimeline::calls() noexcept {
  static std::size_t v = 0;
  return v;
}

std::chrono::steady_clock::time_point& EvalImplTimeline::last_exit() noexcept {
  static std::chrono::steady_clock::time_point t{};
  return t;
}

bool& EvalImplTimeline::have_last() noexcept {
  static bool b = false;
  return b;
}

std::uint64_t& EvalImplTimeline::max_gap_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::uint64_t& EvalImplTimeline::prod_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::uint64_t& EvalImplTimeline::phase_ns() noexcept {
  static std::uint64_t v = 0;
  return v;
}

std::unordered_map<std::size_t, std::size_t>&
BuildMeter::build_count() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
  return *m;
}

std::size_t& BuildMeter::total() noexcept {
  static std::size_t t = 0;
  return t;
}

std::unordered_map<std::size_t, std::string>& BuildMeter::label_of() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::string>();
  return *m;
}

std::unordered_map<std::size_t, std::size_t>&
BuildMeter::read_count() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
  return *m;
}

std::unordered_map<std::size_t, std::size_t>&
BuildMeter::max_life_of() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
  return *m;
}

std::size_t& DefUseMeter::op_clock() noexcept {
  static std::size_t c = 0;
  return c;
}

std::unordered_map<std::size_t, std::size_t>&
DefUseMeter::produce_at() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
  return *m;
}

std::unordered_map<std::size_t, std::size_t>&
DefUseMeter::store_count() noexcept {
  static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
  return *m;
}

std::size_t& DefUseMeter::reads() noexcept {
  static std::size_t r = 0;
  return r;
}

std::uint64_t& DefUseMeter::sum_dist() noexcept {
  static std::uint64_t s = 0;
  return s;
}

long double& DefUseMeter::sum_dist_bytes() noexcept {
  static long double s = 0;
  return s;
}

long double& DefUseMeter::sum_bytes() noexcept {
  static long double s = 0;
  return s;
}

std::size_t& DefUseMeter::max_dist() noexcept {
  static std::size_t m = 0;
  return m;
}

namespace detail {
std::size_t& peak_compose_max() noexcept {
  static std::size_t m = 0;
  return m;
}
}  // namespace detail

namespace {
// report the meters' totals at exit (if enabled), once per process
const LookupMeter::Reporter lookup_meter_reporter;
const PhaseTimer::Reporter phase_timer_reporter;
const EvalImplTimeline::Reporter eval_impl_timeline_reporter;
const BuildMeter::Reporter build_meter_reporter;
const DefUseMeter::Reporter def_use_meter_reporter;
}  // namespace

}  // namespace sequant::eval

namespace sequant {

namespace {

template <meta::eval_node Node, typename N, bool F>
void max_cache(Node const& node,        //
               CacheManager<N, F>& cm,  //
               AsyCost& curr,           //
               AsyCost& max) {
  if (auto ptr = cm.access(node); ptr) {
    if (cm.life(node) == 0) {
      curr -= Memory{}(node);
    }
    return;
  }
  if (!node.leaf()) {
    max_cache(node.left(), cm, curr, max);
    max_cache(node.right(), cm, curr, max);
    if (cm.exists(node)) {
      curr += Memory{}(node);
      max = std::max(curr, max);
      // simulate cache store
      auto s = cm.store_and_access(node, nullptr);
    }
  }
}

}  // namespace

AsyCost peak_cache(Sum const& expr, std::optional<size_t> min_repeats) {
  // Materialize into a vector so that nodes have stable addresses for
  // pointer-based scanning in cache_manager(). per-summand binarize for
  // cache-cost analysis only; the positional head doesn't escape.
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto const nodes_vec = expr | ranges::views::transform([](const auto& expr) {
                           return binarize<EvalExpr>(expr);
                         }) |
                         ranges::to_vector;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  auto cm = min_repeats ? cache_manager(nodes_vec, min_repeats.value())
                        : cache_manager(nodes_vec);
  auto max = AsyCost::zero();
  auto curr = AsyCost::zero();
  for (auto const& n : nodes_vec) max_cache(n, cm, curr, max);
  return max;
}

}  // namespace sequant
