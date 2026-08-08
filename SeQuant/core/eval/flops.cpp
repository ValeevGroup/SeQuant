#include <SeQuant/core/eval/flops.hpp>

#include <memory>
#include <mutex>
#include <vector>

namespace sequant::eval {

namespace {

/// One accumulator per thread that has recorded anything. The registry OWNS
/// the bins, so a bin outlives the thread that wrote it (a thread_local
/// object would be destroyed at thread exit, losing its contribution and
/// leaving the registry with a dangling pointer).
std::mutex& registry_mutex() {
  static std::mutex m;
  return m;
}

std::vector<std::unique_ptr<FlopReport>>& registry() {
  static std::vector<std::unique_ptr<FlopReport>> bins;
  return bins;
}

FlopReport& local_bin() {
  thread_local FlopReport* bin = [] {
    auto owned = std::make_unique<FlopReport>();
    FlopReport* raw = owned.get();
    std::lock_guard<std::mutex> g(registry_mutex());
    registry().push_back(std::move(owned));
    return raw;
  }();
  return *bin;
}

}  // namespace

double FlopReport::flops() const noexcept {
  double s = 0;
  for (auto const& t : by_category) s += t.flops;
  return s;
}

double FlopReport::elements() const noexcept {
  double s = 0;
  for (auto const& t : by_category) s += t.elements;
  return s;
}

double FlopReport::ops() const noexcept {
  double s = 0;
  for (auto const& t : by_category) s += t.ops;
  return s;
}

void FlopCounter::record(FlopCategory category, double flops, double elements,
                         double ops) noexcept {
  auto& t = local_bin().by_category[static_cast<std::size_t>(category)];
  t.flops += flops;
  t.elements += elements;
  t.ops += ops;
}

void FlopCounter::record_kernels(double kernels, double m_sum, double n_sum,
                                 double k_sum) noexcept {
  auto& b = local_bin();
  b.kernels += kernels;
  b.m_sum += m_sum;
  b.n_sum += n_sum;
  b.k_sum += k_sum;
}

void FlopCounter::record_unsourced() noexcept { local_bin().unsourced_ops += 1; }

void FlopCounter::reset() {
  std::lock_guard<std::mutex> g(registry_mutex());
  for (auto& bin : registry()) *bin = FlopReport{};
}

FlopReport FlopCounter::read() {
  FlopReport out;
  std::lock_guard<std::mutex> g(registry_mutex());
  for (auto const& bin : registry()) {
    for (std::size_t c = 0; c < flop_category_count; ++c) {
      out.by_category[c].flops += bin->by_category[c].flops;
      out.by_category[c].elements += bin->by_category[c].elements;
      out.by_category[c].ops += bin->by_category[c].ops;
    }
    out.kernels += bin->kernels;
    out.m_sum += bin->m_sum;
    out.n_sum += bin->n_sum;
    out.k_sum += bin->k_sum;
    out.unsourced_ops += bin->unsourced_ops;
  }
  return out;
}

}  // namespace sequant::eval
