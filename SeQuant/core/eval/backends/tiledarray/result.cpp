#include <SeQuant/core/eval/backends/tiledarray/result.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/backends/tiledarray/tiledarray.hpp>
#include <SeQuant/core/logger.hpp>

#include <atomic>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace sequant {

namespace detail {

std::atomic<std::size_t>& fence_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
std::atomic<std::size_t>& wait_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
std::atomic<std::size_t>& slice_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
std::atomic<long long>& slice_ns() {
  static std::atomic<long long> c{0};
  return c;
}

namespace {
struct FenceReporter {
  ~FenceReporter() {
    if (std::getenv("SEQUANT_SYNC_STATS"))
      std::cerr << "TOTAL gop.fence() calls = " << fence_counter().load()
                << " ; wait_for_lazy_cleanup calls = " << wait_counter().load()
                << " ; slice_mode calls = " << slice_counter().load()
                << " ; slice_mode total = " << (slice_ns().load() / 1e9) << " s"
                << "\n";
  }
};
const FenceReporter fence_reporter;
}  // namespace

}  // namespace detail

void log_ta_tensor_host_memory_use([[maybe_unused]] madness::World& world,
                                   [[maybe_unused]] std::string_view label) {
#if defined(SEQUANT_EVAL_TRACE) && defined(TA_TENSOR_MEM_PROFILE)
  auto logger = Logger::instance();
  if (logger.eval.level < 3) return;
  std::vector<std::uint64_t> hwsize(world.size(), 0);
  std::vector<std::uint64_t> currsize(world.size(), 0);
  std::vector<std::uint64_t> actsize(world.size(), 0);
  hwsize[world.rank()] =
      TA::hostEnv::instance()->host_allocator_getActualHighWatermark();
  currsize[world.rank()] =
      TA::hostEnv::instance()->host_allocator().getCurrentSize();
  actsize[world.rank()] =
      TA::hostEnv::instance()->host_allocator().getActualSize();
  world.gop.sum(hwsize.data(), hwsize.size());
  world.gop.sum(currsize.data(), currsize.size());
  world.gop.sum(actsize.data(), actsize.size());

  std::ostringstream oss;
  oss << label << ": TA_TENSOR_MEM_PROFILE allocation statistics (MiB):\n";
  oss << std::setw(5) << "rank"  //
      << std::setw(12) << "hw"   //
      << std::setw(12) << "cur"  //
      << std::setw(12) << "act"  //
      << '\n';                   //
  oss << "--------------------------------------------\n";
  std::uint64_t total = 0;
  for (auto rank = 0; rank != world.size(); ++rank) {
    oss << std::setw(5) << rank                         //
        << std::setw(12) << hwsize[rank] / (1 << 20)    //
        << std::setw(12) << currsize[rank] / (1 << 20)  //
        << std::setw(12) << actsize[rank] / (1 << 20)   //
        << '\n';
    total += currsize[rank] / (1 << 20);
  }
  oss << std::setw(5) << "total"  //
      << std::setw(12) << ""      //
      << std::setw(12) << total   //
      << std::setw(12) << ""      //
      << '\n';
  oss << "--------------------------------------------" << std::endl;
  write_log(logger, oss.str());
#endif
}

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY
