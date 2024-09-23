#include <SeQuant/domain/eval/eval_result.hpp>

namespace sequant {

EvalResult::id_t EvalResult::next_id() noexcept {
  static std::atomic<id_t> grand_type_id = 0;
  return ++grand_type_id;
}

bool EvalResult::has_value() const noexcept { return value_.has_value(); }

void log_ta_tensor_host_memory_use(madness::World& world,
                                   std::string_view label) {
#if defined(TA_TENSOR_MEM_PROFILE) && defined(SEQUANT_EVAL_TRACE)
  auto logger = Logger::instance();
  if (logger.log_level_eval < 3) return;
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
