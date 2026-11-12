#include <SeQuant/core/logger.hpp>

namespace sequant {

// N.B. defined out-of-line so that the Singleton state lives (only) in
// SeQuant-symb, not in every library/executable that uses the Logger

Logger& Logger::instance() { return Singleton<Logger>::instance(); }

Logger* Logger::instance_ptr() { return Singleton<Logger>::instance_ptr(); }

Logger& Logger::set_instance(int log_level) {
  return Singleton<Logger>::set_instance(log_level);
}

}  // namespace sequant
