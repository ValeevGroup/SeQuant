//
// Created by Eduard Valeyev on 2019-02-26.
//

#ifndef SEQUANT_RUNTIME_HPP
#define SEQUANT_RUNTIME_HPP

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include <SeQuant/core/ranges.hpp>

#ifdef SEQUANT_HAS_EXECUTION_HEADER
#include <execution>
#else
#include <atomic>
#include <mutex>
#endif

namespace sequant {

namespace detail {
inline int& nthreads_accessor() {
  auto init_nthreads = []() {
    const auto nthreads_cstr = std::getenv("SEQUANT_NUM_THREADS");
    int nthreads = nthreads_cstr ? std::atoi(nthreads_cstr)
                                 : (std::thread::hardware_concurrency() > 0
                                        ? std::thread::hardware_concurrency()
                                        : 1);
    return nthreads;
  };
  static int nthreads = init_nthreads();
  return nthreads;
}
}  // namespace detail

/// sets the number of threads to use for concurrent work
inline void set_num_threads(int nt) {
  if (nt < 1)
    throw std::invalid_argument("set_num_threads(nthreads): invalid nthreads");
  detail::nthreads_accessor() = nt;
}

/// @return the number of threads to use for concurrent work
/// @note by default use the value returned std::thread::hardware_concurrency()
/// if available, otherwise 1
/// @sa set_num_threads()
inline int num_threads() { return detail::nthreads_accessor(); }

/// Fires off @c nthreads instances of lambda in parallel, each in its own
/// thread (thus @c nthreads-1 std::thread objects are created), where @c
/// nthreads is the value returned by get_num_threads() .
/// @tparam Lambda a function type for which @c Lambda(int) is valid
/// @param lambda the function object to execute, each will be invoked as @c
/// lambda(thread_id) where @c thread_id is an integer in
///        @c [0,nthreads) .
/// @sa get_num_threads()
template <typename Lambda>
void parallel_do(Lambda&& lambda) {
  std::vector<std::thread> threads;
  const auto nthreads = num_threads();
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(std::forward<Lambda>(lambda), thread_id));
    else
      std::forward<Lambda>(lambda)(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
}

/// Parallel version of std::for_each , using either parallel C++ algorithms
/// or manual threaded implementation with at most
/// @c nthreads instances executing concurrently, where @c nthreads is the value
/// returned by get_num_threads() .
/// @tparam SizedRange a sied range
/// @tparam UnaryOp a function type for which @c Lambda(int) is valid
/// @param rng the \p SizedRange object
/// @param op the function object to execute, each will be invoked as
/// @c op(std::advance(begin(rng) + task_id)) where @c task_id is an integer in
///        @c [0,size(rng)) . @c op(t1) will be commenced not
/// after @c op(t2) if @c t1<t2 .
/// @note The load is balanced dynamically.
/// @sa get_num_threads()
template <typename SizedRange, typename UnaryOp>
void for_each(SizedRange& rng, const UnaryOp& op) {
  using ranges::begin;
  using ranges::end;
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  std::for_each(std::execution::par_unseq, begin(rng), end(rng), op);
#else
  std::atomic<size_t> work = 0;
  auto task = [&work, &op, &rng, ntasks = ranges::size(rng)]() {
    auto it = ranges::begin(rng);
    size_t prev_task_id = 0;
    size_t task_id = work.fetch_add(1);
    while (task_id < ntasks) {
      std::advance(it, task_id - prev_task_id);
      op(*it);
      prev_task_id = task_id;
      task_id = work.fetch_add(1);
    }
  };

  const auto nthreads = num_threads();
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(task));
    else
      task();
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}

/// Does map+reduce (i.e., std::transform_reduce) on a range
/// using up to get_num_threads() threads.
/// @tparam SizedRange a sized range
/// @tparam BinaryReductionOp a function type such that
/// `reduce(identity,map(*begin(rng)))`, where `reduce` and `identity` are
/// objects of type `ReduceLambda` and `Identity`, respectively, is valid
/// @tparam UnaryMapOp a function type such that `map(*begin(rng))`, where `map`
/// is an object of type `MapLambda`, is valid
/// @tparam T a result type of \p ReduceLambda
/// @param rng the \p Range object
/// @param init the initial value for reduction
/// @param reduce the \p ReduceLambda object
/// @param map the \p MapLambda object
/// @sa get_num_threads()
template <typename SizedRange, typename T, typename BinaryReductionOp,
          typename UnaryMapOp>
T transform_reduce(SizedRange&& rng, T init, const BinaryReductionOp& reduce,
                   const UnaryMapOp& map) {
  using ranges::begin;
  using ranges::end;
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  return std::transform_reduce(std::execution::par_unseq, begin(rng), end(rng),
                               init, reduce, map);
#else
  std::atomic<size_t> work = 0;
  std::mutex mtx;
  T result = init;
  auto task = [&work, &map, &reduce, &rng, &mtx, &result,
               ntasks = ranges::size(rng)]() {
    size_t task_id = work.fetch_add(1);
    while (task_id < ntasks) {
      const auto& item = rng[task_id];
      auto mapped_item = map(item);
      {  // critical section
        std::scoped_lock<std::mutex> lock(mtx);
        result = reduce(result, mapped_item);
      }
      task_id = work.fetch_add(1);
    }
  };

  const auto nthreads = num_threads();
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(task));
    else
      task();
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();

  return result;
#endif
}

/// Set up the global locale settings
void set_locale();

/// @brief disables Thousands separators in the current locale
/// @note Need this in case of JSON/CSV output
void disable_thousands_separator();

}  // namespace sequant

#endif  // SEQUANT_RUNTIME_HPP
