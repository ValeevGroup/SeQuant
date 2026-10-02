//
// Created by Eduard Valeyev on 2019-02-26.
//

#ifndef SEQUANT_RUNTIME_HPP
#define SEQUANT_RUNTIME_HPP

#include <cstdlib>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#include <SeQuant/core/ranges.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/exception.hpp>

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
    int nthreads = nthreads_cstr ? string_to<unsigned int>(nthreads_cstr)
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
  if (nt < 1) throw Exception("set_num_threads(nthreads): invalid nthreads");
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
/// @note each instance is invoked on its own copy of @p lambda and sees the
/// scoped contexts of the calling thread
template <typename Lambda>
void parallel_do(Lambda&& lambda) {
  std::vector<std::thread> threads;
  const auto nthreads = num_threads();
  const auto overlays = detail::implicit_context_overlays();
  auto task = [&overlays, lambda](int thread_id) mutable {
    detail::ImplicitContextOverlaysScope scope(overlays);
    lambda(thread_id);
  };
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(task, thread_id));
    else
      task(thread_id);
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
/// @note each invocation of @p op sees the scoped contexts of the calling
/// thread
/// @sa get_num_threads()
template <typename SizedRange, typename UnaryOp>
void for_each(SizedRange& rng, const UnaryOp& op) {
  using ranges::begin;
  using ranges::end;
  const auto overlays = detail::implicit_context_overlays();
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  // the parallel STL's TBB backend partitions over tbb::blocked_range, which
  // rejects iterators whose difference is not a plain integer (e.g. range-v3
  // views), so the algorithm runs over a vector of the range's iterators
  std::vector<decltype(begin(rng))> iterators;
  iterators.reserve(static_cast<std::size_t>(ranges::size(rng)));
  for (auto it = begin(rng); it != end(rng); ++it) iterators.push_back(it);
  std::for_each(std::execution::par, iterators.begin(), iterators.end(),
                [&overlays, &op](auto it) {
                  detail::ImplicitContextOverlaysScope scope(overlays);
                  op(*it);
                });
#else
  std::atomic<size_t> work = 0;
  auto task = [&work, &op, &rng, &overlays, ntasks = ranges::size(rng)]() {
    detail::ImplicitContextOverlaysScope scope(overlays);
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
/// @note each invocation of @p map sees the scoped contexts of the calling
/// thread
/// @sa get_num_threads()
template <typename SizedRange, typename T, typename BinaryReductionOp,
          typename UnaryMapOp>
T transform_reduce(SizedRange&& rng, T init, const BinaryReductionOp& reduce,
                   const UnaryMapOp& map) {
  using ranges::begin;
  using ranges::end;
  const auto overlays = detail::implicit_context_overlays();
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  // see for_each() for why the algorithm runs over the range's iterators
  std::vector<decltype(begin(rng))> iterators;
  iterators.reserve(static_cast<std::size_t>(ranges::size(rng)));
  for (auto it = begin(rng); it != end(rng); ++it) iterators.push_back(it);
  return std::transform_reduce(
      std::execution::par, iterators.begin(), iterators.end(), init, reduce,
      [&overlays, &map](auto it) {
        detail::ImplicitContextOverlaysScope scope(overlays);
        return map(*it);
      });
#else
  std::atomic<size_t> work = 0;
  std::mutex mtx;
  T result = init;
  auto task = [&work, &map, &reduce, &rng, &mtx, &result, &overlays,
               ntasks = ranges::size(rng)]() {
    detail::ImplicitContextOverlaysScope scope(overlays);
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

void set_locale();

}  // namespace sequant

#endif  // SEQUANT_RUNTIME_HPP
