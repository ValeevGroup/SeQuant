//
// Created by Eduard Valeyev on 2019-02-26.
//

#ifndef SEQUANT_RUNTIME_HPP
#define SEQUANT_RUNTIME_HPP

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <SeQuant/core/ranges.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#ifdef SEQUANT_HAS_EXECUTION_HEADER
#include <execution>
#else
#include <atomic>
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

/// thrown by the parallel primitives (sequant::parallel_do, sequant::for_each,
/// sequant::transform_reduce) when more than one invocation of the work
/// function throws; a single exception is rethrown as is
class ParallelExceptions : public Exception {
 public:
  /// the ordinal of a work item (the thread id for parallel_do) and the
  /// exception its invocation threw
  using Entry = std::pair<std::size_t, std::exception_ptr>;

  /// @param exceptions at least two entries, ordered by item ordinal
  explicit ParallelExceptions(std::vector<Entry> exceptions)
      : Exception(make_message(exceptions)),
        exceptions_(std::move(exceptions)) {}

  /// @return the exceptions, ordered by the ordinal of the item that threw
  /// each; std::rethrow_exception recovers their types
  const std::vector<Entry>& exceptions() const { return exceptions_; }

 private:
  std::vector<Entry> exceptions_;

  static std::string make_message(const std::vector<Entry>& exceptions) {
    std::string first;
    try {
      std::rethrow_exception(exceptions.front().second);
    } catch (const std::exception& e) {
      first = e.what();
    } catch (...) {
      first = "an exception not derived from std::exception";
    }
    return std::to_string(exceptions.size()) +
           " exceptions thrown by parallel work; the first, by item " +
           std::to_string(exceptions.front().first) + ": " + first;
  }
};

namespace detail {

/// collects the exceptions thrown by the invocations of the work function of
/// a parallel primitive, which must not escape them: an exception escaping a
/// parallel algorithm's element function or a std::thread calls
/// std::terminate
class ParallelExceptionCollector {
 public:
  /// invokes @p f, recording an exception it throws as that of item @p item
  template <typename F>
  void run(std::size_t item, F&& f) noexcept {
    try {
      std::forward<F>(f)();
    } catch (...) {
      add(item, std::current_exception());
    }
  }

  /// once the work is done: rethrows the recorded exception if there is one,
  /// throws ParallelExceptions if there are several, else does nothing
  void rethrow() {
    if (exceptions_.empty()) return;
    if (exceptions_.size() == 1)
      std::rethrow_exception(exceptions_.front().second);
    std::stable_sort(
        exceptions_.begin(), exceptions_.end(),
        [](const auto& e1, const auto& e2) { return e1.first < e2.first; });
    throw ParallelExceptions(std::move(exceptions_));
  }

 private:
  std::mutex mtx_;
  std::vector<ParallelExceptions::Entry> exceptions_;

  // the exceptions of a nested primitive are recorded individually, as
  // exceptions of the item whose invocation ran it
  void add(std::size_t item, std::exception_ptr exception) {
    std::vector<ParallelExceptions::Entry> entries;
    try {
      std::rethrow_exception(exception);
    } catch (const ParallelExceptions& nested) {
      for (const auto& [nested_item, nested_exception] : nested.exceptions())
        entries.emplace_back(item, nested_exception);
    } catch (...) {
      entries.emplace_back(item, std::move(exception));
    }
    std::scoped_lock lock(mtx_);
    exceptions_.insert(exceptions_.end(), entries.begin(), entries.end());
  }
};

}  // namespace detail

/// sets the number of threads to use for concurrent work
inline void set_num_threads(int nt) {
  SEQUANT_ENFORCE(nt >= 1, "set_num_threads(nthreads): invalid nthreads");
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
/// @note each spawned thread invokes its own copy of @p lambda, the calling
/// thread @p lambda itself; every instance sees the scoped contexts of the
/// calling thread
/// @throw the exception thrown by an instance, once all have finished;
/// ParallelExceptions if several throw
template <typename Lambda>
void parallel_do(Lambda&& lambda) {
  std::vector<std::thread> threads;
  const auto nthreads = num_threads();
  const auto overlays = detail::implicit_context_overlays();
  detail::ParallelExceptionCollector exceptions;
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(
          std::thread([&overlays, &exceptions, lambda, thread_id]() mutable {
            detail::ImplicitContextOverlaysScope scope(overlays);
            exceptions.run(thread_id, [&] { lambda(thread_id); });
          }));
    else
      exceptions.run(thread_id, [&] { lambda(thread_id); });
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
  exceptions.rethrow();
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
/// @throw the exception thrown by an invocation of @p op, once all have
/// finished; ParallelExceptions if several throw
/// @sa get_num_threads()
template <typename SizedRange, typename UnaryOp>
void for_each(SizedRange& rng, const UnaryOp& op) {
  using ranges::begin;
  using ranges::end;
  const auto overlays = detail::implicit_context_overlays();
  detail::ParallelExceptionCollector exceptions;
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  // the parallel STL's TBB backend partitions over tbb::blocked_range, which
  // rejects iterators whose difference is not a plain integer (e.g. range-v3
  // views), so the algorithm runs over a vector of the range's iterators
  std::vector<decltype(begin(rng))> iterators;
  iterators.reserve(static_cast<std::size_t>(ranges::size(rng)));
  for (auto it = begin(rng); it != end(rng); ++it) iterators.push_back(it);
  std::for_each(std::execution::par, iterators.begin(), iterators.end(),
                [&overlays, &op, &exceptions, &iterators](const auto& it) {
                  detail::ImplicitContextOverlaysScope scope(overlays);
                  exceptions.run(
                      static_cast<std::size_t>(&it - iterators.data()),
                      [&] { op(*it); });
                });
#else
  std::atomic<size_t> work = 0;
  auto task = [&work, &op, &rng, &overlays, &exceptions,
               ntasks = ranges::size(rng)]() {
    detail::ImplicitContextOverlaysScope scope(overlays);
    auto it = ranges::begin(rng);
    size_t prev_task_id = 0;
    size_t task_id = work.fetch_add(1);
    while (task_id < ntasks) {
      std::advance(it, task_id - prev_task_id);
      exceptions.run(task_id, [&] { op(*it); });
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
  exceptions.rethrow();
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
/// @note @p reduce must not throw
/// @throw the exception thrown by an invocation of @p map, once all have
/// finished; ParallelExceptions if several throw
/// @sa get_num_threads()
template <typename SizedRange, typename T, typename BinaryReductionOp,
          typename UnaryMapOp>
T transform_reduce(SizedRange&& rng, T init, const BinaryReductionOp& reduce,
                   const UnaryMapOp& map) {
  using ranges::begin;
  using ranges::end;
  const auto overlays = detail::implicit_context_overlays();
  detail::ParallelExceptionCollector exceptions;
#ifdef SEQUANT_HAS_EXECUTION_HEADER
  // see for_each() for why the algorithm runs over the range's iterators
  std::vector<decltype(begin(rng))> iterators;
  iterators.reserve(static_cast<std::size_t>(ranges::size(rng)));
  for (auto it = begin(rng); it != end(rng); ++it) iterators.push_back(it);
  // an item whose map throws contributes nothing
  using MaybeT = std::optional<T>;
  MaybeT result = std::transform_reduce(
      std::execution::par, iterators.begin(), iterators.end(),
      MaybeT(std::move(init)),
      [&reduce](MaybeT a, MaybeT b) -> MaybeT {
        if (!a) return b;
        if (!b) return a;
        return reduce(std::move(*a), std::move(*b));
      },
      [&overlays, &map, &exceptions, &iterators](const auto& it) -> MaybeT {
        detail::ImplicitContextOverlaysScope scope(overlays);
        MaybeT mapped;
        exceptions.run(static_cast<std::size_t>(&it - iterators.data()),
                       [&] { mapped = map(*it); });
        return mapped;
      });
  exceptions.rethrow();
  return std::move(*result);
#else
  std::atomic<size_t> work = 0;
  std::mutex mtx;
  T result = init;
  auto task = [&work, &map, &reduce, &rng, &mtx, &result, &overlays,
               &exceptions, ntasks = ranges::size(rng)]() {
    detail::ImplicitContextOverlaysScope scope(overlays);
    size_t task_id = work.fetch_add(1);
    while (task_id < ntasks) {
      exceptions.run(task_id, [&] {
        const auto& item = rng[task_id];
        auto mapped_item = map(item);
        {  // critical section
          std::scoped_lock<std::mutex> lock(mtx);
          result = reduce(result, mapped_item);
        }
      });
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

  exceptions.rethrow();
  return result;
#endif
}

void set_locale();

}  // namespace sequant

#endif  // SEQUANT_RUNTIME_HPP
