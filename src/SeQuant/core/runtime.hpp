//
// Created by Eduard Valeyev on 2019-02-26.
//

#ifndef SEQUANT_RUNTIME_HPP
#define SEQUANT_RUNTIME_HPP

#include <mutex>
#include <thread>

namespace sequant {

namespace detail {
inline int& nthreads_accessor() {
  static int nthreads = std::thread::hardware_concurrency() > 0
                            ? std::thread::hardware_concurrency()
                            : 1;
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
/// @note by default use the value returned std::thread::hardware_concurrency() if available, otherwise 1
/// @sa set_num_threads()
inline int num_threads() {
  return detail::nthreads_accessor();
}

/// Fires off @c nthreads instances of lambda in parallel, each in its own thread (thus @c nthreads-1 std::thread objects are created),
/// where @c nthreads is the value returned by get_num_threads() .
/// @tparam Lambda a function type for which @c Lambda(int) is valid
/// @param lambda the function object to execute, each will be invoked as @c lambda(thread_id) where @c thread_id is an integer in
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

/// Fires off @c ntasks instances of lambda in parallel, with at most @c nthreads instances executing concurrently,
/// where @c nthreads is the value returned by get_num_threads() .
/// @tparam Lambda a function type for which @c Lambda(int) is valid
/// @param lambda the function object to execute, each will be invoked as @c lambda(task_id) where @c task_id is an integer in
///        @c [0,ntasks) . @c lambda(t1) will be commenced not after @c lambda(t2) if @c t1<t2 .
/// @node The load is balanced dynamically.
/// @sa get_num_threads()
template <typename Lambda>
void parallel_for_each(Lambda&& lambda, const size_t ntasks) {
  std::atomic<size_t> work = 0;
  auto task = [&work, &lambda, ntasks](int thread_id) {
    size_t task_id = work.fetch_add(1);
    while(task_id < ntasks) {
      std::forward<Lambda>(lambda)(task_id);
      task_id = work.fetch_add(1);
    }
  };

  const auto nthreads = num_threads();
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(task, thread_id));
    else
      task(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
}

}


#endif //SEQUANT_RUNTIME_HPP
