//
// Created by Bimal Gaudel on 9/26/23.
//

#include <SeQuant/domain/eval/cache_manager.hpp>

#include <SeQuant/core/eval_node.hpp>

namespace sequant {

CacheManager::entry::entry(size_t count) noexcept
    : max_life{count},  //
      life_c{count},    //
      data_p{nullptr} {}

ERPtr CacheManager::entry::access() noexcept {
  if (!data_p) return nullptr;

  return decay() == 0 ? std::move(data_p) : data_p;
}

void CacheManager::entry::store(ERPtr data) noexcept { data_p = data; }

void CacheManager::entry::reset() noexcept {
  life_c = max_life;
  data_p = nullptr;
}

[[nodiscard]] size_t CacheManager::entry::life_count() const noexcept {
  return life_c;
}

[[nodiscard]] size_t CacheManager::entry::max_life_count() const noexcept {
  return max_life;
}

[[nodiscard]] int CacheManager::entry::decay() noexcept {
  return life_c > 0 ? static_cast<int>(--life_c) : 0;
}

ERPtr CacheManager::store(CacheManager::entry& entry, ERPtr data) noexcept {
  entry.store(data);
  return entry.access();
}

void CacheManager::reset() noexcept {
  for (auto&& [k, v] : cache_map_) v.reset();
}

ERPtr CacheManager::access(key_type key) noexcept {
  if (auto&& found = cache_map_.find(key); found != cache_map_.end())
    return found->second.access();

  return nullptr;
}

ERPtr CacheManager::store(key_type key, ERPtr data) noexcept {
  if (auto&& found = cache_map_.find(key); found != cache_map_.end())
    return store(found->second, data);
  return data;
}

bool CacheManager::exists(key_type key) const noexcept {
  return cache_map_.find(key) != cache_map_.end();
}

int CacheManager::life(key_type key) const noexcept {
  auto iter = cache_map_.find(key);
  auto end = cache_map_.end();
  return iter == end ? -1 : static_cast<int>(iter->second.life_count());
}

int CacheManager::max_life(key_type key) const noexcept {
  auto iter = cache_map_.find(key);
  auto end = cache_map_.end();
  return iter == end ? -1 : static_cast<int>(iter->second.max_life_count());
}

container::set<size_t> CacheManager::keys() const noexcept {
  return ranges::views::keys(cache_map_) | ranges::to<container::set<size_t>>;
}

CacheManager CacheManager::empty() noexcept { return CacheManager{{}}; }

template <typename NodeT,
          typename = std::enable_if_t<meta::is_eval_node<NodeT>>>
void max_cache(NodeT const& node,  //
               CacheManager& cm,   //
               AsyCost& curr,      //
               AsyCost& max) {
  auto const k = hash::value(*node);
  if (auto ptr = cm.access(k); ptr) {
    // std::cout << "[ACCESS][" << k << "]\n";
    if (cm.life(k) == 0) {
      curr -= Memory{}(node);
      // std::cout << "[RELEASE][" << k << "]\n";
    }
    return;
  }
  if (!node.leaf()) {
    max_cache(node.left(), cm, curr, max);
    max_cache(node.right(), cm, curr, max);
    if (cm.exists(k)) {
      curr += Memory{}(node);
      max = std::max(curr, max);
      // simulate cache store
      auto s = cm.store(k, eval_result<EvalScalar<double>>(0));
      //      std::cout << "[STORE][" << k << "]\n";
      //      std::cout << "[ACCESS][" << k << "]\n";
    }
  }
}

AsyCost peak_cache(Sum const& expr) {
  auto const nodes = expr | ranges::views::transform([](const auto& expr) {
                       return binarize<EvalExpr>(expr);
                     });
  auto cm = cache_manager(nodes);
  auto max = AsyCost::zero();
  auto curr = AsyCost::zero();
  for (auto const& n : nodes) max_cache(n, cm, curr, max);
  return max;
}

}  // namespace sequant
