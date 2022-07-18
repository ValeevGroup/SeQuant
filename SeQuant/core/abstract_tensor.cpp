//
// Created by Eduard Valeyev on 2019-03-24.
//

#include "abstract_tensor.hpp"

namespace sequant {

TensorCanonicalizer::~TensorCanonicalizer() = default;

std::pair<container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>*,
          std::unique_lock<std::recursive_mutex>>
TensorCanonicalizer::instance_map_accessor() {
  static container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
      map_;
  static std::recursive_mutex mtx_;
  return std::make_pair(&map_, std::unique_lock<std::recursive_mutex>{mtx_});
}

container::vector<std::wstring>&
TensorCanonicalizer::cardinal_tensor_labels_accessor() {
  static container::vector<std::wstring> ctlabels_;
  return ctlabels_;
}

std::shared_ptr<TensorCanonicalizer>
TensorCanonicalizer::nondefault_instance_ptr(std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  // look for label-specific canonicalizer
  auto it = map_ptr->find(std::wstring{label});
  if (it != map_ptr->end()) {
    return it->second;
  } else
    return {};
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance_ptr(
    std::wstring_view label) {
  auto result = nondefault_instance_ptr(label);
  if (!result)  // not found? look for default
    result = nondefault_instance_ptr(L"");
  return result;
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance(
    std::wstring_view label) {
  auto inst_ptr = instance_ptr(label);
  if (!inst_ptr)
    throw std::runtime_error(
        "must first register canonicalizer via "
        "TensorCanonicalizer::register_instance(...)");
  return inst_ptr;
}

void TensorCanonicalizer::register_instance(
    std::shared_ptr<TensorCanonicalizer> can, std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  (*map_ptr)[std::wstring{label}] = can;
}

bool TensorCanonicalizer::try_register_instance(
    std::shared_ptr<TensorCanonicalizer> can, std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  if (!map_ptr->contains(std::wstring{label})) {
    (*map_ptr)[std::wstring{label}] = can;
    return true;
  } else
    return false;
}

void TensorCanonicalizer::deregister_instance(std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  auto it = map_ptr->find(std::wstring{label});
  if (it != map_ptr->end()) {
    map_ptr->erase(it);
  }
}

ExprPtr DefaultTensorCanonicalizer::apply(AbstractTensor& t) {
  // tag all indices as ext->true/ind->false
  auto braket_view = braket(t);
  ranges::for_each(braket_view, [this](auto& idx) {
    auto it = external_indices_.find(std::wstring(idx.label()));
    auto is_ext = it != external_indices_.end();
    idx.tag().assign(
        is_ext ? 0 : 1);  // ext -> 0, int -> 1, so ext will come before
  });

  auto result = this->apply(t, std::less<Index>{});
  reset_tags(t);

  return result;
}

}  // namespace sequant
