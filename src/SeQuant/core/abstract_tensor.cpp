//
// Created by Eduard Valeyev on 2019-03-24.
//

#include "abstract_tensor.hpp"

namespace sequant {

TensorCanonicalizer::~TensorCanonicalizer() = default;

container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> &TensorCanonicalizer::instance_map_accessor() {
  static container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> map_;
  return map_;
}

container::vector<std::wstring>
&TensorCanonicalizer::cardinal_tensor_labels_accessor() {
  static container::vector<std::wstring> ctlabels_;
  return ctlabels_;
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance(
    std::wstring_view label) {
  auto &map = instance_map_accessor();
  // look for label-specific canonicalizer
  auto it = map.find(std::wstring{label});
  if (it != map.end()) {
    return it->second;
  } else {  // look for the default canonicalizer
    auto it = map.find(L"");
    if (it != map.end()) {
      return it->second;
    }
  }
  throw std::runtime_error("must first register canonicalizer via TensorCanonicalizer::register_instance(...)");
}

void TensorCanonicalizer::register_instance(std::shared_ptr<TensorCanonicalizer> can, std::wstring_view label) {
  auto &map = instance_map_accessor();
  map[std::wstring{label}] = can;
}


std::shared_ptr<Expr> DefaultTensorCanonicalizer::apply(AbstractTensor &t) {
  // tag all indices as ext->true/ind->false
  auto braket_view = braket(t);
  ranges::for_each(braket_view, [this](auto &idx) {
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
