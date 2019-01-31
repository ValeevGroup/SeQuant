//
// Created by Eduard Valeyev on 2019-01-30.
//

#include "tensor.hpp"

namespace sequant2 {

Tensor::~Tensor() = default;
TensorCanonicalizer::~TensorCanonicalizer() = default;

std::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> &TensorCanonicalizer::instance_map_accessor() {
  static std::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> map_;
  return map_;
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance(std::wstring_view label) {
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

std::shared_ptr<Expr> DefaultTensorCanonicalizer::apply(Tensor &t) {
  auto symmetry = t.symmetry();
  auto is_antisymm = symmetry == Symmetry::antisymm;

  // can only handle anisymmetric case so far
  assert(is_antisymm);

  // tag all indices as ext->true/ind->false
  auto braket_view = this->braket(t);
  ranges::for_each(braket_view, [this](auto &idx) {
    auto it = external_indices_.find(std::wstring(idx.label()));
    auto is_ext = it != external_indices_.end();
    if (is_ext) {
      idx.tag(is_ext);
    }
  });

  //
  auto sort_swappables = [](const auto &begin, const auto &end, const auto &compare) {
    const auto len = end - begin;
    switch (len) {
      case 0: return;

      case 1: return;

      case 2: {
        auto &elem1 = *begin;
        auto &elem2 = *(begin + 1);
        if (compare(elem1, elem2))
          return;
        else {
          swap(elem1, elem2);
          return;
        }
      }
        break;

      case 3: {
        // bubble sort
        // {2,3} -> [2,3]
        {
          auto &elem2 = *(begin + 1);
          auto &elem3 = *(begin + 2);
          const auto lt23 = compare(elem2, elem3);
          if (!lt23)
            swap(elem2, elem3);
        }
        // sort {1,[2,3]} -> [1, 2, 3]
        {
          auto &elem1 = *(begin);
          auto &elem2 = *(begin + 1);
          auto &elem3 = *(begin + 2);
          const auto lt12 = compare(elem1, elem2);
          const auto lt13 = compare(elem1, elem3);
          if (!lt12)
            swap(elem1, elem2);
          if (!lt13)
            swap(elem2, elem3);
        }
        return;
      }
        break;

      default:abort();  // not yet implemented
    }
  };

  bool even = true;
  switch (symmetry) {
    case Symmetry::antisymm: {
      auto comp = [](const auto &idx1, const auto &idx2) {
        const auto idx1_is_ext = idx1.has_tag();
        const auto idx2_is_ext = idx2.has_tag();
        if (idx1_is_ext == idx2_is_ext)
          return idx1 < idx2;
        else if (idx1_is_ext) {
          return true;
        } else {
          return false;
        }
        abort();
      };
      auto &_bra = this->bra(t);
      auto &_ket = this->ket(t);
      using std::begin;
      using std::end;
//      std::wcout << "canonicalizing " << to_latex(t);
      IndexSwapper::thread_instance().reset();
      // std::{stable_}sort does not necessarily use swap! so much implement sort outselves .. thankfully ranks will be low so can stick with bubble
      sort_swappables(begin(_bra), end(_bra), comp);
      sort_swappables(begin(_ket), end(_ket), comp);
      even = IndexSwapper::thread_instance().even_num_of_swaps();
//      std::wcout << " is " << (even ? "even" : "odd") << " and produces " << to_latex(t) << std::endl;
    }
      break;

    case Symmetry::symm: {

    }
      break;

    case Symmetry::nonsymm: {

    }
      break;

    default:abort();
  }

  std::shared_ptr<Expr> result = is_antisymm ? (even == false ? make<Constant>(-1) : nullptr) : nullptr;
  return result;
}

std::shared_ptr<Expr> Tensor::canonicalize() {
  const auto &canonicalizer = TensorCanonicalizer::instance(label_);
  return canonicalizer->apply(*this);
}

}  // namespace sequant2