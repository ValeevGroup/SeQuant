//
// Created by Eduard Valeyev on 2019-04-01.
//

#include "convention.hpp"
#include "SeQuant/core/index.hpp"
#include "SeQuant/core/tensor.hpp"
#include "op.hpp"

namespace sequant {
namespace mbpt {

namespace {
void register_index(IndexRegistry& reg, const Index& idx, long size) {
  IndexSpace space = idx.space();
  reg.make(idx, [size,space](const Index& index) -> long {
    assert(index.space() == space);
    const auto& proto_indices = index.proto_indices();
    if (proto_indices.empty())
      return size;
    else {
      if (proto_indices.size() == 1) {
        assert(proto_indices[0].space().type() == IndexSpace::active_occupied);
        if (space == IndexSpace::active_occupied)    // 1-index-specific occupied
          return 20;
        else if (space == IndexSpace::active_unoccupied) // OSV or PAO?
          return 500;
        else if (space == IndexSpace::all)
          return 520;
      }
      else if (proto_indices.size() == 2) {
        assert(proto_indices[1].space().type() == IndexSpace::active_occupied);
        if (space == IndexSpace::active_occupied)
          return 40;
        else if (space == IndexSpace::occupied)
          return 50;
        else if (space == IndexSpace::active_unoccupied) // CVS, e.g. PNO
          return 50;
        else if (space == IndexSpace::all)
          return 100;
      }
      abort();  // what to return here by default?
    }
  });
}
}  // anonymous namespace

namespace qcifs {

namespace {
enum class qns { none = 0, alpha = 1, beta = 2 };
auto qndecorate(qns qn, std::wstring_view label) {
  switch (static_cast<int>(qn)) {
    case 0: return std::wstring(label);
    case 1: return std::wstring(label) + L"⁺";
    case 2: return std::wstring(label) + L"⁻";
    default: assert(false && "invalid quantum number");
  }
};
};

/// @brief registers standard instances of IndexSpace objects
void register_standard_instances() {
  const bool do_not_throw = true;
  for(int s=0; s<=2; ++s) {
    auto qnattr = s==0 ? IndexSpace::nullqns : (s==1 ? IndexSpace::alpha : IndexSpace::beta);
    auto declab = [s](auto&& label) {
      return qndecorate(static_cast<qns>(s), label);
    };
    IndexSpace::register_instance(declab(L"i"), IndexSpace::active_occupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"m"), IndexSpace::occupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"a"), IndexSpace::active_unoccupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"e"), IndexSpace::unoccupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"p"), IndexSpace::all, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"⍺'"), IndexSpace::other_unoccupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"⍺"), IndexSpace::complete_unoccupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"κ"), IndexSpace::complete, qnattr, do_not_throw);
  }
}

/// @brief creates an IndexRegistry
void make_default_indexregistry() {
  auto idxreg = std::make_shared<IndexRegistry>();
  auto& idxreg_ref = *idxreg;

  for(int s=0; s<=2; ++s) {
    auto declab = [s](auto &&label) {
      return qndecorate(static_cast<qns>(s), label);
    };
    register_index(idxreg_ref, Index{declab(L"i")}, 100);
    register_index(idxreg_ref, Index{declab(L"m")}, 110);
    register_index(idxreg_ref, Index{declab(L"a")}, 1000);
    register_index(idxreg_ref, Index{declab(L"e")}, 1000);
    register_index(idxreg_ref, Index{declab(L"p")}, 1110);
    register_index(idxreg_ref, Index{declab(L"⍺'")}, 3000);
    register_index(idxreg_ref, Index{declab(L"⍺")}, 4000);
    register_index(idxreg_ref, Index{declab(L"κ")}, 4110);
  }
}

}  // namespace qcifs

/// Loads defaults for Convention @c conv
void set_default_convention(Convention conv) {
  switch (conv) {
    case Convention::QCiFS: {
      using namespace qcifs;
      register_standard_instances();
      make_default_indexregistry();
      TensorCanonicalizer::set_cardinal_tensor_labels(mbpt::cardinal_tensor_labels);
    }
  }
}

}
}

