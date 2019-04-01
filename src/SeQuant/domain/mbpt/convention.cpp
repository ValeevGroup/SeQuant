//
// Created by Eduard Valeyev on 2019-04-01.
//

#include "convention.hpp"
#include "../../core/index.hpp"
#include "../../core/tensor.hpp"
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
        else if (space == IndexSpace::active_unoccupied) // PNO
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

/// @brief registers standard instances of IndexSpace objects
void register_standard_instances() {
  const bool do_not_throw = false;
  IndexSpace::register_instance(L"i", IndexSpace::active_occupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"m", IndexSpace::occupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"a", IndexSpace::active_unoccupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"e", IndexSpace::unoccupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"p", IndexSpace::all, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"⍺'", IndexSpace::other_unoccupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"⍺", IndexSpace::complete_unoccupied, IndexSpace::nullqns, do_not_throw);
  IndexSpace::register_instance(L"κ", IndexSpace::complete, IndexSpace::nullqns, do_not_throw);
}

/// @brief creates an IndexRegistry
void make_default_indexregistry() {
  auto idxreg = std::make_shared<IndexRegistry>();
  auto& idxreg_ref = *idxreg;

  register_index(idxreg_ref, Index{L"i"}, 100);
  register_index(idxreg_ref, Index{L"m"}, 110);
  register_index(idxreg_ref, Index{L"a"}, 1000);
  register_index(idxreg_ref, Index{L"e"}, 1000);
  register_index(idxreg_ref, Index{L"p"}, 1110);
  register_index(idxreg_ref, Index{L"⍺'"}, 3000);
  register_index(idxreg_ref, Index{L"⍺"}, 4000);
  register_index(idxreg_ref, Index{L"κ"}, 4110);
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

