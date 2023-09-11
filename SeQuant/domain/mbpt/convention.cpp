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
  reg.make(idx, [size, space](const Index& index) -> long {
    assert(index.space() == space);
    const auto& proto_indices = index.proto_indices();
    if (proto_indices.empty())
      return size;
    else {
      if (proto_indices.size() == 1) {
        assert(proto_indices[0].space().type() == IndexSpace::active_occupied);
        if (space == IndexSpace::active_occupied)  // 1-index-specific occupied
          return 20;
        else if (space == IndexSpace::active_unoccupied)  // OSV or PAO?
          return 500;
        else if (space == IndexSpace::all)
          return 520;
      } else if (proto_indices.size() == 2) {
        assert(proto_indices[1].space().type() == IndexSpace::active_occupied);
        if (space == IndexSpace::active_occupied)
          return 40;
        else if (space == IndexSpace::occupied)
          return 50;
        else if (space == IndexSpace::active_unoccupied)  // CVS, e.g. PNO
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
    case 0:
      return std::wstring(label);
    case 1:
      return std::wstring(label) + L"↑";
    case 2:
      return std::wstring(label) + L"↓";
    default:
      assert(false && "invalid quantum number");
  }
  abort();  // unreachable
};
};  // namespace

/// @brief registers "standard" instances of IndexSpace objects
void register_standard_instances() {
  const bool do_not_throw = true;
  for (int s = 0; s <= 2; ++s) {
    auto qnattr = s == 0 ? IndexSpace::nullqns
                         : (s == 1 ? IndexSpace::alpha : IndexSpace::beta);
    auto declab = [s](auto&& label) {
      return qndecorate(static_cast<qns>(s), label);
    };
    // these spaces are used in Fock space methods
    // based on single-determinant references
    // p,q,r... for OBS spstates introduced in DOI 10.1063/1.444231 (QCiFS I)
    IndexSpace::register_instance(declab(L"p"), IndexSpace::all, qnattr,
                                  do_not_throw);
    // {i,j,k.../a,b,c...} for active {occupied/unoccupied} spstates introduced
    // in DOI 10.1063/1.446736 (QCiFS III)
    IndexSpace::register_instance(declab(L"i"), IndexSpace::active_occupied,
                                  qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"a"), IndexSpace::active_unoccupied,
                                  qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"g"), IndexSpace::inactive_unoccupied,
                                 qnattr, do_not_throw);
    // {α,β.../κ,𝛌...} for complete {unoccupied/any} spstates introduced in
    // DOI 10.1063/1.459921 (MP2-R12 I)
    IndexSpace::register_instance(declab(L"α"), IndexSpace::complete_unoccupied,
                                  qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"κ"), IndexSpace::complete, qnattr,
                                  do_not_throw);
    // for orthogonal complement to p introduced in
    // DOI 10.1016/j.cplett.2004.07.061 (CABS)
    IndexSpace::register_instance(declab(L"α'"), IndexSpace::other_unoccupied,
                                  qnattr, do_not_throw);
    // m,n... for all occupied (including inactive/frozen orbitals) de facto
    // introduced in [DOI 10.1016/j.cplett.2004.07.061
    // (CABS)](https://dx.doi.org/10.1016/j.cplett.2004.07.061), though formally
    // not explicitly defined so
    IndexSpace::register_instance(declab(L"m"), IndexSpace::occupied, qnattr,
                                  do_not_throw);
    // introduced in MPQC LCAOWavefunction
    IndexSpace::register_instance(declab(L"e"), IndexSpace::unoccupied, qnattr,
                                  do_not_throw);
    // introduced in MPQC for GF, CT-F12, and other ad hoc uses
    IndexSpace::register_instance(declab(L"x"), IndexSpace::all_active, qnattr,
                                  do_not_throw);
    // used in virtual projection CT-F12 theory methods in MPQC.
    IndexSpace::register_instance(declab(L"c"), IndexSpace::complete_inactive_unoccupied, qnattr,
                                  do_not_throw);
    // e.g. see DOI 10.1063/5.0067511
    IndexSpace::register_instance(declab(L"u"), IndexSpace::active, qnattr,
                                  do_not_throw);
    // DOI 10.1063/5.0067511 uses I,J,K... and A,B,C... for these
    // although QCiFS uses capital letters for spin-free indices, using I/A this
    // way seems preferable
    IndexSpace::register_instance(
        declab(L"I"), IndexSpace::active_maybe_occupied, qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"A"),
                                  IndexSpace::active_maybe_unoccupied, qnattr,
                                  do_not_throw);
    // introduced here
    IndexSpace::register_instance(declab(L"M"), IndexSpace::maybe_occupied,
                                  qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"E"), IndexSpace::maybe_unoccupied,
                                  qnattr, do_not_throw);
    IndexSpace::register_instance(declab(L"Δ"),
                                  IndexSpace::complete_maybe_unoccupied, qnattr,
                                  do_not_throw);
  }
}

/// @brief creates an IndexRegistry
void make_default_indexregistry() {
  auto idxreg = std::make_shared<IndexRegistry>();
  auto& idxreg_ref = *idxreg;

  for (int s = 0; s <= 2; ++s) {
    auto declab = [s](auto&& label) {
      return qndecorate(static_cast<qns>(s), label);
    };
    register_index(idxreg_ref, Index{declab(L"u")}, 20);
    register_index(idxreg_ref, Index{declab(L"i")}, 100);
    register_index(idxreg_ref, Index{declab(L"I")}, 120);
    register_index(idxreg_ref, Index{declab(L"m")}, 110);
    register_index(idxreg_ref, Index{declab(L"a")}, 200);
    register_index(idxreg_ref, Index{declab(L"g")}, 800);
    register_index(idxreg_ref, Index{declab(L"e")}, 1000);
    register_index(idxreg_ref, Index{declab(L"A")}, 1020);
    register_index(idxreg_ref, Index{declab(L"x")}, 320);
    register_index(idxreg_ref, Index{declab(L"p")}, 1130);
    register_index(idxreg_ref, Index{declab(L"α'")}, 3000);
    register_index(idxreg_ref, Index{declab(L"c")}, 3800);
    register_index(idxreg_ref, Index{declab(L"α")}, 4000);
    register_index(idxreg_ref, Index{declab(L"κ")}, 4130);
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
      TensorCanonicalizer::set_cardinal_tensor_labels(
          mbpt::cardinal_tensor_labels());
    }
  }
}

}  // namespace mbpt
}  // namespace sequant
