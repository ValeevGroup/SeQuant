#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

int main() {
  using namespace sequant;

  // start-snippet-1
  // the one-line way to configure SeQuant for standard single-reference
  // quantum chemistry: populates the default IndexBasisRegistry with the
  // usual occupied/virtual partitioning and switches the default Context to
  // the single-product (quasiparticle) vacuum
  mbpt::load(mbpt::Convention::SR, mbpt::SpinConvention::None);

  const Context& ctx = get_default_context();
  SEQUANT_ASSERT(ctx.vacuum() == Vacuum::SingleProduct);
  SEQUANT_ASSERT(ctx.index_basis_registry() != nullptr);
  // end-snippet-1

  // start-snippet-2
  // most of SeQuant/domain/mbpt additionally needs an OpRegistry, which maps
  // operator labels ("t", "f", "g", ...) to their OpClass (excitation,
  // de-excitation, or general); this lives in a separate, mbpt-specific
  // context that must be configured on top of the core one
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_legacy_registry()});

  SEQUANT_ASSERT(mbpt::to_op_class(L"t") == mbpt::OpClass::Ex);
  // end-snippet-2

  // start-snippet-4
  // CSV controls whether OpMaker builds excitation/de-excitation operators
  // (e.g. cluster amplitudes "t") with cluster-specific (index-dependent)
  // virtuals or plain, independent ones; it defaults to CSV::No
  auto csv_resetter = mbpt::set_scoped_default_mbpt_context(
      {.csv = mbpt::CSV::Yes, .op_registry_ptr = mbpt::make_legacy_registry()});
  SEQUANT_ASSERT(mbpt::get_default_mbpt_context().csv() == mbpt::CSV::Yes);
  // end-snippet-4

  // start-snippet-6
  // several bases of one space can meet in one expression, e.g. canonical and
  // localized occupied orbitals. An Index made from an IndexSpace runs over
  // the space's own basis; one made from an IndexBasis with a basis instance
  // runs over that basis of the space. The instance is written after the
  // proto indices: i_1<;1>
  const auto& isr = get_default_context().index_basis_registry();
  const IndexSpace occ = isr->retrieve(L"i");
  const Index i_canonical(occ, 1);
  const Index i_localized(IndexBasis(occ, 1), 1);
  SEQUANT_ASSERT(i_canonical.space() == i_localized.space());
  SEQUANT_ASSERT(i_canonical != i_localized);
  SEQUANT_ASSERT(i_localized.full_label() == L"i_1<;1>");

  // a (de)excitation operator's legs take an instance through its registry
  // entry: here the unoccupied legs of t are in basis instance 1, so the
  // amplitudes OpMaker builds, and the projectors of the equations solved for
  // them, carry it
  auto granted_registry = mbpt::make_legacy_registry();
  granted_registry->grant_basis(L"t", isr->retrieve(L"a"), 1);
  auto grant_resetter = mbpt::set_scoped_default_mbpt_context(
      {.op_registry_ptr = granted_registry});
  const ExprPtr t1 = mbpt::op::tensor::t(1);
  for (const Index& idx : get_used_indices(t1))
    SEQUANT_ASSERT(idx.basis().basis_instance() ==
                   (isr->is_pure_unoccupied(idx.space())
                        ? IndexBasis::optional_instance{1}
                        : IndexBasis::optional_instance{}));
  // end-snippet-6

  // start-snippet-3
  // temporarily switching context (e.g. to a spin-free basis for a single
  // calculation) without disturbing the enclosing default: the RAII
  // resetter restores the previous context when it goes out of scope
  {
    auto resetter = set_scoped_default_context({.spbasis = SPBasis::Spinfree});
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinfree);
  }
  SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);
  // end-snippet-3

  // start-snippet-5
  // the Context also owns the canonicalizer configuration: within the scope,
  // tensors labeled "A" are canonicalized by the null canonicalizer, so their
  // slots keep the order the graph-based canonicalization of the product gives
  // them, here the input order (hence no phase is produced); a lone
  // Tensor::canonicalize() would use the canonicalizer of the empty label
  // instead
  auto make_product = [] {
    return ex<Tensor>(L"A", bra{L"a_1", L"i_1"}, ket{L"i_2", L"a_2"},
                      Symmetry::Antisymm) *
           ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_3"});
  };
  {
    auto resetter = set_scoped_modified_default_context([](Context& scoped) {
      scoped.set_tensor_canonicalizer(
          L"A", std::make_shared<NullTensorCanonicalizer>());
    });
    auto product = make_product();
    product->canonicalize();
    SEQUANT_ASSERT(product.as<Product>().scalar() == 1);
  }

  // outside of the scope the default canonicalizer reorders the indices of A,
  // picking up a sign
  auto product = make_product();
  product->canonicalize();
  SEQUANT_ASSERT(product.as<Product>().scalar() == -1);
  // end-snippet-5

  // start-snippet-7
  {
    auto normalization_resetter = mbpt::set_scoped_default_mbpt_context(
        mbpt::Context{mbpt::get_default_mbpt_context()}.set(
            mbpt::NormalizationConvention::Symmetric));
    SEQUANT_ASSERT(
        mbpt::get_default_mbpt_context().normalization_convention() ==
        mbpt::NormalizationConvention::Symmetric);
  }
  SEQUANT_ASSERT(mbpt::get_default_mbpt_context().normalization_convention() ==
                 mbpt::NormalizationConvention::Default);
  // end-snippet-7

  (void)ctx;

  return 0;
}
