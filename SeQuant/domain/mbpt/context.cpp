#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/context.hpp>

#include <range/v3/algorithm/contains.hpp>

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
#include <mutex>
#endif

namespace sequant::mbpt {

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
static std::recursive_mutex mbpt_ctx_mtx;  // used to protect the MBPT context
#endif

Context::Context(Options options)
    : csv_(options.csv),
      op_registry_(options.op_registry_ptr
                       ? std::move(options.op_registry_ptr)
                       : (options.op_registry
                              ? std::make_shared<OpRegistry>(
                                    std::move(options.op_registry.value()))
                              : nullptr)) {}

Context Context::clone() const {
  Context ctx(*this);
  if (op_registry_) {
    ctx.op_registry_ = std::make_shared<OpRegistry>(op_registry_->clone());
  }
  return ctx;
}

CSV Context::csv() const { return csv_; }

std::shared_ptr<const OpRegistry> Context::op_registry() const {
  SEQUANT_ASSERT(op_registry_, "mbpt::Context has null OpRegistry");
  return op_registry_;
}

std::shared_ptr<OpRegistry> Context::mutable_op_registry() const {
  SEQUANT_ASSERT(op_registry_, "mbpt::Context has null OpRegistry");
  return op_registry_;
}

Context& Context::set(const OpRegistry& op_registry) {
  op_registry_ = std::make_shared<OpRegistry>(op_registry);
  return *this;
}

Context& Context::set(std::shared_ptr<OpRegistry> op_registry) {
  op_registry_ = std::move(op_registry);
  return *this;
}

Context& Context::set(CSV csv) {
  csv_ = csv;
  return *this;
}

bool operator==(Context const& left, Context const& right) {
  if (left.csv() != right.csv()) return false;

  // both null -> equal; one null -> not equal
  if (!left.op_registry_ && !right.op_registry_) return true;
  if (!left.op_registry_ || !right.op_registry_) return false;

  return *left.op_registry_ == *right.op_registry_;
}

bool operator!=(Context const& left, Context const& right) {
  return !(left == right);
}

const Context& get_default_mbpt_context() {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(mbpt_ctx_mtx);
#endif
  return sequant::detail::get_implicit_context<Context>();
}

void set_default_mbpt_context(const Context& ctx) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(mbpt_ctx_mtx);
#endif
  sequant::detail::set_implicit_context(ctx);
}

void set_default_mbpt_context(const Context::Options& options) {
  return set_default_mbpt_context(Context(options));
}

void reset_default_mbpt_context() {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(mbpt_ctx_mtx);
#endif
  sequant::detail::reset_implicit_context<Context>();
}

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context& f) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(mbpt_ctx_mtx);
#endif
  return sequant::detail::set_scoped_implicit_context(f);
}

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context::Options& f) {
  return set_scoped_default_mbpt_context(Context(f));
}

std::shared_ptr<OpRegistry> make_minimal_registry() {
  auto registry = std::make_shared<OpRegistry>();

  registry
      ->add(L"h", OpClass::Gen)   /// 1-body Hamiltonian
      .add(L"g", OpClass::Gen)    /// 2-body Coulomb
      .add(L"f", OpClass::Gen)    /// Fock operator
      .add(L"θ", OpClass::Gen)    /// general fock space operator
      .add(L"t", OpClass::Ex)     /// cluster operator
      .add(L"λ", OpClass::Deex)   /// deexcitation cluster operator
      .add(L"R", OpClass::Ex)     /// right-hand eigenstate
      .add(L"L", OpClass::Deex);  /// left-hand eigenstate

  return registry;
}

std::shared_ptr<OpRegistry> make_legacy_registry() {
  auto registry = std::make_shared<OpRegistry>();

  registry->add(L"h", OpClass::Gen)
      .add(L"f", OpClass::Gen)
      /// closed Fock operator (i.e. Fock operator due to fully-occupied
      /// orbitals)
      .add(L"f̃", OpClass::Gen)
      .add(L"g", OpClass::Gen)
      .add(L"θ", OpClass::Gen)
      .add(L"t", OpClass::Ex)
      .add(L"λ", OpClass::Deex)
      .add(L"R", OpClass::Ex)
      .add(L"L", OpClass::Deex)
      /// R12
      .add(L"F", OpClass::Gen)
      .add(L"GR", OpClass::Gen)
      .add(L"C", OpClass::Gen)
      /// RDM and RDM Cumulant
      .add(L"γ", OpClass::Gen)
      .add(L"κ", OpClass::Gen);

  return registry;
}

OpClass to_op_class(const std::wstring& op) {
  // check reserved labels first
  if (ranges::contains(reserved::labels(), op)) {
    return OpClass::Gen;  // all reserved labels are Gen
  } else {
    return get_default_mbpt_context().op_registry()->to_class(op);
  }
}

Hermiticity op_hermiticity(const std::wstring& op) {
  // reserved labels are OpClass::Gen, hence Hermitian by default
  if (ranges::contains(reserved::labels(), op)) {
    return default_hermiticity(OpClass::Gen);
  } else {
    return get_default_mbpt_context().op_registry()->hermiticity(op);
  }
}

}  // namespace sequant::mbpt
