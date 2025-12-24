#include <SeQuant/domain/mbpt/context.hpp>

namespace sequant::mbpt {

/// forward declaration of OpClass
enum class OpClass;

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
  ctx.op_registry_ = std::make_shared<OpRegistry>(op_registry_->clone());
  return ctx;
}

CSV Context::csv() const { return csv_; }

std::shared_ptr<const OpRegistry> Context::op_registry() const {
  return op_registry_;
}

std::shared_ptr<OpRegistry> Context::mutable_op_registry() {
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
  return left.csv() == right.csv() &&
         *(left.op_registry()) == *(right.op_registry());
}

bool operator!=(Context const& left, Context const& right) {
  return !(left == right);
}

const Context& get_default_mbpt_context() {
  return sequant::detail::get_implicit_context<Context>();
}

void set_default_mbpt_context(const Context& ctx) {
  sequant::detail::set_implicit_context(ctx);
}

void set_default_mbpt_context(const Context::Options& options) {
  return set_default_mbpt_context(Context(options));
}

void reset_default_mbpt_context() {
  sequant::detail::reset_implicit_context<Context>();
}

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context& f) {
  return sequant::detail::set_scoped_implicit_context(f);
}

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context::Options& f) {
  return sequant::detail::set_scoped_implicit_context(Context(f));
}

std::shared_ptr<OpRegistry> make_minimal_registry() {
  auto registry = std::make_shared<OpRegistry>();

  registry
      ->add(L"h", OpClass::gen)   /// 1-body Hamiltonian
      .add(L"g", OpClass::gen)    /// 2-body Coulomb
      .add(L"f", OpClass::gen)    /// Fock operator
      .add(L"θ", OpClass::gen)    /// general fock space operator
      .add(L"t", OpClass::ex)     /// cluster operator
      .add(L"λ", OpClass::deex)   /// deexcitation cluster operator
      .add(L"R", OpClass::ex)     /// right-hand eigenstate
      .add(L"L", OpClass::deex);  /// left-hand eigenstate

  return registry;
}

std::shared_ptr<OpRegistry> make_legacy_registry() {
  auto registry = std::make_shared<OpRegistry>();

  registry->add(L"h", OpClass::gen)
      .add(L"f", OpClass::gen)
      /// closed Fock operator (i.e. Fock operator due to fully-occupied
      /// orbitals)
      .add(L"f̃", OpClass::gen)
      .add(L"g", OpClass::gen)
      .add(L"θ", OpClass::gen)
      .add(L"t", OpClass::ex)
      .add(L"λ", OpClass::deex)
      .add(L"R", OpClass::ex)
      .add(L"L", OpClass::deex)
      /// R12
      .add(L"F", OpClass::gen)
      .add(L"GR", OpClass::gen)
      .add(L"C", OpClass::gen)
      /// RDM and RDM Cumulant
      .add(L"γ", OpClass::gen)
      .add(L"κ", OpClass::gen);

  return registry;
}
}  // namespace sequant::mbpt
