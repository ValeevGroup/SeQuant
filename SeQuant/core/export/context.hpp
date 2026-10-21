#ifndef SEQUANT_CORE_EXPORT_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_CONTEXT_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace sequant {

/// Different options for what to do when loading a tensor in the sense of
/// bringing it into memory to make it available for subsequent operations
enum class LoadStrategy {
  Load,
  Create,
};

/// Different options for when to (re)set a tensor's value to zero
enum class ZeroStrategy {
  NeverZero = 0,
  ZeroOnLoad = 1,
  ZeroOnCreate = 1 << 1,
  ZeroOnReuse = 1 << 2,
  AlwaysZero = ZeroOnLoad | ZeroOnCreate | ZeroOnReuse,
};

ZeroStrategy operator|(ZeroStrategy lhs, ZeroStrategy rhs);
ZeroStrategy operator&(ZeroStrategy lhs, ZeroStrategy rhs);
bool operator!(ZeroStrategy strategy);

/// A combination of a Load and ZeroStrategy
///
/// The two strategies are tracked independently: an entry only carries the
/// strategy that was explicitly set for it. This matters because setting one
/// strategy must not implicitly assert a value for the other -- e.g. requesting
/// a zero strategy for a result must not also pin its load strategy (and
/// thereby shadow a load strategy set elsewhere for the same object).
/// @sa LoadStrategy
/// @sa ZeroStrategy
struct StrategyPair {
  StrategyPair() = default;
  StrategyPair(LoadStrategy load, ZeroStrategy zero);
  StrategyPair(LoadStrategy load);
  StrategyPair(ZeroStrategy zero);

  std::optional<LoadStrategy> load;
  std::optional<ZeroStrategy> zero;
};

class AnyIrrepContext;
class AnyIrrepSet;

/// The origin of the IRREPs of a tensor (block) registered during export
enum class IrrepOrigin {
  /// IRREPs of a leaf tensor as obtained from the IRREP context. Repeated
  /// registrations must agree.
  Input,
  /// IRREPs of an intermediate as computed during export. Repeated
  /// registrations must agree.
  Intermediate,
  /// IRREPs of a result as computed during export. Repeated registrations
  /// (e.g. from different expressions contributing to the same result) are
  /// merged.
  Result,
};

namespace detail {

struct ExportIrrepData;

/// Value-semantic, opaque holder of the IRREP-related data of an
/// ExportContext (keeps the IRREP machinery out of this header)
class ExportIrrepDataHolder {
 public:
  ExportIrrepDataHolder();
  ExportIrrepDataHolder(const ExportIrrepDataHolder &other);
  ExportIrrepDataHolder(ExportIrrepDataHolder &&other) noexcept;
  ExportIrrepDataHolder &operator=(const ExportIrrepDataHolder &other);
  ExportIrrepDataHolder &operator=(ExportIrrepDataHolder &&other) noexcept;
  ~ExportIrrepDataHolder();

  /// @returns The held data or nullptr if there is none
  [[nodiscard]] const ExportIrrepData *get() const;
  /// @returns The held data (created if necessary)
  ExportIrrepData &get_or_create();

 private:
  std::unique_ptr<ExportIrrepData> m_data;
};

}  // namespace detail

/// Different ways a given tensor/variable can be used in a given code.
enum class Usage {
  None = 0,
  Result = 0b001,
  Intermediate = 0b010,
  Terminal = 0b100,
};

/// A set of Usage scenarios for a given tensor/variable. Contrary to Usage,
/// this can hold a combination of usage scenarios.
/// @sa Usage
class UsageSet {
 public:
  using Value = std::underlying_type_t<Usage>;
  UsageSet() = default;

  UsageSet &operator=(Usage usage);

  UsageSet &operator|=(UsageSet usage);
  UsageSet &operator&=(UsageSet usage);
  UsageSet &operator|=(Usage usage);
  UsageSet &operator&=(Usage usage);

  friend bool operator&(UsageSet set, Usage usage);
  friend bool operator&(Usage usage, UsageSet set);

  friend bool operator==(UsageSet set, Usage usage);
  friend bool operator==(Usage usage, UsageSet set);

  friend bool operator!=(UsageSet set, Usage usage);
  friend bool operator!=(Usage usage, UsageSet set);

 private:
  Value m_val = static_cast<Value>(Usage::None);
};

/// The ExportContext is intended to provide context during export operations.
/// This task is achieved by two core tasks:
/// 1. Acting as a meta-information storage. Before export starts, certain meta
///    information are collected and stored on the context object.
/// 2. Providing a customization point for users to influence any export
///    behavior that is not universal for a given backend. In this sense, the
///    context object can be seen as a configuration storage system.
/// The context object has to be user-provided when initiating an export and
/// is passed down to every single function during the export context. Hence,
/// it acts as a communication channel down into the export process.
class ExportContext {
 public:
  using TensorStrategyMap =
      std::map<Tensor, StrategyPair, TensorBlockLessThanComparator>;
  using VariableStrategyMap = std::map<Variable, StrategyPair>;

  ExportContext() = default;
  ExportContext(TensorStrategyMap tensorMap,
                VariableStrategyMap variableMap = {});
  ExportContext(VariableStrategyMap map);
  virtual ~ExportContext();

  /// @returns The LoadStrategy for the given tensor
  /// @sa LoadStrategy
  [[nodiscard]] virtual LoadStrategy loadStrategy(const Tensor &tensor) const;
  /// @returns The LoadStrategy for the given variable
  /// @sa LoadStrategy
  [[nodiscard]] virtual LoadStrategy loadStrategy(
      const Variable &variable) const;
  /// Sets the LoadStrategy for the given tensor
  /// @param tensor The Tensor to set the strategy for
  /// @param strategy The LoadStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa LoadStrategy
  /// @sa ExportNode
  virtual void setLoadStrategy(
      const Tensor &tensor, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  /// Sets the LoadStrategy for the given variable
  /// @param variable The Variable to set the strategy for
  /// @param strategy The LoadStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa LoadStrategy
  /// @sa ExportNode
  virtual void setLoadStrategy(
      const Variable &variable, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  /// @returns The ZeroStrategy for the given tensor
  /// @sa ZeroStrategy
  [[nodiscard]] virtual ZeroStrategy zeroStrategy(const Tensor &tensor) const;
  /// @returns The ZeroStrategy for the given variable
  /// @sa ZeroStrategy
  [[nodiscard]] virtual ZeroStrategy zeroStrategy(
      const Variable &variable) const;
  /// Sets the ZeroStrategy for the given tensor
  /// @param tensor The Tensor to set the strategy for
  /// @param strategy The ZeroStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa ZeroStrategy
  /// @sa ExportNode
  virtual void setZeroStrategy(
      const Tensor &tensor, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  /// Sets the ZeroStrategy for the given variable
  /// @param variable The Variable to set the strategy for
  /// @param strategy The ZeroStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa ZeroStrategy
  /// @sa ExportNode
  virtual void setZeroStrategy(
      const Variable &variable, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  /// Performs necessary (backend-specific) modifications (if any) on the
  /// given Tensor in-place
  /// @param tensor The Tensor to be modified
  /// @returns Whether any modification has been performed
  virtual bool rewrite(Tensor &tensor) const;
  /// Performs necessary (backend-specific) modifications (if any) on the
  /// given Variable in-place
  /// @param variable The Variable to be modified
  /// @returns Whether any modification has been performed
  virtual bool rewrite(Variable &variable) const;

  /// @returns Whether the current export workflow is generating code inside a
  /// named section
  [[nodiscard]] virtual bool inside_named_section() const;
  /// @returns The name of the current section. If not inside a named section,
  /// this function will throw.
  [[nodiscard]] virtual const std::string &current_section_name() const;
  /// Sets the name of the section that subsequent export shall be for.
  /// Implicitly notifies the context that from now on we're in a named section.
  /// @param name The name of the current section
  virtual void set_current_section_name(std::string name);
  /// Resets the current section name. Implicitly notifies the context that
  /// we're no longer inside a named section.
  virtual void clear_current_section_name();

  /// @returns Whether or not the ID of the current expression is available
  [[nodiscard]] virtual bool has_current_expression_id() const;
  /// @returns The ID of the current expression. If none is available, this
  /// function will throw.
  [[nodiscard]] virtual std::size_t current_expression_id() const;
  /// Sets the current expression ID
  /// @param id The ID of the current expression
  virtual void set_current_expression_id(std::size_t id);
  /// Resets the ID of the current expression
  virtual void clear_current_expression_id();

  /// @param id The ID of the relevant expression
  /// @returns The set of indices that the given expression should be batched
  /// over
  virtual std::span<const Index> batch_indices(
      std::optional<std::size_t> id = {}) const;
  /// Sets the batch indices for the given expression
  /// @param indices The batch indices
  /// @param id ID of the target expression
  virtual void set_batch_indices(std::span<const Index> indices,
                                 std::optional<std::size_t> id = {});
  /// @param id ID of the relevant expression
  /// @returns Whether this expression is computed in index batches
  virtual bool is_batched(std::optional<std::size_t> id = {}) const;
  /// @param idx The index to check
  /// @param id ID of the relevant expression
  /// @returns Whether this expression uses idx as a batch index
  virtual bool is_index_batched(const Index &idx,
                                std::optional<std::size_t> id = {}) const;

  /// Sets the context used to determine the IRREPs of all tensors occurring
  /// during export. Without such a context, no IRREP information is
  /// determined.
  virtual void set_irrep_context(AnyIrrepContext ctx);
  /// @returns Whether an IRREP context has been set
  [[nodiscard]] virtual bool has_irrep_context() const;
  /// @returns The IRREP context. Throws if none has been set.
  [[nodiscard]] virtual const AnyIrrepContext &irrep_context() const;
  /// @returns The IRREPs of the given tensor (block) as registered during the
  /// export or nullptr if they are unknown
  [[nodiscard]] virtual const AnyIrrepSet *irreps_of(
      const Tensor &tensor) const;
  /// @returns The IRREPs of the given tensor (block) if it has been computed
  /// (as an intermediate or result) during the export. Otherwise, nullptr.
  [[nodiscard]] virtual const AnyIrrepSet *computed_irreps_of(
      const Tensor &tensor) const;
  /// Registers the IRREPs of the given tensor (block)
  /// @param tensor The tensor
  /// @param irreps The IRREPs the tensor spans
  /// @param origin The origin of the IRREPs, which determines how repeated
  ///        registrations of the same tensor block are handled
  /// @throws Exception if the IRREPs conflict with a previous registration
  /// @sa IrrepOrigin
  virtual void register_irreps(const Tensor &tensor, AnyIrrepSet irreps,
                               IrrepOrigin origin);

 protected:
  static constexpr std::size_t ID_GLOBAL =
      std::numeric_limits<std::size_t>::max();

 private:
  std::map<std::size_t, TensorStrategyMap> m_tensorStrategies;
  std::map<std::size_t, VariableStrategyMap> m_variableStrategies;
  std::optional<std::string> m_currentSection;
  std::optional<std::size_t> m_currentExpressionID;
  detail::ExportIrrepDataHolder m_irrepData;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_CONTEXT_HPP
