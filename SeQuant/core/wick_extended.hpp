#ifndef SEQUANT_CORE_WICK_EXTENDED_HPP
#define SEQUANT_CORE_WICK_EXTENDED_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>

#include <cstddef>
#include <optional>
#include <utility>

namespace sequant {

/// controls extended_wick() and cumulant_expand()
struct ExtendedWickOptions {
  /// keep only terms with no surviving operators
  bool full_contractions = true;
  /// largest cumulant rank to form; nullopt = no bound; 0 or 1 = no cumulants
  std::optional<std::size_t> max_cumulant_rank;
  /// rewrite every η as δ - γ
  bool eta_as_delta_minus_gamma = false;
  /// forwarded to WickTheorem::use_topology
  bool use_topology = true;
  /// pairs of input NormalOperator ordinals that must end up connected
  container::svector<std::pair<std::size_t, std::size_t>> nop_connections;
  /// pairs of input NormalOperator ordinals that must not be connected
  container::svector<std::pair<std::size_t, std::size_t>>
      nop_avoided_connections;
};

/// maps the index of a surviving Op to the ordinal of the input
/// NormalOperator it came from
using OpProvenance = container::map<Index, std::size_t>;

namespace detail {

/// the value of a cumulant block; a spin-free variant would override this
template <Statistics S>
ExprPtr block_value(const NormalOperator<S> &block) {
  return density::make_cumulant(block);
}

/// a per-term scalar weight (1 for spin-orbital); a spin-free variant would
/// put its cycle factor here
template <Statistics S>
rational term_weight(
    const NormalOperator<S> & /*survivors*/,
    const container::svector<container::svector<std::size_t>> & /*blocks*/) {
  return 1;
}

}  // namespace detail

/// expands the surviving operators of each term of a MultiProduct-vacuum
/// WickTheorem output (partial contractions) into cumulant blocks
/// @param wick_output the WickTheorem output
/// @param provenance the input NormalOperator ordinal of every surviving Op
/// @param opts only `full_contractions` and `max_cumulant_rank` are used
/// @note a block has k creators and k annihilators, 2 <= k <=
///       `max_cumulant_rank`, all active, and legs from at least two input
///       NormalOperators; its sign is the parity of moving each block's legs,
///       in order, to the front of the surviving operator string, and the
///       remaining operators (if `!opts.full_contractions`) are kept as a
///       MultiProduct-vacuum NormalOperator
/// @pre every γ/η index is pure-active, and every surviving index is
///      pure-active, pure core or pure virtual
template <Statistics S>
ExprPtr cumulant_expand(const ExprPtr &wick_output,
                        const OpProvenance &provenance,
                        const ExtendedWickOptions &opts);

extern template ExprPtr cumulant_expand<Statistics::FermiDirac>(
    const ExprPtr &, const OpProvenance &, const ExtendedWickOptions &);

}  // namespace sequant

#endif  // SEQUANT_CORE_WICK_EXTENDED_HPP
