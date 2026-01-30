#include "utils.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cassert>
#include <optional>

#include <range/v3/view/zip.hpp>

using namespace sequant;

std::size_t IndexSpaceMeta::getSize(const IndexSpace &space) const {
  return space.approximate_size();
}

std::size_t IndexSpaceMeta::getSize(const Index &index) const {
  return getSize(index.space());
}

std::string IndexSpaceMeta::getName(const IndexSpace &space) const {
  auto iter = m_entries.find(space);
  if (iter == m_entries.end()) {
    throw std::runtime_error("No known name for index space " +
                             toUtf8(space.base_key()));
  }

  return iter->second.name;
}

std::string IndexSpaceMeta::getTag(const IndexSpace &space) const {
  auto iter = m_entries.find(space);
  if (iter == m_entries.end()) {
    throw std::runtime_error("No known tag for index space " +
                             toUtf8(space.base_key()));
  }

  return iter->second.tag;
}

void IndexSpaceMeta::registerSpace(IndexSpace space, Entry entry) {
  m_entries.insert({std::move(space), std::move(entry)});
}

container::svector<container::svector<Index> > getExternalIndexPairs(
    const Tensor &tensor) {
  container::svector<container::svector<Index> > externalPairs;

  auto zipped = ranges::views::zip(tensor.bra(), tensor.ket());

  for (const auto &current : zipped) {
    externalPairs.push_back({current.first, current.second});
  }

  for (std::size_t i = zipped.size(); i < tensor.bra_rank(); ++i) {
    externalPairs.push_back({tensor.bra()[i]});
  }

  for (std::size_t i = zipped.size(); i < tensor.ket_rank(); ++i) {
    externalPairs.push_back({tensor.ket()[i]});
  }

  for (const Index &current : tensor.aux()) {
    externalPairs.push_back({current});
  }

  return externalPairs;
}

container::svector<container::svector<Index> > getExternalIndexPairs(
    const ExprPtr &expression) {
  if (expression->is<Tensor>()) {
    return getExternalIndexPairs(expression->as<Tensor>());
  }

  if (expression->is<Sum>()) {
    // External indices are expected to be the same across summands
    return getExternalIndexPairs(expression->as<Sum>().summand(0));
  }

  const Product &product = expression->as<Product>();

  if (product.at(0)->is<Tensor>() &&
      (product.at(0).as<Tensor>().label() == reserved::antisymm_label() ||
       product.at(0).as<Tensor>().label() == reserved::symm_label())) {
    // Special case in the presence of an (anti)symmetrizer operator
    return getExternalIndexPairs(product.at(0).as<Tensor>());
  }

  IndexGroups groups = get_unique_indices(expression);
  container::svector<container::svector<Index> > externalPairs;

  auto visitor = [&groups, &externalPairs](const ExprPtr &expr) {
    SEQUANT_ASSERT(expr.is<Tensor>());

    const Tensor &tensor = expr.as<Tensor>();

    for (std::size_t i = 0; i < std::min(tensor.bra_rank(), tensor.ket_rank());
         ++i) {
      const Index &braIdx = tensor.bra()[i];
      const Index &ketIdx = tensor.ket()[i];

      const bool braIsExternal = std::find(groups.bra.begin(), groups.bra.end(),
                                           braIdx) != groups.bra.end();
      const bool ketIsExternal = std::find(groups.ket.begin(), groups.ket.end(),
                                           ketIdx) != groups.ket.end();

      if (braIsExternal && ketIsExternal) {
        // Bra and ket index are paired -> keep them paired in external index
        // list
        externalPairs.push_back({braIdx, ketIdx});
      } else if (braIsExternal) {
        externalPairs.push_back({braIdx});
      } else if (ketIsExternal) {
        externalPairs.push_back({ketIdx});
      }
    }
  };

  expression->visit(visitor, true);

  for (Index &current : groups.aux) {
    externalPairs.push_back({std::move(current)});
  }

  return externalPairs;
}

bool needsSymmetrization(const sequant::ExprPtr &expression) {
  bool containsSymmetrizer = false;
  // Note: the assumption is that we never encounter cases where only part of
  // the expression requires symmetrization
  expression->visit(
      [&containsSymmetrizer](const ExprPtr &expr) {
        // Right now symmetrizer operators are represented as tensor objects
        // with name "S"
        if (expr.is<Tensor>() &&
            expr.as<Tensor>().label() == reserved::symm_label()) {
          containsSymmetrizer = true;
        }
      },
      true);

  return containsSymmetrizer;
}

sequant::ExprPtr generateResultSymmetrization(const ResultExpr &result,
                                              std::wstring_view precursorName) {
  IndexGroups<std::vector<Index> > externals;

  for (const auto &[bra, ket] :
       result.index_particle_grouping<std::pair<Index, Index> >()) {
    externals.bra.push_back(bra);
    externals.ket.push_back(ket);
  }

  return generateResultSymmetrization(precursorName, externals);
}

sequant::ExprPtr generateResultSymmetrization(const Tensor &result,
                                              std::wstring_view precursorName) {
  SEQUANT_ASSERT(result.bra_rank() == result.ket_rank());

  IndexGroups<std::vector<Index> > externals;

  for (const auto [bra, ket] : ranges::views::zip(result.bra(), result.ket())) {
    externals.bra.push_back(bra);
    externals.ket.push_back(ket);
  }

  return generateResultSymmetrization(precursorName, externals);
}

sequant::ExprPtr generateResultSymmetrization(
    std::wstring_view precursorName,
    const sequant::IndexGroups<std::vector<sequant::Index> > &externals) {
  SEQUANT_ASSERT(externals.bra.size() == externals.ket.size());

  // Note: we're only symmetrizing over bra-ket, not over auxiliary indices
  ExprPtr symmetrization = ex<Sum>(ExprPtrList{});
  for (std::size_t i = 0; i < externals.bra.size(); ++i) {
    std::vector<Index> symBra;
    std::vector<Index> symKet;

    for (std::size_t j = 0; j < externals.ket.size(); ++j) {
      symBra.push_back(externals.bra[(i + j) % externals.bra.size()]);
      symKet.push_back(externals.ket[(i + j) % externals.ket.size()]);
    }

    symmetrization += ex<Tensor>(precursorName, bra(std::move(symBra)),
                                 ket(std::move(symKet)), aux(externals.aux));
  }

  return symmetrization;
}

std::optional<ExprPtr> pop_symmetrizer(ResultExpr &expr) {
  std::optional<ExprPtr> symmetrizer =
      pop_tensor(expr.expression(), reserved::symm_label());

  if (!symmetrizer.has_value()) {
    symmetrizer = pop_tensor(expr.expression(), reserved::antisymm_label());
  }

  return symmetrizer;
}
