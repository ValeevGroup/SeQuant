//
// Created by Robert Adam on 2023-09-04
//

#include "itf.hpp"

#include <SeQuant/core/utility/indices.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cwchar>
#include <functional>
#include <map>
#include <unordered_map>

namespace sequant {

std::wstring to_itf(const itf::CodeBlock &block, const itf::Context &ctx) {
  itf::detail::ITFGenerator generator(ctx);
  generator.addBlock(block);
  return generator.generate();
}

namespace itf {

Context::~Context() {}

Result::Result(ExprPtr expression, Tensor resultTensor, bool importResultTensor)
    : expression(std::move(expression)),
      resultTensor(std::move(resultTensor)),
      importResultTensor(importResultTensor) {}

Tensor generateResultTensor(ExprPtr expr) {
  // This (same as ITF itself) assumes that all indices that are repeated in the
  // expression are also contracted (summed) over and therefore don't contribute
  // to the result of the expression
  IndexGroups externals = get_unique_indices(expr);

  return Tensor(L"Result", bra(std::move(externals.bra)),
                ket(std::move(externals.ket)), aux(std::move(externals.aux)));
}

Result::Result(ExprPtr expression, bool importResultTensor)
    : expression(std::move(expression)),
      resultTensor(generateResultTensor(expression)),
      importResultTensor(importResultTensor) {}

CodeBlock::CodeBlock(std::wstring blockName, Result result)
    : CodeBlock(std::move(blockName), std::vector<Result>{std::move(result)}) {}

CodeBlock::CodeBlock(std::wstring blockName, std::vector<Result> results)
    : name(std::move(blockName)), results(std::move(results)) {}

namespace detail {

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor,
                                         const Context &ctx);

Tensor to_tensor(const Expr &expr) {
  if (expr.is<Tensor>()) {
    return expr.as<Tensor>();
  }

  if (expr.is<Variable>()) {
    // A variable can be seen as a tensor without indices
    return Tensor(std::wstring(expr.as<Variable>().label()), {}, {}, {});
  }

  throw std::runtime_error(
      "Invalid Expr type encountered (can't be converted to a Tensor) in "
      "to_tensor");
}

std::vector<Contraction> to_contractions(const Product &product,
                                         const Tensor &resultTensor,
                                         const Context &ctx) {
  static std::size_t intermediateCounter = 1;

  if (product.factors().size() == 1) {
    assert(product.scalar().imag() == 0);

    return {Contraction{product.scalar().real(),
                        resultTensor,
                        to_tensor(*product.factor(0)),
                        {}}};
  }

  // We assume that we're dealing with a binary tree
  assert(product.factors().size() == 2);

  std::unordered_map<const ExprPtr::element_type *, Tensor> intermediates;
  std::vector<Contraction> contractions;

  // Handle intermediates
  for (const ExprPtr &factor : product.factors()) {
    if (factor.is<Product>()) {
      // Create intermediate that computes this nested product
      // This (same as ITF itself) assumes that all indices that are repeated in
      // the expression are also contracted (summed) over and therefore don't
      // contribute to the result of the expression
      IndexGroups intermediateIndexGroups = get_unique_indices(factor);

      // Collect all intermediate indices and sort them such that the order of
      // index spaces is the canonical one within ITF
      // This is possible, because the index ordering of intermediates is
      // arbitrary
      std::vector<Index> intermediateIndices;
      intermediateIndices.reserve(intermediateIndexGroups.bra.size() +
                                  intermediateIndexGroups.ket.size() +
                                  intermediateIndexGroups.aux.size());
      intermediateIndices.insert(intermediateIndices.end(),
                                 intermediateIndexGroups.bra.begin(),
                                 intermediateIndexGroups.bra.end());
      intermediateIndices.insert(intermediateIndices.end(),
                                 intermediateIndexGroups.ket.begin(),
                                 intermediateIndexGroups.ket.end());
      intermediateIndices.insert(intermediateIndices.end(),
                                 intermediateIndexGroups.aux.begin(),
                                 intermediateIndexGroups.aux.end());
      std::sort(intermediateIndices.begin(), intermediateIndices.end(),
                [&ctx](const Index &lhs, const Index &rhs) {
                  int result = ctx.compare(lhs, rhs);
                  return result == 0 ? lhs < rhs : result < 0;
                });

      std::array<wchar_t, 64> intermediateName;
      swprintf(intermediateName.data(), intermediateName.size(), L"INTER%06u",
               intermediateCounter++);

      // There is no notion of bra and ket for intermediates, so we dump all
      // indices in the bra for now
      Tensor intermediate(intermediateName.data(),
                          bra(std::move(intermediateIndices)),
                          ket(std::vector<Index>{}));

      std::vector<Contraction> intermediateContractions =
          to_contractions(factor, intermediate, ctx);
      contractions.reserve(contractions.size() +
                           intermediateContractions.size());
      contractions.insert(
          contractions.end(),
          std::make_move_iterator(intermediateContractions.begin()),
          std::make_move_iterator(intermediateContractions.end()));

      intermediates.insert({factor.get(), std::move(intermediate)});
    } else if (factor.is<Sum>()) {
      // TODO: Handle on-the-fly antisymmetrization (K[abij] - K[baij])
      throw std::invalid_argument(
          "Products of sums can not yet be translated to ITF");
    }
  }

  // Now create the contraction for the two factors
  auto lhsIntermediate = intermediates.find(product.factor(0).get());
  auto rhsIntermediate = intermediates.find(product.factor(1).get());

  assert(product.scalar().imag() == 0);
  contractions.push_back(Contraction{
      product.scalar().real(), resultTensor,
      lhsIntermediate == intermediates.end() ? to_tensor(*product.factor(0))
                                             : lhsIntermediate->second,
      rhsIntermediate == intermediates.end() ? to_tensor(*product.factor(1))
                                             : rhsIntermediate->second});

  return contractions;
}

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor,
                                         const Context &ctx) {
  std::wstring itfCode;

  if (expression.is<Constant>()) {
    // Make use of special One[] tensor to represent adding constants
    return {Contraction{expression.as<Constant>().value().real(), resultTensor,
                        Tensor(L"One", {}, {}, {})}};
  } else if (expression.is<Tensor>()) {
    return {Contraction{1, resultTensor, expression.as<Tensor>(), {}}};
  } else if (expression.is<Product>()) {
    // Separate into binary contractions
    return to_contractions(expression.as<Product>(), resultTensor, ctx);
  } else if (expression.is<Sum>()) {
    // Process each summand
    std::vector<Contraction> contractions;

    for (const ExprPtr &summand : expression.as<Sum>().summands()) {
      std::vector<Contraction> currentContractions =
          to_contractions(summand, resultTensor, ctx);

      contractions.reserve(contractions.size() + currentContractions.size());
      contractions.insert(contractions.end(),
                          std::make_move_iterator(currentContractions.begin()),
                          std::make_move_iterator(currentContractions.end()));
    }

    return contractions;
  } else {
    throw std::invalid_argument(
        "Unhandled expression type in to_contractions function");
  }
}

template <typename IndexContainer, typename SpaceTypeContainer>
bool isSpacePattern(const IndexContainer &indices,
                    const SpaceTypeContainer &pattern) {
  static_assert(std::is_same_v<typename IndexContainer::value_type, Index>);
  static_assert(std::is_same_v<typename SpaceTypeContainer::value_type,
                               IndexSpace::Type>);
  assert(indices.size() == pattern.size());

  auto indexIter = indices.cbegin();
  auto typeIter = pattern.cbegin();

  while (indexIter != indices.cend() && typeIter != pattern.cend()) {
    if (indexIter->space().type() != *typeIter) {
      return false;
    }

    ++indexIter;
    ++typeIter;
  }

  return true;
}

void one_electron_integral_remapper(ExprPtr &expr,
                                    const std::wstring_view integralTensorLabel,
                                    const Context &ctx) {
  if (!expr.is<Tensor>()) {
    return;
  }

  const Tensor &tensor = expr.as<Tensor>();

  if (tensor.label() != integralTensorLabel || tensor.bra().size() != 1 ||
      tensor.ket().size() != 1) {
    return;
  }

  assert(tensor.bra().size() == 1);
  assert(tensor.ket().size() == 1);

  auto braIndices = tensor.bra();
  auto ketIndices = tensor.ket();

  // Use the bra-ket (hermitian) symmetry of the integrals to exchange creators
  // and annihilators such that the larger index space is on the left (in the
  // bra)
  int res = ctx.compare(ketIndices[0], braIndices[0]);
  if (res < 0) {
    std::swap(braIndices[0], ketIndices[0]);
  } else if (res == 0 && ketIndices[0] < braIndices[0]) {
    // Cosmetic exchange to arrive at a more canonical index ordering
    std::swap(braIndices[0], ketIndices[0]);
  }

  expr = ex<Tensor>(tensor.label(), bra(std::move(braIndices)),
                    ket(std::move(ketIndices)), aux(tensor.aux()));
}

template <typename BraContainer, typename KetContainer>
bool isExceptionalJ(const BraContainer &braIndices,
                    const KetContainer &ketIndices, const Context &ctx) {
  assert(braIndices.size() == 2);
  assert(ketIndices.size() == 2);
  // integrals with 3 external (virtual) indices ought to be converted to
  // J-integrals
  // Here, we generalize this to all integrals for which the bra indices and
  // the first ket index are of the same space and that space compares less
  // than the space of the bras.
  // TODO: outsource determining of exceptional Js to context
  const bool bras_are_same = ctx.compare(braIndices[0], braIndices[1]) == 0;
  const bool kets_are_ordered = ctx.compare(ketIndices[0], ketIndices[1]) < 0;
  const bool first_ket_same_as_bra =
      ctx.compare(ketIndices[0], braIndices[0]) == 0;
  return bras_are_same && kets_are_ordered && first_ket_same_as_bra;
}

void two_electron_integral_remapper(ExprPtr &expr,
                                    const std::wstring_view integralTensorLabel,
                                    const Context &ctx) {
  if (!expr.is<Tensor>()) {
    return;
  }

  const Tensor &tensor = expr.as<Tensor>();

  if (tensor.label() != integralTensorLabel || tensor.bra().size() != 2 ||
      tensor.ket().size() != 2) {
    return;
  }

  assert(tensor.bra().size() == 2);
  assert(tensor.ket().size() == 2);

  // Copy indices as we might have to mutate them
  auto braIndices = tensor.bra();
  auto ketIndices = tensor.ket();
  assert(tensor.aux().empty());

  // Step 1: Use 8-fold permutational symmetry of spin-summed integrals
  // to bring indices into a canonical order in terms of the index
  // spaces they belong to. Note: This symmetry is generated by the two
  // individual bra-ket symmetries for indices for particle one and two
  // as well as the particle-1,2-symmetry (column-symmetry)
  //
  // The final goal is to order the indices in descending index space
  // size, where the assumed relative sizes are
  // occ < virt

  // Step 1a: Particle-intern bra-ket symmetry
  for (std::size_t i = 0; i < braIndices.size(); ++i) {
    if (ctx.compare(ketIndices.at(i), braIndices.at(i)) < 0) {
      // This bra index belongs to a smaller space than the ket index ->
      // swap them
      std::swap(braIndices[i], ketIndices[i]);
    }
  }

  // Step 1b: Particle-1,2-symmetry
  bool switchColumns = false;
  if (int res = ctx.compare(braIndices[1], braIndices[0]); res != 0) {
    switchColumns = res < 0;
  } else if (int res = ctx.compare(ketIndices[1], ketIndices[0]); res != 0) {
    switchColumns = res < 0;
  }

  if (switchColumns) {
    std::swap(braIndices[0], braIndices[1]);
    std::swap(ketIndices[0], ketIndices[1]);
  }

  // Step 2: Look at the index space patterns to figure out whether
  // this is a K or a J integral. If the previously attempted sorting
  // of index spaces can be improved by switching the second and third
  // index, do that and thereby produce a J tensor. Otherwise, we retain
  // the index sequence as-is and thereby produce a K tensor.
  // There are some explicit exceptions for J-tensors though.
  std::wstring tensorLabel;
  Index *particle1_1 = nullptr;
  Index *particle1_2 = nullptr;
  Index *particle2_1 = nullptr;
  Index *particle2_2 = nullptr;

  if (isExceptionalJ(braIndices, ketIndices, ctx) ||
      ctx.compare(ketIndices[0], braIndices[1]) < 0) {
    std::swap(braIndices[1], ketIndices[0]);
    tensorLabel = L"J";

    particle1_1 = &braIndices[0];
    particle1_2 = &braIndices[1];
    particle2_1 = &ketIndices[0];
    particle2_2 = &ketIndices[1];
  } else {
    tensorLabel = L"K";

    particle1_1 = &braIndices[0];
    particle1_2 = &ketIndices[0];
    particle2_1 = &braIndices[1];
    particle2_2 = &ketIndices[1];
  }

  // Go through the symmetries again to try and produce the most canonical
  // index ordering possible without breaking the index-space ordering
  // established up to this point.
  // This is a purely cosmetic change, but it is very useful for testing
  // purposes to have  unique representation of the integrals.
  if (particle1_1->space().type() == particle1_2->space().type() &&
      *particle1_2 < *particle1_1) {
    std::swap(*particle1_1, *particle1_2);
  }
  if (particle2_1->space().type() == particle2_2->space().type() &&
      *particle2_2 < *particle2_1) {
    std::swap(*particle2_1, *particle2_2);
  }
  if (particle1_1->space().type() == particle2_1->space().type() &&
      particle1_2->space().type() == particle2_2->space().type()) {
    if (*particle2_1 < *particle1_1 ||
        (*particle1_1 == *particle2_1 && *particle2_2 < *particle1_2)) {
      std::swap(*particle1_1, *particle2_1);
      std::swap(*particle1_2, *particle2_2);
    }
  }

  expr = ex<Tensor>(std::move(tensorLabel), bra(std::move(braIndices)),
                    ket(std::move(ketIndices)), tensor.aux());
}

void integral_remapper(ExprPtr &expr, const Context &ctx,
                       std::wstring_view oneElectronIntegralName,
                       std::wstring_view twoElectronIntegralName,
                       std::wstring_view dfTensorName) {
  two_electron_integral_remapper(expr, twoElectronIntegralName, ctx);
  one_electron_integral_remapper(expr, oneElectronIntegralName, ctx);
  // We can reuse the same logic for the DF tensors
  one_electron_integral_remapper(expr, dfTensorName, ctx);
}

void remap_integrals(ExprPtr &expr, const Context &ctx,
                     std::wstring_view oneElectronIntegralName,
                     std::wstring_view twoElectronIntegralName,
                     std::wstring_view dfTensorName) {
  auto remapper = [&](ExprPtr &expr) {
    return integral_remapper(expr, ctx, oneElectronIntegralName,
                             twoElectronIntegralName, dfTensorName);
  };

  const bool visitedRoot = expr->visit(remapper, true);

  if (!visitedRoot) {
    remapper(expr);
  }
}

ITFGenerator::ITFGenerator(const Context &ctx) : m_ctx(&ctx) {}

void ITFGenerator::addBlock(const itf::CodeBlock &block) {
  m_codes.reserve(m_codes.size() + block.results.size());

  std::vector<std::vector<Contraction>> contractionBlocks;

  for (const Result &currentResult : block.results) {
    ExprPtr expression = currentResult.expression;
    remap_integrals(expression, *m_ctx);

    contractionBlocks.push_back(
        to_contractions(expression, currentResult.resultTensor, *m_ctx));

    if (currentResult.importResultTensor) {
      m_importedTensors.insert(currentResult.resultTensor);
    } else {
      m_createdTensors.insert(currentResult.resultTensor);
    }

    // If we encounter a tensor in an expression that we have not yet seen
    // before, it must be an imported tensor (otherwise the expression would be
    // invalid).
    // Additionally, all result tensors that we find in the produced
    // contractions that is not imported must be created in order for the
    // expression to be valid.
    for (const Contraction &currentContraction : contractionBlocks.back()) {
      if (m_createdTensors.find(currentContraction.lhs) ==
          m_createdTensors.end()) {
        m_importedTensors.insert(currentContraction.lhs);
      }
      if (currentContraction.rhs.has_value()) {
        if (m_createdTensors.find(currentContraction.rhs.value()) ==
            m_createdTensors.end()) {
          m_importedTensors.insert(currentContraction.rhs.value());
        }
      }
      if (m_importedTensors.find(currentContraction.result) ==
          m_importedTensors.end()) {
        m_createdTensors.insert(currentContraction.result);
      }
    }
  }

  m_codes.push_back(CodeSection{block.name, std::move(contractionBlocks)});
}

struct IndexComponents {
  IndexSpace space;
  std::size_t id;
};

IndexComponents decomposeIndex(const Index &idx) {
  // The labels are of the form <letter>_<number> and we want <number>
  // Note that SeQuant uses 1-based indexing, but we want 0-based
  int num =
      std::stoi(std::wstring(idx.label().substr(idx.label().find('_') + 1))) -
      1;
  assert(num >= 0 && num <= 6);

  return {idx.space(), static_cast<std::size_t>(num)};
}

std::map<IndexSpace, std::set<std::size_t>> indicesBySpace(
    const std::set<Index> &indices) {
  std::map<IndexSpace, std::set<std::size_t>> indexMap;

  for (const Index &current : indices) {
    IndexComponents components = decomposeIndex(current);

    indexMap[components.space].insert(components.id);
  }

  return indexMap;
}

std::wstring to_itf(const Tensor &tensor, const Context &ctx,
                    bool includeIndexing = true) {
  std::wstring tags;
  std::wstring indices;

  // Note that it is important to iterate over the auxiliary indices first as
  // ITF expects those the be listed first
  for (const Index &current :
       ranges::views::concat(tensor.aux(), tensor.bra(), tensor.ket())) {
    IndexComponents components = decomposeIndex(current);

    assert(components.id <= 7);
    tags += ctx.get_tag(current.space());
    assert(ctx.get_base_label(current.space()).size() == 1);
    indices += static_cast<wchar_t>(ctx.get_base_label(current.space())[0] +
                                    components.id);
  }

  return std::wstring(tensor.label()) + (tags.empty() ? L"" : L":" + tags) +
         (includeIndexing ? L"[" + indices + L"]" : L"");
}

std::wstring ITFGenerator::generate() const {
  std::wstring itf =
      L"// This ITF algo file has been generated via SeQuant's ITF export\n\n";

  itf += L"---- decl\n";

  // Index declarations
  std::map<IndexSpace, std::set<std::size_t>> indexGroups =
      indicesBySpace(m_encounteredIndices);
  for (auto iter = indexGroups.begin(); iter != indexGroups.end(); ++iter) {
    const IndexSpace &space = iter->first;

    std::wstring label = m_ctx->get_base_label(space);
    if (label.size() != 1) {
      throw std::runtime_error(
          "Base labels are restricted to a size of 1 (at the moment)");
    }
    wchar_t baseLabel = label[0];

    std::wstring spaceName = m_ctx->get_name(space);
    std::wstring spaceTag = m_ctx->get_tag(space);

    itf += L"index-space: ";
    for (std::size_t i : iter->second) {
      assert(i <= 7);
      itf += static_cast<wchar_t>(baseLabel + i);
    }
    itf += L", " + spaceName + L", " + spaceTag + L"\n";
  }

  itf += L"\n";

  // Tensor declarations
  for (const Tensor &current : m_importedTensors) {
    if (current.indices().size() == 0 && current.label() == L"One") {
      // The One[] tensor exists implicitly
      continue;
    }

    itf += L"tensor: " + to_itf(current, *m_ctx) + L", " +
           to_itf(current, *m_ctx, false) + L"\n";
  }
  itf += L"\n";
  for (const Tensor &current : m_createdTensors) {
    if (current.indices().size() == 0 && current.label() == L"One") {
      // The One[] tensor exists implicitly
      continue;
    }

    itf += L"tensor: " + to_itf(current, *m_ctx);
    if (current.indices().size() > 0) {
      itf += +L", !Create{type:disk}\n";
    } else {
      itf += +L", !Create{type:scalar}\n";
    }
  }
  itf += L"\n\n";

  // Actual code
  for (const CodeSection &currentSection : m_codes) {
    itf += L"---- code(\"" + currentSection.name + L"\")\n";

    std::set<Tensor, TensorBlockCompare> allocatedTensors;

    for (const std::vector<Contraction> &currentBlock :
         currentSection.contractionBlocks) {
      for (const Contraction &currentContraction : currentBlock) {
        if (currentContraction.factor == 0) {
          if (currentContraction.lhs.label() == L"One" &&
              !currentContraction.rhs.has_value()) {
            // This is likely the only contraction belonging to the given result
            // meaning that we simply want to explicitly set it to zero
            itf +=
                L"alloc " + to_itf(currentContraction.result, *m_ctx) + L"\n";
            itf +=
                L"store " + to_itf(currentContraction.result, *m_ctx) + L"\n";
          }

          continue;
        }

        // For now we'll do a really silly contribution-by-contribution
        // load-process-store strategy
        if (allocatedTensors.find(currentContraction.result) ==
            allocatedTensors.end()) {
          itf += L"alloc " + to_itf(currentContraction.result, *m_ctx) + L"\n";
          allocatedTensors.insert(currentContraction.result);
        } else {
          itf += L"load " + to_itf(currentContraction.result, *m_ctx) + L"\n";
        }
        itf += L"load " + to_itf(currentContraction.lhs, *m_ctx) + L"\n";
        if (currentContraction.rhs.has_value()) {
          itf +=
              L"load " + to_itf(currentContraction.rhs.value(), *m_ctx) + L"\n";
        }

        itf += L"." + to_itf(currentContraction.result, *m_ctx) + L" ";
        int sign = currentContraction.factor < 0 ? -1 : 1;

        itf += (sign < 0 ? L"-= " : L"+= ");
        if (currentContraction.factor * sign != 1) {
          itf += to_wstring(currentContraction.factor * sign) + L" * ";
        }
        itf += to_itf(currentContraction.lhs, *m_ctx) +
               (currentContraction.rhs.has_value()
                    ? L" " + to_itf(currentContraction.rhs.value(), *m_ctx)
                    : L"") +
               L"\n";

        if (currentContraction.rhs.has_value()) {
          itf +=
              L"drop " + to_itf(currentContraction.rhs.value(), *m_ctx) + L"\n";
        }
        itf += L"drop " + to_itf(currentContraction.lhs, *m_ctx) + L"\n";
        itf += L"store " + to_itf(currentContraction.result, *m_ctx) + L"\n";
      }

      itf += L"\n";
    }
  }

  itf += L"\n---- end\n";

  return itf;
}

}  // namespace detail

}  // namespace itf

}  // namespace sequant
