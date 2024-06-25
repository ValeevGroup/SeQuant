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

std::wstring to_itf(const itf::CodeBlock &block) {
  itf::detail::ITFGenerator generator;
  generator.addBlock(block);
  return generator.generate();
}

namespace itf {

struct IndexTypeComparer {
  bool operator()(const IndexSpace::Type &lhs,
                  const IndexSpace::Type &rhs) const {
    assert(get_default_context().index_space_registry()->retrieve("i").type() <
           get_default_context().index_space_registry()->retrieve("a").type());
    return lhs < rhs;
  }
};

Result::Result(ExprPtr expression, Tensor resultTensor, bool importResultTensor)
    : expression(std::move(expression)),
      resultTensor(std::move(resultTensor)),
      importResultTensor(importResultTensor) {}

Tensor generateResultTensor(ExprPtr expr) {
  // This (same as ITF itself) assumes that all indices that are repeated in the
  // expression are also contracted (summed) over and therefore don't contribute
  // to the result of the expression
  IndexGroups externals = get_unique_indices(expr);

  return Tensor(L"Result", std::move(externals.bra), std::move(externals.ket),
                std::move(externals.aux));
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

bool TensorBlockCompare::operator()(const Tensor &lhs,
                                    const Tensor &rhs) const {
  if (lhs.label() != rhs.label()) {
    return lhs.label() < rhs.label();
  }
  if (lhs.braket().size() != rhs.braket().size()) {
    return lhs.braket().size() < rhs.braket().size();
  }
  auto lhsBraket = lhs.braket();
  auto rhsBraket = rhs.braket();

  for (std::size_t i = 0; i < lhsBraket.size(); ++i) {
    if (lhsBraket.at(i).space() != rhsBraket.at(i).space()) {
      return lhsBraket.at(i).space() < rhsBraket.at(i).space();
    }
  }

  return false;
}

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor);

std::vector<Contraction> to_contractions(const Product &product,
                                         const Tensor &resultTensor) {
  static std::size_t intermediateCounter = 1;

  if (product.factors().size() == 1) {
    assert(product.factor(0).is<Tensor>());
    assert(product.scalar().imag() == 0);

    return {Contraction{product.scalar().real(),
                        resultTensor,
                        product.factor(0).as<Tensor>(),
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
      // index spaces is the canonical one within ITF (largest space leftmost).
      // This is possible, because the index ordering of intermediates is
      // arbitrary
      std::vector<Index> intermediateIndices;
      intermediateIndices.reserve(intermediateIndexGroups.bra.size() +
                                  intermediateIndexGroups.ket.size());
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
                [](const Index &lhs, const Index &rhs) {
                  IndexTypeComparer cmp;
                  if (cmp(lhs.space().type(), rhs.space().type())) {
                    return false;
                  } else if (cmp(rhs.space().type(), lhs.space().type())) {
                    return true;
                  } else {
                    // Indices are of same space
                    return lhs < rhs;
                  }
                });

      std::array<wchar_t, 64> intermediateName;
      swprintf(intermediateName.data(), intermediateName.size(), L"INTER%06u",
               intermediateCounter++);

      // There is no notion of bra and ket for intermediates, so we dump all
      // indices in the bra for now
      Tensor intermediate(intermediateName.data(),
                          std::move(intermediateIndices), std::vector<Index>{},
                          std::vector<Index>{});

      std::vector<Contraction> intermediateContractions =
          to_contractions(factor, intermediate);
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
      lhsIntermediate == intermediates.end() ? product.factor(0).as<Tensor>()
                                             : lhsIntermediate->second,
      rhsIntermediate == intermediates.end() ? product.factor(1).as<Tensor>()
                                             : rhsIntermediate->second});

  return contractions;
}

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor) {
  std::wstring itfCode;

  if (expression.is<Constant>()) {
    throw std::invalid_argument("Can't transform constants into contractions");
  } else if (expression.is<Tensor>()) {
    return {Contraction{1, resultTensor, expression.as<Tensor>(), {}}};
  } else if (expression.is<Product>()) {
    // Separate into binary contractions
    return to_contractions(expression.as<Product>(), resultTensor);
  } else if (expression.is<Sum>()) {
    // Process each summand
    std::vector<Contraction> contractions;

    for (const ExprPtr &summand : expression.as<Sum>().summands()) {
      std::vector<Contraction> currentContractions =
          to_contractions(summand, resultTensor);

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

void one_electron_integral_remapper(
    ExprPtr &expr, const std::wstring_view integralTensorLabel) {
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
  assert(tensor.auxiliary().empty());

  IndexTypeComparer cmp;

  // Use the bra-ket (hermitian) symmetry of the integrals to exchange creators
  // and annihilators such that the larger index space is on the left (in the
  // bra)
  if (cmp(braIndices[0].space().type(), ketIndices[0].space().type())) {
    std::swap(braIndices[0], ketIndices[0]);
  } else if (braIndices[0].space().type() == ketIndices[0].space().type() &&
             ketIndices[0] < braIndices[0]) {
    // Cosmetic exchange to arrive at a more canonical index ordering
    std::swap(braIndices[0], ketIndices[0]);
  }

  expr = ex<Tensor>(tensor.label(), std::move(braIndices),
                    std::move(ketIndices), tensor.auxiliary());
}

template <typename Container>
bool isExceptionalJ(const Container &braIndices, const Container &ketIndices) {
  assert(braIndices.size() == 2);
  assert(ketIndices.size() == 2);
  // integrals with 3 external (virtual) indices ought to be converted to
  // J-integrals
  return braIndices[0].space().type() == get_default_context()
                                             .index_space_registry()
                                             ->retrieve("a")
                                             .type() &&
         braIndices[1].space().type() == get_default_context()
                                             .index_space_registry()
                                             ->retrieve("a")
                                             .type() &&
         ketIndices[0].space().type() == get_default_context()
                                             .index_space_registry()
                                             ->retrieve("a")
                                             .type() &&
         ketIndices[1].space().type() !=
             get_default_context().index_space_registry()->retrieve("a").type();
}

void two_electron_integral_remapper(
    ExprPtr &expr, const std::wstring_view integralTensorLabel) {
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
  assert(tensor.auxiliary().empty());

  IndexTypeComparer cmp;

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
    if (cmp(braIndices.at(i).space().type(), ketIndices.at(i).space().type())) {
      // This bra index belongs to a smaller space than the ket index ->
      // swap them
      std::swap(braIndices[i], ketIndices[i]);
    }
  }

  // Step 1b: Particle-1,2-symmetry
  bool switchColumns = false;
  if (braIndices[0].space().type() != braIndices[1].space().type()) {
    switchColumns =
        cmp(braIndices[0].space().type(), braIndices[1].space().type());
  } else if (ketIndices[0].space().type() != ketIndices[1].space().type()) {
    switchColumns =
        cmp(ketIndices[0].space().type(), ketIndices[1].space().type());
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

  if (isExceptionalJ(braIndices, ketIndices) ||
      cmp(braIndices[1].space().type(), ketIndices[0].space().type())) {
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

  expr = ex<Tensor>(std::move(tensorLabel), std::move(braIndices),
                    std::move(ketIndices), tensor.auxiliary());
}

void integral_remapper(ExprPtr &expr, std::wstring_view oneElectronIntegralName,
                       std::wstring_view twoElectronIntegralName) {
  two_electron_integral_remapper(expr, twoElectronIntegralName);
  one_electron_integral_remapper(expr, oneElectronIntegralName);
}

void remap_integrals(ExprPtr &expr, std::wstring_view oneElectronIntegralName,
                     std::wstring_view twoElectronIntegralName) {
  auto remapper = std::bind(integral_remapper, std::placeholders::_1,
                            oneElectronIntegralName, twoElectronIntegralName);

  const bool visitedRoot = expr->visit(remapper, true);

  if (!visitedRoot) {
    remapper(expr);
  }
}

void ITFGenerator::addBlock(const itf::CodeBlock &block) {
  m_codes.reserve(m_codes.size() + block.results.size());

  std::vector<std::vector<Contraction>> contractionBlocks;

  for (const Result &currentResult : block.results) {
    ExprPtr expression = currentResult.expression;
    remap_integrals(expression);

    contractionBlocks.push_back(
        to_contractions(expression, currentResult.resultTensor));

    if (currentResult.importResultTensor) {
      m_importedTensors.insert(currentResult.resultTensor);
    } else {
      m_createdTensors.insert(currentResult.resultTensor);
    }

    // If we encounter a tensor in an expression that we have not yet seen
    // before, it must be an imported tensor (otherwise the expression would be
    // invalid)
    expression->visit(
        [this](const ExprPtr &expr) {
          if (expr.is<Tensor>()) {
            const Tensor &tensor = expr.as<Tensor>();
            if (m_createdTensors.find(tensor) == m_createdTensors.end()) {
              m_importedTensors.insert(tensor);
            }
            m_encounteredIndices.insert(tensor.braket().begin(),
                                        tensor.braket().end());
          }
        },
        true);

    // Now go through all result tensors of the contractions that we have
    // produced and add all new tensors to the set of created tensors
    for (const Contraction &currentContraction : contractionBlocks.back()) {
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

std::wstring to_itf(const Tensor &tensor, bool includeIndexing = true) {
  std::wstring tags;
  std::wstring indices;

  for (const Index &current : tensor.braket()) {
    IndexComponents components = decomposeIndex(current);

    assert(components.id <= 7);

    if (components.space.type() ==
        get_default_context().index_space_registry()->retrieve("i").type()) {
      tags += L"c";
      indices += static_cast<wchar_t>(L'i' + components.id);
    } else if (components.space.type() == get_default_context()
                                              .index_space_registry()
                                              ->retrieve("a")
                                              .type()) {
      tags += L"e";
      indices += static_cast<wchar_t>(L'a' + components.id);
    } else {
      throw std::runtime_error("Encountered unhandled index space type");
    }
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
    wchar_t baseLabel;
    std::wstring spaceLabel;
    std::wstring spaceTag;
    if (iter->first.type() ==
        get_default_context().index_space_registry()->retrieve("i").type()) {
      baseLabel = L'i';
      spaceLabel = L"Closed";
      spaceTag = L"c";
    } else if (iter->first.type() == get_default_context()
                                         .index_space_registry()
                                         ->retrieve("a")
                                         .type()) {
      baseLabel = L'a';
      spaceLabel = L"External";
      spaceTag = L"e";
    } else {
      throw std::runtime_error("Encountered unhandled index space type");
    }

    itf += L"index-space: ";
    for (std::size_t i : iter->second) {
      assert(i <= 7);
      itf += static_cast<wchar_t>(baseLabel + i);
    }
    itf += L", " + spaceLabel + L", " + spaceTag + L"\n";
  }

  itf += L"\n";

  // Tensor declarations
  for (const Tensor &current : m_importedTensors) {
    itf +=
        L"tensor: " + to_itf(current) + L", " + to_itf(current, false) + L"\n";
  }
  itf += L"\n";
  for (const Tensor &current : m_createdTensors) {
    itf += L"tensor: " + to_itf(current) + L", !Create{type:disk}\n";
  }
  itf += L"\n\n";

  // Actual code
  for (const CodeSection &currentSection : m_codes) {
    itf += L"---- code(\"" + currentSection.name + L"\")\n";

    std::set<Tensor, TensorBlockCompare> allocatedTensors;

    for (const std::vector<Contraction> &currentBlock :
         currentSection.contractionBlocks) {
      for (const Contraction &currentContraction : currentBlock) {
        // For now we'll do a really silly contribution-by-contribution
        // load-process-store strategy
        if (allocatedTensors.find(currentContraction.result) ==
            allocatedTensors.end()) {
          itf += L"alloc " + to_itf(currentContraction.result) + L"\n";
          allocatedTensors.insert(currentContraction.result);
        } else {
          itf += L"load " + to_itf(currentContraction.result) + L"\n";
        }
        itf += L"load " + to_itf(currentContraction.lhs) + L"\n";
        if (currentContraction.rhs.has_value()) {
          itf += L"load " + to_itf(currentContraction.rhs.value()) + L"\n";
        }

        itf += L"." + to_itf(currentContraction.result) + L" ";
        int sign = currentContraction.factor < 0 ? -1 : 1;

        itf += (sign < 0 ? L"-= " : L"+= ");
        if (currentContraction.factor * sign != 1) {
          itf += to_wstring(currentContraction.factor * sign) + L" * ";
        }
        itf += to_itf(currentContraction.lhs) +
               (currentContraction.rhs.has_value()
                    ? L" " + to_itf(currentContraction.rhs.value())
                    : L"") +
               L"\n";

        if (currentContraction.rhs.has_value()) {
          itf += L"drop " + to_itf(currentContraction.rhs.value()) + L"\n";
        }
        itf += L"drop " + to_itf(currentContraction.lhs) + L"\n";
        itf += L"store " + to_itf(currentContraction.result) + L"\n";
      }

      itf += L"\n";
    }

    itf += L"\n---- end\n";
  }

  return itf;
}

}  // namespace detail

}  // namespace itf

}  // namespace sequant
