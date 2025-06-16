//
// Created by Robert Adam on 2023-09-04
//

#ifndef SEQUANT_CORE_EXPORT_ITF_HPP
#define SEQUANT_CORE_EXPORT_ITF_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <functional>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace sequant {

namespace itf {

/// A result consists of a single result tensor along with an expression tree
/// that contains the expressions needed to build the aforementioned result
/// tensor
struct Result {
  ExprPtr expression;
  Tensor resultTensor;
  bool importResultTensor;

  Result(ExprPtr expression, Tensor resultTensor,
         bool importResultTensor = true);
  Result(ExprPtr expression, bool importResultTensor = false);
};

/// A code block groups one or multiple Results together. Upon ITF code
/// generation a CodeBlock object will be mapped to a single '---- code("xyz")'
/// block in ITF. The individual results will be computed serially inside the
/// code block and the computation of individual results will be separated in
/// different "paragraphs" in the generated code.
struct CodeBlock {
  std::wstring name;
  std::vector<Result> results;

  CodeBlock(std::wstring blockName, Result result);
  CodeBlock(std::wstring blockName, std::vector<Result> results);
};

class Context {
 public:
  virtual ~Context();

  virtual int compare(const Index &lhs, const Index &rhs) const = 0;
  virtual std::wstring get_base_label(const IndexSpace &space) const = 0;
  virtual std::wstring get_tag(const IndexSpace &space) const = 0;
  virtual std::wstring get_name(const IndexSpace &space) const = 0;
};

namespace detail {

/// Comparator that identifies Tensors only by their "block", which is defined
/// by its name, the amount of its indices as well as the space these indices
/// belong to. Note that it explicitly does not depend on the explicit index
/// labelling.
struct TensorBlockCompare {
  bool operator()(const Tensor &lhs, const Tensor &rhs) const;
};

/// Replaces one- and two-electron integrals in the given expression with
/// versions that use the index ordering and tensor naming as expected in ITF
void remap_integrals(ExprPtr &expr, const Context &ctx,
                     std::wstring_view oneElectronIntegralName = L"f",
                     std::wstring_view twoElectronIntegralName = L"g",
                     std::wstring_view dfTensorName = L"DF");

/// Represents a single contraction where the contraction of lhs with rhs,
/// multiplied by the given factor contributes to the given result. If rhs is
/// not given, this implies that the contribution is only lhs scaled with the
/// provided factor.
struct Contraction {
  rational factor;

  Tensor result;
  Tensor lhs;
  std::optional<Tensor> rhs;
};

/// This is essentially a low-level representation of a CodeBlock object, where
/// all expressions have been broken down into binary contractions, which leads
/// to the creation of intermediate result tensors.
struct CodeSection {
  std::wstring name;
  std::vector<std::vector<Contraction>> contractionBlocks;
};

/// This is stateful object that is used to implement the actual ITF code
/// generating capabilities
class ITFGenerator {
 public:
  ITFGenerator(const itf::Context &ctx);

  void addBlock(const itf::CodeBlock &block);

  std::wstring generate() const;

 private:
  std::set<Index> m_encounteredIndices;
  std::set<Tensor, TensorBlockCompare> m_importedTensors;
  std::set<Tensor, TensorBlockCompare> m_createdTensors;
  std::vector<CodeSection> m_codes;
  const Context *m_ctx;
};

}  // namespace detail

}  // namespace itf

/// Translates the given ITF CodeBlock to executable ITF code
std::wstring to_itf(const itf::CodeBlock &block, const itf::Context &ctx);

/// Translates the given collection/range of ITF CodeBlocks or Results to
/// executable ITF code
template <typename Container,
          typename = std::enable_if_t<!std::is_same_v<
              std::remove_const_t<std::remove_reference_t<Container>>,
              itf::CodeBlock>>>
std::wstring to_itf(Container &&container, const itf::Context &ctx) {
  using ContainerType = std::remove_const_t<std::remove_reference_t<Container>>;
  static_assert(
      std::is_same_v<typename ContainerType::value_type, itf::CodeBlock> ||
      std::is_same_v<typename ContainerType::value_type, ExprPtr>);
  itf::detail::ITFGenerator generator(ctx);

  if constexpr (std::is_same_v<typename ContainerType::value_type,
                               itf::CodeBlock>) {
    for (const itf::CodeBlock &current : container) {
      generator.addBlock(current);
    }
  } else {
    static_assert(
        std::is_same_v<typename ContainerType::value_type, itf::Result>,
        "Container::value_type must either be itf::CodeBlock or itf::Result");

    itf::CodeBlock block(L"Generate_Results",
                         std::forward<Container>(container));

    generator.addBlock(block);
  }

  return generator.generate();
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_ITF_HPP
