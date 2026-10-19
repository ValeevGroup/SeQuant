#ifndef SEQUANT_EXTERNAL_INTERFACE_TESTS_TEST_EXTINT_HPP
#define SEQUANT_EXTERNAL_INTERFACE_TESTS_TEST_EXTINT_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "to_export_tree_step.hpp"

#include <SeQuant/core/io/shorthands.hpp>

#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint::test_support {

/// Redirects std::cout and std::wcout into throwaway buffers for its lifetime
/// so that the pipeline's progress chatter and the output step's dumps do not
/// clutter the test report. Restores the original buffers on destruction.
class ScopedStreamSilencer {
 public:
  ScopedStreamSilencer()
      : cout_(std::cout.rdbuf(sink_.rdbuf())),
        wcout_(std::wcout.rdbuf(wsink_.rdbuf())) {}

  ~ScopedStreamSilencer() {
    std::cout.rdbuf(cout_);
    std::wcout.rdbuf(wcout_);
  }

  ScopedStreamSilencer(const ScopedStreamSilencer &) = delete;
  ScopedStreamSilencer &operator=(const ScopedStreamSilencer &) = delete;

 private:
  std::ostringstream sink_;
  std::wostringstream wsink_;
  std::streambuf *cout_;
  std::wstreambuf *wcout_;
};

/// Absolute path to a file under tests/inputs/, independent of the working
/// directory the test binary is launched from.
inline std::string input_path(std::string_view name) {
  return std::string(SEQUANT_EXTINT_TESTS_SOURCE_DIR) + "/inputs/" +
         std::string(name);
}

/// Builds an ExpressionData from serialized ResultExpr definitions.
inline ExpressionData expr_data(
    std::initializer_list<std::wstring_view> definitions) {
  ExpressionData data;
  data.expressions.reserve(definitions.size());
  for (std::wstring_view def : definitions) {
    data.expressions.push_back(deserialize<ResultExpr>(def));
  }
  return data;
}

/// Converts serialized ResultExpr definitions into ExportTreeData (one entry
/// per definition) via the to_export_tree step.
inline ExportTreeData tree_data(
    std::initializer_list<std::wstring_view> definitions) {
  ExecutionContext ctx;
  ctx.set_data("in", 0, ProcessingData{expr_data(definitions)});
  ToExportTreeStep step;
  step.run("tree", ctx, {"in.0"});
  return convert_data<ExportTreeData>(
      ctx.get_data("tree.0").front().data.get());
}

}  // namespace sequant::util::extint::test_support

#endif  // SEQUANT_EXTERNAL_INTERFACE_TESTS_TEST_EXTINT_HPP
