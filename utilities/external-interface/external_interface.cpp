#include "format_support.hpp"
#include "processing.hpp"
#include "utils.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generation_optimizer.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>  // for remove_tensor

#include <CLI/CLI.hpp>

#include <nlohmann/json.hpp>

#include <spdlog/spdlog.h>

#include <boost/algorithm/string.hpp>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <variant>

using nlohmann::json;
using namespace sequant;

template <>
struct std::hash<Tensor> {
  std::size_t operator()(const Tensor &tensor) const {
    return tensor.hash_value();
  }
};

class ItfExportContext : public ItfContext {
 public:
  ItfExportContext(const IndexSpaceMeta &meta) : m_meta(&meta) {}

  std::string get_tag(const IndexSpace &space) const override {
    assert(m_meta);
    return m_meta->getTag(space);
  }

  std::string get_name(const IndexSpace &space) const override {
    assert(m_meta);
    return m_meta->getName(space);
  }

 private:
  const IndexSpaceMeta *m_meta = nullptr;
};

ProcessingOptions extractProcessingOptions(
    const json &details, const ProcessingOptions &defaults = {}) {
  ProcessingOptions options = defaults;

  if (details.contains("density_fitting")) {
    options.density_fitting = details.at("density_fitting").get<bool>();
  }

  if (details.contains("spintracing")) {
    const std::string spintrace = details.at("spintracing").get<std::string>();

    if (spintrace == "none") {
      options.spintrace = SpinTracing::None;
    } else if (spintrace == "closed_shell") {
      options.spintrace = SpinTracing::ClosedShell;
    } else if (spintrace == "rigorous") {
      options.spintrace = SpinTracing::Rigorous;
    } else {
      throw std::runtime_error("Invalid spintracing option '" + spintrace +
                               "'");
    }
  }

  if (details.contains("projection")) {
    const std::string projection = details.at("projection").get<std::string>();

    if (projection == "primitive") {
      options.transform = ProjectionTransformation::None;
    } else if (projection == "biorthogonal") {
      options.transform = ProjectionTransformation::Biorthogonal;
    } else {
      throw std::runtime_error("Invalid projection option '" + projection +
                               "'");
    }
  }

  if (details.contains("optimize")) {
    options.factorize_to_binary = details.at("optimize").get<bool>();
  }

  if (details.contains("expand_symmetrizer")) {
    options.expand_symmetrizer = details.at("expand_symmetrizer").get<bool>();
  }

  if (details.contains("term_by_term")) {
    options.term_by_term = details.at("term_by_term").get<bool>();
  }

  if (details.contains("subexpression_elimination")) {
    options.subexpression_elimination =
        details.at("subexpression_elimination").get<bool>();
  }

  return options;
}

ExportNode<> prepareForExport(const ResultExpr &result,
                              const ItfGenerator<ItfExportContext> &generator,
                              ItfExportContext &ctx, bool importResult,
                              bool createResult) {
  ExportNode<> tree = to_export_tree(result);

  if (result.produces_tensor()) {
    assert(result.has_label());
    Tensor result_tensor = result.result_as_tensor();
    ctx.rewrite(result_tensor);
    if (importResult) {
      ctx.set_import_name(result_tensor,
                          generator.get_name(result_tensor, ctx));
    }
    if (createResult) {
      ctx.setLoadStrategy(result_tensor, LoadStrategy::Create, tree->id());
    } else {
      ctx.setLoadStrategy(result_tensor, LoadStrategy::Load, tree->id());
    }
  } else {
    if (importResult) {
      assert(result.has_label());
      ctx.set_import_name(result.result_as_variable(), toUtf8(result.label()));
    }
    if (createResult) {
      ctx.setLoadStrategy(result.result_as_variable(), LoadStrategy::Create,
                          tree->id());
    } else {
      ctx.setLoadStrategy(result.result_as_variable(), LoadStrategy::Load,
                          tree->id());
    }
  }

  return tree;
}

std::vector<ResultExpr> splitContributions(const ResultExpr &result) {
  if (!result.expression()->is<Sum>()) {
    return {result};
  }

  assert(result.expression()->is<Sum>());

  Tensor resultTensor(result.label(), bra(result.bra()), ket(result.ket()),
                      aux(result.aux()), result.symmetry(),
                      result.braket_symmetry(), result.particle_symmetry());

  std::vector<ResultExpr> contributions;
  contributions.reserve(result.expression()->size());

  for (const ExprPtr &subExpr : result.expression()) {
    contributions.emplace_back(resultTensor, subExpr);
  }

  return contributions;
}

void generateITF(const json &blocks, std::string_view out_file,
                 const ProcessingOptions &defaults,
                 const IndexSpaceMeta &spaceMeta) {
  ItfExportContext context(spaceMeta);
  // We assume index IDs start at 1
  context.set_index_id_offset(1);

  ItfGenerator<ItfExportContext> itfgen;
  GenerationOptimizer<ItfGenerator<ItfExportContext>> generator(itfgen);

  container::svector<ExpressionGroup<>> groups;
  groups.reserve(blocks.size());

  for (const json &current_block : blocks) {
    const std::string block_name = current_block.at("name");

    spdlog::debug("Processing ITF code block '{}'", block_name);

    container::svector<ExportNode<>> results;

    for (const json &current_result : current_block.at("results")) {
      const std::string result_name = current_result.at("name");
      const std::string input_file = current_result.at("equation_file");

      spdlog::debug("Processing equations from '{}' to result '{}'", input_file,
                    result_name);

      if (!std::filesystem::exists(input_file)) {
        throw std::runtime_error("Specified input file '" + input_file +
                                 "' does not exist");
      }

      // Read input file
      std::ifstream in(input_file);
      const std::string input(std::istreambuf_iterator<char>(in), {});

      sequant::ResultExpr result =
          sequant::parse_result_expr(toUtf16(input), Symmetry::antisymm);

      if (current_result.contains("name")) {
        result.set_label(toUtf16(current_result.at("name").get<std::string>()));
      }

      spdlog::debug("Initial equation is:\n{}", result);

      ProcessingOptions options =
          extractProcessingOptions(current_result, defaults);

      std::vector<ResultExpr> resultParts =
          options.term_by_term ? splitContributions(result)
                               : std::vector<ResultExpr>{result};

      std::unordered_set<Tensor> tensorsToSymmetrize;

      std::set<std::variant<Tensor, Variable>> createdResults;
      for (const ResultExpr &contribution : resultParts) {
        if (resultParts.size() > 1) {
          spdlog::debug("Current contribution:\n{}", contribution);
        }

        for (ResultExpr &current :
             postProcess(contribution, spaceMeta, options)) {
          spdlog::debug("Fully processed equation is:\n{}", current);

          if (*current.expression() == Constant(0)) {
            continue;
          }

          const bool createResult = [&]() {
            // We only want to create a given result once to not overwrite
            // previous contributions
            if (contribution.produces_tensor()) {
              if (createdResults.find(current.result_as_tensor()) ==
                  createdResults.end()) {
                createdResults.insert(current.result_as_tensor());
                return true;
              }
            } else if (createdResults.find(current.result_as_variable()) ==
                       createdResults.end()) {
              createdResults.insert(current.result_as_variable());
              return true;
            }

            return false;
          }();

          if (needsSymmetrization(current.expression())) {
            std::optional<ExprPtr> symmetrizer =
                pop_tensor(current.expression(), L"S");
            assert(symmetrizer.has_value());

            Tensor resultTensor(current.label(), bra(current.bra()),
                                ket(current.ket()), aux(current.aux()),
                                current.symmetry(), current.braket_symmetry(),
                                current.particle_symmetry());
            tensorsToSymmetrize.insert(resultTensor);

            current.set_label(current.label() + L"u");

            spdlog::debug("After popping S tensor:\n{}", current);

            results.push_back(prepareForExport(current, itfgen, context, false,
                                               createResult));
          } else {
            results.push_back(prepareForExport(
                current, itfgen, context, current_result.value("import", true),
                createResult));
          }
        }
      }

      for (const Tensor &current : tensorsToSymmetrize) {
        ExprPtr symmetrization = generateResultSymmetrization(
            current, std::wstring(current.label()) + L"u");

        ResultExpr symmetrizedResult(current, std::move(symmetrization));

        spdlog::debug("Result symmetrization via\n{}", symmetrizedResult);

        results.push_back(prepareForExport(symmetrizedResult, itfgen, context,
                                           current_result.value("import", true),
                                           true));
      }
    }

    if (current_block.value("subexpression_elimination",
                            defaults.subexpression_elimination)) {
      eliminate_common_subexpressions(results, [](const auto &expr) {
        // Note: the lambda is needed to make the callable usable for
        // ExprPtr as well as ResultExpr objects
        return to_export_tree(expr);
      });
    }

    groups.emplace_back(std::move(results),
                        current_block.at("name").get<std::string>());
  }

  export_groups(groups, generator, context);
  std::string itfCode = generator.get_generated_code();

  std::ofstream output(out_file.data());
  output << itfCode;
}

void generateCode(const json &details, const IndexSpaceMeta &spaceMeta) {
  const std::string format = details.at("output_format");
  const std::string out_path = details.at("output_path");

  ProcessingOptions defaultOptions;
  if (details.contains("default_options")) {
    defaultOptions = extractProcessingOptions(details.at("default_options"));
  }

  if (boost::iequals(format, "itf")) {
    generateITF(details.at("code_blocks"), out_path, defaultOptions, spaceMeta);
  } else {
    throw std::runtime_error("Unknown code generation target format '" +
                             std::string(format) + "'");
  }
}

void registerIndexSpaces(const json &spaces, IndexSpaceMeta &meta) {
  IndexSpaceRegistry &registry =
      *get_default_context().mutable_index_space_registry();

  std::vector<std::pair<std::wstring, IndexSpaceMeta::Entry>> spaceList;
  spaceList.reserve(spaces.size());

  for (std::size_t i = 0; i < spaces.size(); ++i) {
    const json &current = spaces.at(i);
    const IndexSpace::Type type(1 << i);
    const int size = current.at("size");

    if (size <= 0) {
      throw std::runtime_error("Index space sizes must be > 0");
    }

    IndexSpaceMeta::Entry entry;
    entry.name = current.at("name").get<std::string>();
    entry.tag = current.at("tag").get<std::string>();

    std::wstring label = toUtf16(current.at("label").get<std::string>());
    registry.add(label, type, size);

    spdlog::debug(
        "Registered index space '{}' with label '{}', tag '{}' and size {}",
        entry.name, toUtf8(label), entry.tag, size);

    spaceList.push_back(std::make_pair(std::move(label), std::move(entry)));
  }

  mbpt::add_fermi_spin(registry);

  for (auto &[label, entry] : spaceList) {
    meta.registerSpace(Index(label + L"_1").space(), std::move(entry));
  }
}

void process(const json &driver, IndexSpaceMeta &spaceMeta) {
  if (!driver.contains("index_spaces")) {
    throw std::runtime_error("Missing index_spaces definition");
  }

  registerIndexSpaces(driver.at("index_spaces"), spaceMeta);

  if (driver.contains("code_generation")) {
    const json &details = driver.at("code_generation");

    generateCode(details, spaceMeta);
  }
}

void generalSetup() {
  TensorCanonicalizer::set_cardinal_tensor_labels(
      mbpt::cardinal_tensor_labels());
}

int main(int argc, char **argv) {
  set_locale();
  set_default_context(Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  generalSetup();

  CLI::App app(
      "Interface for reading in equations generated outside of SeQuant");
  argv = app.ensure_utf8(argv);

  std::filesystem::path driver;
  app.add_option("--driver", driver,
                 "Path to the JSON file used to drive the processing")
      ->required();
  bool verbose = false;
  app.add_flag("--verbose", verbose, "Whether to enable verbose output");

  CLI11_PARSE(app, argc, argv);

  if (!std::filesystem::exists(driver)) {
    throw std::runtime_error("Specified driver file '" + driver.string() +
                             "' does not exist");
  } else if (!std::filesystem::is_regular_file(driver)) {
    throw std::runtime_error("Specified driver file '" + driver.string() +
                             "' is not a file");
  }

  // Change directory to where the driver file is located so that all relative
  // paths specified in it resolve to be relative to the driver file.
  // Before that, we have to get the absolute driver path though (or else the
  // relative path will no longer work)
  driver = std::filesystem::absolute(driver);
  std::filesystem::current_path(driver.parent_path());

  if (verbose) {
    spdlog::set_level(spdlog::level::debug);
  }

  IndexSpaceMeta spaceMeta;

  try {
    std::ifstream in(driver);
    json driver_info;
    driver_info =
        json::parse(in, /*callback*/ nullptr, /*allow_exceptions*/ true,
                    /*skip_comments*/ true);

    process(driver_info, spaceMeta);
  } catch (const std::exception &e) {
    spdlog::error("Unexpected error: {}", e.what());
    return 1;
  }

  return 0;
}
