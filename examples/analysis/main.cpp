#include "algorithm.hpp"
#include "cmdline.hpp"

#include <clipp.h>
#include <fmt/format.h>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <iostream>
#include <locale>

void set_locale() {
  // set global C++ locale (and C locale)
  std::locale target_locale{};
  try {  // use en_US.UTF-8, if supported
    target_locale = std::locale{"en_US.UTF-8"};
  } catch (std::exception&) {  // use default if en_US.UTF-8 not available
  }
  std::locale::global(target_locale);
  // set C++ streams locale to target
  std::ios_base::sync_with_stdio(false);
  std::cout.imbue(target_locale);
  std::cerr.imbue(target_locale);
  std::clog.imbue(target_locale);
  std::wcout.imbue(target_locale);
  std::wcerr.imbue(target_locale);
  std::wclog.imbue(target_locale);
  std::ios_base::sync_with_stdio(true);
}

void load_sequant_defaults() {
  using namespace sequant;
  set_locale();
  detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
}

int main(int argc, char* argv[]) {
  using std::cerr;
  using std::cout;
  using std::endl;
  using namespace sequant;

  // initialization
  load_sequant_defaults();
  initialize_db();
  //

  auto doc_fmt = clipp::doc_formatting{};
  doc_fmt.first_column(2);
  doc_fmt.indent_size(2);
  doc_fmt.last_column(120);
  doc_fmt.max_flags_per_param_in_usage(1);

  enum struct Command { Graph, Latex, Help, Invalid };

  auto cmnd = Command::Invalid;
  clipp::group sqntly;

  auto latex_opts = LatexOpts{};
  auto graph_opts = GraphOpts{};
  auto cck_opts = CCkOpts{};
  container::vector<std::string> unknowns;

  {
    clipp::group latex_group, graph_group, cck_group;

    latex_opts.add_options(latex_group);
    graph_opts.add_options(graph_group);

    auto latex_command = clipp::command("latex")
                             .set(cmnd, Command::Latex)
                             .doc("Generate latex expressions") &
                         latex_group;
    auto graph_command =
        clipp::command("graph").set(cmnd, Command::Graph).doc("Generate graphs")
        /* & graph_group //disabled for now */
        ;
    auto help_command = clipp::required("--help", "-h")
                            .set(cmnd, Command::Help)
                            .doc("Show help and exit");

    cck_opts.add_options(cck_group);
    cck_group.doc("Select a coupled-cluster equation");
    sqntly.push_back(help_command |
                     ((latex_command | graph_command), cck_group));
  }

  sqntly.push_back(clipp::any_other(unknowns));

  auto parse_result = clipp::parse(argc, argv, sqntly);
  if (cmnd == Command::Help) {
    cout << clipp::make_man_page(sqntly, "sqntly", doc_fmt) << endl;
    return 0;
  }

  if (!(parse_result && unknowns.empty())) {
    cerr << "Error parsing command line arguments\n";
    if (!unknowns.empty())
      cerr << fmt::format("Unknown arguments: {}\n", fmt::join(unknowns, ", "));

    cout << "Usage:\n"
         << clipp::usage_lines(sqntly, "sqntly", doc_fmt) << "\n\n"
         << "Run with '--help' for more." << endl;

    return 1;
  }

  assert(cmnd != Command::Invalid);
  cck_opts.sanity_check();

  auto expr = read_cck(cck_opts, database_dir());
  if (cmnd == Command::Graph) {
    for (auto const& vec : cse_graph(expr.as<Sum>()))
      cout << fmt::format("{}\n", fmt::join(vec, ","));
  } else if (cmnd == Command::Latex) {
    cout << to_string(to_latex_align(expr, latex_opts.lines_per_block(),
                                     latex_opts.terms_per_line()))
         << endl;
  }

  return 0;
}