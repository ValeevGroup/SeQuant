#include <SeQuant/core/optimize.hpp>
#include "algorithm.hpp"
#include "utils.hpp"

#include "cmdline.hpp"

#include <clipp.h>
#include <fmt/format.h>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <fstream>
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

void load_sequant() {
  using namespace sequant;
  set_locale();
  detail::OpIdRegistrar op_id_registrar;
  mbpt::set_default_convention();
}

int main(int argc, char* argv[]) {
  using std::cerr;
  using std::cout;
  using std::endl;
  using namespace sequant;

  // initialization
  load_sequant();
  initialize_db();
  //

  auto doc_fmt = clipp::doc_formatting{};
  doc_fmt.first_column(2);
  doc_fmt.indent_size(2);
  doc_fmt.last_column(120);

  container::vector<std::string> all_args;
  all_args.reserve(argc - 1);
  for (auto i = 1; i < argc; ++i) all_args.emplace_back(argv[i]);

  // ///////////////////////////////////////
  // Expression actions
  // ///////////////////////////////////////
  auto graph_cmd = cmdl::Graph{};
  auto latex_cmd = cmdl::Latex{};
  auto expr_cmd = cmdl::Expr{};
  container::vector<cmdl::Command*> expr_actions{&expr_cmd, &latex_cmd,
                                                 &graph_cmd};
  //

  // ///////////////////////////////////////
  // Expression sources
  // ///////////////////////////////////////
  auto cck_cmd = cmdl::CCk{};
  auto read_cmd = cmdl::Read{};
  container::vector<cmdl::Command*> expr_sources{&cck_cmd, &read_cmd};
  //

  // when asked for a command help, show help and exit
  for (auto ptr : ranges::views::concat(expr_actions, expr_sources))
    ptr->parse(argc, argv);
  // done single command help

  bool show_help;
  std::string outfile;
  container::vector<std::string> action_args, source_args;

  clipp::group cli =
      (clipp::option("-h", "--help")
           .set(show_help)
           .doc("Show this screen\n"
                "run 'sqntly COMMAND --help' for more"),
       clipp::option("--out").doc("Write output to this file") &
           clipp::value("FILE", outfile),
       (cmdl::one_of_cmds(expr_actions)
                .doc("COMMAND: Action to perform on a sequant::Sum") &
            clipp::opt_values("ARGS", action_args),
        cmdl::one_of_cmds(expr_sources)
                .doc("COMMAND: A source of sequant::Sum") &
            clipp::opt_values("ARGS", source_args)));

  assert(cli.flags_are_prefix_free());

  auto pr = clipp::parse(argc, argv, cli);
  if (show_help) {
    cout << clipp::make_man_page(cli, "sqntly", doc_fmt) << endl;
    return 0;
  }

  clipp::parsing_result pr_action, pr_source;
  for (auto ptr : expr_actions)
    if (ptr->set()) pr_action = ptr->parse(action_args);
  for (auto ptr : expr_sources)
    if (ptr->set()) pr_source = ptr->parse(source_args);
  if (!(pr && pr_action && pr_source)) {
    cout << "Usage:\n";
    cout << clipp::usage_lines(cli, "sqntly", doc_fmt) << endl;
    return 1;
  }

  ExprPtr expr;
  for (auto ptr : expr_sources)
    if (ptr->set()) expr = dynamic_cast<cmdl::ExprSource const&>(*ptr).expr();

  for (auto ptr : expr_actions)
    if (ptr->set()) {
      std::ofstream ofs{outfile};
      dynamic_cast<cmdl::ExprAction const&>(*ptr).run(
          outfile.empty() ? std::cout : ofs, expr);
    }
  return 0;
}