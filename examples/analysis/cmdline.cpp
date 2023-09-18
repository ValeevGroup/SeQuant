#include "cmdline.hpp"
#include "algorithm.hpp"
#include "utils.hpp"

#include <fmt/format.h>
#include <SeQuant/core/optimize.hpp>

#include <fstream>
#include <iostream>

namespace sequant::cmdl {
void err(std::string_view msg = "Invalid command-line") {
  std::cerr << msg << std::endl;
}

void err(container::vector<std::string> const& args) {
  std::cerr << fmt::format(
      "Invalid command-line{}{}\n",
      args.empty() ? "" : ", unknown args: ", fmt::join(args, ", "));
}

template <typename T>
std::string doc_default(T& tgt, std::string_view doc) {
  return fmt::format("{}, default '{}'", doc, tgt);
}

// /////////////////////////////////////////////////////////////////////////////
// Command implementation
// /////////////////////////////////////////////////////////////////////////////
Command::Command(std::string_view name, std::string_view doc)
    : is_set_{}, name_{name}, doc_{doc} {}

std::string_view Command::name() const noexcept { return name_; }

std::string_view Command::doc() const noexcept { return doc_; }

void Command::parse(int argc, char** argv) {
  bool show_help;

  clipp::parse(argc, argv,
               clipp::command(name().data()) &
                   clipp::option("-h", "--help").set(show_help));
  if (show_help) {
    std::cout << clipp::make_man_page(group().empty()
                                          ? clipp::group(clipp::command(""))
                                          : group(),
                                      name().data())
                     .prepend_section("DESCRIPTION", doc().data())
              << std::endl;

    std::exit(0);
  }
}

clipp::parsing_result Command::parse(args_type const& args) {
  args_type unknown;
  auto ugroup = clipp::group(clipp::any_other(unknown));
  auto pr = clipp::parse(args, group().empty() ? ugroup : (group(), ugroup));
  if (!(pr && unknown.empty())) {
    err(unknown);
    std::exit(1);
  }
  return pr;
}

bool Command::set() const noexcept { return is_set_; }

void Command::set(bool yn) noexcept { is_set_ = yn; }

// /////////////////////////////////////////////////////////////////////////////
// ExprAction derived class implementations follow
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// Latex implementation
// /////////////////////////////////////////////////////////////////////////////
Latex::Latex() : ExprAction{"latex", "Generate LaTeX expression"} {}

clipp::group Latex::group() noexcept {
  return (clipp::option("--tpl").doc(
              doc_default(tpl_, "Number of terms per line")) &
              clipp::value("NUM", tpl_),
          clipp::option("--lpb").doc(
              doc_default(lpb_,
                          "Number of lines per align env.\n"
                          "(0 implies all lines in a single block)")) &
              clipp::value("NUM", lpb_));
}

void Latex::run(std::ostream& os, ExprPtr const& expr) const {
  os << to_string(to_latex_align(expr, lpb_, tpl_)) << std::endl;
}

// /////////////////////////////////////////////////////////////////////////////
// Graph implementation
// /////////////////////////////////////////////////////////////////////////////

Graph::Graph() : ExprAction{"graph", "Generate graph of terms"} {}

void Graph::run(std::ostream& os, ExprPtr const& expr) const {
  auto prev_locale = os.getloc();
  os.imbue(std::locale::classic());

  assert(expr.is<Sum>());
  if (format() != Format::AdjSet) {
    std::cerr << "Only 'set' format is supported so far\n";
    std::exit(1);
  }
  auto graph = adjacency_set(intermediates(expr.as<Sum>()));
  for (auto const& [k, pos] : graph) {
    os << k;
    for (auto p : pos) os << "," << p;
    os << '\n';
  }

  os.imbue(prev_locale);
}

clipp::group Graph::group() noexcept {
  clipp::group cmds;

  for (auto const& tpl : format_labels) {
    auto const& lbl = std::get<0>(tpl);
    auto const& doc = std::get<1>(tpl);
    auto val = std::get<2>(tpl);
    auto cm = clipp::group(
        clipp::command(lbl.data()).doc(doc.data()).call([this, val]() {
          format_ = val;
        }));
    cmds = cmds.empty() ? cm : cmds | cm;
  }

  std::string_view dflt;
  for (auto&& [l, _, f] : format_labels)
    if (format_ == f) dflt = l;

  return clipp::option("--fmt").doc(doc_default(dflt, "Graph formats")) & cmds;
}

Graph::Format Graph::format() const noexcept { return format_; }

// /////////////////////////////////////////////////////////////////////////////
// Expr implementation
// /////////////////////////////////////////////////////////////////////////////

Expr::Expr() : ExprAction{"expr", "Generate re-parsable text"} {}

clipp::group Expr::group() noexcept { return clipp::group{}; }

void Expr::run(std::ostream& os, ExprPtr const& expr) const {
  assert(expr.is<Sum>());
  write_sum(os, expr.as<Sum>());
}

// /////////////////////////////////////////////////////////////////////////////
// ExprSource derived classes implementations follow
// /////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// CCk implementation
// /////////////////////////////////////////////////////////////////////////////

CCk::CCk() : ExprSource{"CCk", "Generate coupled cluster equations"} {}

size_t CCk::excit() const noexcept { return excit_; }

size_t CCk::r() const noexcept { return r_; }

size_t CCk::nalpha() const noexcept { return nalpha_; }

bool CCk::optimize() const noexcept { return optimize_; }

OrbitalType CCk::orbital_type() const noexcept { return orbital_type_; }

SpinTraceAlgorithm CCk::st_algorithm() const noexcept { return st_algo_; }

clipp::group CCk::group() noexcept {
  using clipp::option;
  using clipp::required;
  using clipp::value;
  clipp::group g;
  g.push_back(required("--ex")
                  .doc(fmt::format("The highest level of excitation"
                                   " in coupled-cluster method, min {}, max {}",
                                   min_excit, max_excit))
                  .if_missing([]() { std::cerr << "--excit is required\n"; }) &
              clipp::integer("NUM", excit_));
  g.push_back(required("--eq")
                  .doc("Select a particular equation (1 to --ex)")
                  .if_missing([]() { std::cerr << "--eqn is required\n"; }) &
              clipp::integer("NUM", r_));

  g.push_back(option("--opt")
                  .doc("Optimize expressions for better flops")
                  .set(optimize_));
  // openshell options
  clipp::group os;
  os.push_back(
      (option("--osh")
           .doc("Expressions with openshell spin-integrated orbitals")
           .set(orbital_type_, OrbitalType::Openshell) &
       (required("--na")
            .doc("No. of particles with implicit alpha spin (0 to --eqn)")
            .if_missing([]() {
              std::cerr << "No. of particles with implicit alpha spin "
                           "required while specifying openshell equation\n";
            }) &
        clipp::integer("NUM", nalpha_)))
          .doc("Openshell options"));

  // closedshell options
  clipp::group cs;
  cs.push_back(
      (option("--csh")
           .doc("Expressions with closedshell spinfree orbitals")
           .set(orbital_type_, OrbitalType::Closedshell)
           .set(st_algo_, SpinTraceAlgorithm::Fast) &
       (option("--fast")
            .doc("Use faster spintracing algorithm (generates more "
                 "terms, default)")
            .set(st_algo_, SpinTraceAlgorithm::Fast) |
        option("--robust")
            .set(st_algo_, SpinTraceAlgorithm::Robust)
            .doc("Use robuster spintracing algorithm (takes longer time)")))
          .doc("Closedshell options"));
  g.push_back(os | cs);
  return g;
}

ExprPtr CCk::expr() const { return read(database_dir()); }

void CCk::write(std::filesystem::path const& odir) const {
  assert(std::filesystem::is_directory(odir));
  auto const E = excit();
  auto fnames = orbital_type() == OrbitalType::Spinorbital
                    ? canonical_fnames_cck(E, false, ".expr")
                    : canonical_fnames_cck(E, st_algorithm(), false, ".expr");
  std::clog << "Following files will be written:\n";
  for (auto const& f : fnames) std::clog << f << '\n';
  std::clog << "Generating equations" << std::endl;
  auto eqs = orbital_type() == OrbitalType::Spinorbital
                 ? cck_equations(E)
                 : cck_equations(E, st_algorithm());
  std::clog << "Done generating equations" << std::endl;
  assert(fnames.size() == eqs.size());
  for (auto&& [f, e] : ranges::views::zip(fnames, eqs)) {
    write_sum(e.as<Sum>(), std::filesystem::path{odir}.append(f));
  }
}

ExprPtr CCk::read(std::filesystem::path const& idir) const {
  assert(std::filesystem::is_directory(idir));
  auto fname =
      orbital_type() == OrbitalType::Openshell
          ? canonical_fname_cck(excit(), r(), nalpha(), false, ".expr")
          : canonical_fname_cck(excit(), r(), st_algorithm(), false, ".expr");

  auto fpath = std::filesystem::path{idir}.append(fname);

  if (!std::filesystem::is_regular_file(fpath)) write(idir);
  assert(std::filesystem::is_regular_file(fpath));
  return optimize()
             ? sequant::optimize(
                   read_sum(fpath),
                   [](Index const& idx) -> size_t {
                     return idx.space() == IndexSpace::active_occupied ? 2 : 20;
                   })
             : read_sum(fpath);
}

// /////////////////////////////////////////////////////////////////////////////
// Read implementation
// /////////////////////////////////////////////////////////////////////////////
Read::Read() : ExprSource{"read", "Read expressions from file"} {}

clipp::group Read::group() noexcept {
  return clipp::group(clipp::value("FILE")
                          .doc("File to read")
                          .call([this](std::string const& f) {
                            ifile_ = std::filesystem::path(f);
                          }));
}

std::filesystem::path const& Read::ifile() const noexcept { return ifile_; }

ExprPtr Read::expr() const { return read_sum(ifile()); }

}  // namespace sequant::cmdl