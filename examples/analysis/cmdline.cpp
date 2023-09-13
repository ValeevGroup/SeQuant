#include "cmdline.hpp"

#include <fmt/format.h>
#include <SeQuant/core/optimize.hpp>

#include <fstream>

namespace sequant {

size_t CCkOpts::excit() const noexcept { return excit_; }

size_t CCkOpts::r() const noexcept { return r_; }

size_t CCkOpts::nalpha() const noexcept { return nalpha_; }

bool CCkOpts::optimize() const noexcept { return optimize_; }

OrbitalType CCkOpts::orbital_type() const noexcept { return orbital_type_; }

SpinTraceAlgorithm CCkOpts::st_algorithm() const noexcept { return st_algo_; }

void CCkOpts::sanity_check() const {
  if (excit() < min_excit || excit() > max_excit) {
    std::cerr << fmt::format("CC excitation {} not in the range [{},{}]\n",
                             excit(), min_excit, max_excit);
    std::exit(1);
  }
  if (r() < 1 || r() > excit()) {
    std::cerr << fmt::format("CC equation index {} not in the range [1,{}]\n",
                             r(), excit());
    std::exit(1);
  }
  if (orbital_type() == OrbitalType::Openshell && nalpha() > excit()) {
    std::cerr << fmt::format("No. of alpha value {} not in the range [0,{}]\n",
                             nalpha(), r());
    std::exit(1);
  }
}

void CCkOpts::add_options(clipp::group& g) {
  using clipp::option;
  using clipp::required;
  using clipp::value;

  g.push_back(
      required("--excit")
          .if_missing([]() { std::cerr << "--excit is required\n"; })
          .doc(fmt::format("The highest level of excitation in coupled-cluster "
                           "method min {}, max {}",
                           CCkOpts::min_excit, CCkOpts::max_excit, excit_)) &
      value("NUM", excit_));
  g.push_back(required("--eqn")
                  .if_missing([]() { std::cerr << "--eqn is required\n"; })
                  .doc("Select a particular equation (1 to --excit)") &
              value("NUM", r_));
  g.push_back(option("--optimize")
                  .doc("Optimize expression for better flops")
                  .set(optimize_));

  // openshell options
  clipp::group os;
  os.push_back(option("--oshell")
                   .doc("Expressions with openshell spinfree orbitals")
                   .set(orbital_type_, OrbitalType::Openshell) &
               (required("--nalpha")
                    .doc("No. of particles with alpha label (0 to --eqn)")
                    .if_missing([]() {
                      std::cerr
                          << "No. of particles with alpha labels required "
                             "while specifying openshell equation\n";
                    }) &
                value("NUM", nalpha_).set(nalpha_)));

  // closedshell options
  clipp::group cs;
  cs.push_back(
      option("--cshell")
          .doc("Expressions with closedshell spinfree orbitals")
          .set(orbital_type_, OrbitalType::Closedshell)
          .set(st_algo_, SpinTraceAlgorithm::Fast) &
      (option("--fast")
           .doc("Use faster spintracing algorithm (generates more "
                "terms, default)")
           .set(st_algo_, SpinTraceAlgorithm::Fast) |
       option("--robust")
           .set(st_algo_, SpinTraceAlgorithm::Robust)
           .doc("Use robuster spintracing algorithm (takes longer time)")));
  g.push_back(os | cs);
}

size_t LatexOpts::terms_per_line() const noexcept { return terms_; }

size_t LatexOpts::lines_per_block() const noexcept { return lines_; }

void LatexOpts::add_options(clipp::group& g) {
  g.push_back((clipp::option("-t", "--terms-per-line")
                       .doc(fmt::format("No. of terms per line, default {}",
                                        terms_per_line())) &
                   clipp::value("NUM", terms_),
               clipp::option("-l", "--lines-per-block")
                       .doc("No. of lines per block, default single block with "
                            "all lines") &
                   clipp::value("NUM", lines_)) %
              "Control latex output");
}

void GraphOpts::add_options(clipp::group& g) {
  // todo
}

void write_cck(CCkOpts const& opts, std::filesystem::path const& odir) {
  assert(std::filesystem::is_directory(odir));
  auto const E = opts.excit();
  auto fnames =
      opts.orbital_type() == OrbitalType::Spinorbital
          ? canonical_fnames_cck(E, false, ".expr")
          : canonical_fnames_cck(E, opts.st_algorithm(), false, ".expr");
  std::clog << "Following files will be written:\n";
  for (auto const& f : fnames) std::clog << f << '\n';
  std::clog << "Generating equations" << std::endl;
  auto eqs = opts.orbital_type() == OrbitalType::Spinorbital
                 ? cck_equations(E)
                 : cck_equations(E, opts.st_algorithm());
  std::clog << "Done generating equations" << std::endl;
  assert(fnames.size() == eqs.size());
  for (auto&& [f, e] : ranges::views::zip(fnames, eqs)) {
    write_sum(e.as<Sum>(), std::filesystem::path{odir}.append(f));
  }
}

ExprPtr read_cck(CCkOpts const& opts, std::filesystem::path const& idir) {
  assert(std::filesystem::is_directory(idir));
  auto fname = opts.orbital_type() == OrbitalType::Openshell
                   ? canonical_fname_cck(opts.excit(), opts.r(), opts.nalpha(),
                                         false, ".expr")
                   : canonical_fname_cck(opts.excit(), opts.r(),
                                         opts.st_algorithm(), false, ".expr");

  auto fpath = std::filesystem::path{idir}.append(fname);

  if (!std::filesystem::is_regular_file(fpath)) write_cck(opts, idir);
  assert(std::filesystem::is_regular_file(fpath));
  return opts.optimize()
             ? optimize(
                   read_sum(fpath),
                   [](Index const& idx) -> size_t {
                     return idx.space() == IndexSpace::active_occupied ? 2 : 20;
                   })
             : read_sum(fpath);
}

}  // namespace sequant
