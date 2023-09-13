#include "utils.hpp"

#include <fmt/format.h>
#include <range/v3/view.hpp>

#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <fstream>

namespace sequant {

std::string canonical_fname_cck(size_t excit,                //
                                size_t r,                    //
                                SpinTraceAlgorithm st_algo,  //
                                bool optimize,               //
                                std::string_view ext) {
  using namespace fmt::literals;
  return fmt::format(
      std::locale::classic(),
      "CCk{ex}_R{r}{st}{opt}{ext}",  //
      "ex"_a = excit,                //
      "r"_a = r,                     //
      "st"_a = (st_algo == SpinTraceAlgorithm::Fast     ? "_fast-st"
                : st_algo == SpinTraceAlgorithm::Robust ? "_robust-st"
                                                        : ""),  //
      "opt"_a = (optimize ? "_opt" : ""),                       //
      "ext"_a = ext);
}

std::string canonical_fname_cck(size_t excit,   //
                                size_t r,       //
                                size_t nalpha,  //
                                bool optimize,  //
                                std::string_view ext) {
  using namespace fmt::literals;

  assert(0 <= nalpha && nalpha <= r);

  return fmt::format(std::locale::classic(),              //
                     "CCk{ex}_R{r}_{a}α{b}β{opt}{ext}",   //
                     "ex"_a = excit,                      //
                     "r"_a = r,                           //
                     "a"_a = nalpha,                      //
                     "b"_a = (r - nalpha),                //
                     "opt"_a = (optimize ? "_opt" : ""),  //
                     "ext"_a = ext);
}

template <typename Os>
void set_uft8(Os& os) {
  try {
    os.imbue(std::locale{"en_US.UTF-8"});
  } catch (std::exception const& ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << "Cannot read/write expressions without en_US.UTF8 encoding"
              << std::endl;
    std::exit(1);
  }
}

void write_sum(Sum const& sum, std::filesystem::path const& file) {
  auto ofs = std::wofstream{file};
  ofs.exceptions(std::ios::badbit | std::ios::failbit);
  set_uft8(ofs);
  try {
    for (auto const& term : sum) ofs << deparse_expr(term) << '\n';
  } catch (std::ios_base::failure const& ex) {
    std::cerr << "Caught exception: " << ex.what()
              << ", error code = " << ex.code() << std::endl;
    std::exit(1);
  }
}

ExprPtr read_sum(std::filesystem::path const& file) {
  auto ifs = std::wifstream{file};
  //
  // todo: cannot set this on MacOS AppleClang 14.0.3, std::getline fails
  // ifs.exceptions(std::ios::badbit | std::ios::failbit);
  //
  set_uft8(ifs);
  try {
    Sum result;
    for (std::wstring term; std::getline(ifs, term);)
      result.append(parse_expr(term));
    return ex<Sum>(std::move(result));
  } catch (std::ios_base::failure const& ex) {
    std::cerr << "Caught exception: " << ex.what()
              << ", error code = " << ex.code() << std::endl;
    std::exit(1);
  }
}

container::vector<ExprPtr> cck_equations(size_t excit) {
  using ranges::views::tail;

  auto cc_rs_ = sequant::mbpt::sr::cceqs(excit).t();
  assert(!cc_rs_[0]);

  // remove leading 'A' tensors
  for (auto& r : tail(cc_rs_)) r = opt::tail_factor(r);

  return tail(cc_rs_) | ranges::to<container::vector<ExprPtr>>;
}

container::vector<ExprPtr> cck_equations(size_t excit,
                                         SpinTraceAlgorithm st_algo) {
  using ranges::views::join;
  using ranges::views::tail;
  using ranges::views::transform;

  auto cc_rs_ = sequant::mbpt::sr::cceqs(excit).t();
  assert(!cc_rs_[0]);

  if (st_algo == SpinTraceAlgorithm::Fast)
    return tail(cc_rs_) | transform(closed_shell_CC_spintrace) |
           transform(opt::tail_factor) | ranges::to<container::vector<ExprPtr>>;

  if (st_algo == SpinTraceAlgorithm::Robust)
    return tail(cc_rs_) | transform(closed_shell_CC_spintrace_rigorous) |
           transform(opt::tail_factor) | ranges::to<container::vector<ExprPtr>>;

  assert(st_algo == SpinTraceAlgorithm::OpenShell);

  return tail(cc_rs_) |
         transform([](auto const& r) { return open_shell_CC_spintrace(r); }) |
         join | ranges::to<container::vector<ExprPtr>>;
}

std::filesystem::path database_dir() {
  return std::filesystem::path{std::getenv("HOME")}.append(database_dir_name);
}

void initialize_db() {
  auto db_dir = database_dir();
  if (!std::filesystem::exists(db_dir)) {
    try {
      std::filesystem::create_directories(db_dir);
    } catch (...) {
      std::cerr << fmt::format("Could not initialize database directory: {}",
                               db_dir.string())
                << std::endl;
      std::exit(1);
    }
  }
}

std::string canonical_name_cck(size_t excit, size_t r, std::string_view ext) {
  using namespace fmt::literals;
  return fmt::format(std::locale::classic(), "CCk{x}_R{r}{t}", "x"_a = excit,
                     "r"_a = r, "t"_a = ext);
}

container::vector<std::string> canonical_fnames_cck(size_t excit, bool optimize,
                                                    std::string_view ext) {
  using ranges::views::iota;
  using ranges::views::take;
  using ranges::views::transform;

  return iota(size_t{1}) | take(excit) |
         transform([excit, optimize, ext](size_t r) {
           return canonical_fname_cck(excit, r, SpinTraceAlgorithm::Invalid,
                                      optimize, ext);
         }) |
         ranges::to<container::vector<std::string>>;
}

container::vector<std::string> canonical_fnames_cck(size_t excit,
                                                    SpinTraceAlgorithm st_algo,
                                                    bool optimize,
                                                    std::string_view ext) {
  using ranges::views::iota;
  using ranges::views::join;
  using ranges::views::take;
  using ranges::views::transform;

  if (st_algo == SpinTraceAlgorithm::OpenShell) {
    return iota(size_t{1}) | take(excit) | transform([=](size_t r) {
             return iota(size_t{0}, r + 1) | transform([=](size_t nbeta) {
                      return canonical_fname_cck(excit, r, (r - nbeta),
                                                 optimize, ext);
                    });
           }) |
           join | ranges::to<container::vector<std::string>>;
  } else {
    return iota(size_t{1}) | take(excit) | transform([=](size_t r) {
             return canonical_fname_cck(excit, r, st_algo, optimize, ext);
           }) |
           ranges::to<container::vector<std::string>>;
  }
}

}  // namespace sequant
