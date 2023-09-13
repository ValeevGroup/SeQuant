#ifndef SEQUANT_ANALYSIS_CMDLINE_HPP
#define SEQUANT_ANALYSIS_CMDLINE_HPP

#include "utils.hpp"

#include <clipp.h>

namespace sequant {

enum struct OrbitalType { Spinorbital, Closedshell, Openshell };

struct Options {
 protected:
  Options() = default;

 public:
  virtual void add_options(clipp::group&) = 0;
};

class CCkOpts : public Options {
 public:
  inline static constexpr size_t min_excit = 2;
  inline static constexpr size_t max_excit = 4;

  CCkOpts() = default;

  void add_options(clipp::group& g) override;

  [[nodiscard]] bool optimize() const noexcept;

  [[nodiscard]] size_t excit() const noexcept;

  [[nodiscard]] size_t r() const noexcept;

  [[nodiscard]] size_t nalpha() const noexcept;

  [[nodiscard]] OrbitalType orbital_type() const noexcept;

  [[nodiscard]] SpinTraceAlgorithm st_algorithm() const noexcept;

  void sanity_check() const;

 private:
  bool optimize_{false};

  size_t excit_{2};

  size_t r_{1};

  size_t nalpha_{0};

  OrbitalType orbital_type_{OrbitalType::Spinorbital};

  SpinTraceAlgorithm st_algo_{SpinTraceAlgorithm::Invalid};
};

struct LatexOpts : public Options {
 public:
  LatexOpts() = default;

  void add_options(clipp::group&) override;

  [[nodiscard]] size_t terms_per_line() const noexcept;

  [[nodiscard]] size_t lines_per_block() const noexcept;

 private:
  size_t terms_{5};
  size_t lines_{0};
};

struct GraphOpts : public Options {
 public:
  GraphOpts() = default;
  void add_options(clipp::group&) override;
};

void write_cck(CCkOpts const& opts, std::filesystem::path const& odir);

ExprPtr read_cck(CCkOpts const& opts, std::filesystem::path const& idir);

}  // namespace sequant

#endif  // SEQUANT_ANALYSIS_CMDLINE_HPP
