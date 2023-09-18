#ifndef SEQUANT_ANALYSIS_CMDLINE_HPP
#define SEQUANT_ANALYSIS_CMDLINE_HPP

#include "utils.hpp"

#include <clipp.h>

namespace sequant::cmdl {

struct Option {
  [[nodiscard]] virtual clipp::group group() noexcept = 0;

 protected:
  Option() = default;
};

class Command : public Option {
 public:
  using args_type = container::vector<std::string>;

  [[nodiscard]] std::string_view name() const noexcept;

  [[nodiscard]] std::string_view doc() const noexcept;

  [[maybe_unused]] clipp::parsing_result parse(args_type const&);

  void parse(int argc, char** argv);

  [[nodiscard]] bool set() const noexcept;

  void set(bool) noexcept;

 protected:
  explicit Command(std::string_view name, std::string_view doc);

 private:
  bool is_set_;
  std::string_view name_;
  std::string_view doc_;
};

template <typename Iterable>
clipp::group one_of_cmds(Iterable& cmds) {
  static_assert(std::is_same_v<decltype(**ranges::begin(cmds)), Command&>);
  clipp::group result;
  for (auto& cm : cmds) {
    auto cm_grp = clipp::group(
        clipp::command(cm->name().data()).doc(cm->doc().data()).call([&cm]() {
          cm->set(true);
        }));
    result = result.empty() ? cm_grp : result | cm_grp;
  }
  return result;
}

struct ExprSource : public Command {
  [[nodiscard]] virtual ExprPtr expr() const = 0;

  using Command::Command;
};

struct ExprAction : public Command {
  virtual void run(std::ostream&, ExprPtr const&) const = 0;

  using Command::Command;
};

// /////////////////////////////////////////////////////////////////////////////
// classes derived from ExprAction follow
// /////////////////////////////////////////////////////////////////////////////
struct Latex : public ExprAction {
 public:
  Latex();

  void run(std::ostream&, ExprPtr const&) const override;

  [[nodiscard]] clipp::group group() noexcept override;

 private:
  size_t tpl_{5};
  size_t lpb_{0};
};

struct Graph : public ExprAction {
 public:
  enum struct Format { AdjSet, AdjMat, Tikz, Dot, Invalid };
  inline static container::set<
      std::tuple<std::string_view, std::string_view, Format>>
      format_labels{{"set", "Adjacency set", Format::AdjSet},
                    {"mat", "Adjacency matrix", Format::AdjMat},
                    {"tikz", "Tikz (LaTeX)", Format::Tikz},
                    {"dot", "Dot language", Format::Dot}};

  Graph();

  void run(std::ostream&, ExprPtr const&) const override;

  [[nodiscard]] clipp::group group() noexcept override;

  [[nodiscard]] Format format() const noexcept;

 private:
  Format format_{Format::AdjSet};
};

struct Expr : public ExprAction {
 public:
  Expr();

  void run(std::ostream&, ExprPtr const&) const override;

  [[nodiscard]] clipp::group group() noexcept override;
};

// /////////////////////////////////////////////////////////////////////////////
// classes derived from ExprSource follow
// /////////////////////////////////////////////////////////////////////////////

enum struct OrbitalType { Spinorbital, Closedshell, Openshell };

class CCk : public ExprSource {
 public:
  CCk();

  [[nodiscard]] bool optimize() const noexcept;

  [[nodiscard]] size_t excit() const noexcept;

  [[nodiscard]] size_t r() const noexcept;

  [[nodiscard]] size_t nalpha() const noexcept;

  [[nodiscard]] OrbitalType orbital_type() const noexcept;

  [[nodiscard]] SpinTraceAlgorithm st_algorithm() const noexcept;

  [[nodiscard]] clipp::group group() noexcept override;

  [[nodiscard]] ExprPtr expr() const override;

  void sanity_check() const;

  inline static constexpr size_t min_excit = 2;
  inline static constexpr size_t max_excit = 4;

 private:
  bool optimize_{false};

  size_t excit_{2};

  size_t r_{1};

  size_t nalpha_{0};

  OrbitalType orbital_type_{OrbitalType::Spinorbital};

  SpinTraceAlgorithm st_algo_{SpinTraceAlgorithm::Invalid};

  [[nodiscard]] ExprPtr read(std::filesystem::path const& idir) const;

  void write(std::filesystem::path const& odir) const;
};

class Read : public ExprSource {
 public:
  Read();

  [[nodiscard]] clipp::group group() noexcept override;

  [[nodiscard]] ExprPtr expr() const override;

  [[nodiscard]] std::filesystem::path const& ifile() const noexcept;

 private:
  std::filesystem::path ifile_{};
};

}  // namespace sequant::cmdl

#endif  // SEQUANT_ANALYSIS_CMDLINE_HPP
