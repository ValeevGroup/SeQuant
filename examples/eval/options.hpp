//
// Created by Bimal Gaudel on 7/28/21.
//

#ifndef SEQUANT_EVAL_OPTIONS_HPP
#define SEQUANT_EVAL_OPTIONS_HPP

#include <SeQuant/core/container.hpp>
#include <fstream>
#include <iomanip>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>
#include <sstream>

namespace sequant::eval {

struct OptionsEquations {
  inline static constexpr size_t MIN_EXCIT = 2;
  inline static constexpr size_t MAX_EXCIT = 6;
  size_t excit = MIN_EXCIT;
  bool spintrace = false;
};

struct OptionsOptimization {
  bool single_term = true;
  bool reuse_imeds = true;
  bool cache_leaves = true;
};

struct OptionsSCF {
  size_t max_iter;
  double conv;
};

inline static const OptionsSCF LooseSCF{100, 1e-10};
inline static const OptionsSCF NormalSCF{200, 1e-12};
inline static const OptionsSCF TightSCF{400, 1e-14};

struct OptionsLog {
  inline static constexpr size_t MIN_LEVEL = 0;
  inline static constexpr size_t MAX_LEVEL = 1;
  size_t level = 1;
  std::string file = "";
};

namespace detail {

struct ErrorArgValInvalid : public std::runtime_error {
  ErrorArgValInvalid() : std::runtime_error{"invalid value for the argument"} {}
};

struct ErrorArgValOutOfRange : public std::runtime_error {
  ErrorArgValOutOfRange()
      : std::runtime_error{"value out of range for the argument"} {}
};

struct ErrorArgNameInvalid : public std::runtime_error {
  ErrorArgNameInvalid(std::string const& name)
      : std::runtime_error{"invalid argument name " + name} {}
};

inline std::string help(std::string_view arg_name,  //
                        int lbound, int ubound) {
  using namespace std::string_literals;
  std::ostringstream oss{};
  oss << std::setw(15) << std::left << arg_name << "(";
  oss << lbound << ".." << ubound << ")";
  return oss.str();
}

template <typename Iterable>
inline std::string help(std::string_view arg_name, Iterable const& container) {
  using namespace std::string_literals;
  std::ostringstream oss{};
  oss << std::setw(15) << std::left << arg_name << "(";
  oss << ranges::front(container);
  for (auto&& v : ranges::views::tail(container)) oss << "|" << v;
  oss << ")";
  return oss.str();
}

struct ArgValBool {
  const container::svector<std::string_view> true_vals;

  const container::svector<std::string_view> false_vals;

  ArgValBool(std::initializer_list<std::string_view> yes = {"yes"},
             std::initializer_list<std::string_view> no = {"no"})
      : true_vals{yes}, false_vals{no} {}

  [[nodiscard]] bool parse(std::string_view v) const;

  [[nodiscard]] std::string help(std::string_view arg_name) const;
};

struct ArgValInt {
  int const lbound;
  int const ubound;

  ArgValInt(int lo, int hi) : lbound{lo}, ubound{hi} {}

  [[nodiscard]] int parse(std::string_view arg_val) const;

  [[nodiscard]] std::string help(std::string_view arg_name) const;
};

struct ArgValList {
  std::vector<std::string> const list;

  ArgValList(std::initializer_list<std::string> const& allowed_vals)
      : list{allowed_vals} {}

  [[nodiscard]] std::string parse(std::string_view arg_value) const;

  [[nodiscard]] std::string help(std::string_view arg_name) const;
};

class ParseOptionsEquations {
 private:
  const ArgValInt excit_parser{OptionsEquations::MIN_EXCIT,
                               OptionsEquations::MAX_EXCIT};

  inline static const auto spintrace_parser = ArgValBool{};

  OptionsEquations opts_{};

  inline static std::string_view const spintrace{"spintrace"};

  inline static std::string_view const excitation{"excitation"};

 public:
  [[nodiscard]] std::string help() const;

  [[nodiscard]] OptionsEquations opts() const;

  void update(std::string_view arg_name, std::string_view value);
};

class ParseOptionsOptimization {
 private:
  inline static const auto bool_parser = ArgValBool{};

  inline static std::string_view const single_term{"single_term"};

  inline static std::string_view const reuse_imeds{"reuse_imeds"};

  inline static std::string_view const cache_leaves{"cache_leaves"};

  OptionsOptimization opts_;

 public:
  [[nodiscard]] std::string help() const;

  [[nodiscard]] OptionsOptimization opts() const;

  void update(std::string_view arg_name, std::string_view value);
};

class ParseOptionsSCF {
 private:
  inline static std::string const Loose{"loose"};

  inline static std::string const Normal{"normal"};

  inline static std::string const Tight{"tight"};

  inline static ArgValList const keyword_parser =
      ArgValList{{Loose, Normal, Tight}};

  OptionsSCF opts_{NormalSCF};

  std::string_view const tightness = "tightness";

  std::string_view const maxiter   = "maxiter";

 public:
  [[nodiscard]] std::string help() const;

  [[nodiscard]] OptionsSCF opts() const;

  void update(std::string_view arg_name, std::string_view value);
};

class ParseOptionsLog {
 private:
  inline static ArgValInt const level_parser{OptionsLog::MIN_LEVEL,
                                             OptionsLog::MAX_LEVEL};

  inline static std::string_view const level{"level"};

  OptionsLog opts_{};

 public:
  [[nodiscard]] std::string help() const;

  [[nodiscard]] OptionsLog opts() const;

  void update(std::string_view arg_name, std::string_view value);
};

}  // namespace detail

class ParseConfigFile {
 private:
  inline static std::string_view const eqs{"equation"};
  inline static std::string_view const optm{"optimization"};
  inline static std::string_view const scf{"scf"};
  inline static std::string_view const log{"log"};

  detail::ParseOptionsEquations parse_eqs{};

  detail::ParseOptionsOptimization parse_optm{};

  detail::ParseOptionsSCF parse_scf{};

  detail::ParseOptionsLog parse_log{};

 public:
  std::string help() const;

  void parse(std::string_view fname);

  OptionsEquations opts_equations() const;

  OptionsOptimization opts_optimization() const;

  OptionsSCF opts_scf() const;

  OptionsLog opts_log() const;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_OPTIONS_HPP
