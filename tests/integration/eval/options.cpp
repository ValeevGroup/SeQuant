//
// Created by Bimal Gaudel on 7/28/21.
//

#include "options.hpp"

namespace sequant::eval {

namespace detail {

bool ArgValBool::parse(std::string_view v) const {
  if (ranges::contains(true_vals, v))
    return true;
  else if (ranges::contains(false_vals, v))
    return false;
  else
    throw detail::ErrorArgValInvalid{};
}
std::string ArgValBool::help(std::string_view arg_name) const {
  return detail::help(arg_name, ranges::views::concat(true_vals, false_vals));
}
int ArgValInt::parse(std::string_view arg_val) const {
  long val;
  try {
    val = std::stol(arg_val.data());
    if (val < lbound || val > ubound) throw detail::ErrorArgValOutOfRange{};
  } catch (std::invalid_argument const& ex) {
    throw detail::ErrorArgValInvalid{};
  } catch (std::out_of_range const& ex) {
    throw detail::ErrorArgValOutOfRange{};
  }
  return static_cast<int>(val);
}
std::string ArgValInt::help(std::string_view arg_name) const {
  return detail::help(arg_name, lbound, ubound);
}
std::string ArgValList::parse(std::string_view arg_value) const {
  if (!ranges::contains(list, arg_value.data()))
    throw detail::ErrorArgValInvalid{};
  return arg_value.data();
}
std::string ArgValList::help(std::string_view arg_name) const {
  return detail::help(arg_name, list);
}
std::string ParseOptionsEquations::help() const {
  return excit_parser.help(excitation) + "\n" +
         spintrace_parser.help(spintrace);
}
OptionsEquations ParseOptionsEquations::opts() const { return opts_; }

void ParseOptionsEquations::update(std::string_view arg_name,
                                   std::string_view value) {
  if (arg_name == excitation)
    opts_.excit = excit_parser.parse(value);
  else if (arg_name == spintrace)
    opts_.spintrace = spintrace_parser.parse(value);
  else
    throw detail::ErrorArgNameInvalid{arg_name.data()};
}

std::string ParseOptionsOptimization::help() const {
  return bool_parser.help(single_term) + "\n" + bool_parser.help(reuse_imeds) +
         "\n" + bool_parser.help(cache_leaves);
}

OptionsOptimization ParseOptionsOptimization::opts() const { return opts_; }

void ParseOptionsOptimization::update(std::string_view arg_name,
                                      std::string_view value) {
  if (arg_name == single_term)
    opts_.single_term = bool_parser.parse(value);
  else if (arg_name == reuse_imeds)
    opts_.reuse_imeds = bool_parser.parse(value);
  else if (arg_name == cache_leaves)
    opts_.cache_leaves = bool_parser.parse(value);
  else
    throw detail::ErrorArgNameInvalid{arg_name.data()};
}

std::string ParseOptionsSCF::help() const {
  return keyword_parser.help(tightness);
}

OptionsSCF ParseOptionsSCF::opts() const { return opts_; }

void ParseOptionsSCF::update(std::string_view arg_name,
                             std::string_view value) {
  if (arg_name == maxiter) {
    opts_.max_iter = std::stoi(value.data());
    return;
  }

  if (arg_name != tightness) throw detail::ErrorArgNameInvalid{arg_name.data()};

  auto const val = keyword_parser.parse(value);
  if (val == Loose)
    opts_ = LooseSCF;
  else if (val == Normal)
    opts_ = NormalSCF;
  else if (val == Tight)
    opts_ = TightSCF;
  else {
  }  // unreachable
}

std::string ParseOptionsLog::help() const { return level_parser.help(level); }

OptionsLog ParseOptionsLog::opts() const { return opts_; }

void ParseOptionsLog::update(std::string_view arg_name,
                             std::string_view value) {
  if (arg_name != level) throw detail::ErrorArgNameInvalid{arg_name.data()};
  opts_.level = level_parser.parse(value);
}

}  // namespace detail

std::string ParseConfigFile::help() const {
  std::ostringstream oss{};
  oss << eqs << "\n"
      << parse_eqs.help() << "\n\n"
      << optm << "\n"
      << parse_optm.help() << "\n\n"
      << scf << "\n"
      << parse_scf.help() << "\n\n"
      << log << "\n"
      << parse_log.help();
  return oss.str();
}

void ParseConfigFile::parse(std::string_view fname) {
  auto ifs = std::ifstream{fname.data()};
  if (!ifs.good()) throw std::runtime_error{"unable to read file"};

  auto const headers = std::vector<std::string_view>{eqs, optm, scf, log};
  std::string token{}, curr_header{};

  while (ifs >> token) {
    if (ranges::contains(headers, token))
      curr_header = token;
    else if (curr_header.empty())
      throw std::runtime_error{"Invalid line in config file"};
    else {
      std::string val{};
      ifs >> val;

      if (curr_header == eqs)
        parse_eqs.update(token, val);
      else if (curr_header == optm)
        parse_optm.update(token, val);
      else if (curr_header == scf)
        parse_scf.update(token, val);
      else if (curr_header == log)
        parse_log.update(token, val);
      else {
        throw std::runtime_error{"unknown token " + token};
      }
    }
  }
}

OptionsEquations ParseConfigFile::opts_equations() const {
  return parse_eqs.opts();
}

OptionsOptimization ParseConfigFile::opts_optimization() const {
  return parse_optm.opts();
}

OptionsSCF ParseConfigFile::opts_scf() const { return parse_scf.opts(); }

OptionsLog ParseConfigFile::opts_log() const { return parse_log.opts(); }
}  // namespace sequant::eval
