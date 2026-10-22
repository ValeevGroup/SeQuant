#ifndef SEQUANT_TESTS_CATCH2_SEQUANT_H
#define SEQUANT_TESTS_CATCH2_SEQUANT_H

#include <catch2/catch_tostring.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <dtl/dtl.hpp>

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>
#include <variant>

namespace sequant::tests {
/// Shared particle-symmetric symmetry pack for the MBPT test TUs, which must
/// spell out MBPT particle symmetry explicitly since programmatic Tensor
/// construction is Context-independent (see Tensor::Defaults). Defined `inline`
/// here rather than once per TU so that unity test builds see one definition,
/// and in a test-only namespace rather than in `sequant` itself so that it
/// cannot collide with a library symbol.
inline constexpr TensorSymmetries particle_symmetric{.column =
                                                         ColumnSymmetry::Symm};

/// Sets environment variable @p name to @p value for the lifetime of this
/// object and restores its previous state (value or absence) on destruction,
/// including when a failed assertion unwinds the test case.
class ScopedEnv {
 public:
  ScopedEnv(std::string name, const char *value) : name_(std::move(name)) {
    if (const char *prev = std::getenv(name_.c_str())) prev_ = prev;
    set(name_.c_str(), value);
  }

  ~ScopedEnv() {
    if (prev_)
      set(name_.c_str(), prev_->c_str());
    else
      unset(name_.c_str());
  }

  ScopedEnv(const ScopedEnv &) = delete;
  ScopedEnv &operator=(const ScopedEnv &) = delete;

 private:
  std::string name_;
  std::optional<std::string> prev_;

  // POSIX has setenv/unsetenv, MSVC has neither and provides _putenv_s
  // instead (whose "" value removes the variable, matching unsetenv)
  static void set(const char *name, const char *value) {
#ifdef _WIN32
    _putenv_s(name, value);
#else
    setenv(name, value, 1);
#endif
  }

  static void unset(const char *name) {
#ifdef _WIN32
    _putenv_s(name, "");
#else
    unsetenv(name);
#endif
  }
};

/// a stream buffer that may be written to from several threads at once:
/// each thread's output is collected line by line, and only complete lines are
/// appended to the shared contents, so that lines do not interleave
class LineAtomicStreamBuf : public std::streambuf {
 public:
  /// @return the complete lines written so far
  std::string str() const {
    std::scoped_lock lock(mtx_);
    return contents_;
  }

 protected:
  int_type overflow(int_type ch) override {
    if (!traits_type::eq_int_type(ch, traits_type::eof()))
      put(traits_type::to_char_type(ch));
    return traits_type::not_eof(ch);
  }

  std::streamsize xsputn(const char_type *s, std::streamsize n) override {
    for (std::streamsize i = 0; i != n; ++i) put(s[i]);
    return n;
  }

 private:
  mutable std::mutex mtx_;
  std::string contents_;
  std::map<std::thread::id, std::string> lines_;

  void put(char_type ch) {
    std::scoped_lock lock(mtx_);
    auto &line = lines_[std::this_thread::get_id()];
    line.push_back(ch);
    if (ch == '\n') {
      contents_ += line;
      line.clear();
    }
  }
};

/// @return the number of Product canonicalizations performed by @p f , as
/// reported by the canonicalization Logger
/// @param methods the canonicalization methods to count
/// @note @p f may canonicalize concurrently
template <typename F>
std::size_t count_product_canonicalizations(
    F &&f, std::initializer_list<CanonicalizationMethod> methods = {
               CanonicalizationMethod::Topological,
               CanonicalizationMethod::Lexicographic,
               CanonicalizationMethod::Complete}) {
  struct Capture {
    bool logged = Logger::instance().canonicalize;
    LineAtomicStreamBuf buf;
    std::streambuf *cout_buf = std::cout.rdbuf(&buf);
    Capture() { Logger::instance().canonicalize = true; }
    ~Capture() {
      std::cout.rdbuf(cout_buf);
      Logger::instance().canonicalize = logged;
    }
  } capture;
  std::forward<F>(f)();
  const auto log = capture.buf.str();
  std::size_t count = 0;
  for (const auto method : methods) {
    const auto key =
        "Product canonicalization(" + to_string(method) + ") input: ";
    for (auto pos = log.find(key); pos != std::string::npos;
         pos = log.find(key, pos + key.size()))
      ++count;
  }
  return count;
}
/// @return a guard that scopes, for every Statistics, a copy of the default
/// context in effect for it whose canonicalization options are @p opts
[[nodiscard]] inline auto scoped_canonicalize_options(
    CanonicalizeOptions opts) {
  return set_scoped_modified_default_context(
      [&opts](Context &ctx) { ctx.set(opts); });
}

}  // namespace sequant::tests

namespace Catch {

// Make sure Catch uses proper string representation for SeQuant types

template <>
struct StringMaker<sequant::Expr> {
  static std::string convert(const sequant::Expr &expr,
                             bool include_canonical = true) {
    std::string str;
    try {
      str = sequant::serialize(expr, {.annot_symm = true});
    } catch (const std::exception &) {
      // serialize doesn't support all kinds of expressions -> fall back to
      // LaTeX representation
      str = sequant::to_latex(expr);
    }

    if (include_canonical) {
      sequant::ExprPtr clone = expr.clone();
      canonicalize(clone);
      simplify(clone);
      std::string canon_str =
          StringMaker<sequant::Expr>::convert(*clone, false);

      if (canon_str != str) {
        str += " (canonicalized: " + canon_str + ")";
      }
    }

    return str;
  }
};
template <>
struct StringMaker<sequant::ExprPtr> {
  static std::string convert(const sequant::ExprPtr &expr) {
    return StringMaker<sequant::Expr>::convert(*expr);
  }
};
template <>
struct StringMaker<sequant::Tensor> {
  static std::string convert(const sequant::Tensor &tensor) {
    return StringMaker<sequant::Expr>::convert(tensor);
  }
};

template <>
struct StringMaker<sequant::Index> {
  static std::string convert(const sequant::Index &idx) {
    return std::string(idx.full_label());
  }
};

template <>
struct StringMaker<sequant::ResultExpr> {
  static std::string convert(const sequant::ResultExpr &res,
                             bool include_canonical = true) {
    std::string str = sequant::serialize(res, {.annot_symm = true});

    if (include_canonical) {
      sequant::ResultExpr clone = res.clone();
      canonicalize(clone);
      simplify(clone);
      std::string canon_str =
          StringMaker<sequant::ResultExpr>::convert(clone, false);

      if (canon_str != str) {
        str += " (canonicalized: " + canon_str + ")";
      }
    }

    return str;
    std::stringstream sstream;
  }
};

}  // namespace Catch

namespace {

/// Matches two strings for equality. However, if they are not equal,
/// this matcher will produce a line-based diff of the strings rather than just
/// dumping the value of the string to the console
class StringDiffMatcher : public Catch::Matchers::MatcherGenericBase {
 public:
  StringDiffMatcher(std::string_view str) : expected_(str) {}

  bool match(const std::string &input) const {
    if (expected_ == input) {
      return true;
    }

    // Compute and save diff
    std::vector<std::string> expected_lines = to_lines(expected_);
    decltype(expected_lines) actual_lines = to_lines(input);

    dtl::Diff<std::string, decltype(expected_lines)> diff(actual_lines,
                                                          expected_lines);
    diff.compose();
    diff.composeUnifiedHunks();

    std::stringstream sstream;
    diff.printUnifiedFormat(sstream);

    diff_ = sstream.str();

    return false;
  }

  std::string describe() const override {
    return "requires modification:\n" + diff_;
  }

 private:
  std::string expected_;
  // We need to work around the fact that the match function is required to be
  // const, so the diff (which we can only compute inside match()) can be
  // updated there.
  mutable std::string diff_;

  static std::vector<std::string> to_lines(const std::string &input) {
    std::vector<std::string> lines;

    std::istringstream in(input);

    for (std::string line; std::getline(in, line);) {
      lines.push_back(std::move(line));
    }

    if (input.ends_with('\n')) {
      // Add terminating empty line
      lines.emplace_back("");
    }

    return lines;
  }
};

using ExprVar = std::variant<sequant::ExprPtr, sequant::ResultExpr>;

/// Converts the given expression-like object into an actual ExprPtr.
/// It accepts either an actual expression object (as Expr & or ExprPtr) or
/// a string-like object which will then be deserialized to yield the actual
/// expression object.
template <typename T>
ExprVar to_expression(T &&expression) {
  using std::begin;
  using std::end;

  using BaseT = std::remove_cvref_t<T>;

  if constexpr (std::is_convertible_v<BaseT, std::string_view>) {
    const std::string_view string = expression;

    if (string.find('=') != std::string_view::npos) {
      return sequant::deserialize<sequant::ResultExpr>(
          string, {.def_perm_symm = sequant::Symmetry::Nonsymm});
    } else {
      return sequant::deserialize<sequant::ExprPtr>(
          string, {.def_perm_symm = sequant::Symmetry::Nonsymm});
    }
  } else if constexpr (std::same_as<BaseT, sequant::ResultExpr>) {
    return expression;
  } else if constexpr (std::same_as<BaseT, sequant::Expr>) {
    // Clone in order to not have to worry about later modification
    return expression.clone();
  } else {
    static_assert(std::same_as<BaseT, sequant::ExprPtr>,
                  "Invalid type for expression");

    // Clone in order to not have to worry about later modification
    return expression->clone();
  }
}

template <typename Subclass>
class ExpressionMatcher : public Catch::Matchers::MatcherGenericBase {
 public:
  template <typename T>
  ExpressionMatcher(T &&expression)
      : m_expr(to_expression(std::forward<T>(expression))) {
    if (std::holds_alternative<sequant::ResultExpr>(m_expr)) {
      Subclass::pre_comparison(
          std::get<sequant::ResultExpr>(m_expr).expression());
    } else {
      Subclass::pre_comparison(std::get<sequant::ExprPtr>(m_expr));
    }
  }

  bool match(const sequant::ResultExpr &res) const {
    if (!std::holds_alternative<sequant::ResultExpr>(m_expr)) {
      return false;
    }

    const sequant::ResultExpr &self = std::get<sequant::ResultExpr>(m_expr);

    return self.has_label() == res.has_label() && self.label() == res.label() &&
           self.symmetry() == res.symmetry() &&
           self.braket_symmetry() == res.braket_symmetry() &&
           self.column_symmetry() == res.column_symmetry() &&
           do_match(*res.expression());
  }

  bool match(const sequant::ExprPtr &expr) const {
    return std::holds_alternative<sequant::ExprPtr>(m_expr) && do_match(*expr);
  }

  bool match(const sequant::Expr &expr) const {
    return std::holds_alternative<sequant::ExprPtr>(m_expr) && do_match(expr);
  }

  std::string describe() const override {
    return Subclass::comparison_requirement() + ": " + stringify(m_expr);
  }

 protected:
  ExprVar m_expr;

  bool do_match(const sequant::Expr &expr) const {
    // Never modify the expression that we are trying to check in order to avoid
    // side-effects
    sequant::ExprPtr clone = expr.clone();

    Subclass::pre_comparison(clone);

    if (std::holds_alternative<sequant::ResultExpr>(m_expr)) {
      return *clone == *std::get<sequant::ResultExpr>(m_expr).expression();
    }

    return *clone == *std::get<sequant::ExprPtr>(m_expr);
  }

  std::string stringify(const ExprVar &expr) const {
    if (std::holds_alternative<sequant::ResultExpr>(expr)) {
      return Catch::Detail::stringify(std::get<sequant::ResultExpr>(expr));
    }

    return Catch::Detail::stringify(std::get<sequant::ExprPtr>(expr));
  }
};

/// Matches that the tested expression is equivalent to the given one. Two
/// expressions are considered equivalent if they both have the same canonical
/// form (i.e. they are the same expression after canonicalization)
struct EquivalentToMatcher : ExpressionMatcher<EquivalentToMatcher> {
  using ExpressionMatcher::ExpressionMatcher;

  static void pre_comparison(sequant::ExprPtr &expr) {
    auto scope = sequant::tests::scoped_canonicalize_options(
        sequant::CanonicalizeOptions::default_options().copy_and_set(
            sequant::CanonicalizeOptions::IgnoreNamedIndexLabel::No));
    sequant::simplify(expr);
  }

  static std::string comparison_requirement() { return "Equivalent to"; }
};

/// Matches that the tested expression simplifies (**without**
/// re-canonicalization!) to the same form as the given one
struct SimplifiesToMatcher : ExpressionMatcher<SimplifiesToMatcher> {
  using ExpressionMatcher::ExpressionMatcher;

  static void pre_comparison(sequant::ExprPtr &expr) {
    sequant::rapid_simplify(expr);
  }

  static std::string comparison_requirement() { return "Simplifies to"; }
};

}  // namespace

template <typename T>
EquivalentToMatcher EquivalentTo(T &&expression) {
  return EquivalentToMatcher(std::forward<T>(expression));
}

template <typename T>
SimplifiesToMatcher SimplifiesTo(T &&expression) {
  return SimplifiesToMatcher(std::forward<T>(expression));
}

template <typename String>
StringDiffMatcher DiffedStringEquals(String &&str) {
  return StringDiffMatcher(std::forward<String>(str));
}

#endif
