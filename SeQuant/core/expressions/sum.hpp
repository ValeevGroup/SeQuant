#ifndef SEQUANT_EXPRESSIONS_SUM_HPP
#define SEQUANT_EXPRESSIONS_SUM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <optional>
#include <type_traits>

namespace sequant {

/// @brief sum of zero or more summands

/// Sum is associative and is flattened automatically.
class Sum : public Expr {
 public:
  using summands_type = container::svector<ExprPtr, 2>;

  Sum() = default;
  virtual ~Sum() = default;
  Sum(const Sum &) = default;
  Sum(Sum &&) = default;
  Sum &operator=(const Sum &) = default;
  Sum &operator=(Sum &&) = default;

  void swap(Sum &other) {
    Sum tmp = std::move(other);
    other = std::move(*this);
    *this = std::move(tmp);
  }

  /// construct a Sum out of zero or more summands
  /// @param summands an initializer list of summands
  Sum(ExprPtrList summands) {
    // use append to flatten out Sum summands
    for (auto &&summand : summands) {
      append(std::forward<decltype(summand)>(summand));
    }
  }

  /// construct a Sum out of a range of summands
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename Iterator>
  Sum(Iterator begin, Iterator end) {
    // use append to flatten out Sum summands
    for (auto it = begin; it != end; ++it) {
      append(*it);
    }
  }

  /// construct a Sum out of a range of summands
  /// @param rng a range
  template <typename Range>
    requires(meta::is_range_v<std::remove_cvref_t<Range>> &&
             !meta::is_same_v<std::remove_cvref_t<Range>, ExprPtrList>)
  explicit Sum(Range &&rng) {
    // N.B. use append to flatten out Sum summands
    constexpr auto rng_is_expr =
        meta::is_base_of_v<Expr, std::remove_cvref_t<Range>>;
    constexpr auto rng_is_exprptr =
        meta::is_same_v<ExprPtr, std::remove_cvref_t<Range>>;
    if constexpr (rng_is_expr || rng_is_exprptr) {
      ExprPtr rng_as_exprptr;
      if constexpr (rng_is_expr) {
        rng_as_exprptr = rng.exprptr_from_this();
      } else {
        rng_as_exprptr = rng;
      }
      this->append(rng_as_exprptr);
    } else {
      for (auto &&v : rng) {
        append(std::forward<decltype(v)>(v));
      }
    }
  }

  /// tags ctor to move the summands directly
  struct move_only_tag {};

  /// construct a Sum by moving in the summands, no flattening is performed,
  /// but zeros will be omitted and constants added up
  /// @param summands the summands to move in
  explicit Sum(summands_type &&summands, move_only_tag)
      : summands_(std::move(summands)) {
    std::size_t pos = 0;
    for (auto it = summands_.begin(); it != summands_.end(); ++it) {
      auto &summand = *it;
      bool do_erase = false;
      if (summand->is_zero()) {
        do_erase = true;
      } else if (summand->is<Constant>()) {
        auto summand_constant = summand.as_shared_ptr<Constant>();
        if (constant_summand_idx_) {  // add up to the existing constant ...
          SEQUANT_ASSERT(summands_.at(*constant_summand_idx_)->is<Constant>());
          *summands_[*constant_summand_idx_] += *summand_constant;
          do_erase = true;
        } else {  // or memorize the position of the constant
          constant_summand_idx_ = pos;
        }
      }

      // erase if needed
      if (do_erase) {
        summands.erase(it);
        it = summands_.begin();
        std::advance(it, pos);
      } else
        ++pos;
    }
  }

  /// append a summand to the sum
  /// @param summand the summand
  Sum &append(ExprPtr summand) {
    SEQUANT_ASSERT(summand);
    if (!summand->is<Sum>()) {
      if (!summand->is_zero()) {        // exclude zeros
        if (summand->is<Constant>()) {  // add up constants
          // immediately, if possible
          auto summand_constant = summand.as_shared_ptr<Constant>();
          if (constant_summand_idx_) {
            SEQUANT_ASSERT(
                summands_.at(*constant_summand_idx_)->is<Constant>());
            *(summands_[*constant_summand_idx_]) += *summand;
          } else {
            summands_.push_back(summand->clone());
            constant_summand_idx_ = summands_.size() - 1;
          }
        } else {
          summands_.push_back(summand->clone());
        }
        reset_hash_value();
      }
    } else {  // this recursively flattens Sum summands
      for (auto &subsummand : *summand) this->append(subsummand);
    }
    return *this;
  }

  /// prepend a summand to the sum
  /// @param summand the summand
  Sum &prepend(ExprPtr summand) {
    SEQUANT_ASSERT(summand);
    if (!summand->is<Sum>()) {
      if (!summand->is_zero()) {
        // exclude zeros
        if (summand->is<Constant>()) {
          auto summand_constant = summand.as_shared_ptr<Constant>();
          if (constant_summand_idx_) {  // add up to the existing constant ...
            SEQUANT_ASSERT(
                summands_.at(*constant_summand_idx_)->is<Constant>());
            *summands_[*constant_summand_idx_] += *summand_constant;
          } else {  // or include the nonzero constant and update
            // constant_summand_idx_
            summands_.insert(summands_.begin(), summand->clone());
            constant_summand_idx_ = 0;
          }
        } else {
          summands_.insert(summands_.begin(), summand->clone());
          if (constant_summand_idx_)  // if have a constant, update its position
            ++*constant_summand_idx_;
        }
        reset_hash_value();
      }
    } else {  // this recursively flattens Sum summands
      for (auto &subsummand : *summand) this->prepend(subsummand);
    }
    return *this;
  }

  /// Summands accessor
  const auto &summands() const { return summands_; }

  /// Summand accessor
  /// @param i summand index
  /// @return ith summand
  const ExprPtr &summand(size_t i) const { return summands_.at(i); }

  /// Takes the first @c count elements of the sum
  ExprPtr take_n(size_t count) const {
    const auto e = (count >= summands_.size() ? summands_.end()
                                              : (summands_.begin() + count));
    return ex<Sum>(summands_.begin(), e);
  }

  /// Takes the first @c count elements of the sum starting with element @c
  /// offset
  ExprPtr take_n(size_t offset, size_t count) const {
    const auto offset_plus_count = offset + count;
    const auto b = (offset >= summands_.size() ? summands_.end()
                                               : (summands_.begin() + offset));
    const auto e = (offset_plus_count >= summands_.size()
                        ? summands_.end()
                        : (summands_.begin() + offset_plus_count));
    return ex<Sum>(b, e);
  }

  /// @tparam Filter a boolean predicate type, such `Filter(const ExprPtr&)`
  /// evaluates to true
  /// @param f an object of Filter type
  /// Selects elements {`e`} for which `f(e)` is true
  template <typename Filter>
  ExprPtr filter(Filter &&f) const {
    return ex<Sum>(summands_ | ranges::views::filter(f));
  }

  /// @return true if the number of factors is zero
  bool empty() const { return summands_.empty(); }

  /// @return the number of summands in a Sum
  std::size_t size() const { return summands_.size(); }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{ \\bigl(";
    std::size_t counter = 0;
    for (const auto &i : summands()) {
      const auto i_is_product = i->is<Product>();
      if (!i_is_product) {
        result += (counter == 0) ? i->to_latex() : (L" + " + i->to_latex());
      } else {  // i_is_product
        const auto i_prod = i->as<Product>();
        const auto scalar = i_prod.scalar();
        if (scalar.real() < 0 || (scalar.real() == 0 && scalar.imag() < 0)) {
          result += L" - " + i_prod.to_latex(true);
        } else {
          result += (counter == 0) ? i->to_latex() : (L" + " + i->to_latex());
        }
      }
      ++counter;
    }
    result += L"\\bigr) }";
    return result;
  }

  std::wstring to_wolfram() const override {
    std::wstring result;
    result = L"Plus[";
    std::size_t counter = 0;
    for (const auto &i : summands()) {
      result += i->to_wolfram();
      ++counter;
      if (counter != summands().size()) result += L",";
    }
    result += L"]";
    return result;
  }

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<Sum>();
  };

  ExprPtr clone() const override {
    auto cloned_summands =
        summands() | ranges::views::transform(
                         [](const ExprPtr &ptr) { return ptr->clone(); });
    return ex<Sum>(ranges::begin(cloned_summands),
                   ranges::end(cloned_summands));
  }

  /// @brief adjoint of a Sum is a sum of adjoints of its factors
  virtual void adjoint() override;

  virtual Expr &operator+=(const Expr &that) override {
    this->append(const_cast<Expr &>(that).shared_from_this());
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    if (that.is<Constant>())
      this->append(ex<Constant>(-that.as<Constant>().value()));
    else
      this->append(ex<Product>(
          -1, ExprPtrList{const_cast<Expr &>(that).shared_from_this()}));
    return *this;
  }

 private:
  summands_type summands_{};
  std::optional<size_t>
      constant_summand_idx_{};  // points to the constant summand, if any; used
                                // to sum up constants in append/prepend

  cursor begin_cursor() override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() override {
    return summands_.empty() ? Expr::end_cursor()
                             : cursor{&summands_[0] + summands_.size()};
  };
  cursor begin_cursor() const override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() const override {
    return summands_.empty() ? Expr::end_cursor()
                             : cursor{&summands_[0] + summands_.size()};
  };

  /// @return the hash of this object
  /// @note this ensures that hash of a Sum of a single summand is
  /// identical to the hash of the summand itself.
  hash_type memoizing_hash() const override {
    auto compute_hash = [this]() {
      if (summands_.size() == 1)
        return summands_[0]->hash_value();
      else {
        auto deref_summands =
            summands() |
            ranges::views::transform(
                [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
        auto value = hash::range(ranges::begin(deref_summands),
                                 ranges::end(deref_summands));
        return value;
      }
    };

    if (!hash_value_) {
      hash_value_ = compute_hash();
    } else {
      SEQUANT_ASSERT(*hash_value_ == compute_hash());
    }

    return *hash_value_;
  }

  /// @param multipass if true, will do a multipass canonicalization, with extra
  /// cleanup pass after the deep canonization pass
  ExprPtr canonicalize_impl(bool multipass, CanonicalizeOptions opt);

  virtual ExprPtr canonicalize(
      CanonicalizeOptions opt =
          CanonicalizeOptions::default_options()) override {
    return canonicalize_impl(true, opt);
  }
  virtual ExprPtr rapid_canonicalize(
      CanonicalizeOptions opts =
          CanonicalizeOptions::default_options().copy_and_set(
              CanonicalizationMethod::Rapid)) override {
    SEQUANT_ASSERT(opts.method == CanonicalizationMethod::Rapid);
    return canonicalize_impl(false, opts);
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const Sum &>(that);
    if (summands().size() == that_cast.summands().size()) {
      if (this->empty()) return true;
      // compare hash values first
      if (this->hash_value() ==
          that.hash_value())  // hash values agree -> do full comparison
        return std::equal(begin_subexpr(), end_subexpr(), that.begin_subexpr(),
                          expr_ptr_comparer);
      else
        return false;
    } else
      return false;
  }
};  // class Sum

/// @brief utility for eagerly accumulating summands in a hash table
/// Intended to accumulate summands that are already in canonical form and
/// produces canonical Sum
class HashingAccumulator {
 public:
  /// @p summand expr to append to the sum
  /// @p flatten if true, and @p summand is a Sum, will flatten the sum
  HashingAccumulator &append(ExprPtr summand, bool flatten = true);

  SumPtr make_sum();

  SumPtr make_canonicalized_sum();

  /// @param canonicalize if true, will sort the summands to canonical order
  /// defined by ExprPtr::operator<
  /// @return summands as a Sum (if have more than 1 summand), Constant (if have
  /// zero summands), or the lone summand itself
  ExprPtr make_expr(bool canonicalize = true);

  bool empty() const { return summands_.empty(); }

 private:
  /// @brief Common implementation for make_sum and make_canonicalized_sum
  /// @param canonicalize if true, sort the summands by hash value
  SumPtr make_sum_impl(bool canonicalize);

  container::unordered_set<ExprPtr, sequant::hash::_<ExprPtr>, proportional_to>
      summands_;
};

struct TransformSumExprOptions {
  bool canonicalize = true;
  bool flatten = true;
};

/// variant of sequant::transform_reduce for eager sum reduction of Expr's
/// @sa HashingAccumulator
template <typename SizedRange, typename UnaryMapOp>
  requires(meta::is_range_v<std::remove_cvref_t<SizedRange>>)
ExprPtr transform_sum_expr(SizedRange &&rng, const UnaryMapOp &map,
                           const TransformSumExprOptions &options = {}) {
  HashingAccumulator result_acc;
  std::mutex result_mtx;  // serializes updates of result

  auto task = [&result_acc, &result_mtx, &map,
               canonicalize = options.canonicalize,
               flatten = options.flatten](const ExprPtr &input) {
    auto task_result = map(input);
    if (task_result) {
      if (canonicalize) {
        auto bp = task_result->canonicalize();
        if (bp) {
          task_result = bp * task_result;
        }
      }

      std::scoped_lock<std::mutex> lock(result_mtx);
      result_acc.append(task_result, flatten);
    }
  };
  sequant::for_each(std::forward<SizedRange>(rng), task);
  return result_acc.make_expr(options.canonicalize);
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_SUM_HPP
