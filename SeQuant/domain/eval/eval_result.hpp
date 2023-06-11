#ifndef SEQUANT_EVAL_RESULT_HPP
#define SEQUANT_EVAL_RESULT_HPP

#include <tiledarray.h>
#include <any>
#include <memory>
#include <utility>

namespace sequant::eval {

namespace {
std::logic_error invalid_operand(
    std::string_view msg = "Invalid operand for binary op") noexcept {
  return std::logic_error{msg.data()};
}
}  // namespace

struct EvalResult;

using ERPtr = std::shared_ptr<EvalResult>;

template <typename T, typename... Args>
ERPtr eval_result(Args&&... args) noexcept {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

class EvalResult {
 public:
  using id_t = size_t;

  virtual ~EvalResult() = default;

  template <typename T>
  [[nodiscard]] bool is() const noexcept {
    return this->type_id() == id_for_type<std::decay_t<T>>();
  }

  template <typename T>
  [[nodiscard]] T const& as() const {
    assert(this->is<std::decay_t<T>()>());
    return static_cast<T const&>(*this);
  }

  [[nodiscard]] virtual ERPtr sum(EvalResult const&,
                                  std::array<std::any, 3> const&) const = 0;

  [[nodiscard]] virtual ERPtr prod(EvalResult const&,
                                   std::array<std::any, 3> const&) const = 0;

  [[nodiscard]] bool has_value() const noexcept;

  template <typename... Args>
  void set(Args&&... args) noexcept {
    value_.emplace(std::forward<Args...>(args...));
  }

  template <typename T>
  [[nodiscard]] T const& get() const {
    auto ptr = std::any_cast<T>(&value_);
    return *ptr;
  }

 protected:
  template <typename... Args>
  explicit EvalResult(Args&&... args) noexcept
      : value_{std::make_any<Args...>(args...)} {}

  [[nodiscard]] virtual id_t type_id() const noexcept = 0;

  template <typename T>
  [[nodiscard]] static id_t id_for_type() noexcept {
    static id_t id = next_id();
    return id;
  }

 private:
  std::any value_;

  [[nodiscard]] static id_t next_id() noexcept;
};

///
/// T is numeric type such as double
///
template <typename T>
class EvalConstant final : public EvalResult {
 public:
  using EvalResult::id_t;

  explicit EvalConstant(T v) : EvalResult(v) {}

  [[nodiscard]] T value() const noexcept { return get<T>(); }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const&) const override {
    if (other.is<EvalConstant>()) {
      auto const& o = other.as<EvalConstant>();
      return eval_result<EvalConstant>(value() + o.value());
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ERPtr prod(
      EvalResult const& other,
      std::array<std::any, 3> const& empty) const override {
    if (other.is<EvalConstant>()) {
      auto const& o = other.as<EvalConstant>();
      return eval_result<EvalConstant>(value() * o.value());
    } else {
      return other.prod(*this, empty);
    }
  }

 private:
  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalConstant>();
  }
};

template <typename T>
class EvalTensorTA final : public EvalResult {
 public:
  using EvalResult::id_t;
  using numeric_type = typename T::numeric_type;

  explicit EvalTensorTA(T const& arr) : EvalResult(arr) {}

 private:
  struct Annot {
    explicit Annot(std::array<std::any, 3> const& annot)
        : lannot(*std::any_cast<std::string>(&annot[0])),
          rannot(*std::any_cast<std::string>(&annot[1])),
          this_annot(*std::any_cast<std::string>(&annot[2])) {}
    std::string const& lannot;
    std::string const& rannot;
    std::string const& this_annot;
  };

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalTensorTA>();
  }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const& annot) const override {
    assert(other.is<EvalTensorTA<T>>());
    auto const a = Annot{annot};

    T result;
    result(a.this_annot) = get<T>()(a.lannot) + other.get<T>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<EvalTensorTA<T>>(std::move(result));
  }

  [[nodiscard]] ERPtr prod(
      EvalResult const& other,
      std::array<std::any, 3> const& annot) const override {
    if (other.is<EvalConstant<numeric_type>>()) {
      auto const a = Annot{annot};

      auto result = get<T>();
      result(a.this_annot) = other.get<double>() * result(a.lannot);
      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<EvalTensorTA<T>>(std::move(result));
    }

    auto const a = Annot{annot};

    if (a.this_annot.empty())
      return eval_result<EvalConstant<numeric_type>>(
          TA::dot(get<T>()(a.lannot), other.get<T>()(a.rannot)));

    T result;
    result(a.this_annot) = get<T>()(a.lannot) * other.get<T>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<EvalTensorTA<T>>(std::move(result));
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_RESULT_HPP
