//
// Created by Eduard Valeyev on 5/17/23.
//

#ifndef SEQUANT_CORE_COMPLEX_HPP
#define SEQUANT_CORE_COMPLEX_HPP

#include <string>

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/wolfram.hpp>

namespace sequant {

/// analog of std::complex<Real> for non-real numeric rings
template <typename T>
struct Complex {
  using value_type = T;

  T re = {}, im = {};

  constexpr Complex() = delete;
  template <typename A, typename B = A,
            typename = std::enable_if_t<
                std::is_arithmetic_v<A> && std::is_constructible_v<T, A> &&
                std::is_arithmetic_v<B> && std::is_constructible_v<T, B>>>
  constexpr Complex(A real, B imag = 0) : re{real}, im{0} {}
  constexpr Complex(T real, T imag = 0) : re(real), im(imag) {}
  constexpr Complex(const Complex&) = default;
  constexpr Complex(Complex&&) = default;
  constexpr Complex& operator=(const Complex&) = default;
  constexpr Complex& operator=(Complex&&) = default;

  constexpr const T& real() const { return re; }
  constexpr const T& imag() const { return im; }

  constexpr bool is_zero() const { return real() == 0 && imag() == 0; }
  constexpr bool is_identity() const { return real() == 1 && imag() == 0; }

  std::wstring to_latex() const {
    using ::sequant::to_latex;
    std::wstring result = L"{";
    result += to_latex(this->real());
    if (this->imag() > 0) {
      result += L" + i " + to_latex(this->imag());
    } else if (this->imag() < 0)
      result += L" - i " + to_latex(-this->imag());
    result += L"}";
    return result;
  }

  std::wstring to_wolfram() const {
    using ::sequant::to_wolfram;
    if (this->imag() == 0)
      return to_wolfram(this->real());
    else
      return std::wstring(L"Complex[") + to_wolfram(this->real()) + L"," +
             to_wolfram(this->imag()) + L"]";
  }

  std::size_t hash_value() const {
    auto v = hash::value(this->real());
    hash::combine(v, this->imag());
    return v;
  }

  constexpr Complex& operator+=(const T& scalar) {
    re += scalar;
    return *this;
  }
  constexpr Complex& operator-=(const T& scalar) {
    re -= scalar;
    return *this;
  }
  constexpr Complex& operator*=(const T& scalar) {
    re *= scalar;
    im *= scalar;
    return *this;
  }

  template <class X>
  constexpr Complex& operator+=(const Complex<X>& other) {
    re += other.real();
    im += other.imag();
    return *this;
  }
  template <class X>
  constexpr Complex& operator-=(const Complex<X>& other) {
    re -= other.real();
    im -= other.imag();
    return *this;
  }
  template <class X>
  constexpr Complex& operator*=(const Complex<X>& other) {
    *this = Complex{re * other.real() - im * other.imag(),
                    re * other.imag() + im * other.real()};
    return *this;
  }

};  // class Complex

template <class T>
constexpr Complex<T> operator+(const Complex<T>& val) {
  return val;
}

template <class T>
constexpr Complex<T> operator-(const Complex<T>& val) {
  return {-val.real(), -val.imag()};
}

template <class T>
constexpr Complex<T> operator+(const Complex<T>& val1, const Complex<T>& val2) {
  return {val1.real() + val2.real(), val1.imag() + val2.imag()};
}

template <class T>
constexpr Complex<T> operator-(const Complex<T>& val1, const Complex<T>& val2) {
  return {val1.real() - val2.real(), val1.imag() - val2.imag()};
}

template <class T>
constexpr Complex<T> operator*(const Complex<T>& val1, const Complex<T>& val2) {
  return {val1.real() * val2.real() - val1.imag() * val2.imag(),
          val1.real() * val2.imag() + val1.imag() * val2.real()};
}

template <typename T>
constexpr T norm_squared(const Complex<T>& val) {
  return val.real() * val.real() + val.imag() * val.imag();
}

template <class T>
constexpr Complex<T> operator/(const Complex<T>& val1, const Complex<T>& val2) {
  return val1 * conj(val2) * (1 / norm_squared(val2));
}

template <class U, class T>
constexpr Complex<T> operator*(const U& val1, const Complex<T>& val2) {
  return {val1 * val2.real(), val1 * val2.imag()};
}

template <class U, class T>
constexpr Complex<T> operator*(const Complex<T>& val2, const U& val1) {
  return {val1 * val2.real(), val1 * val2.imag()};
}

template <class U, class T>
constexpr Complex<T> operator/(const U& val1, const Complex<T>& val2) {
  return val1 * conj(val2) * (1 / norm_squared(val2));
}

template <class T>
constexpr bool operator==(const Complex<T>& val1, const Complex<T>& val2) {
  return val1.real() == val2.real() && val1.imag() == val2.imag();
}

template <class T>
constexpr bool operator!=(const Complex<T>& val1, const Complex<T>& val2) {
  return !(val1 == val2);
}

template <class T, class X,
          typename = std::enable_if_t<std::is_arithmetic_v<X>>>
constexpr bool operator==(const Complex<T>& val1, const X& val2) {
  return val1.real() == val2 && val1.imag() == 0;
}

template <class T, class X,
          typename = std::enable_if_t<std::is_arithmetic_v<X>>>
constexpr bool operator==(const X& val2, const Complex<T>& val1) {
  return val1.real() == val2 && val1.imag() == 0;
}

template <class T, class X,
          typename = std::enable_if_t<std::is_arithmetic_v<X>>>
constexpr bool operator!=(const Complex<T>& val1, const X& val2) {
  return !(val1 == val2);
}

template <class T, class X,
          typename = std::enable_if_t<std::is_arithmetic_v<X>>>
constexpr bool operator!=(const X& val2, const Complex<T>& val1) {
  return !(val1 == val2);
}

template <class T>
constexpr bool operator==(const Complex<T>& val1, const T& val2) {
  return val1.real() == val2 && val1.imag() == 0;
}

template <class T>
constexpr bool operator==(const T& val2, const Complex<T>& val1) {
  return val1.real() == val2 && val1.imag() == 0;
}

template <class T>
constexpr bool operator!=(const Complex<T>& val1, const T& val2) {
  return !(val1 == val2);
}

template <class T>
constexpr bool operator!=(const T& val2, const Complex<T>& val1) {
  return !(val1 == val2);
}

template <typename T>
constexpr auto abs(const Complex<T>& val) {
  if constexpr (std::is_floating_point_v<T>) {
    using std::sqrt;
    return sqrt(val.real() * val.real() + val.imag() * val.imag());
  } else {
    assert(val.real() == 0 || val.imag() == 0);
    return val.imag() == 0 ? abs(val.real()) : abs(val.imag());
  }
}

template <typename T>
constexpr auto conj(const Complex<T>& val) {
  return Complex<T>{val.real(), -val.imag()};
}

}  // namespace sequant

#endif  // SEQUANT_CORE_COMPLEX_HPP
