//
// Created by Eduard Valeyev on 10/16/25.
//

#ifndef SEQUANT_CORE_UTILITY_EXCEPTION_HPP
#define SEQUANT_CORE_UTILITY_EXCEPTION_HPP

#include <exception>
#include <string>

namespace sequant {

/// basic SeQuant exception
/// @sa SEQUANT_ASSERT
class Exception : public std::exception {
 public:
  Exception(const std::string& str) : msg_(str) {}
  virtual const char* what() const noexcept { return msg_.data(); }

 private:
  std::string msg_;
};  // class Exception

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_EXCEPTION_HPP
