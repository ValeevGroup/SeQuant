//
// Created by Bimal Gaudel on 7/13/21.
//

#ifndef SEQUANT_EVAL_DATA_INFO_HPP
#define SEQUANT_EVAL_DATA_INFO_HPP

#include <string>

namespace sequant::eval {

class DataInfo {
 private:
  size_t nocc_;

  size_t nvirt_;

  std::string fock_file_;

  std::string eri_file_;

 public:
  inline static size_t constexpr fock_rank = 2;

  inline static size_t constexpr eri_rank = 4;

  /// \param fname1 Input file name with Fock(or ERI) tensor data.
  /// \param fname2 Input file name with ERI(or Fock) tensor data.
  DataInfo(std::string_view fname1, std::string_view fname2);

  ///
  /// \return number of occupied orbitals.
  size_t nocc() const;

  ///
  /// \return number of virtual orbitals.
  size_t nvirt() const;

  ///
  /// \return Fock tensor input file name.
  std::string_view fock_file() const;

  ///
  /// \return ERI tensor input file name.
  std::string_view eri_file() const;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_DATA_INFO_HPP
