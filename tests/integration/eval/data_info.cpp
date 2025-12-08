//
// Created by Bimal Gaudel on 7/14/21.
//

#include "data_info.hpp"

#include <SeQuant/core/utility/macros.hpp>

#include <fstream>
#include <sstream>

namespace sequant::eval {

size_t DataInfo::nocc() const { return nocc_; }

size_t DataInfo::nvirt() const { return nvirt_; }

std::string_view DataInfo::fock_file() const { return fock_file_; }

std::string_view DataInfo::eri_file() const { return eri_file_; }

DataInfo::DataInfo(std::string_view fname1, std::string_view fname2) {
  //-----------------------------------------------------------------
  // To read the header line (the first line) from each file
  //-----------------------------------------------------------------
  struct header_line {
    std::string_view fname;
    size_t nocc;
    size_t nvirt;
    size_t rank;
  };
  //
  auto read_header = [](std::string_view fname) -> header_line {
    auto ifs = std::ifstream{fname.data()};
    std::string header{};
    std::getline(ifs, header);
    auto hs = std::istringstream{header};

    auto result = header_line{};
    result.fname = fname;
    hs >> result.rank;
    hs >> result.nocc;
    hs >> result.nvirt;

    return result;
  };
  //-----------------------------------------------------------------

  // initially let's assume the first file is the Fock file
  // and the second file is the ERI file
  auto fock_header = read_header(fname1);
  auto eri_header = read_header(fname2);

  // the one with smaller rank is the Fock file
  if (fock_header.rank > eri_header.rank) std::swap(fock_header, eri_header);

  SEQUANT_ASSERT(fock_header.rank == 2 && "Fock tensor must have rank 2");
  SEQUANT_ASSERT(eri_header.rank == 4 && "ERI tensor must have rank 4");
  SEQUANT_ASSERT(fock_header.nocc == eri_header.nocc &&
                 fock_header.nvirt == eri_header.nvirt &&
                 "Number of occupied and virtual orbitals"
                 " in Fock and ERI tensor do not match.");

  fock_file_ = fock_header.fname;
  eri_file_ = eri_header.fname;
  nocc_ = fock_header.nocc;    // or eri_header.nocc
  nvirt_ = fock_header.nvirt;  // or eri_header.nvirt
}

}  // namespace sequant::eval
