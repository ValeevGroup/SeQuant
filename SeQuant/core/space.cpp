//
// Created by Eduard Valeyev on 3/27/18.
//

#include "space.hpp"

sequant::container::map<sequant::IndexSpace::Attr, std::wstring>
    sequant::IndexSpace::attr2basekey_{};
sequant::container::map<std::wstring, sequant::IndexSpace::Attr,
                        sequant::IndexSpace::KeyCompare>
    sequant::IndexSpace::key2attr_{};
sequant::container::map<sequant::IndexSpace::Attr, sequant::IndexSpace>
    sequant::IndexSpace::instances_{};
sequant::IndexSpace sequant::IndexSpace::null_instance_{
    sequant::IndexSpace::Attr::null()};

namespace sequant {

IndexSpace::Type IndexSpace::nulltype = Type{0};
IndexSpace::Type IndexSpace::nonnulltype = Type{0x7fffffff};
IndexSpace::Type IndexSpace::frozen_occupied = Type{0b000001};
IndexSpace::Type IndexSpace::inactive_occupied = Type{0b000010};
IndexSpace::Type IndexSpace::active_occupied = Type{0b000100};
IndexSpace::Type IndexSpace::occupied = Type{0b000111};
IndexSpace::Type IndexSpace::active_unoccupied = Type{0b001000};
IndexSpace::Type IndexSpace::inactive_unoccupied = Type{0b010000};
IndexSpace::Type IndexSpace::unoccupied = Type{0b011000};
IndexSpace::Type IndexSpace::all_active = Type{0b001100};
IndexSpace::Type IndexSpace::all = Type{0b011111};
IndexSpace::Type IndexSpace::other_unoccupied = Type{0b100000};
IndexSpace::Type IndexSpace::complete_unoccupied = Type{0b111000};
IndexSpace::Type IndexSpace::complete = Type{0b111111};

IndexSpace::QuantumNumbers IndexSpace::nullqns =
    IndexSpace::QuantumNumbers{0b000000};  //!< no quantum numbers
IndexSpace::QuantumNumbers IndexSpace::alpha =
    IndexSpace::QuantumNumbers{0b000001};  //!< spin-up
IndexSpace::QuantumNumbers IndexSpace::beta =
    IndexSpace::QuantumNumbers{0b000010};  //!< spin-down

std::wstring to_wolfram(const IndexSpace& space) {
  std::wstring result = L"particleSpace[";

  // this is a hack due to partial representation of spaces in SeQuant
  if (space.type() == IndexSpace::active_occupied)
    result += L"occupied";
  else if (space.type() == IndexSpace::active_unoccupied)
    result += L"virtual";
  else if (space.type() == IndexSpace::all)
    result += L"occupied,virtual";
  else if (space.type() == IndexSpace::other_unoccupied)
    result += L"othervirtual";
  else if (space.type() == IndexSpace::complete_unoccupied)
    result += L"virtual,othervirtual";
  else if (space.type() == IndexSpace::complete)
    result += L"occupied,virtual,othervirtual";
  else
    throw std::invalid_argument(
        "to_wolfram(IndexSpace) received a nonstandard space");

  result += L"]";

  return result;
}

}  // namespace sequant
