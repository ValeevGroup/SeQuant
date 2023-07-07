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
