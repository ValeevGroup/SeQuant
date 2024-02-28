//
// Created by Eduard Valeyev on 3/27/18.
//

#include "space.hpp"


sequant::IndexSpace sequant::IndexSpace::null_instance_{
    sequant::IndexSpace::Attr::null()};

namespace sequant {

std::wstring to_wolfram(const IndexSpace& space) {
  std::wstring result = L"particleSpace[";

  // this is a hack due to partial representation of spaces in SeQuant
  if (space.get_base_key() == L"M")
    result += L"occupied";
  else if (space.get_base_key() == L"e")
    result += L"virtual";
  else if (space.get_base_key() == L"p")
    result += L"occupied,virtual";
  else if (space.get_base_key() == L"α'")
    result += L"othervirtual";
  else if (space.get_base_key() == L"α")
    result += L"virtual,othervirtual";
  else if (space.get_base_key() == L"κ")
    result += L"occupied,virtual,othervirtual";
  else
    throw std::invalid_argument(
        "to_wolfram(IndexSpace) received a nonstandard space");

  result += L"]";

  return result;
}

}  // namespace sequant
