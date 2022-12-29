//
// Created by Eduard Valeyev on 12/29/22.
//

#include "SeQuant/core/wstring.hpp"

#include <boost/locale/encoding_utf.hpp>

namespace sequant {

std::string to_string(const std::wstring& wstr_utf8) {
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(wstr_utf8.c_str(),
                          wstr_utf8.c_str() + wstr_utf8.size());
}

}  // namespace sequant
