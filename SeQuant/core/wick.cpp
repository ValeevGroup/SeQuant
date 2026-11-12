//
// Created by Eduard Valeyev on 10/9/25.
//

#include <SeQuant/core/wick.hpp>

#include <SeQuant/core/wick.impl.hpp>

namespace sequant {

template class SEQUANT_SYMB_TEMPLATE_DEF_API
    WickTheorem<Statistics::BoseEinstein>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API
    WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant
