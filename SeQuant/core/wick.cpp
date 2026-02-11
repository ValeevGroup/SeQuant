//
// Created by Eduard Valeyev on 10/9/25.
//

#include <SeQuant/core/wick.hpp>

#include <SeQuant/core/wick.impl.hpp>

namespace sequant {

template class WickTheorem<Statistics::BoseEinstein>;
template class WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant
