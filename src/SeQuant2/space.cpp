//
// Created by Eduard Valeyev on 3/27/18.
//

#include "space.hpp"

std::map<sequant2::IndexSpace::Attr, std::wstring> sequant2::IndexSpace::keys_{};
std::map<sequant2::IndexSpace::Attr, sequant2::IndexSpace> sequant2::IndexSpace::instances_{};
sequant2::IndexSpace sequant2::IndexSpace::null_instance_{sequant2::IndexSpace::Attr::null()};

namespace sequant2 {

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

IndexSpace::QuantumNumbers IndexSpace::nullqns = IndexSpace::QuantumNumbers{0b000000};  //!< no quantum numbers
IndexSpace::QuantumNumbers IndexSpace::alpha = IndexSpace::QuantumNumbers{0b000001};  //!< spin-up
IndexSpace::QuantumNumbers IndexSpace::beta = IndexSpace::QuantumNumbers{0b000010};  //!< spin-down

}