#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP

#include <SeQuant/core/math.hpp>

#include <libperm/Permutation.hpp>
#include <libperm/Rank.hpp>
#include <libperm/Utils.hpp>

#include <Eigen/Core>
#include <cstddef>
#include <optional>

namespace sequant {

// hardcoded coefficients for biorthogonalization
std::vector<sequant::rational> get_first_row_biorth_coeffs_rational(
    std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<sequant::rational>{ratio(1, 2)};

    case 2:
      return std::vector<sequant::rational>{ratio(1, 3), ratio(1, 6)};

    case 3:
      return std::vector<sequant::rational>{ratio(17, 120), ratio(-7, 120),
                                            ratio(-1, 120), ratio(-1, 120),
                                            ratio(-1, 120), ratio(-7, 120)};

    case 4:
      return std::vector<sequant::rational>{
          ratio(43, 840), ratio(-19, 1680), ratio(-19, 1680),
          ratio(-1, 105), ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(1, 120),    ratio(-1, 105),
          ratio(1, 120),  ratio(-1, 105),   ratio(-19, 1680),
          ratio(-1, 105), ratio(1, 120),    ratio(1, 120),
          ratio(13, 840), ratio(-1, 105),   ratio(-1, 105),
          ratio(1, 120),  ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(-19, 1680), ratio(1, 120)};

    case 5:
      return std::vector<sequant::rational>{
          ratio(59, 3780),   ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-31, 7560),  ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-23, 30240), ratio(19, 7560),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-5, 3024),   ratio(19, 7560),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-31, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-23, 30240),
          ratio(-13, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(37, 15120),
          ratio(19, 7560),   ratio(-23, 30240), ratio(19, 7560),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-23, 30240), ratio(-31, 7560),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(19, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-31, 7560),
          ratio(37, 15120),  ratio(-31, 7560),  ratio(37, 15120),
          ratio(-5, 3024),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(19, 7560),   ratio(37, 15120),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-13, 7560),
          ratio(-23, 30240), ratio(37, 15120),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-23, 30240),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-23, 30240),
          ratio(19, 7560),   ratio(-13, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(19, 7560),   ratio(-31, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-5, 3024),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-5, 3024),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-13, 7560),
          ratio(19, 7560),   ratio(19, 7560),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(-5, 3024),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(37, 15120),  ratio(-23, 30240)};

    case 6:
      return std::vector<sequant::rational>{
          ratio(43, 10395),   ratio(-13, 83160),  ratio(-13, 83160),
          ratio(-13, 83160),  ratio(-13, 83160),  ratio(-53, 41580),
          ratio(-13, 83160),  ratio(-13, 83160),  ratio(-47, 166320),
          ratio(-47, 166320), ratio(5, 8316),     ratio(41, 83160),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-13, 83160),
          ratio(-47, 166320), ratio(5, 8316),     ratio(41, 83160),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(-13, 83160),  ratio(5, 8316),     ratio(41, 83160),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(41, 83160),
          ratio(41, 83160),   ratio(-53, 41580),  ratio(-13, 83160),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(79, 166320),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(-13, 83160),  ratio(-13, 83160),  ratio(-47, 166320),
          ratio(5, 8316),     ratio(-47, 166320), ratio(41, 83160),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-47, 166320),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(17, 332640),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-1, 33264),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-1, 3780),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(17, 332640),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(-47, 166320),
          ratio(-13, 83160),  ratio(79, 166320),  ratio(-47, 166320),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(-47, 166320), ratio(-47, 166320), ratio(79, 166320),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(17, 332640),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-13, 83160),
          ratio(5, 8316),     ratio(-47, 166320), ratio(41, 83160),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-47, 166320),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-1, 3780),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-1, 3780),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(-47, 166320),
          ratio(-53, 41580),  ratio(-1, 3780),    ratio(-1, 3780),
          ratio(41, 83160),   ratio(41, 83160),   ratio(5, 8316),
          ratio(41, 83160),   ratio(41, 83160),   ratio(17, 332640),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(-47, 166320),
          ratio(41, 83160),   ratio(17, 332640),  ratio(41, 83160),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(41, 83160),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(-13, 83160),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(5, 8316),     ratio(-13, 83160),  ratio(41, 83160),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 3780),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 33264),
          ratio(79, 166320),  ratio(-1, 33264),   ratio(17, 332640),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-47, 166320),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(17, 332640),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(79, 166320),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(79, 166320),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(79, 166320),  ratio(-47, 166320), ratio(-1, 3780),
          ratio(-13, 83160),  ratio(-13, 83160),  ratio(5, 8316),
          ratio(-47, 166320), ratio(-47, 166320), ratio(41, 83160),
          ratio(5, 8316),     ratio(-47, 166320), ratio(5, 8316),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(-1, 3780),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 33264),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(17, 332640),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-1, 3780),
          ratio(17, 332640),  ratio(41, 83160),   ratio(-47, 166320),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-1, 33264),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(-1, 83160),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 33264),
          ratio(-1, 33264),   ratio(79, 166320),  ratio(17, 332640),
          ratio(5, 8316),     ratio(-47, 166320), ratio(5, 8316),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 3780),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 33264),
          ratio(-47, 166320), ratio(79, 166320),  ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 83160),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-1, 33264),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(-1, 83160),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 3780),
          ratio(17, 332640),  ratio(17, 332640),  ratio(79, 166320),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-1, 3780),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-47, 166320),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 33264),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(17, 332640),
          ratio(-53, 41580),  ratio(-1, 3780),    ratio(41, 83160),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(5, 8316),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(17, 332640),
          ratio(17, 332640),  ratio(-1, 83160),   ratio(-1, 33264),
          ratio(41, 83160),   ratio(17, 332640),  ratio(41, 83160),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(-47, 166320),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-1, 3780),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-47, 166320),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-47, 166320),
          ratio(-47, 166320), ratio(5, 8316),     ratio(-1, 3780),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-1, 3780),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-47, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 3780),
          ratio(17, 332640),  ratio(17, 332640),  ratio(79, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 83160),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 33264),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 33264),
          ratio(79, 166320),  ratio(-47, 166320), ratio(17, 332640),
          ratio(-53, 41580),  ratio(-1, 3780),    ratio(41, 83160),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(5, 8316),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-1, 3780),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(-53, 41580),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-13, 83160),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-1, 3780),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(5, 8316),
          ratio(-47, 166320), ratio(-13, 83160),  ratio(41, 83160),
          ratio(-13, 83160),  ratio(79, 166320),  ratio(-47, 166320),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 33264),
          ratio(79, 166320),  ratio(-1, 33264),   ratio(17, 332640),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(5, 8316),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(5, 8316),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(41, 83160),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-1, 3780),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(-47, 166320),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(41, 83160),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(5, 8316),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-1, 3780),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-47, 166320),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(41, 83160),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(-1, 3780),    ratio(-1, 83160),   ratio(-1, 33264),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-47, 166320),
          ratio(-1, 33264),   ratio(5, 8316),     ratio(-1, 3780),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-1, 3780),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-47, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 3780),
          ratio(-1, 83160),   ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(79, 166320),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(79, 166320),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(17, 332640),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(17, 332640),
          ratio(-1, 83160),   ratio(17, 332640),  ratio(-1, 33264),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(41, 83160),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-47, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(79, 166320),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 33264),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(17, 332640),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-1, 33264),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(-1, 83160),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 33264),
          ratio(-1, 33264),   ratio(79, 166320),  ratio(17, 332640),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-47, 166320),
          ratio(-1, 33264),   ratio(79, 166320),  ratio(17, 332640),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-47, 166320),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-1, 3780),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(-1, 83160),   ratio(-1, 3780),    ratio(-1, 33264),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(-1, 3780),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(5, 8316),
          ratio(41, 83160),   ratio(41, 83160),   ratio(17, 332640),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 3780),
          ratio(17, 332640),  ratio(-1, 83160),   ratio(-1, 33264),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(17, 332640),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-47, 166320),
          ratio(5, 8316),     ratio(-47, 166320), ratio(-1, 33264),
          ratio(-47, 166320), ratio(5, 8316),     ratio(-1, 3780),
          ratio(-53, 41580),  ratio(41, 83160),   ratio(-1, 3780),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(5, 8316),
          ratio(-53, 41580),  ratio(-53, 41580),  ratio(41, 83160),
          ratio(41, 83160),   ratio(41, 83160),   ratio(-13, 83160),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(41, 83160),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-47, 166320),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(17, 332640),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-47, 166320),
          ratio(-13, 83160),  ratio(5, 8316),     ratio(-47, 166320),
          ratio(-47, 166320), ratio(-13, 83160),  ratio(41, 83160),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(79, 166320),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(41, 83160),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-47, 166320),
          ratio(-1, 3780),    ratio(-1, 83160),   ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(-1, 33264),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 33264),
          ratio(79, 166320),  ratio(-47, 166320), ratio(17, 332640),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(79, 166320),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(5, 8316),     ratio(5, 8316),     ratio(-1, 33264),
          ratio(-47, 166320), ratio(-47, 166320), ratio(-1, 3780),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-47, 166320),
          ratio(79, 166320),  ratio(-1, 33264),   ratio(17, 332640),
          ratio(-13, 83160),  ratio(5, 8316),     ratio(-47, 166320),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(41, 83160),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(17, 332640),
          ratio(-1, 3780),    ratio(41, 83160),   ratio(-47, 166320),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 3780),
          ratio(-47, 166320), ratio(-47, 166320), ratio(79, 166320),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(17, 332640),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-47, 166320),
          ratio(79, 166320),  ratio(-1, 33264),   ratio(17, 332640),
          ratio(-47, 166320), ratio(79, 166320),  ratio(-1, 33264),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(17, 332640),
          ratio(-1, 3780),    ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(79, 166320),
          ratio(5, 8316),     ratio(-1, 33264),   ratio(-1, 33264),
          ratio(-1, 33264),   ratio(-1, 33264),   ratio(-1, 83160),
          ratio(5, 8316),     ratio(5, 8316),     ratio(-47, 166320),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 3780),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-47, 166320),
          ratio(-1, 33264),   ratio(79, 166320),  ratio(17, 332640),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(-1, 33264),
          ratio(-47, 166320), ratio(79, 166320),  ratio(17, 332640),
          ratio(-1, 3780),    ratio(-1, 83160),   ratio(17, 332640),
          ratio(17, 332640),  ratio(-1, 3780),    ratio(-1, 33264),
          ratio(-13, 83160),  ratio(-47, 166320), ratio(-47, 166320),
          ratio(79, 166320),  ratio(-47, 166320), ratio(-1, 3780),
          ratio(5, 8316),     ratio(5, 8316),     ratio(-47, 166320),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(-1, 3780),
          ratio(-13, 83160),  ratio(5, 8316),     ratio(-13, 83160),
          ratio(-47, 166320), ratio(-47, 166320), ratio(41, 83160),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(79, 166320),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(17, 332640),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-1, 3780),
          ratio(17, 332640),  ratio(41, 83160),   ratio(-47, 166320),
          ratio(41, 83160),   ratio(17, 332640),  ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(-1, 3780),    ratio(-1, 3780),    ratio(17, 332640),
          ratio(17, 332640),  ratio(17, 332640),  ratio(79, 166320),
          ratio(-1, 3780),    ratio(-1, 83160),   ratio(-1, 3780),
          ratio(17, 332640),  ratio(17, 332640),  ratio(-1, 33264),
          ratio(41, 83160),   ratio(-1, 3780),    ratio(-1, 3780),
          ratio(41, 83160),   ratio(17, 332640),  ratio(-47, 166320),
          ratio(-47, 166320), ratio(-1, 33264),   ratio(79, 166320),
          ratio(-1, 33264),   ratio(-47, 166320), ratio(17, 332640)};

    default:
      return std::vector<sequant::rational>{};
  }
}

// constructs the full coeff matrix from the first row
Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
biorth_coeffs_from_first_row_rational(
    const std::vector<sequant::rational>& first_row, std::size_t n_particles) {
  const auto n = first_row.size();
  Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic> M(n, n);

  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      // using (n - 1 - row) to reverse the row order
      perm::Permutation row_perm = perm::unrank(n - 1 - row, n_particles);
      perm::Permutation col_perm = perm::unrank(col, n_particles);

      // row_perm * col_perm
      col_perm->preMultiply(row_perm);

      std::size_t source_idx = perm::rank(col_perm, n_particles);
      M(row, col) = first_row[source_idx];
    }
  }
  return M;
}

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
get_hardcoded_biorth_coeffs_rational(std::size_t n_particles) {
  auto first_row = get_first_row_biorth_coeffs_rational(n_particles);
  if (!first_row.empty()) {
    return biorth_coeffs_from_first_row_rational(first_row, n_particles);
  }
  // empty matrix for unsupported ranks
  return Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>(0, 0);
}

}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
