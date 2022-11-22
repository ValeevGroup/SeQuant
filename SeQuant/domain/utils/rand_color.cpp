#include "rand_color.hpp"

#include <array>
#include <ctime>
#include <random>

namespace sequant::utils {

RandColor::RandColor()
    : randEngine{[]() {
        std::random_device seeder;
        const auto seed = seeder.entropy() ? seeder() : std::time(nullptr);
        return static_cast<std::mt19937_64::result_type>(seed);
      }()} {}

RandColor::RandColor(int seed)
    : randEngine{static_cast<std::mt19937_64::result_type>(seed)} {}

std::array<size_t, 3> RandColor::rand_rgb(double sat, double brit) {
  return RandColor::hsv_to_rgb(rand_hue(), sat, brit);
}

double RandColor::rand_hue() {
  auto attempt = [this]() {
    auto hue = GOLDEN_RATIO_CONJ + uniRealDist(randEngine);
    return hue_cache_.emplace(hue > 1 ? hue - 1 : hue);
  };

  auto result = attempt();
  while (!result.second) {
    // hue already exists in the cache registry
    result = attempt();
  }
  return *result.first;
}

std::array<size_t, 3> RandColor::hsv_to_rgb(double h, double s, double v) {
  // https://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
  size_t h_i = (size_t)(h * 6);
  double f = h * 6 - h_i;
  double p = v * (1 - s);
  double q = v * (1 - f * s);
  double t = v * (1 - (1 - f) * s);

  double r, g, b;
  r = g = b = -1.;

  if (h_i == 0) {
    r = v;
    g = t;
    b = p;
  } else if (h_i == 1) {
    r = q;
    g = v;
    b = p;
  } else if (h_i == 2) {
    r = p;
    g = v;
    b = t;
  } else if (h_i == 3) {
    r = p;
    g = q;
    b = v;
  } else if (h_i == 4) {
    r = t;
    g = p;
    b = v;
  } else if (h_i == 5) {
    r = v;
    g = p;
    b = q;
  }

  return std::array<size_t, 3>{
      (size_t)(256 * r),  //
      (size_t)(256 * g),  //
      (size_t)(256 * b)   //
  };
};

}  // namespace sequant::utils
