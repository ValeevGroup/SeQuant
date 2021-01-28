#include "rand_color.hpp"

#include <array>
#include <ctime>
#include <random>
#include <sstream>
#include <string>

namespace sequant::domain::util {

rand_color::rand_color()
    : randEngine{[]() {
        std::random_device seeder;
        const auto seed = seeder.entropy() ? seeder() : std::time(nullptr);
        return static_cast<std::mt19937_64::result_type>(seed);
      }()} {}

std::wstring rand_color::rand_rgb(double saturation, double brightness) {
  std::wostringstream oss;
  oss << "#" << std::hex;
  for (auto cc : hsv_to_rgb(rand_hue(), saturation, brightness)) oss << cc;
  return oss.str();
}

double rand_color::rand_hue() {
  auto hue = GOLDEN_RATIO_CONJ + uniRealDist(randEngine);
  return hue > 1 ? hue - 1 : hue;
}

std::array<size_t, 3> rand_color::hsv_to_rgb(double h, double s, double v) {
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

}  // namespace sequant::domain::util
