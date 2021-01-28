#ifndef SEQUANT_UTIL_RAND_COLOR_HPP
#define SEQUANT_UTIL_RAND_COLOR_HPP

#include <array>
#include <random>
#include <string>

namespace sequant::domain::util {

/**
 * Random color generator.
 *
 * Reference:
 * https://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
 *
 * @author Bimal Gaudel
 * @version 29 Sep 2020
 */
class rand_color {
 private:
  static constexpr double GOLDEN_RATIO_CONJ = 0.618033988749895;

  std::uniform_real_distribution<double> uniRealDist;

  std::mt19937_64 randEngine;

 public:
  rand_color();

  /**
   * Get a random color RGB hexcode for a given saturation level
   * and a brightness level.
   *
   * @param saturation Saturation level.
   * @param brightness Brightness level.
   *
   * @return String of the pattern '#<R><G><B>' where <X> is a hex number.
   */
  std::wstring rand_rgb(double saturation, double brightness);

  /**
   * Get a random hue.
   *
   * @return A double in the range (0.0, 1.0).
   */
  double rand_hue();

 private:
  /**
   * Convert hue, saturation, brightness value to hex number.
   *
   * @param hue Hue value.
   * @param saturation Saturation value.
   * @param brightness Brightness value.
   *
   * @return Size three array of non-negative integers that represent the RGB
   * color values in that order.
   */
  static std::array<size_t, 3> hsv_to_rgb(double hue, double saturation,
                                          double brightness);
};

}  // namespace sequant::domain::util

#endif  // SEQUANT_UTIL_RAND_COLOR_HPP
