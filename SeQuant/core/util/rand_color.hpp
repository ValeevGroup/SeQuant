#ifndef SEQUANT_UTILS_RAND_COLOR_HPP
#define SEQUANT_UTILS_RAND_COLOR_HPP

#include <array>
#include <random>
#include <set>

namespace sequant {

/**
 * Random color generator.
 *
 * Reference:
 * https://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
 *
 * @author Bimal Gaudel
 * @version 29 Sep 2020
 */

class RandColor {
 private:
  /// Used for better color distribution
  static constexpr double GOLDEN_RATIO_CONJ = 0.618033988749895;

  std::uniform_real_distribution<double> uniRealDist;

  std::mt19937_64 randEngine;

  /// Do not return the same hue again from this object.
  std::set<double> hue_cache_;

 public:
  RandColor();

  explicit RandColor(int seed);

  /**
   * Get a random color RGB hexcode for a given saturation level
   * and a brightness level.
   *
   * @param sat Saturation level.
   * @param brit Brightness level.
   *
   * @return std::array<size_t, 3> for red, blue and green
   *                               levels in the range [0,255]
   */
  std::array<size_t, 3> rand_rgb(double sat, double brit);

  /**
   * Get a random hue.
   *
   * @return A double in the range (0.0, 1.0).
   */
  double rand_hue();

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

}  // namespace sequant

#endif  // SEQUANT_UTILS_RAND_COLOR_HPP
