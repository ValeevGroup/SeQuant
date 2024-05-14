//
// Created by Conner Masteran on 5/1/24.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;
  IndexSpaceRegistry sample_ISR;

  IndexSpace x_space(L"x", 0b001);
  sample_ISR.add(x_space);

  IndexSpace y_space(L"y", 0b010);
  sample_ISR.add(y_space);

  IndexSpace z_space(L"z", 0b100);
  sample_ISR.add(z_space);

  IndexSpace xy_space(L"xy", 0b011);
  sample_ISR.add(xy_space);

  IndexSpace yz_space(L"yz", 0b110);
  sample_ISR.add(yz_space);

  IndexSpace xyz_space(L"xyz", 0b111);
  sample_ISR.add(xyz_space);

  assert(sample_ISR.unIon(z_space, y_space) == yz_space);
  assert(sample_ISR.intersection(xy_space, x_space) == x_space);

  Context new_cxt(sample_ISR, Vacuum::SingleProduct);
  set_default_context(new_cxt);
}
