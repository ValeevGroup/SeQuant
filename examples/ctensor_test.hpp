#include <clocale>
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <algorithm>
#include <btas/btas.h>
#include <btas/tensorview.h>
#include <btas/tensor_func.h>

#include <boost/numeric/interval.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>

#include "contractable_tensor.hpp"
/* #include "contract.hpp" */

using STensor = sequant::Tensor;
using CTensor = sequant::contractable::Tensor;
using BTensor = btas::Tensor<double>;

using std::wcout;
using std::endl;

void print_index(const index_vec& v, std::string&& s="") {
  if (! s.empty())
    wcout << endl << s.c_str() << ":" << endl;

  for (const auto& i: v)
    wcout << i.label() << "  ";

  /* wcout << endl; */
}

void print_norm(const BTensor& bt, std::string&& s="") {
  if (! s.empty())
    wcout << endl << s.c_str() << ": ";

  wcout << sqrt(btas::dot(bt, bt)) << endl;
}

void print_norm(const CTensor& ct, std::string&& s="") {
  if (! s.empty())
    wcout << endl << s.c_str() << ": ";

  wcout << sqrt(btas::dot(ct.bt(), ct.bt())) << endl;
}

void print_ctensor(const CTensor& ct, std::string&& s="") {
  if (! s.empty())
    wcout << endl << s.c_str() << ": ";
  wcout << ct.label();
  wcout << "__";
  print_index(ct.occs());
  wcout << "^^";
  print_index(ct.virts());
  wcout << endl;
}
