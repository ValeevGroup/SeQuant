/*
 *  Copyright (C) 2004-2017 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

#include "hartree-fock.h"

// this reads the geometry in the standard xyz format supported by most chemistry software
std::vector<Atom> read_dotxyz(std::istream& is)
{
  // line 1 = # of atoms
  size_t natom;
  is >> natom;
  // read off the rest of line 1 and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // line 2 = comment (possibly empty)
  std::string comment;
  std::getline(is, comment);

  std::vector<Atom> atoms(natom);
  for (auto i = 0; i < natom; i++) {
    std::string element_label;
    double x, y, z;
    is >> element_label >> x >> y >> z;

    // .xyz files report element labels, hence convert to atomic numbers
    int Z;
    if (element_label == "H")
      Z = 1;
    else if (element_label == "C")
      Z = 6;
    else if (element_label == "N")
      Z = 7;
    else if (element_label == "O")
      Z = 8;
    else if (element_label == "F")
      Z = 9;
    else if (element_label == "S")
      Z = 16;
    else if (element_label == "Cl")
      Z = 17;
    else {
      std::cerr << "read_dotxyz: element label \"" << element_label << "\" is not recognized" << std::endl;
      throw "Did not recognize element label in .xyz file";
    }

    atoms[i].atomic_number = Z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value
    /* const auto angstrom_to_bohr = 1 / 0.5291772083; // 2010 CODATA value */
    atoms[i].x = x * angstrom_to_bohr;
    atoms[i].y = y * angstrom_to_bohr;
    atoms[i].z = z * angstrom_to_bohr;
  }

  return atoms;
}

std::vector<Atom> read_geometry(const std::string& filename)
{

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms)
{

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
            {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
            {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
            {0, false, {0.1543289707, 0.5353281424, 0.4446345420}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
            );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
            {71.616837000, 13.045096000, 3.530512200},
            {
            {0, false, {0.15432897, 0.53532814, 0.44463454}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {2.941249400, 0.683483100, 0.222289900},
            {
            {0, false, {-0.09996723, 0.39951283, 0.70011547}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {2.941249400, 0.683483100, 0.222289900},
            { // contraction 0: p shell (l=1), spherical=false
            {1, false, {0.15591627, 0.60768372, 0.39195739}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
            {99.106169000, 18.052312000, 4.885660200},
            {
            {0, false, {0.15432897, 0.53532814, 0.44463454}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {3.780455900, 0.878496600, 0.285714400},
            {
            {0, false, {-0.09996723, 0.39951283, 0.70011547}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {3.780455900, 0.878496600, 0.285714400},
            { // contraction 0: p shell (l=1), spherical=false
            {1, false, {0.15591627, 0.60768372, 0.39195739}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
            {130.709320000, 23.808861000, 6.443608300},
            {
            {0, false, {0.1543289687, 0.5353281356, 0.4446345363}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {5.033151300, 1.169596100, 0.380389000},
            {
            {0, false, {-0.0999672287, 0.3995128246, 0.7001154606}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        shells.push_back(
            {
            {5.033151300, 1.169596100, 0.380389000},
            { // contraction 0: p shell (l=1), spherical=false
            {1, false, {0.1559162685, 0.6076837141, 0.3919573862}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
            );
        break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  return shells;
}

size_t nbasis(const std::vector<libint2::Shell>& shells)
{
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells)
{
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& shells)
{
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells)
{
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

// computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
// in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<Atom>& atoms)
{

  // compute number of atomic orbitals
  size_t nao = 0;
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) // H, He
      nao += 1;
    else if (Z <= 10) // Li - Ne
      nao += 5;
    else
      throw "SOAD with Z > 10 is not yet supported";
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0; // first AO of this atom
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) { // H, He
      D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
      ao_offset += 1;
    }
    else if (Z <= 10) {
      D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
      D(ao_offset+1, ao_offset+1) = (Z == 3) ? 1 : 2; // Li? only 1 electron in 2s, else 2 electrons
      // smear the remaining electrons in 2p orbitals
      const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4)/3 : 0;
      for(auto xyz=0; xyz!=3; ++xyz)
        D(ao_offset+2+xyz, ao_offset+2+xyz) = num_electrons_per_2p;
      ao_offset += 5;
    }
  }

  return D * 0.5; // we use densities normalized to # of electrons/2
}

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
    libint2::Operator obtype,
    const std::vector<Atom>& atoms)
{
  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = nbasis(shells);
  Matrix result(n,n);

  // construct the overlap integrals engine
  Engine engine(obtype, max_nprim(shells), max_l(shells), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function(shells);

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;
}

Matrix compute_2body_fock_simple(const std::vector<libint2::Shell>& shells,
    const Matrix& D)
{

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(auto s4=0; s4!=shells.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * buf_1234[f1234];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
          const auto* buf_1324 = buf[0];

          for(auto f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) -= D(bf3,bf4) * buf_1324[f1324];
                }
              }
            }
          }

        }
      }
    }
  }

  return G;
}

Matrix compute_2body_fock(const std::vector<libint2::Shell>& shells,
    const Matrix& D)
{

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  std::chrono::duration<double> time_elapsed = std::chrono::duration<double>::zero();

  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  const auto& buf = engine.results();

  // The problem with the simple Fock builder is that permutational symmetries of the Fock,
  // density, and two-electron integrals are not taken into account to reduce the cost.
  // To make the simple Fock builder efficient we must rearrange our computation.
  // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
  // hence we must minimize the number of computed integrals by taking advantage of their permutational
  // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
  // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
  //
  // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
  //
  // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
  // electron 1, and c and d -- to electron 2).
  //
  // It is easy to verify that the following set of nested loops produces a permutationally-unique
  // set of integrals:
  // foreach a = 0 .. n-1
  //   foreach b = 0 .. a
  //     foreach c = 0 .. a
  //       foreach d = 0 .. (a == c ? b : c)
  //         compute (ab|cd)
  //
  // The only complication is that we must compute integrals over shells. But it's not that complicated ...
  //
  // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
  // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3<=s1; ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for(auto s4=0; s4<=s4_max; ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          const auto tstart = std::chrono::high_resolution_clock::now();

          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          const auto tstop = std::chrono::high_resolution_clock::now();
          time_elapsed += tstop - tstart;

          // ANSWER
          // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
          //    F(a,b) += (ab|cd) * D(c,d)
          //    F(c,d) += (ab|cd) * D(a,b)
          //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
          //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
          //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
          //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
          // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
          //    i.e. the number of the integrals/sets equivalent to it
          // 3) the end result must be symmetrized
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                  G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                  G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                  G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                  G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                  G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}

#ifdef SEQUANT_HAS_BTAS
btas::Tensor<double> compute_mo_ints(const Matrix& coff_mat,
    const std::vector<libint2::Shell>& shells) {
  // returns 2e integrals on MO basis in physicist's notation

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;
  using BTensor = btas::Tensor<double>;

  const auto nao = nbasis(shells);

  // *************************************
  // first calculate integrals on AO basis
  // *************************************

  BTensor ints_ao(nao, nao, nao, nao);

  libint2::initialize();

  // construct the electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(auto s4=0; s4!=shells.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  ints_ao(bf1, bf2, bf3, bf4) = buf_1234[f1234];
                }
              }
            }
          }
        }
      }
    }
  } // computed ints_ao on AO basis
  libint2::finalize(); // done with libint

  // *********************
  // transform to MO basis
  // *********************
  //

  // See 'The Smarter Algorithm' at http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4

  BTensor coff_tensor(nao, nao);
  std::copy(coff_mat.data(), coff_mat.data()+coff_mat.size(), coff_tensor.begin());

  //first bracket (from right to left in the equation linked above)
  BTensor s_contract(nao, nao, nao, nao);
  btas::contract(1.0, ints_ao, {1, 2, 3, 4}, coff_tensor, {4, 5}, 0.0, s_contract, {5, 1, 2, 3});

  //second bracket
  BTensor rs_contract(nao, nao, nao, nao);
  btas::contract(1.0,     s_contract, {1, 2, 3, 4}, coff_tensor, {4, 5}, 0.0, rs_contract, {5, 1, 2, 3});

  //third bracket
  BTensor qrs_contract(nao, nao, nao, nao);
  btas::contract(1.0,    rs_contract, {1, 2, 3, 4}, coff_tensor, {4, 5}, 0.0, qrs_contract, {5, 1, 2, 3});

  //final loop
  BTensor pqrs_contract(nao, nao, nao, nao);
  btas::contract(1.0,   qrs_contract, {1, 2, 3, 4}, coff_tensor, {4, 5}, 0.0, pqrs_contract, {5, 1, 2, 3});

  // into Physicist's notation
  return BTensor(btas::permute(pqrs_contract, {0, 2, 1, 3}));
}
#endif

#ifdef SEQUANT_HAS_TILEDARRAY
TA::TArrayD compute_mo_ints(const Matrix& coff_mat,
                     const std::vector<libint2::Shell>& shells,
                     madness::World& world) {
  // returns 2e integrals on MO basis in physicist's notation

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto nao = nbasis(shells);

  // *************************************
  // first calculate integrals on AO basis
  // *************************************

  auto fill_eri_tensor = [&](const TA::Range& range){

    // allocate the tensor to return
    TA::Tensor<double> tile(range);

    libint2::initialize();
    // construct the electron repulsion integrals engine
    Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

    auto shell2bf = map_shell_to_basis_function(shells);

    // buf[0] points to the target shell set after every call  to engine.compute()
    const auto& buf = engine.results();

    // loop over shell pairs of the Fock matrix, {s1,s2}
    // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
    for(auto s1=0; s1!=shells.size(); ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();

      for(auto s2=0; s2!=shells.size(); ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        // loop over shell pairs of the density matrix, {s3,s4}
        // again symmetry is not used for simplicity
        for(auto s3=0; s3!=shells.size(); ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          for(auto s4=0; s4!=shells.size(); ++s4) {

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
            engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

            // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
            // hence some manual labor here:
            // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
            // and 2) add contribution from each integral
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;
                    tile(bf1, bf2, bf3, bf4) = buf_1234[f1234];
                  }
                }
              }
            }
          }
        }
      }
    } // computed ints on AO basis
    libint2::finalize(); // done with libint
    return tile;
  };

  std::vector<size_t> tile_boundaries = {0, nao};
  std::vector<TA::TiledRange1> ranges(
      4, TA::TiledRange1(tile_boundaries.begin(), tile_boundaries.end()));

  // create ints_ao
  TA::TiledRange trange_ao(ranges.begin(), ranges.end());
  TA::TArrayD    ints_ao(world, trange_ao);
  auto tile = ints_ao.world().taskq.add(fill_eri_tensor,
                                ints_ao.trange().make_tile_range(0));
  *(ints_ao.begin()) = tile;

  // *********************
  // transform to MO basis
  // *********************
  //

  auto mat2tensor = [&](const TA::Range& range){
    TA::Tensor<double> tile(range);
    std::copy(coff_mat.data(),
        coff_mat.data()+coff_mat.size(),
        tile.begin());
    return tile;
  };

  std::vector<TA::TiledRange1> ranges2(
      2, TA::TiledRange1(tile_boundaries.begin(), tile_boundaries.end()));

  TA::TiledRange trange_ao2(ranges2.begin(), ranges2.end());

  // transform the coefficient Eigen::MatrixD to TA::Tensor<double> type
  TA::TArrayD coeff_tensor(world, trange_ao2);
  tile = coeff_tensor.world().taskq.add(mat2tensor,
                                      coeff_tensor.trange().make_tile_range(0));
  *(coeff_tensor.begin()) = tile;

  TA::TArrayD ints_mo(world, trange_ao);
  ints_mo("5,1,2,3") = ints_ao("1,2,3,4") * coeff_tensor("4,5"); // s contract
  ints_mo("5,1,2,3") = ints_mo("1,2,3,4") * coeff_tensor("4,5"); // rs contract
  ints_mo("5,1,2,3") = ints_mo("1,2,3,4") * coeff_tensor("4,5"); // qrs contract
  ints_mo("5,1,2,3") = ints_mo("1,2,3,4") * coeff_tensor("4,5"); // pqrs contract

  // into Physicist's notation
  ints_mo("1,2,3,4") = ints_mo("1,3,2,4");
  return ints_mo;
}
#endif
