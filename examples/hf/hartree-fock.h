#ifndef MATRICES_H
#define MATRICES_H
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <libint2.hpp>
#include <vector>

// BTAS
#ifdef SEQUANT_HAS_BTAS
#include <btas/btas.h>
#include <btas/tensor_func.h>
#endif
// TiledArray
#ifdef SEQUANT_HAS_TILEDARRAY
#include <tiledarray.h>
#endif

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;  // import dense, dynamically sized Matrix type from Eigen;
             // this is a matrix with row-major storage
             // (http://en.wikipedia.org/wiki/Row-major_order) to meet the
             // layout of the integrals returned by the Libint integral library

struct Atom {
  int atomic_number;
  double x, y, z;
};

std::vector<Atom> read_geometry(const std::string& filename);

std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms);

size_t nbasis(const std::vector<libint2::Shell>& shells);

std::vector<size_t> map_shell_to_basis_function(
    const std::vector<libint2::Shell>& shells);

size_t max_nprim(const std::vector<libint2::Shell>& shells);

int max_l(const std::vector<libint2::Shell>& shells);

Matrix compute_soad(const std::vector<Atom>& atoms);

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::Operator t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());

// simple-to-read, but inefficient Fock builder; computes ~16 times as many ints
// as possible
Matrix compute_2body_fock_simple(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);

// an efficient Fock builder; *integral-driven* hence computes
// permutationally-unique ints once
Matrix compute_2body_fock(const std::vector<libint2::Shell>& shells,
                          const Matrix& D);

#ifdef SEQUANT_HAS_BTAS
btas::Tensor<double> compute_mo_ints(const Matrix& coff_mat,
                                     const std::vector<libint2::Shell>& shells);
#endif
#ifdef SEQUANT_HAS_TILEDARRAY
TA::TArrayD compute_mo_ints(const Matrix&, const std::vector<libint2::Shell>&,
                            madness::World&);
#endif

#endif /* ifndef MATRICES_H */
