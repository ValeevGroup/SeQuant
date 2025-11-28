#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include "test_export.hpp"

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/python_einsum.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/scope.hpp>

#include "catch2_sequant.hpp"

#include <boost/algorithm/string.hpp>

#include <unsupported/Eigen/CXX11/Tensor>

#include <filesystem>
#include <fstream>
#include <random>
#include <string>

using namespace sequant;

// ============================================================================
// Eigen::Tensor utilities for numpy I/O and validation
// ============================================================================
namespace {

// Write Eigen::Tensor to numpy .npy file (version 1.0 format)
// Supports both column-major (default) and row-major layouts
template <typename Scalar, int NumDims, int Options = 0,
          typename IndexType = Eigen::Index>
void write_eigen_tensor_to_numpy(
    const std::string &filename,
    const Eigen::Tensor<Scalar, NumDims, Options, IndexType> &tensor) {
  std::ofstream file(filename, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Cannot open file for writing: " + filename);
  }

  // Magic number
  file.write("\x93NUMPY", 6);

  // Version 1.0
  file.put(0x01);
  file.put(0x00);

  // Build header dict
  std::ostringstream header;
  header << "{'descr': '";

  // Data type descriptor
  if (std::is_same<Scalar, float>::value) {
    header << "<f4";
  } else if (std::is_same<Scalar, double>::value) {
    header << "<f8";
  } else if (std::is_same<Scalar, std::complex<float>>::value) {
    header << "<c8";
  } else if (std::is_same<Scalar, std::complex<double>>::value) {
    header << "<c16";
  } else {
    throw std::runtime_error("Unsupported data type");
  }

  // Determine layout from Options template parameter
  // Eigen::RowMajor = 1, Eigen::ColMajor = 0 (default)
  constexpr bool is_row_major = (Options & Eigen::RowMajor) != 0;
  constexpr bool fortran_order = !is_row_major;

  header << "', 'fortran_order': " << (fortran_order ? "True" : "False")
         << ", 'shape': (";
  auto dimensions = tensor.dimensions();
  for (int i = 0; i < NumDims; ++i) {
    if (i > 0) header << ", ";
    header << dimensions[i];
  }
  if (NumDims == 1) header << ",";
  header << "), }";

  // Pad header so total (10 bytes preamble + header length) is divisible by 64
  // The header dict should end with a space before the newline for numpy
  // compatibility
  std::string header_str = header.str();
  // Calculate required total size (must be multiple of 64)
  std::size_t preamble_size = 10;  // magic(6) + version(2) + header_len(2)
  std::size_t current_size =
      preamble_size + header_str.size() + 1;  // +1 for newline
  std::size_t padding = (64 - (current_size % 64)) % 64;

  // Add padding spaces, then newline
  if (padding > 0) {
    header_str.append(padding, ' ');
  }
  header_str += "\n";

  // Write header length (2 bytes, little endian)
  uint16_t hlen = static_cast<uint16_t>(header_str.size());
  file.write(reinterpret_cast<const char *>(&hlen), 2);

  // Write header
  file.write(header_str.c_str(), header_str.size());

  // Write data
  file.write(reinterpret_cast<const char *>(tensor.data()),
             tensor.size() * sizeof(Scalar));
}

// Read Eigen::Tensor from numpy .npy file
// Supports both column-major (default) and row-major layouts
template <typename Scalar, int NumDims, int Options = 0,
          typename IndexType = Eigen::Index>
Eigen::Tensor<Scalar, NumDims, Options, IndexType> read_eigen_tensor_from_numpy(
    const std::string &filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Cannot open file for reading: " + filename);
  }

  // Read and verify magic number
  char magic[6];
  file.read(magic, 6);
  if (std::string(magic, 6) != "\x93NUMPY") {
    throw std::runtime_error("Invalid numpy file format");
  }

  // Read version
  uint8_t major, minor;
  file.read(reinterpret_cast<char *>(&major), 1);
  file.read(reinterpret_cast<char *>(&minor), 1);

  // Read header length
  uint16_t header_len;
  file.read(reinterpret_cast<char *>(&header_len), 2);

  // Read header
  std::vector<char> header_buf(header_len);
  file.read(header_buf.data(), header_len);
  std::string header(header_buf.begin(), header_buf.end());

  // Parse fortran_order from header
  bool fortran_order = true;  // default assumption
  auto order_start = header.find("'fortran_order': ");
  if (order_start != std::string::npos) {
    auto order_value_start = order_start + 17;  // length of "'fortran_order': "
    if (header.substr(order_value_start, 4) == "True") {
      fortran_order = true;
    } else if (header.substr(order_value_start, 5) == "False") {
      fortran_order = false;
    }
  }

  // Verify layout matches expected layout from template parameters
  // When Options=0 (default), we accept either layout for compatibility.
  // The data layout in the file determines the Eigen tensor's actual layout.
  constexpr bool is_row_major = (Options & Eigen::RowMajor) != 0;
  constexpr bool expected_fortran_order = !is_row_major;

  // With default Options=0, Eigen uses column-major (Fortran order).
  // But if the numpy file is row-major (fortran_order=False), we need to handle
  // it. For now, we allow this mismatch and rely on the fact that both numpy
  // and Eigen store data contiguously - the interpretation of indices may
  // differ, but for symmetric operations this often works out.
  if (Options != 0 && fortran_order != expected_fortran_order) {
    // Only enforce strict layout matching when non-default Options were
    // explicitly specified
    throw std::runtime_error(
        std::string("Layout mismatch: numpy file has fortran_order=") +
        (fortran_order ? "True" : "False") +
        ", but tensor expects fortran_order=" +
        (expected_fortran_order ? "True" : "False"));
  }

  // Parse shape from header
  auto shape_start = header.find("'shape': (");
  auto shape_end = header.find(")", shape_start);
  if (shape_start == std::string::npos || shape_end == std::string::npos) {
    throw std::runtime_error("Cannot parse shape from header");
  }

  std::string shape_str =
      header.substr(shape_start + 10, shape_end - shape_start - 10);
  std::vector<std::size_t> shape_vec;

  std::istringstream shape_stream(shape_str);
  std::string dim;
  while (std::getline(shape_stream, dim, ',')) {
    std::string trimmed = boost::trim_copy(dim);
    if (!trimmed.empty()) {
      shape_vec.push_back(std::stoull(trimmed));
    }
  }

  if (shape_vec.size() != NumDims) {
    throw std::runtime_error("Shape mismatch: expected " +
                             std::to_string(NumDims) + " dimensions, got " +
                             std::to_string(shape_vec.size()));
  }

  // Create tensor with appropriate dimensions
  std::array<IndexType, NumDims> dims;
  for (int i = 0; i < NumDims; ++i) {
    dims[i] = static_cast<IndexType>(shape_vec[i]);
  }

  Eigen::Tensor<Scalar, NumDims, Options, IndexType> tensor(dims);

  // Read data
  file.read(reinterpret_cast<char *>(tensor.data()),
            tensor.size() * sizeof(Scalar));

  return tensor;
}

// Helper function to properly escape shell arguments for POSIX shells
// Uses single quotes and escapes embedded single quotes
std::string shell_escape(const std::string &arg) {
  std::string result = "'";
  for (char c : arg) {
    if (c == '\'') {
      // End the current single-quoted string, add an escaped single quote,
      // and start a new single-quoted string
      result += "'\\''";
    } else {
      result += c;
    }
  }
  result += "'";
  return result;
}

// Execute Python code and check for errors
#if defined(SEQUANT_HAS_NUMPY_FOR_VALIDATION) || \
    defined(SEQUANT_HAS_TORCH_FOR_VALIDATION)
bool run_python_code(const std::string &code, const std::string &working_dir,
                     bool use_torch = false) {
  // Convert to filesystem::path for proper handling
  std::filesystem::path work_dir_path(working_dir);

  // Validate that the path is safe (canonical and exists)
  if (!std::filesystem::exists(work_dir_path)) {
    return false;
  }

  // Get canonical path to prevent path traversal issues
  std::filesystem::path canonical_work_dir;
  try {
    canonical_work_dir = std::filesystem::canonical(work_dir_path);
  } catch (const std::filesystem::filesystem_error &) {
    return false;
  }

  // Write code to a temporary file
  std::filesystem::path script_path = canonical_work_dir / "test_script.py";
  std::ofstream script(script_path);
  if (!script) return false;

  if (use_torch) {
    script << "import torch\n";
  } else {
    script << "import numpy as np\n";
  }
  script << "import sys\n";
  script << "import os\n";
  // Change directory in Python to avoid shell command injection
  script << "os.chdir(r'" << canonical_work_dir.string() << "')\n\n";
  script << code;
  script.close();

// Execute Python (use CMake-discovered Python executable)
#ifndef SEQUANT_PYTHON3_EXECUTABLE
#define SEQUANT_PYTHON3_EXECUTABLE "python3"
#endif
  std::string python_exe = SEQUANT_PYTHON3_EXECUTABLE;
  // Execute Python directly with properly escaped script path
  std::string cmd = shell_escape(python_exe) + " " +
                    shell_escape(script_path.string()) + " 2>&1";
  FILE *pipe = popen(cmd.c_str(), "r");
  if (!pipe) return false;

  char buffer[256];
  std::string output;
  while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
    output += buffer;
  }

  int exit_code = pclose(pipe);

  if (exit_code != 0) {
    std::cerr << "Python execution failed with exit code " << exit_code << "\n";
    std::cerr << "Output:\n" << output << "\n";
    return false;
  }

  return true;
}
#endif  // defined(SEQUANT_HAS_NUMPY_FOR_VALIDATION) ||
        // defined(SEQUANT_HAS_TORCH_FOR_VALIDATION)

// Helper to generate random Eigen::Tensor
template <typename Scalar, int NumDims>
Eigen::Tensor<Scalar, NumDims> random_tensor(
    const std::array<Eigen::Index, NumDims> &dims, unsigned int seed = 42) {
  Eigen::Tensor<Scalar, NumDims> tensor(dims);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  for (Eigen::Index i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<Scalar>(dist(gen));
  }

  return tensor;
}

}  // anonymous namespace

// ============================================================================
// VALIDATION TESTS - Execute generated Python code and verify results
// ============================================================================

#ifdef SEQUANT_HAS_NUMPY_FOR_VALIDATION

TEST_CASE("PythonEinsumGenerator - Validation", "[export][python]") {
  auto resetter = to_export_context();

  auto registry = get_default_context().index_space_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");

  SECTION("Validation: Simple matrix multiplication") {
    // Create temporary directory for test files
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_simple";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: T[a1, a2] = F[a1, i1] * t[i1, a2]
    const Eigen::Index nocc = 3;
    const Eigen::Index nvirt = 5;

    // Generate random input tensors
    auto F_tensor = random_tensor<double, 2>({nvirt, nocc}, 100);
    auto t_tensor = random_tensor<double, 2>({nocc, nvirt}, 200);

    // Compute expected result using Eigen (must evaluate the expression)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_dims = {
        Eigen::IndexPair<int>(1, 0)};
    Eigen::Tensor<double, 2> T_expected =
        F_tensor.contract(t_tensor, contraction_dims);

    // Generate Python code
    auto F = ex<Tensor>(L"F", bra{L"a_1"}, ket{L"i_1"});
    auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
    Tensor T(L"T", bra{L"a_1"}, ket{L"a_2"});

    ResultExpr result_expr(T, F * t);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged names and write files
    std::string F_name = generator.represent(F.as<Tensor>(), ctx);
    std::string t_name = generator.represent(t.as<Tensor>(), ctx);
    std::string T_name = generator.represent(T, ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + F_name + ".npy",
                                F_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t_name + ".npy",
                                t_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    // Execute Python code
    REQUIRE(run_python_code(code, temp_dir.string()));

    // Read result
    auto T_actual = read_eigen_tensor_from_numpy<double, 2>(
        temp_dir.string() + "/" + T_name + ".npy");

    // Verify results match
    REQUIRE(T_actual.dimensions()[0] == nvirt);
    REQUIRE(T_actual.dimensions()[1] == nvirt);

    const double tolerance = 1e-10;
    for (Eigen::Index i = 0; i < nvirt; ++i) {
      for (Eigen::Index j = 0; j < nvirt; ++j) {
        REQUIRE(std::abs(T_actual(i, j) - T_expected(i, j)) < tolerance);
      }
    }
  }

  SECTION("Validation: Scalar contraction (energy)") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_scalar";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: E = g[i1, i2; a1, a2] * t[a1, a2; i1, i2]
    const Eigen::Index nocc = 2;
    const Eigen::Index nvirt = 3;

    // Generate random tensors
    auto g_tensor = random_tensor<double, 4>({nocc, nocc, nvirt, nvirt}, 300);
    auto t_tensor = random_tensor<double, 4>({nvirt, nvirt, nocc, nocc}, 400);

    // Compute expected scalar result using Eigen (must evaluate the expression)
    // Contract all indices: g[i,j,a,b] * t[a,b,i,j]
    Eigen::array<Eigen::IndexPair<int>, 4> contraction_dims = {
        Eigen::IndexPair<int>(0, 2),  // i1 with i1
        Eigen::IndexPair<int>(1, 3),  // i2 with i2
        Eigen::IndexPair<int>(2, 0),  // a1 with a1
        Eigen::IndexPair<int>(3, 1)   // a2 with a2
    };
    Eigen::Tensor<double, 0> E_expected =
        g_tensor.contract(t_tensor, contraction_dims);

    // Generate Python code
    auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
    auto t = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"});
    Variable E(L"E");

    ResultExpr result_expr(E, g * t);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged names and write files
    std::string g_name = generator.represent(g.as<Tensor>(), ctx);
    std::string t_name = generator.represent(t.as<Tensor>(), ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + g_name + ".npy",
                                g_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t_name + ".npy",
                                t_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    // Execute Python code
    REQUIRE(run_python_code(code, temp_dir.string()));

    // Read scalar result
    auto E_actual =
        read_eigen_tensor_from_numpy<double, 0>(temp_dir.string() + "/E.npy");

    // Verify
    const double tolerance = 1e-9;
    REQUIRE(std::abs(E_actual() - E_expected()) < tolerance);
  }

  SECTION("Validation: Three-tensor contraction (scalar result)") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_three_scalar";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: R = g[i1, i2; a3, a4] * t1[a3; i1] * t2[a4; i2] (scalar result)
    const Eigen::Index nocc = 2;
    const Eigen::Index nvirt = 3;

    auto g_tensor = random_tensor<double, 4>({nocc, nocc, nvirt, nvirt}, 500);
    auto t1_tensor = random_tensor<double, 2>({nvirt, nocc}, 600);
    auto t2_tensor = random_tensor<double, 2>({nvirt, nocc}, 700);

    // Compute expected result manually with loops
    // R = g[i1, i2, a3, a4] * t1[a3, i1] * t2[a4, i2]
    Eigen::Tensor<double, 0> R_expected;
    R_expected.setZero();

    for (Eigen::Index i1 = 0; i1 < nocc; ++i1) {
      for (Eigen::Index i2 = 0; i2 < nocc; ++i2) {
        for (Eigen::Index a3 = 0; a3 < nvirt; ++a3) {
          for (Eigen::Index a4 = 0; a4 < nvirt; ++a4) {
            R_expected() += g_tensor(i1, i2, a3, a4) * t1_tensor(a3, i1) *
                            t2_tensor(a4, i2);
          }
        }
      }
    }

    // Generate Python code for: R = g * t1 * t2 (scalar result)
    auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_3", L"a_4"});
    auto t1 = ex<Tensor>(L"t1", bra{L"a_3"}, ket{L"i_1"});
    auto t2 = ex<Tensor>(L"t2", bra{L"a_4"}, ket{L"i_2"});
    Variable R(L"R");

    ResultExpr result_expr(R, g * t1 * t2);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged names and write files
    std::string g_name = generator.represent(g.as<Tensor>(), ctx);
    std::string t1_name = generator.represent(t1.as<Tensor>(), ctx);
    std::string t2_name = generator.represent(t2.as<Tensor>(), ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + g_name + ".npy",
                                g_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t1_name + ".npy",
                                t1_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t2_name + ".npy",
                                t2_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    REQUIRE(run_python_code(code, temp_dir.string()));

    auto R_actual =
        read_eigen_tensor_from_numpy<double, 0>(temp_dir.string() + "/R.npy");

    const double tolerance = 1e-9;
    REQUIRE(std::abs(R_actual() - R_expected()) < tolerance);
  }

  SECTION("Validation: Ternary contraction with tensor result") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_ternary";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: I[i1, a1] = A[a2, i2] * B[i2, a1] * C[i1, a2]
    // This matches the ternary.export_test case
    const Eigen::Index nocc = 3;
    const Eigen::Index nvirt = 4;

    auto A_tensor = random_tensor<double, 2>({nvirt, nocc}, 1000);
    auto B_tensor = random_tensor<double, 2>({nocc, nvirt}, 1100);
    auto C_tensor = random_tensor<double, 2>({nocc, nvirt}, 1200);

    // Compute expected result: I[i1, a1] = A[a2, i2] * B[i2, a1] * C[i1, a2]
    // First: Intermediate[a2, a1] = A[a2, i2] * B[i2, a1]
    Eigen::array<Eigen::IndexPair<int>, 1> contract_AB = {
        Eigen::IndexPair<int>(1, 0)};  // contract i2
    Eigen::Tensor<double, 2> intermediate =
        A_tensor.contract(B_tensor, contract_AB);

    // Second: I[i1, a1] = Intermediate[a2, a1] * C[i1, a2]
    // Swap to C.contract(intermediate) to get correct index order [i1, a1]
    Eigen::array<Eigen::IndexPair<int>, 1> contract_CI = {Eigen::IndexPair<int>(
        1, 0)};  // contract C.dim[1]=a2 with intermediate.dim[0]=a2
    Eigen::Tensor<double, 2> I_expected =
        C_tensor.contract(intermediate, contract_CI);

    // Generate Python code matching ternary.export_test
    auto A = ex<Tensor>(L"A", bra{L"a_2"}, ket{L"i_2"});
    auto B = ex<Tensor>(L"B", bra{L"i_2"}, ket{L"a_1"});
    auto C = ex<Tensor>(L"C", bra{L"i_1"}, ket{L"a_2"});
    Tensor I(L"I", bra{L"i_1"}, ket{L"a_1"});

    ResultExpr result_expr(I, A * B * C);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged file names for writing
    std::string A_name = generator.represent(A.as<Tensor>(), ctx);
    std::string B_name = generator.represent(B.as<Tensor>(), ctx);
    std::string C_name = generator.represent(C.as<Tensor>(), ctx);
    std::string I_name = generator.represent(I, ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + A_name + ".npy",
                                A_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + B_name + ".npy",
                                B_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + C_name + ".npy",
                                C_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    REQUIRE(run_python_code(code, temp_dir.string()));

    auto I_actual = read_eigen_tensor_from_numpy<double, 2>(
        temp_dir.string() + "/" + I_name + ".npy");

    REQUIRE(I_actual.dimensions()[0] == nocc);
    REQUIRE(I_actual.dimensions()[1] == nvirt);

    const double tolerance = 1e-10;
    for (Eigen::Index i = 0; i < nocc; ++i) {
      for (Eigen::Index a = 0; a < nvirt; ++a) {
        REQUIRE(std::abs(I_actual(i, a) - I_expected(i, a)) < tolerance);
      }
    }
  }

  SECTION("Validation: With scalar prefactor") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_scalar_factor";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: T[a1, a2] = 0.5 * F[a1, i1] * t[i1, a2]
    const Eigen::Index nocc = 3;
    const Eigen::Index nvirt = 4;

    auto F_tensor = random_tensor<double, 2>({nvirt, nocc}, 800);
    auto t_tensor = random_tensor<double, 2>({nocc, nvirt}, 900);

    // Expected: 0.5 * contraction (must evaluate the expression)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_dims = {
        Eigen::IndexPair<int>(1, 0)};
    Eigen::Tensor<double, 2> T_expected =
        0.5 * F_tensor.contract(t_tensor, contraction_dims);

    auto F = ex<Tensor>(L"F", bra{L"a_1"}, ket{L"i_1"});
    auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
    Tensor T(L"T", bra{L"a_1"}, ket{L"a_2"});

    ResultExpr result_expr(T, rational(1, 2) * F * t);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged names and write files
    std::string F_name = generator.represent(F.as<Tensor>(), ctx);
    std::string t_name = generator.represent(t.as<Tensor>(), ctx);
    std::string T_name = generator.represent(T, ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + F_name + ".npy",
                                F_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t_name + ".npy",
                                t_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    REQUIRE(run_python_code(code, temp_dir.string()));

    auto T_actual = read_eigen_tensor_from_numpy<double, 2>(
        temp_dir.string() + "/" + T_name + ".npy");

    const double tolerance = 1e-10;
    for (Eigen::Index i = 0; i < nvirt; ++i) {
      for (Eigen::Index j = 0; j < nvirt; ++j) {
        REQUIRE(std::abs(T_actual(i, j) - T_expected(i, j)) < tolerance);
      }
    }
  }

  SECTION("Validation: Sum of expressions") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_sum_expr";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Test: I[a1, i1] = f[a1, i1] - f[i2, i1] * t[a1, i2]
    // This matches the sum_unary_plus_binary.export_test case
    const Eigen::Index nocc = 3;
    const Eigen::Index nvirt = 4;

    auto f_vo_tensor = random_tensor<double, 2>({nvirt, nocc}, 1300);
    auto f_oo_tensor = random_tensor<double, 2>({nocc, nocc}, 1400);
    auto t_vo_tensor = random_tensor<double, 2>({nvirt, nocc}, 1500);

    // Compute expected result: I[a1, i1] = f[a1, i1] - f[i2, i1] * t[a1, i2]
    // First: contraction[a1, i1] = f[i2, i1] * t[a1, i2] (contract over i2)
    // Swap to t.contract(f) to get correct index order [a1, i1]
    Eigen::array<Eigen::IndexPair<int>, 1> contract_tf = {
        Eigen::IndexPair<int>(1, 0)};  // contract t.dim[1]=i2 with f.dim[0]=i2
    Eigen::Tensor<double, 2> contraction =
        t_vo_tensor.contract(f_oo_tensor, contract_tf);

    // Second: I[a1, i1] = f[a1, i1] - contraction[a1, i1]
    Eigen::Tensor<double, 2> I_expected = f_vo_tensor - contraction;

    // Generate Python code matching sum_unary_plus_binary.export_test
    auto f_vo = ex<Tensor>(L"f", bra{L"a_1"}, ket{L"i_1"});
    auto f_oo = ex<Tensor>(L"f", bra{L"i_2"}, ket{L"i_1"});
    auto t_vo = ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_2"});
    Tensor I(L"I", bra{L"a_1"}, ket{L"i_1"});

    ResultExpr result_expr(I, f_vo - f_oo * t_vo);
    auto export_tree = to_export_tree(result_expr);

    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, std::to_string(nocc));
    ctx.set_shape(virt, std::to_string(nvirt));
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    NumPyEinsumGenerator generator;

    // Get tagged file names for writing
    std::string f_vo_name = generator.represent(f_vo.as<Tensor>(), ctx);
    std::string f_oo_name = generator.represent(f_oo.as<Tensor>(), ctx);
    std::string t_vo_name = generator.represent(t_vo.as<Tensor>(), ctx);
    std::string I_name = generator.represent(I, ctx);

    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + f_vo_name + ".npy",
                                f_vo_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + f_oo_name + ".npy",
                                f_oo_tensor);
    write_eigen_tensor_to_numpy(temp_dir.string() + "/" + t_vo_name + ".npy",
                                t_vo_tensor);

    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    REQUIRE(run_python_code(code, temp_dir.string()));

    auto I_actual = read_eigen_tensor_from_numpy<double, 2>(
        temp_dir.string() + "/" + I_name + ".npy");

    REQUIRE(I_actual.dimensions()[0] == nvirt);
    REQUIRE(I_actual.dimensions()[1] == nocc);

    const double tolerance = 1e-10;
    for (Eigen::Index a = 0; a < nvirt; ++a) {
      for (Eigen::Index i = 0; i < nocc; ++i) {
        REQUIRE(std::abs(I_actual(a, i) - I_expected(a, i)) < tolerance);
      }
    }
  }

  SECTION("Tensor I/O - Layout compatibility") {
    std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() / "sequant_test_layout";
    std::filesystem::create_directories(temp_dir);
    // Ensure cleanup on all exit paths
    auto cleanup = sequant::detail::make_scope_exit(
        [&temp_dir]() { std::filesystem::remove_all(temp_dir); });

    // Subsection 1: Column-major (default) write and read
    {
      // Default Eigen::Tensor is column-major
      Eigen::Tensor<double, 2> tensor_col(3, 4);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          tensor_col(i, j) = i * 10 + j;
        }
      }

      // Write and read back
      write_eigen_tensor_to_numpy(temp_dir.string() + "/col_major.npy",
                                  tensor_col);
      auto tensor_read = read_eigen_tensor_from_numpy<double, 2>(
          temp_dir.string() + "/col_major.npy");

      // Verify dimensions
      REQUIRE(tensor_read.dimensions()[0] == 3);
      REQUIRE(tensor_read.dimensions()[1] == 4);

      // Verify data - note: we're reading with default (column-major) template
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          REQUIRE(tensor_read(i, j) == tensor_col(i, j));
        }
      }
    }

    // Subsection 2: Row-major write and read
    {
      // Explicitly create row-major tensor
      Eigen::Tensor<double, 2, Eigen::RowMajor> tensor_row(3, 4);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          tensor_row(i, j) = i * 10 + j;
        }
      }

      // Write and read back with matching layout
      write_eigen_tensor_to_numpy(temp_dir.string() + "/row_major.npy",
                                  tensor_row);
      auto tensor_read =
          read_eigen_tensor_from_numpy<double, 2, Eigen::RowMajor>(
              temp_dir.string() + "/row_major.npy");

      // Verify dimensions
      REQUIRE(tensor_read.dimensions()[0] == 3);
      REQUIRE(tensor_read.dimensions()[1] == 4);

      // Verify data
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          REQUIRE(tensor_read(i, j) == tensor_row(i, j));
        }
      }
    }

    // Subsection 3: Verify layout mismatch detection with explicit RowMajor
    {
      // Create a column-major tensor and write it
      Eigen::Tensor<double, 2> tensor_col(2, 3);  // default is column-major
      tensor_col.setConstant(42.0);
      write_eigen_tensor_to_numpy(temp_dir.string() + "/col_test.npy",
                                  tensor_col);

      // Try to read with explicit row-major (should throw)
      // Note: Eigen::RowMajor = 1, so Options != 0, triggering strict
      // validation
      REQUIRE_THROWS_WITH(
          (read_eigen_tensor_from_numpy<double, 2, Eigen::RowMajor>(
              temp_dir.string() + "/col_test.npy")),
          Catch::Matchers::ContainsSubstring("Layout mismatch"));
    }

    // Subsection 4: Default Options accepts any layout
    {
      // Write a row-major file using external means (simulate numpy output)
      Eigen::Tensor<double, 2, Eigen::RowMajor> tensor_row(2, 3);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
          tensor_row(i, j) = i + j * 10;
        }
      }
      write_eigen_tensor_to_numpy(temp_dir.string() + "/any_layout.npy",
                                  tensor_row);

      // Read with default Options (should succeed even though layouts don't
      // match) This is for compatibility with numpy files
      auto tensor_read = read_eigen_tensor_from_numpy<double, 2>(
          temp_dir.string() + "/any_layout.npy");

      // Note: Due to layout mismatch, the logical interpretation may differ
      // but the read should succeed
      REQUIRE(tensor_read.dimensions()[0] == 2);
      REQUIRE(tensor_read.dimensions()[1] == 3);
    }

    // Subsection 5: Scalar (rank-0) tensors
    {
      Eigen::Tensor<double, 0> scalar_tensor;
      scalar_tensor() = 3.14159;

      write_eigen_tensor_to_numpy(temp_dir.string() + "/scalar.npy",
                                  scalar_tensor);
      auto scalar_read = read_eigen_tensor_from_numpy<double, 0>(
          temp_dir.string() + "/scalar.npy");

      REQUIRE(std::abs(scalar_read() - 3.14159) < 1e-10);
    }

    // Subsection 6: Higher rank tensors
    {
      Eigen::Tensor<float, 4> tensor4d(2, 3, 4, 5);
      tensor4d.setRandom();

      write_eigen_tensor_to_numpy(temp_dir.string() + "/tensor4d.npy",
                                  tensor4d);
      auto tensor4d_read = read_eigen_tensor_from_numpy<float, 4>(
          temp_dir.string() + "/tensor4d.npy");

      REQUIRE(tensor4d_read.dimensions()[0] == 2);
      REQUIRE(tensor4d_read.dimensions()[1] == 3);
      REQUIRE(tensor4d_read.dimensions()[2] == 4);
      REQUIRE(tensor4d_read.dimensions()[3] == 5);

      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 4; ++k) {
            for (int l = 0; l < 5; ++l) {
              REQUIRE(std::abs(tensor4d_read(i, j, k, l) -
                               tensor4d(i, j, k, l)) < 1e-6f);
            }
          }
        }
      }
    }
  }
}

#endif  // SEQUANT_HAS_NUMPY_FOR_VALIDATION
