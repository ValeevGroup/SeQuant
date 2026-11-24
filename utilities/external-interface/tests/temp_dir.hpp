#ifndef SEQUANT_EXTERNAL_INTERFACE_TESTS_TEMP_DIR_HPP
#define SEQUANT_EXTERNAL_INTERFACE_TESTS_TEMP_DIR_HPP

#include <filesystem>
#include <format>
#include <fstream>
#include <random>
#include <string>
#include <system_error>

namespace sequant::util::extint::test_support {

/// A fresh directory under the system's temporary directory, removed with
/// everything in it on destruction; its random name keeps concurrent test runs
/// apart
struct TempDir {
  std::filesystem::path path;

  TempDir() {
    std::random_device rd;
    path = std::filesystem::temp_directory_path() /
           std::format("sequant_extint_{:08x}{:08x}", rd(), rd());
    std::filesystem::create_directories(path);
  }
  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
  }
  TempDir(const TempDir &) = delete;
  TempDir &operator=(const TempDir &) = delete;

  /// Writes @p content to the file @p name (relative to the directory),
  /// creating its parent directories
  void write(const std::string &name, const std::string &content) const {
    std::filesystem::create_directories((path / name).parent_path());
    std::ofstream(path / name) << content;
  }
};

}  // namespace sequant::util::extint::test_support

#endif  // SEQUANT_EXTERNAL_INTERFACE_TESTS_TEMP_DIR_HPP
