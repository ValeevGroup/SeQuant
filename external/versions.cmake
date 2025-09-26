# for each dependency track both current and previous id (the variable for the latter must contain PREVIOUS)
# to be able to auto-update them

set(SEQUANT_TRACKED_VGCMAKEKIT_TAG 748f2499684c453109098aa06ccd03e10e915336)

set(SEQUANT_TRACKED_RANGEV3_TAG 0.12.0)

set(SEQUANT_TRACKED_TILEDARRAY_TAG 5477174cd8e857f76d29cf6c0dcc80277393a9eb)

set(SEQUANT_TRACKED_LIBPERM_TAG cada3e185549896203cf4d0c7f26ea22c7de428f)

set(SEQUANT_TRACKED_POLYMORPHICVARIANT_TAG 010c69786104c07c5faccffe0e99f99de5a69fd8)

set(SEQUANT_TRACKED_UTFCPP_TAG v4.0.6)

set(SEQUANT_TRACKED_CLI11_TAG v2.5.0)
set(SEQUANT_OLDEST_CLI11_VERSION 2)

set(SEQUANT_TRACKED_SPDLOG_TAG v1.15.3)

set(SEQUANT_TRACKED_JSON_TAG v3.12.0)
set(SEQUANT_OLDEST_JSON_VERSION 3)

# oldest Boost we can tolerate ... with some compilers should be able to use an earlier version, but:
# - Boost.ContainerHash <1.81 uses unary_function that has been deprecated in C++11 and removed in C++17:
#   https://github.com/boostorg/container_hash/blob/boost-1.80.0/include/boost/container_hash/hash.hpp#L132
#   Recent stdlibs have removed this class when using C++17 (e.g. Apple Clang 15)
set(SEQUANT_OLDEST_BOOST_VERSION 1.81)

set(SEQUANT_TRACKED_CATCH2_TAG v3.9.1)
set(SEQUANT_OLDEST_CATCH2_VERSION 3.3)

set(SEQUANT_TRACKED_GOOGLEBENCHMARK_TAG v1.9.4)
set(SEQUANT_OLDEST_GOOGLEBENCHMARK_VERSION 1.9.3)

set(SEQUANT_TRACKED_PYBIND11_TAG v3.0.1)
set(SEQUANT_OLDEST_PYBIND11_VERSION 3)
