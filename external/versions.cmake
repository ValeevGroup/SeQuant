# for each dependency track both current and previous id (the variable for the latter must contain PREVIOUS)
# to be able to auto-update them

set(SEQUANT_TRACKED_VGCMAKEKIT_TAG 6ecd3689f3d33d4426b47f8b68ba81b0efb7c80b)

set(SEQUANT_TRACKED_RANGEV3_TAG 0.12.0)
set(SEQUANT_TRACKED_RANGEV3_PREVIOUS_TAG d800a032132512a54c291ce55a2a43e0460591c7)

set(SEQUANT_TRACKED_TILEDARRAY_TAG 5477174cd8e857f76d29cf6c0dcc80277393a9eb)
set(SEQUANT_TRACKED_TILEDARRAY_PREVIOUS_TAG 131ed04c989a5c6022d13e0dfe848d7112ffa9e4)

set(SEQUANT_TRACKED_LIBPERM_TAG cada3e185549896203cf4d0c7f26ea22c7de428f)
set(SEQUANT_TRACKED_LIBPERM_PREVIOUS_TAG 8e4afd1461baffa5d829c8fed059f5a172a3b060)

set(SEQUANT_TRACKED_POLYMORPHICVARIANT_TAG 010c69786104c07c5faccffe0e99f99de5a69fd8)
set(SEQUANT_TRACKED_POLYMORPHICVARIANT_PREVIOUS_TAG 010c69786104c07c5faccffe0e99f99de5a69fd8)

# oldest Boost we can tolerate ... with some compilers should be able to use an earlier version, but:
# - Boost.ContainerHash <1.81 uses unary_function that has been deprecated in C++11 and removed in C++17:
#   https://github.com/boostorg/container_hash/blob/boost-1.80.0/include/boost/container_hash/hash.hpp#L132
#   Recent stdlibs have removed this class when using C++17 (e.g. Apple Clang 15)
set(SEQUANT_OLDEST_BOOST_VERSION 1.81)

set(SEQUANT_TRACKED_CATCH2_TAG v3.7.1)

set(SEQUANT_TRACKED_GOOGLEBENCHMARK_TAG v1.9.1)
