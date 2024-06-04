# for each dependency track both current and previous id (the variable for the latter must contain PREVIOUS)
# to be able to auto-update them

set(SEQUANT_TRACKED_VGCMAKEKIT_TAG 4f40440dcbda2a5a005fd15f27c582a5587ee779)

set(SEQUANT_TRACKED_RANGEV3_TAG 0.12.0)
set(SEQUANT_TRACKED_RANGEV3_PREVIOUS_TAG d800a032132512a54c291ce55a2a43e0460591c7)

set(SEQUANT_TRACKED_TILEDARRAY_TAG 179e1d39f0b8635f37d72be6d3c91f315a99eaea)
set(SEQUANT_TRACKED_TILEDARRAY_PREVIOUS_TAG 6a134f31f38a1e643583d5c8f19ed85fb61b04ae)

# oldest Boost we can tolerate ... with some compilers should be able to use an earlier version, but:
# - Boost.ContainerHash <1.81 uses unary_function that has been deprecated in C++11 and removed in C++17:
#   https://github.com/boostorg/container_hash/blob/boost-1.80.0/include/boost/container_hash/hash.hpp#L132
#   Recent stdlibs have removed this class when using C++17 (e.g. Apple Clang 15)
set(SEQUANT_OLDEST_BOOST_VERSION 1.81)

set(SEQUANT_TRACKED_CATCH2_TAG v3.5.4)
