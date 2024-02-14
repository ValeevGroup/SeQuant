# for each dependency track both current and previous id (the variable for the latter must contain PREVIOUS)
# to be able to auto-update them

set(SEQUANT_TRACKED_VGCMAKEKIT_TAG 1d1a13673c9acf6b5a58f19504dd1ad399a20f54)

set(SEQUANT_TRACKED_RANGEV3_TAG 0.12.0)
set(SEQUANT_TRACKED_RANGEV3_PREVIOUS_TAG d800a032132512a54c291ce55a2a43e0460591c7)

set(SEQUANT_TRACKED_TILEDARRAY_TAG 3f6629db047417e814b75ad5069b7f4ce26428e7)
set(SEQUANT_TRACKED_TILEDARRAY_PREVIOUS_TAG 954d861f553e938c3cfc4892fce9234bf4bf7193)

# oldest Boost we can tolerate ... with some compilers should be able to use an earlier version, but:
# - Boost.ContainerHash <1.81 uses unary_function that has been deprecated in C++11 and removed in C++17:
#   https://github.com/boostorg/container_hash/blob/boost-1.80.0/include/boost/container_hash/hash.hpp#L132
#   Recent stdlibs have removed this class when using C++17 (e.g. Apple Clang 15)
set(SEQUANT_OLDEST_BOOST_VERSION 1.81)
