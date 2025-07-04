cmake_minimum_required(VERSION 3.27) # 3.27 required for CONFIG_FILE option in doxygen_add_docs

# Preload versions/tags of all dependencies ====================================
include(external/versions.cmake)

###############################################################################
# CMake defaults to address key pain points
###############################################################################

# safety net for dev workflow: accidental install will not affect FindOrFetch*
if (NOT DEFINED CACHE{CMAKE_FIND_NO_INSTALL_PREFIX})
    set(CMAKE_FIND_NO_INSTALL_PREFIX ON CACHE BOOL "Whether find_* commands will search CMAKE_INSTALL_PREFIX and CMAKE_STAGING_PREFIX; see https://cmake.org/cmake/help/latest/variable/CMAKE_FIND_NO_INSTALL_PREFIX.html#variable:CMAKE_FIND_NO_INSTALL_PREFIX")
endif()

###############################################################################
# Bring ValeevGroup cmake toolkit
###############################################################################
include(FetchContent)
if (DEFINED PROJECT_BINARY_DIR)
    set(VG_CMAKE_KIT_PREFIX_DIR PROJECT_BINARY_DIR)
else ()
    set(VG_CMAKE_KIT_PREFIX_DIR CMAKE_CURRENT_BINARY_DIR)
endif ()
FetchContent_Declare(
        vg_cmake_kit
        QUIET
        GIT_REPOSITORY https://github.com/ValeevGroup/kit-cmake.git
        GIT_TAG ${SEQUANT_TRACKED_VGCMAKEKIT_TAG}
        SOURCE_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg
        BINARY_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg-build
        SUBBUILD_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg-subbuild
)
FetchContent_MakeAvailable(vg_cmake_kit)
list(APPEND CMAKE_MODULE_PATH "${vg_cmake_kit_SOURCE_DIR}/modules")

# Set SeQuant version
set(SEQUANT_MAJOR_VERSION 2)
set(SEQUANT_MINOR_VERSION 0)
set(SEQUANT_MICRO_VERSION 0)
set(SEQUANT_PRERELEASE_ID alpha.1)
set(SEQUANT_VERSION "${SEQUANT_MAJOR_VERSION}.${SEQUANT_MINOR_VERSION}.${SEQUANT_MICRO_VERSION}")
if (SEQUANT_PRERELEASE_ID)
    set(SEQUANT_EXT_VERSION "${SEQUANT_VERSION}-${SEQUANT_PRERELEASE_ID}")
else (SEQUANT_PRERELEASE_ID)
    set(SEQUANT_EXT_VERSION "${SEQUANT_VERSION}")
endif (SEQUANT_PRERELEASE_ID)

# extract git metadata
include(GetGitMetadata)
vgkit_cmake_git_metadata()

# make SeQuant project
project(SeQuant LANGUAGES CXX VERSION "${SEQUANT_VERSION}" DESCRIPTION "SEcond QUANTization toolkit")

# need C++20, insist on strict standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ ISO Standard version")
if (CMAKE_CXX_STANDARD LESS 20)
    message(FATAL_ERROR "C++ 2020 ISO Standard or higher is required to compile SeQuant")
endif ()
# C++20 is only configurable via compile features with cmake 3.12 and older
if (CMAKE_CXX_STANDARD EQUAL 20 AND CMAKE_VERSION VERSION_LESS 3.12.0)
    cmake_minimum_required(VERSION 3.12.0)
endif ()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Whether to use extensions of C++ ISO Standard version")

# Enable ccache if not already enabled by symlink masquerading
if (NOT CMAKE_CXX_COMPILER MATCHES ".*/ccache$")
    find_program(CCACHE_EXECUTABLE ccache)
    if (CCACHE_EXECUTABLE)
        mark_as_advanced(CCACHE_EXECUTABLE)
        message(STATUS "Found ccache: ${CCACHE_EXECUTABLE}")
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXECUTABLE}" CACHE STRING "Compiler launcher to use for compiling C++")
        set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_EXECUTABLE}" CACHE STRING "Compiler launcher to use for compiling C")
    else ()
        set(CCACHE_EXECUTABLE "")
    endif ()
endif ()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")

set(FETCHCONTENT_UPDATES_DISCONNECTED OFF CACHE BOOL "Enables UPDATE_DISCONNECTED behavior for all content population")
include(CTest)
include(FetchContent)
include(AddCustomTargetSubproject)
include(FeatureSummary)
include(CMakePackageConfigHelpers)

##########################
# Standard build variables
##########################
set(SEQUANT_INSTALL_BINDIR "bin"
        CACHE PATH "SeQuant BIN install directory")
set(SEQUANT_INSTALL_INCLUDEDIR "include"
        CACHE PATH "SeQuant INCLUDE install directory")
set(SEQUANT_INSTALL_LIBDIR "lib"
        CACHE PATH "SeQuant LIB install directory")
set(SEQUANT_INSTALL_SHAREDIR "share/sequant/${SEQUANT_MAJOR_VERSION}.${SEQUANT_MINOR_VERSION}.${SEQUANT_MICRO_VERSION}"
        CACHE PATH "SeQuant SHARE install directory")
set(SEQUANT_INSTALL_DATADIR "${SEQUANT_INSTALL_SHAREDIR}/data"
        CACHE PATH "SeQuant DATA install directory")
set(SEQUANT_INSTALL_DOCDIR "${SEQUANT_INSTALL_SHAREDIR}/doc"
        CACHE PATH "SeQuant DOC install directory")
set(SEQUANT_INSTALL_CMAKEDIR "lib/cmake/sequant"
        CACHE PATH "SeQuant CMAKE install directory")
option(SEQUANT_BENCHMARKS "Enable SeQuant benchmarks" OFF)

############################
# Additional build variables
############################
option(SEQUANT_IWYU "Whether to use the include-what-you-use tool (if found)" OFF)

##########################
# SeQuant package options
##########################

option(ENABLE_TBB "Enable TBB as an optional prerequisite for C++'s parallel STL (if PSTL is not supported will use of standard C++11 concurrency)" OFF)
add_feature_info(TBB ENABLE_TBB "Intel Thread-Building Blocks is an optional prerequisite for C++'s parallel STL")

option(SEQUANT_PYTHON "Build SeQuant python module" OFF)
add_feature_info(PYTHON SEQUANT_PYTHON "PySeQuant: Python bindings to SeQuant")

option(SEQUANT_USE_SYSTEM_BOOST_HASH "Use system Boost for hashing? Set to OFF to make hashing independent of Boost, thus value-portable" ON)
add_feature_info(SEQUANT_USE_SYSTEM_BOOST_HASH SEQUANT_USE_SYSTEM_BOOST_HASH "SeQuant uses system Boost for hashing (thus results depend on the Boost version)")

##########################
# Prerequisites
##########################

# std::thread requires (on some platforms?) thread support which is not
# provided automatically
find_package(Threads REQUIRED)

# PSTL (used by g++-9 and clang++-8 in c++17+ mode) needs TBB
if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
        ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 9) OR
(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND
        ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 8))
    if (ENABLE_TBB)  # but user must enable the search for TBB since this is an additional source of build entropy
        find_package(TBB REQUIRED)
        # TBB::tbb by default is not GLOBAL, so to allow users of LINALG_LIBRARIES to safely use it we need to make it global
        # more discussion here: https://gitlab.kitware.com/cmake/cmake/-/issues/17256
        set_target_properties(TBB::tbb PROPERTIES IMPORTED_GLOBAL TRUE)
    endif (ENABLE_TBB)
endif ()
# check of <execution> header is usable
include(CheckCXXFeatures)
check_cxx_execution_header(SEQUANT)

# Ranges-V3
include(FindOrFetchRangeV3)

# Boost will be added after defining SeQuant
include(external/boost.cmake)

# LibPerm
include(FindOrFetchLibPerm)

# embedded bliss-0.73
add_library(SeQuant-bliss
        SeQuant/external/bliss/defs.cc
        SeQuant/external/bliss/defs.hh
        SeQuant/external/bliss/graph.cc
        SeQuant/external/bliss/graph.hh
        SeQuant/external/bliss/partition.cc
        SeQuant/external/bliss/partition.hh
        SeQuant/external/bliss/orbit.cc
        SeQuant/external/bliss/orbit.hh
        SeQuant/external/bliss/uintseqhash.cc
        SeQuant/external/bliss/uintseqhash.hh
        SeQuant/external/bliss/heap.cc
        SeQuant/external/bliss/heap.hh
        SeQuant/external/bliss/timer.cc
        SeQuant/external/bliss/timer.hh
        SeQuant/external/bliss/utils.cc
        SeQuant/external/bliss/utils.hh
        SeQuant/external/bliss/bliss_C.cc
        SeQuant/external/bliss/bliss_C.h
        )
target_link_libraries(SeQuant-bliss PUBLIC range-v3::range-v3)

set(SeQuant_src
        ${PROJECT_BINARY_DIR}/SeQuant/version.hpp
        SeQuant/version.cpp
        SeQuant/core/abstract_tensor.hpp
        SeQuant/core/algorithm.hpp
        SeQuant/core/any.hpp
        SeQuant/core/asy_cost.cpp
        SeQuant/core/asy_cost.hpp
        SeQuant/core/attr.hpp
        SeQuant/core/binary_node.hpp
        SeQuant/core/biorthogonalization.cpp
        SeQuant/core/biorthogonalization.hpp
        SeQuant/core/bliss.hpp
        SeQuant/core/complex.hpp
        SeQuant/core/container.hpp
        SeQuant/core/context.cpp
        SeQuant/core/context.hpp
        SeQuant/core/eval_expr.cpp
        SeQuant/core/eval_expr.hpp
        SeQuant/core/eval_node.hpp
        SeQuant/core/export/itf.cpp
        SeQuant/core/export/itf.hpp
        SeQuant/core/expr.cpp
        SeQuant/core/expr.hpp
        SeQuant/core/expr_algorithm.hpp
        SeQuant/core/expr_operator.hpp
        SeQuant/core/hash.cpp
        SeQuant/core/hash.hpp
        SeQuant/core/hugenholtz.hpp
        SeQuant/core/index.cpp
        SeQuant/core/index.hpp
        SeQuant/core/index_space_registry.cpp
        SeQuant/core/index_space_registry.hpp
        SeQuant/core/interval.hpp
        SeQuant/core/latex.cpp
        SeQuant/core/latex.ipp
        SeQuant/core/latex.hpp
        SeQuant/core/logger.hpp
        SeQuant/core/math.hpp
        SeQuant/core/meta.hpp
        SeQuant/core/op.cpp
        SeQuant/core/op.hpp
        SeQuant/core/optimize.hpp
        SeQuant/core/optimize/fusion.cpp
        SeQuant/core/optimize/fusion.hpp
        SeQuant/core/optimize/optimize.cpp
        SeQuant/core/parse/deparse.cpp
        SeQuant/core/parse/parse.cpp
        SeQuant/core/parse.hpp
        SeQuant/core/ranges.hpp
        SeQuant/core/rational.hpp
        SeQuant/core/result_expr.cpp
        SeQuant/core/result_expr.hpp
        SeQuant/core/runtime.cpp
        SeQuant/core/runtime.hpp
        SeQuant/core/space.cpp
        SeQuant/core/space.hpp
        SeQuant/core/tag.hpp
        SeQuant/core/tensor.cpp
        SeQuant/core/tensor.hpp
        SeQuant/core/tensor_canonicalizer.cpp
        SeQuant/core/tensor_canonicalizer.hpp
        SeQuant/core/tensor_network.cpp
        SeQuant/core/tensor_network.hpp
        SeQuant/core/tensor_network_v2.cpp
        SeQuant/core/tensor_network_v2.hpp
        SeQuant/core/tensor_network/canonicals.hpp
        SeQuant/core/tensor_network/slot.hpp
        SeQuant/core/tensor_network/vertex.hpp
        SeQuant/core/tensor_network/vertex_painter.cpp
        SeQuant/core/tensor_network/vertex_painter.hpp
        SeQuant/core/timer.hpp
        SeQuant/core/utility/context.hpp
        SeQuant/core/utility/expr.cpp
        SeQuant/core/utility/expr.hpp
        SeQuant/core/utility/indices.hpp
        SeQuant/core/utility/macros.hpp
        SeQuant/core/utility/nodiscard.hpp
        SeQuant/core/utility/permutation.hpp
        SeQuant/core/utility/singleton.hpp
        SeQuant/core/utility/string.hpp
        SeQuant/core/utility/string.cpp
        SeQuant/core/utility/tuple.hpp
        SeQuant/core/utility/swap.hpp
        SeQuant/core/wick.hpp
        SeQuant/core/wick.impl.hpp
        SeQuant/core/wolfram.hpp
        SeQuant/core/wstring.hpp
        SeQuant/domain/mbpt/antisymmetrizer.cpp
        SeQuant/domain/mbpt/antisymmetrizer.hpp
        SeQuant/domain/mbpt/context.hpp
        SeQuant/domain/mbpt/context.cpp
        SeQuant/domain/mbpt/convention.cpp
        SeQuant/domain/mbpt/convention.hpp
        SeQuant/domain/mbpt/op.cpp
        SeQuant/domain/mbpt/op.hpp
        SeQuant/domain/mbpt/op.ipp
        SeQuant/domain/mbpt/rdm.cpp
        SeQuant/domain/mbpt/rdm.hpp
        SeQuant/domain/mbpt/space_qns.hpp
        SeQuant/domain/mbpt/spin.cpp
        SeQuant/domain/mbpt/spin.hpp
        SeQuant/domain/mbpt/vac_av.ipp
        SeQuant/domain/mbpt/models/cc.cpp
        SeQuant/domain/mbpt/models/cc.hpp
        SeQuant/domain/mbpt/rules/df.cpp
        SeQuant/domain/mbpt/rules/df.hpp
        SeQuant/domain/mbpt/rules/csv.cpp
        SeQuant/domain/mbpt/rules/csv.hpp
)

### optional prereqs
if (SEQUANT_EVAL_TESTS)
    include(FindOrFetchTiledArray)
endif (SEQUANT_EVAL_TESTS)
if (NOT TARGET Eigen3::Eigen)
    # use TA's Eigen, if available
    if (TARGET TiledArray_Eigen)
        add_library(Eigen3::Eigen ALIAS TiledArray_Eigen)
    else()
      # re:NO_CMAKE_PACKAGE_REGISTRY: eigen3 registers its *build* tree with the user package registry ...
      #                               to avoid issues with wiped build directory look for installed eigen
      find_package(Eigen3 3.0 NO_MODULE NO_CMAKE_PACKAGE_REGISTRY)
    endif()
endif()
if (TARGET Eigen3::Eigen)
    set(SEQUANT_HAS_EIGEN ON)
    endif()

if (TARGET tiledarray)
    set(SEQUANT_HAS_TILEDARRAY ON)
    list(APPEND SeQuant_src
            SeQuant/domain/eval/cache_manager.cpp
            SeQuant/domain/eval/cache_manager.hpp
            SeQuant/domain/eval/eval.hpp
            SeQuant/domain/eval/result.cpp
            SeQuant/domain/eval/result.hpp
            SeQuant/domain/eval/eval_fwd.hpp
            )
endif ()

add_library(SeQuant
        ${SeQuant_src}
        )
# feed SEQUANT_GIT_REVISION and SEQUANT_GIT_DESCRIPTION to SeQuant/version.cpp only to avoid recompiling everything
set_source_files_properties(
        SeQuant/version.cpp
        PROPERTIES COMPILE_DEFINITIONS
        "SEQUANT_GIT_REVISION=\"${SEQUANT_GIT_REVISION}\";SEQUANT_GIT_DESCRIPTION=\"${SEQUANT_GIT_DESCRIPTION}\""
)

target_link_libraries(SeQuant PUBLIC range-v3::range-v3 Boost::regex Boost::locale Boost::headers SeQuant-bliss Threads::Threads)
target_link_libraries(SeQuant PRIVATE $<BUILD_INTERFACE:libperm::libperm>)
# modularized Boost has finer grained targets than just Boost::headers
if (Boost_IS_MODULARIZED)
    target_link_libraries(SeQuant PUBLIC
            Boost::container
            Boost::container_hash
            Boost::hana
            Boost::multiprecision
            Boost::numeric_conversion
            Boost::numeric_interval
            Boost::range
            Boost::spirit
    )
endif()
if (TARGET tiledarray)
    target_link_libraries(SeQuant PUBLIC tiledarray)
endif ()
if (TARGET Eigen3::Eigen)
    target_link_libraries(SeQuant PUBLIC Eigen3::Eigen)
    target_compile_definitions(SeQuant PUBLIC SEQUANT_HAS_EIGEN=1)
endif()
if (SEQUANT_HAS_EXECUTION_HEADER_STANDALONE OR SEQUANT_HAS_EXECUTION_HEADER_WITH_TBB)
    target_compile_definitions(SeQuant PUBLIC SEQUANT_HAS_EXECUTION_HEADER)
    if (SEQUANT_HAS_EXECUTION_HEADER_WITH_TBB)
        target_link_libraries(SeQuant PUBLIC ${TBB_LIBRARIES})
        target_include_directories(SeQuant PUBLIC ${TBB_INCLUDE_DIRS})
    endif ()
endif ()
target_include_directories(SeQuant PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_compile_features(SeQuant INTERFACE "cxx_std_20")
if (SEQUANT_EVAL_TRACE)
    target_compile_definitions(SeQuant PUBLIC SEQUANT_EVAL_TRACE=1)
endif ()
if (SEQUANT_USE_SYSTEM_BOOST_HASH)
    target_compile_definitions(SeQuant PUBLIC SEQUANT_USE_SYSTEM_BOOST_HASH=1)
else()
    # if not to use Boost Hash, unpack bundled Boost ContainerHash 1.81
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/SeQuant/external)
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/SeQuant/external/boost.tar.gz
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/SeQuant/external
            RESULT_VARIABLE UNPACK_BOOST_RESULT
            OUTPUT_VARIABLE UNPACK_BOOST_OUTPUT
            ERROR_VARIABLE UNPACK_BOOST_OUTPUT
    )
    if (NOT UNPACK_BOOST_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to unpack the bundled Boost! The tar command output:\n${UNPACK_BOOST_OUTPUT}")
    endif()
    install(DIRECTORY ${PROJECT_BINARY_DIR}/SeQuant
            COMPONENT sequant
            DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}"
            FILES_MATCHING PATTERN "*.hpp"
            )
endif ()
install(TARGETS SeQuant-bliss EXPORT sequant COMPONENT sequant LIBRARY DESTINATION ${SEQUANT_INSTALL_LIBDIR})
install(TARGETS SeQuant EXPORT sequant COMPONENT sequant LIBRARY DESTINATION ${SEQUANT_INSTALL_LIBDIR})
install(DIRECTORY SeQuant
        COMPONENT sequant
        DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.hh"
        )
add_library(SeQuant::SeQuant ALIAS SeQuant)  # to be able to use as subproject

if (SEQUANT_MIMALLOC)
    find_package(mimalloc REQUIRED)
    target_link_libraries(SeQuant PUBLIC mimalloc)
    target_compile_definitions(SeQuant PUBLIC SEQUANT_HAS_MIMALLOC=1)
endif ()

# build all of boost before SeQuant, including parts it does not use
if (Boost_BUILT_FROM_SOURCE AND TARGET build-boost-in-SeQuant)
    add_dependencies(SeQuant build-boost-in-SeQuant)
endif()

if (SEQUANT_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu)

    if (iwyu_path)
        set(iwyu_options_and_path
            "${iwyu_path}"
            -Xiwyu --cxx20ns
            -Xiwyu --no_comments
        )
        set_property(TARGET SeQuant PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_options_and_path})
    endif()
endif()


add_subdirectory(utilities)
add_subdirectory(tests)
add_subdirectory(doc)

if (SEQUANT_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()


####### Python ########

if (SEQUANT_PYTHON)
    if (NOT CMAKE_POSITION_INDEPENDENT_CODE)
        message(FATAL_ERROR "Python module requires CMAKE_POSITION_INDEPENDENT_CODE=ON")
    endif ()
    add_subdirectory(python)
endif ()


##########################
# export SeQuant
##########################

configure_file(
        ${PROJECT_SOURCE_DIR}/SeQuant/version.hpp.in
        ${PROJECT_BINARY_DIR}/SeQuant/version.hpp
)
install(FILES ${PROJECT_BINARY_DIR}/SeQuant/version.hpp
        DESTINATION "${SEQUANT_INSTALL_INCLUDEDIR}/SeQuant")

# Create the version file
write_basic_package_version_file(sequant-config-version.cmake
        VERSION ${SEQUANT_VERSION} COMPATIBILITY AnyNewerVersion)

# Make SeQuant consumable from the build directory without installation
export(EXPORT sequant
        NAMESPACE SeQuant::
        FILE "${PROJECT_BINARY_DIR}/sequant-targets.cmake")

## Create the configure file
configure_package_config_file(cmake/sequant-config.cmake.in
        "${PROJECT_BINARY_DIR}/sequant-config.cmake"
        INSTALL_DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
        PATH_VARS CMAKE_INSTALL_PREFIX SEQUANT_INSTALL_BINDIR
        SEQUANT_INSTALL_INCLUDEDIR SEQUANT_INSTALL_LIBDIR
        SEQUANT_INSTALL_DOCDIR SEQUANT_INSTALL_CMAKEDIR)

## Install config, version, and target files
install(EXPORT sequant
        FILE "sequant-targets.cmake"
        DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
        NAMESPACE SeQuant::
        COMPONENT sequant-config)
install(FILES
        "${PROJECT_BINARY_DIR}/sequant-config.cmake"
        "${PROJECT_BINARY_DIR}/sequant-config-version.cmake"
        DESTINATION "${SEQUANT_INSTALL_CMAKEDIR}"
        COMPONENT sequant-config)
add_custom_target_subproject(sequant install-config
        COMMAND ${CMAKE_COMMAND} -DCOMPONENT=sequant-config -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
        COMMENT "Installing SeQuant config components")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    feature_summary(WHAT ALL
                    DESCRIPTION "=== SeQuant Package/Feature Info ===")
    feature_summary(FILENAME ${CMAKE_CURRENT_BINARY_DIR}/features.log WHAT ALL)
endif()

###############################################################################
# appendix: misc details
###############################################################################
SET(CMAKE_COLOR_MAKEFILE ON)
