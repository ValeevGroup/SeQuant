cmake_minimum_required(VERSION 3.15.0) # need list(PREPEND for toolchains

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
include(CTest)

option(ENABLE_TBB "Enable TBB as an optional prerequisite for C++'s parallel STL (if PSTL is not supported will use of standard C++11 concurrency)" OFF)
add_feature_info(TBB ENABLE_TBB "Intel Thread-Building Blocks is an optional prerequisite for C++'s parallel STL")

option(SEQUANT_TESTS "Enable building of test cases" ${BUILD_TESTING})
add_feature_info(TESTS SEQUANT_TESTS "Enable building of test cases")

option(SEQUANT_EVAL_TESTS "Enable building of evaluation tests (if true, will look for and/or build TiledArray)" OFF)
add_feature_info(EVAL_TESTS SEQUANT_EVAL_TESTS "Build evaluation tests (if true, will look for and/or build TiledArray)")

option(SEQUANT_PYTHON "Build SeQuant python module" OFF)
add_feature_info(PYTHON SEQUANT_PYTHON "PySeQuant: Python bindings to SeQuant")

option(SEQUANT_EVAL_TRACE "Enable tracing of the expression interpretation (configure TiledArray with TA_TENSOR_MEM_PROFILE=ON to trace memory use)" OFF)
add_feature_info(EVAL_TRACE SEQUANT_EVAL_TRACE "Trace expression interpretation")

option(SEQUANT_USE_SYSTEM_BOOST_HASH "Use system Boost for hashing? Set to OFF to make hashing independent of Boost, thus value-portable" ON)
add_feature_info(SEQUANT_USE_SYSTEM_BOOST_HASH SEQUANT_USE_SYSTEM_BOOST_HASH "SeQuant uses system Boost for hashing (thus results depend on the Boost version)")

option(SEQUANT_BUILD_DOCS "Build SeQuant documentation using Doxygen and Sphinx" OFF)
add_feature_info(BUILD_DOCS SEQUANT_BUILD_DOCS "Build SeQuant documentation. Doxygen and Sphinx are required.")

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
        SeQuant/domain/mbpt/models/cc.cpp
        SeQuant/domain/mbpt/models/cc.hpp
        SeQuant/domain/mbpt/op.cpp
        SeQuant/domain/mbpt/op.hpp
        SeQuant/domain/mbpt/op.ipp
        SeQuant/domain/mbpt/rdm.cpp
        SeQuant/domain/mbpt/rdm.hpp
        SeQuant/domain/mbpt/spin.cpp
        SeQuant/domain/mbpt/spin.hpp
        SeQuant/domain/mbpt/vac_av.ipp)

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
            SeQuant/domain/eval/eval_result.cpp
            SeQuant/domain/eval/eval_result.hpp
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

if (SEQUANT_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

### unit tests
include(FindOrFetchCatch2)
include(Catch)

if (NOT SEQUANT_TESTS)
	set(BUILD_BY_DEFAULT "EXCLUDE_FROM_ALL")
endif()

function(build_test_as_needed TARGET_NAME TEST_PREFIX DEPENDENT_TESTS)
	if (SEQUANT_TESTS)
		# Tests are built by default -> no need to add extra build steps
		return()
	endif()

	add_test("${TEST_PREFIX}/build" "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${TARGET_NAME}")

	foreach(current IN LISTS ${DEPENDENT_TESTS})
		set_tests_properties("${current}" PROPERTIES DEPENDS "${TEST_PREFIX}/build")
	endforeach()
endfunction()


set(utests_src
		tests/unit/catch2_sequant.hpp
		tests/unit/test_biorthogonalization.cpp
		tests/unit/test_space.cpp
		tests/unit/test_index.cpp
		tests/unit/test_op.cpp
		tests/unit/test_wick.cpp
		tests/unit/test_tensor.cpp
		tests/unit/test_bliss.cpp
		tests/unit/test_expr.cpp
		tests/unit/test_iterator.cpp
		tests/unit/test_mbpt.cpp
		tests/unit/test_mbpt_cc.cpp
		tests/unit/test_tensor_network.cpp
		tests/unit/test_spin.cpp
		tests/unit/test_canonicalize.cpp
		tests/unit/test_expr.cpp
		tests/unit/test_parse.cpp
		tests/unit/test_eval_expr.cpp
		tests/unit/test_binary_node.cpp
		tests/unit/test_asy_cost.cpp
		tests/unit/test_eval_node.cpp
		tests/unit/test_optimize.cpp
		tests/unit/test_fusion.cpp
		tests/unit/test_runtime.cpp
		tests/unit/test_math.cpp
		tests/unit/test_string.cpp
		tests/unit/test_latex.cpp
		tests/unit/test_utilities.cpp
		tests/unit/test_export.cpp
		tests/unit/test_meta.cpp
		)

if (TARGET tiledarray)
	list(APPEND utests_src
			tests/unit/test_eval_ta.cpp
			tests/unit/test_eval_btas.cpp
			tests/unit/test_cache_manager.cpp
			)
	set_source_files_properties(
			tests/unit/test_eval_btas.cpp tests/unit/test_eval_ta.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif (TARGET tiledarray)

set(utests_deps SeQuant Catch2::Catch2)

set(unit_test_executable unit_tests-sequant)
add_executable(${unit_test_executable} ${BUILD_BY_DEFAULT}
		tests/unit/test_main.cpp
		tests/unit/test_config.hpp
		${utests_src})
target_link_libraries(${unit_test_executable} ${utests_deps})
if (TARGET tiledarray)
	target_link_libraries(${unit_test_executable} tiledarray)
	target_compile_definitions(${unit_test_executable} PRIVATE SEQUANT_HAS_TILEDARRAY)
endif (TARGET tiledarray)

if (SEQUANT_TESTS)
	catch_discover_tests(
		${unit_test_executable}
		TEST_PREFIX "sequant/"
		WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/unit"
	)
else()
	add_test(
		NAME "sequant/unit"
		COMMAND ${unit_test_executable}
		WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/unit"
	)
	set(TEST_NAMES "sequant/unit")
	build_test_as_needed(${unit_test_executable} "sequant/unit" TEST_NAMES)
endif()


####### Integration tests aka "examples" ########

# if target examples-sequant does not exist yet, create it
if (NOT TARGET examples-sequant)
	add_custom_target_subproject(sequant examples)
endif()

# Single-Reference Coupled-Cluster equation generation
# (spin-orbital)
set(example0 srcc)
add_executable(${example0} ${BUILD_BY_DEFAULT}
		examples/${example0}/${example0}.cpp)
target_link_libraries(${example0} SeQuant)
# N.B. empty space, to appease foreach
set(${example0}_cmdline_args " ;3 t std so;3 t csv so;3 lambda std so;3 lambda csv so;2 t std sf;2 t csv sf;2 lambda std sf;2 lambda csv sf")

if (TARGET Eigen3::Eigen)  # these examples require Eigen for full functionality
  # Single-Reference closed-shell Coupled-Cluster equation
  # generation (fast spin-traced)
  set(example1 stcc)
  add_executable(${example1} ${BUILD_BY_DEFAULT}
		  examples/${example1}/${example1}.cpp)
  target_link_libraries(${example1} SeQuant)

  # Single-Reference closed-shell Coupled-Cluster equation
  # generation (rigorous spin-traced)
  set(example2 stcc_rigorous)
  add_executable(${example2} ${BUILD_BY_DEFAULT}
		examples/${example2}/${example2}.cpp)
  target_link_libraries(${example2} SeQuant)
endif()

# Single-Reference open-shell equation generation (spin-traced)
set(example3 osstcc)
add_executable(${example3} ${BUILD_BY_DEFAULT}
		examples/${example3}/${example3}.cpp)
target_link_libraries(${example3} SeQuant)

# uccf12 example moved to MPQC

set(example_eval_src
		examples/eval/eval_utils.hpp
		examples/eval/options.hpp
		examples/eval/options.cpp
		examples/eval/data_info.hpp
		examples/eval/data_info.cpp
		examples/eval/calc_info.hpp
		examples/eval/calc_info.cpp
		examples/eval/scf.hpp
		examples/eval/scf.cpp
		)

if (TARGET tiledarray)
	set(example5 eval_ta)
	add_executable(${example5} ${BUILD_BY_DEFAULT}
			${example_eval_src}
			examples/eval/ta/data_world_ta.hpp
			examples/eval/ta/scf_ta.hpp
			examples/eval/ta/main.cpp)
	target_link_libraries(${example5} SeQuant tiledarray)
endif (TARGET tiledarray)

set(example6 antisymmetrizer_test)
add_executable(${example6} ${BUILD_BY_DEFAULT}
		examples/${example6}/${example6}.cpp)
target_link_libraries(${example6} SeQuant)
target_compile_definitions(${example6} PRIVATE SEQUANT_HAS_TILEDARRAY)

if (BTAS_SOURCE_DIR)
	set(example7 eval_btas)
	add_executable(${example7} ${BUILD_BY_DEFAULT}
			${example_eval_src}
			examples/eval/btas/data_world_btas.hpp
			examples/eval/btas/scf_btas.hpp
			examples/eval/btas/main.cpp)
	target_include_directories(${example7} PUBLIC ${BTAS_SOURCE_DIR})
	target_link_libraries(${example7} SeQuant)
endif (BTAS_SOURCE_DIR)

set(lastexample 13)
foreach (i RANGE 8 ${lastexample})
  math(EXPR s "${i} - 7") # map 8..13 -> 1..6
  set(example${i} synopsis${s})
  add_executable(${example${i}} ${BUILD_BY_DEFAULT}
		  examples/synopsis/${example${i}}.cpp)
  target_link_libraries(${example${i}} SeQuant)
endforeach ()

set(example12 "tensor_network_graphs")
add_executable(${example12} ${BUILD_BY_DEFAULT}
		examples/${example12}/${example12}.cpp)
target_link_libraries(${example12} SeQuant)

set (example13 "eomcc")
add_executable(${example13} ${BUILD_BY_DEFAULT}
		examples/${example13}/${example13}.cpp)
target_link_libraries(${example13} SeQuant)

# add tests for running examples
foreach (i RANGE ${lastexample})
	if (TARGET ${example${i}})
		add_dependencies(examples-sequant ${example${i}})
		if (NOT DEFINED ${example${i}}_cmdline_args)
			set(${example${i}}_cmdline_args " ") # N.B. empty space, to appease foreach
		endif()
		set(${example${i}}_run_counter 0)
		set(VARIANTS "")
		foreach (args ${${example${i}}_cmdline_args})
			string(REPLACE " " ";" args "${args}")
			add_test(NAME sequant/example/${example${i}}/${${example${i}}_run_counter}
					 COMMAND ${example${i}} ${args}
					 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples)
			list(APPEND VARIANTS "sequant/example/${example${i}}/${${example${i}}_run_counter}")
			math(EXPR ${example${i}}_run_counter "${${example${i}}_run_counter} + 1")
		endforeach()

		build_test_as_needed(${example${i}} "sequant/example/${example${i}}" VARIANTS)
	endif ()
endforeach ()

# target for running ALL tests
add_custom_target_subproject(sequant check
		USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "^sequant")

####### Python ########

if (SEQUANT_PYTHON)
    if (NOT CMAKE_POSITION_INDEPENDENT_CODE)
        message(FATAL_ERROR "Python module requires CMAKE_POSITION_INDEPENDENT_CODE=ON")
    endif ()
    add_subdirectory(python)
endif ()

####### DOCS ########
if (SEQUANT_BUILD_DOCS)
    add_subdirectory(doc)
endif ()

####### External interface ########

add_subdirectory(external-interface)


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
