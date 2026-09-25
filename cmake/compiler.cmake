macro(__check_gnu_like_compiler)
    if (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        set(IS_GNU_LIKE_COMPILER TRUE)
    else()
        set(IS_GNU_LIKE_COMPILER FALSE)
    endif()
endmacro()

function(target_warnings_as_errors TARGET)
    __check_gnu_like_compiler()

    if (IS_GNU_LIKE_COMPILER)
        target_compile_options("${TARGET}" PRIVATE "-Werror")
    elseif(MSVC)
        target_compile_options("${TARGET}" PRIVATE "/WX")
    else()
        message(DEBUG "Warnings-as-errors not supported for compiler '${CMAKE_CXX_COMPILER_ID}' - disabling…")
    endif()
endfunction()

function(target_set_compiler_flags TARGET)
    __check_gnu_like_compiler()

    if (MSVC)
        # By default MSVC is not standard-compliant in its preprocessor implementation
        # but we need it to be (partially in headers, which is why this is a public option)
        target_compile_options("${TARGET}" PUBLIC "/Zc:preprocessor")
        # By default MSVC does not set the __cplusplus macro to the correct value
        # That breaks any code that tries to be compatible with different C++ standards
        target_compile_options("${TARGET}" PUBLIC "/Zc:__cplusplus")
        # Don't error due to object files being too big
        target_compile_options("${TARGET}" PRIVATE "/bigobj")
        # Make MSVC use and understand UTF-8 encoding in source files
        target_compile_options("${TARGET}" PUBLIC "/utf-8")

        # Increase the available stack memory to what appears to be the default on Linux/macOS (8MB)
        target_link_options("${TARGET}" PRIVATE "/STACK:8388608")
    endif()

    if (NOT PROJECT_IS_TOP_LEVEL)
        # Disable compiler warnings
        if (IS_GNU_LIKE_COMPILER)
            target_compile_options("${TARGET}" PRIVATE "-w")
        elseif(MSVC)
            target_compile_options("${TARGET}" PRIVATE "/w")
        endif()

        return()
    endif()

    if (SEQUANT_WARNINGS_AS_ERRORS)
        target_warnings_as_errors("${TARGET}")
    endif()

    if (IS_GNU_LIKE_COMPILER)
        target_compile_options("${TARGET}" PRIVATE "-Wall" "-Wpedantic" "-Wextra" "-Wno-sign-conversion" "-Wno-sign-compare" "-Wno-parentheses")
    endif()

    if (CMAKE_COMPILER_IS_GNUCXX)
        # Certain kinds of warnings are no longer suppressed inside system headers (under all circumstances) when using GCC 12+
        # Hence, we have to ensure we're not causing a compile error for those warnings as the warning might
        # be in a dependency which we can't fix.
        # Which warnings belong into this category is unclear as of writing this, so consider the below an incomplete list
        #
        # See also:
        # - https://gcc.gnu.org/bugzilla/show_bug.cgi?id=119388
        # - https://gcc.gnu.org/cgit/gcc/commit/?id=6feb628a706e86eb3f303aff388c74bdb29e7381
        # - https://stackoverflow.com/q/79742311

        # With maybe-uninitialized in particular, it appears as if it creates a lot of false positives
        # causing it to effectively only create noise. Hence, we disable it entirely.
        target_compile_options("${TARGET}" PRIVATE "-Wno-maybe-uninitialized")
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")
        # This warning can be a bit odd in that it seems like some Clang versions emit it incorrectly,
        # others don't emit it and some emit it correctly but in places where fixing the code causes
        # it to no longer be compilable with other compilers (in particular GCC) because support for
        # the exact semantics of when a lambda capture is required seems to be quite lacking across
        # different compilers.
        # See also https://github.com/llvm/llvm-project/issues/35017
        target_compile_options("${TARGET}" PRIVATE "-Wno-unused-lambda-capture")
    endif()
endfunction()


include(CheckIPOSupported)
include(CheckCXXCompilerFlag)

check_ipo_supported(RESULT SEQUANT_CAN_RELY_ON_CMAKE_LTO LANGUAGES CXX)

# Note: these are probed at directory scope (rather than inside target_set_optimization_flags)
# because the archiver selection below has to happen at directory scope in order to be inherited
# by the subdirectories that define our targets. The results are cached, so this costs a single
# round of try-compiles for the entire project.
# The compile-only probe type is scoped to these three checks: a toolchain file or an encompassing
# project may have set CMAKE_TRY_COMPILE_TARGET_TYPE for a reason (typically cross-compiling without
# a usable link step), and since this runs at directory scope every later try_compile in SeQuant
# would inherit whatever we leave behind.
if (DEFINED CMAKE_TRY_COMPILE_TARGET_TYPE)
	set(_sequant_saved_try_compile_target_type "${CMAKE_TRY_COMPILE_TARGET_TYPE}")
endif()
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
check_cxx_compiler_flag("-flto" SEQUANT_LTO_FLAG_SUPPORTED)
check_cxx_compiler_flag("-flto=auto" SEQUANT_LTO_AUTO_SUPPORTED)
check_cxx_compiler_flag("-flto;-ffat-lto-objects" SEQUANT_FAT_LTO_FLAG_SUPPORTED)
if (DEFINED _sequant_saved_try_compile_target_type)
	set(CMAKE_TRY_COMPILE_TARGET_TYPE "${_sequant_saved_try_compile_target_type}")
	unset(_sequant_saved_try_compile_target_type)
else()
	unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
endif()

# Tri-state on purpose: ON/OFF force LTO on/off for all SeQuant targets, whereas the default
# (empty) lets target_set_optimization_flags decide per target type. Declared as a cache
# variable so that it shows up in cmake -LH, ccmake, etc. just like our other knobs.
set(SEQUANT_LTO "" CACHE STRING
	"Whether to build SeQuant's targets with link-time optimization (LTO); leave empty to decide automatically per target type")
set_property(CACHE SEQUANT_LTO PROPERTY STRINGS "" ON OFF)

# Whether target_set_optimization_flags may end up applying LTO to any of our targets. An explicit
# SEQUANT_LTO is honored either way, whereas the automatic per-target defaults only kick in when
# SeQuant is the top-level project: as a subproject, whether to pay LTO's build-time and link-memory
# cost is the encompassing project's call to make, not ours.
if (NOT SEQUANT_LTO STREQUAL "")
	set(SEQUANT_MAY_USE_LTO ${SEQUANT_LTO})
else()
	set(SEQUANT_MAY_USE_LTO ${PROJECT_IS_TOP_LEVEL})
endif()

# Without "fat" objects, a static library built with LTO holds IR rather than machine code, and
# an archiver that doesn't understand that IR produces an archive without a usable symbol index
# ("archive has no index" at link time). CMake substitutes the compiler's LTO-aware archiver
# when it drives LTO itself (via INTERPROCEDURAL_OPTIMIZATION), but we set the LTO flags by hand
# (see below), so we have to make sure the archiver matches. In practice CMake's own CMAKE_AR
# detection already picks e.g. llvm-ar next to clang++, so the substitution below only kicks in
# on toolchains where it doesn't - it is a safety net rather than the common path.
set(SEQUANT_LTO_AWARE_ARCHIVER FALSE)
if (NOT SEQUANT_MAY_USE_LTO)
	# We won't be applying LTO flags to anything (see target_set_optimization_flags below), so
	# there is no reason to touch the archiver either
elseif (CMAKE_CXX_COMPILER_AR)
	set(SEQUANT_LTO_AWARE_ARCHIVER TRUE)
	if (NOT CMAKE_AR STREQUAL CMAKE_CXX_COMPILER_AR)
		# Announced rather than done silently: CMAKE_AR may well have been set deliberately
		message(STATUS "Using the compiler's LTO-aware archiver (${CMAKE_CXX_COMPILER_AR}) for SeQuant's static libraries")
		set(CMAKE_AR "${CMAKE_CXX_COMPILER_AR}")
		if (CMAKE_CXX_COMPILER_RANLIB)
			set(CMAKE_RANLIB "${CMAKE_CXX_COMPILER_RANLIB}")
		endif()
	endif()
elseif (APPLE)
	# On Apple platforms LTO objects are bitcode wrapped in a Mach-O container, which cctools'
	# ar and ranlib index just fine - no special archiver needed
	set(SEQUANT_LTO_AWARE_ARCHIVER TRUE)
endif()

function(target_set_optimization_flags TARGET)
	# CMAKE_BUILD_TYPE only reflects the active configuration for a single-config
	# generator (Ninja, Makefiles); for a multi-config generator (Visual Studio,
	# Ninja Multi-Config) it is always empty at configure time, since one configure
	# step services every configuration at once. Bail out only if EVERY configuration
	# this generator will ever build is Debug (matching the single-config check this
	# replaces); otherwise proceed and let the per-config guards further down (generator
	# expressions for hand-set flags, INTERPROCEDURAL_OPTIMIZATION_<CONFIG> for the
	# CMake-driven path) keep Debug itself LTO-free among the configurations that do want it.
	if (CMAKE_CONFIGURATION_TYPES)
		set(_seq_any_nondebug_config FALSE)
		foreach(_seq_config IN LISTS CMAKE_CONFIGURATION_TYPES)
			if (NOT _seq_config STREQUAL "Debug")
				set(_seq_any_nondebug_config TRUE)
			endif()
		endforeach()
	else()
		if (CMAKE_BUILD_TYPE STREQUAL "Debug")
			set(_seq_any_nondebug_config FALSE)
		else()
			set(_seq_any_nondebug_config TRUE)
		endif()
	endif()
	if (NOT _seq_any_nondebug_config)
		return()
	endif()

	get_target_property(TARGET_TYPE "${TARGET}" TYPE)

	if (TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
		message(WARNING "target_set_optimization_flags is not intended to be used on interface targets")
		return()
	endif()

	# Only static libraries leave the build as-is and may get linked by someone else (with or without
	# LTO); everything else - executables, shared/module libraries, but also object libraries, whose
	# objects only ever end up in a link step of ours - is linked by us.
	if (TARGET_TYPE STREQUAL "STATIC_LIBRARY")
		set(IS_ARCHIVE_LIKE_TARGET TRUE)
	else()
		set(IS_ARCHIVE_LIKE_TARGET FALSE)
	endif()

	if (NOT SEQUANT_LTO STREQUAL "")
		# Always honor explicit user choice - including when we are consumed as a subproject
		set(ENABLE_LTO ${SEQUANT_LTO})
	elseif(NOT PROJECT_IS_TOP_LEVEL)
		# Absent an explicit choice, how SeQuant's targets are optimized is the encompassing
		# project's call when we are consumed as a subproject, not ours
		set(ENABLE_LTO OFF)
	elseif(IS_ARCHIVE_LIKE_TARGET)
		# For static libraries we only want to enable LTO by default, if we can create
		# "fat" object files. Those can still be linked without LTO and hence shouldn't
		# break any downstream use.
		set(ENABLE_LTO ${SEQUANT_FAT_LTO_FLAG_SUPPORTED})
	elseif(SEQUANT_LTO_FLAG_SUPPORTED OR SEQUANT_CAN_RELY_ON_CMAKE_LTO)
		# Anything but static libraries is also linked by us and
		# hence enabling LTO doesn't affect downstream compatibility
		set(ENABLE_LTO ON)
	endif()

	if (ENABLE_LTO)
		if (SEQUANT_LTO_FLAG_SUPPORTED)
			# The archiver concern below only exists because we drive LTO by hand: when CMake does it
			# (the INTERPROCEDURAL_OPTIMIZATION branch) it also substitutes a matching archiver itself
			if (NOT SEQUANT_FAT_LTO_FLAG_SUPPORTED AND NOT SEQUANT_LTO_AWARE_ARCHIVER
					AND TARGET_TYPE STREQUAL "STATIC_LIBRARY")
				message(FATAL_ERROR "Requested LTO for static library '${TARGET}' but your compiler supports neither "
					"\"fat\" LTO objects nor an LTO-aware archiver (CMAKE_CXX_COMPILER_AR). The resulting archive "
					"would not be linkable - point CMAKE_AR and CMAKE_RANLIB at an LTO-aware ar/ranlib (e.g. llvm-ar "
					"and llvm-ranlib) or use SEQUANT_LTO=OFF")
			endif()

			# We prefer to manually set the LTO flag(s) rather than CMake doing it for us
			# due to https://gitlab.kitware.com/cmake/cmake/-/work_items/23136
			# On some compilers, the thin LTO type requested by CMake is incompatible
			# with explicitly asking for fat LTO object files.
			# Besides, it seems like full LTO achieves quite a bit better optimizations
			# with Clang.
			# Wrapped in $<CONFIG:Debug>'s negation (rather than relying on the early
			# return above alone) so that under a multi-config generator, a Debug build
			# stays LTO-free even though some other configuration here wants it.
			if (SEQUANT_LTO_AUTO_SUPPORTED)
				target_compile_options("${TARGET}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-flto=auto>")
				target_link_options("${TARGET}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-flto=auto>")
			else()
				target_compile_options("${TARGET}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-flto>")
				target_link_options("${TARGET}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-flto>")
			endif()

			# Only static libraries benefit from "fat" objects (see above): everything else
			# is linked by us, so the extra machine-code copy is never used - while producing it
			# costs an entire second compilation of every translation unit on GCC.
			if (SEQUANT_FAT_LTO_FLAG_SUPPORTED AND IS_ARCHIVE_LIKE_TARGET)
				target_compile_options("${TARGET}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-ffat-lto-objects>")
			endif()
		else()
			if (NOT SEQUANT_CAN_RELY_ON_CMAKE_LTO)
				message(FATAL_ERROR "Requested LTO but CMake doesn't know how to enable it for your compiler - Use SEQUANT_LTO=OFF")
			endif()
			# INTERPROCEDURAL_OPTIMIZATION (unsuffixed) is a fallback CMake applies to
			# EVERY configuration under a multi-config generator, Debug included; set the
			# per-configuration property for each non-Debug configuration instead so Debug
			# is left alone there. A single-config generator has no per-config properties
			# to set (CMAKE_CONFIGURATION_TYPES is empty) and takes the plain property,
			# exactly as before -- this function already returned above whenever that one
			# configuration is Debug.
			if (CMAKE_CONFIGURATION_TYPES)
				foreach(_seq_config IN LISTS CMAKE_CONFIGURATION_TYPES)
					if (NOT _seq_config STREQUAL "Debug")
						string(TOUPPER "${_seq_config}" _seq_config_upper)
						set_target_properties("${TARGET}" PROPERTIES
							"INTERPROCEDURAL_OPTIMIZATION_${_seq_config_upper}" ON)
					endif()
				endforeach()
			else()
				set_target_properties("${TARGET}" PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
			endif()
		endif()
	endif()
endfunction()
