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
    else()
        message(DEBUG "Warnings-as-errors not supported for compiler '${CMAKE_CXX_COMPILER_ID}' - disabling…")
    endif()
endfunction()

function(target_set_warning_flags TARGET)
    __check_gnu_like_compiler()

    if (NOT PROJECT_IS_TOP_LEVEL)
        if (IS_GNU_LIKE_COMPILER)
            # Disable compiler warnings
            target_compile_options("${TARGET}" PRIVATE "-w")
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
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
check_cxx_compiler_flag("-flto" SEQUANT_LTO_FLAG_SUPPORTED)
check_cxx_compiler_flag("-flto=auto" SEQUANT_LTO_AUTO_SUPPORTED)
check_cxx_compiler_flag("-flto;-ffat-lto-objects" SEQUANT_FAT_LTO_FLAG_SUPPORTED)
unset(CMAKE_TRY_COMPILE_TARGET_TYPE)

# Tri-state on purpose: ON/OFF force LTO on/off for all SeQuant targets, whereas the default
# (empty) lets target_set_optimization_flags decide per target type. Declared as a cache
# variable so that it shows up in cmake -LH, ccmake, etc. just like our other knobs.
set(SEQUANT_LTO "" CACHE STRING
	"Whether to build SeQuant's targets with link-time optimization (LTO); leave empty to decide automatically per target type")
set_property(CACHE SEQUANT_LTO PROPERTY STRINGS "" ON OFF)

# Without "fat" objects, a static library built with LTO holds IR rather than machine code, and
# an archiver that doesn't understand that IR produces an archive without a usable symbol index
# ("archive has no index" at link time). CMake substitutes the compiler's LTO-aware archiver
# when it drives LTO itself (via INTERPROCEDURAL_OPTIMIZATION), but we set the LTO flags by hand
# (see below), so we have to make sure the archiver matches. In practice CMake's own CMAKE_AR
# detection already picks e.g. llvm-ar next to clang++, so the substitution below only kicks in
# on toolchains where it doesn't - it is a safety net rather than the common path.
set(SEQUANT_LTO_AWARE_ARCHIVER FALSE)
if (NOT PROJECT_IS_TOP_LEVEL)
	# We don't apply LTO flags as a subproject (see target_set_optimization_flags below), so
	# there is no reason to touch the encompassing project's archiver either
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
	if (NOT PROJECT_IS_TOP_LEVEL)
		# When SeQuant is consumed as a subproject, how its targets are optimized is the
		# encompassing project's call, not ours (cf. target_set_warning_flags)
		return()
	endif()

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		return()
	endif()

	get_target_property(TARGET_TYPE "${TARGET}" TYPE)

	if (TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
		message(WARNING "target_set_optimization_flags is not intended to be used on interface targets")
		return()
	endif()

	if (TARGET_TYPE STREQUAL "STATIC_LIBRARY" OR TARGET_TYPE STREQUAL "OBJECT_LIBRARY")
		set(IS_ARCHIVE_LIKE_TARGET TRUE)
	else()
		set(IS_ARCHIVE_LIKE_TARGET FALSE)
	endif()

	if (NOT SEQUANT_LTO STREQUAL "")
		# Always honor explicit user choice
		set(ENABLE_LTO ${SEQUANT_LTO})
	elseif(IS_ARCHIVE_LIKE_TARGET)
		# For static/object libraries we only want to enable LTO by default, if we can create
		# "fat" object files. Those can still be linked without LTO and hence shouldn't
		# break any downstream use.
		set(ENABLE_LTO ${SEQUANT_FAT_LTO_FLAG_SUPPORTED})
	elseif(SEQUANT_LTO_FLAG_SUPPORTED OR SEQUANT_CAN_RELY_ON_CMAKE_LTO)
		# Anything but static/object libraries is also linked by us and
		# hence enabling LTO doesn't affect downstream compatibility
		set(ENABLE_LTO ON)
	endif()

	if (ENABLE_LTO)
		if (NOT SEQUANT_FAT_LTO_FLAG_SUPPORTED AND NOT SEQUANT_LTO_AWARE_ARCHIVER
				AND TARGET_TYPE STREQUAL "STATIC_LIBRARY")
			message(FATAL_ERROR "Requested LTO for static library '${TARGET}' but your compiler supports neither "
				"\"fat\" LTO objects nor an LTO-aware archiver (CMAKE_CXX_COMPILER_AR). The resulting archive "
				"would not be linkable - point CMAKE_AR and CMAKE_RANLIB at an LTO-aware ar/ranlib (e.g. llvm-ar "
				"and llvm-ranlib) or use SEQUANT_LTO=OFF")
		endif()

		if (SEQUANT_LTO_FLAG_SUPPORTED)
			# We prefer to manually set the LTO flag(s) rather than CMake doing it for us
			# due to https://gitlab.kitware.com/cmake/cmake/-/work_items/23136
			# On some compilers, the thin LTO type requested by CMake is incompatible
			# with explicitly asking for fat LTO object files.
			# Besides, it seems like full LTO achieves quite a bit better optimizations
			# with Clang.
			if (SEQUANT_LTO_AUTO_SUPPORTED)
				target_compile_options("${TARGET}" PRIVATE -flto=auto)
				target_link_options("${TARGET}" PRIVATE -flto=auto)
			else()
				target_compile_options("${TARGET}" PRIVATE -flto)
				target_link_options("${TARGET}" PRIVATE -flto)
			endif()

			# Only static/object libraries benefit from "fat" objects (see above): everything else
			# is linked by us, so the extra machine-code copy is never used - while producing it
			# costs an entire second compilation of every translation unit on GCC.
			if (SEQUANT_FAT_LTO_FLAG_SUPPORTED AND IS_ARCHIVE_LIKE_TARGET)
				target_compile_options("${TARGET}" PRIVATE -ffat-lto-objects)
			endif()
		else()
			if (NOT SEQUANT_CAN_RELY_ON_CMAKE_LTO)
				message(FATAL_ERROR "Requested LTO but CMake doesn't know how to enable it for your compiler - Use SEQUANT_LTO=OFF")
			endif()
			set_target_properties("${TARGET}" PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
		endif()
	endif()
endfunction()
