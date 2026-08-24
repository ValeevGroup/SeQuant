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



function(target_set_optimization_flags TARGET)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		return()
	endif()

	get_target_property(TARGET_TYPE "${TARGET}" TYPE)

	if (TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
		message(WARNING "target_set_optimization_flags is not intended to be used on interface targets")
		return()
	endif()

	include(CheckCXXCompilerFlag)
	include(CheckIPOSupported)

	check_ipo_supported(RESULT CMAKE_SUPPORTS_COMPILER_LTO LANGUAGES CXX)

	set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
	check_cxx_compiler_flag("-flto" SEQUANT_LTO_FLAG_SUPPORTED)
	check_cxx_compiler_flag("-flto=auto" SEQUANT_LTO_AUTO_SUPPORTED)
	check_cxx_compiler_flag("-flto;-ffat-lto-objects" SEQUANT_FAT_LTO_FLAG_SUPPORTED)

	if (DEFINED SEQUANT_LTO)
		# Always honor explicit user choice
		set(ENABLE_LTO ${SEQUANT_LTO})
	elseif(TARGET_TYPE STREQUAL "STATIC_LIBRARY" OR TARGET_TYPE STREQUAL "OBJECT_LIBRARY")
		# For static/object libraries we only want to enable LTO by default, if we can create
		# "fat" object files. Those can still be linked without LTO and hence shouldn't
		# break any downstream use.
		set(ENABLE_LTO ${SEQUANT_FAT_LTO_FLAG_SUPPORTED})
	elseif(SEQUANT_LTO_FLAG_SUPPORTED OR CMAKE_SUPPORTS_COMPILER_LTO)
		# Anything but static/object libraries is also linked by us and
		# hence enabling LTO doesn't affect downstream compatibility
		set(ENABLE_LTO ON)
	endif()

	if (ENABLE_LTO)
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

			if (SEQUANT_FAT_LTO_FLAG_SUPPORTED)
				target_compile_options("${TARGET}" PRIVATE -ffat-lto-objects)
			endif()
		else()
			set_target_properties("${TARGET}" PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
		endif()
	endif()
endfunction()
