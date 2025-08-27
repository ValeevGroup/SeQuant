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
        target_compile_options("${TARGET}" PRIVATE "-Wno-error=maybe-uninitialized")
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
