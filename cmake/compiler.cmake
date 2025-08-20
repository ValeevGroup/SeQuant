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
endfunction()
