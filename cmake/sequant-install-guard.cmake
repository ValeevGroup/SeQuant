# Install-time guard (CMake script mode, included from cmake_install.cmake): refuses to overwrite
# a different version of a package already installed in the destination.
#
# The packages SeQuant is about to install are enumerated, not hardwired: SeQuant's own install
# script is first run into a staging directory (DESTDIR), and every CMake package version file
# found there identifies a package and its version. This covers whatever the fetched
# dependencies install themselves (e.g. MADNESS with TiledArray).
#
# Expects SEQUANT_INSTALL_SCRIPT, SEQUANT_INSTALL_STAGING_DIR and
# SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES; the staging run is marked by the
# SEQUANT_INSTALL_STAGING environment variable so that the guard does not recurse.

if ("$ENV{SEQUANT_INSTALL_STAGING}")
    return()
endif ()

# <out-var> = <path> as file(INSTALL) places it under <destdir>: a Windows drive letter is dropped
function(_sequant_destdir_path destdir path out_var)
    if (destdir AND path MATCHES "^[A-Za-z]:")
        string(SUBSTRING "${path}" 2 -1 path)
    endif ()
    set(${out_var} "${destdir}${path}" PARENT_SCOPE)
endfunction()

# <out-var> = the version a CMake package version file reports in PACKAGE_VERSION
function(_sequant_package_version file out_var)
    unset(PACKAGE_FIND_VERSION)
    unset(PACKAGE_VERSION)
    include("${file}")
    set(${out_var} "${PACKAGE_VERSION}" PARENT_SCOPE)
endfunction()

# a function keeps the guard's variables out of the cmake_install.cmake scope
function(_sequant_guard_installed_dependencies)
    set(_stage "${SEQUANT_INSTALL_STAGING_DIR}")
    if (NOT _stage OR NOT IS_ABSOLUTE "${_stage}")
        message(FATAL_ERROR "SeQuant: SEQUANT_INSTALL_STAGING_DIR must be an absolute path, got '${_stage}'")
    endif ()
    file(REMOVE_RECURSE "${_stage}")

    message(STATUS "SeQuant: staging the install under ${_stage} to check for conflicts")
    set(_args "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
    foreach (_var CMAKE_INSTALL_CONFIG_NAME CMAKE_INSTALL_COMPONENT CMAKE_INSTALL_LOCAL_ONLY)
        if (DEFINED ${_var})
            list(APPEND _args "-D${_var}=${${_var}}")
        endif ()
    endforeach ()
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "DESTDIR=${_stage}" SEQUANT_INSTALL_STAGING=1
                "${CMAKE_COMMAND}" ${_args} -P "${SEQUANT_INSTALL_SCRIPT}"
        RESULT_VARIABLE _rc
        OUTPUT_QUIET
        ERROR_VARIABLE _err)
    if (NOT _rc EQUAL 0)
        message(FATAL_ERROR "SeQuant: staging the install failed (${_rc}):\n${_err}")
    endif ()

    _sequant_destdir_path("${_stage}" "${CMAKE_INSTALL_PREFIX}" _staged_prefix)
    _sequant_destdir_path("$ENV{DESTDIR}" "${CMAKE_INSTALL_PREFIX}" _dest_prefix)
    file(GLOB_RECURSE _version_files LIST_DIRECTORIES false
         "${_staged_prefix}/*-config-version.cmake" "${_staged_prefix}/*ConfigVersion.cmake")

    set(_summary "")
    set(_conflicts "")
    foreach (_staged IN LISTS _version_files)
        get_filename_component(_name "${_staged}" NAME)
        # SeQuant itself is what is being installed or upgraded
        if (_name STREQUAL "sequant-config-version.cmake")
            continue()
        endif ()
        string(REGEX REPLACE "(-config-version|ConfigVersion)\\.cmake$" "" _pkg "${_name}")
        file(RELATIVE_PATH _rel "${_staged_prefix}" "${_staged}")
        _sequant_package_version("${_staged}" _new)
        list(APPEND _summary "${_pkg} ${_new}")
        set(_installed "${_dest_prefix}/${_rel}")
        if (NOT EXISTS "${_installed}")
            continue()
        endif ()
        _sequant_package_version("${_installed}" _old)
        if (NOT _old VERSION_EQUAL _new)
            get_filename_component(_dir "${_installed}" DIRECTORY)
            list(APPEND _conflicts "${_pkg} ${_old} at ${_dir} with ${_pkg} ${_new}")
        endif ()
    endforeach ()
    file(REMOVE_RECURSE "${_stage}")

    if (_summary)
        list(JOIN _summary ", " _summary)
        message(STATUS "SeQuant installs these packages into ${_dest_prefix}: ${_summary}")
    else ()
        message(STATUS "SeQuant installs no packages besides itself into ${_dest_prefix}")
    endif ()

    if (_conflicts)
        list(JOIN _conflicts ";\n  " _conflicts)
        set(_message "SeQuant would overwrite\n  ${_conflicts}.\nSeQuant fetched its own copy because it did not find the installed one: the install prefix is not searched for dependencies unless it is on CMAKE_PREFIX_PATH (CMAKE_FIND_NO_INSTALL_PREFIX is ON by default), and a copy older than the required minimum is rejected. To proceed: upgrade or remove the installed copy, or make SeQuant use it instead of fetching (put its prefix on CMAKE_PREFIX_PATH, or set <pkg>_DIR) and reconfigure; or install SeQuant into a different prefix; or reconfigure with -DSEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES=ON.")
        if (SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES)
            message(WARNING "${_message}")
        else ()
            message(FATAL_ERROR "${_message}")
        endif ()
    endif ()
endfunction()

_sequant_guard_installed_dependencies()
