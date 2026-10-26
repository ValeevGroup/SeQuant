# Options and helpers for dependencies that SeQuant fetches (FetchContent) and may install.
# Included from the top-level CMakeLists.txt before any FindOrFetch* module.

option(SEQUANT_FETCH_MISSING_DEPENDENCIES "Fetch (via FetchContent) dependencies that are not found; if OFF, missing dependencies are an error" ON)
add_feature_info(SEQUANT_FETCH_MISSING_DEPENDENCIES SEQUANT_FETCH_MISSING_DEPENDENCIES "Fetch missing dependencies")
option(SEQUANT_INSTALL_FETCHED_DEPENDENCIES "Install the fetched primary dependencies (Range-v3, libPerm, polymorphic_variant, ...) with SeQuant" ON)
add_feature_info(SEQUANT_INSTALL_FETCHED_DEPENDENCIES SEQUANT_INSTALL_FETCHED_DEPENDENCIES "Install fetched primary dependencies")
set(SEQUANT_INSTALL_FETCHED_BACKENDS "" CACHE STRING "Fetched eval backends to install with SeQuant; subset of TiledArray;BTAS;TAPP (case-insensitive)")
option(SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES "Allow installing a fetched dependency over a different version already installed in the prefix" OFF)
add_feature_info(SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES "Overwrite different versions of installed dependencies")

# validate SEQUANT_INSTALL_FETCHED_BACKENDS, normalizing the spelling of each entry
set(_sequant_known_backends TiledArray BTAS TAPP)
# the SeQuant module of each backend
set(_sequant_backend_module_TiledArray SeQuant-eval-ta)
set(_sequant_backend_module_BTAS SeQuant-eval-btas)
set(_sequant_backend_module_TAPP SeQuant-eval-tapp)
set(SEQUANT_INSTALL_FETCHED_BACKENDS_NORMALIZED "")
foreach (_backend IN LISTS SEQUANT_INSTALL_FETCHED_BACKENDS)
    set(_normalized "")
    foreach (_known IN LISTS _sequant_known_backends)
        string(TOLOWER "${_backend}" _backend_lc)
        string(TOLOWER "${_known}" _known_lc)
        if (_backend_lc STREQUAL _known_lc)
            set(_normalized ${_known})
        endif ()
    endforeach ()
    if (NOT _normalized)
        message(FATAL_ERROR "SEQUANT_INSTALL_FETCHED_BACKENDS: unknown backend '${_backend}'; allowed: TiledArray;BTAS;TAPP")
    endif ()
    list(APPEND SEQUANT_INSTALL_FETCHED_BACKENDS_NORMALIZED ${_normalized})
endforeach ()
list(REMOVE_DUPLICATES SEQUANT_INSTALL_FETCHED_BACKENDS_NORMALIZED)

# if fetching is forbidden, require <pkg> to be found; call before FetchContent_Declare
macro(sequant_find_unless_fetching pkg)
    if (NOT SEQUANT_FETCH_MISSING_DEPENDENCIES)
        find_package(${pkg} ${ARGN} REQUIRED)
    endif ()
endmacro()

# for dependencies that cannot be found, only fetched
macro(sequant_fetch_only_dependency pkg advice)
    if (NOT SEQUANT_FETCH_MISSING_DEPENDENCIES)
        message(FATAL_ERROR "${pkg} can only be fetched, but SEQUANT_FETCH_MISSING_DEPENDENCIES=OFF; ${advice}")
    endif ()
endmacro()

# sets <out-var> to TRUE iff <backend> is to be installed when fetched
function(sequant_backend_installable backend out_var)
    if (backend IN_LIST SEQUANT_INSTALL_FETCHED_BACKENDS_NORMALIZED)
        set(${out_var} TRUE PARENT_SCOPE)
    else ()
        set(${out_var} FALSE PARENT_SCOPE)
    endif ()
endfunction()

# records a fetched dependency that SeQuant will install, for the configure-time summary;
# what that dependency installs in turn (e.g. MADNESS with TiledArray) is enumerated at
# install time by cmake/sequant-install-guard.cmake
function(sequant_record_installed_dependency pkg)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "VERSION" "")
    if (NOT ARG_VERSION)
        message(FATAL_ERROR "sequant_record_installed_dependency(${pkg}): VERSION must be non-empty")
    endif ()
    set_property(GLOBAL APPEND PROPERTY SEQUANT_INSTALLED_DEPENDENCIES "${pkg} ${ARG_VERSION}")
endfunction()

# summarizes which fetched dependencies and backends will (not) be installed
function(sequant_report_fetched_dependencies)
    get_property(installed GLOBAL PROPERTY SEQUANT_INSTALLED_DEPENDENCIES)
    if (installed)
        list(JOIN installed ", " installed)
        message(STATUS "SeQuant will install these fetched dependencies (and what they install themselves) into ${CMAKE_INSTALL_PREFIX}: ${installed}; the install refuses to overwrite a different version of any package unless SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES=ON")
    else ()
        message(STATUS "SeQuant will install no fetched dependencies")
    endif ()

    set(not_installed "")
    foreach (backend IN LISTS _sequant_known_backends)
        string(TOUPPER ${backend} backend_uc)
        if (SEQUANT_HAS_${backend_uc} AND SEQUANT_FETCHED_${backend_uc} AND NOT SEQUANT_EXPORTS_${backend_uc})
            list(APPEND not_installed ${backend})
        endif ()
    endforeach ()
    if (not_installed)
        list(JOIN not_installed ", " not_installed)
        message(STATUS "SeQuant will not install or export these fetched backends: ${not_installed} (add them to SEQUANT_INSTALL_FETCHED_BACKENDS to install them)")
    endif ()
endfunction()
