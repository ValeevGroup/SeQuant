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

# records a fetched dependency that will be installed; CONFIG_VERSION_FILE is relative to the install prefix
function(sequant_record_installed_dependency pkg)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "VERSION;CONFIG_VERSION_FILE" "")
    if (NOT ARG_VERSION OR NOT ARG_CONFIG_VERSION_FILE)
        message(FATAL_ERROR "sequant_record_installed_dependency(${pkg}): VERSION and CONFIG_VERSION_FILE must be non-empty")
    endif ()
    set_property(GLOBAL APPEND PROPERTY SEQUANT_INSTALLED_DEPENDENCY_NAMES "${pkg}")
    set_property(GLOBAL APPEND PROPERTY SEQUANT_INSTALLED_DEPENDENCY_VERSIONS "${ARG_VERSION}")
    set_property(GLOBAL APPEND PROPERTY SEQUANT_INSTALLED_DEPENDENCY_CONFIG_VERSION_FILES "${ARG_CONFIG_VERSION_FILE}")
endfunction()

# reads the version recorded in <file> (a CMake package version file), sets <out-var> to it in the caller
function(_sequant_installed_package_version file out_var)
    unset(PACKAGE_FIND_VERSION)
    unset(PACKAGE_VERSION)
    include("${file}")
    set(${out_var} "${PACKAGE_VERSION}" PARENT_SCOPE)
endfunction()

# summarizes what will be installed, warns about different versions already installed,
# and writes the manifest consumed at install time
function(sequant_report_fetched_dependencies)
    get_property(names GLOBAL PROPERTY SEQUANT_INSTALLED_DEPENDENCY_NAMES)
    get_property(versions GLOBAL PROPERTY SEQUANT_INSTALLED_DEPENDENCY_VERSIONS)
    get_property(files GLOBAL PROPERTY SEQUANT_INSTALLED_DEPENDENCY_CONFIG_VERSION_FILES)

    if (names)
        set(summary "")
        set(i 0)
        foreach (name IN LISTS names)
            list(GET versions ${i} version)
            list(GET files ${i} file)
            string(APPEND summary "${name} ${version}, ")
            if (EXISTS "${CMAKE_INSTALL_PREFIX}/${file}")
                _sequant_installed_package_version("${CMAKE_INSTALL_PREFIX}/${file}" installed)
                if (NOT installed VERSION_EQUAL version)
                    message(WARNING "${name} ${installed} is already installed at ${CMAKE_INSTALL_PREFIX}/${file} and differs from the ${version} that will be installed; the install will refuse to overwrite it unless SEQUANT_OVERWRITE_INSTALLED_DEPENDENCIES=ON")
                endif ()
            endif ()
            math(EXPR i "${i} + 1")
        endforeach ()
        string(REGEX REPLACE ", $" "" summary "${summary}")
        message(STATUS "SeQuant will install these fetched dependencies into ${CMAKE_INSTALL_PREFIX}: ${summary}")
    else ()
        message(STATUS "SeQuant will install no fetched dependencies")
    endif ()

    file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/sequant-fetched-dependencies.cmake"
         CONTENT "set(SEQUANT_FETCHED_DEPENDENCY_NAMES \"${names}\")
set(SEQUANT_FETCHED_DEPENDENCY_VERSIONS \"${versions}\")
set(SEQUANT_FETCHED_DEPENDENCY_CONFIG_VERSION_FILES \"${files}\")
")
endfunction()
