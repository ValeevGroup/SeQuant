# Standalone project that consumes an installed SeQuant via find_package;
# not part of the SeQuant build, see the consumer steps in .github/workflows/cmake.yml
cmake_minimum_required(VERSION 3.28)
project(sequant-consumer CXX)

set(SEQUANT_CONSUMER_BACKENDS "" CACHE STRING "Exported SeQuant eval backends to link as well; subset of ta;btas;tapp")

find_package(SeQuant CONFIG REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE SeQuant::SeQuant)
foreach (_backend IN LISTS SEQUANT_CONSUMER_BACKENDS)
    if (NOT TARGET SeQuant::eval::${_backend})
        message(FATAL_ERROR "the installed SeQuant does not export SeQuant::eval::${_backend}")
    endif ()
    string(TOUPPER ${_backend} _backend_uc)
    target_link_libraries(consumer PRIVATE SeQuant::eval::${_backend})
    target_compile_definitions(consumer PRIVATE SEQUANT_CONSUMER_${_backend_uc})
endforeach ()

# Windows has no RPATH: an executable finds DLLs in its own directory (or via PATH), so copy the ones it needs
# (SeQuant's and its dependencies') next to it
get_target_property(_sequant_symb_type SeQuant::symb TYPE)
if (WIN32 AND _sequant_symb_type STREQUAL "SHARED_LIBRARY")
    add_custom_command(TARGET consumer POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:consumer> $<TARGET_FILE_DIR:consumer>
            COMMAND_EXPAND_LISTS)
endif ()
