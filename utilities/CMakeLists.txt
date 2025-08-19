set(EXECUTABLE_SOURCES
    "tensor_network_graphs.cpp"
)

foreach(current IN LISTS EXECUTABLE_SOURCES)
    string(REGEX REPLACE "\\.cpp$" "" name "${current}")
    add_executable(${name} ${current})
    target_link_libraries(${name} SeQuant)
    
    if (SEQUANT_WARNINGS_AS_ERRORS)
        target_warnings_as_errors(${name})
    endif()
endforeach()

add_subdirectory(external-interface)
