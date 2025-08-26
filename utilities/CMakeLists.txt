set(EXECUTABLE_SOURCES
    "tensor_network_graphs.cpp"
)

foreach(current IN LISTS EXECUTABLE_SOURCES)
    string(REGEX REPLACE "\\.cpp$" "" name "${current}")
    add_executable(${name} ${current})
    target_link_libraries(${name} SeQuant)
    
    target_set_warning_flags("${name}")
endforeach()

add_subdirectory(external-interface)
