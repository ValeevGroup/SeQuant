set(EXECUTABLE_SOURCES
    "tensor_network_graphs.cpp"
)

foreach(current IN LISTS EXECUTABLE_SOURCES)
    string(REGEX REPLACE "\\.cpp$" "" name "${current}")
    add_executable(${name} ${current})
    set_target_properties(${name} PROPERTIES CXX_SCAN_FOR_MODULES OFF)
    target_link_libraries(${name} SeQuant)
    
    target_set_warning_flags("${name}")
endforeach()

add_subdirectory(external-interface)
