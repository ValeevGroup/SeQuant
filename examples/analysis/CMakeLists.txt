FetchContent_Declare(
        clipp
        QUIET
        GIT_REPOSITORY https://github.com/GerHobbelt/clipp.git
        GIT_TAG 02783b6782ebedbb2bebc2e6ceda738ee51c7df2
)

FetchContent_Declare(
        fmt
        QUIET
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.1.1
)

FetchContent_MakeAvailable(clipp)
FetchContent_MakeAvailable(fmt)

add_executable(sqntly EXCLUDE_FROM_ALL
        main.cpp
        utils.hpp
        utils.cpp
        cmdline.hpp
        cmdline.cpp
        imed.cpp
        imed.hpp
        algorithm.hpp
        algorithm.cpp
)

target_link_libraries(sqntly PRIVATE
        SeQuant::SeQuant
        clipp::clipp
        fmt::fmt
)