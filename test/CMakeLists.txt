cmake_minimum_required(VERSION 3.0)

project(wbc_test)

set(CMAKE_BUILD_TYPE "Debug") # Debug Release

set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

# include_directories("src")
add_executable(wqp_test wqp_test.cpp)
target_link_libraries(wqp_test
    wbckits
)

add_executable(rhp_test rhp_test.cpp)
target_link_libraries(rhp_test
    wbckits
)

add_executable(hqp_test hqp_test.cpp)
target_link_libraries(hqp_test
    wbckits
)

add_executable(ghc_test ghc_test.cpp)
target_link_libraries(ghc_test
    wbckits
)