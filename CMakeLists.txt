cmake_minimum_required(VERSION 3.5.0)
project(hw4-numerical-modeling VERSION 0.1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 17)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

add_executable(hw4-numerical-modeling main.cpp)
link_directories("/opt/homebrew/lib")
include_directories("/opt/homebrew/include")

target_link_libraries(hw4-numerical-modeling wave)

