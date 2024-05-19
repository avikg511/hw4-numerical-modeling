cmake_minimum_required(VERSION 3.5.0)
project(hw4-numerical-modeling VERSION 0.1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 17)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 

    set(SOURCES
        src/models.cpp
        src/main.cpp
    )
    add_executable(${PROJECT_NAME} ${SOURCES})

    target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${PROJECT_SOURCE_DIR}/inc
    )   
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()


if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(my_executable PRIVATE -Wno-deprecated-declarations)
endif()
 

