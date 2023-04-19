cmake_minimum_required(VERSION 3.14)

set(PROJECT_MAIN "${PROJECT_NAME}_main")
project(${PROJECT_MAIN} LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#add_executable(${PROJECT_MAIN} couette-solver-lib_HLLCmain.cpp)
add_executable(${PROJECT_MAIN} couette-solver-lib_Godunovmain.cpp)
target_link_libraries(${PROJECT_MAIN} couette-solver-lib)
