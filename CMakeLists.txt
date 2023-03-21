cmake_minimum_required(VERSION 3.14)

project(couette-solver-lib LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRC_FOLDER "${PROJECT_SOURCE_DIR}/src")
set(HDR_FOLDER "${PROJECT_SOURCE_DIR}/hdr")


set(${PROJECT_NAME}_SRC  ${SRC_FOLDER}/AbstractSolver.cpp
                         ${SRC_FOLDER}/AdditionalSolver.cpp
                         ${SRC_FOLDER}/Mixture.cpp
)

set(${PROJECT_NAME}_HDR  ${HDR_FOLDER}/AbstractSolver.h
                         ${HDR_FOLDER}/AdditionalSolver.h
                         ${HDR_FOLDER}/global.h    
			 ${HDR_FOLDER}/Mixture.h                     
)

add_subdirectory(example)

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SRC} ${${PROJECT_NAME}_HDR})

target_compile_definitions(${PROJECT_NAME} PRIVATE -D_USE_MATH_DEFINES)

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/hdr>"
                                                  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

