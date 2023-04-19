cmake_minimum_required(VERSION 3.14)

project(couette-solver-lib LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRC_FOLDER "${PROJECT_SOURCE_DIR}/src")
set(HDR_FOLDER "${PROJECT_SOURCE_DIR}/hdr")


set(${PROJECT_NAME}_SRC  ${SRC_FOLDER}/abstractsolver.cpp
                         ${SRC_FOLDER}/hllcsolver.cpp
                         ${SRC_FOLDER}/godunovsolver.cpp
			 ${SRC_FOLDER}/coeffsolver.cpp
                         ${SRC_FOLDER}/additionalsolver.cpp
                         ${SRC_FOLDER}/mixture.cpp
                         ${SRC_FOLDER}/bordercondition.cpp
                         ${SRC_FOLDER}/datawriter.cpp
)

set(${PROJECT_NAME}_HDR  ${HDR_FOLDER}/abstractsolver.h
                         ${HDR_FOLDER}/hllcsolver.h
                         ${HDR_FOLDER}/godunovsolver.h
                         ${HDR_FOLDER}/coeffsolver.h
                         ${HDR_FOLDER}/additionalsolver.h
                         ${HDR_FOLDER}/global.h    
                         ${HDR_FOLDER}/mixture.h
                         ${HDR_FOLDER}/bordercondition.h
                         ${HDR_FOLDER}/datawriter.h
)

add_subdirectory(example)

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SRC} ${${PROJECT_NAME}_HDR})

target_compile_definitions(${PROJECT_NAME} PRIVATE -D_USE_MATH_DEFINES)

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")


target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/hdr>"
                                                  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

