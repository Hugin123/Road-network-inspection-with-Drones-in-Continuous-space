cmake_minimum_required(VERSION 3.28)
project(bishe)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(bishe main.cpp
        model.cpp
        model.h
        init_data.cpp
        init_data.h)

set(GUROBI_LIB_DIR "$ENV{GUROBI_HOME}\\lib")
set(GUROBI_INCLUDE_DIR "$ENV{GUROBI_HOME}\\include")

message(STATUS "GUROBI_LIB_DIR: ${GUROBI_LIB_DIR}")
message(STATUS "GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")

link_directories(${GUROBI_LIB_DIR})
include_directories(${GUROBI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} "${GUROBI_LIB_DIR}\\gurobi110.lib" "${GUROBI_LIB_DIR}\\gurobi_c++mdd2017.lib")