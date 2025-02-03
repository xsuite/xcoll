cmake_minimum_required(VERSION 3.17)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O0")

project(g4interface)

add_subdirectory(lib/pybind11)

set(SOURCE_DIR "src/collimasim")

# The headers are in the source directory
include_directories(${SOURCE_DIR})

set(SOURCES ${SOURCE_DIR}/BDSXtrackInterface.cpp
            ${SOURCE_DIR}/BDSXtrackInterface.hh
            ${SOURCE_DIR}/bindings.cpp)

find_package(BDSIM REQUIRED HINTS $ENV{BDSIM})

pybind11_add_module(g4interface ${SOURCES} "${SOURCE_DIR}/bindings.cpp")


target_link_libraries(g4interface PRIVATE ${BDSIM_LIBRARIES})
target_include_directories(g4interface PUBLIC ${BDSIM_INCLUDE_DIR})
