cmake_minimum_required(VERSION 3.17)
project(g4interface LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)
find_package(BDSIM REQUIRED HINTS $ENV{BDSIM})

set(SOURCES
    BDSXtrackInterface.cpp
    bindings.cpp
)

pybind11_add_module(g4interface ${SOURCES})

target_include_directories(g4interface
    PRIVATE
        ${BDSIM_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(g4interface PRIVATE ${BDSIM_LIBRARIES})
