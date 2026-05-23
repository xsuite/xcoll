cmake_minimum_required(VERSION 3.21)  # 3.21+ recommended for Python3_AddLibrary helpers
project(g4interface LANGUAGES CXX)

message(STATUS "----------------------------------------")
message(STATUS "Compilers being used for ${PROJECT_NAME}:")
message(STATUS "  C compiler:      ${CMAKE_C_COMPILER}")
message(STATUS "  C++ compiler:    ${CMAKE_CXX_COMPILER}")
if(CMAKE_Fortran_COMPILER)
  message(STATUS "  Fortran compiler: ${CMAKE_Fortran_COMPILER}")
endif()
message(STATUS "----------------------------------------")


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure the same Python that builds also provides headers & module helper
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

# ---------- pybind11: find or vendor ----------
# Controls (cache so users can override at configure time)
option(PYBIND11_TRY_PYTHON_HINT "Try Python to discover pybind11 cmake dir" ON)
option(PYBIND11_USE_FETCHCONTENT_FALLBACK "Vendor pybind11 if not found" ON)
set(PYBIND11_FC_TAG "v3.0.1" CACHE STRING "pybind11 tag/commit for FetchContent fallback")

# 1) Try a normal config package
find_package(pybind11 CONFIG QUIET)

# 2) If not found, ask the *active* Python for pybind11's cmake dir
if (NOT pybind11_FOUND AND PYBIND11_TRY_PYTHON_HINT)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if (Python3_Interpreter_FOUND)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import pybind11, sys; sys.stdout.write(pybind11.get_cmake_dir())"
      OUTPUT_VARIABLE _pybind11_cmake_dir
      RESULT_VARIABLE _pybind11_dir_status
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if (_pybind11_dir_status EQUAL 0 AND _pybind11_cmake_dir)
      list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmake_dir}")
      find_package(pybind11 CONFIG QUIET)
    endif()
  endif()
endif()

# 3) If still not found, vendor with FetchContent
if (NOT pybind11_FOUND AND PYBIND11_USE_FETCHCONTENT_FALLBACK)
  include(FetchContent)
  message(STATUS "pybind11 not found; vendoring via FetchContent (${PYBIND11_FC_TAG})")
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        ${PYBIND11_FC_TAG}
    FIND_PACKAGE_ARGS # enables find_package(pybind11) semantics if you need them later
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# Defensive check: we should now have pybind11
if (NOT COMMAND pybind11_add_module)
  message(FATAL_ERROR "pybind11 was not found and vendoring disabled. "
                      "Enable PYBIND11_USE_FETCHCONTENT_FALLBACK=ON or provide pybind11.")
endif()

# ---------- Xcoll BDSIM ----------
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
target_include_directories(g4interface PUBLIC ${BDSIM_INCLUDE_DIR})
