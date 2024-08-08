cmake_minimum_required(VERSION 3.5.0)

# set(CMAKE_BUILD_TYPE DEBUG)

file(READ "${CMAKE_SOURCE_DIR}/VERSION.txt" PROJECT_VERSION)
string(STRIP "${PROJECT_VERSION}" PROJECT_VERSION)
message(STATUS "Project version: ${PROJECT_VERSION}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message("cmake build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(segy VERSION ${PROJECT_VERSION})

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${segy_SOURCE_DIR}/Dist CACHE PATH "" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message("Install Prefix: " ${CMAKE_INSTALL_PREFIX})

set(BUILD_PYTHON ON)
set(BUILD_TOOLS OFF)
set(ENABLE_OPENMP OFF)
set(BUILD_TEST OFF)

if (BUILD_PYTHON) 
  if (MSVC)
    # TODO: set this
    set(PYTHON_EXECUTABLE "D:\\Program\\Anaconda\\envs\\py38\\python.exe")
    set(PYTHON_LIBRARIES "D:\\Program\\Anaconda\\envs\\py38\\Lib")
    set(pybind11_ROOT "D:\\Program\\Anaconda\\Lib\\site-packages\\pybind11")
  endif()
  find_package(pybind11 REQUIRED)
endif()

if (APPLE)
  set(CMAKE_BUILD_RPATH "@executable_path")
endif()

# openmp
if (ENABLE_OPENMP)
  if (APPLE)
    set(OpenMP_ROOT /opt/homebrew/opt/libomp/)
  endif()
  find_package(OpenMP)
endif()
if (OpenMP_CXX_FOUND)
  message("Use OpenMP")
  set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
endif()

# fmt
if (APPLE)
  set(fmt_ROOT /opt/homebrew/Cellar/fmt/9.1.0/)  
endif()
if (MSVC)
  set(fmt_ROOT D:\\fmt\\out\\install\\x64-Release)
endif()
find_package(fmt REQUIRED)

if (MSVC)
add_compile_options(/wd4244 /wd4996 /wd4819 /Ox)
endif()

add_subdirectory(src)

if (BUILD_TOOLS)
  add_subdirectory(tools)
endif()

if (BUILD_PYTHON)
  add_subdirectory(python)
endif()

if (BUILD_TEST)
  add_subdirectory(test)
endif()




