###################################################################################
#
# 1D flow solver
#
# Copyright (c) 2023 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

cmake_minimum_required(VERSION 3.16)

option(enable_OPENMP "Enable OpenMP" ON)
option(enable_GLOG "Enable GLOG" OFF)

#-----compiler options
if(NOT DEFINED compiler)
  set(base_compile_options "-Wall -Wextra -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -fopenmp")
   endif()
endif()

if(compiler STREQUAL intel)
  set(CMAKE_C_COMPILER "icx")
  set(CMAKE_CXX_COMPILER "icpx")
  set(base_compile_options "-Wall")
  if(enable_GLOG)
    set(base_compile_options ${base_compile_options} " -D_GLOG")
  endif()
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -qopenmp")
  endif()
elseif(DEFINED compiler)
  message(FATAL ERROR "invalid compiler")
endif()
string(REPLACE " " ";" compile_options ${base_compile_options})
#-----------------------

project(Taylor_Galerkin C CXX)
set(CMAKE_CXX_STANDARD 17)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

#change default install prefix
set(INSTALL_DIR "/Users/t-enomoto/research/1d_blood/src/Flow1D")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "..." FORCE)
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

#######
# Display options
#######
message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

#######
# find Libraries to depend on
#######

# Google-glog
if(enable_GLOG)
  set(GLOG_DIR ${GLOG_DIR} CACHE PATH "...")
  find_package(GLOG REQUIRED)
  IF(GLOG_FOUND)
  message( STATUS "GLOG include dir  : "      ${GLOG_INCLUDE_DIR})
  MESSAGE("find GLOG")
  include_directories(${GLOG_INCLUDE_DIR})
  ENDIF(GLOG_FOUND)
endif()

# Eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

#######
# Sub Directories
#######

# add_subdirectory(include)
add_subdirectory(src)