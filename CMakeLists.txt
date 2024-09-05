# ##############################################################################
# CMakeLists.txt
#
# Copyright (C) 2021-2024 Florian Kurpicz <florian@kurpicz.org>
#
# pasta::bit_vector is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# pasta::bit_vector is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# pasta::bit_vector.  If not, see <http://www.gnu.org/licenses/>.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.25)

project(pasta_bit_vector
  DESCRIPTION "Bit Vector with Compact and Fast Rank and Select Support"
  LANGUAGES CXX
  VERSION 0.1.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

if (PROJECT_IS_TOP_LEVEL)
  # Enable easy formatting
  FetchContent_Declare(
    Format.cmake
    GIT_REPOSITORY https://github.com/TheLartians/Format.cmake
    GIT_TAG v1.8.1
  )
  FetchContent_MakeAvailable(Format.cmake)

  add_subdirectory(docs)
endif ()

# Options when compiling pasta::bit_vector Build tests
option(PASTA_BIT_VECTOR_BUILD_TESTS "Build pasta::bit_vector's tests." OFF)
# Build benchmark tools
option(PASTA_BIT_VECTOR_BUILD_BENCHMARKS
       "Build pasta::bit_vector's benchmarks." OFF
)
# Build pasta::bit_vector with code coverage options
option(PASTA_BIT_VECTOR_COVERAGE_REPORTING
       "Enable coverage reporting for pasta::bit_vector" OFF
)

# Optional code coverage (library compile options are only set if coverage
# reporting is enabled
add_library(pasta_bit_vector_coverage_config INTERFACE)
if (PASTA_BIT_VECTOR_COVERAGE_REPORTING)
  target_compile_options(
    pasta_bit_vector_coverage_config INTERFACE -fprofile-arcs -ftest-coverage
  )
  target_link_libraries(pasta_bit_vector_coverage_config INTERFACE gcov)
endif ()

# pasta::bit_vector interface definitions
add_library(pasta_bit_vector INTERFACE)
target_include_directories(
  pasta_bit_vector INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
  pasta_bit_vector INTERFACE pasta_utils pasta_bit_vector_coverage_config
)

# Use FetchContent to load dependencies
FetchContent_Declare(
  pasta_utils
  GIT_REPOSITORY https://github.com/pasta-toolbox/utils.git
  GIT_TAG origin/main
)
FetchContent_MakeAvailable(pasta_utils)

# Optional test
if (PASTA_BIT_VECTOR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif ()

# Optional benchmarks
if (PASTA_BIT_VECTOR_BUILD_BENCHMARKS)
  FetchContent_Declare(
    tlx
    GIT_REPOSITORY https://github.com/tlx/tlx.git
    GIT_TAG origin/main
  )
  FetchContent_MakeAvailable(tlx)

  add_executable(bit_vector_benchmark benchmarks/bit_vector_benchmark.cpp)

  target_link_libraries(
    bit_vector_benchmark
    PUBLIC pasta_bit_vector
           tlx
           pasta_utils
           optimized
           pasta_memory_monitor
  )
endif ()

# ##############################################################################
