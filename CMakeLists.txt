################################################################################
# CMakeLists.txt
#
# Copyright (C) 2021 Florian Kurpicz <florian@kurpicz.org>
#
# pasta::bit_vector is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pasta::bit_vector is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pasta::bit_vector.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(pasta_bit_vector)

# Options when compiling pasta::bit_vector
option(PASTA_BIT_VECTOR_BUILD_TESTS
  "Build pasta::bit_vector's tests." OFF)

# Add tlx subdirectory only, if it has not been added by another project already
if(NOT TARGET tlx)
  add_subdirectory(extlib/tlx)
endif()

# 
add_library(pasta_bit_vector INTERFACE)
target_include_directories(pasta_bit_vector INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${TLX_INCLUDE_DIRS})

# Optional test
if(PASTA_BIT_VECTOR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

################################################################################
