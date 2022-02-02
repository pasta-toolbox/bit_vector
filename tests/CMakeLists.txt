################################################################################
# tests/CMakeLists.txt
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

macro(pasta_build_test TESTNAME)
  string(REPLACE "/" "_" TESTNAME_REPLACED "${TESTNAME}")
  add_executable(${TESTNAME_REPLACED}
    ${TESTNAME}.cpp)
  target_link_libraries(${TESTNAME_REPLACED} PRIVATE
    pasta_bit_vector
    tlx)
  add_test(
    NAME ${TESTNAME_REPLACED}
    COMMAND ${TESTNAME_REPLACED} ${ARGN})
endmacro(pasta_build_test)

pasta_build_test(bit_vector/bit_vector_test)
pasta_build_test(bit_vector/support/bit_vector_rank_test)
pasta_build_test(bit_vector/support/bit_vector_flat_rank_test)
pasta_build_test(bit_vector/support/bit_vector_rank_select_test)
pasta_build_test(bit_vector/support/bit_vector_flat_rank_select_test)
pasta_build_test(bit_vector/support/bit_vector_wide_rank_test)
pasta_build_test(bit_vector/support/bit_vector_wide_rank_select_test)

################################################################################
