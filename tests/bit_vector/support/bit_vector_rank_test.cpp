/*******************************************************************************
 * tests/container/support/bit_vector_rank_test.cpp
 *
 * Copyright (C) 2019-2021 Florian Kurpicz <florian@kurpicz.org>
 *
 * PaStA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PaStA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PaStA.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <cstdint>

#include <tlx/die.hpp>

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/bit_vector_rank.hpp"

template<typename TestFunction>
void run_test(TestFunction test_config) {
  std::vector<size_t> offsets = { 0, 7, 723, 1347 };
  for (size_t n = 2; n <= 32; n += 10) {
    for (auto const offset : offsets) {
      size_t const vector_size = (1ULL << n) + offset;
      if (n == 32) {
        test_config(vector_size, 1ULL << 2);
        continue;
      }
      for (size_t k = 0; k <= 4; ++k) {
        size_t const set_every_kth = 1ULL << k;
        if (k < n) { // if k > n this testing doesn't make any sense
          test_config(vector_size, set_every_kth);
        }
      }
    }
  }
}

int32_t main() {
  run_test([](size_t N, size_t K){
    pasta::BitVector bv(N, 0);
    size_t set_ones = 0;
    for (size_t i = 0; i < N; i += K) {
      ++set_ones;
      bv[i] = 1;
    }
    // Test optimized for one queries
    {
      pasta::BitVectorRank<pasta::OptimizedFor::ONE_QUERIES> bvr(bv);

      die_unequal(set_ones, bvr.rank1(N));
      for (size_t i = 1; i <= N/K; ++i) {
        die_unequal(i, bvr.rank1((K * i)));
      }

      die_unequal((N - set_ones), bvr.rank0(N));
      for (size_t i = 1; i <= N/K; ++i) {
        die_unequal((K - 1) * i, bvr.rank0((K * i)));
      }
    }
    // Test optimized for zero queries
    {
      pasta::BitVectorRank<pasta::OptimizedFor::ZERO_QUERIES> bvr(bv);

      die_unequal(set_ones, bvr.rank1(N));
      for (size_t i = 1; i <= N/K; ++i) {
        die_unequal(i, bvr.rank1((K * i)));
      }

      die_unequal((N - set_ones), bvr.rank0(N));
      for (size_t i = 1; i <= N/K; ++i) {
        die_unequal((K - 1) * i, bvr.rank0((K * i)));
      }
    }
  });

  return 0;
}

/******************************************************************************/
