/*******************************************************************************
 * tests/bit_vector/support/bit_vector_flattened_rank_test.cpp
 *
 * Copyright (C) 2021 Florian Kurpicz <florian@kurpicz.org>
 *
 * pasta::bit_vector is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pasta::bit_vector is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::bit_vector.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <cstdint>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank.hpp>
#include <tlx/die.hpp>
#include <vector>

template <typename TestFunction>
void run_test(TestFunction test_config) {
  std::vector<size_t> offsets = {0, 7, 723, 1347};
  for (size_t n = 2; n <= 32; n += 10) {
    for (auto const offset : offsets) {
      size_t const vector_size = (1ULL << n) + offset;
      if (n == 32) {
        test_config(vector_size, 1ULL << 2);
        continue;
      }
      for (size_t k = 0; k <= 4; ++k) {
        size_t const set_every_kth = 1 + (1ULL << k);
        if (k < n) { // if k > n this testing doesn't make any sense
          test_config(vector_size, set_every_kth);
        }
      }
    }
  }
}

int32_t main() {
  run_test([](size_t N, size_t K) {
    pasta::BitVector bv(N, 0);
    size_t set_ones = 0;

    auto bv_data = bv.data();
    for (size_t i = 0; i < bv_data.size(); ++i) {
      uint64_t word = 0ULL;
      for (size_t j = 0; j < 64; ++j) {
        word >>= 1;
        if (size_t bit_pos = ((i * 64) + j); bit_pos >= N) {
          word >>= (63 - j);
          break;
        } else if (bit_pos % K == 0) {
          ++set_ones;
          word |= (1ULL << 63);
        }
      }
      bv_data.data()[i] = word;
    }

    size_t const query_pos_offset = (N > (1ULL << 30)) ? 101 : 1;

    // Test optimized for one queries
    {
      pasta::FlatRank<pasta::OptimizedFor::ONE_QUERIES> bvr(bv);

      die_unequal(set_ones, bvr.rank1(N));
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(i, bvr.rank1((K * i)));
      }

      die_unequal((N - set_ones), bvr.rank0(N));
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal((K - 1) * i, bvr.rank0((K * i)));
      }
    }
    // Test optimized for zero queries
    {
      pasta::FlatRank<pasta::OptimizedFor::ZERO_QUERIES> bvr(bv);

      die_unequal(set_ones, bvr.rank1(N));
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(i, bvr.rank1((K * i)));
      }

      die_unequal((N - set_ones), bvr.rank0(N));
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal((K - 1) * i, bvr.rank0((K * i)));
      }
    }
  });

  return 0;
}

/******************************************************************************/
