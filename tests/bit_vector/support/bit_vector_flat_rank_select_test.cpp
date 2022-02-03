/*******************************************************************************
 * tests/container/support/bit_vector_flatrank_select_test.cpp
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
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <tlx/die.hpp>

template <typename TestFunction>
void run_test(TestFunction test_config) {
  std::vector<size_t> offsets = {0, 723};
  std::vector<size_t> bit_sizes = {1ULL << 2,
                                   1ULL << 12,
                                   1ULL << 32,
                                   (1ULL << 32) + (1ULL << 12)};
  // for (size_t n = 2; n <= 32; n += 10) {
  for (auto const& bit_size : bit_sizes) {
    for (auto const offset : offsets) {
      size_t const vector_size = bit_size + offset;
      if (bit_size >= (1ULL << 32)) {
        test_config(vector_size, 3);
        continue;
      }
      for (size_t k = 0; k <= 4; ++k) {
        size_t const set_every_kth = 1ULL << k;
        // if k > bit_size this testing doesn't make any sense
        if (k < bit_size) {
          test_config(vector_size, set_every_kth);
        }
      }
    }
  }
}

int32_t main() {
  // Test select
  run_test([](size_t N, size_t K) {
    pasta::BitVector bv(N, 0);
    auto bv_data = bv.data();
    for (size_t i = 0; i < bv_data.size(); ++i) {
      uint64_t word = 0ULL;
      for (size_t j = 0; j < 64; ++j) {
        word >>= 1;
        if (size_t bit_pos = ((i * 64) + j); bit_pos >= N) {
          word >>= (63 - j);
          break;
        } else if (bit_pos % K == 0) {
          word |= (1ULL << 63);
        }
      }
      bv_data.data()[i] = word;
    }

    size_t const query_pos_offset = (N > (1ULL << 30)) ? 101 : 1;

    // Test optimized for one queries with linear search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::LINEAR_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }
    // Test optimized for zero queries with linear search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::LINEAR_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }

    // Test optimized for one queries with binary search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }
    // Test optimized for zero queries with binary search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }

    // Test optimized for one queries with intrinsics
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::INTRINSICS>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }
    // Test optimized for zero queries with intrinsics
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::INTRINSICS>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += query_pos_offset) {
        die_unequal(K * (i - 1), bvrs.select1(i));
      }
    }

    for (size_t i = 0; i < bv_data.size(); ++i) {
      bv_data.data()[i] = ~bv_data.data()[i];
    }

    // Test optimized for one queries with linear search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::LINEAR_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }
    // Test optimized for zero queries with linear search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::LINEAR_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }

    // Test optimized for one queries with binary search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }
    // Test optimized for zero queries with binary search
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }

    // Test optimized for one queries with intrinsics
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::INTRINSICS>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }
    // Test optimized for zero queries with intrinsics
    {
      pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                            pasta::FindL2FlatWith::INTRINSICS>
          bvrs(bv);
      for (size_t i = 1; i <= N / K; i += (std::max<size_t>(1, N / 100) + 1)) {
        die_unequal(K * (i - 1), bvrs.select0(i));
      }
    }
  });

  return 0;
}

/******************************************************************************/
