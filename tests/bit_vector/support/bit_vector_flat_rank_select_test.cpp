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

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/bit_vector_flat_rank_select.hpp"
#include "bit_vector/support/find_l2_flat_with.hpp"

#include <cstdint>
#include <tlx/die.hpp>

template <typename TestFunction>
void run_test(TestFunction test_config) {
  std::vector<size_t> offsets = {0, 723};
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
  // Test rank
  run_test([](size_t N, size_t K) {
    pasta::BitVector bv(N, 0);

    size_t set_ones = 0;
    for (size_t i = 0; i < N; i += K) {
      ++set_ones;
      bv[i] = 1;
    }

    // We test both configurations for rank in the rank test
    pasta::BitVectorFlatRankSelect bvrs(bv);

    die_unequal(set_ones, bvrs.rank1(N));
    for (size_t i = 1; i <= N / K; i += std::max<size_t>(1, N / 1000)) {
      die_unequal(i, bvrs.rank1((K * i)));
    }

    die_unequal((N - set_ones), bvrs.rank0(N));
    for (size_t i = 1; i <= N / K; i += std::max<size_t>(1, N / 1000)) {
      die_unequal((K - 1) * i, bvrs.rank0((K * i)));
    }
  });

  // Test select
  run_test([](size_t N, size_t K) {
    {
      pasta::BitVector bv(N, 0);
      for (size_t i = 0; i < N; i += K) {
        bv[i] = 1;
      }
      // Test optimized for one queries with linear search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      // Test optimized for zero queries with linear search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }

      // Test optimized for one queries with binary search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      // Test optimized for zero queries with binary search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }

      // Test optimized for one queries with intrinsics
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::INTRINSICS>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      // Test optimized for zero queries with intrinsics
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::INTRINSICS>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
    }
    {
      pasta::BitVector bv(N, 1);
      for (size_t i = 0; i < N; i += K) {
        bv[i] = 0;
      }

      // Test optimized for one queries with linear search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      // Test optimized for zero queries with linear search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }

      // Test optimized for one queries with binary search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      // Test optimized for zero queries with binary search
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }

      // Test optimized for one queries with intrinsics
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                                       pasta::FindL2FlatWith::INTRINSICS>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      // Test optimized for zero queries with intrinsics
      {
        pasta::BitVectorFlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                                       pasta::FindL2FlatWith::INTRINSICS>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
    }
  });

  return 0;
}

/******************************************************************************/
