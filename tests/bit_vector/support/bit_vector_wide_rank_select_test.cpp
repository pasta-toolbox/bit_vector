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
#include <pasta/bit_vector/support/find_l2_wide_with.hpp>
#include <pasta/bit_vector/support/wide_rank_select.hpp>
#include <tlx/die.hpp>
#include <vector>

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
    {
      pasta::BitVector bv(N, 0);
      for (size_t i = 0; i < N; i += K) {
        bv[i] = 1;
      }

      // Test optimized for one queries without intrinsics
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                              pasta::FindL2WideWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                              pasta::FindL2WideWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                              pasta::FindL2WideWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select1(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                              pasta::FindL2WideWith::BINARY_SEARCH>
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

      // Test optimized for one queries without intrinsics
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                              pasta::FindL2WideWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                              pasta::FindL2WideWith::BINARY_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                              pasta::FindL2WideWith::LINEAR_SEARCH>
            bvrs(bv);
        for (size_t i = 1; i <= N / K;
             i += (std::max<size_t>(1, N / 100) + 1)) {
          die_unequal(K * (i - 1), bvrs.select0(i));
        }
      }
      {
        pasta::WideRankSelect<pasta::OptimizedFor::ZERO_QUERIES,
                              pasta::FindL2WideWith::BINARY_SEARCH>
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
