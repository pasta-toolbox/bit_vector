/*******************************************************************************
 * bit_vector/flattened_rank_support.hpp
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

#pragma once

#include <immintrin.h>
#include <emmintrin.h>

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/popcount.hpp"

namespace pasta {

  struct L12Type {
    uint64_t l1 : 40;
    std::array<uint16_t, 7> l2;
  }; // struct L12Type

  /*
   * We store the L1-information of 8 L2-blocks in 40 bits and the
   * corresponding seven L2-information in 12 bits each. Using 12 bits
   * allows us to store the prefix sum of the popcounts in the 512-bit
   * L2-blocks. All this can be stored in 128 bits. To be precise 124 bits
   * would suffice, but the additional 4 bits allow for faster access and
   * result in exactly the same overhead the non-flat popcount rank and
   * select data structure has.
   *
   * +------+------+------+------+------+------+---+------+------+--------+
   * | 8bit | 8bit | 8bit | 8bit | 8bit | 8bit |...| 8bit | 8bit | 40 bit |
   * +------+------+------+------+------+------+---+------+------+--------+
   * | 3 * 12-bit integer | 3 * 12 bit integer |...| 12 bit int. | L1-info|
   * +------+------+------+------+------+------+---+------+------+--------+
   * |   8  | 4/4  |   8  |   8  | 4/4  |  8   |...|  8   | 4/0  | 40 Bit |
   * +------+------+------+------+------+------+---+------+------+--------+
   *
   * The order of the 12-bit integers should remain the same (as they occur
   * in the bit vector). This helps us to determine the correct block later
   * on. To this end, we have to split
   */
  struct L12TypeI {
    __uint128_t data;

    L12TypeI() = default;

    L12TypeI(uint64_t const _l1, std::array<uint16_t, 7>& _l2)
      : data(((__uint128_t{0b111111111111} & _l2[6]) << 116) |
	     ((__uint128_t{0b111111111111} & _l2[5]) << 104) |
	     ((__uint128_t{0b111111111111} & _l2[4]) << 92) |
	     ((__uint128_t{0b111111111111} & _l2[3]) << 80) |
	     ((__uint128_t{0b111111111111} & _l2[2]) << 68) |
	     ((__uint128_t{0b111111111111} & _l2[1]) << 56) |
	     ((__uint128_t{0b111111111111} & _l2[0]) << 44) |
	     ((__uint128_t{0xFFFFFFFFFF} & _l1)))  { }
    inline uint16_t operator[](size_t const index) const {
      return (data >> ((12 * index) + 44) &  uint16_t(0b111111111111));
    }

    inline uint64_t l1() const {
      return uint64_t{0xFFFFFFFFFF} & data;
    }

  };

  static_assert(sizeof(L12TypeI) == 16);

  struct FlattenedRankSelectConfig {
    static constexpr size_t L2_BIT_SIZE = 512;
    static constexpr size_t L1_BIT_SIZE = 8 * L2_BIT_SIZE;

    static constexpr size_t L1_WORD_SIZE = L1_BIT_SIZE / (sizeof(uint64_t) * 8);
    static constexpr size_t L2_WORD_SIZE = L2_BIT_SIZE / (sizeof(uint64_t) * 8);
  }; // struct FlattenedRankSelectConfig

  class FlattenedBitVectorRank{

    size_t data_size_;
    uint64_t const* data_;

    tlx::SimpleVector<L12TypeI, tlx::SimpleVectorMode::NoInitNoDestroy> l12_;

  public:
    FlattenedBitVectorRank() = default;

    FlattenedBitVectorRank(BitVector const& bv)
      : data_size_(bv.size_),
	data_(bv.data_.data()),
	l12_((data_size_ / FlattenedRankSelectConfig::L1_WORD_SIZE) + 1) {
      init();
    }

    size_t rank0(size_t const index) const {
      return index - rank1(index);
    }

    size_t rank1(size_t index) const {
      size_t const l1_pos = index / FlattenedRankSelectConfig::L1_BIT_SIZE;
      __builtin_prefetch(&l12_[l1_pos], 0, 0);
      int32_t l2_pos = ((index % FlattenedRankSelectConfig::L1_BIT_SIZE) /
			FlattenedRankSelectConfig::L2_BIT_SIZE);
      size_t offset = (l1_pos * FlattenedRankSelectConfig::L1_WORD_SIZE) +
	(l2_pos * FlattenedRankSelectConfig::L2_WORD_SIZE);
      __builtin_prefetch(&data_[offset], 0, 0);
      --l2_pos;

      size_t result = l12_[l1_pos].l1() +
	((l2_pos >= 0) ? l12_[l1_pos][l2_pos] : 0);
      index %= FlattenedRankSelectConfig::L2_BIT_SIZE;
      for (size_t i = 0; i < index / 64; ++i) {
	result += std::popcount(data_[offset++]);
      }
      if (index %= 64; index > 0) [[likely]] {
	uint64_t const remaining = data_[offset] << (64 - index);
	result += std::popcount(remaining);
      }
      return result;
    }

  private:

    void init() {
      size_t l12_pos = 0;
      uint64_t l1_entry = 0ULL;

      uint64_t const * data = data_;
      uint64_t const * const data_end = data_ + data_size_;

      std::array<uint16_t, 7> l2_entries = {0, 0, 0, 0, 0, 0, 0};
      while (data + 64 <= data_end) {
	l2_entries[0] = popcount<8>(data);
	data += 8;
	for (size_t i = 1; i < 7; ++i) {
	  l2_entries[i] = l2_entries[i - 1] + popcount<8>(data);
	  data += 8;
	}
	l12_[l12_pos++] = L12TypeI(l1_entry, l2_entries);
	l1_entry += l2_entries.back() + popcount<8>(data);
	data += 8;
      }
      size_t l2_pos = 0;
      while (data + 8 <= data_end) {
	l2_entries[l2_pos++] = popcount<8>(data);
	data += 8;
      }
      while (data < data_end) {
	l2_entries[l2_pos] += popcount<1>(data++);
      }
      l12_[l12_pos++] = L12TypeI(l1_entry, l2_entries);
    }
  }; // class FlattenedBitVectorRank

} // namespace pasta

/******************************************************************************/
