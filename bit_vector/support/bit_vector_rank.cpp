/*******************************************************************************
 * pasta/container/support/bit_vector_rank.cpp
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

#include "container/bit_vector.hpp"
#include "container/support/bit_vector_rank.hpp"
#include "container/support/popcount.hpp"

namespace pasta {

  BitVectorRank::BitVectorRank(BitVector const& bv)
    : data_size_(bv.size_),
      data_(bv.data_.data()),
      l0_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1),
      l12_((data_size_ / PopcntRankSelectConfig::L1_WORD_SIZE) + 1) {
    init();
  }

  BitVectorRank::~BitVectorRank() { }

  void BitVectorRank::init() {
    l0_[0] = 0;

    size_t l0_pos = 1;
    size_t l12_pos = 0;
    uint32_t l1_entry = 0UL;

    uint64_t const * data = data_;
    uint64_t const * const data_end = data_ + data_size_;

    // For each full L12-Block
    std::array<uint16_t, 3> l2_entries = { 0, 0, 0 };
    while (data + 32 <= data_end) {
      uint32_t new_l1_entry = l1_entry;
      for (size_t i = 0; i < 3; ++i) {
        l2_entries[i] = popcount<8>(data);
        data += 8;
        new_l1_entry += l2_entries[i];
      }
      l12_[l12_pos++] = L12Entry(l1_entry, l2_entries);
      new_l1_entry += popcount<8>(data);
      data += 8;
      l1_entry = new_l1_entry;

      if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
		     PopcntRankSelectConfig::L1_WORD_SIZE) == 0) [[unlikely]] {
        l0_[l0_pos] = (l0_[l0_pos - 1] + l1_entry);
        ++l0_pos;
        l1_entry = 0;
      }
    }

    // For the last not full L12-Block
    size_t l2_pos = 0;
    while (data + 8 <= data_end) {
      l2_entries[l2_pos++] = popcount<8>(data);
      data += 8;
    }
    while (data < data_end) {
      l2_entries[l2_pos] += popcount<1>(data++);
    }
    l12_[l12_pos++] = L12Entry(l1_entry, l2_entries);

    if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
		   PopcntRankSelectConfig::L1_WORD_SIZE) == 0) [[unlikely]] {
      l0_[l0_pos] += (l0_[l0_pos - 1] + l1_entry);
      ++l0_pos;
      l1_entry = 0;
    }
  }

  size_t BitVectorRank::rank0(size_t index) const {
    return index - rank1(index);
  }

  size_t BitVectorRank::rank1(size_t index) const {
    size_t const l1_pos = index / PopcntRankSelectConfig::L1_BIT_SIZE;
    __builtin_prefetch(&l12_[l1_pos], 0, 0);
    size_t const l2_pos = (index % PopcntRankSelectConfig::L1_BIT_SIZE) /
      PopcntRankSelectConfig::L2_BIT_SIZE;
    size_t offset = (l1_pos * PopcntRankSelectConfig::L1_WORD_SIZE) +
      (l2_pos * PopcntRankSelectConfig::L2_WORD_SIZE);
    __builtin_prefetch(&data_[offset], 0, 0);

    size_t result = l0_[index / PopcntRankSelectConfig::L0_BIT_SIZE] +
      l12_[l1_pos].l1;

    for (size_t i = 0; i < l2_pos; ++i) {
      result += l12_[l1_pos][i];
    }

    index %= PopcntRankSelectConfig::L2_BIT_SIZE;
    for (size_t i = 0; i < index / 64; ++i) {
      result += std::popcount(data_[offset++]);
    }
    if (index %= 64; index > 0) [[likely]] {
      uint64_t const remaining = data_[offset] << (64 - index);
      result += std::popcount(remaining);
    }
    return result;
  }

  size_t BitVectorRank::space_usage() const {
    return l0_.size() * sizeof(uint64_t)
        + l12_.size() * sizeof(L12Entry)
        + sizeof(*this);
  }
} // pasta

/******************************************************************************/
