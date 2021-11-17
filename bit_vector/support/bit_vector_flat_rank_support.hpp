/*******************************************************************************
 * bit_vector/bit_vector_flat_rank_support.hpp
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

#include <numeric>

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/l12_type.hpp"
#include "bit_vector/support/popcount.hpp"

namespace pasta {

  struct FlattenedRankSelectConfig {
    static constexpr size_t L2_BIT_SIZE = 512;
    static constexpr size_t L1_BIT_SIZE = 8 * L2_BIT_SIZE;

    static constexpr size_t L1_WORD_SIZE = L1_BIT_SIZE / (sizeof(uint64_t) * 8);
    static constexpr size_t L2_WORD_SIZE = L2_BIT_SIZE / (sizeof(uint64_t) * 8);
  }; // struct FlattenedRankSelectConfig

  class BitVectorFlatRank{

    size_t data_size_;
    uint64_t const* data_;

    tlx::SimpleVector<BigL12Type, tlx::SimpleVectorMode::NoInitNoDestroy> l12_;

  public:
    //! Default constructor w/o parameter.
    BitVectorFlatRank() = default;

    /*!
     * \brief Constructor. Creates the auxiliary information for efficient rank
     * queries.
     * \param bv \c BitVector the rank structure is created for.
     */
    BitVectorFlatRank(BitVector const& bv)
      : data_size_(bv.size_),
	data_(bv.data_.data()),
	l12_((data_size_ / FlattenedRankSelectConfig::L1_WORD_SIZE) + 1) {
      init();
    }

    /*!
     * \brief Computes rank of zeros.
     * \param index Index the rank of zeros is computed for.
     * \return Number of zeros (rank) before position \c index.
     */
    [[nodiscard("rank0 computed but not used")]]
    size_t rank0(size_t const index) const {
      return index - rank1(index);
    }

    /*!
     * \brief Computes rank of ones.
     * \param index Index the rank of ones is computed for.
     * \return Number of ones (rank) before position \c index.
     */
    [[nodiscard("rank1 computed but not used")]]
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

    /*!
     * \brief Estimate for the space usage.
     * \return Number of bytes used by this data structure.
     */
    [[nodiscard("space useage computed but not used")]]
    size_t space_useage() const {
      return l12_.size() * sizeof(BigL12Type) + sizeof(*this);
    }

  private:

    //! Function used for initializing data structure to reduce LOCs of
    //! constructor.
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
	l12_[l12_pos++] = BigL12Type(l1_entry, l2_entries);
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
      std::partial_sum(l2_entries.begin(), l2_entries.end(),
		       l2_entries.begin());
      l12_[l12_pos++] = BigL12Type(l1_entry, l2_entries);
    }
  }; // class BitVectorFlatRank

} // namespace pasta

/******************************************************************************/
