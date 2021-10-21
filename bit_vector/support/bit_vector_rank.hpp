/*******************************************************************************
 * pasta/container/support/bit_vector_rank.hpp
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

/*
 * Based on
 *
 * @inproceedings{ZhouAK2013PopcountRankSelect,
 *    author    = {Dong Zhou and David G. Andersen and Michael Kaminsky},
 *    title     = {Space-Efficient, High-Performance Rank and Select Structures
 *                 on Uncompressed Bit Sequences},
 *    booktitle = {12th International Symposium on Experimental Algorithms
 *                 ({SEA})},
 *    series    = {LNCS},
 *    volume    = {7933},
 *    pages     = {151--163},
 *    publisher = {Springer},
 *    year      = {2013},
 *    doi       = {10.1007/978-3-642-38527-8\_15},
 * }
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

#include <tlx/container/simple_vector.hpp>

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/l12_type.hpp"
#include "bit_vector/support/popcount.hpp"

namespace pasta {

  /*!
   * \brief Static configuration for \c BitVectorRank and
   * \c BitVectorRankSelect.
   */
  struct PopcntRankSelectConfig {
    //! Bits covered by an L2-block.
    static constexpr size_t L2_BIT_SIZE = 512;
    //! Bits covered by an L1-block.
    static constexpr size_t L1_BIT_SIZE = 4 * L2_BIT_SIZE;
    //! Bits covered by an L0-block.
    static constexpr size_t L0_BIT_SIZE =
      static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) + 1;

    //! Number of 64-bit words covered by an L2-block.
    static constexpr size_t L2_WORD_SIZE = L2_BIT_SIZE / (sizeof(uint64_t) * 8);
    //! Number of 64-bit words covered by an L1-block.
    static constexpr size_t L1_WORD_SIZE = L1_BIT_SIZE / (sizeof(uint64_t) * 8);
    //! Number of 64-bit words covered by an L0-block.
    static constexpr size_t L0_WORD_SIZE = L0_BIT_SIZE / (sizeof(uint64_t) * 8);
    //! Sample rate of positions for faster select queries.
    static constexpr size_t SELECT_SAMPLE_RATE = 8192;
  }; // struct PopcountRankSelectConfiguration


  //! \addtogroup pasta_bit_vectors
  //! \{

  /*!
   * \brief Rank support for \c BitVector.
   *
   * The rank support is based on popcount and described in detail by
   * Zhou et al. \cite ZhouAK2013PopcountRankSelect. The data structure consists
   * of three levels (L0, L1, and L2) that contain different information
   * regarding the popcount of the current block or all previous blocks.
   *
   * *Note* that rank support is also provided in addition to select support by
   * \c BitVectorRankSelect, which uses this rank support implementation
   * internally.
   */
  class BitVectorRank {

    //! Friend class, using internal information l0_ and l12_, too.
    friend class BitVectorRankSelect;

    //! Size of the bit vector the rank support is constructed for.
    size_t data_size_;
    //! Pointer to the data of the bit vector.
    uint64_t const * data_;

    //! Array containing the number of set bits in the L0-blocks.
    tlx::SimpleVector<uint64_t, tlx::SimpleVectorMode::NoInitNoDestroy> l0_;

    //! Array containing the information about the L1- and L2-blocks.
    tlx::SimpleVector<L12Entry, tlx::SimpleVectorMode::NoInitNoDestroy> l12_;

  public:
    //! Default constructor w/o parameter.
    BitVectorRank() = default;

    /*!
     * \brief Constructor. Creates the auxiliary information for efficient rank
     * queries.
     * \param bv \c BitVector the rank structure is created for.
     */
    BitVectorRank(BitVector const& bv)
      : data_size_(bv.size_),
	data_(bv.data_.data()),
	l0_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1),
	l12_((data_size_ / PopcntRankSelectConfig::L1_WORD_SIZE) + 1) {
      init();
    }

    // Default move constructor.
    BitVectorRank(BitVectorRank&& other) = default;

    // Default move assignment.
    BitVectorRank& operator = (BitVectorRank&& other) = default;

    //BitVectorRank(BitVector::Iterator begin, BitVector::Iterator end);

    //! Destructor. Deleting manually created arrays.
    ~BitVectorRank() = default;

    /*!
     * \brief Computes rank of zeros.
     * \param index Index the rank of zeros is computed for.
     * \return Numbers of zeros (rank) before position \c index.
     */
    [[nodiscard("rank0 computed but not used")]]
    size_t rank0(size_t const index) const {
      return index - rank1(index);
    }

    /*!
     * \brief Computes rank of ones.
     * \param index Index the rank of ones is computed for.
     * \return Numbers of ones (rank) before position \c index.
     */
    [[nodiscard("rank1 computed but not used")]]
    size_t rank1(size_t index) const {
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

    /*!
     * \brief Estimate for the space usage.
     * \return Number of bytes used by this data structure.
     */
    [[nodiscard("space usage computed but not used")]]
    size_t space_usage() const {
      return l0_.size() * sizeof(uint64_t)
        + l12_.size() * sizeof(L12Entry)
        + sizeof(*this);
    }

  private:
    //! Function used initializing data structure to reduce LOCs of constructor.
    void init() {
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
  }; // class BitVectorRank

  //! \}
  
} // namespace pasta

/******************************************************************************/
