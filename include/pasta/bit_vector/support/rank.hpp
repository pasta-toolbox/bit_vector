/*******************************************************************************
 * This file is part of pasta::bit_vector.
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

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/l12_type.hpp"
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/bit_vector/support/popcount.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <pasta/utils/container/aligned_vector.hpp>
#include <pasta/utils/debug_asserts.hpp>
#include <tlx/container/simple_vector.hpp>

namespace pasta {

/*!
 * \ingroup pasta_bit_vector_configuration
 * \brief Static configuration for \c Rank and
 * \c RankSelect.
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

//! \addtogroup pasta_bit_vector_rank
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
 * \c RankSelect, which uses this rank support implementation
 * internally.
 *
 * \tparam OptimizedFor Compile time option to optimize data structure for
 * either 0, 1, or no specific type of query.
 * \tparam VectorType Type of the vector the rank data structure is constructed
 * for, e.g., plain \c BitVector or a compressed bit vector.
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE,
          typename VectorType = BitVector>
class Rank {
public:
  //! Size of the bit vector the rank support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  VectorType::RawDataConstAccess data_;
  //! Size of the bit vector in bits (only used for debug asserts)
  size_t const bit_size_;

  //! Array containing the number of set bits in the L0-blocks.
  tlx::SimpleVector<uint64_t, tlx::SimpleVectorMode::NoInitNoDestroy> l0_;

  //! Array containing the information about the L1- and L2-blocks.
  tlx::SimpleVector<L12Type, tlx::SimpleVectorMode::NoInitNoDestroy> l12_;

public:
  //! Default constructor w/o parameter.
  Rank() = default;

  /*!
   * \brief Constructor. Creates the auxiliary information for efficient rank
   * queries.
   * \param bv \c Vector of type \c VectorType the rank structure is created
   * for.
   */
  Rank(VectorType& bv)
      : data_size_(bv.size_),
        data_(bv.data_.data()),
        bit_size_(bv.size()),
        l0_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 2),
        l12_((data_size_ / PopcntRankSelectConfig::L1_WORD_SIZE) + 1) {
    init();
  }

  /*!
   * \brief Default move constructor.
   * \param other Other rank data structure.
   */
  Rank(Rank&& other) = default;

  /*!
   * \brief Default move assignment.
   * \param other Other rank data structure.
   */
  Rank& operator=(Rank&& other) = default;

  //! Destructor. Deleting manually created arrays.
  ~Rank() = default;

  /*!
   * \brief Computes rank of zeros.
   * \param index Index the rank of zeros is computed for.
   * \return Number of zeros (rank) before position \c index.
   */
  [[nodiscard("rank0 computed but not used")]] inline size_t
  rank0(size_t index) const {
    PASTA_ASSERT(index <= bit_size_, "Index outside of bit vector");
    return index - rank1(index);
  }

  /*!
   * \brief Computes rank of ones.
   * \param index Index the rank of ones is computed for.
   * \return Numbers of ones (rank) before position \c index.
   */
  [[nodiscard("rank1 computed but not used")]] inline size_t
  rank1(size_t index) const {
    PASTA_ASSERT(index <= bit_size_, "Index outside of bit vector");
    size_t offset = ((index / PopcntRankSelectConfig::L2_BIT_SIZE) * 8);
    size_t const l1_pos = index / PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l2_pos = (index % PopcntRankSelectConfig::L1_BIT_SIZE) /
                          PopcntRankSelectConfig::L2_BIT_SIZE;
    size_t result =
        l0_[index / PopcntRankSelectConfig::L0_BIT_SIZE] + l12_[l1_pos].l1;

    auto l2 = l12_[l1_pos].l2_values;
    for (size_t i = 0; i < l2_pos; ++i) {
      result += (l2 & uint16_t(0b1111111111));
      l2 >>= 10;
    }

    // It is faster to not have a specialized rank0 function when
    // optimized for zero queries, because there is no popcount for
    // zero equivalent and for all popcounts in this code, the words
    // would have to be bit-wise negated, which is more expensive than
    // the computation below.
    if constexpr (!optimize_one_or_dont_care(optimized_for)) {
      result = ((l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) +
                (l2_pos * PopcntRankSelectConfig::L2_BIT_SIZE)) -
               result;
    }
    index %= PopcntRankSelectConfig::L2_BIT_SIZE;
    PASTA_ASSERT(index < PopcntRankSelectConfig::L2_BIT_SIZE,
                 "Trying to access bits that should be "
                 "covered in an L1-block");
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
  [[nodiscard("space usage computed but not used")]] virtual size_t
  space_usage() const {
    return l0_.size() * sizeof(uint64_t) + l12_.size() * sizeof(L12Type) +
           sizeof(*this);
  }

private:
  //! Function used for initializing data structure to reduce LOCs of
  //! constructor.
  void init() {
    l0_[0] = 0;

    size_t l0_pos = 1;
    size_t l12_pos = 0;
    uint32_t l1_entry = 0UL;

    uint64_t const* data = data_;
    uint64_t const* const data_end = data_ + data_size_;

    // For each full L12-Block
    std::array<uint16_t, 3> l2_entries = {0, 0, 0};
    while (data + 32 <= data_end) {
      uint32_t new_l1_entry = l1_entry;
      for (size_t i = 0; i < 3; ++i) {
        if constexpr (optimize_one_or_dont_care(optimized_for)) {
          l2_entries[i] = popcount<8>(data);
        } else {
          l2_entries[i] = popcount_zeros<8>(data);
        }
        data += 8;
        new_l1_entry += l2_entries[i];
      }
      l12_[l12_pos++] = L12Type(l1_entry, l2_entries);
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        new_l1_entry += popcount<8>(data);
      } else {
        new_l1_entry += popcount_zeros<8>(data);
      }
      data += 8;
      l1_entry = new_l1_entry;

      if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
                     PopcntRankSelectConfig::L1_WORD_SIZE) ==
          0) [[unlikely]] {
        l0_[l0_pos] = (l0_[l0_pos - 1] + l1_entry);
        ++l0_pos;
        l1_entry = 0;
      }
    }
    // For the last not full L12-Block
    l2_entries = {0, 0, 0};
    size_t l2_pos = 0;
    while (data + 8 < data_end) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        l2_entries[l2_pos++] = popcount<8>(data);
      } else {
        l2_entries[l2_pos++] = popcount_zeros<8>(data);
      }
      data += 8;
    }
    if (l2_pos < 3) {
      while (data < data_end) {
        if constexpr (optimize_one_or_dont_care(optimized_for)) {
          l2_entries[l2_pos] += popcount<1>(data++);
        } else {
          l2_entries[l2_pos] += popcount_zeros<1>(data++);
        }
      }
    }
    l12_[l12_pos] = L12Type(l1_entry, l2_entries);
    if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
                   PopcntRankSelectConfig::L1_WORD_SIZE) ==
        0) [[unlikely]] {
      l0_[l0_pos] += (l0_[l0_pos - 1] + l1_entry);
      ++l0_pos;
      l1_entry = 0;
    } else {
      // Append sentinel (max uint64_t value) to l0_, as this makes some
      // loop-conditions in during select queries easier
      l0_[l0_pos] = std::numeric_limits<uint64_t>::max();
    }
  }
}; // class Rank

//! \}

} // namespace pasta

/******************************************************************************/
