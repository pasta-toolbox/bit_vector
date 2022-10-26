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

#pragma once

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/find_l2_wide_with.hpp"
#include "pasta/bit_vector/support/l12_type.hpp"
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/bit_vector/support/popcount.hpp"
#include "pasta/utils/container/aligned_vector.hpp"

#include <tlx/container/simple_vector.hpp>

namespace pasta {

/*!
 * \ingroup pasta_bit_vector_configuration
 * \brief Static configuration for \c WideRank and
 * \c WideRankSelect
 */
struct WideRankSelectConfig {
  //! Bits covered by an L2-block.
  static constexpr size_t L2_BIT_SIZE = 512;
  //! Bits covered by an L1-block.
  static constexpr size_t L1_BIT_SIZE = 128 * L2_BIT_SIZE;

  //! Number of 64-bit words covered by an L2-block.
  static constexpr size_t L2_WORD_SIZE = L2_BIT_SIZE / (sizeof(uint64_t) * 8);
  //! Number of 64-bit words covered by an L1-block.
  static constexpr size_t L1_WORD_SIZE = L1_BIT_SIZE / (sizeof(uint64_t) * 8);

  //! Sample rate of positions for faster select queries.
  static constexpr size_t SELECT_SAMPLE_RATE = 8192;
}; // struct WideRankSelectConfig

//! \addtogroup pasta_bit_vector_rank
//! \{

/*!
 * \brief %Rank support for \ref BitVector that can be used as an alternative
 * to \ref Rank for bit vectors up to length 2^40.
 *
 * The rank support is an extended and engineered version of the popcount rank
 * support by Zhou et al. \cite ZhouAK2013PopcountRankSelect. This flat rank
 * support hoverer removes the highest utility array (L0) and also groups
 * twice as many L2-blocks in a single L1-block. This allows us to store more
 * information---most importantly the number of ones w.r.t. the beginning of
 * the L1-block---in the L2-blocks. For the L1/2-block layout see
 * \ref BigL12Type.
 *
 * \tparam OptimizedFor Compile time option to optimize data
 * structure for either 0, 1, or no specific type of query.
 * \tparam VectorType Type of the vector the rank data structure is constructed
 * for, e.g., plain \c BitVector or a compressed bit vector.
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE,
          typename VectorType = BitVector>
class WideRank {
  //! Friend class, using internal information l12_.
  template <OptimizedFor o, FindL2WideWith f, typename v>
  friend class WideRankSelect;

  //! Size of the bit vector the rank support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  VectorType::RawDataConstAccess data_;

  //! Array containing the information about the L1-blocks.
  tlx::SimpleVector<uint64_t, tlx::SimpleVectorMode::NoInitNoDestroy> l1_;
  //! Array containing the information about the L2-blocks.
  tlx::SimpleVector<uint16_t, tlx::SimpleVectorMode::NoInitNoDestroy> l2_;

public:
  //! Default constructor w/o parameter.
  WideRank() = default;

  /*!
   * \brief Constructor. Creates the auxiliary information for efficient rank
   * queries.
   * \param bv Vector of type \c VectorType the rank structure is created for.
   */
  WideRank(VectorType& bv)
      : data_size_(bv.size_),
        data_(bv.data_.data()),
        l1_((data_size_ / WideRankSelectConfig::L1_WORD_SIZE) + 1),
        l2_((data_size_ / WideRankSelectConfig::L2_WORD_SIZE) + 1) {
    init();
  }

  /*!
   * \brief Computes rank of zeros.
   * \param index Index the rank of zeros is computed for.
   * \return Number of zeros (rank) before position \c index.
   */
  [[nodiscard("rank0 computed but not used")]] size_t
  rank0(size_t index) const {
    return index - rank1(index);
  }

  /*!
   * \brief Computes rank of ones.
   * \param index Index the rank of ones is computed for.
   * \return Number of ones (rank) before position \c index.
   */
  [[nodiscard("rank1 computed but not used")]] size_t
  rank1(size_t index) const {
    size_t const l1_pos = index / WideRankSelectConfig::L1_BIT_SIZE;
    size_t const l2_pos = index / WideRankSelectConfig::L2_BIT_SIZE;
    size_t result = l1_[l1_pos] + l2_[l2_pos];
    if constexpr (!optimize_one_or_dont_care(optimized_for)) {
      result = (l2_pos * WideRankSelectConfig::L2_BIT_SIZE) - result;
    }

    size_t offset = l2_pos * WideRankSelectConfig::L2_WORD_SIZE;
    size_t const full_words = (index % WideRankSelectConfig::L2_BIT_SIZE) / 64;

    for (size_t i = 0; i < full_words; ++i) {
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
  [[nodiscard("space useage computed but not used")]] virtual size_t
  space_usage() const {
    return (l1_.size() * sizeof(uint64_t)) + (l2_.size() * sizeof(uint16_t)) +
           sizeof(*this);
  }

private:
  //! Function used for initializing data structure to reduce LOCs of
  //! constructor.
  void init() {
    uint64_t const* data = data_;
    uint64_t const* const data_end = data_ + data_size_;
    l1_[0] = 0;
    l2_[0] = 0;

    size_t l1_pos = 0;
    size_t l2_pos = 0;
    size_t l2_entry = 0;
    while (data + 8 < data_end) {
      l2_[l2_pos++] = l2_entry;
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        l2_entry += popcount<8>(data);
      } else {
        l2_entry += popcount_zeros<8>(data);
      }
      data += 8;
      if (l2_pos % 128 == 0 && l2_pos > 0) [[unlikely]] {
        ++l1_pos;
        l1_[l1_pos] = l1_[l1_pos - 1] + l2_entry;
        l2_entry = 0;
      }
    }
    l2_[l2_pos++] = l2_entry;
  }
}; // class BitVectorFlatRank

//! \}

} // namespace pasta

/******************************************************************************/
