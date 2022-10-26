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
#include "pasta/bit_vector/support/find_l2_flat_with.hpp"
#include "pasta/bit_vector/support/l12_type.hpp"
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/bit_vector/support/popcount.hpp"

#include <numeric>
#include <pasta/utils/debug_asserts.hpp>
#include <tlx/container/simple_vector.hpp>

namespace pasta {

/*!
 * \ingroup pasta_bit_vector_configuration
 * \brief Static configuration for \c FlatRank and
 * \c FlatRankSelect
 */
struct FlatRankSelectConfig {
  //! Bits covered by an L2-block.
  static constexpr size_t L2_BIT_SIZE = 512;
  //! Bits covered by an L1-block.
  static constexpr size_t L1_BIT_SIZE = 8 * L2_BIT_SIZE;
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
}; // struct FlatRankSelectConfig

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
 * \tparam OptimizedFor Compile time option to optimize data structure for
 * either 0, 1, or no specific type of query.
 * \tparam VectorType Type of the vector the rank data structure is constructed
 * for, e.g., plain \c BitVector or a compressed bit vector.
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE,
          typename VectorType = BitVector>
class FlatRank {
  //! Friend class, using internal information l12_.
  template <OptimizedFor o, FindL2FlatWith f, typename v>
  friend class FlatRankSelect;

protected:
  //! Size of the bit vector the rank support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  VectorType::RawDataConstAccess data_;

  //! Array containing the information about the L1- and L2-blocks.
  tlx::SimpleVector<BigL12Type, tlx::SimpleVectorMode::NoInitNoDestroy> l12_;
  //! Number of actual existing BigL12-blocks (important for scanning)
  size_t l12_end_ = 0;

public:
  //! Default constructor w/o parameter.
  FlatRank() = default;

  /*!
   * \brief Constructor. Creates the auxiliary information for efficient rank
   * queries.
   * \param bv Vector of \c VectorType the rank structure is created for.
   */
  FlatRank(VectorType& bv)
      : data_size_(bv.size_),
        data_(bv.data_.data()),
        l12_((data_size_ / FlatRankSelectConfig::L1_WORD_SIZE) + 1) {
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
    size_t offset = ((index / 512) * 8);
    size_t const l1_pos = index / FlatRankSelectConfig::L1_BIT_SIZE;
    size_t const l2_pos = ((index % FlatRankSelectConfig::L1_BIT_SIZE) /
                           FlatRankSelectConfig::L2_BIT_SIZE);
    size_t result = l12_[l1_pos].l1() + l12_[l1_pos][l2_pos];

    // It is faster to not have a specialized rank0 function when
    // optimized for zero queries, because there is no popcount for
    // zero equivalent and for all popcounts in this code, the words
    // would have to be bit-wise negated, which is more expensive than
    // the computation below.
    if constexpr (!optimize_one_or_dont_care(optimized_for)) {
      result = ((l1_pos * FlatRankSelectConfig::L1_BIT_SIZE) +
                (l2_pos * FlatRankSelectConfig::L2_BIT_SIZE)) -
               result;
    }

    index %= FlatRankSelectConfig::L2_BIT_SIZE;
    PASTA_ASSERT(index < 512,
                 "Trying to access bits that should be "
                 "covered in an L1-block");
    for (size_t i = 0; i < index / 64; ++i) {
      result += std::popcount(data_[offset++]);
    }
    if (index %= 64; index > 0) [[likely]] {
      uint64_t const remaining = (data_[offset]) << (64 - index);
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
    return l12_.size() * sizeof(BigL12Type) + sizeof(*this);
  }

private:
  //! Function used for initializing data structure to reduce LOCs of
  //! constructor.
  void init() {
    uint64_t const* data = data_;
    uint64_t const* const data_end = data_ + data_size_;

    uint64_t l1_entry = 0ULL;
    std::array<uint16_t, 7> l2_entries = {0, 0, 0, 0, 0, 0, 0};

    while (data + 64 < data_end) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        l2_entries[0] = popcount<8>(data);
      } else {
        l2_entries[0] = popcount_zeros<8>(data);
      }
      data += 8;
      for (size_t i = 1; i < 7; ++i) {
        if constexpr (optimize_one_or_dont_care(optimized_for)) {
          l2_entries[i] = l2_entries[i - 1] + popcount<8>(data);
        } else {
          l2_entries[i] = l2_entries[i - 1] + popcount_zeros<8>(data);
        }
        data += 8;
      }
      l12_[l12_end_++] = BigL12Type(l1_entry, l2_entries);
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        l1_entry += l2_entries.back() + popcount<8>(data);
      } else {
        l1_entry += l2_entries.back() + popcount_zeros<8>(data);
      }
      data += 8;
    }
    size_t l2_pos = 0;
    l2_entries = {0, 0, 0, 0, 0, 0, 0};
    while (data + 8 < data_end) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        l2_entries[l2_pos++] = popcount<8>(data);
      } else {
        l2_entries[l2_pos++] = popcount_zeros<8>(data);
      }
      data += 8;
    }
    if (l2_pos < 7) {
      while (data < data_end) {
        if constexpr (optimize_one_or_dont_care(optimized_for)) {
          l2_entries[l2_pos] += popcount<1>(data++);
        } else {
          l2_entries[l2_pos] += popcount_zeros<1>(data++);
        }
      }
    }
    std::partial_sum(l2_entries.begin(), l2_entries.end(), l2_entries.begin());
    l12_[l12_end_++] = BigL12Type(l1_entry, l2_entries);
  }
}; // class FlatRank

//! \}

} // namespace pasta

/******************************************************************************/
