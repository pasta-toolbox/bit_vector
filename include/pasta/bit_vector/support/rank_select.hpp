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
 * @inproceedings{Zhou2013RankSelect,
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
 * }
 */

#pragma once

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/l12_type.hpp"
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/bit_vector/support/popcount.hpp"
#include "pasta/bit_vector/support/rank.hpp"
#include "pasta/bit_vector/support/select.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <tlx/container/simple_vector.hpp>
#include <vector>

namespace pasta {

//! \addtogroup pasta_bit_vector_rank_select
//! \{

/*!
 * \brief Select support for \c BitVector
 *
 * The rank and select support is based on popcount and described in detail by
 * Zhou et al. \cite ZhouAK2013PopcountRankSelect. The data structure consists
 * of three levels (L0, L1, and L2) that contain different information
 * regarding the popcount of the current block or all previous blocks.
 *
 * Since the rank data structures used by rank are a strict subset of the data
 * structures required by select, we provide a (ever so slightly) more
 * space-efficient rank only data structure in \c BitVectorRank.
 *
 * \tparam OptimizedFor Compile time option to optimize data structure for
 * either 0, 1, or no specific type of query.
 * \tparam VectorType Type of the vector the rank and select data structure is
 * constructed for, e.g., plain \c BitVector or a compressed bit vector.
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE,
          typename VectorType = BitVector>
class RankSelect final : public Rank<optimized_for, VectorType> {
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using Rank<optimized_for>::data_size_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using Rank<optimized_for>::data_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using Rank<optimized_for>::l0_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using Rank<optimized_for>::l12_;

  template <typename T>
  using Array = tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy>;

  // Members for the structure (needed only for select)
  //! Staring positions of the samples of zeros (w.r.t. L0-blocks)
  Array<uint64_t> samples0_pos_;
  //! Staring positions of the samples of ones (w.r.t. L0-blocks)
  Array<uint64_t> samples1_pos_;
  //! Positions of every \c SELECT_SAMPLE_RATE zero.
  std::vector<uint32_t> samples0_;
  //! Positions of every \c SELECT_SAMPLE_RATE one.
  std::vector<uint32_t> samples1_;

public:
  //! Default constructor w/o parameter.
  RankSelect() = default;

  /*!
   * \brief Constructor. Creates the auxiliary information for
   * efficient rank and select queries.
   *
   * \param bv Vector of \c VectorType the rank and select structure is created
   * for.
   */
  RankSelect(VectorType& bv)
      : Rank<optimized_for>(bv),
        samples0_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1),
        samples1_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1) {
    init();
  }

  //! Default move constructor.
  RankSelect(RankSelect&& other) = default;

  //! Default move assignment.
  RankSelect& operator=(RankSelect&& other) = default;

  //! Destructor. Deleting manually created arrays.
  ~RankSelect() = default;

  /*!
   * \brief Get position of specific zero, i.e., select.
   * \param rank Rank of zero the position is searched for.
   * \return Position of the rank-th zero.
   */
  [[nodiscard("select0 computed but not used")]] size_t
  select0(size_t rank) const {
    size_t const l0_end = l0_.size();
    size_t const l12_end = l12_.size();
    size_t l0_pos = 0;

    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l0_pos + 1 < l0_end &&
             ((l0_pos + 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                     l0_[l0_pos + 1] <
                 rank) {
        ++l0_pos;
      }
    } else {
      while (l0_pos + 1 < l0_end && l0_[l0_pos + 1] < rank) {
        ++l0_pos;
      }
    }
    if (l0_pos == l0_end) [[unlikely]] {
      return data_size_;
    }
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      rank -= (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l0_[l0_pos];
    } else {
      rank -= l0_[l0_pos];
    }

    size_t const sample_pos =
        ((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
        samples0_pos_[l0_pos];

    size_t l1_pos = samples0_[sample_pos];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
              PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
        std::min<size_t>(
            ((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
                             PopcntRankSelectConfig::L1_WORD_SIZE)),
            l12_end) -
        1;
    l1_pos = std::min<size_t>(l1_pos, l0_block_end);
    size_t l2_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l0_block_end &&
             ((l1_pos + 1) * PopcntRankSelectConfig::L1_BIT_SIZE) -
                     l12_[l1_pos + 1].l1 <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
              (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l12_[l1_pos].l1;

      auto l2 = l12_[l1_pos].l2_values;
      while (l2_pos < 3 && PopcntRankSelectConfig::L2_BIT_SIZE -
                                   (l2 & uint16_t(0b1111111111)) <
                               rank) {
        rank -=
            PopcntRankSelectConfig::L2_BIT_SIZE - (l2 & uint16_t(0b1111111111));
        l2 >>= 10;
        ++l2_pos;
      }
    } else {
      while (l1_pos + 1 < l0_block_end && l12_[l1_pos + 1].l1 < rank) {
        ++l1_pos;
      }
      rank -= l12_[l1_pos].l1;
      auto l2 = l12_[l1_pos].l2_values;
      while (l2_pos < 3 && (l2 & uint16_t(0b1111111111)) < rank) {
        rank -= (l2 & uint16_t(0b1111111111));
        l2 >>= 10;
        ++l2_pos;
      }
    }

    size_t last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                      (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t popcount = 0;

    while ((popcount = pasta::popcount_zeros<1>(data_ + last_pos)) < rank) {
      ++last_pos;
      rank -= popcount;
    }

    return (last_pos * 64) + select(~data_[last_pos], rank - 1);
  }

  /*!
   * \brief Get position of specific one, i.e., select.
   * \param rank Rank of one the position is searched for.
   * \return Position of the rank-th one.
   */
  [[nodiscard("select1 computed but not used")]] size_t
  select1(size_t rank) const {
    size_t const l0_end = l0_.size();
    size_t const l12_end = l12_.size();
    size_t l0_pos = 0;

    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l0_pos + 1 < l0_end && l0_[l0_pos + 1] < rank) {
        ++l0_pos;
      }
    } else {
      while (l0_pos + 1 < l0_end &&
             ((l0_pos + 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                     l0_[l0_pos + 1] <
                 rank) {
        ++l0_pos;
      }
    }
    if (l0_pos == l0_end) [[unlikely]] {
      return data_size_;
    }
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      rank -= l0_[l0_pos];
    } else {
      rank -= (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l0_[l0_pos];
    }

    size_t const sample_pos =
        ((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
        samples1_pos_[l0_pos];

    size_t l1_pos = samples1_[sample_pos];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
              PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
        std::min<size_t>(
            ((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
                             PopcntRankSelectConfig::L1_WORD_SIZE)),
            l12_end) -
        1;
    l1_pos = std::min(l1_pos, l0_block_end);
    size_t l2_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l0_block_end && l12_[l1_pos + 1].l1 < rank) {
        ++l1_pos;
      }
      rank -= l12_[l1_pos].l1;
      auto l2 = l12_[l1_pos].l2_values;
      while (l2_pos < 3 && (l2 & uint16_t(0b1111111111)) < rank) {
        rank -= (l2 & uint16_t(0b1111111111));
        l2 >>= 10;
        ++l2_pos;
      }
    } else {
      while (l1_pos + 1 < l0_block_end &&
             ((l1_pos + 1) * PopcntRankSelectConfig::L1_BIT_SIZE) -
                     l12_[l1_pos + 1].l1 <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
              (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l12_[l1_pos].l1;

      auto l2 = l12_[l1_pos].l2_values;
      while (l2_pos < 3 && PopcntRankSelectConfig::L2_BIT_SIZE -
                                   (l2 & uint16_t(0b1111111111)) <
                               rank) {
        rank -=
            PopcntRankSelectConfig::L2_BIT_SIZE - (l2 & uint16_t(0b1111111111));
        l2 >>= 10;
        ++l2_pos;
      }
    }

    size_t last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                      (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t popcount = 0;

    while ((popcount = pasta::popcount<1>(data_ + last_pos)) < rank) {
      ++last_pos;
      rank -= popcount;
    }
    return (last_pos * 64) + select(data_[last_pos], rank - 1);
  }

  /*!
   * \brief Estimate for the space usage.
   * \return Number of bytes used by this data structure.
   */
  [[nodiscard("space usage computed but not used")]] size_t
  space_usage() const final {
    return samples0_.size() * sizeof(uint32_t) +
           samples1_.size() * sizeof(uint32_t) +
           samples0_pos_.size() * sizeof(uint64_t) +
           samples1_pos_.size() * sizeof(uint64_t) + sizeof(*this);
  }

private:
  //! Function used initializing data structure to reduce LOCs of constructor.
  void init() {
    size_t const l12_end = l12_.size();

    size_t next_sample0_value = 1;
    size_t next_sample1_value = 1;
    for (size_t l0_pos = 0, l12_pos = 0; l12_pos < l12_end; ++l12_pos) {
      if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
                     PopcntRankSelectConfig::L1_WORD_SIZE) ==
          0) [[unlikely]] {
        samples0_pos_[l0_pos] = samples0_.size();
        samples1_pos_[l0_pos++] = samples1_.size();
        next_sample0_value = 1;
        next_sample1_value = 1;
      }
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        if ((l12_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
                ((l0_pos - 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                l12_[l12_pos].l1 >=
            next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if (l12_[l12_pos].l1 >= next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      } else {
        if (l12_[l12_pos].l1 >= next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if ((l12_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
                ((l0_pos - 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                l12_[l12_pos].l1 >=
            next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      }
    }
    // Add at least one entry.
    if (samples0_.size() == 0) [[unlikely]] {
      samples0_.push_back(0);
    }
    if (samples1_.size() == 0) [[unlikely]] {
      samples1_.push_back(0);
    }
  }
}; // class RankSelect

//! \}

} // namespace pasta

/******************************************************************************/
