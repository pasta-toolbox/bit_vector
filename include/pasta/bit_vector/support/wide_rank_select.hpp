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
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/bit_vector/support/popcount.hpp"
#include "pasta/bit_vector/support/select.hpp"
#include "pasta/bit_vector/support/wide_rank.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <tlx/container/simple_vector.hpp>
#include <tlx/math.hpp>
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
          FindL2WideWith find_with = FindL2WideWith::LINEAR_SEARCH,
          typename VectorType = BitVector>
class WideRankSelect : public WideRank<optimized_for> {
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using WideRank<optimized_for>::data_size_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using WideRank<optimized_for>::data_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using WideRank<optimized_for>::l1_;
  //! Get access to protected members of base class, as dependent
  //! names are not considered.
  using WideRank<optimized_for>::l2_;

  template <typename T>
  using Array = tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy>;

  // Members for the structure (needed only for select)
  //! Positions of every \c SELECT_SAMPLE_RATE zero.
  std::vector<uint32_t> samples0_;
  //! Positions of every \c SELECT_SAMPLE_RATE one.
  std::vector<uint32_t> samples1_;

public:
  //! Default constructor w/o parameter.
  WideRankSelect() = default;
  /*!
   * \brief Constructor. Creates the auxiliary information for
   * efficient rank and select queries.
   *
   * \param bv Vector of type \c VectorType the rank and select structure is
   * created for.
   */
  WideRankSelect(VectorType& bv) : WideRank<optimized_for, VectorType>(bv) {
    init();
  }

  //! Default move constructor.
  WideRankSelect(WideRankSelect&& other) = default;

  //! Default move assignment.
  WideRankSelect& operator=(WideRankSelect&& other) = default;

  //! Destructor. Deleting manually created arrays.
  ~WideRankSelect() = default;

  /*!
   * \brief Get position of specific zero, i.e., select.
   * \param rank Rank of zero the position is searched for.
   * \return Position of the rank-th zero.
   */
  [[nodiscard("select0 computed but not used")]] size_t
  select0(size_t rank) const {
    size_t const l1_end = l1_.size();
    size_t const l2_end = l2_.size();

    size_t l2_pos = ((rank - 1) / WideRankSelectConfig::SELECT_SAMPLE_RATE);
    size_t l1_pos = l2_pos / 128;

    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l1_end &&
             ((l1_pos + 1) * WideRankSelectConfig::L1_BIT_SIZE) -
                     l1_[l1_pos + 1] <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * WideRankSelectConfig::L1_BIT_SIZE) - l1_[l1_pos];
    } else {
      while (l1_pos + 1 < l1_end && l1_[l1_pos + 1] < rank) {
        ++l1_pos;
      }
      rank -= l1_[l1_pos];
    }

    l2_pos = std::max(l1_pos * 128, l2_pos);

    if constexpr (use_linear_search(find_with)) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        size_t added = 0;
        while (l2_pos + 1 < l2_end &&
               ((added + 1) * WideRankSelectConfig::L2_BIT_SIZE) -
                       l2_[l2_pos + 1] <
                   rank) {
          ++l2_pos;
          ++added;
        }
        rank -= (added * WideRankSelectConfig::L2_BIT_SIZE) - l2_[l2_pos];
      } else {
        while (l2_pos + 1 < l2_end && l2_[l2_pos + 1] < rank) {
          ++l2_pos;
        }
        rank -= l2_[l2_pos];
      }
    } else if constexpr (use_binary_search(find_with)) {
      size_t const end = std::min((l1_pos + 1) * 128, l2_end);
      size_t const iterations = tlx::integer_log2_ceil(end - l2_pos + 1);
      size_t size = 1ULL << (iterations - 1);
      size_t mid = end - size;
      size >>= 1;
      // Starting positions of the next two intervals which we are
      // going to use to prefetch the potential next comparison positions (mid).
      size_t start_next_left = l2_pos;
      size_t start_next_right = mid + 1;

      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        size_t const offset = (128 * l1_pos);
        while (size > 0) {
          // The search space does not fit into one cache line, so we do
          // some prefetching. We assume that a cache line has size 64
          // bytes (until \c
          // std::hardware_constructive_interference_size is finally
          // available in gcc) and we have an array containing 2 byte
          // elements (l2).
          if (size > 16) {
            __builtin_prefetch(&l2_[start_next_left + size], 0, 0);
            __builtin_prefetch(&l2_[start_next_right + size], 0, 0);
          }
          start_next_left =
              (rank > ((mid - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                          l2_[mid]) ?
                  start_next_right :
                  start_next_left;
          start_next_right = start_next_left + size;
          mid = start_next_left + size - 1;
          size >>= 1;
        }
        l2_pos = (rank > ((mid - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                             l2_[mid]) ?
                     mid :
                     start_next_left - 1;
        rank -= ((l2_pos - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                l2_[l2_pos];
      } else {
        while (size > 0) {
          // The search space does not fit into one cache line, so we do
          // some prefetching. We assume that a cache line has size 64
          // bytes (until \c
          // std::hardware_constructive_interference_size is finally
          // available in gcc) and we have an array containing 2 byte
          // elements (l2).
          if (size > 16) {
            __builtin_prefetch(&l2_[start_next_left + size], 0, 0);
            __builtin_prefetch(&l2_[start_next_right + size], 0, 0);
          }
          start_next_left =
              (rank > l2_[mid]) ? start_next_right : start_next_left;
          start_next_right = start_next_left + size;
          mid = start_next_left + size - 1;
          size >>= 1;
        }
        l2_pos = (rank > l2_[mid]) ? mid : start_next_left - 1;
        rank -= l2_[l2_pos];
      }
    }

    size_t last_pos = l2_pos * WideRankSelectConfig::L2_WORD_SIZE;
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
    size_t const l1_end = l1_.size();
    size_t const l2_end = l2_.size();

    size_t l2_pos = ((rank - 1) / WideRankSelectConfig::SELECT_SAMPLE_RATE);
    size_t l1_pos = l2_pos / 128;

    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l1_end && l1_[l1_pos + 1] < rank) {
        ++l1_pos;
      }
      rank -= l1_[l1_pos];
    } else {
      while (l1_pos + 1 < l1_end &&
             ((l1_pos + 1) * WideRankSelectConfig::L1_BIT_SIZE) -
                     l1_[l1_pos + 1] <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * WideRankSelectConfig::L1_BIT_SIZE) - l1_[l1_pos];
    }

    l2_pos = std::max(l1_pos * 128, l2_pos);

    if constexpr (use_linear_search(find_with)) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        while (l2_pos + 1 < l2_end && l2_[l2_pos + 1] < rank) {
          ++l2_pos;
        }
        rank -= l2_[l2_pos];
      } else {
        size_t added = 0;
        while (l2_pos + 1 < l2_end &&
               ((added + 1) * WideRankSelectConfig::L2_BIT_SIZE) -
                       l2_[l2_pos + 1] <
                   rank) {
          ++l2_pos;
          ++added;
        }
        rank -= (added * WideRankSelectConfig::L2_BIT_SIZE) - l2_[l2_pos];
      }
    } else if constexpr (use_binary_search(find_with)) {
      size_t const end = std::min((l1_pos + 1) * 128, l2_end);
      size_t const iterations = tlx::integer_log2_ceil(end - l2_pos + 1);
      size_t size = 1ULL << (iterations - 1);
      size_t mid = end - size;
      size >>= 1;
      // Starting positions of the next two intervals which we are
      // going to use to prefetch the potential next comparison positions (mid).
      size_t start_next_left = l2_pos;
      size_t start_next_right = mid + 1;

      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        while (size > 0) {
          // The search space does not fit into one cache line, so we do
          // some prefetching. We assume that a cache line has size 64
          // bytes (until \c
          // std::hardware_constructive_interference_size is finally
          // available in gcc) and we have an array containing 2 byte
          // elements (l2).
          if (size > 16) {
            __builtin_prefetch(&l2_[start_next_left + size], 0, 0);
            __builtin_prefetch(&l2_[start_next_right + size], 0, 0);
          }
          start_next_left =
              (rank > l2_[mid]) ? start_next_right : start_next_left;
          start_next_right = start_next_left + size;
          mid = start_next_left + size - 1;
          size >>= 1;
        }
        l2_pos = (rank > l2_[mid]) ? mid : start_next_left - 1;
        rank -= l2_[l2_pos];
      } else {
        size_t const offset = (128 * l1_pos);
        while (size > 0) {
          // The search space does not fit into one cache line, so we do
          // some prefetching. We assume that a cache line has size 64
          // bytes (until \c
          // std::hardware_constructive_interference_size is finally
          // available in gcc) and we have an array containing 2 byte
          // elements (l2).
          if (size > 16) {
            __builtin_prefetch(&l2_[start_next_left + size], 0, 0);
            __builtin_prefetch(&l2_[start_next_right + size], 0, 0);
          }
          start_next_left =
              (rank > ((mid - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                          l2_[mid]) ?
                  start_next_right :
                  start_next_left;
          start_next_right = start_next_left + size;
          mid = start_next_left + size - 1;
          size >>= 1;
        }
        l2_pos = (rank > ((mid - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                             l2_[mid]) ?
                     mid :
                     start_next_left - 1;
        rank -= ((l2_pos - offset) * WideRankSelectConfig::L2_BIT_SIZE) -
                l2_[l2_pos];
      }
    }

    size_t last_pos = l2_pos * WideRankSelectConfig::L2_WORD_SIZE;
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
  space_usage() const override {
    return samples0_.size() * sizeof(uint32_t) +
           samples1_.size() * sizeof(uint32_t) +
           sizeof(*this); // included in sizeof(*this)
  }

private:
  //! Function used initializing data structure to reduce LOCs of constructor.
  void init() {
    size_t const l2_end = l2_.size();
    size_t next_sample0_value = 1;
    size_t next_sample1_value = 1;
    for (size_t l2_pos = 0, offset = 0; l2_pos < l2_end; ++l2_pos) {
      offset = l1_[l2_pos / 128];
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        if ((l2_pos * WideRankSelectConfig::L2_BIT_SIZE) -
            (offset + l2_[l2_pos] >= next_sample1_value)) {
          samples0_.push_back(l2_pos - 1);
          next_sample0_value += WideRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if (offset + l2_[l2_pos] >= next_sample1_value) {
          samples1_.push_back(l2_pos - 1);
          next_sample1_value += WideRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      } else {
        if (offset + l2_[l2_pos] >= next_sample1_value) {
          samples0_.push_back(l2_pos - 1);
          next_sample0_value += WideRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if ((l2_pos * WideRankSelectConfig::L2_BIT_SIZE) -
            (offset + l2_[l2_pos] >= next_sample1_value)) {
          samples1_.push_back(l2_pos - 1);
          next_sample1_value += WideRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      }
    }
  }
}; // class   WideRankSelect

//! \}

} // namespace pasta

/******************************************************************************/
