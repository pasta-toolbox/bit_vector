/*******************************************************************************
 * pasta/container/support/bit_vector_rank_select.hpp
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

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/bit_vector_rank.hpp"
#include "bit_vector/support/l12_type.hpp"
#include "bit_vector/support/optimized_for.hpp"
#include "bit_vector/support/popcount.hpp"
#include "bit_vector/support/select.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <tlx/container/simple_vector.hpp>
#include <vector>

namespace pasta {

//! \addtogroup pasta_bit_vectors
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
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE>
class BitVectorRankSelect {
  template <typename T>
  using Array = tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy>;

  //! Rank structure (requires a strict subset of data structures of select).
  BitVectorRank<optimized_for> rank_;

  //! Size of the bit vector the select support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  uint64_t const* data_;

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
  BitVectorRankSelect() = default;
  /*!
   * \brief Constructor. Creates the auxiliary information for
   * efficient rank and select queries.
   *
   * \param bv \c BitVector the rank and select structure is created for.
   */
  BitVectorRankSelect(BitVector const& bv)
      : rank_(bv),
        data_size_(bv.size_),
        data_(bv.data_.data()),
        samples0_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1),
        samples1_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1) {
    init();
  }

  //! Default move constructor.
  BitVectorRankSelect(BitVectorRankSelect&& other) = default;

  //! Default move assignment.
  BitVectorRankSelect& operator=(BitVectorRankSelect&& other) = default;

  //! Destructor. Deleting manually created arrays.
  ~BitVectorRankSelect() = default;

  /*!
   * \brief Computes rank of zeros.
   * \param index Index the rank of zeros is computed for.
   * \return Numbers of zeros (rank) before position \c index.
   */
  [[nodiscard("rank0 computed but not used")]] size_t
  rank0(size_t const index) const {
    return rank_.rank0(index);
  }

  /*!
   * \brief Computes rank of ones.
   * \param index Index the rank of ones is computed for.
   * \return Numbers of ones (rank) before position \c index.
   */
  [[nodiscard("rank1 computed but not used")]] size_t
  rank1(size_t const index) const {
    return rank_.rank1(index);
  }

  /*!
   * \brief Get position of specific zero, i.e., select.
   * \param index Rank of zero the position is searched for.
   * \return Position of the rank-th zero.
   */
  [[nodiscard("select0 computed but not used")]] size_t
  select0(size_t rank) const {
    Array<uint64_t> const& l0 = rank_.l0_;
    Array<L12Type> const& l12 = rank_.l12_;

    size_t l0_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l0_pos + 1 < l0.size() &&
             ((l0_pos + 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                     l0[l0_pos + 1] <
                 rank) {
        ++l0_pos;
      }
    } else {
      while (l0_pos + 1 < l0.size() && l0[l0_pos + 1] < rank) {
        ++l0_pos;
      }
    }
    if (l0_pos == l0.size()) [[unlikely]] {
      return data_size_;
    }
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      rank -= (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l0[l0_pos];
    } else {
      rank -= l0[l0_pos];
    }

    size_t const sample_pos =
        ((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
        samples0_pos_[l0_pos];
    size_t l1_pos = (sample_pos >= samples0_.size()) ?
                        (l0_pos * (PopcntRankSelectConfig::L0_WORD_SIZE /
                                   PopcntRankSelectConfig::L1_WORD_SIZE)) :
                        samples0_[sample_pos];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
              PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
        std::min<size_t>(
            ((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
                             PopcntRankSelectConfig::L1_WORD_SIZE)),
            l12.size()) -
        1;

    size_t l2_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos < l0_block_end &&
             ((l1_pos + 1) * PopcntRankSelectConfig::L1_BIT_SIZE) -
                     l12[l1_pos + 1].l1 <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
              (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l12[l1_pos].l1;

      while (l2_pos < 3 &&
             PopcntRankSelectConfig::L2_BIT_SIZE - l12[l1_pos][l2_pos] < rank) {
        rank -= PopcntRankSelectConfig::L2_BIT_SIZE - l12[l1_pos][l2_pos++];
      }
    } else {
      while (l1_pos + 1 < l0_block_end && l12[l1_pos + 1].l1 < rank) {
        ++l1_pos;
      }
      rank -= l12[l1_pos].l1;
      while (l2_pos < 3 && l12[l1_pos][l2_pos] < rank) {
        rank -= l12[l1_pos][l2_pos++];
      }
    }

    size_t const last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                            (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t additional_words = 0;
    size_t popcount = 0;

    while ((popcount = pasta::popcount_zeros<1>(data_ + last_pos +
                                                additional_words)) < rank) {
      ++additional_words;
      rank -= popcount;
    }

    return (PopcntRankSelectConfig::L2_BIT_SIZE * l2_pos) +
           (PopcntRankSelectConfig::L1_BIT_SIZE * l1_pos) +
           (additional_words * 64) +
           pasta::select1_reverse(~data_[last_pos + additional_words],
                                  popcount - rank + 1);
  }

  /*!
   * \brief Get position of specific one, i.e., select.
   * \param index Rank of one the position is searched for.
   * \return Position of the rank-th one.
   */
  [[nodiscard("select1 computed but not used")]] size_t
  select1(size_t rank) const {
    Array<uint64_t> const& l0 = rank_.l0_;
    Array<L12Type> const& l12 = rank_.l12_;

    size_t l0_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l0_pos + 1 < l0.size() && l0[l0_pos + 1] < rank) {
        ++l0_pos;
      }
    } else {
      while (l0_pos + 1 < l0.size() &&
             ((l0_pos + 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                     l0[l0_pos + 1] <
                 rank) {
        ++l0_pos;
      }
    }
    if (l0_pos == l0.size()) [[unlikely]] {
      return data_size_;
    }
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      rank -= l0[l0_pos];
    } else {
      rank -= (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l0[l0_pos];
    }

    size_t const sample_pos =
        ((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
        samples1_pos_[l0_pos];
    size_t l1_pos = (sample_pos >= samples1_.size()) ?
                        (l0_pos * (PopcntRankSelectConfig::L0_WORD_SIZE /
                                   PopcntRankSelectConfig::L1_WORD_SIZE)) :
                        samples1_[sample_pos];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
              PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
        std::min<size_t>(
            ((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
                             PopcntRankSelectConfig::L1_WORD_SIZE)),
            l12.size()) -
        1;

    size_t l2_pos = 0;
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l0_block_end && l12[l1_pos + 1].l1 < rank) {
        ++l1_pos;
      }
      rank -= l12[l1_pos].l1;
      while (l2_pos < 3 && l12[l1_pos][l2_pos] < rank) {
        rank -= l12[l1_pos][l2_pos++];
      }
    } else {
      while (l1_pos < l0_block_end &&
             ((l1_pos + 1) * PopcntRankSelectConfig::L1_BIT_SIZE) -
                     l12[l1_pos + 1].l1 <
                 rank) {
        ++l1_pos;
      }
      rank -= (l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
              (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l12[l1_pos].l1;

      while (l2_pos < 3 &&
             PopcntRankSelectConfig::L2_BIT_SIZE - l12[l1_pos][l2_pos] < rank) {
        rank -= PopcntRankSelectConfig::L2_BIT_SIZE - l12[l1_pos][l2_pos++];
      }
    }

    size_t const last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                            (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t additional_words = 0;
    size_t popcount = 0;

    while ((popcount = pasta::popcount<1>(data_ + last_pos +
                                          additional_words)) < rank) {
      ++additional_words;
      rank -= popcount;
    }
    return //(PopcntRankSelectConfig::L0_BIT_SIZE * l0_pos) +
        (PopcntRankSelectConfig::L2_BIT_SIZE * l2_pos) +
        (PopcntRankSelectConfig::L1_BIT_SIZE * l1_pos) +
        (additional_words * 64) +
        pasta::select1_reverse(data_[last_pos + additional_words],
                               popcount - rank + 1);
  }

  /*!
   * \brief Estimate for the space usage.
   * \return Number of bytes used by this data structure.
   */
  [[nodiscard("space usage computed but not used")]] size_t
  space_usage() const {
    return samples0_.size() * sizeof(uint32_t) +
           samples1_.size() * sizeof(uint32_t) +
           samples0_pos_.size() * sizeof(uint64_t) +
           samples1_pos_.size() * sizeof(uint64_t) + rank_.space_usage() +
           sizeof(*this) - sizeof(rank_); // included in sizeof(*this)
  }

private:
  //! Function used initializing data structure to reduce LOCs of constructor.
  void init() {
    Array<L12Type> const& l12 = rank_.l12_;

    size_t const l12_end = rank_.l12_.size();
    size_t next_sample0_value = 1;
    size_t next_sample1_value = 1;
    size_t l12_pos = 0;
    size_t l0_pos = 0;
    for (; l12_pos < l12_end; ++l12_pos) {
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
                l12[l12_pos].l1 >=
            next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if (l12[l12_pos].l1 >= next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      } else {
        if (l12[l12_pos].l1 >= next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if ((l12_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
                ((l0_pos - 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
                l12[l12_pos].l1 >=
            next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      }
    }
    // Add at least one entry.
    if (samples0_.size() == 0) {
      samples0_.push_back(0);
    }
    if (samples1_.size() == 0) {
      samples1_.push_back(0);
    }
    samples0_pos_[0] = 0;
    samples1_pos_[0] = 0;
  }
}; // class BitVectorRankSelect

//! \}

} // namespace pasta

/******************************************************************************/
