/*******************************************************************************
 * pasta/container/support/bit_vector_flat_rank_select.hpp
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

#pragma once

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/bit_vector_flat_rank.hpp"
#include "bit_vector/support/find_l2_flat_with.hpp"
#include "bit_vector/support/l12_type.hpp"
#include "bit_vector/support/optimized_for.hpp"
#include "bit_vector/support/popcount.hpp"
#include "bit_vector/support/select.hpp"

#include <cstddef>
#include <cstdint>
#include <emmintrin.h>
#include <immintrin.h>
#include <limits>
#include <tlx/container/simple_vector.hpp>
#include <utils/debug_asserts.hpp>
#include <vector>

namespace pasta {

//! \addtogroup pasta_bit_vectors
//! \{

/*!
 * \brief Select support for \ref BitVector that can be used as an alternative
 * to \ref BitVectorRankSelect for bit vectors up to length 2^40
 *
 * The select support is an extended and engineered version of the popcount
 * select support by Zhou et al. \cite ZhouAK2013PopcountRankSelect. Similar
 * to the \ref BitVectorFlatRank support, the highest utility array (L0) is
 * removed. For more details see \ref BitVectorFlatRank and \ref BigL12Type.
 *
 * \tparam use_intrinsic Set \c true if intrinsic functions should be used to
 * find L2-block where the select query has to search the last 512 bits.
 * Currently slower than simple loop.
 *
 * \tparam OptimizedFor Compile time option to optimize data structure for
 * either 0, 1, or neither type of query. Default is \c neither.
 * \tparam use_intrinsic \c bool flag that specifies whether intrinsics should
 * be used during select queries (currently using them is slower). Default is
 * \c false.
 */
template <OptimizedFor optimized_for = OptimizedFor::DONT_CARE,
          FindL2FlatWith find_with = FindL2FlatWith::LINEAR_SEARCH>
class BitVectorFlatRankSelect {
  template <typename T>
  using Array = tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy>;

  //! Rank structure (requires a strict subset of data structures of select).
  BitVectorFlatRank<optimized_for> rank_;

  //! Size of the bit vector the select support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  uint64_t const* data_;

  // Members for the structure (needed only for select)
  //! Positions of every \c SELECT_SAMPLE_RATE zero.
  std::vector<uint32_t> samples0_;
  //! Positions of every \c SELECT_SAMPLE_RATE one.
  std::vector<uint32_t> samples1_;

public:
  //! Default constructor w/o parameter.
  BitVectorFlatRankSelect() = default;
  /*!
   * \brief Constructor. Creates the auxiliary information for
   * efficient rank and select queries.
   *
   * \param bv \c BitVector the rank and select structure is created for.
   */
  BitVectorFlatRankSelect(BitVector const& bv)
      : rank_(bv),
        data_size_(bv.size_),
        data_(bv.data_.data()) {
    init();
  }

  //! Default move constructor.
  BitVectorFlatRankSelect(BitVectorFlatRankSelect&&) = default;

  //! Default move assignment.
  BitVectorFlatRankSelect& operator=(BitVectorFlatRankSelect&&) = default;

  //! Destructor. Deleting manually created arrays.
  ~BitVectorFlatRankSelect() = default;

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
    Array<BigL12Type> const& l12 = rank_.l12_;
    size_t const l12_end = l12.size();

    size_t const sample_pos =
        ((rank - 1) / FlattenedRankSelectConfig::SELECT_SAMPLE_RATE);
    size_t l1_pos = (sample_pos >= samples0_.size()) ? samples0_.back() :
                                                       samples0_[sample_pos];
    l1_pos += ((rank - 1) % FlattenedRankSelectConfig::SELECT_SAMPLE_RATE) /
              FlattenedRankSelectConfig::L1_BIT_SIZE;

    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l1_pos + 1 < l12_end &&
             ((l1_pos + 1) * FlattenedRankSelectConfig::L1_BIT_SIZE) -
                     l12[l1_pos + 1].l1() <
                 rank) {
        ++l1_pos;
      }
      rank -=
          (l1_pos * FlattenedRankSelectConfig::L1_BIT_SIZE) - l12[l1_pos].l1();
    } else {
      while (l1_pos + 1 < l12_end && l12[l1_pos + 1].l1() < rank) {
        ++l1_pos;
      }
      rank -= l12[l1_pos].l1();
    }

    size_t l2_pos = 0;
    if constexpr (use_intrinsics(find_with)) {
      __m128i value =
          _mm_loadu_si128(reinterpret_cast<__m128i const*>(&l12[l1_pos]));
      __m128i const shuffle_mask = _mm_setr_epi8(10,
                                                 11,
                                                 8,
                                                 9,
                                                 7,
                                                 8,
                                                 5,
                                                 6,
                                                 -1,
                                                 1,
                                                 14,
                                                 15,
                                                 13,
                                                 14,
                                                 11,
                                                 12);
      value = _mm_shuffle_epi8(value, shuffle_mask);
      // The values consisting of a complete upper byte and half a lower byte,
      // which have to be shifted to the right to obtain the correct value.
      __m128i const upper_values = _mm_srli_epi16(value, 4);
      // Mask that covers the last 12 bits of a 16 bit word
      __m128i const lower_mask = _mm_set1_epi16(uint16_t{0b0000111111111111});
      // The values consisting of a half upper byte and a complete lower byte,
      // where we have to mask the lower 12 bytes to obtain the correct value.
      __m128i const lower_values = _mm_and_si128(value, lower_mask);
      // Both [upper|lower]_values contain half of the values we want. We
      // blend them together to obtain all required values in a 128 bit word.
      value = _mm_blend_epi16(upper_values, lower_values, 0b01010101);

      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        __m128i const max_ones = _mm_setr_epi16(
            uint16_t{5 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{4 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{3 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{2 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            std::numeric_limits<int16_t>::max(), // Sentinel
            uint16_t{8 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{7 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{6 * FlattenedRankSelectConfig::L2_BIT_SIZE});

        value = _mm_sub_epi16(max_ones, value);
      } else {
        // To circumvent that the last value is a zero and thus the comparison
        // fails in the next step, we add a maximum value to this. As
        // intrinsics only consider signed integers, we have to add a signed
        // 16 bit max!
        value = _mm_insert_epi16(value, std::numeric_limits<int16_t>::max(), 4);
      }

      // We want to compare the L2-values with the remaining number of bits
      // (rank) that are remaining
      __m128i cmp_value;
      PASTA_ASSERT(rank <= std::numeric_limits<uint16_t>::max(),
                   "Rank is too large. This should not occur because in this "
                   "block the number of previous bits should reduce the "
                   "local rank further.");
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        cmp_value = _mm_set1_epi16(rank);
      } else {
        cmp_value = _mm_set1_epi16(rank - 1);
      }
      // We now have a 128 bit word, where all consecutive 16 bit words are
      // either 0 (if values is less equal) or 16_BIT_MAX (if values is
      // greater than)
      __m128i cmp_result = _mm_cmpgt_epi16(value, cmp_value);

      // Obtain the most significant bit of each 8 bit word in the
      // result of the comparison. Note that the 16 MSBs will be 0.
      // Within the other 16 bits, we have 2 zero bits for each
      // element that is less than the rank.
      uint32_t const result = _mm_movemask_epi8(cmp_result);

      // Compute the number of entries that are less than the rank
      // based on the movemask-operation above.
      l2_pos = (16 - std::popcount(result)) / 2;
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        rank -= ((l2_pos * FlattenedRankSelectConfig::L2_BIT_SIZE) -
                 l12[l1_pos][l2_pos]);
      } else {
        rank -= l12[l1_pos][l2_pos];
      }
    } else if constexpr (use_linear_search(find_with)) {
      auto tmp = l12[l1_pos].data >> 32;
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        while ((l2_pos + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                       ((tmp >> 12) & uint16_t(0b111111111111)) <
                   rank &&
               l2_pos < 7) {
          tmp >>= 12;
          ++l2_pos;
        }
      } else {
        while (((tmp >> 12) & uint16_t(0b111111111111)) < rank && l2_pos < 7) {
          tmp >>= 12;
          ++l2_pos;
        }
      }
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        rank -= (l2_pos * FlattenedRankSelectConfig::L2_BIT_SIZE) -
                (tmp & uint16_t(0b111111111111));
      } else {
        rank -= (tmp & uint16_t(0b111111111111));
      }
    } else if (use_binary_search(find_with)) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        auto tmp = l12[l1_pos].data >> 44;
        if (uint16_t const mid =
                (3 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                ((tmp >> 36) & uint16_t(0b111111111111));
            mid < rank) {
          if (uint16_t const right =
                  (5 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                  ((tmp >> 60) & uint16_t(0b111111111111));
              right < rank) {
            if (uint16_t const leaf =
                    (6 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 72) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 7;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 6;
              rank -= (right - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          } else {
            if (uint16_t const leaf =
                    (4 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 48) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 5;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 4;
              rank -= (mid - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          }
        } else {
          if (uint16_t const left =
                  (1 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                  ((tmp >> 12) & uint16_t(0b111111111111));
              left < rank) {
            if (uint16_t const leaf =
                    (2 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 24) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 3;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 2;
              rank -= (left - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          } else {
            if (uint16_t const leaf =
                    (0 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    (tmp & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 1;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          }
        }
      } else {
        auto tmp = l12[l1_pos].data >> 44;
        if (uint16_t const mid = ((tmp >> 36) & uint16_t(0b111111111111));
            mid < rank) {
          if (uint16_t const right = ((tmp >> 60) & uint16_t(0b111111111111));
              right < rank) {
            if (uint16_t const leaf = ((tmp >> 72) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 7;
              rank -= leaf;
            } else {
              l2_pos = 6;
              rank -= right;
            }
          } else {
            if (uint16_t const leaf = ((tmp >> 48) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 5;
              rank -= leaf;
            } else {
              l2_pos = 4;
              rank -= mid;
            }
          }
        } else {
          if (uint16_t const left = ((tmp >> 12) & uint16_t(0b111111111111));
              left < rank) {
            if (uint16_t const leaf = ((tmp >> 24) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 3;
              rank -= leaf;
            } else {
              l2_pos = 2;
              rank -= left;
            }
          } else {
            if (uint16_t const leaf = (tmp & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 1;
              rank -= leaf;
            }
          }
        }
      }
    } else {
      static_assert(use_linear_search(find_with) ||
                        use_binary_search(find_with) ||
                        use_intrinsics(find_with),
                    "Using unsupported search method for l2 entries");
    }

    size_t last_pos = (FlattenedRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                      (FlattenedRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t popcount = 0;

    while ((popcount = pasta::popcount_zeros<1>(data_ + last_pos)) < rank) {
      ++last_pos;
      rank -= popcount;
    }
    return (last_pos * 64) + select(~data_[last_pos], rank - 1);
  }

  /*!
   * \brief Get position of specific one, i.e., select.
   * \param index Rank of one the position is searched for.
   * \return Position of the rank-th one.
   */
  [[nodiscard("select1 computed but not used")]] size_t
  select1(size_t rank) const {
    Array<BigL12Type> const& l12 = rank_.l12_;
    Array<uint64_t> const& l0 = rank_.l0_;
    size_t const l12_end = l12.size();

    size_t const sample_pos =
        ((rank - 1) / FlattenedRankSelectConfig::SELECT_SAMPLE_RATE);
    size_t l1_pos = (sample_pos >= samples1_.size()) ? samples1_.back() :
                                                       samples1_[sample_pos];
    // std::cout << "b: l1_pos " << l1_pos << '\n';
    //  l1_pos += ((rank - 1) % FlattenedRankSelectConfig::SELECT_SAMPLE_RATE) /
    //            FlattenedRankSelectConfig::L1_BIT_SIZE;
    // std::cout << "a: l1_pos " << l1_pos << '\n';
    if constexpr (optimize_one_or_dont_care(optimized_for)) {
      while (l12[l1_pos + 1].l1() < rank && (l1_pos + 1) < l12_end) {
        ++l1_pos;
      }
      // std::cout << "f: l1_pos " << l1_pos << '\n';
      rank -= (l12[l1_pos].l1() + l0[l1_pos >> 20]);
    } else {
      while (l1_pos + 1 < l12_end &&
             ((l1_pos + 1) * FlattenedRankSelectConfig::L1_BIT_SIZE) -
                     l12[l1_pos + 1].l1() <
                 rank) {
        ++l1_pos;
      }
      rank -=
          (l1_pos * FlattenedRankSelectConfig::L1_BIT_SIZE) - l12[l1_pos].l1();
    }
    size_t l2_pos = 0;
    if constexpr (use_intrinsics(find_with)) {
      __m128i value =
          _mm_loadu_si128(reinterpret_cast<__m128i const*>(&l12[l1_pos]));
      __m128i const shuffle_mask = _mm_setr_epi8(10,
                                                 11,
                                                 8,
                                                 9,
                                                 7,
                                                 8,
                                                 5,
                                                 6,
                                                 -1,
                                                 1,
                                                 14,
                                                 15,
                                                 13,
                                                 14,
                                                 11,
                                                 12);
      value = _mm_shuffle_epi8(value, shuffle_mask);
      // The values consisting of a complete upper byte and half a lower byte,
      // which have to be shifted to the right to obtain the correct value.
      __m128i const upper_values = _mm_srli_epi16(value, 4);
      // Mask that covers the last 12 bits of a 16 bit word
      __m128i const lower_mask = _mm_set1_epi16(uint16_t{0b0000111111111111});
      // The values consisting of a half upper byte and a complete lower byte,
      // where we have to mask the lower 12 bytes to obtain the correct value.
      __m128i const lower_values = _mm_and_si128(value, lower_mask);
      // Both [upper|lower]_values contain half of the values we want. We
      // blend them together to obtain all required values in a 128 bit word.
      value = _mm_blend_epi16(upper_values, lower_values, 0b01010101);

      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        // To circumvent that the last value is a zero and thus the comparison
        // fails in the next step, we add a maximum value to this. As
        // intrinsics only consider signed integers, we have to add a signed
        // 16 bit max!
        value = _mm_insert_epi16(value, std::numeric_limits<int16_t>::max(), 4);
      } else {
        __m128i const max_ones = _mm_setr_epi16(
            uint16_t{5 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{4 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{3 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{2 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            std::numeric_limits<int16_t>::max(), // Sentinel
            uint16_t{8 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{7 * FlattenedRankSelectConfig::L2_BIT_SIZE},
            uint16_t{6 * FlattenedRankSelectConfig::L2_BIT_SIZE});

        value = _mm_sub_epi16(max_ones, value);
      }

      // We want to compare the L2-values with the remaining number of bits
      // (rank) that are remaining
      PASTA_ASSERT(rank <= std::numeric_limits<uint16_t>::max(),
                   "Rank is too large. This should not occur because in this "
                   "block the number of previous bits should reduce the "
                   "local rank further.");
      __m128i const cmp_value = _mm_set1_epi16(rank - 1);
      // We now have a 128 bit word, where all consecutive 16 bit words are
      // either 0 (if values is less equal) or 16_BIT_MAX (if values is
      // greater than)
      __m128i cmp_result = _mm_cmpgt_epi16(value, cmp_value);

      // Obtain the most significant bit of each 8 bit word in the
      // result of the comparison. Note that the 16 MSBs will be 0.
      // Within the other 16 bits, we have 2 zero bits for each
      // element that is less than the rank.
      uint32_t const result = _mm_movemask_epi8(cmp_result);

      // Compute the number of entries that are less than the rank
      // based on the movemask-operation above.
      l2_pos = (16 - std::popcount(result)) / 2;
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        rank -= l12[l1_pos][l2_pos];
      } else {
        rank -= ((l2_pos * FlattenedRankSelectConfig::L2_BIT_SIZE) -
                 l12[l1_pos][l2_pos]);
      }
    } else if constexpr (use_linear_search(find_with)) {
      auto tmp = l12[l1_pos].data >> 32;
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        while (((tmp >> 12) & uint16_t(0b111111111111)) < rank && l2_pos < 7) {
          tmp >>= 12;
          ++l2_pos;
        }
        rank -= (tmp & uint16_t(0b111111111111));
      } else {
        while ((l2_pos + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                       ((tmp >> 12) & uint16_t(0b111111111111)) <
                   rank &&
               l2_pos < 7) {
          tmp >>= 12;
          ++l2_pos;
        }
        rank -= (l2_pos * FlattenedRankSelectConfig::L2_BIT_SIZE) -
                (tmp & uint16_t(0b111111111111));
      }
    } else if constexpr (use_binary_search(find_with)) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        auto tmp = l12[l1_pos].data >> 44;
        if (uint16_t const mid = ((tmp >> 36) & uint16_t(0b111111111111));
            mid < rank) {
          if (uint16_t const right = ((tmp >> 60) & uint16_t(0b111111111111));
              right < rank) {
            if (uint16_t const leaf = ((tmp >> 72) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 7;
              rank -= leaf;
            } else {
              l2_pos = 6;
              rank -= right;
            }
          } else {
            if (uint16_t const leaf = ((tmp >> 48) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 5;
              rank -= leaf;
            } else {
              l2_pos = 4;
              rank -= mid;
            }
          }
        } else {
          if (uint16_t const left = ((tmp >> 12) & uint16_t(0b111111111111));
              left < rank) {
            if (uint16_t const leaf = ((tmp >> 24) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 3;
              rank -= leaf;
            } else {
              l2_pos = 2;
              rank -= left;
            }
          } else {
            if (uint16_t const leaf = (tmp & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 1;
              rank -= leaf;
            }
          }
        }
      } else {
        auto tmp = l12[l1_pos].data >> 44;
        if (uint16_t const mid =
                (3 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                ((tmp >> 36) & uint16_t(0b111111111111));
            mid < rank) {
          if (uint16_t const right =
                  (5 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                  ((tmp >> 60) & uint16_t(0b111111111111));
              right < rank) {
            if (uint16_t const leaf =
                    (6 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 72) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 7;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 6;
              rank -= (right - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          } else {
            if (uint16_t const leaf =
                    (4 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 48) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 5;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 4;
              rank -= (mid - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          }
        } else {
          if (uint16_t const left =
                  (1 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                  ((tmp >> 12) & uint16_t(0b111111111111));
              left < rank) {
            if (uint16_t const leaf =
                    (2 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    ((tmp >> 24) & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 3;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            } else {
              l2_pos = 2;
              rank -= (left - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          } else {
            if (uint16_t const leaf =
                    (0 + 2) * FlattenedRankSelectConfig::L2_BIT_SIZE -
                    (tmp & uint16_t(0b111111111111));
                leaf < rank) {
              l2_pos = 1;
              rank -= (leaf - FlattenedRankSelectConfig::L2_BIT_SIZE);
            }
          }
        }
      }
    } else {
      static_assert(use_linear_search(find_with) ||
                        use_binary_search(find_with) ||
                        use_intrinsics(find_with),
                    "Using unsupported search method for l2 entries");
    }

    size_t last_pos = (FlattenedRankSelectConfig::L2_WORD_SIZE * l2_pos) +
                      (FlattenedRankSelectConfig::L1_WORD_SIZE * l1_pos);
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
  space_usage() const {
    return samples0_.size() * sizeof(uint32_t) +
           samples1_.size() * sizeof(uint32_t) + rank_.space_usage() +
           sizeof(*this) - sizeof(rank_); // included in sizeof(*this)
  }

private:
  //! Function used initializing data structure to reduce LOCs of constructor.
  void init() {
    Array<BigL12Type> const& l12 = rank_.l12_;
    size_t const l12_end = rank_.l12_.size();
    size_t next_sample0_value = 1;
    size_t next_sample1_value = 1;
    for (size_t l12_pos = 0; l12_pos < l12_end; ++l12_pos) {
      if constexpr (optimize_one_or_dont_care(optimized_for)) {
        if ((l12_pos * FlattenedRankSelectConfig::L1_BIT_SIZE) -
                l12[l12_pos].l1() >=
            next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += FlattenedRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if (l12[l12_pos].l1() >= next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += FlattenedRankSelectConfig::SELECT_SAMPLE_RATE;
        }
      } else {
        if (l12[l12_pos].l1() >= next_sample0_value) {
          samples0_.push_back(l12_pos - 1);
          next_sample0_value += FlattenedRankSelectConfig::SELECT_SAMPLE_RATE;
        }
        if ((l12_pos * FlattenedRankSelectConfig::L1_BIT_SIZE) -
                l12[l12_pos].l1() >=
            next_sample1_value) {
          samples1_.push_back(l12_pos - 1);
          next_sample1_value += FlattenedRankSelectConfig::SELECT_SAMPLE_RATE;
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
  }
}; // class BitVectorFlatRankSelect

//! \}

} // namespace pasta

/******************************************************************************/
