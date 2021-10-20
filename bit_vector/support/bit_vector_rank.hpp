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

#include "container/support/l12_type.hpp"

namespace pasta {

  //! Forward Declaration
  class BitVector;

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
    BitVectorRank(BitVector const& bv);

    // Default move constructor.
    BitVectorRank(BitVectorRank&& other) = default;

    // Default move assignment.
    BitVectorRank& operator = (BitVectorRank&& other) = default;

    //BitVectorRank(BitVector::Iterator begin, BitVector::Iterator end);

    //! Destructor. Deleting manually created arrays.
    ~BitVectorRank();

    /*!
     * \brief Computes rank of zeros.
     * \param index Index the rank of zeros is computed for.
     * \return Numbers of zeros (rank) before position \c index.
     */
    [[nodiscard("rank0 computed but not used")]]
    size_t rank0(size_t const index) const;

    /*!
     * \brief Computes rank of ones.
     * \param index Index the rank of ones is computed for.
     * \return Numbers of ones (rank) before position \c index.
     */
    [[nodiscard("rank1 computed but not used")]]
    size_t rank1(size_t const index) const;

    /*!
     * \brief Estimate for the space usage.
     * \return Number of bytes used by this data structure.
     */
    [[nodiscard("space usage computed but not used")]]
    size_t space_usage() const;

  private:
    //! Function used initializing data structure to reduce LOCs of constructor.
    void init();    
  }; // class BitVectorRank

  //! \}
  
} // namespace pasta

/******************************************************************************/
