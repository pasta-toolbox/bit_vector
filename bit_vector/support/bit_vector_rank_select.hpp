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

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <tlx/container/simple_vector.hpp>

#include "container/support/bit_vector_rank.hpp"

namespace pasta {

  //! Forward Declaration
  class BitVector;
  //! Forward Declaration
  struct L12Entry;

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
   */
  class BitVectorRankSelect {
    template <typename T>
    using Array = tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy>;

    //! Rank structure (requires a strict subset of data structures of select).
    BitVectorRank rank_;

    //! Size of the bit vector the select support is constructed for.
    size_t data_size_;
    //! Pointer to the data of the bit vector.
    uint64_t const * data_;

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
    BitVectorRankSelect(BitVector const& bv);

    //! Default move constructor.
    BitVectorRankSelect(BitVectorRankSelect&& other) = default;

    //! Default move assignment.
    BitVectorRankSelect& operator = (BitVectorRankSelect&& other) = default;

    //! Destructor. Deleting manually created arrays.
    ~BitVectorRankSelect();

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
     * \brief Get position of specific zero, i.e., select.
     * \param index Rank of zero the position is searched for.
     * \return Position of the rank-th zero.
     */
    [[nodiscard("select0 computed but not used")]]
    size_t select0(size_t rank) const;

    /*!
     * \brief Get position of specific one, i.e., select.
     * \param index Rank of one the position is searched for.
     * \return Position of the rank-th one.
     */
    [[nodiscard("select1 computed but not used")]]
    size_t select1(size_t rank) const;

    /*!
     * \brief Estimate for the space usage.
     * \return Number of bytes used by this data structure.
     */
    [[nodiscard("space usage computed but not used")]]
    size_t space_usage() const;

  private:
    //! Function used initializing data structure to reduce LOCs of constructor.
    void init();
  }; // class bitVectorRankSelect

  //! \}
  
} // namespace pasta

/******************************************************************************/
