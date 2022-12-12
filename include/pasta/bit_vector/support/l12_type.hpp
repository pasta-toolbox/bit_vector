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

#include <array>
#include <cstdint>
#include <tlx/define.hpp>

namespace pasta {

/*!
 * \brief Struct used to store L1- and L2-blocks for \c BitVectorRank and
 * \c BitVectorSelect.
 *
 * In \c L12Entry, the 32-bit L1-block entry is stored, as well as the three
 * corresponding L2-block entries (each requiring 10 bits). Note that 2 bits
 * are still unused.
 */
struct L12Type {
  //! Constructor. Empty constructor required for \c tlx::SimpleVector.
  L12Type() = default;

  /*!
   * \brief Constructor. Setting all values and packing the L2-block entries.
   * \param _l1 Value of the L1-block entry.
   * \param _l2 Values of the three L2-block entries ( as\c std::array).
   */
  L12Type(uint32_t const _l1, std::array<uint16_t, 3> const _l2)
      : l1(_l1),
        l2_values(((uint32_t(0b1111111111) & _l2[2]) << 20) |
                  ((uint32_t(0b1111111111) & _l2[1]) << 10) |
                  (uint32_t(0b1111111111) & _l2[0])) {}

  /*!
   * \brief Access operator used to access the L2-block entries individually.
   * \param index The index (0, 1, or 3) of the L2-block.
   * \return Popcount of the corresponding L2-block.
   */
  inline uint16_t operator[](size_t const index) const {
    return (l2_values >> (10 * index)) & uint16_t(0b1111111111);
  }

  //! L1-block value.
  uint32_t l1;
  //! Packed L2-block values.
  uint32_t l2_values;
} TLX_ATTRIBUTE_PACKED; // struct L12Entry

//! Check that L12Type requires only 64 bits
static_assert(sizeof(L12Type) == 8);

/*!
 * \brief Struct used to store L1- and L2-blocks for \c BitVectorFlatRank and
 * \c BitVectorFlatRankSelect.
 *
 * We store the L1-information of 8 L2-blocks in 40 bits and the
 * corresponding seven L2-information in 12 bits each. Using 12 bits
 * allows us to store the prefix sum of the popcounts in the 512-bit
 * L2-blocks. All this can be stored in 128 bits. To be precise 124 bits
 * would suffice, but the additional 4 bits allow for faster access and
 * result in exactly the same overhead the non-flat popcount rank and
 * select data structure has.
 *
 * +------+------+------+------+------+------+---+------+------+--------+
 * | 8bit | 8bit | 8bit | 8bit | 8bit | 8bit |...| 8bit | 8bit | 40 bit |
 * +------+------+------+------+------+------+---+------+------+--------+
 * | 3 * 12-bit integer | 3 * 12 bit integer |...| 12 bit int. | L1-info|
 * +------+------+------+------+------+------+---+------+------+--------+
 * |   8  | 4/4  |   8  |   8  | 4/4  |  8   |...|  8   | 4/4  + 40 Bit |
 * +------+------+------+------+------+------+---+------+------+--------+
 *
 * The order of the 12-bit integers should remain the same (as they occur
 * in the bit vector). This helps us to determine the correct block later
 * on. To this end, we have to split
 */
struct BigL12Type {
  //! Constructor. Empty constructor required for \c tlx::SimpleVector.
  BigL12Type() = default;

  /*!
   * \brief Constructor. Setting all values and packing the L2-block entries.
   * \param _l1 Value of the L1-block entry.
   * \param _l2 Values of the three L2-block entries ( as\c std::array).
   */
  BigL12Type(uint64_t const _l1, std::array<uint16_t, 7>& _l2)
      : data(((__uint128_t{0b111111111111} & _l2[6]) << 116) |
             ((__uint128_t{0b111111111111} & _l2[5]) << 104) |
             ((__uint128_t{0b111111111111} & _l2[4]) << 92) |
             ((__uint128_t{0b111111111111} & _l2[3]) << 80) |
             ((__uint128_t{0b111111111111} & _l2[2]) << 68) |
             ((__uint128_t{0b111111111111} & _l2[1]) << 56) |
             ((__uint128_t{0b111111111111} & _l2[0]) << 44) |
             ((__uint128_t{0xFFFFFFFFFFF} & _l1))) {}

  /*!
   * \brief Access operator used to access the L2-block entries individually.
   * \param index The index (0 to 6) of the L2-block.
   * \return Popcount of the corresponding L2-block.
   */
  inline uint64_t operator[](size_t const index) const {
    return (index == 0) ?
               0 :
               ((data >> ((12 * index) + 32)) & uint64_t(0b111111111111));
    // Here, shifting by 32 bits + X is sufficient, because index > 1.
  }

  /*!
   * \brief Get the L1-value of the L12-block.
   * \returns L1-value of the L12-block.
   */
  inline uint64_t l1() const {
    return uint64_t{0xFFFFFFFFFFF} & data;
  }

  //! All data of the \c BigL12Type packed into 128 bits.
  __uint128_t data;
} TLX_ATTRIBUTE_PACKED; // struct BigL12Type

//! Check that BigL12Type requires only 128 bits
static_assert(sizeof(BigL12Type) == 16);

} // namespace pasta

/******************************************************************************/
