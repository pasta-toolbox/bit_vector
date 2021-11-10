/*******************************************************************************
 * pasta/container/support/l12_type.hpp
 *
 * Copyright (C) 2021 Florian Kurpicz <florian@kurpicz.org>
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
  struct L12Entry {
    /*!
     *\brief Constructor. Empty constructor required for \c tlx::SimpleVector.
     */
    L12Entry() = default;

    /*!
     * \brief Constructor. Setting all values and packing the L2-block entries.
     * \param _l1 Value of the L1-block entry.
     * \param _l2 Values of the three L2-block entries ( as\c std::array).
     */
    L12Entry(uint32_t const _l1, std::array<uint16_t, 3> const _l2)
      : l1(_l1), l2_values(((uint32_t(0b1111111111) & _l2[2]) << 20) |
			   ((uint32_t(0b1111111111) & _l2[1]) << 10) |
			   (uint32_t(0b1111111111) & _l2[0])) { }

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
  
} // namespace pasta

/******************************************************************************/
