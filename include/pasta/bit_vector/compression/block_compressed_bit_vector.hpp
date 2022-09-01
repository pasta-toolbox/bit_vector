/*******************************************************************************
 * This file is part of pasta::bit_vector.
 *
 * Copyright (C) 2019-2021 Florian Kurpicz <florian@kurpicz.org>
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

namespace pasta {

/*! \file */

/*!
 * \brief Interface for compressed bit vectors that use pasta's rank and select
 * interface.
 *
 * This interface is currently work in progress and should not yet be used in
 * production.
 */
class BlockCompressedBitVector {
public:
  /*!
   * \brief Constructor for a block compressed bit vector that takes a \c
   * BitVector and prepares the compression.
   *
   * To easier work with the rank and select data structures, the bit vector is
   * only compressed as soon as compress() is called. If a rank and select data
   * structure is computed for this compressed bit vector, this happens
   * automatically, e.g.,
   * ```cpp
   * BitVector bv(1000);
   * // Fill bv
   * BlockCompressedBitVector bcbv(bv);
   * FlatRankSelect(std::move(bcbv));
   * ```
   * results in a compressed bit vector `bcbv`.
   *
   * \param bv Bit vector for which the block compressed bit vector is computed.
   */
  BlockCompressedBitVector(BitVector&& bv) {}

  /*!
   * \brief Creates a compressed version of the bit vector passed in the
   * constructor and discards the uncompressed version.
   */
  void compress() {}

}; // class BlockCompressedBitVector

} // namespace pasta

/******************************************************************************/
