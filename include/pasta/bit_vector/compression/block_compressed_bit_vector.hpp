/*******************************************************************************
 * pasta/container/support/bit_vector_wide_rank_select.hpp
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

namespace pasta {

class BlockCompressedBitVector {
public:
  /*! Constructor for a block compressed bit vector that takes a \c BitVector
   * and prepares the compression.
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

  /*! Creates a compressed version of the bit vector passed in the constructor
   * and discards the uncompressed version.
   */
  void compress() {}

}; // class BlockCompressedBitVector

} // namespace pasta

/******************************************************************************/
