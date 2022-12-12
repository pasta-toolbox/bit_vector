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

#include <bit>
#include <cstdint>
#include <iostream>

namespace pasta {

/*! \file */

/*!
 * \brief Compute popcount of a specific number of 64-bit words.
 *
 * Note that there are no bound checks.
 *
 * \tparam Words Number of 64-bit words the popcount is computed for.
 * \param buffer Pointer to the beginning of the 64-bit words.
 * \return Popcount of the \c Words * 64 bits starting at \c buffer.
 */
template <size_t Words>
[[nodiscard]] uint64_t popcount(uint64_t const* const buffer) {
  uint64_t popcount = 0;
  for (size_t i = 0; i < Words; ++i) {
    popcount += std::popcount(buffer[i]);
  }
  return popcount;
}

/*!
 * \brief Counts the number of zero bits in a specific number of 64-bit words.
 *
 * Note that there are no bound checks.
 *
 * \tparam Words Number of 64-bit words the zeros are counted in.
 * \param buffer Pointer to the beginning of the 64-bit words.
 * \return Number of zeros in the \c Words * 64 bits starting at \c buffer.
 */
template <size_t Words>
[[nodiscard]] uint64_t popcount_zeros(uint64_t const* const buffer) {
  uint64_t popcount = 0;
  for (size_t i = 0; i < Words; ++i) {
    popcount += std::popcount(~buffer[i]);
  }
  return popcount;
}

} // namespace pasta

/******************************************************************************/
