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

namespace pasta {

//! \addtogroup pasta_bit_vector_configuration
//! \{

/*!
 * \brief Enum used to specify whether intrinsic functions should be used.
 *
 * Note that this does not necessarily mean that intrinsic functions are
 * faster. Please refer to the benchmarks for real world practical results
 * obtained in experiments.
 */
enum class FindL2FlatWith {
  LINEAR_SEARCH,
  BINARY_SEARCH,
  INTRINSICS
}; // enum class FindL2FlatWith

/*! \brief Helper function indicating whether a linear search
 * function should be used.
 *
 * \param find_with FindL2FlatWith indicating whether intrinsics
 * should be used.
 * \return \c true if intrinsics should be used and \c false otherwise.
 */
constexpr bool use_linear_search(FindL2FlatWith const find_with) {
  return find_with == FindL2FlatWith::LINEAR_SEARCH;
}

/*! \brief Helper function indicating whether a binary search
 * function should be used.
 *
 * \param find_with FindL2FlatWith indicating whether intrinsics
 * should be used.
 * \return \c true if intrinsics should be used and \c false otherwise.
 */
constexpr bool use_binary_search(FindL2FlatWith const find_with) {
  return find_with == FindL2FlatWith::BINARY_SEARCH;
}

/*! \brief Helper function indicating whether intrinsic function should be
 * used.
 *
 * \param find_with FindL2FlatWith indicating whether intrinsics
 * should be used.
 * \return \c true if intrinsics should be used and \c false otherwise.
 */
constexpr bool use_intrinsics(FindL2FlatWith const find_with) {
#if defined(__x86_64__)
  return find_with == FindL2FlatWith::INTRINSICS;
#else
  return false;
#endif
}

//! \)

} // namespace pasta

/******************************************************************************/
