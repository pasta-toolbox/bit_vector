/*******************************************************************************
 * pasta/bit_vector/support/find_l2_wide_with.hpp
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

namespace pasta {

/*!
 * \brief Enum used to specify whether intrinsic functions should be used.
 *
 * Note that this does not necessarily mean that intrinsic functions are
 * faster. Please refer to the benchmarks for real world practical results
 * obtained in experiments.
 */
enum class FindL2WideWith {
  LINEAR_SEARCH,
  BINARY_SEARCH,
}; // enum class FindL2WideWith

/*! \brief Helper function indicating whether a linear search
 * function should be used.
 *
 * \param find_with \ref FindL2WideWith indicating whether intrinsics
 * should be used.
 * \return \c true if intrinsics should be used and \c false otherwise.
 */
constexpr bool use_linear_search(FindL2WideWith const find_with) {
  return find_with == FindL2WideWith::LINEAR_SEARCH;
}

/*! \brief Helper function indicating whether a binary search
 * function should be used.
 *
 * \param find_with \ref FindL2WideWith indicating whether intrinsics
 * should be used.
 * \return \c true if intrinsics should be used and \c false otherwise.
 */
constexpr bool use_binary_search(FindL2WideWith const find_with) {
  return find_with == FindL2WideWith::BINARY_SEARCH;
}

} // namespace pasta

/******************************************************************************/
