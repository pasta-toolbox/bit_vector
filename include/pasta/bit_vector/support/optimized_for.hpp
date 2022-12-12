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

/*! \file */

//! \addtogroup pasta_bit_vector_configuration
//! \{

/*!
 * \brief Enum used to specify which queries the rank (and select) data
 * structures should be optimized for.
 *
 * Throughout the rank and select data structures, we store information about
 * the number of ones or zeros in 512 bit blocks. If we store the information
 * for ones, the queries for zeros are more complicated and vice versa. This
 * option allows to specify which query the data structure should be optimized
 * for.
 */
enum class OptimizedFor {
  //! It does not matter (both types are called equally often).
  DONT_CARE,
  //! rank_1 and select_1 queries should be optimized.
  ONE_QUERIES,
  //! rank_0 and select_9 queries should be optimized.
  ZERO_QUERIES,
}; // enum class OptimizedFor

/*!
 * \brief Helper function indicating if queries should be optimized for one
 * queries or the caller does not care.
 *
 * \param optimized_for Parameter indicating for what type of queries
 * the data structure should be optimized.
 * \return \c true if the data structure should be optimized for one queries
 * or the caller does not care for what queries the data structure is
 * optimized for. `false` otherwise.
 */
constexpr bool optimize_one_or_dont_care(OptimizedFor const optimized_for) {
  return ((optimized_for == OptimizedFor::DONT_CARE) ||
          (optimized_for == OptimizedFor::ONE_QUERIES));
}

//! \}
} // namespace pasta

/******************************************************************************/
