/*******************************************************************************
 * pasta/container/support/select.hpp
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

/*! \file */

/*!
 * \brief Select set bit in 64-bit word and return its position starting from
 * the LSB (starting from the left).
 *
 * Based on:
 * https://graphics.stanford.edu/~seander/bithacks.html#SelectPosFromMSBRank
 *
 * Starting from the left means that rank 3 in 0111, returns 0 and rank 1 in
 * the same word returns 2.
 *
 * \param data 64-bit word the bit is selected in.
 * \param rank Rank of the bit that is selected, i.e., 1st to 64-th set bit.
 * \return Position of the rank-th bit starting from the LSB (starting from
 * the left).
 */
[[nodiscard]] uint32_t select1_reverse(uint64_t const data, uint32_t rank) {
  // Do a normal parallel bit count for a 64-bit integer,
  // but store all intermediate steps.
  // a = (data & 0x5555...) + ((data >> 1) & 0x5555...);
  uint64_t const a = data - ((data >> 1) & ~0UL / 3);
  // b = (a & 0x3333...) + ((a >> 2) & 0x3333...);
  uint64_t const b = (a & ~0UL / 5) + ((a >> 2) & ~0UL / 5);
  // c = (b & 0x0f0f...) + ((b >> 4) & 0x0f0f...);
  uint64_t const c = (b + (b >> 4)) & ~0UL / 0x11;
  // d = (c & 0x00ff...) + ((c >> 8) & 0x00ff...);
  uint64_t const d = (c + (c >> 8)) & ~0UL / 0x101;
  uint32_t t = (d >> 32) + (d >> 48);
  // Now do branchless select!
  uint32_t s = 64;
  // if (rank > t) {s -= 32; rank -= t;}
  s -= ((t - rank) & 256) >> 3;
  rank -= (t & ((t - rank) >> 8));
  t = (d >> (s - 16)) & 0xff;
  // if (rank > t) {s -= 16; rank -= t;}
  s -= ((t - rank) & 256) >> 4;
  rank -= (t & ((t - rank) >> 8));
  t = (c >> (s - 8)) & 0xf;
  // if (rank > t) {s -= 8; rank -= t;}
  s -= ((t - rank) & 256) >> 5;
  rank -= (t & ((t - rank) >> 8));
  t = (b >> (s - 4)) & 0x7;
  // if (rank > t) {s -= 4; rank -= t;}
  s -= ((t - rank) & 256) >> 6;
  rank -= (t & ((t - rank) >> 8));
  t = (a >> (s - 2)) & 0x3;
  // if (rank > t) {s -= 2; rank -= t;}
  s -= ((t - rank) & 256) >> 7;
  rank -= (t & ((t - rank) >> 8));
  t = (data >> (s - 1)) & 0x1;
  // if (rank > t) s--;
  s -= ((t - rank) & 256) >> 8;
  return s - 1;
}

} // namespace pasta

/******************************************************************************/
