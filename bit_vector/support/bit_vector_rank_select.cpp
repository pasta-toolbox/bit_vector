/*******************************************************************************
 * pasta/container/support/bit_vector_rank_select.cpp
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

#include <tlx/logger.hpp>

#include "container/bit_vector.hpp"
#include "container/support/bit_vector_rank_select.hpp"
#include "container/support/l12_type.hpp"
#include "container/support/popcount.hpp"
#include "container/support/select.hpp"

namespace pasta {

  BitVectorRankSelect::BitVectorRankSelect(BitVector const& bv)
    : rank_(bv),
      data_size_(bv.size_),
      data_(bv.data_.data()),
      samples0_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1),
      samples1_pos_((data_size_ / PopcntRankSelectConfig::L0_WORD_SIZE) + 1) {
    init();
  }

  BitVectorRankSelect::~BitVectorRankSelect() { }

  void BitVectorRankSelect::init() {
    Array<L12Entry> const& l12 = rank_.l12_;

    size_t const l12_end = rank_.l12_.size();
    size_t next_sample0_value = 1;
    size_t next_sample1_value = 1;
    for (size_t l12_pos = 0, l0_pos = 0; l12_pos < l12_end; ++l12_pos) {
      if (l12_pos % (PopcntRankSelectConfig::L0_WORD_SIZE /
		     PopcntRankSelectConfig::L1_WORD_SIZE) == 0) [[unlikely]] {
	samples0_pos_[l0_pos] = samples0_.size();
	samples1_pos_[l0_pos++] = samples1_.size();
	next_sample0_value = 1;
	next_sample1_value = 1;
      }
      if ((l12_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
	  ((l0_pos - 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
	  l12[l12_pos].l1 >= next_sample0_value) {
	samples0_.push_back(l12_pos - 1);
	next_sample0_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
      }
      if (l12[l12_pos].l1 >= next_sample1_value) {
        samples1_.push_back(l12_pos - 1);
	next_sample1_value += PopcntRankSelectConfig::SELECT_SAMPLE_RATE;
      }
    }
    // Add at least one entry.
    if (samples0_.size() == 0) {
      samples0_.push_back(0);
    }
    if (samples1_.size() == 0) {
      samples1_.push_back(0);
    }
    samples0_pos_[0] = 0;
    samples1_pos_[0] = 0;
  }

  size_t BitVectorRankSelect::rank0(size_t const index) const {
    return rank_.rank0(index);
  }

  size_t BitVectorRankSelect::rank1(size_t const index) const {
    return rank_.rank1(index);
  }

  size_t BitVectorRankSelect::select0(size_t rank) const {
    Array<uint64_t> const& l0 = rank_.l0_;
    Array<L12Entry> const& l12 = rank_.l12_;

    size_t l0_pos = 0;
    while (l0_pos + 1 < l0.size() &&
	   ((l0_pos + 1) * PopcntRankSelectConfig::L0_BIT_SIZE) -
	   l0[l0_pos + 1] < rank) {
      ++l0_pos;
    }
    if (l0_pos == l0.size()) [[unlikely]] {
      return data_size_;
    }
    rank -= (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l0[l0_pos];

    size_t l1_pos =
      samples0_[((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
		samples0_pos_[l0_pos]];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
      PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
      std::min<size_t>(((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
					PopcntRankSelectConfig::L1_WORD_SIZE)),
		       l12.size()) - 1;
    while (l1_pos < l0_block_end &&
	   ((l1_pos + 1) * PopcntRankSelectConfig::L1_BIT_SIZE) -
	   l12[l1_pos + 1].l1 < rank) {
      ++l1_pos;
    }
    rank -= (l1_pos * PopcntRankSelectConfig::L1_BIT_SIZE) -
      (l0_pos * PopcntRankSelectConfig::L0_BIT_SIZE) - l12[l1_pos].l1;

    size_t l2_pos = 0;
    while (l2_pos < 3 && PopcntRankSelectConfig::L2_BIT_SIZE -
	    l12[l1_pos][l2_pos] < rank) {
      rank -= PopcntRankSelectConfig::L2_BIT_SIZE - l12[l1_pos][l2_pos++];
    }

    size_t const last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
      (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t additional_words = 0;
    size_t popcount = 0;

    while ((popcount = pasta::popcount_zeros<1>(data_ + last_pos +
						additional_words)) < rank) {
      ++additional_words;
      rank -= popcount;
    }

    return (PopcntRankSelectConfig::L2_BIT_SIZE * l2_pos) +
      (PopcntRankSelectConfig::L1_BIT_SIZE * l1_pos) + (additional_words * 64) +
      pasta::select1_reverse(~data_[last_pos + additional_words],
			     popcount - rank + 1);
  }

  size_t BitVectorRankSelect::select1(size_t rank) const {
    Array<uint64_t> const& l0 = rank_.l0_;
    Array<L12Entry> const& l12 = rank_.l12_;

    size_t l0_pos = 0;
    while (l0_pos + 1 < l0.size() && l0[l0_pos + 1] < rank) {
      ++l0_pos;
    }
    if (l0_pos == l0.size()) [[unlikely]] {
      return data_size_;
    }
    rank -= l0[l0_pos];

    size_t l1_pos =
      samples1_[((rank - 1) / PopcntRankSelectConfig::SELECT_SAMPLE_RATE) +
		samples1_pos_[l0_pos]];
    l1_pos += ((rank - 1) % PopcntRankSelectConfig::SELECT_SAMPLE_RATE) /
      PopcntRankSelectConfig::L1_BIT_SIZE;
    size_t const l0_block_end =
      std::min<size_t>(((l0_pos + 1) * (PopcntRankSelectConfig::L0_WORD_SIZE /
					PopcntRankSelectConfig::L1_WORD_SIZE)),
		       l12.size()) - 1;
    while (l1_pos < l0_block_end && l12[l1_pos + 1].l1 < rank) {
      ++l1_pos;
    }
    rank -= l12[l1_pos].l1;

    size_t l2_pos = 0;
    while (l2_pos < 3 && l12[l1_pos][l2_pos] < rank) {
      rank -= l12[l1_pos][l2_pos++];
    }

    size_t const last_pos = (PopcntRankSelectConfig::L2_WORD_SIZE * l2_pos) +
      (PopcntRankSelectConfig::L1_WORD_SIZE * l1_pos);
    size_t additional_words = 0;
    size_t popcount = 0;

    while ((popcount = pasta::popcount<1>(data_ + last_pos +
					  additional_words)) < rank) {
      ++additional_words;
      rank -= popcount;
    }
    return (PopcntRankSelectConfig::L2_BIT_SIZE * l2_pos) +
      (PopcntRankSelectConfig::L1_BIT_SIZE * l1_pos) + (additional_words * 64) +
      pasta::select1_reverse(data_[last_pos + additional_words],
			     popcount - rank + 1);
  }

  size_t BitVectorRankSelect::space_usage() const {
    return samples0_.size() * sizeof(uint32_t)
        + samples1_.size() * sizeof(uint32_t)
        + samples0_pos_.size() * sizeof(uint64_t)
        + samples1_pos_.size() * sizeof(uint64_t)
        + rank_.space_usage() + sizeof(*this)
        - sizeof(rank_); // included in sizeof(*this)
  }

} // namespace pasta

/******************************************************************************/
