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

#include <limits>
#include <pasta/utils/drop.hpp>
#include <queue>
#include <unordered_map>

namespace pasta {

struct FrequencyItem {
  size_t frequency;
  std::vector<uint64_t> words;

  bool operator<(FrequencyItem const& other) {
    return frequency < other.frequency;
  }

  bool operator>(FrequencyItem const& other) {
    return frequency > other.frequency;
  }
}; // struct FrequencyItem

bool operator>(FrequencyItem const& lhs, FrequencyItem const& rhs) {
  return lhs.frequency > rhs.frequency;
}

template <size_t SampleRate = 64>
class BlockCompressedBlockAccess {
  size_t block_width_;
  size_t const* sampled_pos_;
  BitVector const* compressed_bv_;

  uint64_t const* last_word_of_length_;
  size_t const* block_ends_;
  uint64_t const* blocks_;
  size_t min_length_;
  uint64_t* bv_data_;

  size_t cached_index_ = 0;
  size_t cached_bit_pos_ = 0;

  size_t tmp_cached_index_ = 0;
  size_t tmp_cached_bit_pos_ = 0;

  uint64_t tmp_cur_word = 0ULL;
  size_t tmp_cur_word_pos = 0;
  size_t tmp_bit_pos_in_word = 0;

public:
  BlockCompressedBlockAccess() = default;

  BlockCompressedBlockAccess(size_t const block_width,
                             size_t const* sampled_pos,
                             BitVector const* compressed_bv,
                             uint64_t const* last_word_of_length,
                             size_t const* block_ends,
                             uint64_t const* blocks,
                             size_t min_length)
      : block_width_(block_width),
        sampled_pos_(sampled_pos),
        compressed_bv_(compressed_bv),
        last_word_of_length_(last_word_of_length),
        block_ends_(block_ends),
        blocks_(blocks),
        min_length_(min_length),
        bv_data_(compressed_bv->data().data()) {}

  uint64_t tmp_op(size_t const index) {
    size_t blocks_till_match = index % SampleRate;
    size_t bit_pos = sampled_pos_[index / SampleRate];
    tmp_cur_word_pos = bit_pos / 64;
    tmp_bit_pos_in_word = bit_pos % 64;
    tmp_cur_word = (bv_data_[tmp_cur_word_pos] >> tmp_bit_pos_in_word);
    for (size_t i = 0; i < blocks_till_match; ++i) {
      for (size_t j = 0; j < 64 / block_width_; ++j) {
        uint64_t code_word = 0ULL;
        size_t code_length = 0;
        for (; code_length < min_length_; ++code_length) {
          code_word <<= 1;
          code_word |= (tmp_cur_word & 1ULL);
          if (++tmp_bit_pos_in_word == 64) [[unlikely]] {
            ++tmp_cur_word_pos;
            tmp_cur_word = bv_data_[tmp_cur_word_pos];
            tmp_bit_pos_in_word = 0;
          } else {
            tmp_cur_word >>= 1;
          }
        }
        while (code_word > last_word_of_length_[code_length]) {
          code_word <<= 1;
          code_word |= (tmp_cur_word & 1ULL);
          if (++tmp_bit_pos_in_word == 64) [[unlikely]] {
            ++tmp_cur_word_pos;
            tmp_cur_word = bv_data_[tmp_cur_word_pos];
            tmp_bit_pos_in_word = 0;
          } else {
            tmp_cur_word >>= 1;
          }
          ++code_length;
        }
      }
    }
    uint64_t decoded = 0ULL;
    for (size_t i = 0; i < 64 / block_width_; ++i) {
      decoded <<= block_width_;
      uint64_t code_word = 0ULL;
      size_t code_length = 0;
      for (; code_length < min_length_; ++code_length) {
        code_word <<= 1;
        code_word |= (tmp_cur_word & 1ULL);
        if (++tmp_bit_pos_in_word == 64) [[unlikely]] {
          ++tmp_cur_word_pos;
          tmp_cur_word = bv_data_[tmp_cur_word_pos];
          tmp_bit_pos_in_word = 0;
        } else {
          tmp_cur_word >>= 1;
        }
      }
      while (code_word > last_word_of_length_[code_length]) {
        code_word <<= 1;
        code_word |= (tmp_cur_word & 1ULL);
        if (++tmp_bit_pos_in_word == 64) [[unlikely]] {
          ++tmp_cur_word_pos;
          tmp_cur_word = bv_data_[tmp_cur_word_pos];
          tmp_bit_pos_in_word = 0;
        } else {
          tmp_cur_word >>= 1;
        }
        ++code_length;
      }
      size_t block_pos = block_ends_[code_length] -
                         (last_word_of_length_[code_length] - (code_word));
      uint64_t block = blocks_[block_pos];
      decoded |= block;
    }
    return decoded;
  }

  uint64_t operator[](size_t const index) {
    size_t bit_pos = 0;
    if (index == cached_index_) {
      bit_pos = cached_bit_pos_;
    } else {
      size_t blocks_till_match = index % SampleRate;
      bit_pos = sampled_pos_[index / SampleRate];
      for (size_t i = 0; i < blocks_till_match; ++i) {
        for (size_t j = 0; j < 64 / block_width_; ++j) {
          uint64_t code_word = 0ULL;
          size_t code_length = 0;
          for (; code_length < min_length_; ++code_length) {
            code_word <<= 1;
            code_word |= (*compressed_bv_)[bit_pos++];
          }
          while (code_word > last_word_of_length_[code_length]) {
            code_word <<= 1;
            code_word |= (*compressed_bv_)[bit_pos++];
            ++code_length;
          }
        }
      }
    }
    uint64_t decoded = 0ULL;
    for (size_t i = 0; i < 64 / block_width_; ++i) {
      decoded <<= block_width_;
      uint64_t code_word = 0ULL;
      size_t code_length = 0;
      for (; code_length < min_length_; ++code_length) {
        code_word <<= 1;
        code_word |= (*compressed_bv_)[bit_pos++];
      }
      while (code_word > last_word_of_length_[code_length]) {
        code_word <<= 1;
        code_word |= (*compressed_bv_)[bit_pos++];
        ++code_length;
      }
      size_t block_pos = block_ends_[code_length] -
                         (last_word_of_length_[code_length] - (code_word));
      uint64_t block = blocks_[block_pos];
      decoded |= block;
    }
    cached_index_ = index + 1;
    cached_bit_pos_ = bit_pos;
    return decoded;
  }
};

/*! \file */

/*!
 * \brief Interface for compressed bit vectors that use pasta's rank and select
 * interface.
 *
 * This interface is currently work in progress and should not yet be used in
 * production.
 */
template <size_t SampleRate = 64>
class BlockCompressedBitVector {
private:
  BitVector bv_;
  tlx::SimpleVector<size_t, tlx::SimpleVectorMode::NoInitNoDestroy>
      sampled_pos_;

  BitVector compressed_bv_;
  size_t block_width_;
  size_t min_length_;
  size_t max_length_;
  std::vector<uint64_t> blocks_;
  std::vector<size_t> block_ends_;
  std::vector<size_t> last_word_of_length_;

public:
  //! Type that can be used to access constant values of the raw data used to
  //! represent the bit vector.
  using RawDataConstAccess = BlockCompressedBlockAccess<SampleRate>;

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
  BlockCompressedBitVector(BitVector&& bv)
      : bv_(std::forward<BitVector>(bv)),
        sampled_pos_((bv_.data().size() / 64) + 1) {}

  /*!
   * \brief Creates a compressed version of the bit vector passed in the
   * constructor and discards the uncompressed version.
   */
  void compress() {
    auto space_requirements = [](auto hist, auto code_lengths) {
      size_t result = 0;
      for (auto const& block : hist) {
        result += block.second * code_lengths[block.first];
      }
      result += (hist.size() * 64);
      return result;
    };

    auto const bv_data = bv_.data();
    block_width_ = 8;
    auto min_hist = compute_block_histogram<uint8_t>(bv_.data());
    auto min_huff_code_lengths = huffman_code_lengths(min_hist);

    for (size_t i = 1; i < 4; ++i) {
      decltype(min_hist) cur_hist;
      if (i == 1) {
        cur_hist = compute_block_histogram<uint16_t>(bv_data);
      } else if (i == 2) {
        cur_hist = compute_block_histogram<uint32_t>(bv_data);
      } else /*if (i == 3)*/ {
        cur_hist = compute_block_histogram<uint64_t>(bv_data);
      }
      auto cur_huff_code_lengths = huffman_code_lengths(cur_hist);

      if (space_requirements(cur_hist, cur_huff_code_lengths) <
          space_requirements(min_hist, min_huff_code_lengths)) {
        block_width_ = 8ULL << i;
        min_hist = cur_hist;
        min_huff_code_lengths = cur_huff_code_lengths;
      }
    }

    size_t compressed_bv_size = 0;
    for (auto const& block : min_hist) {
      compressed_bv_size += block.second * min_huff_code_lengths[block.first];
    }

    auto chc = canonical_huffman_code(min_huff_code_lengths);
    compressed_bv_ = std::move(BitVector(compressed_bv_size));
    size_t compressed_bv_pos = 0;

    uint64_t bit_mask = std::numeric_limits<uint64_t>::max();
    if (block_width_ < 64) {
      bit_mask = ((1ULL << block_width_) - 1);
    }
    size_t small_blocks = 64 / block_width_;
    for (size_t i = 0; i < bv_data.size(); ++i) {
      if (i % SampleRate == 0) {
        sampled_pos_[i / SampleRate] = compressed_bv_pos;
      }
      uint64_t cur_val = bv_data[i];
      for (size_t j = 0; j != small_blocks; ++j) {
        uint64_t small_val = uint64_t{
            (cur_val >> (block_width_ * (small_blocks - j - 1))) & bit_mask};
        auto code_word = chc.encoding_map[small_val];
        for (size_t k = 0; k < code_word.length; ++k) {
          compressed_bv_[compressed_bv_pos++] =
              1ULL & (code_word.code_word >> (code_word.length - 1 - k));
        }
      }
    }

    min_length_ = chc.min_length;
    max_length_ = chc.max_length;
    last_word_of_length_ = chc.max_value_per_length;

    blocks_ = chc.blocks;
    block_ends_ = chc.block_ends;

    drop(std::move(bv_));
  }

  /*!
   * \brief Direct access to the raw data of the uncompressed bit vector.
   *
   * Note that the raw data does not contain the bits from left to right. A
   * detailed description can be found at the top of this file.
   * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
   */
  std::span<uint64_t> data() noexcept {
    return bv_.data();
  }

  /*!
   * \brief Direct access to the raw data of the uncompressed bit vector.
   *
   * Note that the raw data does not contain the bits from left to right. A
   * detailed description can be found at the top of this file.
   * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
   */
  std::span<uint64_t> data() const noexcept {
    return bv_.data();
  }

  BlockCompressedBlockAccess<SampleRate> compresed_data() const {
    return BlockCompressedBlockAccess<SampleRate>(block_width_,
                                                  sampled_pos_.data(),
                                                  &compressed_bv_,
                                                  last_word_of_length_.data(),
                                                  block_ends_.data(),
                                                  blocks_.data(),
                                                  min_length_);
  }

  /*!
   * \brief Compute space usage of the block compressed bit vector in bytes.
   * \return Space usage of block compressed bit vector in bytes.
   */
  [[nodiscard("space usage computed but not used")]] size_t
  space_usage() const {
    return (sampled_pos_.size() * sizeof(size_t)) +
           (blocks_.size() * sizeof(uint64_t)) +
           (block_ends_.size() * sizeof(size_t)) +
           (last_word_of_length_.size() * sizeof(size_t)) +
           (compressed_bv_.data().size() * sizeof(uint64_t)) + sizeof(*this);
  }

private:
  /*!
   * \brief Computes histogram of blocks of the (uncompressed bit vector). The
   * blocks are defined by \c BlockWordType.
   *
   * The resulting histogram contains the number of occurrences of consecutive
   * (non-overlapping) blocks of type \c BlockWord. The histogram does *not*
   * contain the number of 0- and 1-bits in the bit vector. \tparam
   * BlockWordType Type of the block that is considered for the histogram. Has
   * to be one of f: uint8_t, uint16_t, uint32_t, or uint64_t. \return Histogram
   * in the form of an \c std::unordered_map
   */
  template <typename BlockWordType>
  std::unordered_map<uint64_t, size_t>
  compute_block_histogram(std::span<uint64_t> bv_data) {
    static_assert(std::is_same_v<BlockWordType, uint8_t> ||
                      std::is_same_v<BlockWordType, uint16_t> ||
                      std::is_same_v<BlockWordType, uint32_t> ||
                      std::is_same_v<BlockWordType, uint64_t>,
                  "BlockWordType has to by one of: uint8_t, uint16_t, "
                  "uint32_t, or uint64_t.");

    std::unordered_map<uint64_t, size_t> hist;

    for (size_t i = 0; i < bv_data.size(); ++i) {
      uint64_t cur_val = bv_data[i];
      for (size_t j = 0; j < 8 / sizeof(BlockWordType); ++j) {
        uint64_t small_val =
            uint64_t{cur_val & std::numeric_limits<BlockWordType>::max()};
        if constexpr (sizeof(BlockWordType) < sizeof(decltype(cur_val))) {
          cur_val >>= 8 * sizeof(BlockWordType);
        }
        if (auto res = hist.find(small_val); res != hist.end()) {
          ++(res->second);
        } else {
          hist.insert({small_val, 1});
        }
      }
    }
    return hist;
  }

  std::unordered_map<uint64_t, size_t>
  huffman_code_lengths(std::unordered_map<uint64_t, size_t> const& hist) {
    // Sort single symbols by number of occurrence
    std::priority_queue<FrequencyItem,
                        std::vector<FrequencyItem>,
                        std::greater<FrequencyItem>>
        frequency_tree;

    std::unordered_map<uint64_t, size_t> code_lengths;
    for (auto const& block : hist) {
      frequency_tree.emplace(block.second, std::vector<uint64_t>{block.first});
      code_lengths.insert({block.first, 0});
    }

    while (frequency_tree.size() > 1) {
      auto ft1 = frequency_tree.top();
      frequency_tree.pop();
      auto ft2 = frequency_tree.top();
      frequency_tree.pop();
      std::copy(ft2.words.begin(),
                ft2.words.end(),
                std::back_inserter(ft1.words));
      for (auto const word : ft1.words) {
        ++code_lengths[word];
      }
      frequency_tree.emplace(ft1.frequency + ft2.frequency, ft1.words);
    }
    return code_lengths;
  }

  struct CodeWord {
    uint64_t code_word;
    size_t length;
  };

  struct CanonicalHuffmanCode {
    std::unordered_map<uint64_t, CodeWord> encoding_map;
    std::vector<uint64_t> max_value_per_length;
    std::vector<CodeWord> code_words;
    std::vector<size_t> block_ends;
    std::vector<uint64_t> blocks;
    size_t min_length;
    size_t max_length;
  }; // struct CanonicalHuffmanCode

  CanonicalHuffmanCode
  canonical_huffman_code(std::unordered_map<uint64_t, size_t> const& huffman) {
    struct BlockToLength {
      uint64_t block;
      size_t length;
    }; // struct BlockToLength

    std::vector<BlockToLength> lengths;
    for (auto const& block : huffman) {
      lengths.emplace_back(block.first, block.second);
    }

    std::sort(lengths.begin(),
              lengths.end(),
              [](BlockToLength const& lhs, BlockToLength const& rhs) {
                return lhs.length < rhs.length;
              });

    size_t min_length = lengths.front().length;
    size_t max_length = lengths.back().length;
    std::vector<CodeWord> code_words;
    std::vector<uint64_t> max_value_per_length(max_length + 1, 0);
    std::unordered_map<uint64_t, CodeWord> encoding_map;

    size_t cur_length = lengths.front().length;
    uint64_t cur_code_word = 0ULL;

    code_words.emplace_back(cur_code_word++, cur_length);
    encoding_map[lengths.front().block] = code_words.back();

    std::vector<uint64_t> blocks;
    std::vector<size_t> block_ends;
    blocks.push_back(lengths.front().block);
    for (size_t i = 0; i < min_length; ++i) {
      block_ends.push_back(0);
    }

    for (size_t i = 1; i < lengths.size(); ++i) {
      if (lengths[i].length > lengths[i - 1].length) {
        size_t length_increase = lengths[i].length - lengths[i - 1].length;
        cur_code_word <<= length_increase;
        cur_length = lengths[i].length;
        for (size_t j = 0; j < length_increase; ++j) {
          block_ends.push_back(blocks.size() - 1);
        }
      }
      // the following is updated for each code word, as it is easier to keep
      // track this way.
      max_value_per_length[cur_length] = cur_code_word;
      code_words.emplace_back(cur_code_word++, cur_length);
      encoding_map[lengths[i].block] = code_words.back();
      blocks.push_back(lengths[i].block);
    }
    block_ends.push_back(blocks.size() - 1);

    return {encoding_map,
            max_value_per_length,
            code_words,
            block_ends,
            blocks,
            min_length,
            max_length};
  }

}; // class BlockCompressedBitVector

} // namespace pasta

/******************************************************************************/
