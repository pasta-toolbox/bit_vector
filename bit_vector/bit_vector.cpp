/*******************************************************************************
 * pasta/container/bit_vector.cpp
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

#include <algorithm>

#include "container/bit_vector.hpp"

namespace pasta {

  BitAccess::BitAccess(uint64_t * const data, size_t const position) noexcept
    : data_(data), position_(position) { }

  BitAccess::operator bool() const noexcept {
    uint8_t const offset = uint64_t(position_) & uint64_t(0b111111);
    return (data_[position_ >> 6] >> offset) & 1ULL;
  }

  BitAccess& BitAccess::operator = (bool const value) noexcept {
    // https://graphics.stanford.edu/~seander/bithacks.html
    // (ConditionalSetOrClearBitsWithoutBranching)
    uint64_t const mask = 1ULL << (uint64_t(position_) & uint64_t(0b111111));
    data_[position_ >> 6] = (data_[position_ >> 6] & ~mask) | (-value & mask);
    return *this;
  }

  std::ptrdiff_t BitAccess::distance(BitAccess const& other) const noexcept{
    return position_ - other.position_;
  }

  bool operator == (BitAccess const& a, BitAccess const& b) noexcept {
    return a.position_ == b.position_ && a.data_ == b.data_;
  }

  bool operator != (BitAccess const& a, BitAccess const& b) noexcept {
    return a.position_ != b.position_ || a.data_ != b.data_;
  }

  BitVector::Iterator::Iterator(uint64_t * const data,
				size_t const position) noexcept
    : bit_access_(data, position) { }

  BitVector::Iterator::reference BitVector::Iterator::operator * () noexcept {
    return bit_access_;
  }

  BitVector::Iterator::pointer BitVector::Iterator::operator -> () noexcept {
    return &bit_access_;
  }

  BitVector::Iterator& BitVector::Iterator::operator ++ () noexcept {
    ++bit_access_.position_;
    return *this;
  }

  BitVector::Iterator BitVector::Iterator::operator ++ (int32_t) noexcept {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator == (BitVector::Iterator const& a,
			   BitVector::Iterator const& b) noexcept {
    return a.bit_access_ == b.bit_access_;
  }

  bool operator != (BitVector::Iterator const& a,
			   BitVector::Iterator const& b) noexcept {
    return a.bit_access_ != b.bit_access_;
  }

  BitVector::Iterator::differece_type operator - (BitVector::Iterator const& a,
						  BitVector::Iterator const& b)
    noexcept {
    return a.bit_access_.distance(b.bit_access_);
  }
  
  BitVector::BitVector(size_t const size) noexcept
    : bit_size_(size),
      size_((bit_size_>> 6) + 1),
      data_(size_),
      raw_data_(data_.data()) { }

  BitVector::BitVector(size_t const size, bool const init_value) noexcept
    : BitVector(size) {
    uint64_t const fill_value = init_value ? ~(0ULL) : 0ULL;
    std::fill_n(raw_data_, size_, fill_value);
  }
  
  BitAccess BitVector::operator [] (size_t const index) noexcept {
    return BitAccess(raw_data_, index);
  }

  BitAccess BitVector::operator [] (size_t const index) const noexcept {
    return BitAccess(raw_data_, index);
  }

  void BitVector::resize(size_t const size) noexcept {
    bit_size_ = size;
    size_ = (bit_size_ >> 6) + 1;
    data_.resize(size_);
    raw_data_ = data_.data();
  }

  TLX_ATTRIBUTE_ALWAYS_INLINE BitVector::Iterator BitVector::begin() noexcept {
    return Iterator(raw_data_, 0);
  }

  TLX_ATTRIBUTE_ALWAYS_INLINE BitVector::Iterator BitVector::end() noexcept {
    return Iterator(raw_data_, bit_size_);
  }
  
  std::span<uint64_t> BitVector::data() noexcept {
    return std::span{raw_data_, size_};
  }

  std::span<uint64_t> BitVector::data() const noexcept {
    return std::span{raw_data_, size_};
  }

  size_t BitVector::size() const noexcept {
    return bit_size_;
  }

  std::ostream& operator << (std::ostream& os, BitVector const& bv) {
    for (size_t i = 0; i < bv.bit_size_; ++i) {
      os << (bv[i] ? "1" : "0");
    }
    return os;
  }
  
} // namespace pasta

/******************************************************************************/
