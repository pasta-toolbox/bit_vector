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

#include "pasta/bit_vector/support/find_l2_flat_with.hpp"
#include "pasta/bit_vector/support/find_l2_wide_with.hpp"
#include "pasta/bit_vector/support/optimized_for.hpp"
#include "pasta/utils/container/aligned_vector.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <span>
#include <tlx/container/simple_vector.hpp>

namespace pasta {

/*!
 * \brief Utility class used for the access operator of the \c BitVector.
 *
 * This utility class is created when the bit vector is accessed through the
 * access operator. It cannot be created explicitly and cannot be assigned. It
 * can only be cast to bool (for read access) or be assigned a bool (for write
 * access).
 */
class BitAccess {
  //! Forward declaration.
  friend class BitVector;

  //! 64-bit word the bit is contained in.
  uint64_t* const data_;
  //! Position of the bit within the 64-bit word.
  size_t position_;

public:
  //! Deleted constructor.
  BitAccess() = delete;
  //! Avoiding implicit generation of the move constructor.
  BitAccess(BitAccess&&) = delete;

  //! Avoiding implicit copy assignment.
  BitAccess& operator=(BitAccess const&) = delete;
  //! Avoiding implicit move assignment.
  BitAccess& operator=(BitAccess&&) = delete;

  /*!
   * \brief Constructor setting 64-bit word and position of the bit in the
   * word.
   *
   * This constructor can be used when a single bit in the bit vector has to
   * be accessed and one does not want to extract the bit from the 64-bit word
   * manually.
   *
   * \param data Pointer to the 64-bit word that contains the bit.
   * \param position Position of the bit within the 64-bit word.
   */
  BitAccess(uint64_t* const data, size_t const position) noexcept
      : data_(data),
        position_(position) {}

  /*!
   * \brief User-defined conversion function to bool.
   *
   * Used for read access to a single bit.
   *
   */
  operator bool() const noexcept {
    uint8_t const offset = uint64_t(position_) & uint64_t(0b111111);
    return (data_[position_ >> 6] >> offset) & 1ULL;
  }

  /*!
   * \brief Assignment operator to set a bit within the bit vector.
   * \param value Value the bit should be written to.
   * \return This after the bit has been written.
   */
  BitAccess& operator=(bool const value) noexcept {
    // https://graphics.stanford.edu/~seander/bithacks.html
    // (ConditionalSetOrClearBitsWithoutBranching)
    uint64_t const mask = 1ULL << (uint64_t(position_) & uint64_t(0b111111));
    data_[position_ >> 6] = (data_[position_ >> 6] & ~mask) | (-value & mask);
    return *this;
  }

  /*!
   * \brief Computes the distance to another \c BitAccess instance.
   *
   * This function is only used within the iterator and usually should not be
   * required in any other situation. It has **undefined** behavior if used
   * for two instances from different \c BitVector.
   * \param other The \c BitAccess instance the distance is computed to.
   * \return The distance between \c this and the other \c BitAccess instance.
   */
  int64_t distance(BitAccess const& other) const noexcept {
    return position_ - other.position_;
  }

  /*!
   * \brief Comparison of two \c BitAccess instances for equality.
   *
   * This does not compare the value the \c BitAccess instance refers to but
   * rather if the two instances refer to the same position in the same
   * \c BitVector.
   * \param a First \c BitAccess instance.
   * \param b Second \c BitAccess instance.
   * \return \c true if both instances refer to the same position in the same
   * \c BitVector and \c false otherwise.
   */
  friend bool operator==(BitAccess const& a, BitAccess const& b) noexcept {
    return a.position_ == b.position_ && a.data_ == b.data_;
  }

  /*!
   * \brief Comparison of two \c BitAccess instances for inequality.
   *
   * This does not compare the value the \c BitAccess instance refers to but
   * rather if the two instances refer to different position or to different
   * \c BitVector.
   * \param a First \c BitAccess instance.
   * \param b Second \c BitAccess instance.
   * \return \c true if both instances refer to different position or
   * different \c BitVector and \c false otherwise.
   */
  friend bool operator!=(BitAccess const& a, BitAccess const& b) noexcept {
    return a.position_ != b.position_ || a.data_ != b.data_;
  }

private:
  /*!
   * \brief The copy constructor is private and should not be used.
   *
   * We just use the copy constructor internally for the \c BitVector. There
   * is no reason for the constructor to be called anywhere else.
   */
  BitAccess(BitAccess const&) = default;
}; // class BitAccess

//! \addtogroup pasta_bit_vector
//! \{

/*!\brief Uncompressed, highly tuned, fixed size bit vector
 *
 * The uncompressed bit vector can be used as a replacement for
 * \c std::vector<bool>, when the size of the vector is known in advance. It
 * provides an access operator.
 *
 * **Important:** If you plan on accessing the data directly, note that the
 * bits are stored in reverse order in the 64-bit words. This saves an
 * subtraction, when shifting the bits for access. To be more precise, the
 * bits are stored as follows (for simplicity, we show it for 8-bit words):
 *
 * | 7 6 5 4 3 2 1 0 | 15 14 13 12 11 10 9 8 | 23 22 21 20 19 18 17 16 | ...
 *
 * \todo Add all required support to make bit vector work as a true (fixed
 * replacement of \c std::vector<bool>).
 * \todo Create dynamic sized bit vector (true replacement of
 * \c std::vector<bool>).
 */
class BitVector {
private:
  //! Forward declaration.
  template <OptimizedFor o, typename v>
  friend class Rank;
  //! Forward declaration.
  template <OptimizedFor o, typename v>
  friend class FlatRank;
  //! Forward declaration.
  template <OptimizedFor o, typename v>
  friend class WideRank;
  //! Forward declaration.
  template <OptimizedFor o, typename v>
  friend class RankSelect;
  //! Forward declaration.
  template <OptimizedFor o, FindL2FlatWith f, typename v>
  friend class FlatRankSelect;
  //! Forward declaration.
  template <OptimizedFor o, FindL2WideWith f, typename v>
  friend class WideRankSelect;

public:
  //! Type that is used to store the raw data of the bit vector.
  using RawDataType = uint64_t;
  //! Pointer to the data type that is used to store the raw data of the bit
  //! vector.
  using RawDataPointer = RawDataType*;
  //! Type that can be used to access constant values of the raw data used to
  //! represent the bit vector.
  using RawDataConstAccess = RawDataType const*;

private:
  //! Size of the bit vector in bits.
  size_t bit_size_ = 0;
  //! Size of the underlying data used to store the bits.
  size_t size_ = 0;
  //! Array of 64-bit words used to store the content of the bit vector.
  tlx::SimpleVector<RawDataType, tlx::SimpleVectorMode::NoInitNoDestroy> data_;
  //! Pointer to the raw data of the bit vector.
  RawDataPointer raw_data_ = nullptr;

public:
  /*!
   * \brief Custom iterator for \c BitVector
   */
  struct Iterator {
    //! Iterator category.
    using iterator_category = std::forward_iterator_tag;
    //! Difference type.
    using differece_type = std::ptrdiff_t;
    //! Value type is not bit but the \c BitAccess object used for access.
    using value_type = BitAccess;
    //! Pointer type.
    using pointer = value_type*;
    //! Reference type.
    using reference = value_type&;

    /*!
     * \brief Constructor. Creates Iterator pointing at a bit and holds
     * pointer to data.
     * \param data Pointer to the beginning of the \c BitVector.
     * \param position Position the iterator is pointing at.
     */
    Iterator(uint64_t* const data, size_t const position) noexcept
        : bit_access_(data, position) {}

    /*!
     * \brief Iterator is dereferenceable. Obtain value it is pointing at.
     * \return Reference to the bit the iterator points at.
     */
    reference operator*() noexcept {
      return bit_access_;
    }

    /*!
     * \brief Iterator is dereferenceable. Obtain value it is pointing at.
     * \return Pointer to the bit the iterator points at.
     */
    pointer operator->() noexcept {
      return &bit_access_;
    }

    //! Prefix increment.
    Iterator& operator++() noexcept {
      ++bit_access_.position_;
      return *this;
    }

    //! Postfix increment.
    Iterator operator++(int32_t) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }

    //! Iterator comparison equality.
    friend bool operator==(Iterator const& a, Iterator const& b) noexcept {
      return a.bit_access_ == b.bit_access_;
    }

    //! Iterator comparison inequality.
    friend bool operator!=(Iterator const& a, Iterator const& b) noexcept {
      return a.bit_access_ != b.bit_access_;
    }

    //! Iterator distance computation.
    friend differece_type operator-(Iterator const& a,
                                    Iterator const& b) noexcept {
      return a.bit_access_.distance(b.bit_access_);
    }

  private:
    BitAccess bit_access_;
  }; // struct BitVector::Iterator

  //! Default empty constructor.
  BitVector() = default;

  //! Deleted copy constructor.
  BitVector(BitVector const&) = delete;

  //! Deleted copy assignment.
  BitVector& operator=(BitVector const&) = delete;

  /*!
   * \brief Constructor. Creates a bit vector that holds a specific, fixed
   * number of bits.
   * \param size Number of bits the bit vector contains.
   */
  BitVector(size_t const size) noexcept
      : bit_size_(size),
        size_((bit_size_ >> 6) + 1),
        data_(size_),
        raw_data_(data_.data()) {}

  /*!
   * \brief Constructor. Creates a bit vector that holds a specific, fixed
   * number of bits set to a default value.
   * \param size Number of bits the bit vector contains.
   * \param init_value Value all bits initially are set to. Either 0
   *  (\c false) or 1 (\c true).
   */
  BitVector(size_t const size, bool const init_value) noexcept
      : BitVector(size) {
    uint64_t const fill_value = init_value ? ~(0ULL) : 0ULL;
    std::fill_n(raw_data_, size_, fill_value);
  }

  /*!
   * \brief Access operator to read/write to a bit of the bit vector.
   * \param index Index of the bit to be read/write to in the bit vector.
   * \return \c BitAccess that allows to access to a single bit.
   */
  BitAccess operator[](size_t const index) noexcept {
    return BitAccess(raw_data_, index);
  }

  /*!
   * \brief Access operator to read to a bit of the bit vector.
   * \param index Index of the bit to be read to in the bit vector.
   * \return \c BitAccess that allows to read access to a single bit.
   */
  BitAccess operator[](size_t const index) const noexcept {
    return BitAccess(raw_data_, index);
  }

  /*!
   * \brief Resize the bit vector to contain \c size bits.
   * \param size Number of bits the resized bit vector contains.
   */
  void resize(size_t const size) noexcept {
    bit_size_ = size;
    size_ = (bit_size_ >> 6) + 1;
    data_.resize(size_);
    raw_data_ = data_.data();
  }

  /*!
   * \brief Resize the bit vector to contain \c size bits and fill new bits with
   * default value. \param size Number of bits the resized bit vector contains.
   * \param init_value Value all bits that are appended to the bit vector (if
   * any) will have.
   */
  void resize(size_t const size, bool const init_value) noexcept {
    size_t const old_bit_size = bit_size_;
    bit_size_ = size;
    size_ = (bit_size_ >> 6) + 1;
    data_.resize(size_);
    raw_data_ = data_.data();

    if (old_bit_size < bit_size_) {
      size_t max_bitwise = std::min(bit_size_, ((old_bit_size + 63) / 64) * 64);
      for (size_t i = old_bit_size; i < max_bitwise; ++i) {
        operator[](i) = init_value;
      }
      size_t const old_size = (old_bit_size + 63) / 64;
      uint64_t const fill_value = init_value ? ~(0ULL) : 0ULL;
      std::fill_n(raw_data_ + old_size, size_ - old_size, fill_value);
    }
  }

  /*!
   * \brief Get iterator representing the first element of the \c BitVector.
   * \return Iterator representing the first element of the \c BitVector.
   */
  Iterator begin() noexcept {
    return Iterator(raw_data_, 0);
  }

  /*!
   * \brief Get iterator representing the end of the \c BitVector.
   *
   * The \c end() method returns an iterator that may refers to an invalid
   * memory address. It should never be accessed directly.
   * \return Iterator representing the end of the \c BitVector.
   */
  Iterator end() noexcept {
    return Iterator(raw_data_, bit_size_);
  }

  /*!
   * \brief Direct access to the raw data of the bit vector.
   *
   * Note that the raw data does not contain the bits from left to right. A
   * detailed description can be found at the top of this file.
   * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
   */
  std::span<uint64_t> data() noexcept {
    return std::span{raw_data_, size_};
  }

  /*!
   * \brief Direct access to the raw data of the bit vector.
   *
   * Note that the raw data does not contain the bits from left to right. A
   * detailed description can be found at the top of this file.
   * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
   */
  std::span<uint64_t> data() const noexcept {
    return std::span{raw_data_, size_};
  }

  /*!
   * \brief Direct access to one 64-bit element of the raw data of the bit
   * vector.
   *
   * Note that the raw data does not contain the bits from left to right. A
   * detailed description can be found at the top of this file.
   * \param index Index of the 64-bit bit word that should be returned.
   * \return index-th 64-bit word of the raw bit vector data.
   */
  uint64_t data(size_t const index) const noexcept {
    return raw_data_[index];
  }

  /*!
   * \brief Estimate for the space usage.
   * \return Number of bytes used by this data structure.
   */
  [[nodiscard("space usage computed but not used")]] size_t
  space_usage() const {
    return (data_.size() * sizeof(RawDataType)) + sizeof(*this);
  }

  /*!
   * \brief Get the size of the bit vector in bits.
   * \return Size of the bit vector in bits.
   */
  size_t size() const noexcept {
    return bit_size_;
  }

  //! formatted output of the \c BitVector
  friend std::ostream& operator<<(std::ostream& os, BitVector const& bv) {
    for (size_t i = 0; i < bv.bit_size_; ++i) {
      os << (bv[i] ? "1" : "0");
    }
    return os;
  }

}; // class BitVector

//! \}

} // namespace pasta

/******************************************************************************/
