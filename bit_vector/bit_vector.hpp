/*******************************************************************************
 * pasta/container/bit_vector.hpp
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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <iostream>
#include <span>

#include <tlx/container/simple_vector.hpp>

#include "container/support/bit_vector_rank_select.hpp"

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
    uint64_t * const data_;
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
    BitAccess(uint64_t * const data, size_t const position) noexcept;
    /*!
     * \brief User-defined conversion function to bool.
     *
     * Used for read access to a single bit.
     *
     */
    operator bool() const noexcept;
    /*!
     * \brief Assignment operator to set a bit within the bit vector.
     * \param value Value the bit should be written to.
     * \return This after the bit has been written.
     */
    BitAccess& operator=(bool const value) noexcept;

    /*!
     * \brief Computes the distance to another \c BitAccess instance.
     *
     * This function is only used within the iterator and usually should not be
     * required in any other situation. It has **undefined** behavior if used
     * for two instances from different \c BitVector.
     * \param other The \c BitAccess instance the distance is computed to.
     * \return The distance between \c this and the other \c BitAccess instance.
     */
    int64_t distance(BitAccess const& other) const noexcept;

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
    friend bool operator == (BitAccess const& a, BitAccess const& b) noexcept;

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
    friend bool operator != (BitAccess const& a, BitAccess const& b) noexcept;

  private:
   /*!
     * \brief The copy constructor is private and should not be used.
     *
     * We just use the copy constructor internally for the \c BitVector. There
     * is no reason for the constructor to be called anywhere else.
     */
    BitAccess(BitAccess const&) = default;
  }; // class BitAccess

  //! \addtogroup pasta_bit_vectors
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

  public:
    using RankSelectSupport = BitVectorRankSelect;
  
  private:
    //! Forward declaration.
    friend class BitVectorRank;
    //! Forward declaration.
    friend class BitVectorRankSelect;

    //! Size of the bit vector in bits.
    size_t bit_size_;
    //! Size of the underlying data used to store the bits.
    size_t size_;
    //! Array of 64-bit words used to store the content of the bit vector.
    tlx::SimpleVector<uint64_t, tlx::SimpleVectorMode::NoInitNoDestroy> data_;
    //! Pointer to the raw data of the bit vector.
    uint64_t * raw_data_;

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
      Iterator(uint64_t * const data, size_t const positon) noexcept;

      /*!
       * \brief Iterator is dereferenceable. Obtain value it is pointing at.
       * \return Reference to the bit the iterator points at.
       */

      reference operator * () noexcept;
      /*!
       * \brief Iterator is dereferenceable. Obtain value it is pointing at.
       * \return Pointer to the bit the iterator points at.
       */
      pointer operator -> () noexcept;

      //! Prefix increment.
      Iterator& operator ++ () noexcept;
      //! Postfix increment.
      Iterator operator ++ (int32_t) noexcept;

      //! Iterator comparison equality.
      friend bool operator == (Iterator const& a, Iterator const& b) noexcept;
      //! Iterator comparison inequality.
      friend bool operator != (Iterator const& a, Iterator const& b) noexcept;
      //! Iterator distance computation.
      friend differece_type operator - (Iterator const& a, Iterator const& b) noexcept;

    private:
      BitAccess bit_access_;
    }; // struct BitVector::Iterator

    //! Default empty constructor.
    BitVector() = default;

    BitVector(BitVector const&) = delete;
    
    BitVector& operator = (BitVector const&) = delete;
    
    /*!
     * \brief Constructor. Creates a bit vector that holds a specific, fixed
     * number of bits.
     * \param size Number of bits the bit vector contains.
     */
    BitVector(size_t const size) noexcept;

    /*!
     * \brief Constructor. Creates a bit vector that holds a specific, fixed
     * number of bits set to a default value.
     * \param size Number of bits the bit vector contains.
     * \param init_value Value all bits initially are set to. Either 0
     *  (\c false) or 1 (\c true).
     */
    BitVector(size_t const size, bool const init_value) noexcept;

    /*!
     * \brief Access operator to read/write to a bit of the bit vector.
     * \param index Index of the bit to be read/write to in the bit vector.
     * \return \c BitAccess that allows to access to a single bit.
     */
    BitAccess operator [] (size_t const index) noexcept;

    /*!
     * \brief Access operator to read to a bit of the bit vector.
     * \param index Index of the bit to be read to in the bit vector.
     * \return \c BitAccess that allows to read access to a single bit.
     */
    BitAccess operator [] (size_t const index) const noexcept;

    /*!
     * \brief Resize the bit vector to contain size bits.
     * \param size Number of bits the resized bit vector contains.
     */
    void resize(size_t const size) noexcept;
    
    /*!
     * \brief Get iterator representing the first element of the \c BitVector.
     * \return Iterator representing the first element of the \c BitVector.
     */
    Iterator begin() noexcept;

    /*!
     * \brief Get iterator representing the end of the \c BitVector.
     *
     * The \c end() method returns an iterator that may refers to an invalid
     * memory address. It should never be accessed directly.
     * \return Iterator representing the end of the \c BitVector.
     */
    Iterator end() noexcept;

    /*!
     * \brief Direct access to the raw data of the bit vector.
     *
     * Note that the raw data does not contain the bits from left to right. A
     * detailed description can be found at the top of this file.
     * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
     */
    std::span<uint64_t> data() noexcept;

    /*!
     * \brief Direct access to the raw data of the bit vector.
     *
     * Note that the raw data does not contain the bits from left to right. A
     * detailed description can be found at the top of this file.
     * \return \c std::span<uint64_t> pointing to the bit vector's raw data.
     */
    std::span<uint64_t> data() const noexcept;

    /*!
     * \brief Get the size of the bit vector in bits.
     * \return Size of the bit vector in bits.
     */
    size_t size() const noexcept;

    //! formatted output of the \c BitVector
    friend std::ostream& operator << (std::ostream& os, BitVector const& bv);


  }; // class BitVector

  //! \}
  
} // namespace pasta
  
/******************************************************************************/
