#pragma once

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/popcount.hpp"
#include "pasta/bit_vector/support/select.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <tlx/container/simple_vector.hpp>
#include <tlx/math.hpp>
#include <vector>

namespace pasta {

template <typename VectorType = BitVector>
class HackathonSelect {

  //! Size of the bit vector the rank support is constructed for.
  size_t data_size_;
  //! Pointer to the data of the bit vector.
  VectorType::RawDataConstAccess data_;
  
public:
  //! Default constructor w/o parameter.
  HackathonSelect() = default;
  /*!
   * \brief Constructor. Creates the auxiliary information for
   * efficient rank and select queries.
   *
   * \param bv Vector of type \c VectorType the rank and select structure is
   * created for.
   */
  HackathonSelect(VectorType& bv)
    : data_size_(bv.size_),
      data_(bv.data_.data()) {
    init();
  }

  //! Default move constructor.
  HackathonSelect(HackathonSelect&& other) = default;

  //! Default move assignment.
  HackathonSelect& operator=(HackathonSelect&& other) = default;

  //! Destructor. Deleting manually created arrays.
  ~HackathonSelect() = default;

  /*!
   * \brief Get position of specific zero, i.e., select.
   * \param rank Rank of zero the position is searched for.
   * \return Position of the rank-th zero.
   */
  [[nodiscard("select0 computed but not used")]] size_t
  select0([[maybe_unused]] size_t rank) const {
    return 0;
  }

  /*!
   * \brief Get position of specific one, i.e., select.
   * \param rank Rank of one the position is searched for.
   * \return Position of the rank-th one.
   */
  [[nodiscard("select1 computed but not used")]] size_t
  select1([[maybe_unused]] size_t rank) const {
    return 0;
  }

  /*!
   * \brief Estimate for the space usage.
   * \return Number of bytes used by this data structure.
   */
  [[nodiscard("space usage computed but not used")]] size_t
  space_usage() const {
    return 0;
  }

private:
  //! Function used initializing data structure to reduce LOCs of constructor.
  void init() {
  }
}; // class HackathonSelect

//! \}

} // namespace pasta

/******************************************************************************/
