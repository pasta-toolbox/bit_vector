/*******************************************************************************
 * bit_vector_test.cpp
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

#include <pasta/bit_vector/bit_vector.hpp>
#include <tlx/die.hpp>
#include <vector>

static constexpr size_t FIB_MAX = 94; // largest 64 bit Fibonacci number

const std::array<uint64_t, FIB_MAX> fib = {
    0UL,
    1UL,
    1UL,
    2UL,
    3UL,
    5UL,
    8UL,
    13UL,
    21UL,
    34UL,
    55UL,
    89UL,
    144UL,
    233UL,
    377UL,
    610UL,
    987UL,
    1597UL,
    2584UL,
    4181UL,
    6765UL,
    10946UL,
    17711UL,
    28657UL,
    46368UL,
    75025UL,
    121393UL,
    196418UL,
    317811UL,
    514229UL,
    832040UL,
    1346269UL,
    2178309UL,
    3524578UL,
    5702887UL,
    9227465UL,
    14930352UL,
    24157817UL,
    39088169UL,
    63245986UL,
    102334155UL,
    165580141UL,
    267914296UL,
    433494437UL,
    701408733UL,
    1134903170UL,
    1836311903UL,
    2971215073UL,
    4807526976UL,
    7778742049UL,
    12586269025UL,
    20365011074UL,
    32951280099UL,
    53316291173UL,
    86267571272UL,
    139583862445UL,
    225851433717UL,
    365435296162UL,
    591286729879UL,
    956722026041UL,
    1548008755920UL,
    2504730781961UL,
    4052739537881UL,
    6557470319842UL,
    10610209857723UL,
    17167680177565UL,
    27777890035288UL,
    44945570212853UL,
    72723460248141UL,
    117669030460994UL,
    190392490709135UL,
    308061521170129UL,
    498454011879264UL,
    806515533049393UL,
    1304969544928657UL,
    2111485077978050UL,
    3416454622906707UL,
    5527939700884757UL,
    8944394323791464UL,
    14472334024676221UL,
    23416728348467685UL,
    37889062373143906UL,
    61305790721611591UL,
    99194853094755497UL,
    160500643816367088UL,
    259695496911122585UL,
    420196140727489673UL,
    679891637638612258UL,
    1100087778366101931UL,
    1779979416004714189UL,
    2880067194370816120UL,
    4660046610375530309UL,
    7540113804746346429UL,
    12200160415121876738UL,
};

void direct_access_test() {
  size_t const N = 1'000'000;

  // Test bit vector initialized with 0s
  {
    pasta::BitVector bv(N, 0);
    for (size_t i = 0; i < N; ++i) {
      die_unequal(bool{bv[i]}, bool{0});
    }
  }

  // Test bit vector initialized with 1s
  {
    pasta::BitVector bv(N, 1);
    for (size_t i = 0; i < N; ++i) {
      die_unequal(bool{bv[i]}, bool{1});
    }
  }

  // Test setting all bits to 0 and then to 1
  {
    pasta::BitVector bv(N);
    for (size_t i = 0; i < N; ++i) {
      bv[i] = 1;
      die_unequal(bool{bv[i]}, bool{1});
    }

    {
      for (size_t i = 0; i < N; ++i) {
        bv[i] = 0;
        die_unequal(bool{bv[i]}, bool{0});
      }
    }
  }

  // Test setting all bits to 1 and then to 0
  {
    pasta::BitVector bv(N);
    for (size_t i = 0; i < N; ++i) {
      bv[i] = 0;
      die_unequal(bool{bv[i]}, bool{0});
    }

    {
      for (size_t i = 0; i < N; ++i) {
        bv[i] = 1;
        die_unequal(bool{bv[i]}, bool{1});
      }
    }
  }

  // Test setting every k-th bit
  {
    pasta::BitVector bv(N);

    for (size_t k = 2; k < 7; ++k) {
      for (size_t i = 0; i < N; ++i) {
        bv[i] = (i % k == 0);
      }
      for (size_t i = 0; i < N; ++i) {
        die_unequal(bool{bv[i]}, (i % k == 0));
      }
    }
  }

  // Test encoding and decoding Fibonacci numbers using a bit vector
  {
    for (size_t k = 0; k < FIB_MAX; ++k) {
      pasta::BitVector bv(64);

      // Write the k-th Fibonacci number to the bit vector
      uint64_t v = fib[k];
      for (size_t i = 0; i < 64; ++i) {
        bv[i] = (v & 1ULL);
        v >>= 1ULL;
      }

      // Read the k-th Fibonacci number from the bit vector
      v = fib[k];
      for (size_t i = 0; i < 64; ++i) {
        die_unequal(bool{bv[i]}, v & 1ULL);
        v >>= 1ULL;
      }
    }
  }
}

void iterator_test() {
  size_t const N = 1'000'000;

  // Test bit vector initialized with 0s
  {
    pasta::BitVector bv(N, 0);
    for (auto it = bv.begin(); it != bv.end(); ++it) {
      die_unequal(bool{*it}, bool{0});
    }
  }

  // Test bit vector initialized with 1s
  {
    pasta::BitVector bv(N, 1);
    for (auto it = bv.begin(); it != bv.end(); ++it) {
      die_unequal(bool{*it}, bool{1});
    }
  }

  // Test setting all bits to 0 and then to 1
  {
    pasta::BitVector bv(N);
    for (auto it = bv.begin(); it != bv.end(); ++it) {
      *it = 1;
      die_unequal(bool{*it}, bool{1});
    }

    {
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        *it = 0;
        die_unequal(bool{*it}, bool{0});
      }
    }
  }

  // Test setting all bits to 1 and then to 0
  {
    pasta::BitVector bv(N);
    for (auto it = bv.begin(); it != bv.end(); ++it) {
      *it = 0;
      die_unequal(bool{*it}, bool{0});
    }

    {
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        *it = 1;
        die_unequal(bool{*it}, bool{1});
      }
    }
  }

  // Test setting every k-th bit
  {
    pasta::BitVector bv(N);

    for (size_t k = 2; k < 7; ++k) {
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        *it = ((bv.end() - it) % k == 0);
      }
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        die_unequal(bool{*it}, ((bv.end() - it) % k == 0));
      }
    }
  }

  // Test encoding and decoding Fibonacci numbers using a bit vector
  {
    for (size_t k = 0; k < FIB_MAX; ++k) {
      pasta::BitVector bv(64);

      // Write the k-th Fibonacci number to the bit vector
      uint64_t v = fib[k];
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        *it = (v & 1ULL);
        v >>= 1ULL;
      }

      // Read the k-th Fibonacci number from the bit vector
      v = fib[k];
      for (auto it = bv.begin(); it != bv.end(); ++it) {
        die_unequal(bool{*it}, v & 1ULL);
        v >>= 1ULL;
      }
    }
  }
}

void resize_test() {
  {
    pasta::BitVector bv(100, 0);
    bv.resize(200, 1);

    for (size_t i = 0; i < 100; ++i) {
      die_unless(bv[i] == false);
    }
    for (size_t i = 100; i < 200; ++i) {
      die_unless(bv[i] == true);
    }
  }
  {
    pasta::BitVector bv(1000, 1);
    bv.resize(2051, 0);

    for (size_t i = 0; i < 1000; ++i) {
      die_unless(bv[i] == true);
    }
    for (size_t i = 1000; i < 2051; ++i) {
      die_unless(bv[i] == false);
    }
  }
  {
    pasta::BitVector bv(221341, 1);
    bv.resize(63, 0);

    for (size_t i = 0; i < 63; ++i) {
      die_unless(bv[i] == true);
    }
  }
  {
    size_t size = 714010;
    pasta::BitVector H(size);
    std::vector<bool> content(size);
    for (size_t i = 0; i < size; i++) {
      H[i] = rand() % 2;
      content[i] = H[i];
    }
    H.resize(2 * size, 0);
    for (size_t i = 0; i < size; i++) {
      die_unless(content[i] == H[i]);
    }
    for (size_t i = size; i < H.size(); ++i) {
      die_unless(H[i] == false);
    }
  }
}

int32_t main() {
  direct_access_test();
  iterator_test();
  resize_test();

  return 0;
}

/******************************************************************************/
