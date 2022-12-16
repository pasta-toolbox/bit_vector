# pasta::bit_vector

<p align="center">
   <img width=250 height=175 src="https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/logo_pasta_bit_vector.svg" />
</p>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/419381885.svg)](https://zenodo.org/badge/latestdoi/419381885)
[![pasta::bit_vector CI](https://github.com/pasta-toolbox/bit_vector/actions/workflows/ctest.yml/badge.svg)](https://github.com/pasta-toolbox/bit_vector/actions/workflows/ctest.yml)
[![codecov](https://codecov.io/gh/pasta-toolbox/bit_vector/branch/main/graph/badge.svg?token=2QD6ME44SU)](https://codecov.io/gh/pasta-toolbox/bit_vector)

This header-only library contains a highly tuned (uncompressed) bit vector implementation with multiple space efficient rank and select support data structures.
Our fastest rank and select support has a space overhead of only ~3.51% and makes use of data level parallelism via SIMD instructions.

If you use this code in a scientific context, please cite our paper.
```bibtex
@inproceedings{Kurpicz2022CompactRankSelect,
  author    = {Florian Kurpicz},
  title     = {Engineering Compact Data Structures for Rank and Select Queries on Bit Vectors},
  booktitle = {{SPIRE}},
  series    = {Lecture Notes in Computer Science},
  volume    = {13617},
  pages     = {257--272},
  publisher = {Springer},
  year      = {2022},
  doi       = {10.1007/978-3-031-20643-6\_19}
}
```

## Contents
This repository contains the following algorithms and data structures.
Our [documentation][] contains in-depth on the usage of all these algorithms and data structures including easy to follow examples.
You can find the example in the screenshot below as text, too.

[![Screenshot Documentation](https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/screenshot_documentation_v1.0.0.png)](https://www.pasta-toolbox.org/bit_vector/)

### Bit Vectors
Bit vectors play an important role in many compressed text indices, e.g., the FM-index.
This repository contains the following bit vector implementations:

- highly tuned [uncompressed bit vector][] with access operator
- compact [rank](include/pasta/bit_vector/support/rank.hpp) and [select](include/pasta/bit_vector/support/rank_select.hpp) support for the uncompressed bit vector based on

> Dong Zhou and David G. Andersen and Michael Kaminsky,
> Space-Efficient, High-Performance Rank and Select Structures on Uncompressed Bit Sequences,
> SEA 2013.

- improved [rank](include/pasta/bit_vector/support/flat_rank.hpp) and [select](include/pasta/bit_vector/support/flat_rank_select.hpp) support requiring the same amount of memory but providing faster rank (up to 8% speedup) and select (up to 16.5% speedup) queries, and
- a very fast [rank](include/pasta/bit_vector/support/wide_rank.hpp) support that can also answer [select](include/pasta/bit_vector/support/wide_rank_select.hpp) queries.

[uncompressed bit vector]: include/pasta/bit_vector/bit_vector.hpp

### Easy to Use

Since this is a header-only library, you have to simply add it to your projects include path to use it.
A small example can be found below.
We refer to the [documentation][] for more information.

  ```cpp
  #include <pasta/bit_vector/bit_vector.hpp>
  #include <pasta/bit_vector/support/flat_rank_select.hpp>

  // Create a bit vector of size 1000 containing only zeros and flip every other bit.
  pasta::BitVector bv(1000, 0);
  for (size_t i = 0; i < bv.size(); ++i) {
    if (i % 2 == 0) {  bv[i] = 1; }
  }
  // Print the bit vector to see that everything worked ;-)
  for (auto it = bv.begin(); it != bv.end(); ++it) {
    std::cout << ((*it == true) ? '1' : '0');
  }
  std::cout << std::endl;

  // Create rank and select support and print the result of some queries.
  pasta::FlatRankSelect rs(bv);
  std::cout << rs.rank0(250) << ", " << rs.rank1(250)
            << ", "
            << rs.select0(250) << ", " << rs.rank1(250)
            << std::endl;
  ```

### Benchmarks and Tests

There exist an easy to use [benchmark][], which helps to compare the implementations in this repository.
To build the benchmark, run the CMake command with `-DPASTA_BIT_VECTOR_BUILD_BENCHMARKS=On`.
Our tests are contained in the folder [tests][].
To build the tests, run the CMake command with `-DPASTA_BIT_VECTOR_BUILD_TESTS=On`.

We also conducted an extensive experimental evaluation.
To this end, we use our [rank and select benchmark][] where we compare our implementations with many other compact rank and select data structures.

We refer to our paper for a full description of the results, i.e., hardware, inputs, and competitors.
Below, you can find some of the figures we present in the paper.

![Screenshot Documentation](https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/space_requirements_v1.0.0.png)

![Screenshot Documentation](https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/rank_times_v1.0.0.png)

![Screenshot Documentation](https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/select_times_v1.0.0.png)

![Screenshot Documentation](https://raw.githubusercontent.com/pasta-toolbox/bit_vector/main/docs/images/select_times_pasta_only_v1.0.0.png)

[benchmark]: benchmarks/bit_vector_benchmark.cpp
[rank and select benchmark]: https://github.com/pasta-toolbox/bit_vector_experiments
[tests]: tests/

## How to Get This
Below, we list all commands that are required to build the code in this repository.
To this end, we provide three CMake presets (_debug_, _release_, and _release with debug information_).

- The debug preset creates a `debug` folder and uses the compiler flags `-DDEBUG -O0 -g -ggdb -fsanitize=address`.
- The release preset creates a `build` folder and uses the compiler flags `-DNDEBUG -march=native -O3`.
- The release with debug information preset creates a `build_with_debug_info` folder and uses the compiler flags `-DDEBUG -g -march=native -O3`.

Per default, we use the following compiler flags: `-Wall -Wextra -Wpedantic -fdiagnostics-color=always`.

### Requirements
pasta::bit_vector is written in C++20 and requires a compiler that [supports][] it.
We use [Ninja][] as build system.
For more information on how to use this library, please refer to our [documentation][].

[supports]: https://en.cppreference.com/w/cpp/compiler_support
[Ninja]: https://ninja-build.org/

### tl;dr
To just clone the source code, use the following.
```bash
git clone git@github.com:pasta-toolbox/bit_vector
cd bit_vector
git submodule update --init --recursive
```
If you also want to build the test, please continue with the following commands.
```bash
cmake --preset=[debug|build|relwithdeb]-DPASTA_BIT_VECTOR_BUILD_TESTS=On
cmake --build --preset=[debug|release|relwithdeb]
ctest --test-dir [debug|build|relwithdeb]
```

[documentation]: https://www.pasta-toolbox.org/bit_vector/
