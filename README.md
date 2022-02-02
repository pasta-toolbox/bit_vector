# pasta::bit_vector

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pasta::bit_vector CI](https://github.com/pasta-toolbox/bit_vector/actions/workflows/ctest.yml/badge.svg)](https://github.com/pasta-toolbox/bit_vector/actions/workflows/ctest.yml)
[![codecov](https://codecov.io/gh/pasta-toolbox/bit_vector/branch/main/graph/badge.svg?token=2QD6ME44SU)](https://codecov.io/gh/pasta-toolbox/bit_vector)

This header-only library contains a highly tuned (uncompressed) bit vector implementation with multiple space efficient rank and select support data structures.
The our fastest rank and select support has a space overhead of only ~3% and makes use of data level parallelism via SIMD instructions.

## Contents
This repository contains the following algorithms and data structures.

### Bit Vectors
Bit vectors play an important role in many compressed text indices, e.g., the FM-index.
This repository contains the following bit vector implementations:

- highly tuned [uncompressed bit vector][] with access operator
- compact [rank][] and [select][] support for the uncompressed bit vector based on

> Dong Zhou and David G. Andersen and Michael Kaminsky,
> Space-Efficient, High-Performance Rank and Select Structures on Uncompressed Bit Sequences,
> SEA 2013.

- improved [rank][] and [select][] support requiring the same amount of memory but providing faster rank (most of the time) and select (always) queries, and
- a very fast [rank][] support that can also answer [select][] queries

TODO FIX THE LINKS
TODO ADD BENCHMARK RESULTS
TODO ADD EXAMPLES
TODO REMOVE SDSL
TODO MORE FOCUS ON THE RANK AND SELECT, LESS FOCUS ON THE BIT VECTOR

[uncompressed bit vector]: bit_vector/bit_vector.hpp
[rank]: bit_vectorsupport/bit_vector_rank.hpp
[select]: bit_vector/bit_vector_rank_select.hpp

### Benchmarks and Tests

There exist an easy to use [benchmark][], which helps to compare the implementations in this repository with the ones contained in the [SDSL][].
To build the benchmark, run the CMake command with `-DPASTA_BIT_VECTOR_BUILD_BENCHMARKS=On`.
Our tests are contained in the folder [tests][].
To build the tests, run the CMake command with `-DPASTA_BIT_VECTOR_BUILD_TESTS=On`.

[benchmark]: benchmarks/bit_vector_benchmark.cpp
[SDSL]: https://github.com/simongog/sdsl-lite
[tests]: tests/

## How to Get This
Below, we list all commands that are required to build the code in this repository.
To this end, we provide three CMake presets (_debug_, _release_, and _release with debug information_).
The debug preset creates a `debug` folder and uses the compiler flags `-DDEBUG -O0 -g -ggdb -fsanitize=address`.
The release preset creates a `build` folder and uses the compiler flags `-DNDEBUG -march=native -O3`.
The release with debug information preset creates a `build_with_debug_info` folder and uses the compiler flags `-DDEBUG -g -march=native -O3`.
Per default, we use the following compiler flags: `-Wall -Wextra -Wpedantic -fdiagnostics-color=always`.

### Requirements
pasta::bit_vector is written in C++20 and requires a compiler that [supports][] it.
We use [Ninja][] as build system.

[supports]: https://en.cppreference.com/w/cpp/compiler_support
[Ninja]: https://ninja-build.org/

### tl;dr
To just clone the source code, use the following.
```
git clone git@github.com:pasta-toolbox/bit_vector
cd bit_vector
git submodule update --init --recursive
```
If you also want to build the test, please continue with the following commands.
```
cmake --preset=[debug|build|relwithdeb]-DPASTA_BIT_VECTOR_BUILD_TESTS=On
cmake --build --preset=[debug|release|relwithdeb]
ctest --test-dir [debug|build|build_with_debug_info]
```

