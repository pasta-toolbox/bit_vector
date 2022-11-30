/*******************************************************************************
 * benchmarks/bit_vector/bit_vector_benchmark.cpp
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
#include <pasta/bit_vector/compression/block_compressed_bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#if defined(DNDEBUG)
#  include <pasta/utils/memory_monitor.hpp>
#endif
#include <cstdint>
#include <fstream>
#include <iostream>
#include <pasta/utils/benchmark/timer.hpp>
#include <queue>
#include <random>
#include <string>
#include <string_view>
#include <tlx/cmdline_parser.hpp>
#include <tlx/die.hpp>
#include <tlx/logger.hpp>
#include <tlx/math.hpp>

template <typename SplitWordType>
float simple_hist(std::span<uint64_t> bv_data,
                  size_t const m,
                  [[maybe_unused]] size_t const n,
                  [[maybe_unused]] size_t const a) {
  std::unordered_map<uint64_t, size_t> hist;

  for (size_t i = 0; i < bv_data.size(); ++i) {
    uint64_t cur_val = bv_data[i];
    for (size_t j = 0; j < 8 / sizeof(SplitWordType); ++j) {
      uint64_t small_val =
          (cur_val & std::numeric_limits<SplitWordType>::max());
      if constexpr (sizeof(SplitWordType) < sizeof(decltype(cur_val))) {
        cur_val >>= 8 * sizeof(SplitWordType);
      }
      if (auto res = hist.find(small_val); res != hist.end()) {
        ++(res->second);
      } else {
        hist.insert({small_val, 1});
      }
    }
  }

  std::vector<size_t> popcounts(8 * sizeof(SplitWordType) + 1, 0);

  for (auto const& block : hist) {
    ++popcounts[std::popcount(block.first)];
  }

  size_t compressed_hist_bit_size = 0;
  for (size_t i = 0; i < popcounts.size(); ++i) {
    if (popcounts[i] > 0) {
      compressed_hist_bit_size +=
          popcounts[i] * i * tlx::integer_log2_ceil(8 * sizeof(SplitWordType));
      // std::cout << i << " " << popcounts[i] << std::endl;
    }
  }
  size_t const bits_per_code_word = tlx::integer_log2_ceil(hist.size());
  size_t const number_codes = bv_data.size() * (8 / sizeof(SplitWordType));

  double const bits_per_block = double(bits_per_code_word * number_codes) / m;

  size_t const dict_size_in_bits = hist.size() * sizeof(SplitWordType) * 8;

  return bits_per_block + ((float(n) * 0.0351) / m) +
         float(dict_size_in_bits) / m;
}

template <typename SplitWordType>
float shifted_hist(std::span<uint64_t> bv_data,
                   size_t const m,
                   [[maybe_unused]] size_t const n,
                   [[maybe_unused]] size_t const a) {
  std::unordered_map<uint64_t, size_t> hist;

  for (size_t i = 0; i < bv_data.size(); ++i) {
    uint64_t cur_val = bv_data[i];
    for (size_t j = 0; j < 8 / sizeof(SplitWordType); ++j) {
      uint64_t small_val =
          (cur_val & std::numeric_limits<SplitWordType>::max());
      if constexpr (sizeof(SplitWordType) < sizeof(decltype(cur_val))) {
        cur_val >>= 8 * sizeof(SplitWordType);
      }
      bool inserted = false;
      for (size_t k = 0; k < (8 * sizeof(SplitWordType)) - 1; ++k) {
        if (auto res = hist.find(small_val >> k); res != hist.end()) {
          ++(res->second);
          inserted = true;
          break;
        }
      }
      if (!inserted) {
        hist.insert({small_val, 1});
      }
    }
  }

  size_t const bits_per_code_word = tlx::integer_log2_ceil(hist.size());
  size_t const number_codes = bv_data.size() * (8 / sizeof(SplitWordType));

  double const bits_per_block = double(bits_per_code_word * number_codes) / m;
  size_t const shift_bits = tlx::integer_log2_ceil(8 * sizeof(SplitWordType));

  size_t const dict_size_in_bits = hist.size() * sizeof(SplitWordType) * 8;

  return bits_per_block + shift_bits + ((float(n) * 0.0351) / m) +
         float(dict_size_in_bits) / m;
}

struct FrequencyItem {
  uint64_t frequency;
  std::vector<uint64_t> words;

  bool operator<(FrequencyItem const& other) {
    return frequency < other.frequency;
  }

  bool operator>(FrequencyItem const& other) {
    return frequency > other.frequency;
  }
}; // struct CodePair

bool operator>(FrequencyItem const& a, FrequencyItem const& b) {
  return a.frequency > b.frequency;
}

struct DetailedHuffSizes {
  double bv_per_block;
  double dictionary_per_block;
  double rs_per_block;

  double total_size() const {
    return bv_per_block + dictionary_per_block + rs_per_block;
  }
}; // struct DetailedHuffSizes

template <typename SplitWordType>
DetailedHuffSizes huffman_hist(std::span<uint64_t> bv_data,
                               size_t const m,
                               [[maybe_unused]] size_t const n,
                               [[maybe_unused]] size_t const a) {
  std::unordered_map<uint64_t, size_t> hist;

  for (size_t i = 0; i < bv_data.size(); ++i) {
    uint64_t cur_val = bv_data[i];
    for (size_t j = 0; j < 8 / sizeof(SplitWordType); ++j) {
      uint64_t small_val =
          (cur_val & std::numeric_limits<SplitWordType>::max());
      if constexpr (sizeof(SplitWordType) < sizeof(decltype(cur_val))) {
        cur_val >>= 8 * sizeof(SplitWordType);
      }
      if (auto res = hist.find(small_val); res != hist.end()) {
        ++(res->second);
      } else {
        hist.insert({small_val, 1});
      }
    }
  }

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

  size_t compressed_bv_size = 0;
  for (auto const& block : hist) {
    compressed_bv_size += block.second * code_lengths[block.first];
  }

  double const bits_per_block = double(compressed_bv_size) / m;

  size_t const dict_size_in_bits = hist.size() * sizeof(SplitWordType) * 8;

  return DetailedHuffSizes{bits_per_block,
                           float(dict_size_in_bits) / m,
                           (float(n) * 0.0475) / m};
}

DetailedHuffSizes huff_dist([[maybe_unused]] std::string path,
                            std::span<uint64_t> bv_data,
                            size_t const m,
                            [[maybe_unused]] size_t const n,
                            [[maybe_unused]] size_t const a) {
  std::unordered_map<size_t, size_t> distances;

  // std::ofstream outfile;
  //  outfile.open(path + std::string(".hist.csv"));

  size_t last_distance = 0;
  size_t max_distance = 0;
  for (size_t i = 0; i < bv_data.size(); ++i) {
    uint64_t cur_val = bv_data[i];
    for (size_t j = 0; j < 64; ++j) {
      while (j < 64 && !((cur_val >> j) & 1ULL)) {
        ++j;
        ++last_distance;
      }
      if (j != 64) {
        // std::cout << text << " " << last_distance << " " << (i * 64) + j <<
        // "\n"; outfile << last_distance << "\n";

        if (auto res = distances.find(last_distance); res != distances.end()) {
          ++(res->second);
        } else {
          distances.insert({last_distance, 1});
          if (last_distance > max_distance) {
            max_distance = last_distance;
          }
        }
        last_distance = 0;
      }
    }
  }

  std::vector<std::pair<size_t, size_t>> tmp;
  for (auto const& block : distances) {
    tmp.emplace_back(block.first, block.second);
  }
  std::sort(tmp.begin(),
            tmp.end(),
            [](std::pair<size_t, size_t> const& lhs,
               std::pair<size_t, size_t> const& rhs) {
              return lhs.first < rhs.first;
            });

  // std::ofstream outfile;
  // outfile.open(path + std::string(".hist.csv"));

  // outfile << "length" << " " << "occ" << "\n";

  size_t max_length = tmp.back().first;

  if (max_length < 256) {
    std::ofstream outfile;
    outfile.open(path + std::string(".entropy-text"));

    last_distance = 0;
    max_distance = 0;
    for (size_t i = 0; i < bv_data.size(); ++i) {
      uint64_t cur_val = bv_data[i];
      for (size_t j = 0; j < 64; ++j) {
        while (j < 64 && !((cur_val >> j) & 1ULL)) {
          ++j;
          ++last_distance;
        }
        if (j != 64) {
          // std::cout << text << " " << last_distance << " " << (i * 64) + j <<
          // "\n";
          outfile << uint8_t(last_distance);

          if (auto res = distances.find(last_distance);
              res != distances.end()) {
            ++(res->second);
          } else {
            distances.insert({last_distance, 1});
            if (last_distance > max_distance) {
              max_distance = last_distance;
            }
          }
          last_distance = 0;
        }
      }
    }
    outfile.close();
  }

  // size_t cur_length = 0;
  // for (size_t i = 0; i < tmp.size(); ++i) {
  //   while (cur_length != max_length && cur_length++ != tmp[i].first) {
  //     outfile << (cur_length - 1) << " " << 0 << "\n";
  //   }
  //   outfile << tmp[i].first << " " << tmp[i].second << "\n";
  // }
  // outfile.close();

  // for (auto const& e : tmp) {
  //   std::cout << e.first << " " << e.second << "\n";
  // }

  // Sort single symbols by number ob occurrence
  std::priority_queue<FrequencyItem,
                      std::vector<FrequencyItem>,
                      std::greater<FrequencyItem>>
      frequency_tree;

  std::unordered_map<uint64_t, size_t> code_lengths;
  for (auto const& block : distances) {
    frequency_tree.emplace(block.second, std::vector<size_t>{block.first});
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

  size_t compressed_bv_size = 0;
  for (auto const& block : distances) {
    compressed_bv_size += block.second * code_lengths[block.first];
  }

  double const bits_per_block = double(compressed_bv_size) / m;

  size_t dict_entry_size = 1;

  if (max_distance > std::numeric_limits<uint8_t>::max()) {
    dict_entry_size = 2;
  }
  if (max_distance > std::numeric_limits<uint16_t>::max()) {
    dict_entry_size = 4;
  }
  if (max_distance > std::numeric_limits<uint32_t>::max()) {
    dict_entry_size = 8;
  }

  size_t const dict_size_in_bits = distances.size() * dict_entry_size * 8;

  return DetailedHuffSizes{bits_per_block,
                           float(dict_size_in_bits) / m,
                           (float(n) * 0.0475) / m};
}

int32_t main(int32_t argc, char const* const argv[]) {
  std::string input_path = "";
  tlx::CmdlineParser cp;
  cp.set_description("Benchmark tool for PaStA's bit vector implementation.");
  cp.set_author("Florian Kurpicz <florian@kurpicz.org>");

  cp.add_param_string("input", input_path, "Path to input file.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  std::ifstream fileIn(input_path, std::ios::in | std::ios::binary);
  size_t size = 0;
  fileIn.read(reinterpret_cast<char*>(&size), sizeof(size_t));
  pasta::BitVector bitVector(size);
  pasta::BitVector bitVectorToCompress(size);
  unsigned long* inputBitData = bitVector.data().data();
  fileIn.read(reinterpret_cast<char*>(inputBitData),
              bitVector.data().size_bytes());
  fileIn.seekg(0, std::ios_base::beg);
  fileIn.read(reinterpret_cast<char*>(&size), sizeof(size_t));
  unsigned long* inputToCompressBitData = bitVectorToCompress.data().data();
  fileIn.read(reinterpret_cast<char*>(inputToCompressBitData),
              bitVectorToCompress.data().size_bytes());
  fileIn.close();

  std::cout << "bitVectorToCompress.data().data()[0] "
            << bitVectorToCompress.data().data()[0] << '\n';
  std::cout << "bitVector.data().data()[0] " << bitVector.data().data()[0]
            << '\n';

  pasta::BlockCompressedBitVector bcbv(std::move(bitVectorToCompress));
  bcbv.compress();

  std::cout << "bitVector.size() " << bitVector.size() << '\n';
  [[maybe_unused]] pasta::FlatRankSelect rs(bcbv);
  // pasta::FlatRankSelect ucrs(bitVector);

  // std::cout << "bitVector.size() " << bitVector.size() << '\n';

  // for (size_t i = 0; i < bitVector.size(); ++i) {
  //   auto const comp_rank = rs.rank1(i);
  //   auto const ucommp_rank = ucrs.rank1(i);
  //   if (comp_rank != ucommp_rank) {
  //     std::cout << "i " << i << '\n';
  //     std::cout << "bitVector.data().data()[i/64] " <<
  //     bitVector.data().data()[i/64] << '\n'; std::cout <<
  //     "bitVector.data().data()[(i/64)+1] " <<
  //     bitVector.data().data()[(i/64)+1] << '\n'; std::cout << "ist: " <<
  //     ucommp_rank << " soll: " << comp_rank << "\n";
  //   }
  // }

  // size_t number_of_ones = rs.rank1(bitVector.size() - 1);

  // std::cout << "number_of_ones " << number_of_ones << '\n';

  // for (size_t i = 1; i < number_of_ones; ++i) {
  //   auto const comp_select = rs.select1(i);
  //   auto const ucommp_select = ucrs.select1(i);
  //   if (comp_select != ucommp_select) {
  //     std::cout << i << "\n";
  //     std::cout << "##WARNING" << "\n";
  //   }
  // }

  // size_t number_of_zeros = rs.rank0(bitVector.size() - 1);
  // std::cout << "number_of_zeros " << number_of_zeros << "\n";

  // for (size_t i = 1; i < number_of_ones; ++i) {
  //   auto const comp_select = rs.select0(i);
  //   auto const ucommp_select = ucrs.select0(i);
  //   if (comp_select != ucommp_select) {
  //     std::cout << i << "\n";
  //     std::cout << "##WARNING" << "\n";
  //   }
  // }

  // pasta::FlatRankSelect rs(bitVector);
  // size_t const m = rs.rank1(bitVector.size());
  // size_t const n = rs.rank0(bitVector.size());
  // size_t const a = n / m;

  // input_path = input_path.substr(input_path.find_last_of("/\\") + 1);

  // std::cout << input_path << " & ";// << a << " & ";

  // float ef = 2 + tlx::integer_log2_ceil(a);

  // std::cout << ef << " & ";

  // float succincter = a * std::log2((float(a) + 1) / float(a)) +
  // tlx::integer_log2_ceil(a) + std::log2(1 + (1/float(a)));

  // std::cout << succincter << " & ";

  // std::vector<DetailedHuffSizes> sizes(4);
  // sizes[0] = huffman_hist<uint8_t>(bitVector.data(), m, n, a);
  // sizes[1] = huffman_hist<uint16_t>(bitVector.data(), m, n, a);
  // sizes[2] = huffman_hist<uint32_t>(bitVector.data(), m, n, a);
  // sizes[3] = huffman_hist<uint64_t>(bitVector.data(), m, n, a);

  // size_t idx = 0;
  // DetailedHuffSizes details = sizes[idx];
  // float min = details.total_size();
  // for (size_t i = 1; i < sizes.size(); ++i) {
  //   if (sizes[i].total_size() < min) {
  //     details = sizes[i];
  //     min = details.total_size();
  //     idx = i;
  //   }
  // }

  // [[maybe_unused]] DetailedHuffSizes dist = huff_dist(input_path,
  // bitVector.data(), m, n, a);

  // std::string pre = "";
  // std::string post = "";

  // pre = min > ef ? "\\textcolor{red}{" : pre;
  // post = min > ef ? "}" : post;

  // if (dist.total_size() > min) {
  //   pre = min < succincter ? "\\textbf{" : "";
  //   post = min < succincter ? "}" : "";
  // }

  // std::cout << pre << min << post << " & (" << std::pow(2, idx) * 8 << ") [ "
  // << sizes[1].total_size() << " (16)] & ";

  // pre = "";
  // post = "";

  // if (dist.total_size() < min) {
  //   pre = dist.total_size() < succincter ? "\\textbf{" : "";
  //   post = dist.total_size() < succincter ? "}" : "";
  // }

  // pre = dist.total_size() > ef ? "\\textcolor{red}{" : pre;
  // post = dist.total_size() > ef ? "}" : post;

  // std::cout << pre << dist.total_size() << post << "\\\\\n";

  return 0;
}

/******************************************************************************/
