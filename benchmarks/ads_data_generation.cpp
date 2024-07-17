/*******************************************************************************
 * benchmarks/bit_vector/asd_data_generation.cpp
 *
 * Copyright (C) 2024 Florian Kurpicz <florian@kurpicz.org>
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

#include <cstdint>
#include <iostream>
#include <fstream>
#include <random>
#include <tlx/cmdline_parser.hpp>

class BenchmarkDataGenerator {

private:
  size_t one_bits_ = 0;

public:

  size_t bit_size_ = 1024 * 1024;
  uint32_t fill_percentage_ = 50;
  size_t query_count_ = 10000;

  void run(std::string name) {
    std::ofstream out_stream(name.c_str());

    std::random_device rd;
    std::mt19937 const gen(rd());

    generate_number_queries(out_stream);
    generate_bit_vector(gen, out_stream);
    generate_queries(gen, out_stream);
    
    out_stream.close();
  }
  
private:
  template <typename OutFunction>
  void generate_number_queries(OutFunction& out) {
    out << query_count_ << '\n';
  }

  template <typename OutFunction>
  void generate_queries(std::mt19937 randomness, OutFunction& out) {
    std::uniform_int_distribution<> query_dist(0, 2);
    std::uniform_int_distribution<> bit_dist(0, 1);
    std::uniform_int_distribution<> pos_dist(0, bit_size_);
    std::uniform_int_distribution<> one_pos_dist(0, one_bits_);
    std::uniform_int_distribution<> zero_pos_dist(0, bit_size_ - one_bits_);

    for (size_t i = 0; i < query_count_; ++i) {    
      if (auto rand = query_dist(randomness); rand == 0) {
        // access query
        out << "access " << pos_dist(randomness) << '\n';
      } else if (rand == 1) {
        // rank query
        out << "rank " << bit_dist(randomness) << " " << pos_dist(randomness) << '\n';
      } else {
        // select query
        out << "select ";
        if (bit_dist(randomness) == 0) {
          out << "0 " << zero_pos_dist(randomness) << '\n';
        } else {
          out << "1 " << one_pos_dist(randomness) << '\n';
        }
      }
    }
  }

  template <typename OutFunction>
  void generate_bit_vector(std::mt19937 randomness, OutFunction& out) {

    std::uniform_int_distribution<> bit_dist(0, 99);
    for (size_t i = 0; i < bit_size_; ++i) {
      if (static_cast<uint32_t>(bit_dist(randomness)) < fill_percentage_) {
        ++one_bits_;
        out << "1";
      } else {
        out << "0";
      }
    }
    out << '\n';
  }

};
  
int32_t main(int32_t argc, char const* const argv[]) {
  tlx::CmdlineParser cp;
  BenchmarkDataGenerator bdg;

  cp.set_description("Benchmark tool for PaStA's bit vector implementation.");
  cp.set_author("Florian Kurpicz <florian@kurpicz.org>");

  cp.add_bytes('b',
               "bit_size",
               bdg.bit_size_,
               "Size of the bit vector in bits "
               "(accepts SI units, default 1024^2.");

  cp.add_uint('f',
              "fill_percentage",
              bdg.fill_percentage_,
              "Percentage of set "
              "bits in the bit vector (default 50%).");

  cp.add_bytes('q',
               "query_count",
               bdg.query_count_,
               "Number of rank and select"
               " queries (accepts SI units, default is 10000)");

  std::string file_name;
  cp.add_string('n',
                "file_name",
                file_name,
                "Name of the benchmark file.");
  
  if (!cp.process(argc, argv)) {
    return -1;
  }

  bdg.run(file_name);

  return 0;
}

/******************************************************************************/
