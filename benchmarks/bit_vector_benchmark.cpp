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

#include <cstdint>
#include <iostream>
#include <random>

#include <tlx/math/aggregate.hpp>
#include <tlx/cmdline_parser.hpp>
#include <tlx/die.hpp>
#include <tlx/logger.hpp>

#include <sdsl/bit_vectors.hpp>

#include "bit_vector/bit_vector.hpp"
#include "bit_vector/support/bit_vector_rank.hpp"
#include "bit_vector/support/bit_vector_rank_select.hpp"
#include "bit_vector/support/flattened_rank_support.hpp"
#include "utils/do_not_optimize.hpp"
// #include "util/git_commit.hpp"
#include "utils/memory_monitor.hpp"
#include "utils/timer.hpp"

class BitVectorBenchmark {

  static constexpr bool debug = true;
  static constexpr auto LOG_PREFIX = "[BitVectorBenchmark] ";

public:

  void run() {

    die_verbose_unless(fill_percentage_ <= 100, "-f [--fill_percentage] must "
		       "be between 0 and 100 inclusive.");

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, bit_size_ - 1);

    {
      LOG << LOG_PREFIX << "Creating PaStA bit vector";

      pasta::Timer timer;
      pasta::MemoryMonitor& mem_monitor = pasta::MemoryMonitor::instance();
      mem_monitor.reset();
      pasta::BitVector bv(bit_size_, 0);

      size_t const bv_construction_time = timer.get_and_reset();
      auto const bv_construction_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Flipping bits with uniform distribution";
      for (size_t i = 0; i < bit_size_; ++i) {
	bv[i] = (static_cast<uint32_t>(rand() % 100) < fill_percentage_);
      }

      size_t const bv_set_bits_time = timer.get_and_reset();
      auto const bv_set_bits_mem = mem_monitor.get_and_reset();

      pasta::BitVectorRankSelect bvrs(bv);

      size_t const rs_construction_time = timer.get_and_reset();
      LOG << LOG_PREFIX << "Preparing queries";
      timer.reset();
      auto const rs_construction_mem = mem_monitor.get_and_reset();

      std::uniform_int_distribution<> rank_dist(0, bit_size_ - 1);
      std::vector<size_t> rank_positions(query_count_);

      tlx::Aggregate<size_t> rank_query_properties;
      for (auto& pos : rank_positions) {
	pos = rank_dist(gen);
	rank_query_properties.add(pos);
      }

      size_t const zero_bits = bvrs.rank0(bit_size_);
      size_t const one_bits = bvrs.rank1(bit_size_);

      std::vector<size_t> select0_positions(query_count_ / 2);
      std::vector<size_t> select1_positions(query_count_ / 2 +
					    ((query_count_ % 2 == 0) ? 0 : 1));
      std::uniform_int_distribution<> select0_dist(1, zero_bits);
      std::uniform_int_distribution<> select1_dist(1, one_bits);

      tlx::Aggregate<size_t> select0_query_properties;
      tlx::Aggregate<size_t> select1_query_properties;
      for (auto& pos : select0_positions) {
	pos = select0_dist(gen);
	select0_query_properties.add(pos);
      }
      for (auto& pos : select1_positions) {
	pos = select1_dist(gen);
	select1_query_properties.add(pos);
      }

      LOG << LOG_PREFIX << "Benchmarking queries";
      timer.reset();
      mem_monitor.reset();

      for (size_t i = 0; i < rank_positions.size() / 2; ++i) {
	[[maybe_unused]]
	size_t const result = bvrs.rank0(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }
      for (size_t i = rank_positions.size() / 2; i < rank_positions.size();
	   ++i) {
	[[maybe_unused]]
	size_t const result = bvrs.rank1(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }

      size_t const rank_query_time = timer.get_and_reset();

      for (auto const pos : select0_positions) {
	[[maybe_unused]]
	size_t const result = bvrs.select0(pos);
	PASTA_DO_NOT_OPTIMIZE(result);
      }
      for (auto const pos : select1_positions) {
	[[maybe_unused]]
	size_t const result = bvrs.select1(pos);
	PASTA_DO_NOT_OPTIMIZE(result);
      }

      size_t const select_query_time = timer.get_and_reset();
      auto const rs_query_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Query stats";
      LOG << LOG_PREFIX << "Rank positions min/max/avg:"
	  << rank_query_properties.min() << " / "
	  << rank_query_properties.max() << " / "
	  << rank_query_properties.avg();

      LOG << LOG_PREFIX << "Select0 rank min/max/avg:"
	  << select0_query_properties.min() << " / "
	  << select0_query_properties.max() << " / "
	  << select0_query_properties.avg();

      LOG << LOG_PREFIX << "Select1 rank min/max/avg:"
	  << select1_query_properties.min() << " / "
	  << select1_query_properties.max() << " / "
	  << select1_query_properties.avg();

      LOG << LOG_PREFIX << "Finished PaStA bit vector benchmark";

      std::cout << "RESULT "
                << "algo=popcount_bv_uncompressed "
                << "bit_size=" << bit_size_ << " "
                << "fill_percentage=" << fill_percentage_ << " "
                << "bv_construction_time=" << bv_construction_time << " "
                << "bv_construction_mem=" << bv_construction_mem.cur_peak << " "
                << "bv_set_bits_time=" << bv_set_bits_time << " "
                << "bv_set_bits_mem=" << bv_set_bits_mem.cur_peak << " "
                << "rs_construction_time=" << rs_construction_time << " "
                << "rs_construction_mem=" << rs_construction_mem.cur_peak << " "
                << "query_count=" << query_count_ << " "
                << "rank_query_time=" << rank_query_time << " "
                << "select_query_time=" << select_query_time << " "
                << "rs_query_mem=" << rs_query_mem.cur_peak << " "
                << "\n";
    }

    {
      LOG << LOG_PREFIX << "Creating PaStA bit vector";

      pasta::Timer timer;
      pasta::MemoryMonitor& mem_monitor = pasta::MemoryMonitor::instance();
      mem_monitor.reset();
      pasta::BitVector bv(bit_size_, 0);

      size_t const bv_construction_time = timer.get_and_reset();
      auto const bv_construction_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Flipping bits with uniform distribution";
      for (size_t i = 0; i < bit_size_; ++i) {
	bv[i] = (static_cast<uint32_t>(rand() % 100) < fill_percentage_);
      }

      size_t const bv_set_bits_time = timer.get_and_reset();
      auto const bv_set_bits_mem = mem_monitor.get_and_reset();

      pasta::FlattenedBitVectorRank bvrs(bv);

      size_t const rs_construction_time = timer.get_and_reset();
      LOG << LOG_PREFIX << "Preparing queries";
      timer.reset();
      auto const rs_construction_mem = mem_monitor.get_and_reset();

      std::uniform_int_distribution<> rank_dist(0, bit_size_ - 1);
      std::vector<size_t> rank_positions(query_count_);

      tlx::Aggregate<size_t> rank_query_properties;
      for (auto& pos : rank_positions) {
	pos = rank_dist(gen);
	rank_query_properties.add(pos);
      }

      size_t const zero_bits = bvrs.rank0(bit_size_);
      size_t const one_bits = bvrs.rank1(bit_size_);

      std::vector<size_t> select0_positions(query_count_ / 2);
      std::vector<size_t> select1_positions(query_count_ / 2 +
					    ((query_count_ % 2 == 0) ? 0 : 1));
      std::uniform_int_distribution<> select0_dist(1, zero_bits);
      std::uniform_int_distribution<> select1_dist(1, one_bits);

      tlx::Aggregate<size_t> select0_query_properties;
      tlx::Aggregate<size_t> select1_query_properties;
      for (auto& pos : select0_positions) {
	pos = select0_dist(gen);
	select0_query_properties.add(pos);
      }
      for (auto& pos : select1_positions) {
	pos = select1_dist(gen);
	select1_query_properties.add(pos);
      }

      LOG << LOG_PREFIX << "Benchmarking queries";
      timer.reset();
      mem_monitor.reset();

      for (size_t i = 0; i < rank_positions.size() / 2; ++i) {
	[[maybe_unused]]
	size_t const result = bvrs.rank0(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }
      for (size_t i = rank_positions.size() / 2; i < rank_positions.size();
	   ++i) {
	[[maybe_unused]]
	size_t const result = bvrs.rank1(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }

      size_t const rank_query_time = timer.get_and_reset();

      // for (auto const pos : select0_positions) {
      // 	[[maybe_unused]]
      // 	size_t const result = bvrs.select0(pos);
      // 	PASTA_DO_NOT_OPTIMIZE(result);
      // }
      // for (auto const pos : select1_positions) {
      // 	[[maybe_unused]]
      // 	size_t const result = bvrs.select1(pos);
      // 	PASTA_DO_NOT_OPTIMIZE(result);
      // }

      size_t const select_query_time = timer.get_and_reset();
      auto const rs_query_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Query stats";
      LOG << LOG_PREFIX << "Rank positions min/max/avg:"
	  << rank_query_properties.min() << " / "
	  << rank_query_properties.max() << " / "
	  << rank_query_properties.avg();

      LOG << LOG_PREFIX << "Select0 rank min/max/avg:"
	  << select0_query_properties.min() << " / "
	  << select0_query_properties.max() << " / "
	  << select0_query_properties.avg();

      LOG << LOG_PREFIX << "Select1 rank min/max/avg:"
	  << select1_query_properties.min() << " / "
	  << select1_query_properties.max() << " / "
	  << select1_query_properties.avg();

      LOG << LOG_PREFIX << "Finished PaStA bit vector benchmark";

      std::cout << "RESULT "
                << "algo=flat_popcount_bv_uncompressed "
                << "bit_size=" << bit_size_ << " "
                << "fill_percentage=" << fill_percentage_ << " "
                << "bv_construction_time=" << bv_construction_time << " "
                << "bv_construction_mem=" << bv_construction_mem.cur_peak << " "
                << "bv_set_bits_time=" << bv_set_bits_time << " "
                << "bv_set_bits_mem=" << bv_set_bits_mem.cur_peak << " "
                << "rs_construction_time=" << rs_construction_time << " "
                << "rs_construction_mem=" << rs_construction_mem.cur_peak << " "
                << "query_count=" << query_count_ << " "
                << "rank_query_time=" << rank_query_time << " "
                << "select_query_time=" << select_query_time << " "
                << "rs_query_mem=" << rs_query_mem.cur_peak << " "
                << "\n";
    }

    {
      LOG << LOG_PREFIX << "Creating SDLS bit vector";

      pasta::Timer timer;
      pasta::MemoryMonitor& mem_monitor = pasta::MemoryMonitor::instance();
      mem_monitor.reset();
      sdsl::bit_vector bv(bit_size_, 0);

      size_t const bv_construction_time = timer.get_and_reset();
      auto const bv_construction_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Flipping bits with uniform distribution";
      for (size_t i = 0; i < bit_size_; ++i) {
	bv[i] = (static_cast<uint32_t>(rand() % 100) < fill_percentage_);
      }

      size_t const bv_set_bits_time = timer.get_and_reset();
      auto const bv_set_bits_mem = mem_monitor.get_and_reset();

      sdsl::bit_vector::select_0_type bvs0(&bv);
      sdsl::bit_vector::select_1_type bvs1(&bv);

      sdsl::bit_vector::rank_0_type bvr0(&bv);
      sdsl::bit_vector::rank_1_type bvr1(&bv);

      size_t const rs_construction_time = timer.get_and_reset();
      LOG << LOG_PREFIX << "Preparing queries";
      timer.reset();
      auto const rs_construction_mem = mem_monitor.get_and_reset();

      std::uniform_int_distribution<> rank_dist(0, bit_size_ - 1);
      std::vector<size_t> rank_positions(query_count_);

      tlx::Aggregate<size_t> rank_query_properties;
      for (auto& pos : rank_positions) {
	pos = rank_dist(gen);
	rank_query_properties.add(pos);
      }

      size_t const zero_bits = bvr0.rank(bit_size_);
      size_t const one_bits = bvr1.rank(bit_size_);

      std::vector<size_t> select0_positions(query_count_ / 2);
      std::vector<size_t> select1_positions(query_count_ / 2 +
					    ((query_count_ % 2 == 0) ? 0 : 1));
      std::uniform_int_distribution<> select0_dist(1, zero_bits);
      std::uniform_int_distribution<> select1_dist(1, one_bits);

      tlx::Aggregate<size_t> select0_query_properties;
      tlx::Aggregate<size_t> select1_query_properties;
      for (auto& pos : select0_positions) {
	pos = select0_dist(gen);
	select0_query_properties.add(pos);
      }
      for (auto& pos : select1_positions) {
	pos = select1_dist(gen);
	select1_query_properties.add(pos);
      }

      LOG << LOG_PREFIX << "Benchmarking queries";
      timer.reset();
      mem_monitor.reset();

      for (size_t i = 0; i < rank_positions.size() / 2; ++i) {
	[[maybe_unused]]
	size_t const result = bvr0.rank(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }
      for (size_t i = rank_positions.size() / 2; i < rank_positions.size();
	   ++i) {
	[[maybe_unused]]
	size_t const result = bvr1.rank(rank_positions[i]);
	PASTA_DO_NOT_OPTIMIZE(result);
      }

      size_t const rank_query_time = timer.get_and_reset();

      for (auto const pos : select0_positions) {
	[[maybe_unused]]
	size_t const result = bvs0.select(pos);
	PASTA_DO_NOT_OPTIMIZE(result);
      }
      for (auto const pos : select1_positions) {
	[[maybe_unused]]
	size_t const result = bvs1.select(pos);
	PASTA_DO_NOT_OPTIMIZE(result);
      }

      size_t const select_query_time = timer.get_and_reset();
      auto const rs_query_mem = mem_monitor.get_and_reset();

      LOG << LOG_PREFIX << "Query stats";
      LOG << LOG_PREFIX << "Rank positions min/max/avg:"
	  << rank_query_properties.min() << " / "
	  << rank_query_properties.max() << " / "
	  << rank_query_properties.avg();

      LOG << LOG_PREFIX << "Select0 rank min/max/avg:"
	  << select0_query_properties.min() << " / "
	  << select0_query_properties.max() << " / "
	  << select0_query_properties.avg();

      LOG << LOG_PREFIX << "Select1 rank min/max/avg:"
	  << select1_query_properties.min() << " / "
	  << select1_query_properties.max() << " / "
	  << select1_query_properties.avg();

      LOG << LOG_PREFIX << "Finished SDSL bit vector benchmark";

      std::cout << "RESULT "
                << "algo=sdsl_bv_uncompressed "
                << "bit_size=" << bit_size_ << " "
                << "fill_percentage=" << fill_percentage_ << " "
                << "bv_construction_time=" << bv_construction_time << " "
                << "bv_construction_mem=" << bv_construction_mem.cur_peak << " "
                << "bv_set_bits_time=" << bv_set_bits_time << " "
                << "bv_set_bits_mem=" << bv_set_bits_mem.cur_peak << " "
                << "rs_construction_time=" << rs_construction_time << " "
                << "rs_construction_mem=" << rs_construction_mem.cur_peak << " "
                << "query_count=" << query_count_ << " "
                << "rank_query_time=" << rank_query_time << " "
                << "select_query_time=" << select_query_time << " "
                << "rs_query_mem=" << rs_query_mem.cur_peak << " "
                << "\n";
    }
  }

  size_t bit_size_ = 1024*1024;
  uint32_t fill_percentage_ = 50;
  size_t query_count_ = 10000;

}; // class BitVectorBenchmark

int32_t main(int32_t argc, char const * const argv[]) {
  BitVectorBenchmark bvb;

  tlx::CmdlineParser cp;

  cp.set_description("Benchmark tool for PaStA's bit vector implementation.");
  cp.set_author("Florian Kurpicz <florian@kurpicz.org>");

  cp.add_bytes('b', "bit_size", bvb.bit_size_, "Size of the bit vector in bits "
	       "(accepts SI units, default 1024^2.");

  cp.add_uint('f', "fill_percentage", bvb.fill_percentage_, "Percentage of set "
	      "bits in the bit vector (default 50%).");

  cp.add_bytes('q', "query_count", bvb.query_count_, "Number of rank and select"
	       " queries (accepts SI units, default is 10000)");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  bvb.run();

  return 0;
}

/******************************************************************************/
