/**
   C++ implementation of the RS algorithm
   Copyright Vera-Licona Research Group (C) 2015
   Author: Andrew Gainer-Dewar, Ph.D. <andrew.gainer.dewar@gmail.com>

   This file is part of MHSGenerationAlgorithms.

   MHSGenerationAlgorithms is free software: you can redistribute it
   and/or modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation, either version 3 of
   the License, or (at your option) any later version.

   MHSGenerationAlgorithms is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.
**/

#ifndef _RS__H
#define _RS__H

#include "hypergraph.hpp"
#include "shd-algorithm.hpp"

#include <atomic>

namespace agdmhs {
    struct RSCounters {
        std::atomic<unsigned> iterations;
        std::atomic<unsigned> violators;
        std::atomic<unsigned> critical_fails;
        std::atomic<unsigned> update_loops;
        std::atomic<unsigned> tasks_waiting;
    };

    class RSAlgorithm: public SHDAlgorithm {
        size_t num_threads;
        size_t cutoff_size;

    public:
        RSAlgorithm (size_t cutoff_size = 0);
        RSAlgorithm (size_t num_threads, size_t cutoff_size = 0);
        Hypergraph transversal (const Hypergraph& H) const override;

    private:
        void extend_or_confirm_set (const Hypergraph& H, const Hypergraph& T, RSCounters& counters, bsqueue& hitting_sets, bitset& S, Hypergraph& crit, bitset& uncov, const bitset& violating_vertices) const;
        static bool any_edge_critical_after_i (const hindex& i, const bitset& S, const Hypergraph& crit);
    };
}

#endif
