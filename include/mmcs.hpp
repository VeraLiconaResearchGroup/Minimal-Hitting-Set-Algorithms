/**
   C++ implementation of the MMCS algorithm
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

#ifndef _MMCS__H
#define _MMCS__H

#include "hypergraph.hpp"
#include "shd-algorithm.hpp"

#include <atomic>

namespace agdmhs {
    struct MMCSCounters {
        std::atomic<unsigned> iterations;
        std::atomic<unsigned> violators;
        std::atomic<unsigned> update_loops;
        std::atomic<unsigned> tasks_waiting;
    };

    class MMCSAlgorithm: public SHDAlgorithm {
        size_t num_threads;
        size_t cutoff_size;

    public:
        MMCSAlgorithm (size_t num_threads, size_t cutoff_size);
        Hypergraph transversal (const Hypergraph& H) const override;

    private:
        void extend_or_confirm_set (const Hypergraph& H, const Hypergraph& T, MMCSCounters& counters, bsqueue& hitting_sets, bitset& S, bitset& CAND, Hypergraph& crit, bitset& uncov) const;
    };
}

#endif
