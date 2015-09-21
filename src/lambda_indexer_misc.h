// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2015, Hannes Hauswedell, FU Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
// ==========================================================================
// lambda_indexer_misc.hpp: misc stuff for indexer
// ==========================================================================

#ifndef LAMBDA_INDEXER_MISC_HPP_
#define LAMBDA_INDEXER_MISC_HPP_

// ----------------------------------------------------------------------------
// Class ComparisonCounter
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct ComparisonCounter;

// no counting
template <typename TText>
struct ComparisonCounter<TText, Nothing>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
    uint64_t _lastPercent = 0;
    ComparisonCounter(TText const &,
                      uint64_t expectedComparisons = 0u)
    {
        (void)expectedComparisons;
    }

    // may be constexpr in c++14
    inline void inc() const
    {}
};

// every thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::false_type>
#else
template <typename TText, typename TSpec>
struct ComparisonCounter
#endif
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.2 * double(l) * std::log(l) / std::log(2);
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
        _comparisons = 0;
        _lastPercent = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        uint64_t comp = ++_comparisons;
        // it is not important that the henceforth _comparisons be actually
        // the same value (might not be due to SMP)

        // progress reporting
        if (comp & _checkEveryNHits)
        {
            uint64_t curPerc = comp * 100 / _expectedComparisons;
            if (curPerc < 100)
                printProgressBar(_lastPercent, curPerc);
        }
    }
};

// only one thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::true_type>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.5 * double(l) * std::log(l) / std::log(2) /
                                   omp_get_max_threads();
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
//         _comparisons = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        if (omp_get_thread_num() == 0) // only one thread counts
        {
            uint64_t comp = ++_comparisons;

            // progress reporting
            if (comp & _checkEveryNHits)
            {
                uint64_t curPerc = comp * 100 / _expectedComparisons;
                if (curPerc < 100)
                    printProgressBar(_lastPercent, curPerc);
            }
        }
    }
};
#endif

#endif // LAMBDA_INDEXER_MISC_HPP_