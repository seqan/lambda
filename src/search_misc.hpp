// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2020, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2020, Knut Reinert and Freie Universität Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is Free Software: you can redistribute it and/or modify it
// under the terms found in the LICENSE[.md|.rst] file distributed
// together with this file.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// ==========================================================================
// search_misc.h: misc stuff for search
// ==========================================================================

#pragma once

#include <vector>

#include <seqan/blast.h>

// ============================================================================
// Exceptions
// ============================================================================

struct IndexException : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

struct QueryException : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

// ============================================================================
// Alignment-related
// ============================================================================

inline int64_t _bandSize(uint64_t const seqLength)
{
    // Currently only used for padding of the subject sequences
    return static_cast<int64_t>(std::sqrt(seqLength)) + 1;
}

// ----------------------------------------------------------------------------
// Function computeEValueThreadSafe
// ----------------------------------------------------------------------------

template <typename TBlastMatch, typename TScore, seqan::BlastProgram p, seqan::BlastTabularSpec h>
inline double computeEValueThreadSafe(TBlastMatch & match, uint64_t ql, seqan::BlastIOContext<TScore, p, h> & context)
{
#if defined(__FreeBSD__)
    // && version < 11 && defined(STDLIB_LLVM) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=192320
    // || version >= 11 && defined(STDLIB_GNU) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=215709
    static std::vector<std::unordered_map<uint64_t, uint64_t>> _cachedLengthAdjustmentsArray(omp_get_num_threads());
    std::unordered_map<uint64_t, uint64_t> &                   _cachedLengthAdjustments =
      _cachedLengthAdjustmentsArray[omp_get_thread_num()];
#else
    static thread_local std::unordered_map<uint64_t, uint64_t> _cachedLengthAdjustments;
#endif

    // convert to 64bit and divide for translated sequences
    ql = ql / (seqan::qIsTranslated(context.blastProgram) ? 3 : 1);
    // length adjustment not yet computed
    if (_cachedLengthAdjustments.find(ql) == _cachedLengthAdjustments.end())
        _cachedLengthAdjustments[ql] = _lengthAdjustment(context.dbTotalLength, ql, context.scoringScheme);

    uint64_t adj = _cachedLengthAdjustments[ql];

    match.eValue =
      _computeEValue(match.alignStats.alignmentScore, ql - adj, context.dbTotalLength - adj, context.scoringScheme);
    return match.eValue;
}

// ----------------------------------------------------------------------------
// compute LCA
// ----------------------------------------------------------------------------

template <typename T, typename T2>
T computeLCA(std::vector<T> const & taxParents, std::vector<T2> const & taxHeights, T n1, T n2)
{
    if (n1 == n2)
        return n1;

    // move up so that nodes are on same height
    for (auto i = taxHeights[n1]; i > taxHeights[n2]; --i)
        n1 = taxParents[n1];

    for (auto i = taxHeights[n2]; i > taxHeights[n1]; --i)
        n2 = taxParents[n2];

    while ((n1 != 0) && (n2 != 0))
    {
        // common ancestor
        if (n1 == n2)
            return n1;

        // move up in parallel
        n1 = taxParents[n1];
        n2 = taxParents[n2];
    }

    throw std::runtime_error{"LCA-computation error: One of the paths didn't lead to root."};
    return 0; // avoid warnings on clang
}
