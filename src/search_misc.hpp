// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2017, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2017, Knut Reinert and Freie Universität Berlin
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
// match.h: Main File for the match class
// ==========================================================================

#ifndef LAMBDA_SEARCH_MISC_H_
#define LAMBDA_SEARCH_MISC_H_

#include <vector>

using namespace seqan;

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
// Seeding related
// ============================================================================

template <typename TGH, typename TAlph>
inline void
myHyperSortSingleIndex(std::vector<Match<TAlph>> & matches,
                       bool const doubleIndexing,
                       TGH const &)
{
    using TId = typename Match<TAlph>::TQId;

    // regular sort
    std::sort(matches.begin(), matches.end());

    //                    trueQryId, begin,    end
    std::vector<std::tuple<TId, TId, TId>> intervals;
    for (TId i = 1; i <= length(matches); ++i)
    {
        if ((i == length(matches))                      ||
            (matches[i-1].qryId != matches[i].qryId)    ||
            (matches[i-1].subjId / sNumFrames(TGH::blastProgram)) !=
             (matches[i].subjId / sNumFrames(TGH::blastProgram)))
        {
            if (length(intervals) == 0) // first interval
                intervals.emplace_back(std::make_tuple(matches[i-1].qryId
                                                       / qNumFrames(TGH::blastProgram),
                                                       0,
                                                       i));
            else
                intervals.emplace_back(std::make_tuple(matches[i-1].qryId
                                                       / qNumFrames(TGH::blastProgram),
                                                       std::get<2>(intervals.back()),
                                                       i));
        }
    }

    if (doubleIndexing)
    {
        // sort by trueQryId, then lengths of interval
        std::sort(intervals.begin(), intervals.end(),
                [] (std::tuple<TId, TId, TId> const & i1,
                    std::tuple<TId, TId, TId> const & i2)
        {
            return (std::get<0>(i1) != std::get<0>(i2))
                    ? (std::get<0>(i1) < std::get<0>(i2))
                    : ((std::get<2>(i1) - std::get<1>(i1))
                     > (std::get<2>(i2) - std::get<1>(i2)));
        });
    } else
    {
        // sort by lengths of interval, since trueQryId is the same anyway
        std::sort(intervals.begin(), intervals.end(),
                [] (std::tuple<TId, TId, TId> const & i1,
                    std::tuple<TId, TId, TId> const & i2)
        {
            return (std::get<2>(i1) - std::get<1>(i1))
                >  (std::get<2>(i2) - std::get<1>(i2));
        });
    }

    std::vector<Match<TAlph>> tmpVector;
    tmpVector.resize(matches.size());

    TId newIndex = 0;
    for (auto const & i : intervals)
    {
        TId limit = std::get<2>(i);
        for (TId j = std::get<1>(i); j < limit; ++j)
        {
            tmpVector[newIndex] = matches[j];
            newIndex++;
        }
    }
    std::swap(tmpVector, matches);
}


// ============================================================================
// Alignment-related
// ============================================================================

template <typename T1, typename T2>
inline uint64_t
quickHamming(T1 const & s1, T2 const & s2)
{
    SEQAN_ASSERT_EQ(length(s1), length(s2));

    uint64_t ret = 0;

    for (uint64_t i = 0; i < length(s1); ++i)
        if (s1[i] != s2[i])
            ++ret;

    return ret;
}


template <typename TSource0, typename TGapsSpec0,
          typename TSource1, typename TGapsSpec1,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignContext>
inline TScoreValue
localAlignment2(Gaps<TSource0, TGapsSpec0> & row0,
                Gaps<TSource1, TGapsSpec1> & row1,
                Score<TScoreValue, TScoreSpec> const & scoringScheme,
                int const lowerDiag,
                int const upperDiag,
                TAlignContext & alignContext)
{
    clear(alignContext.traceSegment);

    typedef FreeEndGaps_<True, True, True, True> TFreeEndGaps;
    typedef AlignConfig2<LocalAlignment_<>,
                         DPBand,
                         TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace,
                                                      GapsLeft> > > TAlignConfig;

    TScoreValue score;
    DPScoutState_<Default> scoutState;
    score = _setUpAndRunAlignment(alignContext.dpContext,
                                  alignContext.traceSegment,
                                  scoutState,
                                  row0,
                                  row1,
                                  scoringScheme,
                                  TAlignConfig(lowerDiag, upperDiag));

    _adaptTraceSegmentsTo(row0, row1, alignContext.traceSegment);
    return score;
}


template <typename TLocalHolder>
inline int
_bandSize(uint64_t const seqLength, TLocalHolder & lH)
{
    switch (lH.options.band)
    {
        case -3: case -2:
        {
            int ret = 0;
            auto fit = lH.bandTable.find(seqLength);
            if (fit != lH.bandTable.end())
            {
                ret = fit->second;
            } else
            {
                if (lH.options.band == -3)
                    ret = ceil(std::log2(seqLength));
                else
                    ret = floor(sqrt(seqLength));
            }
            lH.bandTable[seqLength] = ret;
            return ret;
        } break;
        case -1:
            return std::numeric_limits<int>::max();
        default:
            return lH.options.band;
    }
}

// ----------------------------------------------------------------------------
// Function computeEValueThreadSafe
// ----------------------------------------------------------------------------

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeEValueThreadSafe(TBlastMatch & match,
                        uint64_t ql,
                        BlastIOContext<TScore, p, h> & context)
{
#if defined(__FreeBSD__)
    // && version < 11 && defined(STDLIB_LLVM) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=192320
    // || version >= 11 && defined(STDLIB_GNU) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=215709
    static std::vector<std::unordered_map<uint64_t, uint64_t>> _cachedLengthAdjustmentsArray(omp_get_num_threads());
    std::unordered_map<uint64_t, uint64_t> & _cachedLengthAdjustments = _cachedLengthAdjustmentsArray[omp_get_thread_num()];
#else
    static thread_local std::unordered_map<uint64_t, uint64_t> _cachedLengthAdjustments;
#endif

    // convert to 64bit and divide for translated sequences
    ql = ql / (qIsTranslated(context.blastProgram) ? 3 : 1);
    // length adjustment not yet computed
    if (_cachedLengthAdjustments.find(ql) == _cachedLengthAdjustments.end())
        _cachedLengthAdjustments[ql] = _lengthAdjustment(context.dbTotalLength, ql, context.scoringScheme);

    uint64_t adj = _cachedLengthAdjustments[ql];

    match.eValue = _computeEValue(match.alignStats.alignmentScore,
                                  ql - adj,
                                  context.dbTotalLength - adj,
                                  context.scoringScheme);
    return match.eValue;
}

// ----------------------------------------------------------------------------
// compute LCA
// ----------------------------------------------------------------------------

template <typename T, typename T2>
T computeLCA(String<T> const & taxParents, String<T2> const & taxHeights, T n1, T n2)
{
    if (n1 == n2)
        return n1;

    // move up so that nodes are on same height
    for (auto i = taxHeights[n1]; i > taxHeights[n2]; --i)
        n1 = taxParents[n1];

    for (auto i = taxHeights[n2]; i > taxHeights[n1]; --i)
        n2 = taxParents[n2];

    while ((n1 != 0) && ( n2 != 0))
    {
        // common ancestor
        if (n1 == n2)
            return n1;

        // move up in parallel
        n1 = taxParents[n1];
        n2 = taxParents[n2];
    }

    SEQAN_FAIL("One of the paths didn't lead to root.");
    return 0; // avoid warnings on clang
}

#endif // header guard
