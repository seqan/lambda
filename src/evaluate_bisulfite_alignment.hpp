// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2019, Sara Hetzel and MPI für Molekulare Genetik
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
// Scoring schemes for bisulfite converted data on dna5 alphabet
// ==========================================================================

#pragma once

namespace seqan
{

template <typename TSource0, typename TGapsSpec0,
          typename TSource1, typename TGapsSpec1,
          typename TScoreVal>
TScoreVal computeAlignmentStats(AlignmentStats & stats,
                                Gaps<TSource0, TGapsSpec0> const & row0,
                                Gaps<TSource1, TGapsSpec1> const & row1,
                                Score<TScoreVal, ScoreMatrix<Dna5, BisulfiteMatrix>> const & scoringScheme)
{
    clear(stats);

    typedef typename Iterator<Gaps<TSource0, TGapsSpec0> const, Standard>::Type TGapsIter0;
    typedef typename Iterator<Gaps<TSource1, TGapsSpec1> const, Standard>::Type TGapsIter1;
    typedef typename Value<TSource0>::Type TAlphabet;

    // Get iterators.
    TGapsIter0 it0      = begin(row0);
    TGapsIter0 itEnd0   = end(row0);
    TGapsIter1 it1      = begin(row1);
    TGapsIter1 itEnd1   = end(row1);

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;

    for (; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1)
    {
        if (isGap(it0))
        {
            if (!isGapOpen0)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            stats.numDeletions += 1;
            isGapOpen0 = true;
        }
        else
        {
            isGapOpen0 = false;
        }

        if (isGap(it1))
        {
            if (!isGapOpen1)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            stats.numInsertions += 1;
            isGapOpen1 = true;
        }
        else
        {
            isGapOpen1 = false;
        }

        if (!isGap(it0) && !isGap(it1))
        {
            // Compute the alignment score and register in stats.
            TAlphabet c0 = convert<TAlphabet>(*it0);
            TAlphabet c1 = convert<TAlphabet>(*it1);
            TScoreVal scoreVal = score(scoringScheme, c0, c1);
            stats.alignmentScore += scoreVal;
            // Register other statistics.
            bool isMatch = score(scoringScheme, c0, c1) == score(scoringScheme, c0, c0);
            bool isPositive = (scoreVal > 0);
            stats.numMatches += isMatch;
            stats.numMismatches += !isMatch;
            stats.numPositiveScores += isPositive;
            stats.numNegativeScores += !isPositive;
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);

    stats.numGaps = stats.numGapOpens + stats.numGapExtensions;

    // Finally, compute the alignment similarity from the various counts
    stats.alignmentLength = length(row0);
    stats.alignmentSimilarity = 100.0 * static_cast<float>(stats.numPositiveScores)
                                / static_cast<float>(stats.alignmentLength);
    stats.alignmentIdentity = 100.0 * static_cast<float>(stats.numMatches)
                              / static_cast<float>(stats.alignmentLength);

    return stats.alignmentScore;
}

} // namespace seqan
