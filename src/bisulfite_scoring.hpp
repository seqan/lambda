// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2019, Sara Hetzel <hetzel @ molgen.mpg.de>
// Copyright (c) 2016-2019, Knut Reinert and Freie Universität Berlin
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

#include <seqan/score.h>
#include <seqan/sequence.h>

#include <algorithm>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

enum class bsDirection
{
    fwd,
    rev
};

namespace seqan
{

struct BisulfiteMatrix
{};

template <>
struct ScoringMatrixData_<int, Dna5, BisulfiteMatrix>
{
    enum
    {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE   = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData()
    {
        // clang-format off
        static int const _data[TAB_SIZE] =
        {
             0, -1, -1, -1, -1,
            -1,  0, -1, -1, -1,
            -1, -1,  0, -1, -1,
            -1,  0, -1,  0, -1,
            -1, -1, -1, -1, -1
        };
        // clang-format on
        return _data;
    }
};

template <typename T>
inline void setScoreBisulfiteMatrix(Score<int, ScoreMatrix<Dna5, BisulfiteMatrix>> & sc,
                                    T                                                matchScore,
                                    T                                                mismatchScore,
                                    bsDirection const                                dir = bsDirection::fwd)
{
    for (size_t i = 0; i < ValueSize<Dna5>::VALUE; ++i)
    {
        for (size_t j = 0; j < ValueSize<Dna5>::VALUE; ++j)
        {
            if (dir == bsDirection::fwd)
            {
                if (((i == j) || (i == 3 && j == 1)) && i != 4)
                    setScore(sc, i, j, matchScore);
                else
                    setScore(sc, i, j, mismatchScore);
            }
            else
            {
                if (((i == j) || (i == 0 && j == 2)) && i != 4)
                    setScore(sc, i, j, matchScore);
                else
                    setScore(sc, i, j, mismatchScore);
            }
        }
    }
}

} // namespace seqan
