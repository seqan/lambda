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

#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <algorithm>

enum class bsDirection
{
    fwd,
    rev
};

template <seqan3::arithmetic score_type = int8_t>
class bisulfite_scoring_scheme : public seqan3::scoring_scheme_base<bisulfite_scoring_scheme<score_type>, seqan3::dna5, score_type>
{
private:
    using base_t = seqan3::scoring_scheme_base<bisulfite_scoring_scheme<score_type>, seqan3::dna5, score_type>;
    using base_t::matrix;
public:
    using base_t::base_t;
    using typename base_t::matrix_type;
    using matrix_size_type = std::remove_const_t<decltype(seqan3::alphabet_size<seqan3::dna5>)>;
    static constexpr matrix_size_type matrix_size = seqan3::alphabet_size<seqan3::dna5>;

    constexpr bisulfite_scoring_scheme() noexcept {}

    template <seqan3::arithmetic score_arg_t>
    constexpr bisulfite_scoring_scheme(seqan3::match_score<score_arg_t> const ms,
        seqan3::mismatch_score<score_arg_t> const mms, bsDirection const dir = bsDirection::fwd)
    {
        set_bisulfite_scheme(ms, mms, dir);
    }

    constexpr bisulfite_scoring_scheme(matrix_type const & _matrix) noexcept {}

    template <seqan3::arithmetic score_arg_t>
    constexpr void set_bisulfite_scheme(seqan3::match_score<score_arg_t> const ms,
        seqan3::mismatch_score<score_arg_t> const mms, bsDirection const dir = bsDirection::fwd)
    {
        std::conditional_t<std::integral<score_type>, int64_t, double> i_ms = static_cast<score_arg_t>(ms);
        std::conditional_t<std::integral<score_type>, int64_t, double> i_mms = static_cast<score_arg_t>(mms);
        if ((i_ms  < std::numeric_limits<score_type>::lowest() || i_ms  > std::numeric_limits<score_type>::max()) ||
            (i_mms < std::numeric_limits<score_type>::lowest() || i_mms > std::numeric_limits<score_type>::max()))
        {
            throw std::invalid_argument{"You passed a score value to set_bisulfite_scheme that is out of range of the "
                                        "scoring scheme's underlying type. Define your scoring scheme with a larger "
                                        "template parameter or down-cast you score value beforehand to prevent "
                                        "this exception."};
        }

        // Assume query is horizontal sequence
        for (matrix_size_type i = 0; i < matrix_size; ++i)
            for (matrix_size_type j = 0; j < matrix_size; ++j)
                matrix[i][j] = (i == j) ? static_cast<score_type>(i_ms) : static_cast<score_type>(i_mms);

        // Explicitly override for C>T or G>A
        if (dir == bsDirection::fwd)
            matrix[4][1] = static_cast<score_type>(i_ms);
        else
            matrix[0][2] = static_cast<score_type>(i_ms);
    }
};

bisulfite_scoring_scheme() -> bisulfite_scoring_scheme<int8_t>;

template <seqan3::arithmetic score_arg_type>
bisulfite_scoring_scheme(seqan3::match_score<score_arg_type>,
                         seqan3::mismatch_score<score_arg_type>) -> bisulfite_scoring_scheme<int8_t>;

template <seqan3::arithmetic score_arg_type>
bisulfite_scoring_scheme(std::array<std::array<score_arg_type, 5>, 5>) -> bisulfite_scoring_scheme<score_arg_type>;

namespace seqan
{

struct BisulfiteMatrix{};

template <>
struct ScoringMatrixData_<int, Dna5, BisulfiteMatrix>
{
    enum
    {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
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
inline void
setScoreBisulfiteMatrix(Score<int, ScoreMatrix<Dna5, BisulfiteMatrix> > & sc, T matchScore, T mismatchScore, bsDirection const dir = bsDirection::fwd)
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
