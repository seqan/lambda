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
// view_dna_n_to_random.hpp: View that converts N into random selection of {A,C,G,T}
// ==========================================================================

#pragma once

#include <random>

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/std/ranges>

struct dna_n_to_random_fn
{
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to dna_n_to_random cannot be a temporary of a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
            "The range parameter to dna_n_to_random must model std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
            "The range parameter to dna_n_to_random must model std::ranges::random_access_range.");
        static_assert(seqan3::nucleotide_alphabet<seqan3::innermost_value_type_t<urng_t>>,
            "The range parameter to dna_n_to_random must be over elements of seqan3::nucleotide_alphabet.");

        std::mt19937 rng(0xDEADBEEF);

        return std::forward<urng_t>(urange) | std::views::transform([rng] (auto const in) mutable
        {
            if (seqan3::to_char(in) == 'N')
            {
                return seqan3::assign_rank_to(rng() % 4, seqan3::dna4{});
            }
            else
            {
                return static_cast<seqan3::dna4>(in);
            }
        });
    }

    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, dna_n_to_random_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

namespace seqan3::views
{

inline constexpr auto dna_n_to_random = seqan3::views::deep{dna_n_to_random_fn{}};

} // namespace seqan3::view
