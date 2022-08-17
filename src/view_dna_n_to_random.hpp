// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2019, Sara Hetzel and MPI für Molekulare Genetik
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
// view_dna_n_to_random.hpp: View that converts N into random selection of
// {A,C,G,T} if target alphabet is dna4
// ==========================================================================

#pragma once

#include <random>
#include <ranges>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/views/deep.hpp>
#include "view_pos_transform.hpp"

// Definition of the range adaptor object type for views::dna_n_to_random.
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
        static_assert(std::same_as<seqan3::dna5, seqan3::range_innermost_value_t<urng_t>>,
                      "The range parameter to dna_n_to_random must be over elements of seqan3::dna5.");

        std::shared_ptr<std::mt19937> rng{new std::mt19937{0xDEADBEEF}};

        return std::forward<urng_t>(urange) | std::views::transform(
                                                [rng](seqan3::dna5 const c)
                                                {
                                                    return (seqan3::to_char(c) == 'N')
                                                             ? seqan3::assign_rank_to((*rng)() % 4, seqan3::dna4{})
                                                             : static_cast<seqan3::dna4>(c);
                                                });
    }

    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, dna_n_to_random_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

// A view that converts every ambiguous nucleotide 'N' randomly into a letter of a new alphabet (e.g. seqan3::dna4).
namespace views
{

inline constexpr auto dna_n_to_random = seqan3::views::deep{dna_n_to_random_fn{}};

} // namespace views
