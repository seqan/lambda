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
// view_reduce_to_bisulfite.hpp: View that converts dna sequences into bisulfite
// forward and reverse sequences
// ==========================================================================

#pragma once

#include <seqan3/alphabet/composite/semialphabet_any.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/std/ranges>
#include "bisulfite_scoring.hpp"
#include "view_pos_transform.hpp"

// Definition of the range adaptor object type for views::to_bisulfite_semialphabet.
struct to_bisulfite_semialphabet_fn
{
    constexpr auto operator()(bsDirection const direction) const
    {
        return seqan3::detail::adaptor_from_functor{*this, direction};
    }

    // Rank
    // 0 = A forward
    // 1 = C/T forward
    // 2 = G forward
    // 3 = A/G reverse
    // 4 = C reverse
    // 5 = T reverse
private:
    static constexpr std::array<uint8_t, 4> dna4_to_rank_bs_fwd = {0, 1, 2, 1};
    static constexpr std::array<uint8_t, 4> dna4_to_rank_bs_rev = {3, 4, 3, 5};

    static auto func_fwd(auto && urange, size_t pos) { return seqan3::assign_rank_to(dna4_to_rank_bs_fwd[seqan3::to_rank(urange[pos])], seqan3::semialphabet_any<6>{}); }
    static auto func_rev(auto && urange, size_t pos) { return seqan3::assign_rank_to(dna4_to_rank_bs_rev[seqan3::to_rank(urange[pos])], seqan3::semialphabet_any<6>{}); }

public:
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, bsDirection const direction) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::to_bisulfite_semialphabet cannot be a temporary of a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
            "The range parameter to views::to_bisulfite_semialphabet must model std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
            "The range parameter to views::to_bisulfite_semialphabet must model std::ranges::random_access_range.");
        static_assert(std::is_same_v<std::remove_const_t<seqan3::range_innermost_value_t<urng_t>>, seqan3::dna4>,
            "The range parameter to views::to_bisulfite_semialphabet must be over elements of seqan3::dna4.");

        auto l = &func_fwd<std::views::all_t<urng_t> const &>;
        if (direction == bsDirection::rev)
            l = &func_rev<std::views::all_t<urng_t> const &>;
        return std::forward<urng_t>(urange) | views::pos_transform(l);
    }
};

// A view that converts elements of the seqan3::dna4 alphabet to a semialphabet that models reduced forward and reverse
// alphabets for the bisulfite mode. Depending on the `direction` parameter, letters get converted to the rank of the bisulfite
// reduced alphabet for forward sequences or reverse sequences which both are included in the same semialphabet of size 6.
namespace views
{

inline constexpr auto to_bisulfite_semialphabet = seqan3::views::deep{to_bisulfite_semialphabet_fn{}};

} // namespace views

// Definition of the range adaptor object type for views::reduce_to_bisulfite.
struct reduce_to_bisulfite_fn
{
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange) const
    {
        static_assert(seqan3::range_dimension_v<urng_t> == 2,
                      "This adaptor only handles range-of-range (two dimensions) as input.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::reduce_to_bisulfite cannot be a temporary of a non-view range.");
        static_assert(std::ranges::viewable_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::reduce_to_bisulfite cannot be a "
                      "temporary of a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
                      "The range parameter to views::reduce_to_bisulfite must model std::ranges::sized_range.");
        static_assert(std::ranges::sized_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::reduce_to_bisulfite must model "
                      "std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
                      "The range parameter to views::reduce_to_bisulfite must model std::ranges::random_access_range.");
        static_assert(std::ranges::random_access_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::reduce_to_bisulfite must model "
                      "std::ranges::random_access_range.");
        static_assert(std::is_same_v<std::remove_const_t<seqan3::range_innermost_value_t<urng_t>>, seqan3::dna4>,
                      "The range parameter to views::reduce_to_bisulfite must be over a range over elements of "
                      "seqan3::dna4.");

        return std::forward<urng_t>(urange) | views::pos_transform([] (auto && urange, size_t pos)
        {
            return urange[pos] | views::to_bisulfite_semialphabet((bsDirection) (pos % 2));
        });
    }

    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, reduce_to_bisulfite_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

// A view that operates on a range-of-ranges where the inner range needs to be of type seqan3::dna4.
// Inner ranges with an even index get reduced to forward bisulfite sequences while inner ranges with an odd index get
// converted to reverse bisulfite sequences. The output range is a range-of-range where the type of the inner ranges is
// a semialphabet of size 6.
namespace views
{

inline constexpr auto reduce_to_bisulfite = reduce_to_bisulfite_fn{};

} // namespace views
