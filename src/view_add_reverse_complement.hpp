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
// view_add_reverse_complement.hpp: View that returns the reverse complement
// of every second range in a range-of-ranges
// ==========================================================================

#pragma once

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/std/ranges>
#include "view_pos_transform.hpp"

// Definition of the range adaptor object type for views::reverse_complement_or_not.
struct reverse_complement_or_not_fn
{
private:
    static auto func_nop(auto && urange, size_t pos) { return urange[pos]; }
    static auto func_revcomp(auto && urange, size_t pos) { return seqan3::complement(urange[std::ranges::size(urange) - pos - 1]); }

public:
    constexpr auto operator()(bool const transform) const
    {
        return seqan3::detail::adaptor_from_functor{*this, transform};
    }

    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, bool const transform) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::reverse_complement_or_not cannot be a temporary of a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
            "The range parameter to views::reverse_complement_or_not must model std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
            "The range parameter to views::reverse_complement_or_not must model std::ranges::random_access_range.");
        static_assert(seqan3::nucleotide_alphabet<seqan3::range_innermost_value_t<urng_t>>,
            "The range parameter to views::reverse_complement_or_not must be over elements of seqan3::nucleotide_alphabet.");

        auto l = &func_nop<std::views::all_t<urng_t> const &>;
        if (transform)
            l = &func_revcomp<std::views::all_t<urng_t> const &>;
        return std::forward<urng_t>(urange) | views::pos_transform(l);
    }
};

// A view that reverse complements elements of a nucleotide alphabet or returns the elements itself depending on the
// `transform` parameter.
namespace views
{

inline constexpr auto reverse_complement_or_not = seqan3::views::deep{reverse_complement_or_not_fn{}};

} // namespace views

// Definition of the range adaptor object type for views::add_reverse_complement.
struct add_reverse_complement_fn
{
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange) const
    {
        static_assert(seqan3::range_dimension_v<urng_t> == 2,
                      "This adaptor only handles range-of-range (two dimensions) as input.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::add_reverse_complement cannot be a temporary of a non-view range.");
        static_assert(std::ranges::viewable_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::add_reverse_complement cannot be a temporary of "
                      "a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
                      "The range parameter to views::add_reverse_complement must model std::ranges::sized_range.");
        static_assert(std::ranges::sized_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::add_reverse_complement must model "
                      "std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
                      "The range parameter to views::add_reverse_complement must model std::ranges::random_access_range.");
        static_assert(std::ranges::random_access_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::add_reverse_complement must model "
                      "std::ranges::random_access_range.");
        static_assert(seqan3::nucleotide_alphabet<seqan3::range_innermost_value_t<urng_t>>,
                      "The range parameter to views::add_reverse_complement must be over elements of seqan3::nucleotide_alphabet.");

        return std::forward<urng_t>(urange) | views::pos_transform([] (auto && urange, size_t pos)
        {
            return urange[pos / 2] | views::reverse_complement_or_not(pos % 2 == 1);
        }, [] (auto && urange) { return std::ranges::size(urange) * 2; });
    }

    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, add_reverse_complement_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

// A view that reverse complements every second inner range in a range-of-ranges.
namespace views
{

inline constexpr auto add_reverse_complement = add_reverse_complement_fn{};

} // namespace views
