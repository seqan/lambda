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
// view_duplicate.hpp: View that duplicates a range in a range-of-ranges
// ==========================================================================

#pragma once

#include <ranges>

#include <bio/ranges/views/transform_by_pos.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/range/concept.hpp>

// Definition of the range adaptor object type for views::duplicate.
struct duplicate_fn
{
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange) const
    {
        static_assert(bio::ranges::range_dimension_v<urng_t> == 2,
                      "This adaptor only handles range-of-range (two dimensions) as input.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::duplicate cannot be a temporary of a non-view range.");
        static_assert(std::ranges::viewable_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::duplicate cannot be a temporary of "
                      "a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
                      "The range parameter to views::duplicate must model std::ranges::sized_range.");
        static_assert(std::ranges::sized_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::duplicate must model "
                      "std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
                      "The range parameter to views::duplicate must model std::ranges::random_access_range.");
        static_assert(std::ranges::random_access_range<std::ranges::range_reference_t<urng_t>>,
                      "The inner range of the range parameter to views::duplicate must model "
                      "std::ranges::random_access_range.");

        return std::forward<urng_t>(urange) |
               bio::views::transform_by_pos([](auto && urange, size_t pos) -> decltype(auto)
                                            { return urange[(pos - (pos % 2)) / 2]; },
                                            [](auto && urange) { return std::ranges::size(urange) * 2; });
    }

    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, duplicate_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

// A view that duplicates the inner ranges in a range-of-ranges. The output range has twice the size of the input range
// where duplicates appear directly after each other.
namespace views
{

inline constexpr auto duplicate = duplicate_fn{};

} // namespace views
