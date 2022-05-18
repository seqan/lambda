// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2019, Hannes Hauswedell <h2 @ fsfe.org>
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
// view_pos_transform.hpp
// ==========================================================================

#pragma once

#include <ranges>

#include <seqan3/core/range/detail/random_access_iterator.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/utility/range/concept.hpp>

// The return type of views::pos_transform
template <std::ranges::view urng_t, typename pos_transform_t, typename size_transform_t>
class view_pos_transform : public std::ranges::view_base
{

private:
    // The data members of view_pos_transform.
    urng_t                                        urange;
    ::ranges::semiregular_box_t<pos_transform_t>  pos_transform;
    ::ranges::semiregular_box_t<size_transform_t> size_transform;

    // Associated types iterator.
    using reference         = decltype(pos_transform(urange, 0));
    using const_reference   = decltype(pos_transform(std::as_const(urange), 0));
    using value_type        = std::remove_cvref_t<reference>;
    using size_type         = std::ranges::range_size_t<urng_t>;
    using difference_type   = std::ranges::range_difference_t<urng_t>;
    using iterator          = seqan3::detail::random_access_iterator<view_pos_transform>;
    using const_iterator    = seqan3::detail::random_access_iterator<view_pos_transform const>;

    // Befriend the following class s.t. iterator and const_iterator can be defined for this type.
    template <typename, typename>
    friend class seqan3::detail::random_access_iterator_base;

public:

    static_assert(std::ranges::sized_range<urng_t>,
                  "The range parameter to views::pos_transform must model std::ranges::sized_range.");
    static_assert(std::ranges::random_access_range<urng_t>,
                  "The range parameter to views::pos_transform must model std::ranges::random_access_range.");
    static_assert(std::ranges::sized_range<urng_t>,
                  "The range parameter to views::pos_transform must model std::ranges::sized_range.");

    // Constructors, destructor and assignment
    view_pos_transform()                                                     noexcept = default; //!< Defaulted.
    constexpr view_pos_transform(view_pos_transform const & rhs)             noexcept = default; //!< Defaulted.
    constexpr view_pos_transform(view_pos_transform && rhs)                  noexcept = default; //!< Defaulted.
    constexpr view_pos_transform & operator=(view_pos_transform const & rhs) noexcept = default; //!< Defaulted.
    constexpr view_pos_transform & operator=(view_pos_transform && rhs)      noexcept = default; //!< Defaulted.
    ~view_pos_transform()                                                    noexcept = default; //!< Defaulted.

    // Construct from another view
    view_pos_transform(urng_t _urange, pos_transform_t _pos_transform, size_transform_t _size_transform)
        : urange{std::move(_urange)},
          pos_transform{std::move(_pos_transform)},
          size_transform{std::move(_size_transform)}
    {}

    // Construct from another range
    template <typename rng_t>
        requires (!std::same_as<std::remove_cvref_t<rng_t>, view_pos_transform> &&
                  (std::ranges::viewable_range<rng_t> &&
                  std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<rng_t>>>))
    view_pos_transform(rng_t && _urange, pos_transform_t _pos_transform, size_transform_t _size_transform)
     : view_pos_transform{std::views::all(std::forward<rng_t>(_urange)),
                          std::move(_pos_transform),
                          std::move(_size_transform)}
    {}

    // Iterators
    iterator begin() noexcept
    {
        return {*this, 0};
    }

    const_iterator begin() const noexcept
        requires seqan3::const_iterable_range<urng_t>
    {
        return {*this, 0};
    }

    iterator end() noexcept
    {
        return {*this, size()};
    }

    const_iterator end() const noexcept
        requires seqan3::const_iterable_range<urng_t>
    {
        return {*this, size()};
    }

    // Returns the number of elements in the view.
    size_type size() noexcept
    {
        return size_transform(urange);
    }

    size_type size() const noexcept
        requires seqan3::const_iterable_range<urng_t>
    {
        return size_transform(urange);
    }

    // Element access
    reference operator[](size_type const n)
    {
        assert(n < size());
        return pos_transform(urange, n);
    }

    const_reference operator[](size_type const n) const
        requires seqan3::const_iterable_range<urng_t>
    {
        assert(n < size());
        return pos_transform(urange, n);
    }
};

// Template argument deduction for view_pos_transform.
template <typename urng_t, typename pos_transform_t, typename size_transform_t>
view_pos_transform(urng_t &&, pos_transform_t, size_transform_t) -> view_pos_transform<std::views::all_t<urng_t>, pos_transform_t, size_transform_t>;

// Definition of the range adaptor object type for views::pos_transform.
struct pos_transform_fn
{
private:
    struct default_size_fn
    {
        template <typename urng_t>
        auto operator()(urng_t && urange) const
        {
            return std::ranges::size(urange);
        }
    };

public:
    template <typename pos_transform_t, typename size_transform_t = default_size_fn>
    constexpr auto operator()(pos_transform_t pos_transform, size_transform_t size_transform = default_size_fn{}) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(pos_transform), std::move(size_transform)};
    }

    template <std::ranges::range urng_t,
              typename pos_transform_t,
              typename size_transform_t = default_size_fn>
    constexpr auto operator()(urng_t && urange,
                              pos_transform_t pos_transform,
                              size_transform_t size_transform = default_size_fn{}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::pos_transform cannot be a temporary of a non-view range.");
        static_assert(std::ranges::sized_range<urng_t>,
                      "The range parameter to views::pos_transform must model std::ranges::sized_range.");
        static_assert(std::ranges::random_access_range<urng_t>,
                      "The range parameter to views::pos_transform must model std::ranges::random_access_range.");

        return view_pos_transform{std::forward<urng_t>(urange), std::move(pos_transform), std::move(size_transform)};
    }
};

// A view that applies a transformation function to an element of a range using its positional index.
namespace views
{

inline constexpr auto pos_transform = pos_transform_fn{};

} // namespace views
