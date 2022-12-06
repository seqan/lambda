// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Tristate Holder Implementation.
// ==========================================================================

#pragma once

#include <seqan/basic.h>
#include <bio/ranges/views/slice.hpp>
#include <bio/ranges/views/translate_single.hpp>

using overload_t = decltype(std::declval<std::vector<bio::alphabet::dna5> &>() | bio::views::translate_single | bio::views::slice(0,1));

template <typename t>
concept overload_c = std::same_as<std::remove_cvref_t<t>, overload_t>;

namespace seqan
{

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
template <overload_c TValue>
struct Holder<TValue, Tristate>
{
    enum EHolderState
    {
        EMPTY = 0,
        OWNER = 1,
        DEPENDENT = 2
    };

    typedef std::optional<TValue> THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    THostValue data_value;
    EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

    Holder() : data_value{}, data_state(EMPTY)
    {
    }

    Holder(Holder const & source_) : data_value{}, data_state(EMPTY)
    {
        data_value = source_;
    }

    explicit
    Holder(THostValue & value_) : data_value{}, data_state(EMPTY)
    {
        data_value = value_;
    }

    explicit
    Holder(THostValue const & value_) : data_value{}, data_state(EMPTY)
    {
        data_value = value_;
    }

    ~Holder()
    {
        clear(*this);
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

    Holder & operator=(Holder const &) = default;

    Holder & operator=(THostValue const & value_)
    {
        data_value = value_;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

    inline operator THostValue()
    {
        return data_value;
    }
};

// ============================================================================
// Functions
// ============================================================================

/// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <overload_c TValue>
inline void
clear(Holder<TValue, Tristate> & me)
{
    me.data_value.reset();
    me.data_state = Holder<TValue, Tristate>::EMPTY;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <overload_c TValue>
inline void
setValue(Holder<TValue, Tristate> & me,
         TValue const & value_)
{
    me.data_value = value_;
    me.data_state = Holder<TValue, Tristate>::OWNER;
}

template <overload_c TValue>
inline void
setValue(Holder<TValue const, Tristate> & me,
         TValue const & value_)
{
    me.data_value = value_;
    me.data_state = Holder<TValue const, Tristate>::OWNER;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <overload_c TValue>
auto & value(Holder<TValue, Tristate> & me)
{
    SEQAN_ASSERT_NOT(empty(me)); // HERE BE DRAGONS!!!
    return * me.data_value;
}

template <overload_c TValue>
auto & value(Holder<TValue, Tristate> const & me)
{
    SEQAN_ASSERT_NOT(empty(me));

    return * me.data_value;
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <overload_c TValue, typename TSource>
inline void
assignValue(Holder<TValue, Tristate> & me,
            TSource const & value_)
{
    setValue(me, value_);
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <overload_c TValue, typename TSource>
inline void
moveValue(Holder<TValue, Tristate> & me,
          TSource const & value_)
{
    setValue(me, value_);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <overload_c TValue>
inline void
assign(Holder<TValue, Tristate> & target_,
       Holder<TValue, Tristate> const & source_)
{
    target_ = source_;
}

template <overload_c TValue>
inline void
assign(Holder<TValue const, Tristate> & target_,
       Holder<TValue const, Tristate> const & source_)
{
    target_ = source_;
}

}  // namespace seqan
