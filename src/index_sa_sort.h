// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
// ==========================================================================

#ifndef LAMBDA_INDEX_SA_SORT_H
#define LAMBDA_INDEX_SA_SORT_H

#include <atomic>
#ifdef _OPENMP
#include <parallel/algorithm>
#endif


namespace SEQAN_NAMESPACE_MAIN
{

template <typename TAlgoSpec>
struct SaAdvancedSort {};

struct QuickSortTag {};

struct MergeSortTag {};

template <typename TIndex>
struct SaAdvancedSortAlgoTag
{
#ifdef _OPENMP
    typedef __gnu_parallel::default_parallel_tag Type;
#else
    typedef void Type;
#endif
};

template <>
struct SaAdvancedSortAlgoTag<MergeSortTag>
{
#ifdef _OPENMP
    typedef __gnu_parallel::multiway_mergesort_tag Type;
#else
    typedef void Type;
#endif
};

template <>
struct SaAdvancedSortAlgoTag<QuickSortTag>
{
#ifdef _OPENMP
    typedef __gnu_parallel::quicksort_tag Type;
#else
    typedef void Type;
#endif
};


// compare two suffices of a given text
template < typename TSAValue, typename TText >
struct AdvancedSuffixLess_ :
    public ::std::binary_function <TSAValue, TSAValue, bool>
{
    typedef typename Iterator<TText const, Standard>::Type TIter;
    typedef typename Size<TText>::Type TOffset;
    TIter _begin, _end;
    std::function<void (void)> const & _lambda;

    AdvancedSuffixLess_(TText & text,
                        TOffset offset = 0,
                        std::function<void (void)> lambda = [] () {} ) :
        _begin(begin(text, Standard()) + offset),
        _end(end(text, Standard())),
        _lambda(lambda)
    {}

    inline bool operator() (TSAValue const a, TSAValue const b)
    {
        _lambda();

        if (a == b) return false;
        TIter itA = _begin + a;
        TIter itB = _begin + b;
        if (a <= b)
        {
            for (; itB != _end; ++itB, ++itA)
            {
                if (ordLess(*itA, *itB))
                    return true;
                if (ordLess(*itB, *itA))
                    return false;
            }
            return false;
        } else
        {
            for (; itA != _end; ++itA, ++itB)
            {
                if (ordLess(*itA, *itB))
                    return true;
                if (ordLess(*itB, *itA))
                    return false;
            }
            return true;
        }
    }
};

// compare two suffices of a given text
template < typename TSAValue, typename TString, typename TSetSpec>
struct AdvancedSuffixLess_<TSAValue, StringSet<TString, TSetSpec> const> :
    public ::std::binary_function <TSAValue, TSAValue, bool>
{
    typedef StringSet<TString, TSetSpec> const TText;
    typedef typename Size<TString>::Type TOffset;

    TText & _text;
    TOffset _offset;
    std::function<void (void)> const & _lambda;

    AdvancedSuffixLess_(TText & text,
                        TOffset offset = 0,
                        std::function<void (void)> lambda = [] () {} ) :
        _text(text), _offset(offset), _lambda(lambda)
    {}

    inline bool operator() (TSAValue const a, TSAValue const b)
    {
        _lambda();

        typedef typename Iterator<TString const, Standard>::Type TIter;
        if (a == b) return false;
        TIter itA = begin(getValue(_text, getSeqNo(a)), Standard()) +
                    getSeqOffset(a) + _offset;
        TIter itB = begin(getValue(_text, getSeqNo(b)), Standard()) +
                    getSeqOffset(b) + _offset;
        TIter itAEnd = end(getValue(_text, getSeqNo(a)), Standard());
        TIter itBEnd = end(getValue(_text, getSeqNo(b)), Standard());
        if (itAEnd - itA < itBEnd - itB)
        {
            for(; itA != itAEnd; ++itA, ++itB)
            {
                if (ordLess(*itA, *itB))
                    return true;
                if (ordLess(*itB, *itA))
                    return false;
            }
            return true;
        } else
        {
            for (; itB != itBEnd; ++itB, ++itA)
            {
                if (ordLess(*itA, *itB))
                    return true;
                if (ordLess(*itB, *itA))
                    return false;
            }
            if (itA != itAEnd)
                return false;
            return getSeqNo(a) > getSeqNo(b);
        }
    }
};

// template < typename TSAValue, typename TString, typename TSetSpec >
// std::atomic<uint64_t> AdvancedSuffixLess_<TSAValue, StringSet<TString, TSetSpec> const >::_comparisons;

template < typename TSA, 
        typename TText,
        typename TSize>
void _sortBucketAdvancedSort(
    TSA &sa,
    TText &text,
    TSize lcp)
{
    // sort bucket with quicksort
    ::std::sort(
        begin(sa, Standard()),
        end(sa, Standard()),
        AdvancedSuffixLess_<typename Value<TSA>::Type, TText>(text, lcp));
}

template <typename TSA, typename TText, typename TAlgo>
inline void createSuffixArray(
    TSA &SA,
    TText const &s,
    SaAdvancedSort<TAlgo> const &)
{
    typedef typename Size<TSA>::Type TSize;
    typedef typename Iterator<TSA, Standard>::Type TIter;

    // 1. Fill suffix array with a permutation (the identity)
    TIter it = begin(SA, Standard());
    TIter itEnd = end(SA, Standard());
    for(TSize i = 0; it != itEnd; ++it, ++i)
        *it = i;

    // 2. Sort suffix array with quicksort
    ::std::sort(
        begin(SA, Standard()),
        end(SA, Standard()),
        SuffixLess_<typename Value<TSA>::Type, TText const>(s));
}

template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgoSpec>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  SaAdvancedSort<TAlgoSpec> const &,
                  std::function<void(void)> progressCallback = [] () {})
{
    typedef StringSet< TString, TSSetSpec > TText;
    typedef typename Size<TSA>::Type TSize;
    typedef typename Iterator<TSA, Standard>::Type TIter;
    typedef typename SaAdvancedSortAlgoTag<TAlgoSpec>::Type TAlgo;

    // 1. Fill suffix array with a permutation (the identity)
    TIter it = begin(SA, Standard());
    for(unsigned j = 0; j < length(s); ++j)
    {
        TSize len = length(s[j]);
        for(TSize i = 0; i < len; ++i, ++it)
            *it = Pair<unsigned, TSize>(j, i);
    }

    // 2. Sort suffix array with quicksort
    if (std::is_same<TAlgoSpec, std::true_type>::value)
    {
        std::cout << "calling std::sort\n";
        std::sort(
        begin(SA, Standard()),
        end(SA, Standard()),
        AdvancedSuffixLess_<typename Value<TSA>::Type, TText const>(s, 0, progressCallback));
    } else if (std::is_same<TAlgoSpec, std::false_type>::value)
    {
        std::sort(
        begin(SA, Standard()),
        end(SA, Standard()),
        SuffixLess_<typename Value<TSA>::Type, TText const>(s));
    } else
#ifdef _OPENMP
    __gnu_parallel::sort(
        begin(SA, Standard()),
        end(SA, Standard()),
        AdvancedSuffixLess_<typename Value<TSA>::Type, TText const>(s, 0, progressCallback),
        TAlgo());
#else
    std::sort(
        begin(SA, Standard()),
        end(SA, Standard()),
        AdvancedSuffixLess_<typename Value<TSA>::Type, TText const>(s, 0, progressCallback));
#endif

}

template <typename TInput, typename TAlgoSpec>
struct Pipe< TInput, SaAdvancedSort<TAlgoSpec> >
{
    typedef typename Value<TInput>::Type    TValue;
    typedef typename SAValue<TInput>::Type    TSAValue;

    typedef String<TValue, Alloc<> >        TText;
    typedef String<TSAValue, Alloc<> >        TSA;
    typedef Pipe<TSA, Source<> >            TSource;

    TSA        sa;
    TSource    in;

    Pipe(TInput &_textIn):
        in(sa)
    {
        TText text;
        text << _textIn;

        resize(sa, length(_textIn), Exact());
        createSuffixArray(sa, text, SaAdvancedSort<TAlgoSpec>());
    }

    inline typename Value<TSource>::Type const & operator*()
    {
        return *in;
    }

    inline Pipe& operator++()
    {
        ++in;
        return *this;
    }
};

}

#endif
