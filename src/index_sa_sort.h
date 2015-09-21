// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
#if defined(_OPENMP) && defined(__GNUC__) && !defined(__clang__)
#define GNUOMP 1
#include <parallel/algorithm>
#endif

namespace seqan
{

// ==========================================================================
// Tags
// ==========================================================================

template <typename TAlgoSpec>
struct SaAdvancedSort {};

struct QuickSortTag {};

struct MergeSortTag {};

struct QuickSortBucketTag {};

template <typename TIndex>
struct SaAdvancedSortAlgoTag
{
#ifdef GNUOMP
    typedef __gnu_parallel::default_parallel_tag Type;
#else
    typedef void Type;
#endif
};

template <>
struct SaAdvancedSortAlgoTag<MergeSortTag>
{
#ifdef GNUOMP
    typedef __gnu_parallel::multiway_mergesort_tag Type;
#else
    typedef void Type;
#endif
};

template <>
struct SaAdvancedSortAlgoTag<QuickSortTag>
{
#ifdef GNUOMP
    typedef __gnu_parallel::quicksort_tag Type;
#else
    typedef void Type;
#endif
};

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TSA, typename TText, typename TAlgo>
inline void
createSuffixArray(TSA &SA, TText const &s, SaAdvancedSort<TAlgo> const &);

} // namespace seqan

// ==========================================================================
// Metafunctions
// ==========================================================================

namespace SEQAN_NAMESPACE_MAIN
{
template <typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSa<SaAdvancedSort<TSpec> > >, FibreTempSA>
{
    typedef Index<TText, IndexSa<SaAdvancedSort<TSpec> > >          TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template <typename TText, typename TSpec>
struct DefaultIndexCreator<Index<TText, IndexSa<SaAdvancedSort<TSpec>> >, FibreSA>
{
    typedef SaAdvancedSort<TSpec> Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >, FibreTempSA>
{
    typedef Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >  TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template < typename TText, typename TSpec, typename TConfig>
struct DefaultIndexCreator<Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >, FibreSA>
{
    typedef SaAdvancedSort<TSpec> Type;
};

} // SEQAN_NAMESPACE_MAIN

// ==========================================================================
// Classes
// ==========================================================================

namespace seqan
{

// ----------------------------------------------------------------------------
// Class AdvancedSuffixLess_
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Class Pipe
// ----------------------------------------------------------------------------

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

// template < typename TSAValue, typename TString, typename TSetSpec >
// std::atomic<uint64_t> AdvancedSuffixLess_<TSAValue, StringSet<TString, TSetSpec> const >::_comparisons;

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Helper Functions
// ----------------------------------------------------------------------------

template <typename TSA,
         typename TText,
         typename TSize>
inline void
_sortBucketAdvancedSort(TSA &sa, TText &text, TSize lcp)
{
    // sort bucket with quicksort
    ::std::sort(
        begin(sa, Standard()),
        end(sa, Standard()),
        AdvancedSuffixLess_<typename Value<TSA>::Type, TText>(text, lcp));
}

template <typename T>
inline void
progressCallbackWrapper(std::function<void(void)>, T const &)
{
}

inline void
progressCallbackWrapper(std::function<void(uint64_t)> const & lambda, uint64_t v)
{
    lambda(v);
}

// ----------------------------------------------------------------------------
// Function createSuffixArray
// ----------------------------------------------------------------------------

// no callback implies that algo is discarded as well (TODO do we want that?)
template <typename TSA, typename TText, typename TAlgo>
inline void
createSuffixArray(TSA &SA,
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

// different Algos with callback
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgoSpec,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  SaAdvancedSort<TAlgoSpec> const &,
                  TLambda progressCallback = [] () {})
{
    typedef StringSet< TString, TSSetSpec > TText;
    typedef typename Size<TSA>::Type TSize;
    typedef typename Iterator<TSA, Standard>::Type TIter;

    // 1. Fill suffix array with a permutation (the identity)
    TIter it = begin(SA, Standard());
    for(unsigned j = 0; j < length(s); ++j)
    {
        TSize len = length(s[j]);
        for(TSize i = 0; i < len; ++i, ++it)
            *it = Pair<unsigned, TSize>(j, i);
    }

    // 2. Sort suffix array with algo
#ifdef GNUOMP
    typedef typename SaAdvancedSortAlgoTag<TAlgoSpec>::Type TAlgo;
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

// QuicksortBucket with callback
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TLambda>
inline void
createSuffixArray(TSA & sa,
                  StringSet<TString, TSSetSpec> const & text,
                  SaAdvancedSort<QuickSortBucketTag> const &,
                  TLambda progressCallback)
{
    typedef StringSet<TString, TSSetSpec>               TText;
//     typedef typename Size<TText>::Type                  TTextSize;
//     typedef uint64_t TTextSize;
//     typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Value<TString>::Type               TIndexAlphabet;
//     typedef typename Size<TText>::Type                     TIndexSize;
//     typedef typename Fibre<TIndex, FibreSA>::Type           TIndexSAFibre;
    typedef typename Value<TSA>::Type                   TIndexSAPos;
    typedef typename Value<TIndexSAPos, 1>::Type        TSASet;
    typedef typename Value<TIndexSAPos, 2>::Type        TSASeq;
    typedef typename Iterator<TSA, Standard>::Type      TIter;
//     typedef typename Iterator<TIndex, TopDown<> >::Type     TIterator;

    // we don't know the maximal seed length in advance as it depends on the query
    // hence we choose a sufficiently large number
//     unsigned maxSeedLength = 10000;

//     // 1. create and sort q-gram buckets
//     Shape<TIndexAlphabet, SimpleShape> shape;
//     String<TIndexSize> dir;
// 
//     unsigned shapeLength;
//     if (ValueSize<TIndexAlphabet>::VALUE <= 5)
//         shapeLength = _min(maxSeedLength, 10u);
//     else
//         shapeLength = _min(maxSeedLength, 3u);
// 
//     TTextSize stepSize = 1;
// //     if (IsSameType<TSpec, Exact>::VALUE)
// //         stepSize = maxSeedLength;
// 
//     resize(shape, shapeLength);
// //     resize(sa, _qgramQGramCount(text, shape, stepSize), Exact());
//     resize(dir, _fullDirLength(shape), Exact());
//     Nothing nothing;
// 
//     createQGramIndex(sa, dir, nothing, text, shape, stepSize);

    TSASeq initialSortLength;
    if (ValueSize<TIndexAlphabet>::VALUE <= 5)
        initialSortLength = 10u;
    else if (ValueSize<TIndexAlphabet>::VALUE < 10)
        initialSortLength = 3u;
    else
        initialSortLength = 2u;

    TIter it = begin(sa, Standard());
    TSASet numSeqs = length(text);
    for(TSASet j = 0; j < numSeqs; ++j)
    {
        TSASeq seqLen = length(text[j]);
        for(TSASeq i = 0; i < seqLen; ++i, ++it)
            *it = Pair<TSASet, TSASeq>(j, i);
    }

//     if (it != end(sa, Standard()))
//     {
//         std::cerr << "Not all SAVs written\n" 
//                   << "lengthSum(text) == " << lengthSum(text) << "\n"
//                   << "length(sa) == " << length(sa) << "\n"
//                   << "it - begin(sa) == " << (it - begin(sa, Standard())) << "\n"
//                   << std::endl;
//     }
//     std::cout << "typeof(TSize) == " << typeid(TTextSize).name() << std::endl;
    // sort up to certain depth
#ifdef GNUOMP
    __gnu_parallel::sort(begin(sa, Standard()),
                         end(sa, Standard()),
                         QGramLess_<TIndexSAPos, TText const>(text,
                                                              initialSortLength),
                         __gnu_parallel::quicksort_tag());
#else
    std::sort(begin(sa, Standard()),
              end(sa, Standard()),
              QGramLess_<TIndexSAPos, TText const>(text, initialSortLength));
#endif

//     std::cout << "First 20 sa" << std::endl;
//     for (unsigned i = 0; i < 20; ++i)
//         std::cout << suffix(text[sa[i].i1], sa[i].i2) << "\n";
    // create dir
    String<uint64_t> dir;
    appendValue(dir, 0u);
    for(uint64_t j = 1; j < length(sa); ++j)
    {
        if (infix(text[sa[j].i1],
                  sa[j].i2,
                  std::min(static_cast<TSASeq>(sa[j].i2 + initialSortLength),
                           static_cast<TSASeq>(length(text[sa[j].i1])))) !=
            infix(text[sa[j-1].i1],
                  sa[j-1].i2,
                  std::min(static_cast<TSASeq>(sa[j-1].i2 + initialSortLength),
                           static_cast<TSASeq>(length(text[sa[j-1].i1])))))
        {
//             if (j < 50)
//             {
//                 std::cout << "last: " << infix(text[sa[j-1].i1],
//                   sa[j-1].i2,
//                   sa[j-1].i2 + std::min(initialSortLength, length(text[sa[j-1].i1])))
//                           << "\ncur:  " << infix(text[sa[j].i1],
//                                                sa[j].i2,
//                 sa[j].i2 + std::min(initialSortLength, length(text[sa[j].i1])))
//                           << std::endl;
//             }
            appendValue(dir, j);
        }
    }
    appendValue(dir, length(sa));

//     std::cout << "First 20 sa[dir]" << std::endl;
//     for (unsigned i = 0; i < 20; ++i)
//         std::cout << suffix(text[sa[dir[i]].i1], sa[dir[i]].i2) << "\n";


    // 2. refine q-gram buckets and sort up to their maxSeedLength prefix

    TIter saBegin = begin(sa, Standard());

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic,1))
    for (uint64_t i = 1; i < length(dir); ++i)
    {
        if (dir[i - 1] + 1 < dir[i])
#ifdef GNUOMP
            __gnu_parallel::sort(
                saBegin + dir[i - 1],
                saBegin + dir[i],
                SuffixLess_<TIndexSAPos, TText const>(text, initialSortLength),
                __gnu_parallel::sequential_tag()); // deactivate internal SMP
#else
            std::sort(
                saBegin + dir[i - 1],
                saBegin + dir[i],
                SuffixLess_<TIndexSAPos, TText const>(text, initialSortLength));
#endif

        progressCallbackWrapper(progressCallback,
                                uint64_t(i * 100 / length(dir)));
    }

//     std::cout << "POST REFINE First 20 sa" << std::endl;
//     for (unsigned i = 0; i < 20; ++i)
//         std::cout << suffix(text[sa[i].i1], sa[i].i2) << "\n";
//     std::cout << "POST REFINE First 20 sa[dir]" << std::endl;
//     for (unsigned i = 0; i < 20; ++i)
//         std::cout << suffix(text[sa[dir[i]].i1], sa[dir[i]].i2) << "\n";
}

// general case discards the callback
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgo,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  TAlgo const &,
                  TLambda const &)
{
    return createSuffixArray(SA, s, TAlgo());
}

// ----------------------------------------------------------------------------
// Function indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TLambda>
inline bool indexCreate(Index<TText, FMIndex<TSpec, TConfig> > & index,
                        TText const & text,
                        FibreSALF const &,
                        TLambda const & progressCallback)
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >      TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type   TTempSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;
    typedef typename Size<TIndex>::Type                 TSize;

    indexText(index) = text;

    if (empty(text))
        return false;

    TTempSA tempSA;

    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    createSuffixArray(tempSA, text, TAlgo(), progressCallback);

    // Create the LF table.
    createLF(indexLF(index), text, tempSA);

    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());

    // Create the compressed SA.
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);

    return true;
}

template <typename TText, typename TSpec, typename TLambda>
inline bool indexCreate(Index<TText, IndexSa<TSpec> > & index,
                        TText const & text,
                        FibreSA const &,
                        TLambda const & progressCallback)
{
    typedef Index<TText, IndexSa<TSpec> >      TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type       TSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    indexText(index) = text;

    if (empty(text))
        return false;

    TSA & sa = getFibre(index, FibreSA());

    // Create the full SA.
    resize(sa, lengthSum(text), Exact());
    createSuffixArray(sa, text, TAlgo(), progressCallback);

    return true;
}

} // namespace seqan

#endif
