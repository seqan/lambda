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

#include "radix_inplace.h"

namespace seqan
{

// ==========================================================================
// Tags
// ==========================================================================

struct RadixSortSACreateTag {};

// ==========================================================================
// Forwards
// ==========================================================================

} // namespace seqan

// ==========================================================================
// Metafunctions
// ==========================================================================

namespace SEQAN_NAMESPACE_MAIN
{
template <typename TText>
struct Fibre<Index<TText, IndexSa<RadixSortSACreateTag > >, FibreTempSA>
{
    typedef Index<TText, IndexSa<RadixSortSACreateTag > >          TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template <typename TText>
struct DefaultIndexCreator<Index<TText, IndexSa<RadixSortSACreateTag> >, FibreSA>
{
    typedef RadixSortSACreateTag Type;
};

template <typename TText, typename TConfig>
struct Fibre<Index<TText, FMIndex<RadixSortSACreateTag, TConfig> >, FibreTempSA>
{
    typedef Index<TText, FMIndex<RadixSortSACreateTag, TConfig> >  TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template < typename TText, typename TConfig>
struct DefaultIndexCreator<Index<TText, FMIndex<RadixSortSACreateTag, TConfig> >, FibreSA>
{
    typedef RadixSortSACreateTag Type;
};

} // SEQAN_NAMESPACE_MAIN

// ==========================================================================
// Classes
// ==========================================================================

namespace seqan
{

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Helper Functions
// ----------------------------------------------------------------------------

inline std::function<void(unsigned)>
progressWrapper2(std::function<void(void)> const &)
{
    return [] (unsigned) {};
}

inline std::function<void(unsigned)>
progressWrapper2(std::function<void(unsigned)> const & l)
{
    return l;
}

// ----------------------------------------------------------------------------
// Function createSuffixArray
// ----------------------------------------------------------------------------

template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  RadixSortSACreateTag const &,
                  TLambda const & progressCallback = [] (unsigned) {})
{
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

    // 2. Sort suffix array with inplace radix Sort
//     std::cerr << "I am really in radix sort" << std::endl;
    inPlaceRadixSort(SA, s, progressWrapper2(progressCallback));
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
