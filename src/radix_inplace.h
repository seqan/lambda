// ==========================================================================
//                              radix_inplace.h
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
// Author: Sascha Meiers <sascha.meiers@embl.de>
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
// ==========================================================================
// The Radix Sort functions are adapted from Martin Frith's "last"
// tool (last.cbrc.jp), but he himself adapted the code from McIlroy, Bostic:
// "Engineering radix sort" as well as Karkkainen, Rantala: "Engineering radix
// sort for strings". Thanks to Martin for showing this to me.
// ============================================================================

#ifndef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_
#define CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_

namespace SEQAN_NAMESPACE_MAIN
{

// ============================================================================
// struct RadixTextAccessor
// ============================================================================

template <
    typename TSAValue,          // input
    typename TString,           // string object that is referenced
    typename TSpec = void,      // Suffix modifier
    typename TSize = unsigned>  // return type (ordValue)
struct RadixTextAccessor;
/*
 * NOTE:
 * These accessors cannot resolve the correct order of out-of-bound-positions,
 * i.e. when suffixes are equal up to their last character.
 * All these cases get collected in a 0 bucket.
 * The InplaceRadixSorter takes care of that by calling a special
 * sort function on the 0 buckets.
 */

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                                            [String]
// ----------------------------------------------------------------------------

template <typename TSAValue, typename TString, typename TSize>
struct RadixTextAccessor<TSAValue, TString, void, TSize> :
    public std::unary_function<TSAValue, TSize>
{
    TString const & text;
    typename Size<TString>::Type const L;

    RadixTextAccessor(TString const &str) : text(str), L(length(str))
    {}

    template <typename TSize2>
    inline TSize operator()(TSAValue const &x, TSize2 depth) const
    {
        typename Size<TString>::Type pos = x + depth;
        if (pos >= L)   return 0;
        TSize ret = ordValue(text[pos]);
        return ret+1;
    }
};

// ----------------------------------------------------------------------------
// struct RadixTextAccessor                                         [StringSet]
// ----------------------------------------------------------------------------

template <typename TSAValue, typename TString, typename TSetSpec, typename TSize>
struct RadixTextAccessor<TSAValue, StringSet<TString, TSetSpec>, void, TSize> :
    public std::unary_function<TSAValue, TSize>

{
    StringSet<TString, TSetSpec> const & text;
    String<typename Size<TString>::Type> L;

    RadixTextAccessor(StringSet<TString, TSetSpec> const &str) : text(str)
    {
        resize(L, length(text), Exact());
        for(typename Size<TString>::Type i = 0; i < length(text); ++i)
        L[i] = length(text[i]);
    }

    template <typename TSize2>
    inline TSize operator()(TSAValue const &x, TSize2 depth) const
    {
        typename Size<TString>::Type pos = getSeqOffset(x) + depth;
        typename Size<TString>::Type seq = getSeqNo(x);
        if (pos >= L[seq])   return 0;
        TSize ret = ordValue(text[seq][pos]);
        return ret+1;
    }
};

// ----------------------------------------------------------------------------
// struct _ZeroBucketComparator                                    [StringSet]
// ----------------------------------------------------------------------------
// Functors to compare suffixes from 0 bucket (suffixes that are lex. equal)
// ----------------------------------------------------------------------------

template <typename TSAValue, typename TLimitsString=Nothing const>
struct _ZeroBucketComparator
{
    TLimitsString const & limits;
    _ZeroBucketComparator(TLimitsString const & lim) : limits(lim)  { /*std::cout << "limits: " << limits << std::endl;*/   }

    inline bool operator()(TSAValue const & a, TSAValue const & b) const
    {
        typename Size<TLimitsString>::Type lena = limits[getSeqNo(a)+1]-limits[getSeqNo(a)] - getSeqOffset(a);
        typename Size<TLimitsString>::Type lenb = limits[getSeqNo(b)+1]-limits[getSeqNo(b)] - getSeqOffset(b);  
        if (lena == lenb)
            return getSeqNo(a) > getSeqNo(b);
        else
            return lena < lenb;
    }
};

// ----------------------------------------------------------------------------
// struct _ZeroBucketComparator                                       [String]
// ----------------------------------------------------------------------------

template <typename TSAValue>
struct _ZeroBucketComparator<TSAValue, Nothing const>
{
    _ZeroBucketComparator(Nothing const &) {}
    _ZeroBucketComparator(Nothing &) {}


    inline bool operator()(TSAValue const & a, TSAValue const & b) const
    {
        return a > b;
    }
};

// ----------------------------------------------------------------------------
// struct RadixSortContext_
// ----------------------------------------------------------------------------

template <typename TAccessFunctor,                // text accessor
          typename TOrderFunctor,                 // For seperate sort of the 0 bucket.
          typename TSize,                         // type of depth and bucketCount a.s.o
          unsigned Q>                             // alph size = ValueSize + 1
struct RadixSortContext_
{
    static_assert(Q < 256, "Alphabet size must be smaller 256!"); //TODO really?
    typedef typename TAccessFunctor::argument_type      TSAValue;
    typedef typename TAccessFunctor::result_type        TOrdValue; // unsigned

    static const unsigned ORACLESIZE = 256;
    TAccessFunctor     textAccess;
    TOrderFunctor      comp;

    TSize bucketSize[Q];
    std::array<TSAValue*,Q> bucketEnd;

    RadixSortContext_(TAccessFunctor const & f, TOrderFunctor const & c) :
        textAccess(f), comp(c)
    {}
};

template <typename TAccessFunctor,
          typename TOrderFunctor,
          typename TSize,
          unsigned Q>
inline void
clear(RadixSortContext_<TAccessFunctor, TOrderFunctor, TSize, Q> & context)
{
    memset(context.bucketSize, 0, sizeof(TSize)*Q);
}

// ----------------------------------------------------------------------------
// Function radixSort()
// ----------------------------------------------------------------------------

template <typename TSAValue, typename TSize,
          typename TAccessFunctor, typename TOrderFunctor, unsigned Q>
inline void
radixSort(std::vector<std::tuple<TSAValue*, TSAValue*, TSize> > & stack,
          RadixSortContext_<TAccessFunctor, TOrderFunctor, TSize, Q> & context,
          std::tuple<TSAValue*, TSAValue*, TSize> const & item)
{
    typedef RadixSortContext_<TAccessFunctor, TOrderFunctor, TSize, Q> TContext;
    typedef typename TContext::TOrdValue TOrdValue;
    static_assert(std::is_same<TSAValue, typename TContext::TSAValue>::value, "TSAValue mismatch!");

    clear(context);

    // get bucket sizes (i.e. letter counts):
    // The intermediate oracle array makes it faster (see "Engineering
    // Radix Sort for Strings" by J Karkkainen & T Rantala)
    for(TSAValue* i = std::get<0>(item); i < std::get<1>(item); /* noop */ )
    {
        // buffer for the next chars
        TOrdValue oracle [TContext::ORACLESIZE];
        TOrdValue* oracleEnd = oracle + std::min(static_cast<std::size_t>(TContext::ORACLESIZE),
                                                 static_cast<std::size_t>(std::get<1>(item) - i));

        for(TOrdValue* j = oracle; j < oracleEnd; ++j )
            *j = context.textAccess(*i++, std::get<2>(item));

        for(TOrdValue* j = oracle; j < oracleEnd; ++j )
            ++context.bucketSize[*j];
    }

    // get bucket std::get<1>(item)s, and put buckets on the stack to sort within them later:
    // EDIT: 0 bucket is not sorted here, but later.
    TSize zeroBucketSize    = context.bucketSize[0];
    TSAValue* pos           = std::get<0>(item) + context.bucketSize[0];
    context.bucketEnd[0]    = pos;

    for(unsigned i = 1; i < Q; ++i )
    {
        TSAValue* nextPos = pos + context.bucketSize[i];
        if (nextPos - pos > 1)
            stack.emplace_back(pos, nextPos, std::get<2>(item)+1);

        pos = nextPos;
        context.bucketEnd[i] = pos;
    }

    // permute items into the correct buckets:
    for(TSAValue* i = std::get<0>(item); i < std::get<1>(item); )
    {
        TOrdValue subset;  // unsigned is faster than uchar!
        TSAValue holdOut = *i;
        while(--context.bucketEnd[subset = context.textAccess(holdOut, std::get<2>(item))] > i )
            std::swap(*context.bucketEnd[subset], holdOut);
        *i = holdOut;
        i += context.bucketSize[subset];
        context.bucketSize[subset] = 0;  // reset it so we can reuse it //TODO check if we need this, since we clear already!
    }

    // sort the 0 bucket using std::sort
    if(zeroBucketSize > 1)
        std::sort(std::get<0>(item), std::get<0>(item) + zeroBucketSize, context.comp);
}

// ----------------------------------------------------------------------------
// Function inplaceFullRadixSort()                                    [default]
// ----------------------------------------------------------------------------

//TODO: play with this value
#define _RADIX_SORT_SWITCH_TO_QUICKSORT_AT 100

#ifdef _OPENMP
#define N_THREADS omp_get_max_threads()
#define I_THREAD omp_get_thread_num()
#else
#define N_THREADS 1
#define I_THREAD 0
#endif

// TODO: serial version
// TODO: possibly do multiple runs of "secondStep" if alphabet size to small
// TODO: possibly quicksort directly on buckets in third steps, if buckets have been made small enough

template <typename TSA, typename TText, typename TLambda>
void inPlaceRadixSort(TSA & sa, TText const & str, TLambda const & progressCallback = [] (unsigned) {})
{
    typedef typename Value<typename Concatenator<TText>::Type>::Type TAlphabet;
    typedef typename Value<TSA>::Type                             TSAValue;
    typedef typename Size<TText>::Type                            TSize;
    typedef std::tuple<TSAValue*, TSAValue*, TSize>               TItem;
    typedef typename StringSetLimits<TText const>::Type           TLimitsString; // "Nothing" for Strings

    typedef RadixTextAccessor<TSAValue, TText>                    TAccessor;
    typedef _ZeroBucketComparator<TSAValue,TLimitsString>         TZeroComp;

    static const unsigned SIGMA = static_cast<unsigned>(ValueSize<TAlphabet>::VALUE) + 1;
    SEQAN_ASSERT_LT_MSG(SIGMA, 1000u, "Attention: inplace radix sort is not suited for large alphabets");

    typedef RadixSortContext_<TAccessor, TZeroComp, TSize, SIGMA>   TContext;

    if (empty(sa)) return; // otherwise access sa[0] fails

    /* stacks */
    std::vector<TItem> firstStack;
    firstStack.reserve(SIGMA);
    std::vector<TItem> secondStack;
    secondStack.reserve(SIGMA*SIGMA);
    std::vector<std::vector<TItem>> lStack(N_THREADS); // one per thread
    // reduce memory allocations in threads by reserving space
    for (auto & stack : lStack)
        stack.reserve(length(sa) / 1000);

    /* contexts */
    TContext firstSecondContext{TAccessor(str), TZeroComp(stringSetLimits(str))};
    std::vector<TContext> lContext(N_THREADS, TContext{TAccessor(str), TZeroComp(stringSetLimits(str))});

    // sort by the first character
    radixSort(firstStack, firstSecondContext, TItem(&sa[0], &sa[0]+length(sa), 0));

    progressCallback(5); // 5% progress guess after first char

    // sort by second character
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (unsigned j = 0; j < length(firstStack); ++j)
    {
        TItem & i = firstStack[j];
        if (std::get<1>(i) - std::get<0>(i) < _RADIX_SORT_SWITCH_TO_QUICKSORT_AT)
            std::sort(std::get<0>(i),
                      std::get<1>(i),
                      SuffixLess_<TSAValue, TText const>(str, std::get<2>(i)));
        else if (std::get<1>(i) - std::get<0>(i) >= 2)
            radixSort(lStack[I_THREAD], lContext[I_THREAD], i);
    }

    // merge local stacks and clear for next round
    for (auto & stack : lStack)
    {
        secondStack.insert(secondStack.end(), stack.begin(), stack.end());
        stack.clear();
    }

    progressCallback(10); // 10% progress guess after second char

    // sort rest
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (unsigned j = 0; j < secondStack.size(); ++j)
    {
        lStack[I_THREAD].push_back(secondStack[j]);

        while (!lStack[I_THREAD].empty())
        {
            TItem i = lStack[I_THREAD].back();
            lStack[I_THREAD].pop_back();

            if (std::get<1>(i) - std::get<0>(i) < _RADIX_SORT_SWITCH_TO_QUICKSORT_AT)
                std::sort(std::get<0>(i),
                          std::get<1>(i),
                          SuffixLess_<TSAValue, TText const>(str, std::get<2>(i)));
            else if (std::get<1>(i) - std::get<0>(i) >= 2)
                radixSort(lStack[I_THREAD], lContext[I_THREAD], i);
        }

        // progressCallback must be thread safe and cope with smaller numbers after big numbers
        // remaining characters alloted 90% of total progress
        progressCallback(10 + (j * 90) / secondStack.size());
    }

    progressCallback(100); // done
}

}

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_RADIX_INPLACE_H_
