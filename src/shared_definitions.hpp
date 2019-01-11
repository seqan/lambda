// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2019, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2019, Knut Reinert and Freie Universität Berlin
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
// options.h: contains the options and argument parser
// ==========================================================================

#ifndef SEQAN_SHARED_DEFINITIONS_H_
#define SEQAN_SHARED_DEFINITIONS_H_

#include <seqan/index.h>
#include <seqan/blast.h>

using namespace seqan;

// ==========================================================================
// Metafunctions
// ==========================================================================

// SIZE TYPES

// Expected Number of Sequences
template <typename TAlph>
using SizeTypeNum_ = uint32_t;

// Expected Lengths of Sequences
template <typename T>
struct SizeTypePosMeta_
{
#ifdef LAMBDA_LONG_PROTEIN_SUBJ_SEQS
    using Type = uint32_t;
#else
    using Type = uint16_t;
#endif
};

template <>
struct SizeTypePosMeta_<Dna5>
{
    // DNA sequences are expected to be longer
    using Type = uint32_t;
};

template <typename TAlph>
using SizeTypePos_ = typename SizeTypePosMeta_<TAlph>::Type;


// suffix array overloads
namespace seqan
{

template<typename TSpec1, typename TSpec2, typename TSpec3>
struct SAValue<StringSet<String<ReducedAminoAcid<TSpec1>, TSpec2>, TSpec3> >
{
    typedef Pair<SizeTypeNum_<TSpec1>, SizeTypePos_<TSpec1>, Pack> Type;
};

template<typename TSpec1, typename TSpec2, typename TSpec3, typename TFunctor>
struct SAValue<StringSet<ModifiedString<String<TSpec1, TSpec2>, TFunctor>, TSpec3> >
{
    typedef Pair<SizeTypeNum_<TSpec1>, SizeTypePos_<TSpec1>, Pack> Type;
};

template<typename TSpec1, typename TSpec2, typename TSpec3, typename TFunctor, typename TFunctor2>
struct SAValue<StringSet<ModifiedString<ModifiedString<String<TSpec1, TSpec2>, TFunctor>, TFunctor2>, TSpec3> >
{
    typedef Pair<SizeTypeNum_<TSpec1>, SizeTypePos_<TSpec1>, Pack> Type;
};

template<typename TSpec1, typename TSpec2, typename TSpec3>
struct SAValue<StringSet<String<TSpec1, TSpec2>, TSpec3> >
{
    typedef Pair<SizeTypeNum_<TSpec1>, SizeTypePos_<TSpec1>, Pack> Type;
};

template <typename TString, typename TSpec>
struct DefaultIndexStringSpec<StringSet<TString, TSpec>>
{
#if !defined(LAMBDA_INDEXER) && defined(LAMBDA_MMAPPED_DB)
    using Type    = MMap<>;
#else
    using Type    = Alloc<>;
#endif
};

// our custom Bam Overload
template <typename TDirection, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Bam, TDirection, BlastTabular>, TStorageSpec>
{
    typedef typename DefaultIndexStringSpec<StringSet<void, void>>::Type TStringSpec; // see above
    typedef StringSet<Segment<String<char, TStringSpec>, InfixSegment> > TNameStore;
    typedef NameStoreCache<TNameStore>                                   TNameStoreCache;
    typedef BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>      Type;
};

}

// Index Specs
struct LambdaFMIndexConfig
{
    using LengthSum = size_t;
#if !defined(LAMBDA_INDEXER) && defined(LAMBDA_MMAPPED_DB)
    using TAlloc    = MMap<>;
#else
    using TAlloc    = Alloc<>;
#endif

    using Bwt       = Levels<void, LevelsRDConfig<LengthSum, TAlloc, 3, 3> >;
    using Sentinels = Levels<void, LevelsRDConfig<LengthSum, TAlloc> >;

    static const unsigned SAMPLING = 10;
};

struct LambdaFMIndexConfigInBi : LambdaFMIndexConfig
{
    using Bwt       = Levels<void, LevelsPrefixRDConfig<LengthSum, TAlloc, 3, 3> >;
};

template <typename TSpec = void>
using TFMIndex = FMIndex<TSpec, LambdaFMIndexConfig>;

template <typename TSpec = void>
using TFMIndexInBi = FMIndex<TSpec, LambdaFMIndexConfigInBi>;

// lazy...
template <typename TString>
using TCDStringSet = StringSet<TString, Owner<ConcatDirect<> > >;

template <BlastProgram p>
using OrigQryAlph = typename std::conditional<
                                           (p == BlastProgram::BLASTN) ||
                                           (p == BlastProgram::BLASTX) ||
                                           (p == BlastProgram::TBLASTX),
                                           Dna5,
                                           AminoAcid>::type;

template <BlastProgram p>
using OrigSubjAlph = typename std::conditional<
                                           (p == BlastProgram::BLASTN) ||
                                           (p == BlastProgram::TBLASTN) ||
                                           (p == BlastProgram::TBLASTX),
                                           Dna5,
                                           AminoAcid>::type;

template <BlastProgram p>
using TransAlph = typename std::conditional<(p == BlastProgram::BLASTN),
                                            Dna5,
                                            AminoAcid>::type;


template <BlastProgram p, typename TRedAlph_>
using RedAlph = typename std::conditional<(p == BlastProgram::BLASTN),
                                          Dna5,
                                          TRedAlph_>::type;



// ==========================================================================
// Global variables
// ==========================================================================

// this is increased after incompatible changes to on-disk format
constexpr uint64_t indexGeneration = 1;

#endif // header guard
