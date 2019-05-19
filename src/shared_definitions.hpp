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

#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

#include "shared_options.hpp"

// ==========================================================================
// Global variables
// ==========================================================================

// this is increased after incompatible changes to on-disk format
inline constexpr uint64_t currentIndexGeneration = 0;

// ==========================================================================
// Metafunctions
// ==========================================================================

// SIZE TYPES
#if 0
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
#endif

template <typename TText = void>
using TFMIndex = seqan3::fm_index<TText>;

template <typename TText = void>
using TFMIndexInBi = seqan3::bi_fm_index<TText>;

// lazy...
template <typename TString>
using TCDStringSet = std::vector<TString>; //TODO seqan3::concatenated_sequences


template <DbIndexType           dbIndexType,
          AlphabetEnum          origAlph,
          AlphabetEnum          transAlph,
          AlphabetEnum          redAlph>    // <- all members of index_file_options that influence types
struct index_file
{
    index_file_options options{};

    TCDStringSet<std::string>                                   ids;
    std::vector<uint64_t>                                       origSeqLengths; // only used when origAlph != transAlph
    TCDStringSet<std::vector<_alphabetEnumToType<transAlph>>>   transSeqs;
    std::vector<std::vector<uint32_t>>                          sTaxIds; //TODO double check int-width

    std::vector<uint32_t>                                       taxonParentIDs;
    std::vector<uint8_t>                                        taxonHeights;
    std::vector<std::string>                                    taxonNames;

    TCDStringSet<std::vector<_alphabetEnumToType<redAlph>>>     redSeqs; // temporary
    std::conditional_t<dbIndexType == DbIndexType::BI_FM_INDEX,
                       seqan3::bi_fm_index<decltype(redSeqs)>,
                       seqan3::fm_index<decltype(redSeqs)>>     index;

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options",          options),
                cereal::make_nvp("ids",              ids),
                cereal::make_nvp("origSeqLengths",   origSeqLengths),
                cereal::make_nvp("transSeqs",        transSeqs),
                cereal::make_nvp("sTaxIds",          sTaxIds),
                cereal::make_nvp("taxonParentIDs",   taxonParentIDs),
                cereal::make_nvp("taxonHeights",     taxonHeights),
                cereal::make_nvp("taxonNames",       taxonNames),
                cereal::make_nvp("redSeqs",          redSeqs),
                cereal::make_nvp("index",            index));
    }
};


struct fake_index_file
{
    index_file_options & options;

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options",          options);
    }
};

#endif // header guard
