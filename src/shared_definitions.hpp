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

#pragma once

#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/deep.hpp>

#include "shared_options.hpp"
#include "view_dna_n_to_random.hpp"

// ==========================================================================
//  Global variables
// ==========================================================================

// this is increased after incompatible changes to on-disk format
inline constexpr uint64_t currentIndexGeneration = 0;

// ==========================================================================
//  Index specs TODO experiment with these
// ==========================================================================

/* default */
template <size_t>
using IndexSpec = seqan3::default_sdsl_index_type;

/* default, spelled out */
// template <size_t>
// using IndexSpec = sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
//                                               sdsl::rank_support_v<>,
//                                               sdsl::select_support_scan<>,
//                                               sdsl::select_support_scan<0> >,
//                                 16,
//                                 10000000,
//                                 sdsl::sa_order_sa_sampling<>,
//                                 sdsl::isa_sampling<>,
//                                 sdsl::plain_byte_alphabet>;

/* huffman, doesn't satisfy concept at the moment */
// template <size_t>
// using IndexSpec = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector,
//                                               sdsl::rank_support_v<>,
//                                               sdsl::select_support_scan<>,
//                                               sdsl::select_support_scan<0> >,
//                                 16,
//                                 10000000,
//                                 sdsl::sa_order_sa_sampling<>,
//                                 sdsl::isa_sampling<>,
//                                 sdsl::plain_byte_alphabet>;

/* epr, in experimental branch */
// template <size_t alph_size>
// using IndexSpec = sdsl::csa_wt<sdsl::wt_epr<alph_size + 2>, // +1 for sentinels, +1 for collection
//                                 16,
//                                 1'0000'000,
//                                 sdsl::sa_order_sa_sampling<>,
//                                 sdsl::isa_sampling<>,
//                                 sdsl::plain_byte_alphabet>;

// ==========================================================================
//  Misc. aliases
// ==========================================================================

template <typename TString>
using TCDStringSet = std::vector<TString>; //TODO seqan3::concatenated_sequences

template <typename TSpec>
using TTransAlphModString =
  decltype(std::declval<TSpec &>() | seqan3::views::translate_join);

template <typename TSpec, typename TAlph>
using TRedAlphModString =
  decltype(std::declval<TSpec &>() | seqan3::views::deep{seqan3::views::convert<TAlph>});

template <typename TSpec, typename TAlph>
using TRedNuclAlphModString =
  decltype(std::declval<TSpec &>() | seqan3::views::dna_n_to_random);

// ==========================================================================
//  The index
// ==========================================================================

/* actual index */
template <DbIndexType           dbIndexType,
          AlphabetEnum          origAlph,
          AlphabetEnum          redAlph>    // <- all members of index_file_options that influence types
struct index_file
{
    index_file_options options{};

    TCDStringSet<std::string>                                   ids;
    TCDStringSet<std::vector<_alphabetEnumToType<origAlph>>>    seqs;
    std::vector<std::vector<uint32_t>>                          sTaxIds; //TODO double check int-width

    std::vector<uint32_t>                                       taxonParentIDs;
    std::vector<uint8_t>                                        taxonHeights;
    std::vector<std::string>                                    taxonNames; //TODO TCDStringSet?

    using TRedAlph      = _alphabetEnumToType<redAlph>;
    using TIndexSpec    = IndexSpec<seqan3::alphabet_size<TRedAlph>>;
    using TIndex        = std::conditional_t<dbIndexType == DbIndexType::BI_FM_INDEX,
        seqan3::bi_fm_index<TRedAlph, seqan3::text_layout::collection, TIndexSpec>,
        seqan3::fm_index<TRedAlph, seqan3::text_layout::collection, TIndexSpec>>;

    TIndex index;

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options",          options),
                cereal::make_nvp("ids",              ids),
                cereal::make_nvp("seqs",             seqs),
                cereal::make_nvp("sTaxIds",          sTaxIds),
                cereal::make_nvp("taxonParentIDs",   taxonParentIDs),
                cereal::make_nvp("taxonHeights",     taxonHeights),
                cereal::make_nvp("taxonNames",       taxonNames),
                cereal::make_nvp("index",            index));
    }
};

/* pseudo-index that can be used to only load the options from disk */
struct fake_index_file
{
    index_file_options & options;

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options",          options));
    }
};

// ==========================================================================
//  Misc. functions
// ==========================================================================

template <typename TTargetAlph, typename TRange, typename TAdaptProt, typename TAdaptNucl>
    requires std::same_as<seqan3::innermost_value_type_t<TRange>, TTargetAlph>
TRange & initHelper(TRange & input, TAdaptProt &&, TAdaptNucl &&)
{
    return input;
}

template <typename TTargetAlph, typename TRange, typename TAdaptProt, typename TAdaptNucl>
    requires seqan3::nucleotide_alphabet<TTargetAlph> &&
             (!std::same_as<seqan3::innermost_value_type_t<TRange>, TTargetAlph>)
auto initHelper(TRange & input, TAdaptProt &&, TAdaptNucl && adaptNucl)
{
    return input | std::forward<TAdaptNucl>(adaptNucl);
}

template <typename TTargetAlph, typename TRange, typename TAdaptProt, typename TAdaptNucl>
    requires !std::same_as<seqan3::innermost_value_type_t<TRange>, TTargetAlph>
auto initHelper(TRange & input, TAdaptProt && adaptProt, TAdaptNucl &&)
{
    return input | std::forward<TAdaptProt>(adaptProt);
}


