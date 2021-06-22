// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2020, Hannes Hauswedell <h2 @ fsfe.org>
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
// shared_definitions.h: commonly used types and functions
// ==========================================================================

#pragma once

#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>
#include <seqan3/alphabet/aminoacid/aa10li.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/alphabet/views/translate_join.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/deep.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

#include "view_add_reverse_complement.hpp"
#include "view_dna_n_to_random.hpp"
#include "view_duplicate.hpp"
#include "view_reduce_to_bisulfite.hpp"

// ==========================================================================
//  DbIndexType
// ==========================================================================

enum class DbIndexType : uint8_t
{
    FM_INDEX,
    BI_FM_INDEX
};

inline std::string
_indexEnumToName(DbIndexType const t)
{
    switch (t)
    {
        case DbIndexType::FM_INDEX:      return "fm_index";
        case DbIndexType::BI_FM_INDEX:   return "bi_fm_index";
    }

    throw std::runtime_error("Error: unknown index type");
    return "";
}

inline DbIndexType
_indexNameToEnum(std::string const t)
{
    if (t == "bi_fm_index")
        return DbIndexType::BI_FM_INDEX;
    else if (t == "fm_index")
        return DbIndexType::FM_INDEX;

    throw std::runtime_error("Error: unknown index type");
    return DbIndexType::FM_INDEX;
}

// ==========================================================================
//  Alphabet stuff
// ==========================================================================

constexpr const char *
_alphTypeToName(seqan3::semialphabet_any<6> const & /**/)
{
    return "dna3bs";
}

constexpr const char *
_alphTypeToName(seqan3::dna4 const & /**/)
{
    return "dna4";
}

constexpr const char *
_alphTypeToName(seqan3::dna5 const & /**/)
{
    return "dna5";
}

constexpr const char *
_alphTypeToName(seqan3::aa27 const & /**/)
{
    return "aminoacid";
}

constexpr const char *
_alphTypeToName(seqan3::aa10murphy const & /**/)
{
    return "murphy10";
}

constexpr const char *
_alphTypeToName(seqan3::aa10li const & /**/)
{
    return "li10";
}

enum class AlphabetEnum : uint8_t
{
    UNDEFINED,
    DNA3BS,
    DNA4,
    DNA5,
    AMINO_ACID,
    MURPHY10,
    LI10,
};

inline std::string
_alphabetEnumToName(AlphabetEnum const t)
{
    switch (t)
    {
        case AlphabetEnum::UNDEFINED:   return "UNDEFINED";
        case AlphabetEnum::DNA3BS:      return _alphTypeToName(seqan3::semialphabet_any<6>{});
        case AlphabetEnum::DNA4:        return _alphTypeToName(seqan3::dna4{});
        case AlphabetEnum::DNA5:        return _alphTypeToName(seqan3::dna5{});
        case AlphabetEnum::AMINO_ACID:  return _alphTypeToName(seqan3::aa27{});
        case AlphabetEnum::MURPHY10:    return _alphTypeToName(seqan3::aa10murphy{});
        case AlphabetEnum::LI10:        return _alphTypeToName(seqan3::aa10li{});
    }

    throw std::runtime_error("Error: unknown alphabet type");
    return "";
}

inline AlphabetEnum
_alphabetNameToEnum(std::string const t)
{
    if ((t == "UNDEFINED") || (t == "auto"))
        return AlphabetEnum::UNDEFINED;
    else if (t == _alphTypeToName(seqan3::semialphabet_any<6>{}))
        return AlphabetEnum::DNA3BS;
    else if (t == _alphTypeToName(seqan3::dna4{}))
        return AlphabetEnum::DNA4;
    else if (t == _alphTypeToName(seqan3::dna5{}))
        return AlphabetEnum::DNA5;
    else if (t == _alphTypeToName(seqan3::aa27{}))
        return AlphabetEnum::AMINO_ACID;
    else if (t == _alphTypeToName(seqan3::aa10murphy{}))
        return AlphabetEnum::MURPHY10;
    else if (t == _alphTypeToName(seqan3::aa10li{}))
        return AlphabetEnum::LI10;

    throw std::runtime_error("Error: unknown alphabet type");
    return AlphabetEnum::DNA4;
}

template <AlphabetEnum e>
struct _alphabetEnumToType_;

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA3BS>
{
    using type = seqan3::semialphabet_any<6>;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA4>
{
    using type = seqan3::dna4;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA5>
{
    using type = seqan3::dna5;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::AMINO_ACID>
{
    using type = seqan3::aa27;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::MURPHY10>
{
    using type = seqan3::aa10murphy;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::LI10>
{
    using type = seqan3::aa10li;
};

template <AlphabetEnum e>
using _alphabetEnumToType = typename _alphabetEnumToType_<e>::type;

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

// Translated subject sequences
// Case: DNA to amino acid translation
template <typename TOrigSeqs, typename TOrigAlph, typename TTransAlph, typename TRedAlph>
struct TTransSbjSeqsImpl
{
    using type = decltype(std::declval<TOrigSeqs&>() | seqan3::views::translate_join);
};

// Case: bisulfite mode
template <typename TOrigSeqs, typename TSameAlph>
struct TTransSbjSeqsImpl<TOrigSeqs, TSameAlph, TSameAlph, seqan3::semialphabet_any<6>>
{
    using type = decltype(std::declval<TOrigSeqs&>() | views::duplicate);
};

// Case: no translation
template <typename TOrigSeqs, typename TSameAlph, typename TRedAlph>
struct TTransSbjSeqsImpl<TOrigSeqs, TSameAlph, TSameAlph, TRedAlph>
{
    using type = decltype(seqan3::views::type_reduce(std::declval<TOrigSeqs &>()));
};

// Translated query sequences
// Case: DNA to amino acid translation
template <typename TOrigSeqs, typename TOrigAlph, typename TTransAlph, typename TRedAlph>
struct TTransQrySeqsImpl
{
    using type = decltype(std::declval<TOrigSeqs&>() | seqan3::views::translate_join);
};

// Case: nucleotide mode (not bisulfite)
template <typename TOrigSeqs, typename TRedAlph>
    requires seqan3::nucleotide_alphabet<TRedAlph>
struct TTransQrySeqsImpl<TOrigSeqs, seqan3::dna5, seqan3::dna5, TRedAlph>
{
    using type = decltype(std::declval<TOrigSeqs&>() | views::add_reverse_complement);
};

// Case: bisulfite mode
template <typename TOrigSeqs>
struct TTransQrySeqsImpl<TOrigSeqs, seqan3::dna5, seqan3::dna5, seqan3::semialphabet_any<6>>
{
    using type = decltype(std::declval<TOrigSeqs&>() | views::add_reverse_complement | views::duplicate);
};

// Case: no translation
template <typename TOrigSeqs, typename TRedAlph>
struct TTransQrySeqsImpl<TOrigSeqs, seqan3::aa27, seqan3::aa27, TRedAlph>
{
    using type = decltype(seqan3::views::type_reduce(std::declval<TOrigSeqs &>()));
};

// Reduced sequences
// Case: reduction of amino acid alphabet
template <typename TTransSeqs, typename TTransAlph, typename TRedAlph>
struct TRedSeqsImpl
{
    using type = decltype(std::declval<TTransSeqs&>() | seqan3::views::deep{seqan3::views::convert<TRedAlph>});
};

// Case: no reduction
template <typename TTransSeqs, typename TSameAlph>
struct TRedSeqsImpl<TTransSeqs, TSameAlph, TSameAlph> // no reduction
{
    using type = decltype(seqan3::views::type_reduce(std::declval<TTransSeqs &>()));
};

// Case: nucleotide reduction to dna4 (replcae Ns with random A,C,G,T)
template <typename TTransSeqs>
struct TRedSeqsImpl<TTransSeqs, seqan3::dna5, seqan3::dna4>
{
    using type = decltype(std::declval<TTransSeqs &>() | views::dna_n_to_random);
};

// Case: bisulfite mode
template <typename TTransSeqs>
struct TRedSeqsImpl<TTransSeqs, seqan3::dna5, seqan3::semialphabet_any<6>>
{
    using type = decltype(std::declval<TTransSeqs &>() | views::dna_n_to_random
                                                       | views::reduce_to_bisulfite);
};

// ==========================================================================
//  The index
// ==========================================================================

struct index_file_options
{
    uint64_t indexGeneration{0}; // bump this on incompatible changes

    DbIndexType indexType{};

    AlphabetEnum origAlph{};
    AlphabetEnum transAlph{};
    AlphabetEnum redAlph{};

    seqan3::genetic_code geneticCode{};

    //TODO reserve space here for more vars?

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("generation", indexGeneration),
                cereal::make_nvp("index type", indexType),
                cereal::make_nvp("orig alph", origAlph),
                cereal::make_nvp("trans alph", transAlph),
                cereal::make_nvp("red alph", redAlph),
                cereal::make_nvp("genetic code", geneticCode));
    }
};

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
