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

#include <cmath>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/details/traits.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <fmindex-collection/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>

#include <bio/alphabet/aminoacid/aa10li.hpp>
#include <bio/alphabet/aminoacid/aa10murphy.hpp>
#include <bio/alphabet/aminoacid/aa27.hpp>
#include <bio/alphabet/aminoacid/translation_genetic_code.hpp>
#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/convert.hpp>
#include <bio/ranges/views/deep.hpp>
#include <bio/ranges/views/translate_join.hpp>
#include <bio/ranges/views/type_reduce.hpp>

#include <bio/ranges/views/add_reverse_complement.hpp>
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

inline std::string _indexEnumToName(DbIndexType const t)
{
    switch (t)
    {
        case DbIndexType::FM_INDEX:
            return "fm_index";
        case DbIndexType::BI_FM_INDEX:
            return "bi_fm_index";
    }

    throw std::runtime_error("Error: unknown index type");
    return "";
}

inline DbIndexType _indexNameToEnum(std::string const t)
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

constexpr char const * _alphTypeToName(bio::alphabet::semialphabet_any<6> const & /**/)
{
    return "dna3bs";
}

constexpr char const * _alphTypeToName(bio::alphabet::dna4 const & /**/)
{
    return "dna4";
}

constexpr char const * _alphTypeToName(bio::alphabet::dna5 const & /**/)
{
    return "dna5";
}

constexpr char const * _alphTypeToName(bio::alphabet::aa27 const & /**/)
{
    return "aminoacid";
}

constexpr char const * _alphTypeToName(bio::alphabet::aa10murphy const & /**/)
{
    return "murphy10";
}

constexpr char const * _alphTypeToName(bio::alphabet::aa10li const & /**/)
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

inline std::string _alphabetEnumToName(AlphabetEnum const t)
{
    switch (t)
    {
        case AlphabetEnum::UNDEFINED:
            return "UNDEFINED";
        case AlphabetEnum::DNA3BS:
            return _alphTypeToName(bio::alphabet::semialphabet_any<6>{});
        case AlphabetEnum::DNA4:
            return _alphTypeToName(bio::alphabet::dna4{});
        case AlphabetEnum::DNA5:
            return _alphTypeToName(bio::alphabet::dna5{});
        case AlphabetEnum::AMINO_ACID:
            return _alphTypeToName(bio::alphabet::aa27{});
        case AlphabetEnum::MURPHY10:
            return _alphTypeToName(bio::alphabet::aa10murphy{});
        case AlphabetEnum::LI10:
            return _alphTypeToName(bio::alphabet::aa10li{});
    }

    throw std::runtime_error("Error: unknown alphabet type");
    return "";
}

inline AlphabetEnum _alphabetNameToEnum(std::string const t)
{
    if ((t == "UNDEFINED") || (t == "auto"))
        return AlphabetEnum::UNDEFINED;
    else if (t == _alphTypeToName(bio::alphabet::semialphabet_any<6>{}))
        return AlphabetEnum::DNA3BS;
    else if (t == _alphTypeToName(bio::alphabet::dna4{}))
        return AlphabetEnum::DNA4;
    else if (t == _alphTypeToName(bio::alphabet::dna5{}))
        return AlphabetEnum::DNA5;
    else if (t == _alphTypeToName(bio::alphabet::aa27{}))
        return AlphabetEnum::AMINO_ACID;
    else if (t == _alphTypeToName(bio::alphabet::aa10murphy{}))
        return AlphabetEnum::MURPHY10;
    else if (t == _alphTypeToName(bio::alphabet::aa10li{}))
        return AlphabetEnum::LI10;

    throw std::runtime_error("Error: unknown alphabet type");
    return AlphabetEnum::DNA4;
}

template <AlphabetEnum e>
struct _alphabetEnumToType_;

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA3BS>
{
    using type = bio::alphabet::semialphabet_any<6>;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA4>
{
    using type = bio::alphabet::dna4;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::DNA5>
{
    using type = bio::alphabet::dna5;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::AMINO_ACID>
{
    using type = bio::alphabet::aa27;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::MURPHY10>
{
    using type = bio::alphabet::aa10murphy;
};

template <>
struct _alphabetEnumToType_<AlphabetEnum::LI10>
{
    using type = bio::alphabet::aa10li;
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

template <size_t AlphabetSize>
using IndexSpec = fmindex_collection::occtable::interleavedEPR32V2::OccTable<AlphabetSize + 1>;

// ==========================================================================
//  Misc. aliases
// ==========================================================================

template <typename TString>
using TCDStringSet = bio::ranges::concatenated_sequences<TString>;

template <AlphabetEnum c_origSbjAlph, AlphabetEnum c_transAlph, AlphabetEnum c_redAlph>
inline constexpr auto sbjTransView = []()
{
    if constexpr (c_redAlph == AlphabetEnum::DNA3BS)
        return views::duplicate;
    else if constexpr (c_origSbjAlph != c_transAlph)
        return bio::views::translate_join;
    else
        return bio::views::type_reduce;
}();

template <AlphabetEnum c_origQryAlph, AlphabetEnum c_transAlph, AlphabetEnum c_redAlph>
constexpr auto qryTransView = []()
{
    if constexpr (c_redAlph == AlphabetEnum::DNA3BS)
        return bio::views::add_reverse_complement | views::duplicate;
    else if constexpr (bio::alphabet::nucleotide_alphabet<_alphabetEnumToType<c_redAlph>>)
        return bio::views::add_reverse_complement;
    else if constexpr (c_origQryAlph == c_transAlph)
        return bio::views::type_reduce;
    else
        return bio::views::translate_join;
}();

template <AlphabetEnum c_transAlph, AlphabetEnum c_redAlph>
constexpr auto redView = []()
{
    if constexpr (c_transAlph == c_redAlph)
        return bio::views::type_reduce;
    else if constexpr (c_transAlph == AlphabetEnum::AMINO_ACID)
        return bio::views::deep{bio::views::convert<_alphabetEnumToType<c_redAlph>>};
    else if constexpr (c_redAlph == AlphabetEnum::DNA3BS)
        return views::dna_n_to_random | views::reduce_to_bisulfite;
    else
        return views::dna_n_to_random;
}();

// ==========================================================================
//  overload serialisation
// ==========================================================================

namespace cereal
{

template <bio::alphabet::alphabet alph_t>
    requires std::is_trivially_copyable_v<alph_t>
void save(cereal::BinaryOutputArchive & archive, std::vector<alph_t> const & vec)
{
    archive(static_cast<uint64_t>(vec.size()));
    archive(cereal::binary_data(vec.data(), vec.size() * sizeof(alph_t)));
}

template <bio::alphabet::alphabet alph_t>
    requires std::is_trivially_copyable_v<alph_t>
void load(cereal::BinaryInputArchive & archive, std::vector<alph_t> & vec)
{
    uint64_t s = 0;
    archive(s);

    vec.resize(s);
    archive(cereal::binary_data(vec.data(), s * sizeof(alph_t)));
}

} // namespace cereal

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

    bio::alphabet::genetic_code geneticCode{};

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
template <DbIndexType  dbIndexType,
          AlphabetEnum origAlph,
          AlphabetEnum redAlph> // <- all members of index_file_options that influence types
struct index_file
{
    index_file_options options{};

    TCDStringSet<std::string>                                ids;
    TCDStringSet<std::vector<_alphabetEnumToType<origAlph>>> seqs;
    TCDStringSet<std::vector<uint32_t>>                      sTaxIds; //TODO double check int-width

    std::vector<uint32_t>     taxonParentIDs;
    std::vector<uint8_t>      taxonHeights;
    TCDStringSet<std::string> taxonNames;

    using TRedAlph   = _alphabetEnumToType<redAlph>;
    using TIndexSpec = IndexSpec<bio::alphabet::size<TRedAlph>>;
    using TIndex     = std::conditional_t<dbIndexType == DbIndexType::BI_FM_INDEX,
                                      fmindex_collection::BiFMIndex<TIndexSpec>,
                                      fmindex_collection::ReverseFMIndex<TIndexSpec>>;

    // Special c'tor that supports 'default' initialization to allow deserialize
    TIndex index{fmindex_collection::cereal_tag{}};

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options", options),
                cereal::make_nvp("ids", ids),
                cereal::make_nvp("seqs", seqs),
                cereal::make_nvp("sTaxIds", sTaxIds),
                cereal::make_nvp("taxonParentIDs", taxonParentIDs),
                cereal::make_nvp("taxonHeights", taxonHeights),
                cereal::make_nvp("taxonNames", taxonNames),
                cereal::make_nvp("index", index));
    }
};

/* pseudo-index that can be used to only load the options from disk */
struct fake_index_file
{
    index_file_options & options;

    template <typename TArchive>
    void serialize(TArchive & archive)
    {
        archive(cereal::make_nvp("options", options));
    }
};
