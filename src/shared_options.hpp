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
// shared_options.h: contains the options and argument parser
// ==========================================================================

#ifndef LAMBDA_SHARED_OPTIONS_H_
#define LAMBDA_SHARED_OPTIONS_H_

#include <bitset>
#include <cstdio>
#include <thread>
#include <unistd.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/std/filesystem>

// #include <seqan/basic.h>
// #include <seqan/modifier.h>
// #include <seqan/arg_parse.h>
// #include <seqan/index.h>
// #include <seqan/blast.h>
// #include <seqan/misc/terminal.h>
// #include <seqan/translation.h>
// #include <seqan/reduced_aminoacid.h>

// using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Enum DbIndexType
// --------------------------------------------------------------------------

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


// --------------------------------------------------------------------------
// Enum AlphabetEnum
// --------------------------------------------------------------------------

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

// constexpr const char *
// _alphTypeToName(ReducedAminoAcid<Murphy10> const & /**/)
// {
//     return "murphy10";
// }

enum class AlphabetEnum : uint8_t
{
    UNDEFINED,
    DNA4,
    DNA5,
    AMINO_ACID,
    MURPHY10,
};

inline std::string
_alphabetEnumToName(AlphabetEnum const t)
{
    switch (t)
    {
        case AlphabetEnum::UNDEFINED:   return "UNDEFINED";
        case AlphabetEnum::DNA4:        return _alphTypeToName(seqan3::dna4{});
        case AlphabetEnum::DNA5:        return _alphTypeToName(seqan3::dna5{});
        case AlphabetEnum::AMINO_ACID:  return _alphTypeToName(seqan3::aa27{});
        case AlphabetEnum::MURPHY10:    return _alphTypeToName(seqan3::aa10murphy{});
    }

    throw std::runtime_error("Error: unknown alphabet type");
    return "";
}

inline AlphabetEnum
_alphabetNameToEnum(std::string const t)
{
    if ((t == "UNDEFINED") || (t == "auto"))
        return AlphabetEnum::UNDEFINED;
    else if (t == _alphTypeToName(seqan3::dna4{}))
        return AlphabetEnum::DNA4;
    else if (t == _alphTypeToName(seqan3::dna5{}))
        return AlphabetEnum::DNA5;
    else if (t == _alphTypeToName(seqan3::aa27{}))
        return AlphabetEnum::AMINO_ACID;
    else if (t == _alphTypeToName(seqan3::aa10murphy{}))
        return AlphabetEnum::MURPHY10;

    throw std::runtime_error("Error: unknown alphabet type");
    return AlphabetEnum::DNA4;
}

template <AlphabetEnum e>
struct _alphabetEnumToType_;

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

template <AlphabetEnum e>
using _alphabetEnumToType = typename _alphabetEnumToType_<e>::type;

// inline uint64_t
// _alphabetEnumToSize(AlphabetEnum const t)
// {
//     switch (t)
//     {
//         case AlphabetEnum::DNA4:        return sizeof(SizeTypePos_<Dna>);
//         case AlphabetEnum::DNA5:        return sizeof(SizeTypePos_<Dna5>);
//         case AlphabetEnum::AMINO_ACID:  return sizeof(SizeTypePos_<AminoAcid>);
//         case AlphabetEnum::MURPHY10:    return sizeof(SizeTypePos_<ReducedAminoAcid<Murphy10>>);
//     }
//
//     throw std::runtime_error("Error: unknown alphabet type");
//     return 0;
// }

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

// --------------------------------------------------------------------------
// Class SharedOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose
    int32_t     verbosity = 1;

    std::string commandLine;

    std::filesystem::path indexFilePath;

    index_file_options indexFileOptions{};

//     seqan::BlastProgram blastProgram   = seqan::BlastProgram::UNKNOWN;
    bool        nucleotide_mode = false;
    bool        need_to_translate = false;

    bool        isTerm          = true;
    unsigned    terminalCols    = 80;

    uint64_t    threads         = std::thread::hardware_concurrency();

    bool        hasSTaxIds      = false;

    SharedOptions()
    {
//         isTerm = seqan::isTerminal();
//         if (isTerm)
//         {
//             unsigned _rows;
//             seqan::getTerminalSize(terminalCols, _rows);
//         }
    }
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function sharedSetup()
// --------------------------------------------------------------------------

void sharedSetup(seqan3::argument_parser & parser)
{
    // Set short description, version, and date
    parser.info.version = SEQAN_APP_VERSION;
    parser.info.date = __DATE__;
    parser.info.citation = "Hauswedell et al (2014); doi: 10.1093/bioinformatics/btu439";
    parser.info.short_copyright = "2013-2019 Hannes Hauswedell, released under the GNU AGPL v3 (or later); "
                                  "2016-2019 Knut Reinert and Freie Universität Berlin, released under the 3-clause-BSDL";
    parser.info.long_copyright = " Copyright (c) 2013-2019, Hannes Hauswedell\n"
                                 " All rights reserved.\n"
                                 "\n"
                                 " This program is free software: you can redistribute it and/or modify\n"
                                 " it under the terms of the GNU Affero General Public License as\n"
                                 " published by the Free Software Foundation, either version 3 of the\n"
                                 " License, or (at your option) any later version.\n"
                                 "\n"
                                 " Lambda is distributed in the hope that it will be useful,\n"
                                 " but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
                                 " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
                                 " GNU General Public License for more details.\n"
                                 "\n"
                                 " You should have received a copy of the GNU Affero General Public License\n"
                                 " along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
                                 "\n"
                                 " Copyright (c) 2016-2019 Knut Reinert and Freie Universität Berlin\n"
                                 " All rights reserved.\n"
                                 "\n"
                                 " Redistribution and use in source and binary forms, with or without\n"
                                 " modification, are permitted provided that the following conditions are met:\n"
                                 "\n"
                                 " * Redistributions of source code must retain the above copyright\n"
                                 "   notice, this list of conditions and the following disclaimer.\n"
                                 " * Redistributions in binary form must reproduce the above copyright\n"
                                 "   notice, this list of conditions and the following disclaimer in the\n"
                                 "   documentation and/or other materials provided with the distribution.\n"
                                 " * Neither the name of Knut Reinert or the FU Berlin nor the names of\n"
                                 "   its contributors may be used to endorse or promote products derived\n"
                                 "   from this software without specific prior written permission.\n"
                                 "\n"
                                 " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n"
                                 " AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n"
                                 " IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n"
                                 " ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE\n"
                                 " FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\n"
                                 " DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
                                 " SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
                                 " CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
                                 " LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY\n"
                                 " OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n"
                                 " DAMAGE.\n";
    parser.info.description.push_back("Lambda is a local aligner optimized for many query "
                                      "sequences and searches in protein space. It is compatible to BLAST, but "
                                      "much faster than BLAST and many other comparable tools.");


    parser.info.description.push_back("Detailed information is available in the wiki: "
                                      "<https://github.com/seqan/lambda/wiki>");
}

#endif // header guard
