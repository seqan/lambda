// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2017, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2017, Knut Reinert and Freie Universität Berlin
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

#ifndef LAMBDA_SHARED_OPTIONS_H_
#define LAMBDA_SHARED_OPTIONS_H_

#include <cstdio>
#include <unistd.h>
#include <bitset>

#include <seqan/basic.h>
#include <seqan/translation.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/blast.h>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Enum DbIndexType
// --------------------------------------------------------------------------

enum class DbIndexType : uint8_t
{
    SUFFIX_ARRAY,
    FM_INDEX,
    BI_FM_INDEX
};

inline std::string
_indexEnumToName(DbIndexType const t)
{
    switch (t)
    {
        case DbIndexType::SUFFIX_ARRAY:  return "suffix_array";
        case DbIndexType::FM_INDEX:      return "fm_index";
        case DbIndexType::BI_FM_INDEX:   return "bi_fm_index";
    }

    throw std::runtime_error("Error: unknown index type");
    return "";
}

inline DbIndexType
_indexNameToEnum(std::string const t)
{
    if (t == "suffix_array")
        return DbIndexType::SUFFIX_ARRAY;
    else if (t == "bi_fm_index")
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
_alphTypeToName(Dna const & /**/)
{
    return "dna4";
}

constexpr const char *
_alphTypeToName(Dna5 const & /**/)
{
    return "dna5";
}

constexpr const char *
_alphTypeToName(AminoAcid const & /**/)
{
    return "aminoacid";
}

constexpr const char *
_alphTypeToName(ReducedAminoAcid<Murphy10> const & /**/)
{
    return "murphy10";
}

enum class AlphabetEnum : uint8_t
{
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
        case AlphabetEnum::DNA4:        return _alphTypeToName(Dna{});
        case AlphabetEnum::DNA5:        return _alphTypeToName(Dna5{});
        case AlphabetEnum::AMINO_ACID:  return _alphTypeToName(AminoAcid{});
        case AlphabetEnum::MURPHY10:    return _alphTypeToName(ReducedAminoAcid<Murphy10>{});
    }

    throw std::runtime_error("Error: unknown alphabet type");
    return "";
}

inline AlphabetEnum
_alphabetNameToEnum(std::string const t)
{
    if (t == _alphTypeToName(Dna{}))
        return AlphabetEnum::DNA4;
    else if (t == _alphTypeToName(Dna5{}))
        return AlphabetEnum::DNA5;
    else if (t == _alphTypeToName(AminoAcid{}))
        return AlphabetEnum::AMINO_ACID;
    else if (t == _alphTypeToName(ReducedAminoAcid<Murphy10>{}))
        return AlphabetEnum::MURPHY10;

    throw std::runtime_error("Error: unknown alphabet type");
    return AlphabetEnum::DNA4;
}

// --------------------------------------------------------------------------
// Class SharedOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity = 1;

    std::string commandLine;

    std::string indexDir;

    DbIndexType dbIndexType;

    AlphabetEnum subjOrigAlphabet;
    AlphabetEnum transAlphabet;
    AlphabetEnum reducedAlphabet;

    GeneticCodeSpec geneticCode = static_cast<GeneticCodeSpec>(0);//CANONICAL;

    BlastProgram blastProgram   = BlastProgram::UNKNOWN;

    bool        isTerm          = true;
    unsigned    terminalCols    = 80;

    unsigned    threads         = 1;
    bool        hasSTaxIds      = false;

    SharedOptions()
    {
        isTerm = isTerminal();
        if (isTerm)
        {
            unsigned _rows;
            getTerminalSize(terminalCols, _rows);
        }
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function sharedSetup()
// --------------------------------------------------------------------------

void
sharedSetup(ArgumentParser & parser)
{
    // Set short description, version, and date.
    std::string versionString = SEQAN_APP_VERSION;
    setVersion(parser, versionString);
    setDate(parser, __DATE__);
    setShortCopyright(parser, "2013-2017 Hannes Hauswedell, released under the GNU AGPL v3 (or later); "
                              "2016-2017 Knut Reinert and Freie Universität Berlin, released under the 3-clause-BSDL");

    setCitation(parser, "Hauswedell et al (2014); doi: 10.1093/bioinformatics/btu439");

    setLongCopyright(parser,
        " Copyright (c) 2013-2017, Hannes Hauswedell\n"
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
        " Copyright (c) 2016-2017 Knut Reinert and Freie Universität Berlin\n"
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
        " DAMAGE.\n");

    addDescription(parser, "Lambda is a local aligner optimized for many query "
        "sequences and searches in protein space. It is compatible to BLAST, but "
        "much faster than BLAST and many other comparable tools.");

    addDescription(parser, "Detailed information is available in the wiki: "
        "<https://github.com/seqan/lambda/wiki>");
}

ArgumentParser::ParseResult
parseCommandLineShared(SharedOptions & options, ArgumentParser & parser)
{
    int buf = 0;

#ifdef _OPENMP
    getOptionValue(options.threads, parser, "threads");
    omp_set_num_threads(options.threads);
#else
    options.threads = 1;
#endif

    getOptionValue(buf, parser, "verbosity");
    switch(buf)
    {
        case 0: options.verbosity = 0; break;
        case 2: options.verbosity = 2; break;
        default: options.verbosity = 1; break;
    }

    return ArgumentParser::PARSE_OK;
}

#endif // header guard
