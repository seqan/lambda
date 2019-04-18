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
// lambda.cpp: Main File for the main application
// ==========================================================================

#include <initializer_list>
#include <iostream>

#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#define LAMBDA_INDEXER 1 // some things are different for the indexer binary

#include "shared_options.hpp"
#include "shared_definitions.hpp"
#include "shared_misc.hpp"

#include "mkindex_misc.hpp"
#include "mkindex_options.hpp"
// #include "mkindex_saca.hpp"
#include "mkindex_algo.hpp"

// using namespace seqan;

// ==========================================================================
// Forwards
// ==========================================================================

void argConv0(LambdaIndexerOptions & options);

template <DbIndexType c_indexType>
void argConv1(LambdaIndexerOptions & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origAlph>
void argConv2(LambdaIndexerOptions const & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origAlph,
          AlphabetEnum  c_transAlph>
void argConv3(LambdaIndexerOptions const & options);

template <DbIndexType           dbIndexType,
          AlphabetEnum          origAlph,
          AlphabetEnum          transAlph,
          AlphabetEnum          redAph>
void realMain(LambdaIndexerOptions     const & options);

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int mkindexMain(int const argc, char const ** argv)
{
    LambdaIndexerOptions options;
    parseCommandLine(options, argc, argv);

#ifdef NDEBUG
    try
    {
        argConv0(options);
    } catch (std::bad_alloc const & e)
    {
        std::cerr << "ERROR: Lambda ran out of memory :(\n"
                     "       You need to split your file into smaller segments.\n";
        return -1;
    } catch (std::exception const & e)
    {
        std::cerr << "\n\nERROR: The following unspecified exception was thrown:\n"
                  <<     "       \"" << e.what() << "\"\n"
                  <<     "       If the problem persists, report an issue at https://github.com/seqan/lambda/issues "
                  << "and include this output, as well as the output of `lambda2 --version`, thanks!\n";
        return -1;
    }
#else
    // In debug mode we don't catch the exceptions so that we get a backtrace from SeqAn's handler
    argConv0(options);
#endif
    return 0;
}

void argConv0(LambdaIndexerOptions & options)
{
    switch (options.indexFileOptions.dbIndexType)
    {
        case DbIndexType::FM_INDEX:         return argConv1<DbIndexType::FM_INDEX>(options);
        case DbIndexType::BI_FM_INDEX:      return argConv1<DbIndexType::BI_FM_INDEX>(options);
        default:                            throw 42;
    }
}

template <DbIndexType c_indexType>
void argConv1(LambdaIndexerOptions & options)
{
    if (options.indexFileOptions.origAlph = AlphabetEnum::UNDEFINED)
    {
        myPrint(options, 1, "Detecting database alphabet... ");
        options.indexFileOptions.origAlph = detectSeqFileAlphabet(options.dbFile);
        myPrint(options, 1, _alphabetEnumToName(options.subjOrigAlphabet), " detected.\n");
    }

    switch (options.indexFileOptions.origAlph)
    {
//         case AlphabetEnum::DNA4:            return argConv2<c_indexType, AlphabetEnum::DNA4>(options);
        case AlphabetEnum::DNA5:            return argConv2<c_indexType, AlphabetEnum::DNA5>(options);
        case AlphabetEnum::AMINO_ACID:      return argConv3<c_indexType, AlphabetEnum::AMINO_ACID, AlphabetEnum::AMINO_ACID>(options);
        default:                            throw 43;
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origAlph>
void argConv2(LambdaIndexerOptions const & options)
{
    switch (options.indexFileOptions.transAlph)
    {
//         case AlphabetEnum::DNA4:            return realMain<c_indexType, c_origAlph, AlphabetEnum::DNA4, AlphabetEnum::DNA4>(options);
        case AlphabetEnum::DNA5:            return realMain<c_indexType, c_origAlph, AlphabetEnum::DNA5, AlphabetEnum::DNA5>(options);
        case AlphabetEnum::AMINO_ACID:      return argConv3<c_indexType, c_origAlph, AlphabetEnum::AMINO_ACID>(options);
        default:                            throw 44;
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origAlph,
          AlphabetEnum  c_transAlph>
void argConv3(LambdaIndexerOptions const & options)
{
    switch (options.indexFileOptions.redAlph)
    {
        case AlphabetEnum::AMINO_ACID:      return realMain<c_indexType, c_origAlph, c_transAlph, AlphabetEnum::AMINO_ACID>(options);
        case AlphabetEnum::MURPHY10:        return realMain<c_indexType, c_origAlph, c_transAlph, AlphabetEnum::MURPHY10>(options);
        //TODO other reduced alphabets
        default:                            throw 45;
    }
}

template <DbIndexType   c_dbIndexType,
          AlphabetEnum  c_origAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAph>
void realMain(LambdaIndexerOptions     const & options)
{
    index_file<c_dbIndexType, c_origAlph, c_transAlph, c_redAlph> f;

    using TOrigSubjAlph = _alphabetEnumToType_<origAlph>;
    using TTransSubjAlph = _alphabetEnumToType_<transAlph>;

    using TOrigSet  = TCDStringSet<std::vector<TOrigSubjAlph>>;
    using TTransSet = TCDStringSet<std::vector<TTransSubjAlph>>;

    TTransSet translatedSeqs;

    {
        TOrigSet originalSeqs;
        std::unordered_map<std::string, uint64_t> accToIdRank;

        // ids get saved to disk again immediately and are not kept in memory
        loadSubjSeqsAndIds(originalSeqs, accToIdRank, options);

        // preserve lengths of untranslated sequences
        if constexpr (need_to_translate)
            _saveOriginalSeqLengths(originalSeqs, options);

        if (options.hasSTaxIds)
        {
            std::vector<bool> taxIdIsPresent;
            taxIdIsPresent.reserve(2'000'000);

            // read the mapping file and save relevant mappings to disk
            mapAndDumpTaxIDs(taxIdIsPresent, accToIdRank, std::ranges::size(originalSeqs), options);

            // read the mapping file and save relevant mappings to disk
            parseAndDumpTaxTree(taxIdIsPresent, options);
        }

        // translate or swap depending on program
        translateOrSwap(translatedSeqs, originalSeqs, options);
    }

    // dump translated and unreduced sequences (except where they are included in index)
    if ((options.alphReduction != 0) || (options.dbIndexType != DbIndexType::SUFFIX_ARRAY))
        dumpTranslatedSeqs(translatedSeqs, options);

    // see if final sequence set actually fits into index
//     checkIndexSize(translatedSeqs, options, seqan::BlastProgramSelector<p>());

    if (options.dbIndexType == DbIndexType::FM_INDEX)
        generateIndexAndDump<false>(translatedSeqs, options, TRedAlph{});
    else if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
        generateIndexAndDump<true>(translatedSeqs, options, TRedAlph{});

    // dump options
    for (auto && s : std::initializer_list<std::pair<std::string, std::string>>
         {
             { options.indexDir + "/option:db_index_type",   std::to_string(static_cast<uint32_t>(options.dbIndexType))},
             { options.indexDir + "/option:alph_original",   std::string(_alphTypeToName(TOrigSubjAlph())) },
             { options.indexDir + "/option:alph_translated", std::string(_alphTypeToName(TTransSubjAlph())) },
             { options.indexDir + "/option:alph_reduced",    std::string(_alphTypeToName(TRedAlph())) },
             { options.indexDir + "/option:genetic_code",    std::to_string(static_cast<uint16_t>(options.geneticCode)) },
//              { options.indexDir + "/option:subj_seq_len_bits", std::to_string(sizeof(SizeTypePos_<TRedAlph>) * 8)},
             { options.indexDir + "/option:generation",      std::to_string(indexGeneration) },
         })
    {
        std::ofstream f{std::get<0>(s).c_str(),  std::ios_base::out | std::ios_base::binary};
        f << std::get<1>(s);
        f.close();
    }
}
