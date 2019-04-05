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

template <bool nucleotide_mode, bool need_to_translate>
void argConv1(LambdaIndexerOptions const & options);

template <bool nucleotide_mode,
          bool need_to_translate,
          typename TRedAlph>
void realMain(LambdaIndexerOptions const & options);

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
    // set blastProgram
    if (!options.nucleotide_mode) // already implies != BLASTN
    {
        //TODO apparently ignores all manually set options
        myPrint(options, 1, "Detecting database alphabet... ");
        options.subjOrigAlphabet = detectSeqFileAlphabet(options.dbFile);
        myPrint(options, 1, _alphabetEnumToName(options.subjOrigAlphabet), " detected.\n");
        if (options.subjOrigAlphabet == AlphabetEnum::DNA5) // needs to be translated
            options.need_to_translate = true;
        else
            options.need_to_translate = true;
    }

    if (options.nucleotide_mode)
        realMain<true, false, seqan3::dna5>(options);
    else if (options.need_to_translate)
        argConv1<false, true>(options);
    else
        argConv1<false, false>(options);
}

/// Alphabet reduction (skipped in case == BLASTN)
template <bool nucleotide_mode, bool need_to_translate>
void argConv1(LambdaIndexerOptions           const & options)
{
    switch (options.reducedAlphabet)
    {
        case AlphabetEnum::AMINO_ACID:
            return realMain<nucleotide_mode, need_to_translate, seqan3::aa27>(options);
//TODO reactivate when we have it in SeqAn3
//         case AlphabetEnum::MURPHY10:
//             return realMain<nucleotide_mode, need_to_translate, seqan3::aa10murphy>(options);
#if 0
        case 10:
            return argConv2(options, ReducedAminoAcid<ClusterReduction<10>>());
        case 1:
            return argConv2(options, AminoAcid10());
        case 8:
            return argConv2(options, ReducedAminoAcid<ClusterReduction<8>>());
        case 12:
            return argConv2(options, ReducedAminoAcid<ClusterReduction<12>>());
#endif
        default:
            break;
    }
    throw std::invalid_argument("ERROR: Could not determine alphabet reduction.\n");
}

template <bool nucleotide_mode,
          bool need_to_translate,
          typename TRedAlph>
void realMain(LambdaIndexerOptions     const & options)
{
    using TOrigSubjAlph = std::conditional_t<nucleotide_mode,
                                             seqan3::dna5,
                                             std::conditional_t<need_to_translate, seqan3::dna5, seqan3::aa27>>;
    using TTransSubjAlph = std::conditional_t<nucleotide_mode, seqan3::dna5, seqan3::aa27>;

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
