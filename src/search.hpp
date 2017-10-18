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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/reduced_aminoacid.h>
#include <seqan/misc/terminal.h>

#include "shared_definitions.hpp"
#include "shared_options.hpp"
#include "shared_misc.hpp"

#include "search_output.hpp"
#include "search_options.hpp"
#include "search_datastructures.hpp"
#include "search_misc.hpp"
#include "search_algo.hpp"

// forwards

inline int
argConv0(LambdaOptions & options);
//-
template <typename TOutFormat,
          BlastTabularSpec h>
inline int
argConv1(LambdaOptions                       & options,
         TOutFormat                    const & /**/,
         BlastTabularSpecSelector<h>   const &);
//-
template <typename TOutFormat,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv2(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &);
//-
template <typename TOutFormat,
          typename TRedAlph,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv3(LambdaOptions                        & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &);
//-
template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreExtension,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv4(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreExtension                const & /**/);
//-
template <typename TIndexSpec,
          typename TRedAlph,
          typename TScoreExtension,
          typename TOutFormat,
          BlastProgram p,
          BlastTabularSpec h>
inline int
realMain(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreExtension                const & /**/);

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int searchMain(int const argc, char const ** argv)
{
    LambdaOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (std::string(CMAKE_BUILD_TYPE) != "Release")
        std::cerr << "WARNING: This binary is not built in release mode and will be much slower than it should be!\n";

#ifdef NDEBUG
    int ret = 0;
    try
    {
        ret = argConv0(options);
    } catch (std::bad_alloc const & e)
    {
        std::cerr << "ERROR: Lambda ran out of memory :(\n"
                     "       You need to split your file into smaller segments or search against a smaller database.\n";
        return -1;
    } catch (std::exception const & e)
    {
        std::cerr << "\n\n**An exception was thrown**\n" << e.what() << '\n';
        return -1;
    }
    return ret;
#else
    // In debug mode we don't catch the exceptions so that we get a backtrace from SeqAn's handler
    return argConv0(options);
#endif
}

// CONVERT Run-time options to compile-time Format-Type
inline int
argConv0(LambdaOptions & options)
{
    myPrint(options, 1, "LAMBDA - the Local Aligner for Massive Biological DatA"
                        "\n======================================================"
                        "\nVersion ", SEQAN_APP_VERSION, "\n\n");

    // Index
    myPrint(options, 1, "Reading index properties... ");
    readIndexOptions(options);
    myPrint(options, 1, "done.\n");

#ifndef LAMBDA_LEGACY_PATHS
    if (options.dbIndexType == DbIndexType::SUFFIX_ARRAY)
        throw std::runtime_error("ERROR: The index is of type suffix array, but support was removed.\n"
                                 "       Either rebuild lambda2 with '-DLAMBDA_LEGACY_PATHS=1' or re-create the index.\n");
#endif // LAMBDA_LEGACY_PATHS

    myPrint(options, 2, "Index properties\n"
                        "  type:                ", _indexEnumToName(options.dbIndexType), "\n",
                        "  original   alphabet: ", _alphabetEnumToName(options.subjOrigAlphabet), "\n");
    if (_alphabetEnumToName(options.subjOrigAlphabet) == _alphabetEnumToName(options.transAlphabet))
    {
        myPrint(options, 2, "  translated alphabet: not translated\n");
    }
    else
    {
        myPrint(options, 2, "  translated alphabet: ", _alphabetEnumToName(options.transAlphabet), "\n");
        myPrint(options, 2, "    translation code:  ", options.geneticCode, "\n");
    }

    if (_alphabetEnumToName(options.transAlphabet) == _alphabetEnumToName(options.reducedAlphabet))
    {
        myPrint(options, 2, "  reduced    alphabet:  not reduced\n");
    }
    else
    {
        myPrint(options, 2, "  reduced    alphabet: ", _alphabetEnumToName(options.reducedAlphabet), "\n\n");
    }

    if ((options.blastProgram == BlastProgram::BLASTN) &&
        (options.reducedAlphabet != AlphabetEnum::DNA4))
    {
        throw std::runtime_error("ERROR: You are attempting a nucleotide search on a protein index.\n");
    }

    // query file
    if (options.qryOrigAlphabet == AlphabetEnum::DNA4) // means "auto", as dna4 not valid as argument to --query-alphabet
    {
        myPrint(options, 1, "Detecting query alphabet... ");
        options.qryOrigAlphabet = detectSeqFileAlphabet(options.queryFile);
        myPrint(options, 1, _alphabetEnumToName(options.qryOrigAlphabet), " detected.\n");
    }

    // set blastProgram
    if (options.blastProgram == BlastProgram::UNKNOWN)
    {
        if ((options.transAlphabet == AlphabetEnum::DNA5) && (options.qryOrigAlphabet == AlphabetEnum::AMINO_ACID))
        {
            throw std::runtime_error("ERROR: Query file is protein, but index is nucleotide. "
                                    "Re-create the index with 'lambda mkindexp'.\n");
        }
        else if ((options.transAlphabet == AlphabetEnum::DNA5) && (options.qryOrigAlphabet == AlphabetEnum::DNA5))
        {
            options.blastProgram = BlastProgram::BLASTN;
        }
        else if ((options.qryOrigAlphabet == AlphabetEnum::DNA5)) // query will be translated
        {
            if (options.subjOrigAlphabet == options.transAlphabet)
                options.blastProgram = BlastProgram::BLASTX;
            else
                options.blastProgram = BlastProgram::TBLASTX;
        }
        else // query is aminoacid already
        {
            if (options.subjOrigAlphabet == options.transAlphabet)
                options.blastProgram = BlastProgram::BLASTP;
            else
                options.blastProgram = BlastProgram::TBLASTN;
        }
    }

    // sizes
    int ret = checkRAM(options);
    if (ret)
        return ret;

    // output
    CharString output = options.output;
    if (endsWith(output, ".gz"))
        output = prefix(output, length(output) - 3);
    else if (endsWith(output, ".bz2"))
        output = prefix(output, length(output) - 4);

    if (endsWith(output, ".m0"))
        return argConv1(options, BlastReport(), BlastTabularSpecSelector<BlastTabularSpec::NO_COMMENTS>());
    else if (endsWith(output, ".m8"))
        return argConv1(options, BlastTabular(), BlastTabularSpecSelector<BlastTabularSpec::NO_COMMENTS>());
    else if (endsWith(output, ".m9"))
        return argConv1(options, BlastTabular(), BlastTabularSpecSelector<BlastTabularSpec::COMMENTS>());
    else if (endsWith(output, ".sam") || endsWith(output, ".bam")) // handled elsewhere
        return argConv1(options, BlastTabular(), BlastTabularSpecSelector<BlastTabularSpec::COMMENTS>());

    std::cerr << "ERROR: Cannot handle output extension.\n";
    return -1;
}

template <typename TOutFormat,
          BlastTabularSpec h>
inline int
argConv1(LambdaOptions                       & options,
         TOutFormat                    const & /**/,
         BlastTabularSpecSelector<h>   const &)
{
    switch(options.blastProgram)
    {
#ifndef FASTBUILD
        case BlastProgram::BLASTN:
            return argConv3(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTN>(),
                            Dna5());
#endif
        case BlastProgram::BLASTP:
            return argConv2(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTP>());
        case BlastProgram::BLASTX:
            return argConv2(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTX>());
#ifndef FASTBUILD
        case BlastProgram::TBLASTN:
            return argConv2(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::TBLASTN>());
        case BlastProgram::TBLASTX:
            return argConv2(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::TBLASTX>());
#endif
        default:
            break;
    }
    std::cerr << "ERROR: Cannot handle program mode (perhaps you are building in FASTMODE?).\n";
    return -1;
}


/// Alphabet reduction
template <typename TOutFormat,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv2(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &)
{
    using Th = BlastTabularSpecSelector<h>;
    using Tp = BlastProgramSelector<p>;

    switch (options.reducedAlphabet)
    {
        case AlphabetEnum::AMINO_ACID:
            return argConv3(options, TOutFormat(), Th(), Tp(), AminoAcid());
        case AlphabetEnum::MURPHY10:
            return argConv3(options, TOutFormat(), Th(), Tp(), ReducedAminoAcid<Murphy10>());
#if 0
        case 10:
            return argConv2(options, TOutFormat(), ReducedAminoAcid<ClusterReduction<10>>());
        case 1:
            return argConv2(options, TOutFormat(), AminoAcid10());
        case 8:
            return argConv2(options, TOutFormat(), ReducedAminoAcid<ClusterReduction<8>>());
        case 12:
            return argConv2(options, TOutFormat(), ReducedAminoAcid<ClusterReduction<12>>());
#endif
        default:
            break;
    }
    std::cerr << "ERROR: Cannot handle the specified alphabet reduction.\n";
    return -1;
}

// extension model
template <typename TOutFormat,
          typename TRedAlph,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv3(LambdaOptions                        & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &)
{

    if (options.gapOpen == 0)
#ifndef LAMBDA_LINGAPS_OPT
        std::cerr << "ATTENTION: You have set the additional gap open cost to 0. If you run LAMBDA "
                     "in this configuration regularly, you might want to rebuild it with "
                     "LAMBDA_LINGAPS_OPT=1 to profit from additional optimizations.\n";
#else
        return argConv4(options,
                        TOutFormat(),
                        BlastTabularSpecSelector<h>(),
                        BlastProgramSelector<p>(),
                        TRedAlph(),
                        LinearGaps());
    else
#endif
        return argConv4(options,
                        TOutFormat(),
                        BlastTabularSpecSelector<h>(),
                        BlastProgramSelector<p>(),
                        TRedAlph(),
                        AffineGaps());
}

template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreExtension,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv4(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreExtension                const & /**/)
{
#ifdef LAMBDA_LEGACY_PATHS
    if (options.dbIndexType == DbIndexType::SUFFIX_ARRAY)
        return realMain<IndexSa<>>(options,
                                   TOutFormat(),
                                   BlastTabularSpecSelector<h>(),
                                   BlastProgramSelector<p>(),
                                   TRedAlph(),
                                   TScoreExtension());
    else
#endif // LAMBDA_LEGACY_PATHS
    if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
        return realMain<BidirectionalIndex<TFMIndexInBi<>>>(options,
                                                            TOutFormat(),
                                                            BlastTabularSpecSelector<h>(),
                                                            BlastProgramSelector<p>(),
                                                            TRedAlph(),
                                                            TScoreExtension());
    else
        return realMain<TFMIndex<>>(options,
                                   TOutFormat(),
                                   BlastTabularSpecSelector<h>(),
                                   BlastProgramSelector<p>(),
                                   TRedAlph(),
                                   TScoreExtension());
}

/// REAL MAIN
#ifdef _OPENMP
#define TID omp_get_thread_num()
#else
#define TID 0
#endif
template <typename TIndexSpec,
          typename TRedAlph,
          typename TScoreExtension,
          typename TOutFormat,
          BlastProgram p,
          BlastTabularSpec h>
inline int
realMain(LambdaOptions                        & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreExtension                const & /**/)
{
    using TGlobalHolder = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;
    using TLocalHolder = LocalDataHolder<TGlobalHolder, TScoreExtension>;

    if (options.verbosity >= 2)
        printOptions<TLocalHolder>(options);

    TGlobalHolder globalHolder;
//     context(globalHolder.outfile).scoringScheme._internalScheme = matr;

    int ret = prepareScoring(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSubjects(globalHolder, options);
    if (ret)
        return ret;

    ret = loadDbIndexFromDisk(globalHolder, options);
    if (ret)
        return ret;

    ret = loadTaxonomy(globalHolder, options);
    if (ret)
        return ret;

    ret = loadQuery(globalHolder, options);
    if (ret)
        return ret;

//     std::cout << "1st Query:\n"
//               << front(globalHolder.qrySeqs) << "\n"
//               << front(globalHolder.redQrySeqs) << "\n";
//
//     std::cout << "last Query:\n"
//               << back(globalHolder.qrySeqs) << "\n"
//               << back(globalHolder.redQrySeqs) << "\n";
//
//     std::cout << "1st Subject:\n"
//               << front(globalHolder.subjSeqs) << "\n"
//               << front(globalHolder.redSubjSeqs) << "\n";
//
//     std::cout << "last Subject:\n"
//               << back(globalHolder.subjSeqs) << "\n"
//               << back(globalHolder.redSubjSeqs) << "\n";

    myWriteHeader(globalHolder, options);

    if (options.doubleIndexing)
    {
        myPrint(options, 1,
                "Searching ",
                options.queryPart,
                " blocks of query with ",
                options.threads,
                " threads...\n");
        if ((options.isTerm) && (options.verbosity >= 1))
        {
            for (unsigned char i=0; i< options.threads+3; ++i)
                std::cout << std::endl;
            std::cout << "\033[" << options.threads+2 << "A";
        }
    } else
    {
        myPrint(options, 1, "Searching and extending hits on-line...progress:\n"
                "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n|");
    }
    double start = sysTime();

    uint64_t lastPercent = 0;

    SEQAN_OMP_PRAGMA(parallel)
    {
        TLocalHolder localHolder(options, globalHolder);

        SEQAN_OMP_PRAGMA(for schedule(dynamic))
        for (uint64_t t = 0; t < localHolder.nBlocks; ++t)
        {
            int res = 0;

            localHolder.init(t);

            // seed
        #ifdef LAMBDA_MICRO_STATS
            double buf = sysTime();
        #endif
            if (options.doubleIndexing)
            {
                res = generateSeeds(localHolder);
                if (res)
                    continue;

                res = generateTrieOverSeeds(localHolder);
                if (res)
                    continue;
            }
        #ifdef LAMBDA_MICRO_STATS
            localHolder.stats.timeGenSeeds += sysTime() - buf;

            // search
            buf = sysTime();
        #endif
            search(localHolder); //TODO seed refining if iterateMatches gives 0 results
        #ifdef LAMBDA_MICRO_STATS
            localHolder.stats.timeSearch += sysTime() - buf;
        #endif

//             // TODO DEBUG
//             for (auto const & m : localHolder.matches)
//                 _printMatch(m);

            // sort
            if (options.filterPutativeAbundant || options.filterPutativeDuplicates || options.mergePutativeSiblings)
            {
            #ifdef LAMBDA_MICRO_STATS
                buf = sysTime();
            #endif

                sortMatches(localHolder);

            #ifdef LAMBDA_MICRO_STATS
                localHolder.stats.timeSort += sysTime() - buf;
            #endif
            }

            // extend
            if (length(localHolder.matches) > 0)
                res = iterateMatches(localHolder);

            if (res)
                continue;

            if ((!options.doubleIndexing) && (TID == 0) &&
                (options.verbosity >= 1))
            {
                unsigned curPercent = ((t * 50) / localHolder.nBlocks) * 2; // round to even
                printProgressBar(lastPercent, curPercent);
            }

        } // implicit thread sync here

        if ((!options.doubleIndexing) && (TID == 0) && (options.verbosity >= 1))
            printProgressBar(lastPercent, 100);

        SEQAN_OMP_PRAGMA(critical(statsAdd))
        {
            globalHolder.stats += localHolder.stats;
        }
    }

    if (ret)
        return ret;

    myWriteFooter(globalHolder, options);

    if (!options.doubleIndexing)
    {
        myPrint(options, 2, "Runtime total: ", sysTime() - start, "s.\n\n");
    }

    printStats(globalHolder.stats, options);

    return 0;
}
