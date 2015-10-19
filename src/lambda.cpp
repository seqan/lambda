// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2015, Hannes Hauswedell, FU Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
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

#include "options.hpp"
#include "match.hpp"
#include "lambda.hpp"
#include "misc.hpp"

using namespace seqan;


// inline BlastFormatFile
// _fileType(LambdaOptions const & options)
// {
//     if (endsWith(options.output, ".m0"))
//         return BlastFormatFile::PAIRWISE;
//     else if (endsWith(options.output, ".m8"))
//         return BlastFormatFile::TABULAR;
//     else if (endsWith(options.output, ".m9"))
//         return BlastFormatFile::TABULAR_WITH_COMMENTS;
//     else
//         return BlastFormatFile::UNKNOWN;
// }

// forwards

inline int
argConv0(LambdaOptions const & options);
//-
template <typename TOutFormat,
          BlastTabularSpec h>
inline int
argConv05(LambdaOptions                 const & options,
          TOutFormat                    const & /**/,
          BlastTabularSpecSelector<h>   const &);
//-
template <typename TOutFormat,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv1(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &);
//-
template <typename TOutFormat,
          typename TRedAlph,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv2(LambdaOptions                  const & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &);
//-
template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreScheme,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv3(LambdaOptions                  const & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &,
         TScoreScheme                   const &);
//-
template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv4(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreScheme                   const & /**/,
         TScoreExtension                const & /**/);
//-
template <typename TIndexSpec,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension,
          typename TOutFormat,
          BlastProgram p,
          BlastTabularSpec h>
inline int
realMain(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreScheme                   const & /**/,
         TScoreExtension                const & /**/);

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    LambdaOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

//     std::cout <<   "Match  sizeof : " << sizeof(Match)
//               << "\n       alignof: " << alignof(Match)
//               << "\n    is_trivial: " << std::is_trivial<Match>::value
// //               << "\ntrivially_copy: " << std::is_trivially_copyable<Match>::value
//               << "\n";

    if (std::string(CMAKE_BUILD_TYPE) != "Release")
        std::cerr << "WARNING: This binary is not built in release mode and will be much slower than it should be!\n";
    return argConv0(options);
}

// CONVERT Run-time options to compile-time Format-Type
inline int
argConv0(LambdaOptions const & options)
{
    CharString output = options.output;
    if (endsWith(output, ".gz"))
        output = prefix(output, length(output) - 3);
    else if (endsWith(output, ".bz2"))
        output = prefix(output, length(output) - 4);

    if (endsWith(output, ".m0"))
        return argConv05(options, BlastReport(), BlastTabularSpecSelector<BlastTabularSpec::NO_COMMENTS>());
    else if (endsWith(output, ".m8"))
        return argConv05(options, BlastTabular(), BlastTabularSpecSelector<BlastTabularSpec::NO_COMMENTS>());
    else if (endsWith(output, ".m9"))
        return argConv05(options, BlastTabular(), BlastTabularSpecSelector<BlastTabularSpec::COMMENTS>());
    return -1;
}

template <typename TOutFormat,
          BlastTabularSpec h>
inline int
argConv05(LambdaOptions                 const & options,
          TOutFormat                    const & /**/,
          BlastTabularSpecSelector<h>   const &)
{
    switch(options.blastProgram)
    {
#ifndef FASTBUILD
        case BlastProgram::BLASTN:
            return argConv1(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTN>());
#endif
        case BlastProgram::BLASTP:
            return argConv1(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTP>());
        case BlastProgram::BLASTX:
            return argConv1(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::BLASTX>());
#ifndef FASTBUILD
        case BlastProgram::TBLASTN:
            return argConv1(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::TBLASTN>());
        case BlastProgram::TBLASTX:
            return argConv1(options,
                            TOutFormat(),
                            BlastTabularSpecSelector<h>(),
                            BlastProgramSelector<BlastProgram::TBLASTX>());
#endif
        default:
            break;
    }
    return -1;
}


/// Alphabet reduction
template <typename TOutFormat,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv1(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &)
{
    using TUnred = typename std::conditional<p == BlastProgram::BLASTN, Dna5, AminoAcid>::type;
    using Th = BlastTabularSpecSelector<h>;
    using Tp = BlastProgramSelector<p>;
    switch (options.alphReduction)
    {
        case 0:
            return argConv2(options, TOutFormat(), Th(), Tp(), TUnred());
        case 2:
            return argConv2(options, TOutFormat(), Th(), Tp(), ReducedAminoAcid<Murphy10>());
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
            return -1;
    }
    return -1;
}

/// scoring scheme

template <typename TOutFormat,
          typename TRedAlph,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv2(LambdaOptions                  const & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &)
{
    using Th = BlastTabularSpecSelector<h>;
    using Tp = BlastProgramSelector<p>;
    switch (options.scoringMethod)
    {
#ifndef FASTBUILD
        case 0:
            return argConv3(options, TOutFormat(), Th(), Tp(), TRedAlph(), Score<int, Simple>());
        case 45:
            return argConv3(options, TOutFormat(), Th(), Tp(), TRedAlph(), Blosum45());
        case 80:
            return argConv3(options, TOutFormat(), Th(), Tp(), TRedAlph(), Blosum80());
#endif
        case 62:
            return argConv3(options, TOutFormat(), Th(), Tp(), TRedAlph(), Blosum62());
        default:
            std::cerr << "Unsupported Scoring Scheme selected.\n";
            break;
    }

    return -1;
}

template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreScheme,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv3(LambdaOptions                  const & options,
         TOutFormat                     const &,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const &,
         TScoreScheme                   const &)
{
#ifndef FASTBUILD
    if (options.gapOpen == 0)
        return argConv4(options,
                        TOutFormat(),
                        BlastTabularSpecSelector<h>(),
                        BlastProgramSelector<p>(),
                        TRedAlph(),
                        TScoreScheme(),
                        LinearGaps());
    else
#endif
        return argConv4(options,
                        TOutFormat(),
                        BlastTabularSpecSelector<h>(),
                        BlastProgramSelector<p>(),
                        TRedAlph(),
                        TScoreScheme(),
                        AffineGaps());
}

template <typename TOutFormat,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension,
          BlastTabularSpec h,
          BlastProgram p>
inline int
argConv4(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreScheme                   const & /**/,
         TScoreExtension                const & /**/)
{
    int indexType = options.dbIndexType;
//     if (indexType == -1) // autodetect
//     {
//         //TODO FIX THIS WITH NEW EXTENSIONS
//         CharString file = options.dbFile;
//         append(file, ".sa");
//         struct stat buffer;
//         if (stat(toCString(file), &buffer) == 0)
//         {
//             indexType = 0;
//         } else
//         {
//             file = options.dbFile;
//             append(file, ".sa.val"); // FM Index
//             struct stat buffer;
//             if (stat(toCString(file), &buffer) == 0)
//             {
//                 indexType = 1;
//             } else
//             {
//                 std::cerr << "No Index file could be found, please make sure paths "
//                         << "are correct and the files are readable.\n" << std::flush;
// 
//                 return -1;
//             }
//         }
//     }

    if (indexType == 0)
        return realMain<IndexSa<>>(options,
                                   TOutFormat(),
                                   BlastTabularSpecSelector<h>(),
                                   BlastProgramSelector<p>(),
                                   TRedAlph(),
                                   TScoreScheme(),
                                   TScoreExtension());
    else
        return realMain<TFMIndex<>>(options,
                                   TOutFormat(),
                                   BlastTabularSpecSelector<h>(),
                                   BlastProgramSelector<p>(),
                                   TRedAlph(),
                                   TScoreScheme(),
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
          typename TScoreScheme,
          typename TScoreExtension,
          typename TOutFormat,
          BlastProgram p,
          BlastTabularSpec h>
inline int
realMain(LambdaOptions                  const & options,
         TOutFormat                     const & /**/,
         BlastTabularSpecSelector<h>    const &,
         BlastProgramSelector<p>        const &,
         TRedAlph                       const & /**/,
         TScoreScheme                   const & /**/,
         TScoreExtension                const & /**/)
{
    using TGlobalHolder = GlobalDataHolder<TRedAlph,
                                           TScoreScheme,
                                           TIndexSpec,
                                           TOutFormat,
                                           p,
                                           h>;
    using TLocalHolder = LocalDataHolder<Match, TGlobalHolder, TScoreExtension>;

    myPrint(options, 1, "LAMBDA - the Local Aligner for Massive Biological DatA"
                      "\n======================================================"
                      "\nVersion ", SEQAN_APP_VERSION, "\n\n");

    if (options.verbosity >= 2)
        printOptions<TLocalHolder>(options);

    TGlobalHolder globalHolder;

    int ret = prepareScoring(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSubjects(globalHolder, options);
    if (ret)
        return ret;

    ret = loadDbIndexFromDisk(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSegintervals(globalHolder, options);
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

    open(globalHolder.outfile, toCString(options.output));
    context(globalHolder.outfile).fields = options.columns;
    writeHeader(globalHolder.outfile);

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

    // at least a block for each thread on double-indexing,
    // otherwise a block for each original query (contains 6 queries if
    // translation is used)
    uint64_t nBlocks = (options.doubleIndexing
                        ? options.queryPart
                        : length(globalHolder.qryIds));

    uint64_t lastPercent = 0;

    SEQAN_OMP_PRAGMA(parallel)
    {
        TLocalHolder localHolder(options, globalHolder);

        SEQAN_OMP_PRAGMA(for schedule(dynamic))
        for (uint64_t t = 0; t < nBlocks; ++t)
        {
            int res = 0;

            localHolder.init(t);

            // seed
            res = generateSeeds(localHolder);
            if (res)
                continue;

            if (options.doubleIndexing)
            {
                res = generateTrieOverSeeds(localHolder);
                if (res)
                    continue;
            }

            // search
            search(localHolder);

            // sort
            sortMatches(localHolder);

            // extend
            res = iterateMatches(localHolder);
            if (res)
                continue;


            if ((!options.doubleIndexing) && (TID == 0) &&
                (options.verbosity >= 1))
            {
                unsigned curPercent = ((t * 50) / nBlocks) * 2; // round to even
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

    writeFooter(globalHolder.outfile);

    if (!options.doubleIndexing)
    {
        myPrint(options, 2, "Runtime: ", sysTime() - start, "s.\n\n");
    }

    printStats(globalHolder.stats, options);

    return 0;
}
