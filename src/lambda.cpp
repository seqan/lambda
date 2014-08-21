// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// lambda.cpp: Main File for Lambda
// ==========================================================================

#undef SEQAN_HAS_ZLIB

#define FASTBUILD

// #define SEQAN_DEBUG_INDEX


// #define PARALLEL_SORT 0
// // 0 = off
// // 1 = GCC
// // 2 = omptl
// 
// #if PARALLEL_SORT == 1
//     #include <parallel/algorithm>
//     #define SORT __gnu_parallel::sort
// #elif PARALLEL_SORT == 2
//     #include <omptl/omptl_algorithm>
//     #define SORT omptl::sort
// #else
//     #define SORT std::sort
// #endif

#define SEQAN_CXX11_STANDARD

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/reduced_aminoacid.h>

#include "options.hpp"
#include "match.hpp"
#include "lambda.hpp"
#include "misc.hpp"

using namespace seqan;


inline BlastFormatFile
_fileType(LambdaOptions const & options)
{
    if (endsWith(options.output, ".m0"))
        return BlastFormatFile::PAIRWISE;
    else if (endsWith(options.output, ".m8"))
        return BlastFormatFile::TABULAR;
    else if (endsWith(options.output, ".m9"))
        return BlastFormatFile::TABULAR_WITH_HEADER;
    else
        return BlastFormatFile::INVALID_File;
}



// forwards

inline int
argConv0(LambdaOptions const & options);
//-
template <BlastFormatFile m,
          BlastFormatGeneration g>
inline int
argConv1(LambdaOptions                               const & options,
         BlastFormat<m,BlastFormatProgram::BLASTN,g> const & /**/);
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename std::enable_if<p != BlastFormatProgram::BLASTN, ns_enabler::enabler>::type...>
inline int
argConv1(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph>
inline int
argConv2(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
argConv3(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension>
inline int
realMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/);

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

    return argConv0(options);
}


// CONVERT Run-time options to compile-time Format-Type
inline int
argConv0(LambdaOptions const & options)
{
    switch (options.blastProg)
    {
#ifndef FASTBUILD
         case BlastFormatProgram::BLASTN :
         {
 //             template <BlastFormatFile m>
 //             using TBF = BlastFormat<m,
 //                                 BlastFormatProgram::BLASTN,
 //                                 BlastFormatGeneration::BLAST>;
             switch (_fileType(options))
             {
                 case BlastFormatFile::PAIRWISE:
                 {
                     typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 case BlastFormatFile::TABULAR:
                 {
                     typedef BlastFormat<BlastFormatFile::TABULAR,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 case BlastFormatFile::TABULAR_WITH_HEADER:
                 {
                     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 default:
                     break;
             }
         } break;

        case BlastFormatProgram::BLASTP :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::BLASTP,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;
#endif
        case BlastFormatProgram::BLASTX :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::BLASTX,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
// #ifndef FASTBUILD
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #endif
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #ifndef FASTBUILD
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #endif
                default:
                    break;
            }
        } break;
#ifndef FASTBUILD
        case BlastFormatProgram::TBLASTN :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::TBLASTN,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;

        case BlastFormatProgram::TBLASTX :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::TBLASTX,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;
#endif
        default:
            break;
    }
    return -1;
}

/// Alphabet reduction

template <BlastFormatFile m,
          BlastFormatGeneration g>
inline int
argConv1(LambdaOptions                               const & options,
         BlastFormat<m,BlastFormatProgram::BLASTN,g> const & /**/)
{
    using TFormat = BlastFormat<m,BlastFormatProgram::BLASTN,g>;
    return argConv2(options, TFormat(), Dna5());
}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename std::enable_if<p != BlastFormatProgram::BLASTN, ns_enabler::enabler>::type...>
inline int
argConv1(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
    switch (options.alphReduction)
    {

        case 0:
            return argConv2(options, TFormat(), AminoAcid());
        case 2:
            return argConv2(options, TFormat(), ReducedAminoAcid<Murphy10>());
        case 10:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<10>>());
#ifndef FASTBUILD
        case 1:
            return argConv2(options, TFormat(), AminoAcid10());
        case 8:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<8>>());
        case 12:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<12>>());
#endif
        default:
            return -1;
    }
    return -1;
}

/// scoring scheme

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph>
inline int
argConv2(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
    switch (options.scoringMethod)
    {
// #ifndef FASTBUILD
        case 0:
            return argConv3(options, TFormat(), TRedAlph(), Score<int, Simple>());
        case 45:
            return argConv3(options, TFormat(), TRedAlph(), Blosum45());
        case 80:
            return argConv3(options, TFormat(), TRedAlph(), Blosum80());
// #endif
        case 62:
            return argConv3(options, TFormat(), TRedAlph(), Blosum62());
        default:
            std::cerr << "Unsupported Scoring Scheme selected.\n";
            break;
    }

    return -1;
}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
argConv3(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
#ifndef FASTBUILD
    if (options.gapOpen == 0)
        return realMain(options,
                        TFormat(),
                        TRedAlph(),
                        TScoreScheme(),
                        LinearGaps());
    else
#endif
        return realMain(options,
                        TFormat(),
                        TRedAlph(),
                        TScoreScheme(),
                        AffineGaps());

}

/// REAL MAIN

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension>
inline int
realMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/)
{
    using TGlobalHolder = GlobalDataHolder<TRedAlph,
                                           TScoreScheme,
                                           m, p, g>;
    using TLocalHolder = LocalDataHolder<Match, TGlobalHolder, TScoreExtension>;

    if (options.verbosity >= 2)
        printOptions(options);

    TGlobalHolder globalHolder;
    globalHolder.stats.clear();

    int ret = prepareScoring(globalHolder, options);
    if (ret)
        return ret;

    ret = loadDbIndexFromDisk(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSubjects(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSegintervals(globalHolder, options);
    if (ret)
        return ret;

    ret = loadQuery(globalHolder, options);
    if (ret)
        return ret;

    std::ofstream stream;
    stream.open(toCString(options.output));
    if (!stream.is_open())
    {
        std::cerr << "Error opening output file for writing." << std::endl;
        return -2;
    }

    ret = writeTop(stream,
                   globalHolder.dbSpecs,
                   typename TGlobalHolder::TFormat());
    if (ret)
        return ret;

    myPrint(options, 1,
            "Searching ",
            options.queryPart,
            " blocks of query with ",
            omp_get_max_threads(),
            " threads...\n");

    if (isatty(fileno(stdin)))
    {
        for (unsigned char i=0; i< omp_get_num_threads()+3; ++i)
            std::cout << std::endl;
        std::cout << "\033[" << omp_get_num_threads()+2 << "A";
    }

    // TODO evaluate localHolder outside of loop and firstprivate
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (unsigned short t = 0; t < options.queryPart; ++t)
    {
        int res = 0;
        TLocalHolder localHolder(options, globalHolder);
        localHolder.init(t);

        // seed and searchPipeline
        res = generateSeeds(localHolder);
        if (res)
            continue;

        res = generateTrieOverSeeds(localHolder);
        if (res)
            continue;

        // search
        res = search(localHolder);
        if (res)
            continue;

        // sort, joing and filter
        joinAndFilterMatches(localHolder);

        // extend
        res = iterateMatches(stream, localHolder);
        if (res)
            continue;

        if (options.verbosity >= 2)
        {
            #pragma omp critical(statsAdd)
            {
                globalHolder.stats += localHolder.stats;
            }
        }
    }

    ret = writeBottom(stream,
                      globalHolder.dbSpecs,
                      globalHolder.blastScoringAdapter,
                      typename TGlobalHolder::TFormat());

    stream.close();

    if (ret)
        return ret;

    printStats(globalHolder.stats, options);

    return 0;
}
