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
// lambda.cpp: Main File for the main application
// ==========================================================================

// why is this neccessary?
// #undef SEQAN_HAS_ZLIB

// #define SEQAN_DEBUG_INDEX

// #define PARALLEL_SORT 0
// 0 = off
// 1 = GCC
// 2 = omptl

// #if PARALLEL_SORT == 1
//     #include <parallel/algorithm>
//     #define SORT __gnu_parallel::sort
// #elif PARALLEL_SORT == 2
//     #include <omptl/omptl_algorithm>
//     #define SORT omptl::sort
// #else
//     #define SORT std::sort
// #endif

#if defined(__GNUG__) && defined(_OPENMP)
    #define _GLIBCXX_PARALLEL
#endif
#define _GLIBCXX_USE_C99 1

#include <seqan/basic.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>

#include "lambda_indexer.hpp"

using namespace seqan;

// ==========================================================================
// Forwards
// ==========================================================================

inline int
argConv0(LambdaIndexerOptions const & options);

template <BlastProgram p>
inline int
argConv1(LambdaIndexerOptions           const & options,
         BlastProgramSelector<p>        const &);

template <BlastProgram p,
          typename TRedAlph>
inline int
argConv2(LambdaIndexerOptions     const & options,
         BlastProgramSelector<p>  const &,
         TRedAlph                 const &);

template <BlastProgram p,
          typename TRedAlph,
          typename TIndexSpecSpec>
inline int
realMain(LambdaIndexerOptions     const & options,
         BlastProgramSelector<p>  const &,
         TRedAlph                 const &,
         TIndexSpecSpec           const &);

// ==========================================================================
// Functions
// ==========================================================================
// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    LambdaIndexerOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    return argConv0(options);
}

inline int
argConv0(LambdaIndexerOptions const & options)
{
    switch(options.blastProgram)
    {
        case BlastProgram::BLASTN:
            return argConv1(options, BlastProgramSelector<BlastProgram::BLASTN>());
        case BlastProgram::BLASTP:
            return argConv1(options, BlastProgramSelector<BlastProgram::BLASTP>());
        case BlastProgram::BLASTX:
            return argConv1(options, BlastProgramSelector<BlastProgram::BLASTX>());
        case BlastProgram::TBLASTN:
            return argConv1(options, BlastProgramSelector<BlastProgram::TBLASTN>());
        case BlastProgram::TBLASTX:
            return argConv1(options, BlastProgramSelector<BlastProgram::TBLASTX>());
        default:
            break;
    }
    return -1;
}

/// Alphabet reduction
template <BlastProgram p>
inline int
argConv1(LambdaIndexerOptions           const & options,
         BlastProgramSelector<p>        const &)
{
    using TUnred = typename std::conditional<p == BlastProgram::BLASTN, Dna5, AminoAcid>::type;
    using Tp = BlastProgramSelector<p>;
    switch (options.alphReduction)
    {
        case 0:
            return argConv2(options, Tp(), TUnred());
        case 2:
            return argConv2(options, Tp(), ReducedAminoAcid<Murphy10>());
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
            return -1;
    }
    return -1;
}

template <BlastProgram p,
          typename TRedAlph>
inline int
argConv2(LambdaIndexerOptions     const & options,
         BlastProgramSelector<p>  const &,
         TRedAlph                 const &)
{

    if (options.algo == "mergesort")
        return realMain(options, BlastProgramSelector<p>(), TRedAlph(), SaAdvancedSort<MergeSortTag>());
    else if (options.algo == "quicksort")
        return realMain(options, BlastProgramSelector<p>(), TRedAlph(), SaAdvancedSort<QuickSortTag>());
    else if (options.algo == "quicksortbuckets")
        return realMain(options, BlastProgramSelector<p>(), TRedAlph(), SaAdvancedSort<QuickSortBucketTag>());
    else
        return realMain(options, BlastProgramSelector<p>(), TRedAlph(), Nothing());
}

template <BlastProgram p,
          typename TRedAlph,
          typename TIndexSpecSpec>
inline int
realMain(LambdaIndexerOptions     const & options,
         BlastProgramSelector<p>  const &,
         TRedAlph                 const &,
         TIndexSpecSpec           const &)
{
    using TOrigSet  = TCDStringSet<String<OrigSubjAlph<p>>>;
    using TTransSet = TCDStringSet<String<TransAlph<p>>>;

    TTransSet translatedSeqs;

    {
        TOrigSet originalSeqs;
        int ret = 0;

        // ids get saved to disk again immediately and are not kept in memory
        ret = loadSubjSeqsAndIds(originalSeqs, options);
        if (ret)
            return ret;

        // preserve lengths of untranslated sequences
        if (sIsTranslated(p))
            _saveOriginalSeqLengths(originalSeqs.limits, options);

        // convert the seg file to seqan binary format
        ret = convertMaskingFile(length(originalSeqs), options);
        if (ret)
            return ret;

        // translate or swap depending on program
        translateOrSwap(translatedSeqs, originalSeqs, options);
    }

    // dump translated and unreduced sequences
    dumpTranslatedSeqs(translatedSeqs, options);

    // see if final sequence set actually fits into index 
    if (!checkIndexSize(translatedSeqs))
        return -1;

    if (options.dbIndexType == 1)
    {
        using TIndexSpec = TFMIndex<TIndexSpecSpec>;
        generateIndexAndDump<TIndexSpec,TIndexSpecSpec>(translatedSeqs,
                                                        options,
                                                        BlastProgramSelector<p>(),
                                                        TRedAlph());
    } else
    {
        using TIndexSpec = IndexSa<TIndexSpecSpec>;
        generateIndexAndDump<TIndexSpec,TIndexSpecSpec>(translatedSeqs,
                                                        options,
                                                        BlastProgramSelector<p>(),
                                                        TRedAlph());
    }

    return 0;
}

