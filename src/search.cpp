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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#include "seqan2seqan3.hpp" // must come first

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

void argConv0(LambdaOptions & options);

template <DbIndexType c_indexType>
void argConv1(LambdaOptions const & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph>
void argCon2a(LambdaOptions const & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph>
void argCon2b(LambdaOptions const & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph>
void argCon3(LambdaOptions     const & options);

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void realMain(LambdaOptions     const & options);

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int searchMain(int const argc, char const ** argv)
{
    LambdaOptions options;
    parseCommandLine(options, argc, argv);

#ifdef _OPENMP
    omp_set_num_threads(options.threads);
#else
    options.threads = 1;
#endif

#ifdef NDEBUG
    try
    {
       argConv0(options);
    } catch (std::bad_alloc const & e)
    {
        std::cerr << "\n\nERROR: Lambda ran out of memory :(\n"
                     "       You need to split your file into smaller segments or search against a smaller database.\n";
        return -1;
    } catch (IndexException const & e)
    {
        std::cerr << "\n\nERROR: The following exception was thrown while reading the index:\n"
                  <<     "       \"" << e.what() << "\"\n"
                  <<     "       Make sure the directory exists and is readable; recreate the index and try again.\n"
                  <<     "       If the problem persists, report an issue at https://github.com/seqan/lambda/issues "
                  << "and include this output, as well as the output of `lambda2 --version`, thanks!\n";
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

// CONVERT Run-time options to compile-time Format-Type
void
argConv0(LambdaOptions & options)
{
    myPrint(options, 1, "LAMBDA - the Local Aligner for Massive Biological DatA"
                        "\n======================================================"
                        "\nVersion ", SEQAN_APP_VERSION, "\n\n");

    // Index
    myPrint(options, 1, "Reading index properties... ");
    readIndexOptions(options);
    myPrint(options, 1, "done.\n");

    myPrint(options, 2, "  type:                ", _indexEnumToName(options.indexFileOptions.indexType), "\n",
                        "  original alphabet:   ", _alphabetEnumToName(options.indexFileOptions.origAlph), "\n");
    if (options.indexFileOptions.origAlph == options.indexFileOptions.transAlph)
    {
        myPrint(options, 2, "  translated alphabet: not translated\n");
        if ((int)options.geneticCodeQry == 0) // use same geneticCode as Index, but index wasn't translated
            options.geneticCodeQry = seqan3::genetic_code::CANONICAL;
    }
    else
    {
        myPrint(options, 2, "  translated alphabet: ", _alphabetEnumToName(options.indexFileOptions.transAlph), "\n");
        myPrint(options, 2, "  ^- genetic code:     ", (int)options.indexFileOptions.geneticCode, "\n");
        if ((int)options.geneticCodeQry == 0) // use same geneticCode as Index
        {
            options.geneticCodeQry = options.indexFileOptions.geneticCode;
        } else if (options.geneticCodeQry != options.indexFileOptions.geneticCode)
        {
            std::cerr << "WARNING: The genetic code used when creating the index: " << (int) options.indexFileOptions.geneticCode
                      << "\n         is not the same as now selected for the query sequences: " << (int)options.geneticCodeQry
                      << "\n         Are you sure this is what you want?\n";
        }
    }

    if (options.indexFileOptions.transAlph == options.indexFileOptions.redAlph)
    {
        myPrint(options, 2, "  reduced alphabet:    not reduced\n");
    }
    else
    {
        myPrint(options, 2, "  reduced alphabet:    ", _alphabetEnumToName(options.indexFileOptions.redAlph), "\n\n");
    }

    if ((options.nucleotide_mode) && (options.indexFileOptions.redAlph != AlphabetEnum::DNA5 && options.indexFileOptions.redAlph != AlphabetEnum::DNA4))
    {
        throw std::runtime_error("You are attempting a nucleotide search on a protein index. "
                                 "Did you want to use 'lambda3 searchp' instead?");
    }

    // query file
    if (options.qryOrigAlphabet == AlphabetEnum::DNA4) // means "auto", as dna4 not valid as argument to --query-alphabet
    {
        myPrint(options, 1, "Detecting query alphabet... ");
        options.qryOrigAlphabet = detectSeqFileAlphabet(options.queryFile);
        myPrint(options, 1, _alphabetEnumToName(options.qryOrigAlphabet), " detected.\n");
        myPrint(options, 2, "\n");
    }

#if 0 // blastProgram set in GlobalDataHolder later
    // set blastProgram
    if (options.blastProgram == seqan::BlastProgram::UNKNOWN)
    {
        if ((options.indexFileOptions.transAlph == AlphabetEnum::DNA5) && (options.qryOrigAlphabet == AlphabetEnum::AMINO_ACID))
        {
            throw IndexException("Query file is protein, but index is nucleotide. "
                                 "Recreate the index with 'lambda mkindexp'.");
        }
        else if ((options.indexFileOptions.transAlph == AlphabetEnum::DNA5) && (options.qryOrigAlphabet == AlphabetEnum::DNA5))
        {
            options.blastProgram = seqan::BlastProgram::BLASTN;
        }
        else if (options.qryOrigAlphabet == AlphabetEnum::DNA5) // query will be translated
        {
            if (options.indexFileOptions.origAlph == options.indexFileOptions.transAlph)
                options.blastProgram = seqan::BlastProgram::BLASTX;
            else
                options.blastProgram = seqan::BlastProgram::TBLASTX;
        }
        else // query is aminoacid already
        {
            if (options.indexFileOptions.origAlph == options.indexFileOptions.transAlph)
                options.blastProgram = seqan::BlastProgram::BLASTP;
            else
                options.blastProgram = seqan::BlastProgram::TBLASTN;
        }
    }

#endif

#if 0 // this needs to be moved
    // some blastProgram-specific "late option modifiers"
    if (((options.blastProgram == seqan::BlastProgram::BLASTP) ||
         (options.blastProgram == seqan::BlastProgram::TBLASTN)) &&
        (!options.samBamTags[SamBamExtraTags<>::Q_AA_CIGAR]))
        options.samBamSeq = 0;
#endif
    // sizes
    checkRAM(options);

    switch (options.indexFileOptions.indexType)
    {
        case DbIndexType::FM_INDEX:     return argConv1<DbIndexType::FM_INDEX>(options);
        case DbIndexType::BI_FM_INDEX:  return argConv1<DbIndexType::BI_FM_INDEX>(options);
        default: throw 52;
    }
}

template <DbIndexType c_indexType>
void argConv1(LambdaOptions const & options)
{
    if (options.nucleotide_mode)
    {
        return argCon2a<c_indexType, AlphabetEnum::DNA5>(options);
    }
    else
    {
        switch (options.indexFileOptions.origAlph)
        {
            case AlphabetEnum::DNA5:        return argCon2b<c_indexType,
                                                            AlphabetEnum::DNA5>(options);
            case AlphabetEnum::AMINO_ACID:  return argCon2b<c_indexType,
                                                            AlphabetEnum::AMINO_ACID>(options);
            default: throw 53;
        }
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph>
void argCon2a(LambdaOptions const & options)
{
    // transalph is always amino acid, unless in nucleotide_mode
    switch (options.indexFileOptions.redAlph)
    {
        case AlphabetEnum::DNA5:      return realMain<c_indexType,
                                                      c_origSbjAlph,
                                                      AlphabetEnum::DNA5,
                                                      AlphabetEnum::DNA5,
                                                      AlphabetEnum::DNA5>(options);
        case AlphabetEnum::DNA4:      return realMain<c_indexType,
                                                      c_origSbjAlph,
                                                      AlphabetEnum::DNA5,
                                                      AlphabetEnum::DNA4,
                                                      AlphabetEnum::DNA5>(options);
        default: throw 555;
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph>
void argCon2b(LambdaOptions const & options)
{
    // transalph is always amino acid, unless in nucleotide_mode
    switch (options.indexFileOptions.redAlph)
    {
        case AlphabetEnum::AMINO_ACID:      return argCon3<c_indexType,
                                                            c_origSbjAlph,
                                                            AlphabetEnum::AMINO_ACID,
                                                            AlphabetEnum::AMINO_ACID>(options);
        case AlphabetEnum::MURPHY10:        return argCon3<c_indexType,
                                                            c_origSbjAlph,
                                                            AlphabetEnum::AMINO_ACID,
                                                            AlphabetEnum::MURPHY10>(options);
        case AlphabetEnum::LI10:            return argCon3<c_indexType,
                                                            c_origSbjAlph,
                                                            AlphabetEnum::AMINO_ACID,
                                                            AlphabetEnum::LI10>(options);
        default: throw 54;
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph>
void argCon3(LambdaOptions     const & options)
{
    // transalph is always amino acid, unless in nucleotide_mode
    switch (options.qryOrigAlphabet)
    {
        case AlphabetEnum::DNA5:            return realMain<c_indexType,
                                                            c_origSbjAlph,
                                                            c_transAlph,
                                                            c_redAlph,
                                                            AlphabetEnum::DNA5>(options);
        case AlphabetEnum::AMINO_ACID:      return realMain<c_indexType,
                                                            c_origSbjAlph,
                                                            c_transAlph,
                                                            c_redAlph,
                                                            AlphabetEnum::AMINO_ACID>(options);
        default: throw 55;
    }
}

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void realMain(LambdaOptions     const & options)
{
    using TGlobalHolder = GlobalDataHolder<c_indexType,
                                           c_origSbjAlph,
                                           c_transAlph,
                                           c_redAlph,
                                           c_origQryAlph>;
    using TLocalHolder = LocalDataHolder<TGlobalHolder>;

    if (options.verbosity >= 2)
        printOptions<TLocalHolder>(options);

    TGlobalHolder globalHolder;

    prepareScoring(globalHolder, options);

    loadDbIndexFromDisk(globalHolder, options);

    loadQuery(globalHolder, options);

    myWriteHeader(globalHolder, options);

    myPrint(options, 1, "Searching and extending hits on-line...progress:\n"
                "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n|");

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
            search(localHolder); //TODO seed refining if iterateMatches gives 0 results
        #ifdef LAMBDA_MICRO_STATS
            localHolder.stats.timeSearch += sysTime() - buf;
        #endif

            // extend
            if (localHolder.matches.size() > 0)
                res = iterateMatches(localHolder);

            if (res)
                continue;

            if ((omp_get_thread_num() == 0) && (options.verbosity >= 1))
            {
                unsigned curPercent = ((t * 50) / localHolder.nBlocks) * 2; // round to even
                printProgressBar(lastPercent, curPercent);
            }

        } // implicit thread sync here

        if ((omp_get_thread_num() == 0) && (options.verbosity >= 1))
            printProgressBar(lastPercent, 100);

        SEQAN_OMP_PRAGMA(critical(statsAdd))
        {
            globalHolder.stats += localHolder.stats;
        }
    }

    myPrint(options, 1, "\n");

    myWriteFooter(globalHolder, options);

    myPrint(options, 2, "Runtime total: ", sysTime() - start, "s.\n\n");

    printStats(globalHolder.stats, options);

}
