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
// options.h: contains the options and argument parser
// ==========================================================================


#ifndef SEQAN_LAMBDA_OPTIONS_H_
#define SEQAN_LAMBDA_OPTIONS_H_

#include <cstdio>
#include <unistd.h>

#include <seqan/basic.h>
#include <seqan/translation.h>
#include <seqan/blast.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>

#include <seqan/blast.h>

// ==========================================================================
// Metafunctions
// ==========================================================================

// suffix array overloads
namespace SEQAN_NAMESPACE_MAIN
{

template<typename TSpec1, typename TSpec2, typename TSpec3>
struct SAValue<StringSet<String<ReducedAminoAcid<TSpec1>, TSpec2>, TSpec3> >
{
    typedef Pair<uint32_t, uint16_t, Pack> Type;
};

// Dna Sequences might be longer
template<typename TSpec1, typename TSpec2>
struct SAValue<StringSet<String<Dna5, TSpec1>, TSpec2> >
{
    typedef Pair<uint32_t, uint32_t, Pack> Type;
};

}

using namespace seqan;

using TCDSpec = Packed<>;
//using TCDSpec = Alloc<>;

template <typename TAlph>
using TCDStringSet = StringSet<String<TAlph, TCDSpec>, Owner<ConcatDirect<> > >;


template <BlastFormatProgram p>
struct UnreducedStringSet
{
    typedef TCDStringSet<AminoAcid> Type;
};

template <>
struct UnreducedStringSet<BlastFormatProgram::BLASTN>
{
    typedef TCDStringSet<Dna5> Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LambdaOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.


struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    CharString  dbFile;

    // for indexer, the file format of database sequences
    // for main app, the file format of query sequences
    // 0 -- fasta, 1 -- fastq
    int         fileFormat; 

    int         alphReduction;

    GeneticCodeSpec geneticCode;

    BlastFormatProgram blastProg;

    bool        isTerminal = false;
    unsigned    terminalCols = 80;

    SharedOptions() :
        verbosity(2),fileFormat(0),alphReduction(0),
        blastProg(BlastFormatProgram::BLASTX)
    {
        isTerminal = isatty(fileno(stdout));
        if (isTerminal)
        {
            int cols = 80;
        #ifdef TIOCGSIZE
            struct ttysize ts;
            ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
            cols = ts.ts_cols;
        #elif defined(TIOCGWINSZ)
            struct winsize ts;
            ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
            cols = ts.ws_col;
        #endif /* TIOCGSIZE */
            terminalCols = cols;
        }
    }
};


struct LambdaOptions : public SharedOptions
{

    CharString      queryFile;
    bool            revComp     = true;

    CharString      output;

    unsigned short  queryPart = 0;

//     bool            semiGlobal;

    unsigned char   seedLength  = 0;
    unsigned char   maxSeedDist = 1;
    bool            hammingOnly = false;

    char            seedGravity     = 0;
    unsigned char   seedOffset      = 0;
    unsigned char   minSeedLength   = 0;
    unsigned char   minSeedScore    = 0;
    unsigned int    minSeedEVal     = 0;
    double          minSeedBitS     = -1;

    // 0 = manual, positive X = blosumX, negative Y = pamY
    int             scoringMethod   = 62;
    // scores
    int             gapOpen         = -11;
    int             gapExtend       = -1;
    int             match           = 0; // only for manual
    int             misMatch        = 0; // only for manual

    int             xDropOff    = 0;
    int             band        = -1;
    double          eCutOff     = 0;

    unsigned        threads     = 1; // not exposed
    LambdaOptions() :
        SharedOptions()
    {
        #if defined(_OPENMP)
        threads = omp_get_max_threads();
        #endif
    }
};

struct LambdaIndexerOptions : public SharedOptions
{
    CharString      segFile;

    LambdaIndexerOptions()
        : SharedOptions(), segFile("")
    {}
};



// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("lambda");
    // Set short description, version, and date.
    setShortDescription(parser, "BLAST compatible local aligner optimized for "
                                "NGS and Metagenomics.");
    setVersion(parser, "0.2");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP "
                         "\\fI-d DATABASE.fasta\\fP "
                         "[\\fI-o output.m8\\fP]");
    addDescription(parser, "Lambda is a local aligner that is faster than "
                           "BLAST, optimized for many query sequences.");

    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("q",
                                            "query",
                                            "Query sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "q");
    setValidValues(parser, "query", "fasta fa fna faa fastq fq");



    addOption(parser, seqan::ArgParseOption("d",
                                            "database",
                                            "Database sequences (fasta), with "
                                             "precomputed index (.sa).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "d");
    setValidValues(parser, "database", "fasta fa fna faa");

    // TODO implement some way to transfer options from index builder
    // maybe write a spec file that is loaded by this app

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "File to hold reports on hits (.m8 "
                                            "is blastall -m8 et cetera)",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "m0 m8 m9");
    setDefaultValue(parser, "output", "output.m8");


    addSection(parser, "Program Options");
    addOption(parser, seqan::ArgParseOption("p",
                                            "program",
                                            "Blast Operation Mode.",
                                            seqan::ArgParseArgument::STRING,
                                            "OUT"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

    addOption(parser, seqan::ArgParseOption("pf",
                                            "partition-factor",
                                            "The query sequences "
                                            "are partioned into pf * n parts "
                                            "and searched en bloc, where n is "
                                            "the number of available threads. "
                                            "A pf of 1 yields the best speed, "
                                            "a higher pf results in slight "
                                            "run-time increases, but "
                                            "dramatically decreases memory "
                                            "usage.",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "partition-factor", "2");

    addSection(parser, "Translation and Alphabet");
    addOption(parser, seqan::ArgParseOption("gc",
                                            "genetic-code",
                                            "The translation table to use "
                                            "for nucl -> amino acid translation"
                                            "(not for BlastN, BlastP). See "
               "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
                                            " for ids (default is generic)",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "Alphabet Reduction for AminoAcid "
                                            "Alphabet (0 -> off; 2 -> Murphy10; "
                                            "8,10,12 -> Lambda*).",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "alph", "2");


    addSection(parser, "Seeding / Filtration");
    addOption(parser, seqan::ArgParseOption("su",
                                            "ungapped-seeds",
                                            "Deacitvate EditDistance-Seeding.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "ungapped-seeds", "1");

    addOption(parser, seqan::ArgParseOption("sl",
                                            "seed-length",
                                            "Length of the seeds (0 -> choose "
                                            "automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-length", "0");

    addOption(parser, seqan::ArgParseOption("so",
                                            "seed-offset",
                                            "Offset for seeding (by default "
                                            "same as length -> non-overlapping).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-offset", "0");

    addOption(parser, seqan::ArgParseOption("sd",
                                            "seed-delta",
                                            "maximum seed distance.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-delta", "1");

    addOption(parser, seqan::ArgParseOption("sg",
                                            "seed-gravity",
                                            "Seeds closer than this are joined"
                                            "(-1 -> choose automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-gravity", "-1");

    addOption(parser, seqan::ArgParseOption("sm",
                                            "seed-min-length",
                                            "after postproc shorter seeds are "
                                            "discarded"
                                            "(0 -> choose automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-length", "0");

    addOption(parser, seqan::ArgParseOption("ss",
                                            "seed-min-score",
                                            "after postproc worse seeds are "
                                            "discarded [raw score].",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-score", "32");

//     addOption(parser, seqan::ArgParseOption("se",
//                                             "seedminevalue",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(0 -> off).",
//                                             seqan::ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "seedminevalue", "100000");

//     addOption(parser, seqan::ArgParseOption("sb",
//                                             "seedminbits",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(-1 -> off).",
//                                             seqan::ArgParseArgument::DOUBLE));
//     setDefaultValue(parser, "seedminbits", "-1");

    addSection(parser, "Scoring");

    addOption(parser, seqan::ArgParseOption("sc",
                                            "scoring-scheme",
                                            "'62' for Blosum62 (default); '50' for "
                                            "Blosum50; '0' for manual (default "
                                            "for BlastN)",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "scoring-scheme", "62");

    addOption(parser, seqan::ArgParseOption("ge",
                                            "score-gap",
                                            "Score per gap character.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap", "-1");

    addOption(parser, seqan::ArgParseOption("go",
                                            "score-gap-open",
                                            "Additional cost for opening gap.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap-open", "-11");

    addOption(parser, seqan::ArgParseOption("ma",
                                            "score-match",
                                            "Match score (manual scoring only)",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-match", "2");

    addOption(parser, seqan::ArgParseOption("mi",
                                            "score-mismatch",
                                            "Mismatch score (manual scoring only)",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-mismatch", "-3");

    addSection(parser, "Extension");

    addOption(parser, seqan::ArgParseOption("x",
                                            "x-drop",
                                            "Stop Banded extension if score "
                                            "x below the maximum seen"
                                            "(-1 means no xdrop).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "x-drop", "30");

    addOption(parser, seqan::ArgParseOption("b",
                                            "band",
                                            "Size of the DP-band used in "
                                            "extension (-3 means log2 "
                                            "of query length; -2 means sqrt of "
                                            "query length; -1 means full dp; "
                                            "n means band of size 2n+1)",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "band", "-3");

    addOption(parser, seqan::ArgParseOption("e",
                                            "e-value",
                                            "Minimum E-Value for Results.",
                                            seqan::ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e-value", "0.1");


    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");
    addListItem(parser, "\\fBOMP_NUM_THREADS\\fP",
                        "number of threads to use, defaults to number of "
                        "CPUs/-cores");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    getOptionValue(options.dbFile, parser, "database");

    CharString blastMode;
    getOptionValue(blastMode, parser, "program");
    if (blastMode == "blastn")
    {
        options.blastProg   = BlastFormatProgram::BLASTN;
        options.seedLength  = 20;
        options.alphReduction = 0;
    }
    else if (blastMode == "blastp")
    {
        options.blastProg = BlastFormatProgram::BLASTP;
        options.seedLength  = 10;
    }
    else if (blastMode == "blastx")
    {
        options.blastProg = BlastFormatProgram::BLASTX;
        options.seedLength  = 10;
    }
    else if (blastMode == "tblastn")
    {
        options.blastProg = BlastFormatProgram::TBLASTN;
        options.seedLength  = 10;
    }
    else if (blastMode == "tblastx")
    {
        options.blastProg = BlastFormatProgram::TBLASTX;
        options.seedLength  = 10;
    }
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    //TODO adapt seedlength if alphReduction

    int buf = 0;
    getOptionValue(buf, parser, "seed-length");
    if (buf != 0) // seed length was specified manually
        options.seedLength = buf;

    buf = 0;
    getOptionValue(buf, parser, "seed-offset");
    if (buf != 0) // seed length was specified manually
        options.seedOffset = buf;
    else
        options.seedOffset = options.seedLength;

    buf = 0;
    getOptionValue(buf, parser, "seed-delta");
    options.maxSeedDist = buf;
    if (buf == 0)
    {
        options.hammingOnly = true;
    } else
    {
        getOptionValue(buf, parser, "ungapped-seeds");
        options.hammingOnly = (buf != 0);
        if (!options.hammingOnly)
        {
            std::cerr << "Edit-Distance seeds not supported in the build.\n";
            return seqan::ArgumentParser::PARSE_ERROR;
        }
    }


    getOptionValue(buf, parser, "seed-gravity");
    if (buf == -1)
        options.seedGravity = options.seedLength;
    else
        options.seedGravity = buf;

    unsigned foo = 0;
    getOptionValue(foo, parser, "seed-min-length");
    options.minSeedLength = foo;
    if (options.minSeedLength == 0)
        options.minSeedLength = options.seedLength;

    getOptionValue(foo, parser, "seed-min-score");
    options.minSeedScore = foo;

//     getOptionValue(foo, parser, "seedminevalue");
//     options.minSeedEVal = foo;

//     getOptionValue(options.minSeedBitS, parser, "seedminbits");

    getOptionValue(options.eCutOff, parser, "e-value");

    getOptionValue(options.xDropOff, parser, "x-drop");
    if (options.xDropOff < -1)
    {
        std::cerr << "xDropOff must be in [ -1 : MAXINT ]" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.band, parser, "band");
    if (options.band < -3)
    {
        std::cerr << "band must be in [ -3 : MAXINT ]" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    //verifyFileFormat(); TODO verify fileformats for program mode

    getOptionValue(options.output, parser, "output");

    getOptionValue(options.alphReduction, parser, "alph");
    switch (options.alphReduction)
    {
        case 1: case 2: case 8: case 10: case 12:
        if (options.blastProg == BlastFormatProgram::BLASTN)
        {
            std::cerr << "Alphabet Reduction only makes sense"
                            " for proteins." << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        } break;
        case 0: break;
        break;
        default:
            std::cerr << "Alphabet Reduction must be one of 2, 8, 10, 12 or off" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(foo, parser, "partition-factor");
#ifdef _OPENMP
    options.queryPart = foo * omp_get_max_threads();
#else
    options.queryPart = 1;
#endif

    getOptionValue(options.scoringMethod, parser, "scoring-scheme");
    switch (options.scoringMethod)
    {
        case 0:
            getOptionValue(options.misMatch, parser, "score-mismatch");
            getOptionValue(options.match, parser, "score-match");
            break;
        case 45: case 62: case 80: break;
        default:
            std::cerr << "Unsupported Scoring Scheme selected.\n";
            return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.gapExtend, parser, "score-gap");
    getOptionValue(options.gapOpen, parser, "score-gap-open");

    getOptionValue(foo, parser, "genetic-code");
    switch (foo)
    {
        case 1: case 2: case 3: case 4: case 5: case 6:
        case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        case 21: case 22: case 23: case 24 : case 25:
            options.geneticCode = static_cast<GeneticCodeSpec>(foo);
            break;
        default:
            std::cerr << "Invalid genetic code. See trans_table vars at "
                      << "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
                      << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }
    return seqan::ArgumentParser::PARSE_OK;
}


seqan::ArgumentParser::ParseResult
parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("lambda_indexer");
    // Set short description, version, and date.
    setShortDescription(parser, "Indexer for Lambda");
    setVersion(parser, "0.2");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-i DATABASE.fasta\\fP [\\-o INDEXFILE.sa\\fP]");
    addDescription(parser, "This is the application skelleton and you should modify this string.");


    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("i",
                                            "input",
                                            "Database sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "i");
    setValidValues(parser, "input", "fasta fa fna faa");

    addOption(parser, seqan::ArgParseOption("s",
                                            "segfile",
                                            "SEG intervals for database"
                                            "(optional).",
                                            seqan::ArgParseArgument::INPUTFILE));

    setValidValues(parser, "segfile", "seg");


    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "Index of database sequences",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "sa fm");

    addSection(parser, "Program Options");
    addOption(parser, seqan::ArgParseOption("p",
                                            "program",
                                            "Blast Operation Mode.",
                                            seqan::ArgParseArgument::STRING,
                                            "OUT"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

    addSection(parser, "Translation and Alphabet");
    addOption(parser, seqan::ArgParseOption("gc",
                                            "genetic-code",
                                            "The translation table to use "
                                            "(not for BlastN, BlastP). See "
               "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
                                            " for ids (default is generic)",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "Alphabet Reduction for AminoAcid "
                                            "Alphabet (0 -> off; 2 -> Murphy10; "
                                            "8,10,12 -> Lambda*).",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "alph", "2");

    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Remarks");
    addText(parser, "Note that the indeces created are binary and not "
                    "necessarily platform independant.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.

    seqan::getOptionValue(options.dbFile, parser, "input");

    CharString blastMode;
    getOptionValue(blastMode, parser, "program");
    if (blastMode == "blastn")
    {
        options.blastProg   = BlastFormatProgram::BLASTN;
        options.alphReduction = 0;

    }
    else if (blastMode == "blastp")
    {
        options.blastProg = BlastFormatProgram::BLASTP;

    }
    else if (blastMode == "blastx")
    {
        options.blastProg = BlastFormatProgram::BLASTX;

    }
    else if (blastMode == "tblastn")
    {
        options.blastProg = BlastFormatProgram::TBLASTN;

    }
    else if (blastMode == "tblastx")
    {
        options.blastProg = BlastFormatProgram::TBLASTX;

    }
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    //verifyFileFormat(); TODO

    getOptionValue(options.alphReduction, parser, "alph");
    switch (options.alphReduction)
    {
        case 0: case 1: case 2: case 8: case 10: case 12:
        if (options.blastProg == BlastFormatProgram::BLASTN)
        {
            std::cerr << "Alphabet Reduction only makes sense"
                            " for proteins." << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }
        break;
        default:
            std::cerr << "Alphabet Reduction must be one of 2, 8, 10, 12 or off" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }


    getOptionValue(options.segFile, parser, "segfile");

    int foo = 0;
    getOptionValue(foo, parser, "genetic-code");
    switch (foo)
    {
        case 1: case 2: case 3: case 4: case 5: case 6:
        case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        case 21: case 22: case 23: case 24 : case 25:
            options.geneticCode = static_cast<GeneticCodeSpec>(foo);
            break;
        default:
            std::cerr << "Invalid genetic code. See trans_table vars at "
                      << "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
                      << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }

    return seqan::ArgumentParser::PARSE_OK;
}



inline void
printOptions(LambdaOptions const & options)
{
    std::cout << "Seed-Length:   " << uint(options.seedLength) << "\n"
              << "Seed-Offset:   " << uint(options.seedOffset) << "\n"
              << "Seed-Delta:    " << uint(options.maxSeedDist) << "\n"
              << "Hammingonly:   " << uint(options.hammingOnly) << "\n"
              << "SeedGravity:   " << uint(options.seedGravity) << "\n"
              << "MinSeedLength: " << uint(options.minSeedLength) << "\n"
              << "MinSeedScore:  " << uint(options.minSeedScore) << "\n"
              << "MinSeedEval:   " << uint(options.minSeedEVal) << "\n"
              << "MinSeedBitS:   " << options.minSeedBitS << "\n"
              << "xDropOff:      " << options.xDropOff << "\n"
              << "band:          " << options.band << "\n"
              << "eCutOff:       " << options.eCutOff << "\n"
              << "alphReduction: " << uint(options.alphReduction) << "\n"
              << "queryParts:    " << uint(options.queryPart) << "\n";
}


#endif // header guard
