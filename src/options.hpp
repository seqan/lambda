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

template<typename TSpec2, typename TSpec3>
struct SAValue<StringSet<String<AminoAcid, TSpec2>, TSpec3> >
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

#if defined LAMBDA_BITCOPMRESSED_STRINGS
using TCDSpec = Packed<>;
#else
using TCDSpec = Alloc<>;
#endif

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


template <BlastFormatProgram p>
using OrigQryAlph = typename std::conditional<
                                           (p == BlastFormatProgram::BLASTN) ||
                                           (p == BlastFormatProgram::BLASTX) ||
                                           (p == BlastFormatProgram::TBLASTX),
                                           Dna5,
                                           AminoAcid>::type;

template <BlastFormatProgram p>
using OrigSubjAlph = typename std::conditional<
                                           (p == BlastFormatProgram::BLASTN) ||
                                           (p == BlastFormatProgram::TBLASTN) ||
                                           (p == BlastFormatProgram::TBLASTX),
                                           Dna5,
                                           AminoAcid>::type;

template <BlastFormatProgram p>
using TransAlph = typename std::conditional<(p == BlastFormatProgram::BLASTN),
                                            Dna5,
                                            AminoAcid>::type;


template <BlastFormatProgram p, typename TRedAlph_>
using RedAlph = typename std::conditional<(p == BlastFormatProgram::BLASTN),
                                          Dna5,
                                          TRedAlph_>::type;



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
    int8_t verbosity = 2;


    CharString  dbFile;

    int8_t      dbIndexType = 0;
    // for indexer, the file format of database sequences
    // for main app, the file format of query sequences
    // 0 -- fasta, 1 -- fastq
    int8_t      fileFormat = 0;

    int      alphReduction = 0;

    GeneticCodeSpec geneticCode = CANONICAL;

    BlastFormatProgram blastProg = BlastFormatProgram::BLASTX;

    bool        isTerm = true;
    unsigned    terminalCols = 80;

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


struct LambdaOptions : public SharedOptions
{

    CharString      queryFile;
    bool            revComp     = true;

    CharString      output;

    unsigned short  queryPart = 0;

//     bool            semiGlobal;

    bool            doubleIndexing = true;

    uint8_t         seedLength  = 0;
    uint8_t         maxSeedDist = 1;
    bool            hammingOnly = false;

    int8_t          seedGravity     = 0;
    uint8_t         seedOffset      = 0;
    uint8_t         minSeedLength   = 0;
    uint8_t         minSeedScore    = 0;
//     unsigned int    minSeedEVal     = 0;
//     double          minSeedBitS     = -1;

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

    unsigned        threads     = 1;
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
    CharString      segFile = "";

    bool            indexIsFM = false;

    LambdaIndexerOptions()
        : SharedOptions()
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
    setVersion(parser, "0.3");
    setDate(parser, "September 2014");

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

    addOption(parser, seqan::ArgParseOption("it",
                                            "index-type",
                                            "database index is in this format "
                                            "(auto means \"try sa first then "
                                            "fm\").",
                                            seqan::ArgParseArgument::STRING,
                                            "STR"));
    setValidValues(parser, "index-type", "auto sa fm");
    setDefaultValue(parser, "index-type", "auto");

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
                                            "STR"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

#ifdef _OPENMP
    addOption(parser, seqan::ArgParseOption("t",
                                            "threads",
                                            "number of threads to run "
                                            "concurrently.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", omp_get_max_threads());
#else
    addOption(parser, seqan::ArgParseOption("t",
                                            "threads",
                                            "LAMBDA BUILT WITHOUT OPENMP; "
                                            "setting this option has no effect.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", 1);
#endif

    addOption(parser, seqan::ArgParseOption("qp",
                                            "query-partitions",
                                            "Divide the query into qp number "
                                            "of blocks before processing; "
                                            "should be a multiple of the number "
                                            "of threads; doubling qp halves "
                                            "the memory used for cached hits, "
                                            "see below.",
                                            seqan::ArgParseArgument::INTEGER));
#ifdef _OPENMP
    setDefaultValue(parser, "query-partitions", 2 * omp_get_max_threads());
#else
     setDefaultValue(parser, "query-partitions", 2);
#endif

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
                                            "allow only mismatches in seeds.",
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

    addOption(parser, seqan::ArgParseOption("si",
                                            "seed-index",
                                            "use double indexing or not.",
                                            seqan::ArgParseArgument::STRING));
    setValidValues(parser, "seed-index", "on off");
    setDefaultValue(parser, "seed-index", "on");

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
                                            "Maximum E-Value for Results.",
                                            seqan::ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e-value", "0.1");


    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Memory requirements VS speed");
    addText(parser, "Lambda has three main points of memory consumption:");
    addText(parser, "1) the cache of the hits. This depends on the amount of "
                    "hits and is difficult to estimate. doubling the -qp value "
                    "will reduce this memory by a half at a modest penalty "
                    "to run-time.");
    addText(parser, "2) the database index. This depends on the size n of the "
                    "database and the type of index. The SA index is 6*n in "
                    "size, while the FM index is < 2*n in size. Using the FM "
                    "index increases total running time by about ~ 10%, though."
                    " Choose this when creating the index, default is SA.");
    addText(parser, "3) the sequence strings. If you define "
                    "LAMBDA_BITCOPMRESSED_STRINGS while compiling lambda, you "
                    "can reduce the size per character from 8 to 4 bit. This "
                    "will also increase running time, by ~ 10%.");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    getOptionValue(options.dbFile, parser, "database");

    CharString buffer;
    getOptionValue(buffer, parser, "index-type");
    if (buffer == "auto")
        options.dbIndexType = -1;
    else if (buffer == "sa")
        options.dbIndexType = 0;
    else // if fm
        options.dbIndexType = 1;

    getOptionValue(buffer, parser, "program");
    if (buffer == "blastn")
    {
        options.blastProg   = BlastFormatProgram::BLASTN;
        options.seedLength  = 20;
        options.alphReduction = 0;
    }
    else if (buffer == "blastp")
    {
        options.blastProg = BlastFormatProgram::BLASTP;
        options.seedLength  = 10;
    }
    else if (buffer == "blastx")
    {
        options.blastProg = BlastFormatProgram::BLASTX;
        options.seedLength  = 10;
    }
    else if (buffer == "tblastn")
    {
        options.blastProg = BlastFormatProgram::TBLASTN;
        options.seedLength  = 10;
    }
    else if (buffer == "tblastx")
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

    CharString doubleIndex;
    getOptionValue(doubleIndex, parser, "seed-index");
    options.doubleIndexing = (doubleIndex == "on");

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

#ifdef _OPENMP
    getOptionValue(foo, parser, "threads");
    options.threads = foo;
    omp_set_num_threads(options.threads);
#else
    options.threads = 1;
#endif

    getOptionValue(foo, parser, "query-partitions");
    options.queryPart = foo;
    if ((options.queryPart % options.threads) != 0)
        std::cout << "-qp not a multiple of -t; expect suboptimal performance.\n";

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
    setVersion(parser, "0.3");
    setDate(parser, "September 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-i DATABASE.fasta\\fP [\\-o INDEXFILE.sa\\fP]");
//     addDescription(parser, "This is the application skelleton and you should modify this string.");


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
//     addOption(parser, seqan::ArgParseOption("o",
//                                             "output",
//                                             "Index of database sequences",
//                                             seqan::ArgParseArgument::OUTPUTFILE,
//                                             "OUT"));
//     setValidValues(parser, "output", "sa fm");

    addOption(parser, seqan::ArgParseOption("it",
                                            "index-type",
                                            "suffix array or full-text minute "
                                            "space",
                                            seqan::ArgParseArgument::STRING,
                                            "type"));
    setValidValues(parser, "index-type", "sa fm");
    setDefaultValue(parser, "index-type", "sa");



    addSection(parser, "Program Options");
    addOption(parser, seqan::ArgParseOption("p",
                                            "program",
                                            "Blast Operation Mode.",
                                            seqan::ArgParseArgument::STRING,
                                            "program"));
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

    addTextSection(parser, "Memory requirements");
    addText(parser, "The FM index is < 2* text size, but ~ 10% slower than "
                    "the SA index which is 6x text size. During construction "
                    "a full suffix array always needs to be built in memory.");
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

    CharString buffer;
    getOptionValue(buffer, parser, "program");
    if (buffer == "blastn")
    {
        options.blastProg   = BlastFormatProgram::BLASTN;
        options.alphReduction = 0;

    }
    else if (buffer == "blastp")
    {
        options.blastProg = BlastFormatProgram::BLASTP;

    }
    else if (buffer == "blastx")
    {
        options.blastProg = BlastFormatProgram::BLASTX;

    }
    else if (buffer == "tblastn")
    {
        options.blastProg = BlastFormatProgram::TBLASTN;

    }
    else if (buffer == "tblastx")
    {
        options.blastProg = BlastFormatProgram::TBLASTX;

    }
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    //verifyFileFormat(); TODO

    CharString dbIndexType;
    getOptionValue(dbIndexType, parser, "index-type");

    if (dbIndexType == "fm")
        options.dbIndexType = 1;
    else
        options.dbIndexType = 0;

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


constexpr const char *
_alphName(AminoAcid const & /**/)
{
    return "AminoAcid";
}

constexpr const char *
_alphName(ReducedAminoAcid<Murphy10> const & /**/)
{
    return "Murphy10";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<8>> const & /**/)
{
    return "Lambda08";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<10>> const & /**/)
{
    return "Lambda10";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<12>> const & /**/)
{
    return "Lambda12";
}

constexpr const char *
_alphName(Dna const & /**/)
{
    return "Dna4";
}

constexpr const char *
_alphName(Dna5 const & /**/)
{
    return "Dna5";
}

template <typename TLH>
inline void
printOptions(LambdaOptions const & options)
{
    using TGH = typename TLH::TGlobalHolder;
    using TFormat = typename TGH::TFormat;
    auto constexpr p = getProgramType(TFormat());
    bool isFM = std::is_same<typename TGH::TIndexSpec,FMIndex<>>::value;


    std::string bandStr;
    switch(options.band)
    {
        case -3: bandStr = "2 * log(queryLength) + 1"; break;
        case -2: bandStr = "2 * sqrt(queryLength) + 1"; break;
        case -1: bandStr = "no band"; break;
        default: bandStr = std::to_string(2 * options.band + 1); break;
    }

    std::cout << "OPTIONS\n"
              << " I/O\n"
              << "  query file:               " << options.queryFile << "\n"
              << "  db file:                  " << options.dbFile << "\n"
              << "  db index type:            " << (isFM
                                                    ? "FM-Index\n"
                                                    : "SA-Index\n")
              << "  output file:              " << options.output << "\n"
              << " PROGRAM\n"
              << "  blast mode:               " << _programTagToString(TFormat())
              << "\n"
              << "  threads:                  " << uint(options.threads) << "\n"
              << "  query partitions:         " << uint(options.queryPart) << "\n"
              << " TRANSLATION AND ALPHABETS\n"
              << "  genetic code:             " << uint(options.geneticCode) << "\n"
              << "  original alphabet (query):" << _alphName(OrigQryAlph<p>())
              << "\n"
              << "  original alphabet (subj): " << _alphName(OrigSubjAlph<p>())
              << "\n"
              << "  translated alphabet:      " << _alphName(TransAlph<p>())
              << "\n"
              << "  reduced alphabet:         " << _alphName(typename TGH::TRedAlph())
              << "\n"
              << " SEEDING\n"
              << "  seed length:              " << uint(options.seedLength) << "\n"
              << "  seed offset:              " << uint(options.seedOffset) << "\n"
              << "  seed delta:               " << uint(options.maxSeedDist) << "\n"
              << "  seeds ungapped:           " << uint(options.hammingOnly) << "\n"
              << "  seed gravity:             " << uint(options.seedGravity) << "\n"
              << "  min seed length:          " << uint(options.minSeedLength) << "\n"
              << "  min seed score:           " << uint(options.minSeedScore) << "\n"
              << "  double indexing:          " << options.doubleIndexing << "\n"
//               << "  min seed e-value:         " << uint(options.minSeedEVal) << "\n"
//               << "MinSeedBitS:   " << options.minSeedBitS << "\n"
              << " EXTENSION\n"
              << "  x-drop:                   " << options.xDropOff << "\n"
              << "  band:                     " << bandStr << "\n"
              << "  maximum e-value:          " << options.eCutOff << "\n"
              << " MISC\n"
              << "  stdout is terminal:       " << options.isTerm << "\n"
              << "  terminal width:           " << options.terminalCols << "\n"
              << "  bit-compressed strings:   "
    #if defined LAMBDA_BITCOPMRESSED_STRINGS
              << "on\n"
    #else
              << "off\n"
    #endif
              << "\n";
}


#endif // header guard
