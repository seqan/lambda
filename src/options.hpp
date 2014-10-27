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

template<typename TSpec2, typename TSpec3, typename TFunctor>
struct SAValue<StringSet<ModifiedString<String<AminoAcid, TSpec2>, TFunctor>, TSpec3> >
{
    typedef Pair<uint32_t, uint16_t, Pack> Type;
};

template<typename TSpec2, typename TSpec3>
struct SAValue<StringSet<String<AminoAcid, TSpec2>, TSpec3> >
{
    typedef Pair<uint32_t, uint16_t, Pack> Type;
};

// Dna Sequences might be longer TODO REALLY?
template<typename TSpec1, typename TSpec2>
struct SAValue<StringSet<String<Dna5, TSpec1>, TSpec2> >
{
    typedef Pair<uint32_t, uint16_t, Pack> Type;
};

}

using namespace seqan;

#if defined LAMBDA_BITCOPMRESSED_STRINGS
using PackSpec = Packed<>;
#else
using PackSpec = Alloc<>;
#endif

template <typename TAlph>
using TCDStringSet = StringSet<String<TAlph, PackSpec>, Owner<ConcatDirect<> > >;


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
    int verbosity = 1;

    std::string dbFile;

    int      dbIndexType = 0;
    // for indexer, the file format of database sequences
    // for main app, the file format of query sequences
    // 0 -- fasta, 1 -- fastq
    int      fileFormat = 0;

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

    std::string     queryFile;
    bool            revComp     = true;

    std::string     output;

    unsigned        queryPart = 0;

//     bool            semiGlobal;

    bool            doubleIndexing = true;

    unsigned        seedLength  = 0;
    unsigned        maxSeedDist = 1;
    bool            hammingOnly = true;

    int             seedGravity     = 0;
    unsigned        seedOffset      = 0;
    unsigned        minSeedLength   = 0;
    int             minSeedScore    = 0;
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
    std::string     segFile = "";

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
    setVersion(parser, "0.4");
    setDate(parser, "October 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP "
                         "\\fI-d DATABASE.fasta\\fP "
                         "[\\fI-o output.m8\\fP]");
    addDescription(parser, "Lambda is a local aligner optimized for many query "
    "sequences and searches in protein space. It is faster than BLAST and many "
    "other comparable tools.");

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("q", "query",
        "Query sequences (fasta).",
        seqan::ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "q");
    setValidValues(parser, "query", "fasta fa fna faa fas fastq fq");

    addOption(parser, ArgParseOption("d", "database",
        "Database sequences (fasta), with precomputed index (.sa).",
        seqan::ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "d");
    setValidValues(parser, "database", "fasta fa fna faa fas");

    addOption(parser, ArgParseOption("di", "db-index-type",
        "database index is in this format.",
//         "(auto means \"try sa first then fm\").",
        seqan::ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "db-index-type", "sa fm");
    setDefaultValue(parser, "db-index-type", "sa");

    addSection(parser, "Output Options");
    addOption(parser, ArgParseOption("o", "output",
        "File to hold reports on hits (.m8 is blastall -m8 et cetera)",
        seqan::ArgParseArgument::OUTPUTFILE,
        "OUT"));
    setValidValues(parser, "output", "m0 m8 m9");
    setDefaultValue(parser, "output", "output.m8");

    addOption(parser, ArgParseOption("v", "verbosity",
        "The amount of terminal output printed; 0 [only errors]; 1 [default]; 2 "
        "[+run-time, options and statistics].",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbosity", "1");

    addSection(parser, "Alphabets and Translation");
    addOption(parser, ArgParseOption("p", "program",
        "Blast Operation Mode.",
        seqan::ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

//     addOption(parser, ArgParseOption("qa", "query-alphabet",
//         "original alphabet of the query sequences",
//         seqan::ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "query-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "query-alphabet", "dna5");
// 
//     addOption(parser, ArgParseOption("da", "db-alphabet",
//         "original alphabet of the subject sequences",
//         seqan::ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "db-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "db-alphabet", "aminoacid");
// 
//     addOption(parser, ArgParseOption("sa", "seeding-alphabet",
//         "alphabet to use during seeding (reduction possible)",
//         seqan::ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "seeding-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "seeding-alphabet", "murphy10");

    addOption(parser, ArgParseOption("g", "genetic-code",
        "The translation table to use for nucl -> amino acid translation"
        "(not for BlastN, BlastP). See "
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
        " for ids (default is generic)",
        seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, ArgParseOption("ar", "alphabet-reduction",
        "Alphabet Reduction for seeding phase (ignored for BLASTN).",
        seqan::ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "alphabet-reduction", "none murphy10");
    setDefaultValue(parser, "alphabet-reduction", "murphy10");

    addSection(parser, "General Options");

#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads",
        "number of threads to run concurrently.",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", omp_get_max_threads());
#else
    addOption(parser, ArgParseOption("t", "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", 1);
#endif

    addOption(parser, ArgParseOption("qi", "query-index-type",
        "setting this to \"none\" deactivates double-indexing.",
        seqan::ArgParseArgument::STRING));
    setValidValues(parser, "query-index-type", "radix none");
    setDefaultValue(parser, "query-index-type", "radix");

    addOption(parser, ArgParseOption("qp", "query-partitions",
        "Divide the query into qp number of blocks before processing; should be"
        " a multiple of the number of threads, defaults to one per thread. "
        "Only used with double-indexing; strong influence on memory, see below.",
        seqan::ArgParseArgument::INTEGER));
#ifdef _OPENMP
    setDefaultValue(parser, "query-partitions", omp_get_max_threads());
#else
    setDefaultValue(parser, "query-partitions", 1);
#endif

    addSection(parser, "Seeding / Filtration");
//     addOption(parser, ArgParseOption("su",
//                                             "ungapped-seeds",
//                                             "allow only mismatches in seeds.",
//                                             seqan::ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "ungapped-seeds", "1");

    addOption(parser, ArgParseOption("sl", "seed-length",
        "Length of the seeds (default = 20 for BLASTN).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-length", "10");

    addOption(parser, ArgParseOption("so", "seed-offset",
        "Offset for seeding (if unset = seed-length; non-overlapping).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-offset", "10");

    addOption(parser, ArgParseOption("sd", "seed-delta",
        "maximum seed distance.",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-delta", "1");

    addOption(parser, ArgParseOption("sg", "seed-gravity",
        "Seeds closer than this are merged into region (if unset = "
        "seed-length).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-gravity", "10");

    addOption(parser, ArgParseOption("sm", "seed-min-length",
        "after postproc shorter seeds are discarded (if unset = seed-length).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-length", "10");

    addOption(parser, ArgParseOption("ss", "seed-min-score",
        "after postproc worse seeds are discarded (raw score). Only used in "
        "combination with Alphabet reduction.",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-score", "32");

//     addOption(parser, ArgParseOption("se",
//                                             "seedminevalue",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(0 -> off).",
//                                             seqan::ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "seedminevalue", "100000");

//     addOption(parser, ArgParseOption("sb",
//                                             "seedminbits",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(-1 -> off).",
//                                             seqan::ArgParseArgument::DOUBLE));
//     setDefaultValue(parser, "seedminbits", "-1");

    addSection(parser, "Scoring");

    addOption(parser, ArgParseOption("sc", "scoring-scheme",
        "'62' for Blosum62 (default); '50' for Blosum50; '0' for manual "
        "(default for BlastN)",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "scoring-scheme", "62");

    addOption(parser, ArgParseOption("ge", "score-gap",
        "Score per gap character (default = -2 for BLASTN).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap", "-1");

    addOption(parser, ArgParseOption("go", "score-gap-open",
        "Additional cost for opening gap (default = -5 for BLASTN).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap-open", "-11");

    addOption(parser, ArgParseOption("ma", "score-match",
        "Match score (BLASTN or manual scoring)",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-match", "2");

    addOption(parser, ArgParseOption("mi", "score-mismatch",
        "Mismatch score (BLASTN or manual scoring)",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-mismatch", "-3");

    addSection(parser, "Extension");

    addOption(parser, ArgParseOption("x", "x-drop",
        "Stop Banded extension if score x below the maximum seen (-1 means no "
        "xdrop).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "x-drop", "30");
    setMinValue(parser, "x-drop", "-1");

    addOption(parser, ArgParseOption("b", "band",
        "Size of the DP-band used in extension (-3 means log2 of query length; "
        "-2 means sqrt of query length; -1 means full dp; n means band of size "
        "2n+1)",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "band", "-3");
    setMinValue(parser, "band", "-3");

    addOption(parser, ArgParseOption("e", "e-value",
        "Maximum E-Value for Results.",
        seqan::ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e-value", "0.1");


    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Speed VS sensitivity");
    addText(parser, "Playing with the seeding and alphabet parameters has high "
                    "influence on both speed and sensitivity. We recommend the "
                    "following alternative profiles for protein searches:");
    addText(parser, "fast: \033[1m-a 0 -sl 8 -ss 26 -sd 0 -so 4\033[0m");
    addText(parser, "sensitive: \033[1m-so 5\033[0m");

    addTextSection(parser, "Speed VS memory requirements");
    addText(parser, "Lambda has three main points of memory consumption:");
    addText(parser, "1) the cache of putative hits. Since version 0.4 this is "
                    "quite small, but grows linearly with input data when "
                    "double-indexing is used; doubling \033[1m-qp\033[0m will "
                    "reduce this memory by a half at a modest penalty "
                    "to run-time; however using single-indexing might be "
                    "better if you are memory-constrained.");
    addText(parser, "2) the database index. This depends on the size n of the "
                    "database and the type of index. The SA index is 6*n in "
                    "size, while the FM index is < 2*n in size. Using the FM "
                    "index increases total running time by about ~ 10%, though."
                    " Choose this when creating the index, default is SA.");
    addText(parser, "3) the sequence strings. If you define "
                    "LAMBDA_BITCOPMRESSED_STRINGS while compiling lambda, you "
                    "can reduce the size per character from 8 to 5/4 bit. This "
                    "will also increase running time, by ~ 10%.");
    addText(parser, "As a rule of thumb, if you run out of memory, add "
                    "\033[1m-di fm\033[0m to the indexer call and \033[1m-di fm"
                    " -qi none\033[0m to the lambda call.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    if (endsWith(options.queryFile, ".fastq") ||
        endsWith(options.queryFile, ".fq"))
        options.fileFormat = 1;
    else
        options.fileFormat = 0;

    getOptionValue(options.dbFile, parser, "database");

    getOptionValue(options.output, parser, "output");

    int buf = 0;
    getOptionValue(buf, parser, "verbosity");
    switch(buf)
    {
        case 0: options.verbosity = 0; break;
        case 2: options.verbosity = 2; break;
        default: options.verbosity = 1; break;
    }

    std::string buffer;
    getOptionValue(buffer, parser, "db-index-type");
    if (buffer == "auto")
        options.dbIndexType = -1;
    else if (buffer == "sa")
        options.dbIndexType = 0;
    else // if fm
        options.dbIndexType = 1;

    getOptionValue(buffer, parser, "program");
    if (buffer == "blastn")
        options.blastProg = BlastFormatProgram::BLASTN;
    else if (buffer == "blastp")
        options.blastProg = BlastFormatProgram::BLASTP;
    else if (buffer == "blastx")
        options.blastProg = BlastFormatProgram::BLASTX;
    else if (buffer == "tblastn")
        options.blastProg = BlastFormatProgram::TBLASTN;
    else if (buffer == "tblastx")
        options.blastProg = BlastFormatProgram::TBLASTX;
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    getOptionValue(buffer, parser, "alphabet-reduction");
    if ((buffer == "murphy10") &&
        (options.blastProg != BlastFormatProgram::BLASTN))
        options.alphReduction = 2;
    else
        options.alphReduction = 0;

    getOptionValue(options.seedLength, parser, "seed-length");
    if ((!isSet(parser, "seed-length")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.seedLength = 15;

    if (isSet(parser, "seed-offset"))
        getOptionValue(options.seedOffset, parser, "seed-offset");
    else
        options.seedOffset = options.seedLength;

    if (isSet(parser, "seed-gravity"))
        getOptionValue(options.seedGravity, parser, "seed-gravity");
    else
        options.seedGravity = options.seedLength;

    if (isSet(parser, "seed-min-length"))
        getOptionValue(options.minSeedLength, parser, "seed-min-length");
    else
        options.minSeedLength = options.seedLength;

    getOptionValue(options.maxSeedDist, parser, "seed-delta");

    getOptionValue(options.minSeedScore, parser, "seed-min-score");
    if (options.alphReduction == 0)
        options.minSeedScore = 0;

    getOptionValue(buffer, parser, "query-index-type");
    options.doubleIndexing = (buffer == "radix");

    getOptionValue(options.eCutOff, parser, "e-value");

    getOptionValue(options.xDropOff, parser, "x-drop");
    if ((!isSet(parser, "x-drop")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.xDropOff = 16;

    getOptionValue(options.band, parser, "band");

#ifdef _OPENMP
    getOptionValue(options.threads, parser, "threads");
    omp_set_num_threads(options.threads);
#else
    options.threads = 1;
#endif

    if (options.doubleIndexing)
    {
        getOptionValue(options.queryPart, parser, "query-partitions");
        if ((options.queryPart % options.threads) != 0)
            std::cout << "-qp not a multiple of -t; expect suboptimal performance.\n";
    } else
    {
        options.queryPart = options.threads;
    }

    getOptionValue(options.scoringMethod, parser, "scoring-scheme");
    if (options.blastProg == BlastFormatProgram::BLASTN)
        options.scoringMethod = 0;
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
    if ((!isSet(parser, "score-gap")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.gapExtend = -2;

    getOptionValue(options.gapOpen, parser, "score-gap-open");
    if ((!isSet(parser, "score-gap-open")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.gapOpen = -5;

    getOptionValue(buf, parser, "genetic-code");
    switch (buf)
    {
        case 1: case 2: case 3: case 4: case 5: case 6:
        case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        case 21: case 22: case 23: case 24 : case 25:
            options.geneticCode = static_cast<GeneticCodeSpec>(buf);
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
    setVersion(parser, "0.4");
    setDate(parser, "October 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-i DATABASE.fasta\\fP");
//     addDescription(parser, "This is the application skelleton and you should modify this string.");


    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("d", "database",
        "Database sequences (fasta).",
        seqan::ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "database");
    setValidValues(parser, "database", "fasta fa fna faa");

    addOption(parser, ArgParseOption("s",
        "segfile",
        "SEG intervals for database"
        "(optional).",
        seqan::ArgParseArgument::INPUTFILE));

    setValidValues(parser, "segfile", "seg");


    addSection(parser, "Output Options");
//     addOption(parser, ArgParseOption("o",
//                                             "output",
//                                             "Index of database sequences",
//                                             seqan::ArgParseArgument::OUTPUTFILE,
//                                             "OUT"));
//     setValidValues(parser, "output", "sa fm");

    addOption(parser, ArgParseOption("di", "db-index-type",
        "suffix array or full-text minute space.",
        seqan::ArgParseArgument::STRING,
        "type"));
    setValidValues(parser, "db-index-type", "sa fm");
    setDefaultValue(parser, "db-index-type", "sa");

    addSection(parser, "Alphabets and Translation");
    addOption(parser, ArgParseOption("p", "program",
        "Blast Operation Mode.",
        seqan::ArgParseArgument::STRING,
        "program"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");
    addOption(parser, ArgParseOption("g", "genetic-code",
        "The translation table to use (not for BlastN, BlastP). See "
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
        " for ids (default is generic).",
        seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, ArgParseOption("ar", "alphabet-reduction",
        "Alphabet Reduction for seeding phase (ignored for BLASTN).",
        seqan::ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "alphabet-reduction", "none murphy10");
    setDefaultValue(parser, "alphabet-reduction", "murphy10");

    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Memory requirements");
    addText(parser, "See \033[1mlambda --help\033[0m on details concerning the "
                    "choice of db-index-type. Note also that during "
                    "construction a full suffix array always needs to be built "
                    "in memory, so if you fail to produce an fm-index, you "
                    "could build it on another machine and still succeed "
                    "at running lambda.");
    addTextSection(parser, "Remarks");
    addText(parser, "Note that the indeces created are binary and not "
                    "necessarily platform independant.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.

    seqan::getOptionValue(options.dbFile, parser, "database");

    std::string buffer;
    getOptionValue(buffer, parser, "program");
    if (buffer == "blastn")
        options.blastProg = BlastFormatProgram::BLASTN;
    else if (buffer == "blastp")
        options.blastProg = BlastFormatProgram::BLASTP;
    else if (buffer == "blastx")
        options.blastProg = BlastFormatProgram::BLASTX;
    else if (buffer == "tblastn")
        options.blastProg = BlastFormatProgram::TBLASTN;
    else if (buffer == "tblastx")
        options.blastProg = BlastFormatProgram::TBLASTX;
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    getOptionValue(buffer, parser, "alphabet-reduction");
    if ((buffer == "murphy10") &&
        (options.blastProg != BlastFormatProgram::BLASTN))
        options.alphReduction = 2;
    else
        options.alphReduction = 0;

    getOptionValue(buffer, parser, "db-index-type");

    if (buffer == "fm")
        options.dbIndexType = 1;
    else
        options.dbIndexType = 0;

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
    return "aminoacid";
}

constexpr const char *
_alphName(ReducedAminoAcid<Murphy10> const & /**/)
{
    return "murphy10";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<8>> const & /**/)
{
    return "lambda08";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<10>> const & /**/)
{
    return "lambda10";
}

constexpr const char *
_alphName(ReducedAminoAcid<ClusterReduction<12>> const & /**/)
{
    return "lambda12";
}

constexpr const char *
_alphName(Dna const & /**/)
{
    return "dna4";
}

constexpr const char *
_alphName(Dna5 const & /**/)
{
    return "dna5";
}

template <typename TLH>
inline void
printOptions(LambdaOptions const & options)
{
    using TGH = typename TLH::TGlobalHolder;
    using TFormat = typename TGH::TFormat;
    auto constexpr p = getProgramType(TFormat());

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
              << "  db index type:            " << (TGH::indexIsFM
                                                    ? "FM-Index\n"
                                                    : "SA-Index\n")
              << "  output file:              " << options.output << "\n"
              << " GENERAL\n"
              << "  double indexing:          " << options.doubleIndexing << "\n"
              << "  threads:                  " << uint(options.threads) << "\n"
              << "  query partitions:         " << (options.doubleIndexing
                                                    ? std::to_string(options.queryPart)
                                                    : std::string("n/a")) << "\n"
              << " TRANSLATION AND ALPHABETS\n"
              << "  genetic code:             " << uint(options.geneticCode) << "\n"
              << "  blast mode:               " << _programTagToString(TFormat())
              << "\n"
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
//               << "  min seed e-value:         " << uint(options.minSeedEVal) << "\n"
//               << "MinSeedBitS:   " << options.minSeedBitS << "\n"
              << " SCORING\n"
              << "  scoring scheme:           " << options.scoringMethod << "\n"
              << "  score-match:              " << (options.scoringMethod
                                                    ? std::string("n/a")
                                                    : std::to_string(options.match)) << "\n"
              << "  score-mismatch:           " << (options.scoringMethod
                                                    ? std::string("n/a")
                                                    : std::to_string(options.misMatch)) << "\n"
              << "  score-gap:                " << options.gapExtend << "\n"
              << "  score-gap-open:           " << options.gapOpen << "\n"
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
