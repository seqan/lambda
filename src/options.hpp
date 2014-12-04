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

#include "index_sa_sort.h"

#define LAMBDA_VERSION "0.4.5"

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

// Dna Sequences might be longer (chromosomes, genomes)
template<typename TSpec1, typename TSpec2>
struct SAValue<StringSet<String<Dna5, TSpec1>, TSpec2> >
{
    typedef Pair<uint32_t, uint32_t, Pack> Type;
};

}

using namespace seqan;

// Index Specs
// template <typename TActualSpec = void>
// struct SAqsSpec {} ;

template <typename TSpec = void>
using TFMIndex = FMIndex<TSpec>;

namespace SEQAN_NAMESPACE_MAIN
{
template <typename TText, typename TSpec>
struct Fibre<Index<TText, IndexSa<SaAdvancedSort<TSpec> > >, FibreTempSA>
{
    typedef Index<TText, IndexSa<SaAdvancedSort<TSpec> > >        TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template <typename TText, typename TSpec>
struct DefaultIndexCreator<Index<TText, IndexSa<SaAdvancedSort<TSpec>> >, FibreSA>
{
    typedef SaAdvancedSort<TSpec> Type;
//     static std::function<void(void)> progressCallback;
};

// template <typename TText, typename TSpec>
// static std::function<void(void)>
// DefaultIndexCreator<Index<TText, IndexSa<SaAdvancedSort<TSpec>> >, FibreSA>::progressCallback = [] () {};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >, FibreTempSA>
{
    typedef Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >        TIndex_;
    typedef typename SAValue<TIndex_>::Type                         TSAValue_;

    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type> Type;
};

template < typename TText, typename TSpec, typename TConfig>
struct DefaultIndexCreator<Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >, FibreSA>
{
    typedef SaAdvancedSort<TSpec> Type;
//     std::function<void(void)> progressCallback;
};

// template <typename TText, typename TSpec>
// static std::function<void(void)>
// DefaultIndexCreator<Index<TText, FMIndex<SaAdvancedSort<TSpec>, TConfig> >, FibreSA>::progressCallback = [] () {};

}

#if defined LAMBDA_BITCOPMRESSED_STRINGS
using PackSpec = Packed<>;
#else
using PackSpec = Alloc<>;
#endif

template <typename TAlph>
using TCDStringSet = StringSet<String<TAlph, PackSpec>, Owner<ConcatDirect<> > >;

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

    unsigned        threads     = 1;

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
    int             idCutOff    = 0;
    unsigned long   maxMatches  = 500;

    bool            filterPutativeDuplicates = true;
    bool            filterPutativeAbundant = true;

    int             preScoring = 0; // 0 = off, 1 = seed, 2 = region (
    double          preScoringThresh    = 0.0;

    LambdaOptions() :
        SharedOptions()
    {
    }
};

struct LambdaIndexerOptions : public SharedOptions
{
    std::string     segFile = "";
    std::string     algo = "";

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

// SHARED
ArgumentParser::ParseResult
parseCommandLineShared(SharedOptions & options, ArgumentParser & parser);

ArgumentParser::ParseResult
parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("lambda");
    // Set short description, version, and date.
    setShortDescription(parser, "the Local Aligner for Massive Biological "
    "DatA");
    setVersion(parser, LAMBDA_VERSION);
    setDate(parser, __DATE__);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP "
                         "\\fI-d DATABASE.fasta\\fP "
                         "[\\fI-o output.m8\\fP]");
    addDescription(parser, "Lambda is a local aligner optimized for many query "
    "sequences and searches in protein space. It is compatible to BLAST, but "
    "much faster than BLAST and many other comparable tools.");

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("q", "query",
        "Query sequences.",
        ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "q");
    setValidValues(parser, "query", "fasta fa fna faa fas fastq fq");

    addOption(parser, ArgParseOption("d", "database",
        "Database sequences (fasta), with precomputed index (.sa or .fm).",
        ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "d");
    setValidValues(parser, "database", "fasta fa fna faa fas");

    addOption(parser, ArgParseOption("di", "db-index-type",
        "database index is in this format.",
//         "(auto means \"try sa first then fm\").",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "db-index-type", "sa fm");
    setDefaultValue(parser, "db-index-type", "fm");

    addSection(parser, "Output Options");
    addOption(parser, ArgParseOption("o", "output",
        "File to hold reports on hits (.m8 is blastall -m8 et cetera)",
        ArgParseArgument::OUTPUTFILE,
        "OUT"));
    setValidValues(parser, "output", "m0 m8 m9");
    setDefaultValue(parser, "output", "output.m8");

    addOption(parser, ArgParseOption("id", "percent-identity",
        "Output only matches above this threshold (checked before e-value "
        "check).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "percent-identity", "0");
    setMinValue(parser, "percent-identity", "0");
    setMaxValue(parser, "percent-identity", "100");

    addOption(parser, ArgParseOption("e", "e-value",
        "Output only matches that score below this threshold.",
        ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e-value", "0.1");
    setMinValue(parser, "e-value", "0");

    addOption(parser, ArgParseOption("nm", "num-matches",
        "Print at most this number of matches per query.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "num-matches", "500");
    setMinValue(parser, "num-matches", "1");

    addOption(parser, ArgParseOption("v", "verbosity",
        "The amount of terminal output printed; 0 [only errors]; 1 [default]; 2 "
        "[+run-time, options and statistics].",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbosity", "1");
    setMinValue(parser, "verbosity", "0");
    setMaxValue(parser, "verbosity", "2");

    addSection(parser, "General Options");
#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads",
        "number of threads to run concurrently.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", omp_get_max_threads());
#else
    addOption(parser, ArgParseOption("t", "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", 1);
#endif

    addOption(parser, ArgParseOption("qi", "query-index-type",
        "controls double-indexing.",
        ArgParseArgument::STRING));
    setValidValues(parser, "query-index-type", "radix none");
    setDefaultValue(parser, "query-index-type", "none");

    addOption(parser, ArgParseOption("qp", "query-partitions",
        "Divide the query into qp number of blocks before processing; should be"
        " a multiple of the number of threads, defaults to one per thread. "
        "Only used with double-indexing; strong influence on memory, see below.",
        ArgParseArgument::INTEGER));
#ifdef _OPENMP
    setDefaultValue(parser, "query-partitions", omp_get_max_threads());
#else
    setDefaultValue(parser, "query-partitions", 1);
#endif
    hideOption(parser, "query-partitions"); // HIDDEN

    addSection(parser, "Alphabets and Translation");
    addOption(parser, ArgParseOption("p", "program",
        "Blast Operation Mode.",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

//     addOption(parser, ArgParseOption("qa", "query-alphabet",
//         "original alphabet of the query sequences",
//         ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "query-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "query-alphabet", "dna5");
// 
//     addOption(parser, ArgParseOption("da", "db-alphabet",
//         "original alphabet of the subject sequences",
//         ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "db-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "db-alphabet", "aminoacid");
// 
//     addOption(parser, ArgParseOption("sa", "seeding-alphabet",
//         "alphabet to use during seeding (reduction possible)",
//         ArgParseArgument::STRING,
//         "STR"));
//     setValidValues(parser, "seeding-alphabet", "dna5 aminoacid");
//     setDefaultValue(parser, "seeding-alphabet", "murphy10");

    addOption(parser, ArgParseOption("g", "genetic-code",
        "The translation table to use for nucl -> amino acid translation"
        "(not for BlastN, BlastP). See "
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
        " for ids (default is generic). Six frames are generated.",
        ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, ArgParseOption("ar", "alphabet-reduction",
        "Alphabet Reduction for seeding phase (ignored for BLASTN).",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "alphabet-reduction", "none murphy10");
    setDefaultValue(parser, "alphabet-reduction", "murphy10");

    addSection(parser, "Seeding / Filtration");
//     addOption(parser, ArgParseOption("su",
//                                             "ungapped-seeds",
//                                             "allow only mismatches in seeds.",
//                                             ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "ungapped-seeds", "1");

    addOption(parser, ArgParseOption("sl", "seed-length",
        "Length of the seeds (default = 14 for BLASTN).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-length", "10");

    addOption(parser, ArgParseOption("so", "seed-offset",
        "Offset for seeding (if unset = seed-length, non-overlapping; "
        "default = 5 for BLASTN).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-offset", "10");

    addOption(parser, ArgParseOption("sd", "seed-delta",
        "maximum seed distance.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-delta", "1");

    addOption(parser, ArgParseOption("sg", "seed-gravity",
        "Seeds closer than this are merged into region (if unset = "
        "seed-length).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-gravity", "10");
    hideOption(parser, "seed-gravity"); // HIDDEN

    addOption(parser, ArgParseOption("sm", "seed-min-length",
        "after postproc shorter seeds are discarded (if unset = seed-length).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-length", "10");
    hideOption(parser, "seed-min-length"); // HIDDEN

    addSection(parser, "Miscellaneous Heuristics");

    addOption(parser, ArgParseOption("ps", "pre-scoring",
        "evaluate score of a region NUM times the size of the seed "
        "before extension (0 -> no pre-scoring, 1 -> evaluate seed, n-> area "
        "around seed, as well; default = 0 when no alphabet reduction is used).",
        ArgParseArgument::INTEGER));
    setMinValue(parser, "pre-scoring", "0");
    setDefaultValue(parser, "pre-scoring", "2");

    addOption(parser, ArgParseOption("pt", "pre-scoring-threshold",
        "minimum average score per position in pre-scoring region.",
        ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "pre-scoring-threshold", "2");

    addOption(parser, ArgParseOption("pd", "filter-putative-duplicates",
        "filter hits that will likely duplicate a match already found.",
        ArgParseArgument::STRING));
    setValidValues(parser, "filter-putative-duplicates", "on off");
    setDefaultValue(parser, "filter-putative-duplicates", "on");

    addOption(parser, ArgParseOption("pa", "filter-putative-abundant",
        "If the maximum number of matches per query are found already, "
        "stop searching if the remaining realm looks unfeasable.",
        ArgParseArgument::STRING));
    setValidValues(parser, "filter-putative-abundant", "on off");
    setDefaultValue(parser, "filter-putative-abundant", "on");
//     addOption(parser, ArgParseOption("se",
//                                             "seedminevalue",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(0 -> off).",
//                                             ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "seedminevalue", "100000");

//     addOption(parser, ArgParseOption("sb",
//                                             "seedminbits",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(-1 -> off).",
//                                             ArgParseArgument::DOUBLE));
//     setDefaultValue(parser, "seedminbits", "-1");

    addSection(parser, "Scoring");

    addOption(parser, ArgParseOption("sc", "scoring-scheme",
        "'62' for Blosum62 (default); '50' for Blosum50; '0' for manual "
        "(default for BlastN)",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "scoring-scheme", "62");

    addOption(parser, ArgParseOption("ge", "score-gap",
        "Score per gap character (default = -2 for BLASTN).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap", "-1");

    addOption(parser, ArgParseOption("go", "score-gap-open",
        "Additional cost for opening gap (default = -5 for BLASTN).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-gap-open", "-11");

    addOption(parser, ArgParseOption("ma", "score-match",
        "Match score (BLASTN or manual scoring)",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-match", "2");

    addOption(parser, ArgParseOption("mi", "score-mismatch",
        "Mismatch score (BLASTN or manual scoring)",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "score-mismatch", "-3");

    addSection(parser, "Extension");

    addOption(parser, ArgParseOption("x", "x-drop",
        "Stop Banded extension if score x below the maximum seen (-1 means no "
        "xdrop).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "x-drop", "30");
    setMinValue(parser, "x-drop", "-1");

    addOption(parser, ArgParseOption("b", "band",
        "Size of the DP-band used in extension (-3 means log2 of query length; "
        "-2 means sqrt of query length; -1 means full dp; n means band of size "
        "2n+1)",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "band", "-3");
    setMinValue(parser, "band", "-3");

    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Speed VS sensitivity");
    addText(parser, "Tuning the seeding parameters and (de)activating alphabet "
                    "reduction has a strong "
                    "influence on both speed and sensitivity. We recommend the "
                    "following alternative profiles for protein searches:");
    addText(parser, "fast (high similarity):       \033[1m-ar 0 -sl 7 -sd 0\033[0m");
    addText(parser, "sensitive (lower similarity): \033[1m-so 5\033[0m");

    addTextSection(parser, "Speed VS memory requirements");
    addText(parser, "Lambda requires approximately the following amount of RAM:"
                    " \033[1msize(queryFile) + size(dbIDs) + 2 * size(dbSeqs)\033[0m. "
                    "If you have more RAM, use double indexing and SA:\n"
                    "\033[1m-di sa -qi radix\033[0m "
                    "which will result in an additional speed-up of up to 30% "
                    "compared to the published version (you need to run the "
                    "indexer with \033[1m-di sa \033[0m, as well). The amount "
                    "of RAM required will be: "
                    "\033[1msize(queryFile) + size(dbIDs) + 7 * size(dbSeqs) + n\033[0m "
                    "where n grows slowly but linearly with input size. "
                    "Note that size(dbSeqs) refers to the total "
                    "sequence length and does not include IDs (so it is less "
                    "than the size of the file).");
    addText(parser, "To save more RAM, you can define "
                    "LAMBDA_BITCOPMRESSED_STRINGS while compiling lambda. "
                    "This will reduce memory usage by about:"
                    " \033[1m0.3 * ( size(queryFile) + size(dbSeqs) )\033[0m,"
                    " but slow down lambda by about 10%.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Options shared by lambda and its indexer
    res = parseCommandLineShared(options, parser);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    std::string buffer;

    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    if (endsWith(options.queryFile, ".fastq") ||
        endsWith(options.queryFile, ".fq"))
        options.fileFormat = 1;
    else
        options.fileFormat = 0;

    getOptionValue(options.output, parser, "output");

    getOptionValue(options.seedLength, parser, "seed-length");
    if ((!isSet(parser, "seed-length")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.seedLength = 14;

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


    getOptionValue(buffer, parser, "query-index-type");
    options.doubleIndexing = (buffer == "radix");

    getOptionValue(options.eCutOff, parser, "e-value");
    getOptionValue(options.idCutOff, parser, "percent-identity");

    getOptionValue(options.xDropOff, parser, "x-drop");
//     if ((!isSet(parser, "x-drop")) &&
//         (options.blastProg == BlastFormatProgram::BLASTN))
//         options.xDropOff = 16;

    getOptionValue(options.band, parser, "band");

    if (options.doubleIndexing)
    {
        if (isSet(parser, "query-partitions"))
            getOptionValue(options.queryPart, parser, "query-partitions");
        else
            options.queryPart = options.threads;
        if ((options.queryPart % options.threads) != 0)
            std::cout << "-qp not a multiple of -t; expect suboptimal performance.\n";
    } else
    {
        options.queryPart = 1;
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
            return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.gapExtend, parser, "score-gap");
    if ((!isSet(parser, "score-gap")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.gapExtend = -2;

    getOptionValue(options.gapOpen, parser, "score-gap-open");
    if ((!isSet(parser, "score-gap-open")) &&
        (options.blastProg == BlastFormatProgram::BLASTN))
        options.gapOpen = -5;

    getOptionValue(buffer, parser, "filter-putative-duplicates");
    options.filterPutativeDuplicates = (buffer == "on");

    getOptionValue(buffer, parser, "filter-putative-abundant");
    options.filterPutativeAbundant = (buffer == "on");

    getOptionValue(options.preScoring, parser, "pre-scoring");
    if ((!isSet(parser, "pre-scoring")) &&
        (options.alphReduction == 0))
        options.preScoring = 0;

    getOptionValue(options.preScoringThresh, parser, "pre-scoring-threshold");
    if (options.preScoring == 0)
        options.preScoringThresh = 0;

    return ArgumentParser::PARSE_OK;
}

// INDEXER
ArgumentParser::ParseResult
parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("lambda_indexer");
    setShortDescription(parser, "the Local Aligner for Massive Biological "
    "DatA");
    // Set short description, version, and date.
    setVersion(parser, LAMBDA_VERSION);
    setDate(parser, __DATE__);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-d DATABASE.fasta\\fP");
    addDescription(parser, "Lambda is a local aligner optimized for many query "
    "sequences and searches in protein space. This is the indexer binary "
    "required to pre-process databases for use with the lambda-binary.");

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("d", "database",
        "Database sequences (fasta).",
        ArgParseArgument::INPUTFILE,
        "IN"));
    setRequired(parser, "database");
    setValidValues(parser, "database", "fasta fa fna faa");

    addOption(parser, ArgParseOption("s",
        "segfile",
        "SEG intervals for database"
        "(optional).",
        ArgParseArgument::INPUTFILE));

    setValidValues(parser, "segfile", "seg");


    addSection(parser, "Output Options");
//     addOption(parser, ArgParseOption("o",
//                                             "output",
//                                             "Index of database sequences",
//                                             ArgParseArgument::OUTPUTFILE,
//                                             "OUT"));
//     setValidValues(parser, "output", "sa fm");

    addOption(parser, ArgParseOption("di", "db-index-type",
        "suffix array or full-text minute space.",
        ArgParseArgument::STRING,
        "type"));
    setValidValues(parser, "db-index-type", "sa fm");
    setDefaultValue(parser, "db-index-type", "fm");

    addOption(parser, ArgParseOption("v", "verbosity",
        "The amount of terminal output printed; 0 [only errors]; 1 [default]; 2 "
        "[+run-time, options and statistics].",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbosity", "1");
    setMinValue(parser, "verbosity", "0");
    setMaxValue(parser, "verbosity", "2");

    addSection(parser, "Alphabets and Translation");
    addOption(parser, ArgParseOption("p", "program",
        "Blast Operation Mode.",
        ArgParseArgument::STRING,
        "program"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");
    addOption(parser, ArgParseOption("g", "genetic-code",
        "The translation table to use (not for BlastN, BlastP). See "
        "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
        " for ids (default is generic).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "genetic-code", "1");

    addOption(parser, ArgParseOption("ar", "alphabet-reduction",
        "Alphabet Reduction for seeding phase (ignored for BLASTN).",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "alphabet-reduction", "none murphy10");
    setDefaultValue(parser, "alphabet-reduction", "murphy10");

    addSection(parser, "Algorithm");
    addOption(parser, ArgParseOption("a", "algorithm",
        "Algorithm for SA construction (also used for FM; see Memory "
        " Requirements below!).",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "algorithm", "mergesort quicksortbuckets quicksort skew7ext");
    setDefaultValue(parser, "algorithm", "mergesort");
#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads",
        "number of threads to run concurrently (ignored if a == skew7ext).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", omp_get_max_threads());
#else
    addOption(parser, ArgParseOption("t", "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", 1);
#endif

    addTextSection(parser, "Memory requirements and Speed");
    addText(parser, "\033[1mmergesort [RAM}:\033[0m"
                    "\t14 * size(dbSeqs)");
    addText(parser, "\033[1mquicksort and quicksortbuckets [RAM}:\033[0m"
                    "\t7 * size(dbSeqs)");
    addText(parser, "\033[1mskew7ext [RAM]:\033[0m"
                    "\t2 * size(dbSeqs)");
    addText(parser, "\033[1mskew7ext [DISK]:\033[0m"
                    "\t30 * size(dbSeqs)");
    addText(parser, "size(dbSeqs) refers to the total "
                    "sequence length and does not include IDs (which can "
                    "account for >50% of the file size for protein databases). "
                    "The space is the maximum obseverved factor, for many "
                    "databases the factor is smaller." );
    addText(parser, "mergesort is fully parallelized and thus the fastest on "
                    "many-core-architectures. quicksort only uses up to 2 "
                    "threads, but uses significantly less memory. "
                    "quicksortbuckets uses up to n threads for most of the "
                    "operations, but is still slower than mergesort for large "
                    "files / many threads. "
                    "skew7ext is not parallelized at all, use it "
                    "only if you are very memory constrained. It also only "
                    "works if size(dbSeqs) <= 4GiB! mergesort and "
                    "quicksortbuckets provide a rough progress estimate.");
    addText(parser, "Disk space required is in TMPDIR which you can set as "
                    "an environment variable.");
    addTextSection(parser, "Remarks");
    addText(parser, "Note that the indeces created are binary and not "
                    "compatible between different CPU endiannesses.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Options shared by lambda and its indexer
    res = parseCommandLineShared(options, parser);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values
    getOptionValue(options.segFile, parser, "segfile");
    getOptionValue(options.algo, parser, "algorithm");

    return ArgumentParser::PARSE_OK;
}

// SHARED
ArgumentParser::ParseResult
parseCommandLineShared(SharedOptions & options, ArgumentParser & parser)
{
    int buf = 0;
    std::string buffer;

    getOptionValue(options.dbFile, parser, "database");

    getOptionValue(buffer, parser, "db-index-type");
    if (buffer == "sa")
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
        return ArgumentParser::PARSE_ERROR;

    getOptionValue(buffer, parser, "alphabet-reduction");
    if ((buffer == "murphy10") &&
        (options.blastProg != BlastFormatProgram::BLASTN))
        options.alphReduction = 2;
    else
        options.alphReduction = 0;

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
            return ArgumentParser::PARSE_ERROR;
    }

#ifdef _OPENMP
    getOptionValue(options.threads, parser, "threads");
    omp_set_num_threads(options.threads);
#else
    options.threads = 1;
#endif

    getOptionValue(buf, parser, "verbosity");
    switch(buf)
    {
        case 0: options.verbosity = 0; break;
        case 2: options.verbosity = 2; break;
        default: options.verbosity = 1; break;
    }

    return ArgumentParser::PARSE_OK;
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
              << " INPUT\n"
              << "  query file:               " << options.queryFile << "\n"
              << "  db file:                  " << options.dbFile << "\n"
              << "  db index type:            " << (TGH::indexIsFM
                                                    ? "FM-Index\n"
                                                    : "SA-Index\n")
              << " OUTPUT (file)\n"
              << "  output file:              " << options.output << "\n"
              << "  minimum % identity:       " << options.idCutOff << "\n"
              << "  maximum e-value:          " << options.eCutOff << "\n"
              << "  max #matches per query:   " << options.maxMatches << "\n"
              << " OUTPUT (stdout)\n"
              << "  stdout is terminal:       " << options.isTerm << "\n"
              << "  terminal width:           " << options.terminalCols << "\n"
              << "  verbosity:                " << options.verbosity << "\n"
              << " GENERAL\n"
              << "  double indexing:          " << options.doubleIndexing << "\n"
              << "  threads:                  " << uint(options.threads) << "\n"
              << "  query partitions:         " << (options.doubleIndexing
                                                    ? std::to_string(options.queryPart)
                                                    : std::string("n/a")) << "\n"
              << " TRANSLATION AND ALPHABETS\n"
              << "  genetic code:             "
              << ((p != BlastFormatProgram::BLASTN) &&
                  (p != BlastFormatProgram::BLASTP)
                 ? std::to_string(options.geneticCode)
                 : std::string("n/a")) << "\n"
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
              << " MISCELLANEOUS HEURISTICS\n"
              << "  pre-scoring:              " << (options.preScoring
                                                    ? std::string("on")
                                                    : std::string("off")) << "\n"
              << "  pre-scoring-region:       " << (options.preScoring
                                                    ? std::to_string(
                                                        options.preScoring *
                                                        options.seedLength)
                                                    : std::string("n/a")) << "\n"
              << "  pre-scoring-threshold:    " << (options.preScoring
                                                    ? std::to_string(
                                                       options.preScoringThresh)
                                                    : std::string("n/a")) << "\n"
              << "  putative-abundancy:       " << (options.filterPutativeAbundant
                                                    ? std::string("on")
                                                    : std::string("off")) << "\n"
              << "  putative-duplicates:      " << (options.filterPutativeDuplicates
                                                    ? std::string("on")
                                                    : std::string("off")) << "\n"
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
              << " MISC\n"
              << "  bit-compressed strings:   "
    #if defined LAMBDA_BITCOPMRESSED_STRINGS
              << "on\n"
    #else
              << "off\n"
    #endif
              << "\n";
}


#endif // header guard
