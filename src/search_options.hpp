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
// options.h: contains the options and argument parser
// ==========================================================================

#ifndef LAMBDA_SEARCH_OPTIONS_H_
#define LAMBDA_SEARCH_OPTIONS_H_

#include <cstdio>
#include <unistd.h>
#include <bitset>

#include <seqan/basic.h>
#include <seqan/translation.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/blast.h>

using namespace seqan;

// ==========================================================================
// Forwards
// ==========================================================================

template <typename T>
struct SamBamExtraTags;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LambdaOptions
// --------------------------------------------------------------------------

struct LambdaOptions : public SharedOptions
{

    std::string     queryFile;

    AlphabetEnum    qryOrigAlphabet;
    bool            revComp     = true;

    int             outFileFormat; // 0 = BLAST, 1 = SAM, 2 = BAM
    std::string     output;
    std::vector<BlastMatchField<>::Enum> columns;
    std::string     outputBam;
    std::bitset<64> samBamTags;
    bool            samWithRefHeader;
    unsigned        samBamSeq;
    bool            samBamHardClip;
    bool            versionInformationToOutputFile;

    unsigned        queryPart = 0;

//     bool            semiGlobal;

    bool            doubleIndexing = false;
    bool            adaptiveSeeding;

    unsigned        seedLength  = 0;
    unsigned        maxSeedDist = 1;
    bool            hammingOnly = true;

    int             seedGravity     = 0;
    unsigned        seedOffset      = 0;
    unsigned        minSeedLength   = 0;
    bool            seedDeltaIncreasesLength = true;

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
    int32_t         minBitScore = -1;
    double          maxEValue   = 1e-04;
    int             idCutOff    = 0;
    unsigned long   maxMatches  = 500;

    bool            computeLCA  = false;
    GeneticCodeSpec geneticCodeIndex;

    enum class ExtensionMode : uint8_t
    {
        AUTO,
        XDROP,
        FULL_SERIAL,
        FULL_SIMD
    };
    ExtensionMode   extensionMode;

    bool            filterPutativeDuplicates = true;
    bool            filterPutativeAbundant = true;
    bool            mergePutativeSiblings = true;
    bool            seedHalfExact = false;

    int             preScoring = 0; // 0 = off, 1 = seed, 2 = region (
    double          preScoringThresh    = 0.0;

    LambdaOptions() :
        SharedOptions()
    {
    }
};

ArgumentParser::ParseResult
parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // save commandLine
    for (int i = 0; i < argc; ++i)
        options.commandLine += std::string(argv[i]) + " ";
    eraseBack(options.commandLine);

    std::string programName = "lambda2 " + std::string(argv[0]);

    // this is important for option handling:
    if (std::string(argv[0]) == "searchn")
        options.blastProgram = BlastProgram::BLASTN;

    ArgumentParser parser(programName);
    // Set short description, version, and date.
    setShortDescription(parser, "the Local Aligner for Massive Biological DatA");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP "
                         "\\fI-i INDEX.lambda\\fP "
                         "[\\fI-o output.m8\\fP]");

    sharedSetup(parser);

#ifndef SEQAN_DISABLE_VERSION_CHECK
    // version checker initiated by top-level arg parser
    setDefaultValue(parser, "version-check", "0");
    hideOption(parser, "version-check");
#endif

    addOption(parser, ArgParseOption("v", "verbosity",
        "Display more/less diagnostic output during operation: 0 [only errors]; 1 [default]; 2 "
        "[+run-time, options and statistics].",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbosity", "1");
    setMinValue(parser, "verbosity", "0");
    setMaxValue(parser, "verbosity", "2");

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("q", "query",
        "Query sequences.",
        ArgParseArgument::INPUT_FILE,
        "IN"));
    setValidValues(parser, "query", getFileExtensions(SeqFileIn()));
    setRequired(parser, "q");

    if (options.blastProgram != BlastProgram::BLASTN)
    {
        addOption(parser, ArgParseOption("a", "input-alphabet",
            "Alphabet of the query sequences (specify to override auto-detection). Dna sequences will be translated.",
            ArgParseArgument::STRING));
        setValidValues(parser, "input-alphabet", "auto dna5 aminoacid");
        setDefaultValue(parser, "input-alphabet", "auto");
        setAdvanced(parser, "input-alphabet");

        addOption(parser, ArgParseOption("g", "genetic-code",
            "The translation table to use if input is Dna. See "
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
            " for ids. Default is to use the same table that was used for the index or 1/CANONICAL if the index "
            "was not translated.",
            ArgParseArgument::INTEGER));
        setDefaultValue(parser, "genetic-code", "0");
        setAdvanced(parser, "genetic-code");
    }

    addOption(parser, ArgParseOption("i", "index",
        std::string{"The database index (created by the 'lambda "} +
        (options.blastProgram == BlastProgram::BLASTN ? "mkindexn" : "mkindexp") +
        "' command).",
        ArgParseArgument::INPUT_DIRECTORY,
        "IN"));
    setRequired(parser, "index");
    setValidValues(parser, "index", ".lambda");

    addSection(parser, "Output Options");
    addOption(parser, ArgParseOption("o", "output",
        "File to hold reports on hits (.m* are blastall -m* formats; .m8 is tab-separated, .m9 is tab-separated with "
        "with comments, .m0 is pairwise format).",
        ArgParseArgument::OUTPUT_FILE,
        "OUT"));
    auto exts = getFileExtensions(BlastTabularFileOut<>());
    append(exts, getFileExtensions(BlastReportFileOut<>()));
    append(exts, getFileExtensions(BamFileOut()));
    CharString extsConcat;
    // remove .sam.bam, .sam.vcf.gz, .sam.tbi
    for (auto const & ext : exts)
    {
        if ((!endsWith(ext, ".bam") || startsWith(ext, ".bam")) &&
            (!endsWith(ext, ".vcf.gz")) &&
            (!endsWith(ext, ".sam.tbi")))
        {
            append(extsConcat, ext);
            appendValue(extsConcat, ' ');
        }
    }
    setValidValues(parser, "output", toCString(extsConcat));
    setDefaultValue(parser, "output", "output.m8");

    addOption(parser, ArgParseOption("", "output-columns",
        "Print specified column combination and/or order (.m8 and .m9 outputs only); call -oc help for more details.",
        ArgParseArgument::STRING,
        "STR"));
    setDefaultValue(parser, "output-columns", "std");
    setAdvanced(parser, "output-columns");

    addOption(parser, ArgParseOption("", "percent-identity",
        "Output only matches above this threshold (checked before e-value "
        "check).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "percent-identity", "0");
    setMinValue(parser, "percent-identity", "0");
    setMaxValue(parser, "percent-identity", "100");

    addOption(parser, ArgParseOption("e", "e-value",
        "Output only matches that score below this threshold.",
        ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e-value", "1e-04");
    setMinValue(parser, "e-value", "0");
    setMaxValue(parser, "e-value", "100");

    addOption(parser, ArgParseOption("n", "num-matches",
        "Print at most this number of matches per query.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "num-matches", "256");
    setMinValue(parser, "num-matches", "1");
    setMaxValue(parser, "num-matches", "10000");

    addOption(parser, ArgParseOption("", "sam-with-refheader",
        "BAM files require all subject names to be written to the header. For SAM this is not required, so Lambda does "
        "not automatically do it to save space (especially for protein database this is a lot!). If you still want "
        "them with SAM, e.g. for better BAM compatibility, use this option.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "sam-with-refheader", "off");
    setAdvanced(parser, "sam-with-refheader");

    std::string samBamSeqDescr;

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        samBamSeqDescr = "Write matching DNA subsequence into SAM/BAM file.";
    }
    else
    {
        samBamSeqDescr = "For BLASTX and TBLASTX the matching protein "
        "sequence is \"untranslated\" and positions retransformed to the original sequence. For BLASTP and TBLASTN "
        "there is no DNA sequence so a \"*\" is written to the SEQ column. The matching protein sequence can be "
        "written as an optional tag, see --sam-bam-tags.";
    }

    addOption(parser, ArgParseOption("", "sam-bam-seq",  samBamSeqDescr + " If set to uniq than "
        "the sequence is omitted iff it is identical to the previous match's subsequence.",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "sam-bam-seq", "always uniq never");
    setDefaultValue(parser, "sam-bam-seq", "uniq");
    setAdvanced(parser, "sam-bam-seq");

    addOption(parser, ArgParseOption("", "sam-bam-tags",
        "Write the specified optional columns to the SAM/BAM file. Call --sam-bam-tags help for more details.",
        ArgParseArgument::STRING,
        "STR"));
    setDefaultValue(parser, "sam-bam-tags", "AS NM ae ai qf");
    setAdvanced(parser, "sam-bam-tags");

    addOption(parser, ArgParseOption("", "sam-bam-clip",
        "Whether to hard-clip or soft-clip the regions beyond the local match. Soft-clipping retains the full sequence "
        "in the output file, but obviously uses more space.",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "sam-bam-clip", "hard soft");
    setDefaultValue(parser, "sam-bam-clip", "hard");
    setAdvanced(parser, "sam-bam-clip");

    addOption(parser, ArgParseOption("", "version-to-outputfile",
        "Write the Lambda program tag and version number to the output file.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "version-to-outputfile", "on");
    hideOption(parser, "version-to-outputfile");

    addSection(parser, "General Options");
#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads",
        "number of threads to run concurrently.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", omp_get_max_threads());
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", std::to_string(omp_get_max_threads() * 10));
#else
    addOption(parser, ArgParseOption("t", "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", "1");
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "1");
#endif
    setAdvanced(parser, "threads");

#ifdef LAMBDA_LEGACY_PATHS
    addOption(parser, ArgParseOption("", "query-index-type",
        "controls double-indexing.",
        ArgParseArgument::STRING));
    setValidValues(parser, "query-index-type", "radix none");
    setDefaultValue(parser, "query-index-type", "none");
    setAdvanced(parser, "query-index-type");

    addOption(parser, ArgParseOption("", "query-partitions",
        "Divide the query into qp number of blocks before processing; should be"
        " a multiple of the number of threads, defaults to one per thread. "
        "Only used with double-indexing; strong influence on memory, see below.",
        ArgParseArgument::INTEGER));
#ifdef _OPENMP
    setDefaultValue(parser, "query-partitions", omp_get_max_threads());
#else
    setDefaultValue(parser, "query-partitions", "1");
#endif // _OPENMP
    hideOption(parser, "query-partitions"); // HIDDEN
#endif // LAMBDA_LEGACY_PATHS


    addSection(parser, "Seeding / Filtration");

    addOption(parser, ArgParseOption("", "adaptive-seeding",
        "Grow the seed if it has too many hits (low complexity filter).",
        ArgParseArgument::BOOL));
    if (options.blastProgram == BlastProgram::BLASTN)
        setDefaultValue(parser, "adaptive-seeding", "off");
    else
        setDefaultValue(parser, "adaptive-seeding", "on");
    setAdvanced(parser, "adaptive-seeding");

    unsigned defaultSeedLength = (options.blastProgram == BlastProgram::BLASTN) ? 14 : 10;
    addOption(parser, ArgParseOption("", "seed-length",
        "Length of the seeds.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-length", std::to_string(defaultSeedLength));
    setMinValue(parser, "seed-length", "3");
    setMaxValue(parser, "seed-length", "50");
    setAdvanced(parser, "seed-length");

    addOption(parser, ArgParseOption("", "seed-offset",
        "Offset for seeding (if unset = seed-length/2).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-offset", std::to_string(defaultSeedLength / 2));
    setAdvanced(parser, "seed-offset");
    setMinValue(parser, "seed-offset", "1");
    setMaxValue(parser, "seed-offset", "50");

    addOption(parser, ArgParseOption("", "seed-delta",
        "maximum seed distance.",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-delta", "1");
    setAdvanced(parser, "seed-delta");
    setMinValue(parser, "seed-delta", "0");
    setMaxValue(parser, "seed-delta", "1");

    addOption(parser, ArgParseOption("", "seed-delta-increases-length",
        "Seed delta increases the min. seed length (for affected seeds).",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "seed-delta-increases-length", "off");
    setAdvanced(parser, "seed-delta-increases-length");

    addOption(parser, ArgParseOption("", "seed-half-exact",
        "Allow errors only in second half of seed.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "seed-half-exact", "on");
    setAdvanced(parser, "seed-half-exact");

    addOption(parser, ArgParseOption("", "seed-gravity",
        "Seeds closer than this are merged into region (if unset = "
        "seed-length).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-gravity", "10");
    hideOption(parser, "seed-gravity"); // HIDDEN

    addOption(parser, ArgParseOption("", "seed-min-length",
        "after postproc shorter seeds are discarded (if unset = seed-length).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seed-min-length", "10");
    hideOption(parser, "seed-min-length"); // HIDDEN

    addSection(parser, "Miscellaneous Heuristics");

    addOption(parser, ArgParseOption("", "pre-scoring",
        "evaluate score of a region NUM times the size of the seed "
        "before extension (0 -> no pre-scoring, 1 -> evaluate seed, n-> area "
        "around seed, as well; default = 1 if no reduction is used).",
        ArgParseArgument::INTEGER));
    setMinValue(parser, "pre-scoring", "1");
    setMaxValue(parser, "pre-scoring", "10");
    setDefaultValue(parser, "pre-scoring", "2");
    setAdvanced(parser, "pre-scoring");

    addOption(parser, ArgParseOption("", "pre-scoring-threshold",
        "minimum average score per position in pre-scoring region.",
        ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "pre-scoring-threshold", "2");
    setMinValue(parser, "pre-scoring-threshold", "0");
    setMaxValue(parser, "pre-scoring-threshold", "20");
    setAdvanced(parser, "pre-scoring-threshold");

    addOption(parser, ArgParseOption("", "filter-putative-duplicates",
        "filter hits that will likely duplicate a match already found.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "filter-putative-duplicates", "on");
    setAdvanced(parser, "filter-putative-duplicates");

    addOption(parser, ArgParseOption("", "filter-putative-abundant",
        "If the maximum number of matches per query are found already, "
        "stop searching if the remaining realm looks unfeasible.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "filter-putative-abundant", "on");
    setAdvanced(parser, "filter-putative-abundant");

    addOption(parser, ArgParseOption("", "merge-putative-siblings",
        "Merge seed from one region, "
        "stop searching if the remaining realm looks unfeasable.",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "merge-putative-siblings", "on");
    setAdvanced(parser, "merge-putative-siblings");

    addSection(parser, "Scoring");

    if (options.blastProgram != BlastProgram::BLASTN)
    {
        addOption(parser, ArgParseOption("s", "scoring-scheme",
            "use '45' for Blosum45; '62' for Blosum62 (default); '80' for Blosum80.",
            ArgParseArgument::INTEGER));
        setDefaultValue(parser, "scoring-scheme", "62");
        setAdvanced(parser, "scoring-scheme");
    }

    addOption(parser, ArgParseOption("", "score-gap",
        "Score per gap character.",
        ArgParseArgument::INTEGER));
    if (options.blastProgram == BlastProgram::BLASTN)
        setDefaultValue(parser, "score-gap", "-2");
    else
        setDefaultValue(parser, "score-gap", "-1");
    setMinValue(parser, "score-gap", "-1000");
    setMaxValue(parser, "score-gap", "1000");
    setAdvanced(parser, "score-gap");

    addOption(parser, ArgParseOption("", "score-gap-open",
        "Additional cost for opening gap.",
        ArgParseArgument::INTEGER));
    if (options.blastProgram == BlastProgram::BLASTN)
        setDefaultValue(parser, "score-gap-open", "-5");
    else
        setDefaultValue(parser, "score-gap-open", "-11");
    setMinValue(parser, "score-gap-open", "-1000");
    setMaxValue(parser, "score-gap-open", "1000");
    setAdvanced(parser, "score-gap-open");

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        addOption(parser, ArgParseOption("", "score-match",
            "Match score [only BLASTN])",
            ArgParseArgument::INTEGER));
        setDefaultValue(parser, "score-match", "2");
        setMinValue(parser, "score-match", "-1000");
        setMaxValue(parser, "score-match", "1000");
        setAdvanced(parser, "score-match");

        addOption(parser, ArgParseOption("", "score-mismatch",
            "Mismatch score [only BLASTN]",
            ArgParseArgument::INTEGER));
        setDefaultValue(parser, "score-mismatch", "-3");
        setMinValue(parser, "score-mismatch", "-1000");
        setMaxValue(parser, "score-mismatch", "1000");
        setAdvanced(parser, "score-mismatch");
    }

    addSection(parser, "Extension");

    addOption(parser, ArgParseOption("x", "x-drop",
        "Stop Banded extension if score x below the maximum seen (-1 means no "
        "xdrop).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "x-drop", "30");
    setMinValue(parser, "x-drop", "-1");
    setMaxValue(parser, "x-drop", "1000");
    setAdvanced(parser, "x-drop");

    addOption(parser, ArgParseOption("b", "band",
        "Size of the DP-band used in extension (-3 means log2 of query length; "
        "-2 means sqrt of query length; -1 means full dp; n means band of size "
        "2n+1)",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "band", "-3");
    setMinValue(parser, "band", "-3");
    setMaxValue(parser, "band", "1000");
    setAdvanced(parser, "band");

    addOption(parser, ArgParseOption("m", "extension-mode",
        "Choice of extension algorithms.",
        ArgParseArgument::STRING));
#ifdef SEQAN_SIMD_ENABLED
    setValidValues(parser, "extension-mode", "auto xdrop fullSerial fullSIMD");
#else
    setValidValues(parser, "extension-mode", "auto xdrop fullSerial");
#endif
    setDefaultValue(parser, "extension-mode", "auto");
    setAdvanced(parser, "extension-mode");

    addTextSection(parser, "Tuning");
    addText(parser, "Tuning the seeding parameters and (de)activating alphabet "
                    "reduction has a strong "
                    "influence on both speed and sensitivity. We recommend the "
                    "following alternative profiles for protein searches:");
    addText(parser, "fast (high similarity):       --seed-delta-increases-length on");
    addText(parser, "sensitive (lower similarity): --seed-offset 3");

    addText(parser, "For further information see the wiki: <https://github.com/seqan/lambda/wiki>");
//         addTextSection(parser, "Speed VS memory requirements");
//         addText(parser, "Lambda requires approximately the following amount of RAM:"
//                         " \033[1msize(queryFile) + size(dbIDs) + 2 * size(dbSeqs)\033[0m. "
//                         "If you have more RAM, use double indexing and SA:\n"
//                         "\033[1m-di sa -qi radix\033[0m "
//                         "which will result in an additional speed-up of up to 30% "
//                         "compared to the published version (you need to run the "
//                         "indexer with \033[1m-di sa \033[0m, as well). The amount "
//                         "of RAM required will be: "
//                         "\033[1msize(queryFile) + size(dbIDs) + 7 * size(dbSeqs) + n\033[0m "
//                         "where n grows slowly but linearly with input size. "
//                         "Note that size(dbSeqs) refers to the total "
//                         "sequence length and does not include IDs (so it is less "
//                         "than the size of the file).");
//         addText(parser, "To save more RAM, you can define "
//                         "LAMBDA_BITCOPMRESSED_STRINGS while compiling lambda. "
//                         "This will reduce memory usage by about:"
//                         " \033[1m0.3 * ( size(queryFile) + size(dbSeqs) )\033[0m,"
//                         " but slow down lambda by about 10%.");

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

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        options.qryOrigAlphabet = AlphabetEnum::DNA5;
    }
    else
    {
        getOptionValue(buffer, parser, "input-alphabet");
        if (buffer == "auto")
            options.qryOrigAlphabet = AlphabetEnum::DNA4;
        else if (buffer == "dna5")
            options.qryOrigAlphabet = AlphabetEnum::DNA5;
        else if (buffer == "aminoacid")
            options.qryOrigAlphabet = AlphabetEnum::AMINO_ACID;
        else
            throw std::invalid_argument("ERROR: Invalid argument to --input-alphabet\n");

        int buf = 0;
        getOptionValue(buf, parser, "genetic-code");
        switch (buf)
        {
            case 0: // take code from index
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
    }

    getOptionValue(options.indexDir, parser, "index");

    getOptionValue(options.output, parser, "output");
    buffer = options.output;
    if (endsWith(buffer, ".gz"))
        buffer.resize(length(buffer) - 3);
    else if (endsWith(buffer, ".bz2"))
        buffer.resize(length(buffer) - 4);

    if (endsWith(buffer, ".sam"))
        options.outFileFormat = 1;
    else if (endsWith(buffer, ".bam"))
        options.outFileFormat = 2;
    else
        options.outFileFormat = 0;

    getOptionValue(options.samWithRefHeader, parser, "sam-with-refheader");

    clear(buffer);
    getOptionValue(buffer, parser, "sam-bam-seq");
    if (buffer == "never")
        options.samBamSeq = 0;
    else if (buffer == "uniq")
        options.samBamSeq = 1;
    else
        options.samBamSeq = 2;

    clear(buffer);
    getOptionValue(buffer, parser, "sam-bam-clip");
    options.samBamHardClip = (buffer == "hard");

    clear(buffer);
    getOptionValue(buffer, parser, "output-columns");
    if (buffer == "help")
    {
        std::cout << "Please specify the columns in this format -oc 'column1 column2', i.e. space-separated and "
                  << "enclosed in single quotes.\nThe specifiers are the same as in NCBI Blast, currently "
                  << "the following are supported:\n";
        for (unsigned i = 0; i < length(BlastMatchField<>::implemented); ++i)
        {
            if (BlastMatchField<>::implemented[i])
            {
                std::cout << "\t" << BlastMatchField<>::optionLabels[i]
                          << (length(BlastMatchField<>::optionLabels[i]) >= 8 ? "\t" : "\t\t")
                          << BlastMatchField<>::descriptions[i] << "\n";
            }
        }
        return ArgumentParser::PARSE_HELP;
    }
    else
    {
        StringSet<CharString> fields;
        strSplit(fields, buffer, IsSpace(), false);
        for (auto str : fields)
        {
            bool resolved = false;
            for (unsigned i = 0; i < length(BlastMatchField<>::optionLabels); ++i)
            {
                if (BlastMatchField<>::optionLabels[i] == str)
                {
                    appendValue(options.columns, static_cast<BlastMatchField<>::Enum>(i));
                    resolved = true;
                    if (static_cast<BlastMatchField<>::Enum>(i) == BlastMatchField<>::Enum::S_TAX_IDS)
                        options.hasSTaxIds = true;
                    else if ((static_cast<BlastMatchField<>::Enum>(i) == BlastMatchField<>::Enum::LCA_ID) ||
                             (static_cast<BlastMatchField<>::Enum>(i) == BlastMatchField<>::Enum::LCA_TAX_ID))
                        options.computeLCA = true;
                    break;
                }
            }
            if (!resolved)
            {
                std::cerr << "Unknown column specifier \"" << str << "\". Please see -oc help for valid options.\n";
                return ArgumentParser::PARSE_ERROR;
            }
        }
    }
    clear(buffer);

    getOptionValue(buffer, parser, "sam-bam-tags");
    if (buffer == "help")
    {
        std::cout << "Please specify the tags in this format -oc 'tag1 tag2', i.e. space-separated and "
                  << "enclosed in quotes. The order of tags is not preserved.\nThe following specifiers are "
                  << "supported:\n";

        for (auto const & c : SamBamExtraTags<>::keyDescPairs)
            std::cout << "\t" << std::get<0>(c) << "\t" << std::get<1>(c) << "\n";

        return ArgumentParser::PARSE_HELP;
    }
    else
    {
        StringSet<CharString> fields;
        strSplit(fields, buffer, IsSpace(), false);
        for (auto str : fields)
        {
            bool resolved = false;
            for (unsigned i = 0; i < length(SamBamExtraTags<>::keyDescPairs); ++i)
            {
                if (std::get<0>(SamBamExtraTags<>::keyDescPairs[i]) == str)
                {
                    options.samBamTags[i] = true;
                    resolved = true;
                    break;
                }
            }
            if (!resolved)
            {
                std::cerr << "Unknown column specifier \"" << str
                          << "\". Please see \"--sam-bam-tags help\" for valid options.\n";
                return ArgumentParser::PARSE_ERROR;
            }
        }
    }

    if (options.samBamTags[SamBamExtraTags<>::S_TAX_IDS])
        options.hasSTaxIds = true;

    if (options.samBamTags[SamBamExtraTags<>::LCA_ID] || options.samBamTags[SamBamExtraTags<>::LCA_TAX_ID])
        options.computeLCA = true;

    // lca computation requires tax ids
    if (options.computeLCA)
        options.hasSTaxIds = true;

    getOptionValue(options.versionInformationToOutputFile, parser, "version-to-outputfile");
    getOptionValue(options.adaptiveSeeding, parser, "adaptive-seeding");

    clear(buffer);
    getOptionValue(options.seedLength, parser, "seed-length");

    getOptionValue(options.seedOffset, parser, "seed-offset");

    if (isSet(parser, "seed-gravity"))
        getOptionValue(options.seedGravity, parser, "seed-gravity");
    else
        options.seedGravity = options.seedLength;

    if (isSet(parser, "seed-min-length"))
        getOptionValue(options.minSeedLength, parser, "seed-min-length");
    else
        options.minSeedLength = options.seedLength;

    getOptionValue(options.maxSeedDist, parser, "seed-delta");

    if (options.maxSeedDist == 0)
    {
        // the whole seed is exact, so it is also half-exact :)
        options.seedHalfExact = true;

        if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
        {
            std::cerr << "WARNING: Exact seeeding doesn't benefit from bi-fm-index, so regular index is used.\n";
            options.dbIndexType = DbIndexType::FM_INDEX;
        }
    }

    getOptionValue(options.seedDeltaIncreasesLength, parser, "seed-delta-increases-length");

    getOptionValue(options.eCutOff, parser, "e-value");
    getOptionValue(options.idCutOff, parser, "percent-identity");

    getOptionValue(options.xDropOff, parser, "x-drop");

    getOptionValue(options.band, parser, "band");

#ifdef LAMBDA_LEGACY_PATHS
    getOptionValue(buffer, parser, "query-index-type");
    options.doubleIndexing = (buffer == "radix");

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
#endif

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        options.scoringMethod = 0;
        getOptionValue(options.misMatch, parser, "score-mismatch");
        getOptionValue(options.match, parser, "score-match");
    }
    else
    {
        getOptionValue(options.scoringMethod, parser, "scoring-scheme");

        switch (options.scoringMethod)
        {
            case 45: case 62: case 80: break;
            default:
                std::cerr << "Unsupported Scoring Scheme selected.\n";
                return ArgumentParser::PARSE_ERROR;
        }
    }

    getOptionValue(options.gapExtend, parser, "score-gap");

    getOptionValue(options.gapOpen, parser, "score-gap-open");

    getOptionValue(options.filterPutativeDuplicates, parser, "filter-putative-duplicates");
    getOptionValue(options.filterPutativeAbundant, parser, "filter-putative-abundant");
    getOptionValue(options.mergePutativeSiblings, parser, "merge-putative-siblings");
    getOptionValue(options.seedHalfExact, parser, "seed-half-exact");

    if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
    {
        if (options.seedHalfExact)
            std::cerr << "WARNING: seedHalfExact is already implied by bidirectional indexes.\n";
        else
            options.seedHalfExact = true;
    }

    // TODO always prescore 1
    getOptionValue(options.preScoring, parser, "pre-scoring");

    if ((!isSet(parser, "pre-scoring")) &&
        (options.reducedAlphabet == options.transAlphabet))
        options.preScoring = 1;

    getOptionValue(options.preScoringThresh, parser, "pre-scoring-threshold");
//     if (options.preScoring == 0)
//         options.preScoringThresh = 4;

    int numbuf;
    getOptionValue(numbuf, parser, "num-matches");
    options.maxMatches = static_cast<unsigned long>(numbuf);

    getOptionValue(buffer, parser, "extension-mode");
    if (buffer == "fullSIMD")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::FULL_SIMD;
        options.filterPutativeAbundant = false;
        options.filterPutativeDuplicates = false;
        options.mergePutativeSiblings = false;
        options.xDropOff = -1;
    }
    else if (buffer == "fullSerial")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::FULL_SERIAL;
        options.filterPutativeAbundant = false;
        options.filterPutativeDuplicates = false;
        options.mergePutativeSiblings = false;
        options.xDropOff = -1;
    }
    else if (buffer == "xdrop")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::XDROP;
    }
    else
    {
        options.extensionMode = LambdaOptions::ExtensionMode::AUTO;
    }

    return ArgumentParser::PARSE_OK;
}


// --------------------------------------------------------------------------
// Function printOptions()
// --------------------------------------------------------------------------

template <typename TLH>
inline void
printOptions(LambdaOptions const & options)
{
    using TGH = typename TLH::TGlobalHolder;

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
              << "  index directory:          " << options.indexDir << "\n"
              << "  db index type:            " << _indexEnumToName(options.dbIndexType) << "\n"
              << " OUTPUT (file)\n"
              << "  output file:              " << options.output << "\n"
              << "  maximum e-value:          " << options.maxEValue << "\n"
              << "  minimum bit-score:        " << options.minBitScore << "\n"
              << "  minimum % identity:       " << options.idCutOff << "\n"
              << "  max #matches per query:   " << options.maxMatches << "\n"
              << "  include subj names in sam:" << options.samWithRefHeader << "\n"
              << "  include seq in sam/bam:   " << options.samBamSeq << "\n"
              << "  with subject tax ids:     " << options.hasSTaxIds << '\n'
              << "  compute LCA:              " << options.computeLCA << '\n'
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
              << ((TGH::blastProgram != BlastProgram::BLASTN) &&
                  (TGH::blastProgram != BlastProgram::BLASTP)
                 ? std::to_string(options.geneticCode)
                 : std::string("n/a")) << "\n"
              << "  blast mode:               " << _programTagToString(TGH::blastProgram)
              << "\n"
              << "  original alphabet (query):" << _alphTypeToName(OrigQryAlph<TGH::blastProgram>())
              << "\n"
              << "  original alphabet (subj): " << _alphTypeToName(OrigSubjAlph<TGH::blastProgram>())
              << "\n"
              << "  translated alphabet:      " << _alphTypeToName(TransAlph<TGH::blastProgram>())
              << "\n"
              << "  reduced alphabet:         " << _alphTypeToName(typename TGH::TRedAlph())
              << "\n"
              << " SEEDING\n"
              << "  seed length:              " << uint(options.seedLength) << "\n"
              << "  seed offset:              " << uint(options.seedOffset) << "\n"
              << "  seed delta:               " << uint(options.maxSeedDist) << "\n"
              << "  seeds ungapped:           " << uint(options.hammingOnly) << "\n"
              << "  seed gravity:             " << uint(options.seedGravity) << "\n"
              << "  min seed length:          " << uint(options.minSeedLength) << "\n"
              << "  seed delta length inc.:   " << (options.seedDeltaIncreasesLength
                                                    ? std::string("on")
                                                    : std::string("off")) << "\n"
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
              << "  seed half exact:          " << (options.seedHalfExact
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
              << " EXTENSION\n";
    switch (options.extensionMode)
    {
        case LambdaOptions::ExtensionMode::AUTO:
            std::cout
              << "  extensionMode:            auto (depends on query length)\n"
              << "  x-drop:                   " << options.xDropOff << "\n"
              << "  band:                     " << bandStr << "\n"
              << "  [depending on the automatically chosen mode x-drop or band might get disabled.\n";
              break;
        case LambdaOptions::ExtensionMode::XDROP:
            std::cout
              << "  extensionMode:            individual\n"
              << "  x-drop:                   " << options.xDropOff << "\n"
              << "  band:                     " << bandStr << "\n";
            break;
        case LambdaOptions::ExtensionMode::FULL_SERIAL:
            std::cout
              << "  extensionMode:            batch, but serialized\n"
              << "  x-drop:                   not used\n"
              << "  band:                     " << bandStr << "\n";
            break;
        case LambdaOptions::ExtensionMode::FULL_SIMD:
            std::cout
              << "  extensionMode:            batch with SIMD\n"
              << "  x-drop:                   not used\n"
              << "  band:                     not used\n";
            break;
    }
    std::cout << " BUILD OPTIONS:\n"
              << "  cmake_build_type:         " << std::string(CMAKE_BUILD_TYPE) << "\n"
              << "  fastbuild:                "
    #if defined(FASTBUILD)
              << "on\n"
    #else
              << "off\n"
    #endif
              << "  native_build:             "
    #if defined(LAMBDA_NATIVE_BUILD)
              << "on\n"
    #else
              << "off\n"
    #endif
              << "  static_build:             "
    #if defined(LAMBDA_STATIC_BUILD)
              << "on\n"
    #else
              << "off\n"
    #endif
              << "  mmapped_db:               "
    #if defined(LAMBDA_MMAPPED_DB)
              << "on\n"
    #else
              << "off\n"
    #endif
              << "  lingaps_opt:              "
    #if defined(LAMBDA_LINGAPS_OPT)
              << "on\n"
    #else
              << "off\n"
    #endif
              << "  seqan_simd:               "
    #if defined(SEQAN_SIMD_ENABLED) && defined(__AVX2__)
              << "avx2\n"
    #elif defined(SEQAN_SIMD_ENABLED) && defined(__SSE4_2__)
              << "sse4\n"
    #else
              << "off\n"
    #endif
              << "\n";
}

#endif // header guard
