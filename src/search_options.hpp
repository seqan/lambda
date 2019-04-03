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

#include <seqan3/std/filesystem>

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

    int32_t         outFileFormat; // 0 = BLAST, 1 = SAM, 2 = BAM
    std::string     output = "output.m8";
    std::vector<BlastMatchField<>::Enum> columns;
    std::string     outputBam;
    std::bitset<64> samBamTags;
    bool            samWithRefHeader = false;
    unsigned        samBamSeq;
    bool            samBamHardClip;
    bool            versionInformationToOutputFile = true;

//     bool            semiGlobal;

    bool            adaptiveSeeding;

    unsigned        seedLength  = 0;
    unsigned        maxSeedDist = 1;
    bool            hammingOnly = true;

    int32_t         seedGravity     = 10;
    unsigned        seedOffset      = 0;
    unsigned        minSeedLength   = 0;
    bool            seedDeltaIncreasesLength = false;

//     unsigned int    minSeedEVal     = 0;
//     double          minSeedBitS     = -1;

    // 0 = manual, positive X = blosumX, negative Y = pamY
    int32_t         scoringMethod   = 62;
    // scores
    int32_t         gapOpen         = -11;
    int32_t         gapExtend       = -1;
    int32_t         match           = 2; // only for manual
    int32_t         misMatch        = -3; // only for manual

    int32_t         xDropOff    = 30;
    int32_t         band        = -3;
    double          eCutOff     = 1e-04;
    int32_t         idCutOff    = 0;
    uint64_t        maxMatches  = 256;

    bool            computeLCA  = false;
    GeneticCodeSpec geneticCodeIndex;

    enum class ExtensionMode : uint8_t
    {
        AUTO,
        XDROP,
        FULL_SERIAL,
        FULL_SIMD
    };
    ExtensionMode   extensionMode = ExtensionMode::AUTO;

    bool            filterPutativeDuplicates = true;
    bool            filterPutativeAbundant = true;
    bool            mergePutativeSiblings = true;
    bool            seedHalfExact = true;

    int32_t             preScoring = 2; // 0 = off, 1 = seed, 2 = region
    double          preScoringThresh    = 2.0;

    LambdaOptions() :
        SharedOptions()
    {
    }
};

void parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // save commandLine
    for (int i = 0; i < argc; ++i)
        options.commandLine += std::string(argv[i]) + " ";
    eraseBack(options.commandLine);

    std::string programName = "lambda3 " + std::string(argv[0]);

    // this is important for option handling:
    if (std::string(argv[0]) == "searchn")
        options.blastProgram = BlastProgram::BLASTN;

    seqan3::argument_parser parser(programName, argc, argv);

    // Set short description, version, and date.
    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back("[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP \\fI-i INDEX.lambda\\fP [\\fI-o output.m8\\fP]");

    sharedSetup(parser);

    // TODO version check

    parser.add_option(options.verbosity, 'v', "verbosity", "Display more/less diagnostic output during operation: "
        "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].",
        seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{0, 2});

    parser.add_section("Input options");

    // TODO Better solution for file extensions
    parser.add_option(options.queryFile, 'q', "query", "Query sequences.", seqan3::option_spec::REQUIRED,
        seqan3::path_existence_validator() | seqan3::file_ext_validator({"fa", "fasta", "fq", "fastq"}));

    std::string inputAlphabetTmp = "auto";
    int32_t geneticCodeTmp = 1;

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        options.qryOrigAlphabet = AlphabetEnum::DNA5;
    }
    else
    {
        parser.add_option(inputAlphabetTmp, 'a', "input-alphabet",
            "Alphabet of the query sequences (specify to override auto-detection). Dna sequences will be translated.",
            seqan3::option_spec::ADVANCED, seqan3::value_list_validator({"auto", "dna5", "aminoacid"}));

        parser.add_option(geneticCodeTmp, 'g', "genetic-code",
            "The translation table to use if input is Dna. See "
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
            " for ids. Default is to use the same table that was used for the index or 1/CANONICAL if the index "
            "was not translated.", seqan3::option_spec::ADVANCED,
            seqan3::value_list_validator({0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25}));
    }

    // TODO Does this input directory structure work?
    parser.add_option(options.indexDir, 'i', "index", std::string{"The database index (created by the 'lambda "} +
        (options.blastProgram == BlastProgram::BLASTN ? "mkindexn" : "mkindexp") + "' command).",
        seqan3::option_spec::REQUIRED, seqan3::path_existence_validator());

    parser.add_section("Output options");

    // TODO Fix outout file requirements
    parser.add_option(options.output, 'o', "output",
        "File to hold reports on hits (.m* are blastall -m* formats; .m8 is tab-seperated, .m9 is tab-seperated with "
        "with comments, .m0 is pairwise format).", seqan3::option_spec::DEFAULT, seqan3::path_existence_validator() |
        seqan3::file_ext_validator({"m0", "m8", "m9", "bam", "sam", "gz", "bz2"}));

    std::string outputColumnsTmp = "std";
    parser.add_option(outputColumnsTmp, '\0', "output-columns",
        "Print specified column combination and/or order (.m8 and .m9 outputs only); call -oc help for more details.",
        seqan3::option_spec::ADVANCED);

    parser.add_option(options.idCutOff, '\0', "percent-identity",
        "Output only matches above this threshold (checked before e-value check).", seqan3::option_spec::DEFAULT,
        seqan3::arithmetic_range_validator{0, 100});

    parser.add_option(options.eCutOff, 'e', "e-value", "Output only matches that score below this threshold.",
        seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{0, 100});

    int32_t numMatchesTmp = 256;
    parser.add_option(numMatchesTmp, 'n', "num-matches", "Print at most this number of matches per query.",
        seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{0, 10000});

    parser.add_option(options.samWithRefHeader, '\0', "sam-with-refheader",
        "BAM files require all subject names to be written to the header. For SAM this is not required, so Lambda does "
        "not automatically do it to save space (especially for protein database this is a lot!). If you still want "
        "them with SAM, e.g. for better BAM compatibility, use this option.", seqan3::option_spec::ADVANCED);

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

    std::string samBamSeqDescrTmp = "uniq";
    parser.add_option(samBamSeqDescrTmp, '\0', "sam-bam-seq", samBamSeqDescr + " If set to uniq than the sequence is "
        "omitted iff it is identical to the previous match's subsequence.", seqan3::option_spec::ADVANCED,
        seqan3::value_list_validator({"always", "uniq", "never"}));

    std::string samBamTagsTmp = "AS NM ae ai qf";
    parser.add_option(samBamTagsTmp, '\0', "sam-bam-tags",
        "Write the specified optional columns to the SAM/BAM file. Call --sam-bam-tags help for more details.",
        seqan3::option_spec::ADVANCED);

    std::string samBamClip = "hard";
    parser.add_option(samBamClip, '\0', "sam-bam-clip",
        "Whether to hard-clip or soft-clip the regions beyond the local match. Soft-clipping retains the full sequence "
        "in the output file, but obviously uses more space.", seqan3::option_spec::ADVANCED,
        seqan3::value_list_validator({"hard", "soft"}));

    parser.add_option(options.versionInformationToOutputFile, '\0', "version-to-outputfile",
        "Write the Lambda program tag and version number to the output file.", seqan3::option_spec::HIDDEN);

    parser.add_section("General Options");

#ifdef _OPENMP
    parser.add_option(options.threads, 't', "threads", "Number of threads to run concurrently.",
        seqan3::option_spec::ADVANCED,
        seqan3::arithmetic_range_validator{1, (double) omp_get_max_threads() * 10});
#else
    parser.add_option(options.threads, 't', "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.", seqan3::option_spec::ADVANCED,
        seqan3::arithmetic_range_validator{1, 1});
#endif

    parser.add_section("Seeding / Filtration");

    if (options.blastProgram == BlastProgram::BLASTN)
        options.adaptiveSeeding = false;
    else
        options.adaptiveSeeding = true;

    parser.add_option(options.adaptiveSeeding, '\0', "adaptive-seeding",
        "Grow the seed if it has too many hits (low complexity filter).", seqan3::option_spec::ADVANCED);

    unsigned defaultSeedLength = (options.blastProgram == BlastProgram::BLASTN) ? 14 : 10;

    options.seedLength = defaultSeedLength;
    parser.add_option(options.seedLength, '\0', "seed-length", "Length of the seeds.", seqan3::option_spec::ADVANCED,
        seqan3::arithmetic_range_validator{3, 50});

    options.seedOffset = defaultSeedLength / 2;
    parser.add_option(options.seedOffset, '\0', "seed-offset", "Offset for seeding (if unset = seed-length/2). "
        "If you set 'seed-length', please consider setting this option to 'seedLength / 2'.",
        seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 50});

    parser.add_option(options.maxSeedDist, '\0', "seed-delta",
        "Maximum seed distance.", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{0, 1});

    parser.add_option(options.seedDeltaIncreasesLength, '\0', "seed-delta-increases-length",
        "Seed delta increases the min. seed length (for affected seeds).", seqan3::option_spec::ADVANCED);

    parser.add_option(options.seedHalfExact, '\0', "seed-half-exact",
        "Allow errors only in second half of seed.", seqan3::option_spec::ADVANCED);

    options.seedGravity = defaultSeedLength;
    parser.add_option(options.seedGravity, '\0', "seed-gravity",
        "Seeds closer than this are merged into region (if unset = seed-length). If you set 'seed-length', please "
        "consider setting this option to the same value.", seqan3::option_spec::HIDDEN);

    options.seedGravity = defaultSeedLength;
    parser.add_option(options.minSeedLength, '\0', "seed-min-length",
        "After postproc shorter seeds are discarded (if unset = seed-length).  If you set 'seed-length', please "
        "consider setting this option to the same value.", seqan3::option_spec::HIDDEN);

    parser.add_section("Miscellaneous Heuristics");

    parser.add_option(options.preScoring, '\0', "Pre-scoring",
        "Evaluate score of a region NUM times the size of the seed "
        "before extension (0 -> no pre-scoring, 1 -> evaluate seed, n-> area "
        "around seed, as well; default = 1 if no reduction is used).",
        seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{1, 10});

    parser.add_option(options.preScoringThresh, '\0', "pre-scoring-threshold",
        "Minimum average score per position in pre-scoring region.", seqan3::option_spec::ADVANCED,
        seqan3::arithmetic_range_validator{0, 20});

    parser.add_option(options.filterPutativeDuplicates, '\0', "filter-putative-duplicates",
        "Filter hits that will likely duplicate a match already found.", seqan3::option_spec::ADVANCED);

    parser.add_option(options.filterPutativeAbundant, '\0', "filter-putative-abundant",
        "If the maximum number of matches per query are found already, "
        "stop searching if the remaining realm looks unfeasible.", seqan3::option_spec::ADVANCED);

    parser.add_option(options.mergePutativeSiblings, '\0', "merge-putative-siblings",
        "Merge seed from one region, "
        "stop searching if the remaining realm looks unfeasable.", seqan3::option_spec::ADVANCED);

    parser.add_section("Scoring");

    if (options.blastProgram == BlastProgram::BLASTN)
    {
        parser.add_option(options.match, '\0', "score-match",
            "Match score [only BLASTN]", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-1000, 1000});

        parser.add_option(options.misMatch, '\0', "score-mismatch",
            "Mismatch score [only BLASTN]", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-1000, 1000});
    }
    else
    {
        parser.add_option(options.scoringMethod,'s', "scoring-scheme",
            "Use '45' for Blosum45; '62' for Blosum62 (default); '80' for Blosum80.", seqan3::option_spec::ADVANCED,
            seqan3::value_list_validator({45, 62, 80}));
    }

    if (options.blastProgram == BlastProgram::BLASTN)
        options.gapExtend = -2;
    else
        options.gapExtend = -1;

    parser.add_option(options.gapExtend, '\0', "score-gap",
        "Score per gap character.", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-1000, 1000});

    if (options.blastProgram == BlastProgram::BLASTN)
        options.gapOpen = -5;
    else
        options.gapOpen = -11;

    parser.add_option(options.gapOpen, '\0', "score-gap-open",
        "Additional cost for opening gap.", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-1000, 1000});

    parser.add_section("Extension");

    parser.add_option(options.xDropOff, 'x', "x-drop", "Stop Banded extension if score x below the maximum seen "
        "(-1 means no xdrop).", seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-1, 1000});

    parser.add_option(options.band, 'b', "band", "Size of the DP-band used in extension (-3 means log2 of query length;"
        " -2 means sqrt of query length; -1 means full dp; n means band of size 2n+1)",
        seqan3::option_spec::ADVANCED, seqan3::arithmetic_range_validator{-3, 1000});

    std::string extensionModeTmp;

#ifdef SEQAN_SIMD_ENABLED
    parser.add_option(extensionModeTmp,'m', "extension-mode",
        "Choice of extension algorithms.", seqan3::option_spec::ADVANCED,
        seqan3::value_list_validator({"auto", "xdrop", "fullSerial", "fullSIMD"}));
#else
    parser.add_option(extensionModeTmp,'m', "extension-mode",
        "Choice of extension algorithms.", seqan3::option_spec::ADVANCED,
        seqan3::value_list_validator({"auto", "xdrop", "fullSerial"}));
#endif

    parser.add_section("Tuning");
    parser.add_line("Tuning the seeding parameters and (de)activating alphabet "
                    "reduction has a strong "
                    "influence on both speed and sensitivity. We recommend the "
                    "following alternative profiles for protein searches:", false);
    parser.add_line("fast (high similarity):    --seed-delta-increases-length on", false);
    parser.add_line("sensitive (lower similarity): --seed-offset 3", false);
    parser.add_line("For further information see the wiki: <https://github.com/seqan/lambda/wiki>", false);

    // parse command line.
    parser.parse();

    // set query alphabet and genetic code depending on options
    if (options.blastProgram != BlastProgram::BLASTN)
    {
        if (inputAlphabetTmp == "auto")
            options.qryOrigAlphabet = AlphabetEnum::DNA4;
        else if (inputAlphabetTmp == "dna5")
            options.qryOrigAlphabet = AlphabetEnum::DNA5;
        else if (inputAlphabetTmp == "aminoacid")
            options.qryOrigAlphabet = AlphabetEnum::AMINO_ACID;
        else
            throw seqan3::parser_invalid_argument("ERROR: Invalid argument to --input-alphabet\n");

        options.geneticCode = static_cast<GeneticCodeSpec>(geneticCodeTmp);
    }

    // set output file format
    std::string outputPath;
    if (std::filesystem::path(outputPath).extension() == ".gz")
        outputPath.resize(length(outputPath) - 3);
    else if (std::filesystem::path(outputPath).extension() ==  ".bz2")
        outputPath.resize(length(outputPath) - 4);

    if (std::filesystem::path(outputPath).extension() == ".sam")
        options.outFileFormat = 1;
    else if (std::filesystem::path(outputPath).extension() == ".bam")
        options.outFileFormat = 2;
    else
        options.outFileFormat = 0;

    // help page for output columns
    if (outputColumnsTmp == "help")
    {
        std::cout << "Please specify the columns in this format -oc 'column1 column2', i.e. space-seperated and "
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
        std::exit(0);
    }
    else
    {
        StringSet<CharString> fields;
        strSplit(fields, outputColumnsTmp, IsSpace(), false);
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
                throw seqan3::parser_invalid_argument(std::string("Unknown column specifier \"") + toCString(str) +
                std::string("\". Please see -oc help for valid options.\n"));
            }
        }
    }

    // set max matches
    options.maxMatches = static_cast<uint64_t>(numMatchesTmp);

    // set SamBamSeq
    if (samBamSeqDescrTmp == "never")
        options.samBamSeq = 0;
    else if (samBamSeqDescrTmp == "uniq")
        options.samBamSeq = 1;
    else
        options.samBamSeq = 2;

    // help page for SAM/BAM tags
    if (samBamTagsTmp == "help")
    {
        std::cout << "Please specify the tags in this format -oc 'tag1 tag2', i.e. space-seperated and "
                  << "enclosed in quotes. The order of tags is not preserved.\nThe following specifiers are "
                  << "supported:\n";

        for (auto const & c : SamBamExtraTags<>::keyDescPairs)
            std::cout << "\t" << std::get<0>(c) << "\t" << std::get<1>(c) << "\n";

        std::exit(0);
    }
    else
    {
        StringSet<CharString> fields;
        strSplit(fields, samBamTagsTmp, IsSpace(), false);
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
                throw seqan3::parser_invalid_argument(std::string("Unknown column specifier \"") + toCString(str) +
                    std::string("\". Please see \"--sam-bam-tags help\" for valid options.\n"));
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

    // set samBamHardClip
    options.samBamHardClip = (samBamClip == "hard");

    // set options depending on maxSeedDist
    if (options.maxSeedDist == 0)
    {
        if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
        {
            std::cerr << "WARNING: Exact seeeding doesn't benefit from bi-fm-index, so regular index is used.\n";
            options.dbIndexType = DbIndexType::FM_INDEX;
        }
    }

    // Set options based on dbIndexType
    if (options.dbIndexType == DbIndexType::BI_FM_INDEX)
    {
        if (options.seedHalfExact)
            std::cerr << "WARNING: seedHalfExact is already implied by bidirectional indexes.\n";
        else
            options.seedHalfExact = true;
    }

    // Set options depending on extension mode
    if (extensionModeTmp == "fullSIMD")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::FULL_SIMD;
        options.filterPutativeAbundant = false;
        options.filterPutativeDuplicates = false;
        options.mergePutativeSiblings = false;
        options.xDropOff = -1;
    }
    else if (extensionModeTmp == "fullSerial")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::FULL_SERIAL;
        options.filterPutativeAbundant = false;
        options.filterPutativeDuplicates = false;
        options.mergePutativeSiblings = false;
        options.xDropOff = -1;
    }
    else if (extensionModeTmp == "xdrop")
    {
        options.extensionMode = LambdaOptions::ExtensionMode::XDROP;
    }
    else
    {
        options.extensionMode = LambdaOptions::ExtensionMode::AUTO;
    }
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
              << "  minimum % identity:       " << options.idCutOff << "\n"
              << "  maximum e-value:          " << options.eCutOff << "\n"
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
              << "  threads:                  " << uint(options.threads) << "\n"
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
