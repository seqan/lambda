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
// search_options.h: contains the options and argument parser for the search
// ==========================================================================

#pragma once

#include <cstdio>
#include <unistd.h>
#include <bitset>

#include <seqan/basic.h>
#include <seqan/translation.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/blast.h>

#include <seqan3/argument_parser/all.hpp>
#include <filesystem>

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

    int32_t         outFileFormat; // -1 = BLAST-Report, 0 = BLAST-Tabular, 1 = SAM, 2 = BAM
    bool            blastTabularWithComments = false;
    std::string     output = "output.m8";
    std::vector<seqan::BlastMatchField<>::Enum> columns;
    std::string     outputBam;
    std::bitset<64> samBamTags;
    bool            samWithRefHeader = false;
    unsigned        samBamSeq;
    bool            samBamHardClip;
    bool            versionInformationToOutputFile = true;
    size_t          maximumQueryBlockSize = 10;

    bool            seedHalfExact = true;
    bool            adaptiveSeeding = true;

    unsigned        seedLength  = 0;
    unsigned        maxSeedDist = 1;

    unsigned        seedOffset      = 0;

    // 0 = manual, positive X = blosumX, negative Y = pamY
    int32_t         scoringMethod   = 62;
    // scores
    int32_t         gapOpen         = -11;
    int32_t         gapExtend       = -1;
    int32_t         match           = 2; // only for manual
    int32_t         misMatch        = -3; // only for manual

    int32_t         band        = -3;
    double          eCutOff     = 1e-04;
    int32_t         idCutOff    = 0;
    uint64_t        maxMatches  = 256;

    bool            computeLCA  = false;
    seqan3::genetic_code geneticCodeQry;

    int32_t         preScoring = 2; // 0 = off, 1 = seed, 2 = region
    double          preScoringThresh    = 2.0;
};

void parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // save commandLine
    for (int i = 0; i < argc; ++i)
        options.commandLine += std::string(argv[i]) + " ";
    seqan::eraseBack(options.commandLine);

    std::string programName = "lambda3-" + std::string(argv[0]);

    // this is important for option handling:
    options.nucleotide_mode = (std::string(argv[0]) == "searchn");

    seqan3::argument_parser parser(programName, argc, argv, seqan3::update_notifications::off);

    // Set short description, version, and date.
    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back("[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP \\fI-i INDEX.lambda\\fP [\\fI-o output.m8\\fP]");

    sharedSetup(parser);

    // TODO version check

    parser.add_option(options.verbosity, 'v', "verbosity", "Display more/less diagnostic output during operation: "
        "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].",
        seqan3::option_spec::standard, seqan3::arithmetic_range_validator{0, 2});

    parser.add_section("Input options");

    // TODO Better solution for file extensions
    parser.add_option(options.queryFile, 'q', "query", "Query sequences.", seqan3::option_spec::required,
        seqan3::input_file_validator{{"fa", "fq", "fasta", "fastq", "gz"}});

    std::string inputAlphabetTmp = "auto";
    int32_t geneticCodeTmp = 1;

    if (options.nucleotide_mode) // seqan::BlastProgram::BLASTN
    {
        options.qryOrigAlphabet = AlphabetEnum::DNA5;
    }
    else
    {
        parser.add_option(inputAlphabetTmp, 'a', "input-alphabet",
            "Alphabet of the query sequences (specify to override auto-detection). Dna sequences will be translated.",
            seqan3::option_spec::advanced, seqan3::value_list_validator{"auto", "dna5", "aminoacid"});

#if 0 // TODO: we currently don't support other code, but we should
        parser.add_option(geneticCodeTmp, 'g', "genetic-code",
            "The translation table to use if input is Dna. See "
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
            " for ids. Default is to use the same table that was used for the index or 1/CANONICAL if the index "
            "was not translated.", seqan3::option_spec::advanced,
            seqan3::value_list_validator{0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25});
#endif
    }

    parser.add_option(options.indexFilePath,
                      'i',
                      "index",
                      std::string{"The database index (created by the 'lambda "} +
                        (options.nucleotide_mode ? "mkindexn" : "mkindexp") +
                        "' command).",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"lba", "lta"}});

    parser.add_section("Output options");

    // TODO Fix outout file requirements
    parser.add_option(options.output,
                      'o',
                      "output",
                      "File to hold reports on hits (.m* are blastall -m* formats; .m8 is tab-seperated "
                      ".m9 is tab-seperated with with comments, .m0 is pairwise format).",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                    {"m0", "m8", "m9", "bam", "sam", "gz", "bz2"}});

    std::string outputColumnsTmp = "std";
    parser.add_option(outputColumnsTmp, '\0', "output-columns",
        "Print specified column combination and/or order (.m8 and .m9 outputs only); call -oc help for more details.",
        seqan3::option_spec::advanced);

    parser.add_option(options.idCutOff, '\0', "percent-identity",
        "Output only matches above this threshold (checked before e-value check).", seqan3::option_spec::standard,
        seqan3::arithmetic_range_validator{0, 100});

    parser.add_option(options.eCutOff, 'e', "e-value", "Output only matches that score below this threshold.",
        seqan3::option_spec::standard, seqan3::arithmetic_range_validator{0, 100});

    int32_t numMatchesTmp = 256;
    parser.add_option(numMatchesTmp, 'n', "num-matches", "Print at most this number of matches per query.",
        seqan3::option_spec::standard, seqan3::arithmetic_range_validator{0, 10000});

    parser.add_option(options.samWithRefHeader, '\0', "sam-with-refheader",
        "BAM files require all subject names to be written to the header. For SAM this is not required, so Lambda does "
        "not automatically do it to save space (especially for protein database this is a lot!). If you still want "
        "them with SAM, e.g. for better BAM compatibility, use this option.", seqan3::option_spec::advanced);

    std::string samBamSeqDescr;
    if (options.nucleotide_mode)
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
        "omitted iff it is identical to the previous match's subsequence.", seqan3::option_spec::advanced,
        seqan3::value_list_validator{"always", "uniq", "never"});

    std::string samBamTagsTmp = "AS NM ae ai qf";
    parser.add_option(samBamTagsTmp, '\0', "sam-bam-tags",
        "Write the specified optional columns to the SAM/BAM file. Call --sam-bam-tags help for more details.",
        seqan3::option_spec::advanced);

    std::string samBamClip = "hard";
    parser.add_option(samBamClip, '\0', "sam-bam-clip",
        "Whether to hard-clip or soft-clip the regions beyond the local match. Soft-clipping retains the full sequence "
        "in the output file, but obviously uses more space.", seqan3::option_spec::advanced,
        seqan3::value_list_validator{"hard", "soft"});

    parser.add_option(options.versionInformationToOutputFile, '\0', "version-to-outputfile",
        "Write the Lambda program tag and version number to the output file.", seqan3::option_spec::hidden);

    parser.add_section("General Options");

#ifdef _OPENMP
    parser.add_option(options.threads, 't', "threads", "Number of threads to run concurrently.",
        seqan3::option_spec::advanced,
        seqan3::arithmetic_range_validator{2, 1000});
#else
    parser.add_option(options.threads, 't', "threads",
        "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.", seqan3::option_spec::advanced,
        seqan3::arithmetic_range_validator{2, 2});
#endif

    parser.add_section("Seeding / Filtration");

    parser.add_option(options.adaptiveSeeding, '\0', "adaptive-seeding",
        "Grow the seed if it has too many hits (low complexity filter).", seqan3::option_spec::advanced);

   parser.add_option(options.seedHalfExact, '\0', "seed-half-exact",
        "Allow errors only in second half of seed.", seqan3::option_spec::advanced);

    unsigned defaultSeedLength = options.nucleotide_mode ? 14 : 10;

    options.seedLength = defaultSeedLength;
    parser.add_option(options.seedLength, '\0', "seed-length", "Length of the seeds.", seqan3::option_spec::advanced,
        seqan3::arithmetic_range_validator{3, 50});

    options.seedOffset = options.seedLength / 2;
    parser.add_option(options.seedOffset, '\0', "seed-offset", "Offset for seeding. "
        "If you set 'seed-length', please consider setting this option to half of that.",
        seqan3::option_spec::standard, seqan3::arithmetic_range_validator{1, 50});

    parser.add_option(options.maxSeedDist, '\0', "seed-delta",
        "Maximum seed distance.", seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{0, 5});

    parser.add_section("Miscellaneous Heuristics");

    parser.add_option(options.preScoring, '\0', "pre-scoring",
        "Evaluate score of a region NUM times the size of the seed "
        "before extension (0 -> no pre-scoring, 1 -> evaluate seed, n-> area "
        "around seed, as well; default = 1 if no reduction is used).",
        seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{1, 10});

    parser.add_option(options.preScoringThresh, '\0', "pre-scoring-threshold",
        "Minimum average score per position in pre-scoring region.", seqan3::option_spec::advanced,
        seqan3::arithmetic_range_validator{0, 20});

    parser.add_section("Scoring");

    if (options.nucleotide_mode)
    {
        parser.add_option(options.match, '\0', "score-match",
            "Match score [only BLASTN]", seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{-1000, 1000});

        parser.add_option(options.misMatch, '\0', "score-mismatch",
            "Mismatch score [only BLASTN]", seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{-1000, 1000});
    }
    else
    {
        parser.add_option(options.scoringMethod,'s', "scoring-scheme",
            "Use '45' for Blosum45; '62' for Blosum62 (default); '80' for Blosum80.", seqan3::option_spec::advanced,
            seqan3::value_list_validator{45, 62, 80});
    }

    if (options.nucleotide_mode)
        options.gapExtend = -2;
    else
        options.gapExtend = -1;

    parser.add_option(options.gapExtend, '\0', "score-gap",
        "Score per gap character.", seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{-1000, 1000});

    if (options.nucleotide_mode)
        options.gapOpen = -5;
    else
        options.gapOpen = -11;

    parser.add_option(options.gapOpen, '\0', "score-gap-open",
        "Additional cost for opening gap.", seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{-1000, 1000});

    parser.add_section("Extension");

    //TODO this is only used in serial mode right now
    parser.add_option(options.band, 'b', "band", "Size of the DP-band used in extension (-3 means log2 of query length;"
        " -2 means sqrt of query length; -1 means full dp; n means band of size 2n+1)",
        seqan3::option_spec::advanced, seqan3::arithmetic_range_validator{-3, 1000});

#if 0 //TODO make new guide
    parser.add_section("Tuning");
    parser.add_line("Tuning the seeding parameters and (de)activating alphabet "
                    "reduction has a strong "
                    "influence on both speed and sensitivity. We recommend the "
                    "following alternative profiles for protein searches:", false);
    parser.add_line("fast (high similarity):    --seed-delta-increases-length on", false);
    parser.add_line("sensitive (lower similarity): --seed-offset 3", false);
    parser.add_line("For further information see the wiki: <https://github.com/seqan/lambda/wiki>", false);
#endif

    // parse command line.
    parser.parse();

    // set query alphabet and genetic code depending on options
    if (!options.nucleotide_mode)
    {
        if (inputAlphabetTmp == "auto")
            options.qryOrigAlphabet = AlphabetEnum::DNA4;
        else if (inputAlphabetTmp == "dna5")
            options.qryOrigAlphabet = AlphabetEnum::DNA5;
        else if (inputAlphabetTmp == "aminoacid")
            options.qryOrigAlphabet = AlphabetEnum::AMINO_ACID;
        else
            throw seqan3::argument_parser_error("ERROR: Invalid argument to --input-alphabet\n");

        options.geneticCodeQry = static_cast<seqan3::genetic_code>(geneticCodeTmp);
    }

    // set output file format
    std::string outputPath = options.output;
    if (std::filesystem::path(outputPath).extension() == ".gz")
        outputPath.resize(seqan::length(outputPath) - 3);
    else if (std::filesystem::path(outputPath).extension() ==  ".bz2")
        outputPath.resize(seqan::length(outputPath) - 4);

    if (std::filesystem::path(outputPath).extension() == ".sam")
        options.outFileFormat = 1;
    else if (std::filesystem::path(outputPath).extension() == ".bam")
        options.outFileFormat = 2;
    else if (std::filesystem::path(outputPath).extension() == ".m0")
        options.outFileFormat = -1;
    else if (std::filesystem::path(outputPath).extension() == ".m8")
    {
        options.outFileFormat = 0;
        options.blastTabularWithComments = false;
    } else if (std::filesystem::path(outputPath).extension() == ".m9")
    {
        options.outFileFormat = 0;
        options.blastTabularWithComments = true;
    } else
    {
        throw 99;
    }

    // help page for output columns
    if (outputColumnsTmp == "help")
    {
        std::cout << "Please specify the columns in this format -oc 'column1 column2', i.e. space-separated and "
                  << "enclosed in single quotes.\nThe specifiers are the same as in NCBI Blast, currently "
                  << "the following are supported:\n";
        for (unsigned i = 0; i < seqan::length(seqan::BlastMatchField<>::implemented); ++i)
        {
            if (seqan::BlastMatchField<>::implemented[i])
            {
                std::cout << "\t" << seqan::BlastMatchField<>::optionLabels[i]
                          << (seqan::length(seqan::BlastMatchField<>::optionLabels[i]) >= 8 ? "\t" : "\t\t")
                          << seqan::BlastMatchField<>::descriptions[i] << "\n";
            }
        }
        std::exit(0);
    }
    else
    {
        seqan::StringSet<seqan::CharString> fields;
        seqan::strSplit(fields, outputColumnsTmp, seqan::IsSpace(), false);
        for (auto str : fields)
        {
            bool resolved = false;
            for (unsigned i = 0; i < seqan::length(seqan::BlastMatchField<>::optionLabels); ++i)
            {
                if (seqan::BlastMatchField<>::optionLabels[i] == str)
                {
                    seqan::appendValue(options.columns, static_cast<seqan::BlastMatchField<>::Enum>(i));
                    resolved = true;
                    if (static_cast<seqan::BlastMatchField<>::Enum>(i) == seqan::BlastMatchField<>::Enum::S_TAX_IDS)
                        options.hasSTaxIds = true;
                    else if ((static_cast<seqan::BlastMatchField<>::Enum>(i) == seqan::BlastMatchField<>::Enum::LCA_ID) ||
                             (static_cast<seqan::BlastMatchField<>::Enum>(i) == seqan::BlastMatchField<>::Enum::LCA_TAX_ID))
                        options.computeLCA = true;
                    break;
                }
            }
            if (!resolved)
            {
                throw seqan3::argument_parser_error(std::string("Unknown column specifier \"") + toCString(str) +
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
        std::cout << "Please specify the tags in this format -oc 'tag1 tag2', i.e. space-separated and "
                  << "enclosed in quotes. The order of tags is not preserved.\nThe following specifiers are "
                  << "supported:\n";

        for (auto const & c : SamBamExtraTags<>::keyDescPairs)
            std::cout << "\t" << std::get<0>(c) << "\t" << std::get<1>(c) << "\n";

        std::exit(0);
    }
    else
    {
        seqan::StringSet<seqan::CharString> fields;
        seqan::strSplit(fields, samBamTagsTmp, seqan::IsSpace(), false);
        for (auto str : fields)
        {
            bool resolved = false;
            for (unsigned i = 0; i < seqan::length(SamBamExtraTags<>::keyDescPairs); ++i)
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
                throw seqan3::argument_parser_error(std::string("Unknown column specifier \"") + seqan::toCString(str) +
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
}

// --------------------------------------------------------------------------
// Function printOptions()
// --------------------------------------------------------------------------

template <typename TLH>
inline void
printOptions(LambdaOptions const & options)
{
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
              << "  index file:               " << options.indexFilePath << "\n"
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
              << "  genetic code:             " << (int)options.geneticCodeQry << "\n"
              << "  nucleotide_mode:          " << options.nucleotide_mode << "\n"
              << "  original alphabet (query):" << _alphabetEnumToName(options.qryOrigAlphabet) << "\n"
              << " SEEDING\n"
              << "  seed length:              " << uint(options.seedLength) << "\n"
              << "  seed offset:              " << uint(options.seedOffset) << "\n"
              << "  seed delta:               " << uint(options.maxSeedDist) << "\n"
              << "  adaptive seeding:         " << (options.adaptiveSeeding
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
              << " BUILD OPTIONS:\n"
              << "  cmake_build_type:         " << std::string(CMAKE_BUILD_TYPE) << "\n"
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
              << "  seqan_simd:               "
    #if defined(__AVX2__)
              << "avx2\n"
    #elif defined(__SSE4_2__)
              << "sse4\n"
    #else
              << "off\n"
    #endif
              << "\n";
}
