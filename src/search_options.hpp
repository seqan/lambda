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

#include <bitset>
#include <cstdio>
#include <filesystem>
#include <unistd.h>

#include <seqan/basic.h>
#include <seqan/blast.h>
#include <seqan/index.h>
#include <seqan/translation.h>

#include <sharg/all.hpp>

#include "search_output.hpp"

// ==========================================================================
// Forwards
// ==========================================================================

// template <typename T>
// struct SamBamExtraTags;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LambdaOptions
// --------------------------------------------------------------------------

struct LambdaOptions : public SharedOptions
{
    std::string queryFile;

    AlphabetEnum qryOrigAlphabet;
    bool         revComp = true;

    int32_t                                     outFileFormat; // -1 = BLAST-Report, 0 = BLAST-Tabular, 1 = SAM, 2 = BAM
    bool                                        blastTabularWithComments = false;
    std::string                                 output                   = "output.m8";
    std::vector<seqan::BlastMatchField<>::Enum> columns;
    std::string                                 outputBam;
    std::bitset<64>                             samBamTags;
    bool                                        samWithRefHeader = false;
    unsigned                                    samBamSeq;
    bool                                        samBamHardClip;
    bool                                        versionInformationToOutputFile = true;
    size_t                                      maximumQueryBlockSize          = 10;

    bool seedHalfExact   = true;
    bool adaptiveSeeding = true;

    struct SearchOpts
    {
        size_t seedLength  = 0;
        size_t maxSeedDist = 1;
        size_t seedOffset  = 0;
    };

    SearchOpts searchOpts0;
    SearchOpts searchOpts;

    // 0 = manual, positive X = blosumX, negative Y = pamY
    int32_t scoringMethod = 62;
    // scores
    int32_t gapOpen       = -11;
    int32_t gapExtend     = -1;
    int32_t match         = 2;  // only for manual
    int32_t misMatch      = -3; // only for manual

    int32_t  minBitScore = 42;
    double   maxEValue   = -1;
    int32_t  idCutOff    = 0;
    uint64_t maxMatches  = 25;

    bool                        computeLCA = false;
    bio::alphabet::genetic_code geneticCodeQry;

    int32_t preScoring       = 2; // 0 = off, 1 = seed, 2 = region
    double  preScoringThresh = 2.0;

    bool        iterativeSearch = true;
    std::string profile         = "none";
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

    sharg::parser parser(programName, argc, argv, sharg::update_notifications::off);

    // Set short description, version, and date.
    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back(
      "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP \\fI-i INDEX.lambda\\fP [\\fI-o output.m8\\fP]");

    sharedSetup(parser);

    // TODO version check
    parser.add_option(options.verbosity,
                      sharg::config{
                        .short_id    = 'v',
                        .long_id     = "verbosity",
                        .description = "Display more/less diagnostic output during operation: "
                                       "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].",

                        .validator = sharg::arithmetic_range_validator{0, 2}
    });

    parser.add_section("Input options");

    // TODO Better solution for file extensions
    parser.add_option(options.queryFile,
                      sharg::config{.short_id    = 'q',
                                    .long_id     = "query",
                                    .description = "Query sequences.",
                                    .required    = true,
                                    .validator   = sharg::input_file_validator{{"fa", "fq", "fasta", "fastq", "gz"}}});

    std::string inputAlphabetTmp = "auto";
    int32_t     geneticCodeTmp   = 1;

    if (options.nucleotide_mode) // seqan::BlastProgram::BLASTN
    {
        options.qryOrigAlphabet = AlphabetEnum::DNA5;
    }
    else
    {
        parser.add_option(inputAlphabetTmp,
                          sharg::config{
                            .short_id    = 'a',
                            .long_id     = "input-alphabet",
                            .description = "Alphabet of the query sequences (specify to override auto-detection). "
                                           "Dna sequences will be translated.",
                            .advanced    = true,
                            .validator   = sharg::value_list_validator{"auto", "dna5", "aminoacid"}
        });
    }

    parser.add_option(options.indexFilePath,
                      sharg::config{.short_id    = 'i',
                                    .long_id     = "index",
                                    .description = std::string{"The database index (created by the 'lambda "} +
                                                   (options.nucleotide_mode ? "mkindexn" : "mkindexp") + "' command).",
                                    .required  = true,
                                    .validator = sharg::input_file_validator{{"lba", "lta"}}});

    parser.add_section("Output options");

    // TODO Fix output file requirements
    parser.add_option(options.output,
                      sharg::config{
                        .short_id = 'o',
                        .long_id  = "output",
                        .description =
                          "File to hold reports on hits (.m* are blastall -m* formats; .m8 is tab-seperated "
                          ".m9 is tab-seperated with with comments, .m0 is pairwise format).",
                        .validator = sharg::output_file_validator{sharg::output_file_open_options::create_new,
                                                                  {"m0", "m8", "m9", "bam", "sam", "gz", "bz2"}}
    });

    std::string outputColumnsTmp = "std";
    parser.add_option(
      outputColumnsTmp,
      sharg::config{.short_id    = '\0',
                    .long_id     = "output-columns",
                    .description = "Print specified column combination and/or order (.m8 and .m9 outputs only); call "
                                   "-oc help for more details.",
                    .advanced    = true});

    parser.add_option(options.idCutOff,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "percent-identity",
                        .description = "Output only matches above this threshold (checked before e-value check).",

                        .validator = sharg::arithmetic_range_validator{0, 100}
    });

    parser.add_option(options.minBitScore,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "bit-score",
                        .description = "Output only matches that score >= this threshold (-1 means no check).",

                        .validator = sharg::arithmetic_range_validator{-1, 1000}
    });

    parser.add_option(options.maxEValue,
                      sharg::config{
                        .short_id    = 'e',
                        .long_id     = "e-value",
                        .description = "Output only matches that score below this threshold (-1 means no check).",

                        .validator = sharg::arithmetic_range_validator{-1, 100}
    });

    int32_t numMatchesTmp = 25;
    parser.add_option(numMatchesTmp,
                      sharg::config{
                        .short_id    = 'n',
                        .long_id     = "num-matches",
                        .description = "Print at most this number of matches per query.",

                        .validator = sharg::arithmetic_range_validator{0, 10000}
    });

    parser.add_option(
      options.samWithRefHeader,
      sharg::config{.short_id    = '\0',
                    .long_id     = "sam-with-refheader",
                    .description = "BAM files require all subject names to be written to the header. For SAM this is "
                                   "not required, so Lambda "
                                   "does "
                                   "not automatically do it to save space (especially for protein database this is a "
                                   "lot!). If you still want "
                                   "them with SAM, e.g. for better BAM compatibility, use this option.",
                    .advanced    = true});

    std::string samBamSeqDescr;
    if (options.nucleotide_mode)
    {
        samBamSeqDescr                  = "Write matching DNA subsequence into SAM/BAM file.";
        options.searchOpts0.seedLength  = 14;
        options.searchOpts0.seedOffset  = 9;
        options.searchOpts0.maxSeedDist = 0;
        options.searchOpts.seedLength   = 14;
        options.searchOpts.seedOffset   = 7;
        options.searchOpts.maxSeedDist  = 1;
        options.preScoringThresh        = 1.4;

        options.gapOpen   = -5;
        options.gapExtend = -2;
    }
    else
    {
        samBamSeqDescr =
          "For BLASTX and TBLASTX the matching protein "
          "sequence is \"untranslated\" and positions retransformed to the original sequence. For BLASTP and "
          "TBLASTN "
          "there is no DNA sequence so a \"*\" is written to the SEQ column. The matching protein sequence "
          "can be "
          "written as an optional tag, see --sam-bam-tags.";

        options.searchOpts0.seedLength  = 10;
        options.searchOpts0.seedOffset  = 5;
        options.searchOpts0.maxSeedDist = 0;
        options.searchOpts.seedLength   = 11;
        options.searchOpts.seedOffset   = 3;
        options.searchOpts.maxSeedDist  = 1;

        options.gapOpen   = -11;
        options.gapExtend = -1;
    }

    std::string samBamSeqDescrTmp = "uniq";
    parser.add_option(samBamSeqDescrTmp,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "sam-bam-seq",
                        .description = samBamSeqDescr +
                                       " If set to uniq than the sequence is "
                                       "omitted iff it is identical to the previous match's subsequence.",
                        .advanced  = true,
                        .validator = sharg::value_list_validator{"always", "uniq", "never"}
    });

    std::string samBamTagsTmp = "AS NM ae ai qf";
    parser.add_option(
      samBamTagsTmp,
      sharg::config{.short_id    = '\0',
                    .long_id     = "sam-bam-tags",
                    .description = "Write the specified optional columns to the SAM/BAM file. Call --sam-bam-tags help "
                                   "for more details.",
                    .advanced    = true});

    std::string samBamClip = "hard";
    parser.add_option(samBamClip,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "sam-bam-clip",
                        .description = "Whether to hard-clip or soft-clip the regions beyond the local match. "
                                       "Soft-clipping retains the full sequence "
                                       "in the output file, but obviously uses more space.",
                        .advanced    = true,
                        .validator   = sharg::value_list_validator{"hard", "soft"}
    });

    parser.add_option(
      options.versionInformationToOutputFile,
      sharg::config{.short_id    = '\0',
                    .long_id     = "version-to-outputfile",
                    .description = "Write the Lambda program tag and version number to the output file.",
                    .hidden      = true});

    parser.add_section("General Options");
#ifdef _OPENMP
    parser.add_option(options.threads,
                      sharg::config{
                        .short_id    = 't',
                        .long_id     = "threads",
                        .description = "Number of threads to run concurrently.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{2, 1000}
    });
#else
    parser.add_option(options.threads,
                      sharg::config{
                        .short_id    = 't',
                        .long_id     = "threads",
                        .description = "LAMBDA BUILT WITHOUT OPENMP; setting this option has no effect.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{2, 2}
    });
#endif

    parser.add_option(options.profile,
                      sharg::config{
                        .short_id    = 'p',
                        .long_id     = "profile",
                        .description = "Profiles are presets of a group of parameters. See below.",

                        .validator = sharg::value_list_validator{"none", "fast", "sensitive"}
    });

    parser.add_section("Seeding / Filtration");

    parser.add_option(options.adaptiveSeeding,
                      sharg::config{.short_id    = '\0',
                                    .long_id     = "adaptive-seeding",
                                    .description = "Grow the seed if it has too many hits (low complexity filter).",
                                    .advanced    = true});

    parser.add_option(options.seedHalfExact,
                      sharg::config{.short_id    = '\0',
                                    .long_id     = "seed-half-exact",
                                    .description = "Allow errors only in second half of seed.",
                                    .advanced    = true});

    parser.add_option(options.searchOpts.seedLength,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-length",
                        .description = "Length of the seeds.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{3, 50}
    });

    parser.add_option(options.searchOpts.seedOffset,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-offset",
                        .description = "Offset for seeding. "
                                       "Distance between seed begin positions.",

                        .validator = sharg::arithmetic_range_validator{1, 50}
    });

    parser.add_option(options.searchOpts.maxSeedDist,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-delta",
                        .description = "Maximum seed distance.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{0, 3}
    });

    /* iterative search parameters */
    parser.add_option(options.iterativeSearch,
                      sharg::config{.short_id    = '\0',
                                    .long_id     = "search0",
                                    .description = "If (cheaper) pre-search yield results, skip regular search.",
                                    .advanced    = true});

    parser.add_option(options.searchOpts0.seedLength,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-length0",
                        .description = "Length of the seeds.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{3, 50}
    });

    parser.add_option(options.searchOpts0.seedOffset,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-offset0",
                        .description = "Offset for seeding. "
                                       "Distance between seed begin positions.",

                        .validator = sharg::arithmetic_range_validator{1, 50}
    });

    parser.add_option(options.searchOpts0.maxSeedDist,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "seed-delta0",
                        .description = "Maximum seed distance.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{0, 5}
    });

    parser.add_section("Miscellaneous Heuristics");

    parser.add_option(options.preScoring,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "pre-scoring",
                        .description = "Evaluate score of a region NUM times the size of the seed "
                                       "before extension (0 -> no pre-scoring, 1 -> evaluate seed, n-> area "
                                       "around seed, as well; default = 1 if no reduction is used).",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{1, 10}
    });

    parser.add_option(options.preScoringThresh,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "pre-scoring-threshold",
                        .description = "Minimum average score per position in pre-scoring region.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{0, 20}
    });

    parser.add_section("Scoring");

    if (options.nucleotide_mode)
    {
        parser.add_option(options.match,
                          sharg::config{
                            .short_id    = '\0',
                            .long_id     = "score-match",
                            .description = "Match score",
                            .advanced    = true,
                            .validator   = sharg::arithmetic_range_validator{-1000, 1000}
        });

        parser.add_option(options.misMatch,
                          sharg::config{
                            .short_id    = '\0',
                            .long_id     = "score-mismatch",
                            .description = "Mismatch score",
                            .advanced    = true,
                            .validator   = sharg::arithmetic_range_validator{-1000, 1000}
        });
    }
    else
    {
        parser.add_option(options.scoringMethod,
                          sharg::config{
                            .short_id    = 's',
                            .long_id     = "scoring-scheme",
                            .description = "Use '45' for Blosum45; '62' for Blosum62 (default); '80' for Blosum80.",
                            .advanced    = true,
                            .validator   = sharg::value_list_validator{45, 62, 80}
        });
    }

    parser.add_option(options.gapExtend,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "score-gap",
                        .description = "Score per gap character.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{-1000, 1000}
    });

    parser.add_option(options.gapOpen,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "score-gap-open",
                        .description = "Additional cost for opening gap.",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{-1000, 1000}
    });

    parser.add_section("Profiles");
    parser.add_line(
      "Setting a profile other than \"none\" always overwrites manually given command line "
      "arguments!",
      true);

    if (options.nucleotide_mode)
    {
        parser.add_line("\"fast\"", false);
        parser.add_line("--seed-length 14 --seed-offset 9 --seed-delta 0", true);
        parser.add_line("\"sensitive\"", false);
        parser.add_line("--seed-length0 14 --seed-offset0 3 --seed-length 14 --seed-offset 3", true);
    }
    else
    {
        parser.add_line("\"fast\"", false);
        parser.add_line("--seed-length0 12 --seed-offset0 8 --seed-length 10 --seed-offset 5 --seed-delta 0", true);
        parser.add_line("\"sensitive\"", false);
        parser.add_line(
          "--seed-length0 9 --seed-offset0 4 --seed-length 8 --seed-offset 3 --pre-scoring 3 "
          "--pre-scoring-threshold 1.9",
          true);
    }
    parser.add_line("For further information see the wiki: <https://github.com/seqan/lambda/wiki>", false);

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
            throw sharg::parser_error("ERROR: Invalid argument to --input-alphabet\n");

        options.geneticCodeQry = static_cast<bio::alphabet::genetic_code>(geneticCodeTmp);
    }

    if (options.profile == "fast")
    {
        if (options.nucleotide_mode)
        {
            options.iterativeSearch        = false;
            options.searchOpts.seedOffset  = 9;
            options.searchOpts.maxSeedDist = 0;
        }
        else
        {
            options.searchOpts0.seedLength = 12;
            options.searchOpts0.seedOffset = 8;
            options.searchOpts.seedLength  = 10;
            options.searchOpts.seedOffset  = 5;
            options.searchOpts.maxSeedDist = 0;
        }
    }
    else if (options.profile == "sensitive")
    {
        if (options.nucleotide_mode)
        {
            options.searchOpts0.seedOffset = 3;
            options.searchOpts.seedOffset  = 3;
        }
        else
        {
            options.searchOpts0.seedLength = 9;
            options.searchOpts0.seedOffset = 4;
            options.searchOpts.seedLength  = 8;
            options.searchOpts.seedOffset  = 3;

            options.preScoring       = 3;
            options.preScoringThresh = 1.9;
        }
    }

    // set output file format
    std::string outputPath = options.output;
    if (std::filesystem::path(outputPath).extension() == ".gz")
        outputPath.resize(seqan::length(outputPath) - 3);
    else if (std::filesystem::path(outputPath).extension() == ".bz2")
        outputPath.resize(seqan::length(outputPath) - 4);

    if (std::filesystem::path(outputPath).extension() == ".sam")
        options.outFileFormat = 1;
    else if (std::filesystem::path(outputPath).extension() == ".bam")
        options.outFileFormat = 2;
    else if (std::filesystem::path(outputPath).extension() == ".m0")
        options.outFileFormat = -1;
    else if (std::filesystem::path(outputPath).extension() == ".m8")
    {
        options.outFileFormat            = 0;
        options.blastTabularWithComments = false;
    }
    else if (std::filesystem::path(outputPath).extension() == ".m9")
    {
        options.outFileFormat            = 0;
        options.blastTabularWithComments = true;
    }
    else
    {
        throw 99;
    }

    // help page for output columns
    if (outputColumnsTmp == "help")
    {
        std::cout << "Please specify the columns in this format -oc 'column1 column2', i.e. "
                     "space-separated "
                     "and "
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
                    else if ((static_cast<seqan::BlastMatchField<>::Enum>(i) ==
                              seqan::BlastMatchField<>::Enum::LCA_ID) ||
                             (static_cast<seqan::BlastMatchField<>::Enum>(i) ==
                              seqan::BlastMatchField<>::Enum::LCA_TAX_ID))
                        options.computeLCA = true;
                    break;
                }
            }
            if (!resolved)
            {
                throw sharg::parser_error(std::string("Unknown column specifier \"") + toCString(str) +
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
                  << "enclosed in quotes. The order of tags is not preserved.\nThe following specifiers "
                     "are "
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
                    resolved              = true;
                    break;
                }
            }
            if (!resolved)
            {
                std::cerr << "Unknown column specifier \"" << str
                          << "\". Please see \"--sam-bam-tags help\" for valid options.\n";
                throw sharg::parser_error(std::string("Unknown column specifier \"") + seqan::toCString(str) +
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
inline void printOptions(LambdaOptions const & options)
{
    std::cout << "OPTIONS\n"
              << " INPUT\n"
              << "  query file:               " << options.queryFile << "\n"
              << "  index file:               " << options.indexFilePath << "\n"
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
              << "  threads:                  " << uint(options.threads) << "\n"
              << " TRANSLATION AND ALPHABETS\n"
              << "  genetic code:             " << (int)options.geneticCodeQry << "\n"
              << "  nucleotide_mode:          " << options.nucleotide_mode << "\n"
              << "  original alphabet (query):" << _alphabetEnumToName(options.qryOrigAlphabet) << "\n"
              << " SEEDING\n"
              << "  seed length:              " << options.searchOpts.seedLength << "\n"
              << "  seed offset:              " << options.searchOpts.seedOffset << "\n"
              << "  seed delta:               " << options.searchOpts.maxSeedDist << "\n"
              << "  adaptive seeding:         " << (options.adaptiveSeeding ? std::string("on") : std::string("off"))
              << "\n"
              << "  pre-search:               " << (options.iterativeSearch ? std::string("on") : std::string("off"))
              << "\n"
              << "  seed length0:             " << options.searchOpts0.seedLength << "\n"
              << "  seed offset0:             " << options.searchOpts0.seedOffset << "\n"
              << "  seed delta0:              " << options.searchOpts0.maxSeedDist << "\n"
              << " MISCELLANEOUS HEURISTICS\n"
              << "  pre-scoring:              " << (options.preScoring ? std::string("on") : std::string("off")) << "\n"
              << "  pre-scoring-region:       "
              << (options.preScoring ? std::to_string(options.preScoring * options.searchOpts.seedLength)
                                     : std::string("n/a"))
              << "\n"
              << "  pre-scoring-threshold:    "
              << (options.preScoring ? std::to_string(options.preScoringThresh) : std::string("n/a")) << "\n"
              << " SCORING\n"
              << "  scoring scheme:           " << options.scoringMethod << "\n"
              << "  score-match:              "
              << (options.scoringMethod ? std::string("n/a") : std::to_string(options.match)) << "\n"
              << "  score-mismatch:           "
              << (options.scoringMethod ? std::string("n/a") : std::to_string(options.misMatch)) << "\n"
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
