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
// mkindex_options.h: contains the options and argument parser for the indexer
// ==========================================================================

#pragma once

#include <bitset>
#include <cstdio>
#include <unistd.h>

#include <sharg/all.hpp>

#include "shared_options.hpp"

// --------------------------------------------------------------------------
// Class LambdaIndexerOptions
// --------------------------------------------------------------------------

struct LambdaIndexerOptions : public SharedOptions
{
    std::string dbFile;
    std::string algo = "default";
    std::string accToTaxMapFile;
    std::string taxDumpDir;

    std::string tmpdir = std::filesystem::current_path();

    bool truncateIDs = true;

    int alphReduction;

    LambdaIndexerOptions() : SharedOptions() {}
};

// ==========================================================================
// Functions
// ==========================================================================

// INDEXER
void parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    std::string const subcommand  = std::string(argv[0]);
    std::string const programName = "lambda3-" + subcommand;

    // this is important for option handling:
    if (subcommand == "mkindexp")
        options.domain = domain_t::protein;
    else if (subcommand == "mkindexn")
        options.domain = domain_t::nucleotide;
    else if (subcommand == "mkindexbs")
        options.domain = domain_t::bisulfite;
    else
        throw std::runtime_error{"Unknown subcommand."};

    sharg::parser parser(programName, argc, argv, sharg::update_notifications::off);

    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back("lambda3 "s + subcommand +
                                   " [\\fIOPTIONS\\fP] \\-d DATABASE.fasta [-i INDEX.lba]\\fP");

    parser.info.description.push_back("This is the indexer command for creating lambda-compatible databases.");

    sharedSetup(parser);

    parser.add_option(options.verbosity,
                      sharg::config{
                        .short_id    = 'v',
                        .long_id     = "verbosity",
                        .description = "Display more/less diagnostic output during operation: "
                                       "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].",
                        .validator   = sharg::arithmetic_range_validator{0, 2}
    });

    parser.add_section("Input Options");

    // TODO Change file extensions, make more generic
    parser.add_option(options.dbFile,
                      sharg::config{.short_id    = 'd',
                                    .long_id     = "database",
                                    .description = "Database sequences.",
                                    .required    = true,
                                    .validator   = sharg::input_file_validator{{"fa", "fq", "fasta", "fastq", "gz"}}});

    std::vector<std::string> taxExtensions{"accession2taxid", "dat"};
#ifdef SEQAN_HAS_ZLIB
    taxExtensions.push_back("accession2taxid.gz");
    taxExtensions.push_back("accession2taxid.bgzf");
    taxExtensions.push_back("dat.gz");
    taxExtensions.push_back("dat.bgzf");
#endif
#ifdef SEQAN_HAS_BZIP2
    taxExtensions.push_back("accession2taxid.bz2");
    taxExtensions.push_back("dat.bz2");
#endif

    parser.add_option(
      options.accToTaxMapFile,
      sharg::config{.short_id = 'm',
                    .long_id  = "acc-tax-map",
                    .description =
                      "An NCBI or UniProt accession-to-taxid mapping file. Download from "
                      "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ or "
                      "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ .",
                    .validator = sharg::input_file_validator(taxExtensions)});

    parser.add_option(options.taxDumpDir,
                      sharg::config{.short_id    = 'x',
                                    .long_id     = "tax-dump-dir",
                                    .description = "A directory that contains nodes.dmp and names.dmp; unzipped from "
                                                   "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                                    .validator   = sharg::input_directory_validator()});

    parser.add_section("Output Options");

    options.indexFilePath = "»INPUT«.lba";
    parser.add_option(options.indexFilePath,
                      sharg::config{
                        .short_id    = 'i',
                        .long_id     = "index",
                        .description = "The output path for the index file.",
                        .validator   = sharg::output_file_validator{sharg::output_file_open_options::create_new,
                                                                    {"lba", "lta", "lba.gz", "lta.gz"}}
    });

    options.threads = std::max<size_t>(2ul, std::min<size_t>(std::thread::hardware_concurrency(), 4ul));
    parser.add_option(options.threads,
                      sharg::config{
                        .short_id    = 't',
                        .long_id     = "threads",
                        .description = "Number of threads (only used for compression).",
                        .advanced    = true,
                        .validator   = sharg::arithmetic_range_validator{2, 1000}
    });

    std::string dbIndexTypeTmp = "fm";
    parser.add_option(dbIndexTypeTmp,
                      sharg::config{
                        .short_id    = '\0',
                        .long_id     = "db-index-type",
                        .description = "FM-Index oder bidirectional FM-Index.",
                        .advanced    = true,
                        .validator   = sharg::value_list_validator{"fm", "bifm"}
    });

    parser.add_option(
      options.truncateIDs,
      sharg::config{.short_id = '\0',
                    .long_id  = "truncate-ids",
                    .description =
                      "Truncate IDs at first whitespace. This saves a lot of space and is irrelevant for all "
                      "LAMBDA output formats other than BLAST Pairwise (.m0)."});

    std::string inputAlphabetTmp = "auto";
    std::string alphabetReductionTmp;
    int         geneticCodeTmp = 1;

    switch (options.domain)
    {
        case domain_t::protein:
            alphabetReductionTmp               = "li10";
            options.indexFileOptions.origAlph  = AlphabetEnum::UNDEFINED;
            options.indexFileOptions.transAlph = AlphabetEnum::AMINO_ACID;
            options.indexFileOptions.redAlph   = AlphabetEnum::LI10;

            parser.add_section("Alphabet and Translation");

            parser.add_option(inputAlphabetTmp,
                              sharg::config{
                                .short_id    = 'a',
                                .long_id     = "input-alphabet",
                                .description = "Alphabet of the database sequences (specify to override "
                                               "auto-detection); if input is Dna, it will be translated.",
                                .advanced    = true,
                                .validator   = sharg::value_list_validator{"auto", "dna5", "aminoacid"}
            });

            parser.add_option(alphabetReductionTmp,
                              sharg::config{
                                .short_id    = 'r',
                                .long_id     = "alphabet-reduction",
                                .description = "Alphabet Reduction for seeding phase.",
                                .advanced    = true,
                                .validator   = sharg::value_list_validator{"none", "murphy10", "li10"}
            });
            break;
        case domain_t::nucleotide:
            options.indexFileOptions.origAlph  = AlphabetEnum::DNA5;
            options.indexFileOptions.transAlph = AlphabetEnum::DNA5;
            options.indexFileOptions.redAlph   = AlphabetEnum::DNA4;
            break;
        case domain_t::bisulfite:
            options.indexFileOptions.origAlph  = AlphabetEnum::DNA5;
            options.indexFileOptions.transAlph = AlphabetEnum::DNA5;
            options.indexFileOptions.redAlph   = AlphabetEnum::DNA3BS;
            break;
    }

    parser.add_section("Remarks");
    parser.add_line("Indexes are *NOT* compatible between different CPU endiannesses.", false);

    // parse command line.
    parser.parse();

    // set db index type
    if (dbIndexTypeTmp == "bifm")
        options.indexFileOptions.indexType = DbIndexType::BI_FM_INDEX;
    else
        options.indexFileOptions.indexType = DbIndexType::FM_INDEX;

    // set options for protein alphabet, genetic code and alphabet reduction
    if (options.domain == domain_t::protein)
    {
        options.indexFileOptions.origAlph = _alphabetNameToEnum(inputAlphabetTmp);
        if (alphabetReductionTmp == "none")
            options.indexFileOptions.redAlph = AlphabetEnum::AMINO_ACID;
        else
            options.indexFileOptions.redAlph = _alphabetNameToEnum(alphabetReductionTmp);
        options.indexFileOptions.geneticCode = static_cast<bio::alphabet::genetic_code>(geneticCodeTmp);
    }

    setEnv("TMPDIR", options.tmpdir);

    // set hasSTaxIds based on taxonomy file
    options.hasSTaxIds = (options.accToTaxMapFile != "");

    if (options.indexFilePath == "»INPUT«.lba")
        options.indexFilePath = options.dbFile + ".lba";

    if (std::filesystem::exists(options.indexFilePath))
    {
        throw sharg::parser_error("ERROR: An output file already exists at " + options.indexFilePath.string() +
                                  "\n       Remove it, or choose a different location.\n");
    }

    if (!options.taxDumpDir.empty())
    {
        if (!options.hasSTaxIds)
        {
            throw sharg::parser_error(
              "ERROR: There is no point in including a taxonomic tree in the index, if\n"
              "       you don't also include taxonomic IDs for your sequences.\n");
        }
    }
}
