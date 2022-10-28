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

#include <seqan3/argument_parser/all.hpp>

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
    std::string programName = "lambda3-" + std::string(argv[0]);

    // this is important for option handling:
    options.nucleotide_mode = (std::string(argv[0]) == "mkindexn");

    seqan3::argument_parser parser(programName, argc, argv, seqan3::update_notifications::off);

    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back("[\\fIOPTIONS\\fP] \\-d DATABASE.fasta [-i INDEX.lba]\\fP");

    parser.info.description.push_back("This is the indexer command for creating lambda-compatible databases.");

    sharedSetup(parser);

    parser.add_option(options.verbosity,
                      'v',
                      "verbosity",
                      "Display more/less diagnostic output during operation: "
                      "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 2});

    parser.add_section("Input Options");

    // TODO Change file extensions, make more generic
    parser.add_option(options.dbFile,
                      'd',
                      "database",
                      "Database sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{
                        {"fa", "fq", "fasta", "fastq", "gz"}
    });

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

    parser.add_option(options.accToTaxMapFile,
                      'm',
                      "acc-tax-map",
                      "An NCBI or UniProt accession-to-taxid mapping file. Download from "
                      "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ or "
                      "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ .",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator(taxExtensions));

    parser.add_option(options.taxDumpDir,
                      'x',
                      "tax-dump-dir",
                      "A directory that contains nodes.dmp and names.dmp; unzipped from "
                      "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                      seqan3::option_spec::standard,
                      seqan3::input_directory_validator());

    parser.add_section("Output Options");

    options.indexFilePath = "»INPUT«.lba";
    parser.add_option(options.indexFilePath,
                      'i',
                      "index",
                      "The output path for the index file.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{
                        seqan3::output_file_open_options::create_new,
                        {"lba", "lta"}
    });

    std::string dbIndexTypeTmp = "fm";
    parser.add_option(dbIndexTypeTmp,
                      '\0',
                      "db-index-type",
                      "FM-Index oder bidirectional FM-Index.",
                      seqan3::option_spec::advanced,
                      seqan3::value_list_validator{"fm", "bifm"});

    parser.add_option(
      options.truncateIDs,
      '\0',
      "truncate-ids",
      "Truncate IDs at first whitespace. This saves a lot of space and is irrelevant for all LAMBDA output formats "
      "other than BLAST Pairwise (.m0).");

    std::string inputAlphabetTmp = "auto";
    std::string alphabetReductionTmp;
    int         geneticCodeTmp = 1;

    if (options.nucleotide_mode)
    {
        alphabetReductionTmp               = "dna4";
        options.indexFileOptions.origAlph  = AlphabetEnum::DNA5;
        options.indexFileOptions.transAlph = AlphabetEnum::DNA5;
        options.indexFileOptions.redAlph   = AlphabetEnum::DNA4;

        parser.add_section("Alphabet reduction");

        parser.add_option(alphabetReductionTmp,
                          'r',
                          "alphabet-reduction",
                          "Alphabet Reduction for seeding phase.",
                          seqan3::option_spec::advanced,
                          seqan3::value_list_validator{"none", "dna4", "dna3bs"});
    }
    else
    {
        alphabetReductionTmp               = "li10";
        options.indexFileOptions.origAlph  = AlphabetEnum::UNDEFINED;
        options.indexFileOptions.transAlph = AlphabetEnum::AMINO_ACID;
        options.indexFileOptions.redAlph   = AlphabetEnum::LI10;

        parser.add_section("Alphabet and Translation");

        parser.add_option(inputAlphabetTmp,
                          'a',
                          "input-alphabet",
                          "Alphabet of the database sequences (specify to override auto-detection); "
                          "if input is Dna, it will be translated.",
                          seqan3::option_spec::advanced,
                          seqan3::value_list_validator{"auto", "dna5", "aminoacid"});

        parser.add_option(alphabetReductionTmp,
                          'r',
                          "alphabet-reduction",
                          "Alphabet Reduction for seeding phase.",
                          seqan3::option_spec::advanced,
                          seqan3::value_list_validator{"none", "murphy10", "li10"});
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
    if (!options.nucleotide_mode)
    {
        options.indexFileOptions.origAlph = _alphabetNameToEnum(inputAlphabetTmp);
        if (alphabetReductionTmp == "none")
            options.indexFileOptions.redAlph = AlphabetEnum::AMINO_ACID;
        else
            options.indexFileOptions.redAlph = _alphabetNameToEnum(alphabetReductionTmp);
        options.indexFileOptions.geneticCode = static_cast<seqan3::genetic_code>(geneticCodeTmp);
    }
    else
    {
        if (alphabetReductionTmp == "none")
            options.indexFileOptions.redAlph = AlphabetEnum::DNA5;
        else
            options.indexFileOptions.redAlph = _alphabetNameToEnum(alphabetReductionTmp);
    }

    setEnv("TMPDIR", options.tmpdir);

    // set hasSTaxIds based on taxonomy file
    options.hasSTaxIds = (options.accToTaxMapFile != "");

    if (options.indexFilePath == "»INPUT«.lba")
        options.indexFilePath = options.dbFile + ".lba";

    if (std::filesystem::exists(options.indexFilePath))
    {
        throw seqan3::argument_parser_error("ERROR: An output file already exists at " +
                                            options.indexFilePath.string() +
                                            "\n       Remove it, or choose a different location.\n");
    }

    if (!options.taxDumpDir.empty())
    {
        if (!options.hasSTaxIds)
        {
            throw seqan3::argument_parser_error(
              "ERROR: There is no point in including a taxonomic tree in the index, if\n"
              "       you don't also include taxonomic IDs for your sequences.\n");
        }
        //TODO check existance of directory
    }
}
