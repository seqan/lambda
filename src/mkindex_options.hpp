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


#ifndef SEQAN_LAMBDA_OPTIONS_H_
#define SEQAN_LAMBDA_OPTIONS_H_

#include <cstdio>
#include <unistd.h>
#include <bitset>

// --------------------------------------------------------------------------
// Class LambdaIndexerOptions
// --------------------------------------------------------------------------

struct LambdaIndexerOptions : public SharedOptions
{
    std::string     dbFile;
//     std::string     segFile = "";
    std::string     algo = "radixsort";
    std::string     accToTaxMapFile;
    std::string     taxDumpDir;

    std::string     tmpdir = std::filesystem::current_path();

    bool            truncateIDs = true;

    int alphReduction;

    LambdaIndexerOptions()
        : SharedOptions()
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// INDEXER
void parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    std::string programName = "lambda3 " + std::string(argv[0]);

    // this is important for option handling:
    options.nucleotide_mode = (std::string(argv[0]) == "mkindexn");

    seqan3::argument_parser parser(programName, argc, argv);

    parser.info.short_description = "the Local Aligner for Massive Biological DatA";

    // Define usage line and long description.
    parser.info.synopsis.push_back("[\\fIOPTIONS\\fP] \\-d DATABASE.fasta [-i INDEX.lambda]\\fP");

    parser.info.description.push_back("This is the indexer command for creating lambda-compatible databases.");


    sharedSetup(parser);

    // TODO add version check

    parser.add_option(options.verbosity, 'v', "verbosity", "Display more/less diagnostic output during operation: "
        "0 [only errors]; 1 [default]; 2 [+run-time, options and statistics].", seqan3::option_spec::DEFAULT,
        seqan3::arithmetic_range_validator{0, 2});

    parser.add_section("Input Options");

    // TODO Change file extensions, make more generic
    parser.add_option(options.dbFile, 'd', "database", "Database sequences.", seqan3::option_spec::REQUIRED,
        seqan3::path_existence_validator() | seqan3::file_ext_validator({"fa", "fq", "fasta", "fastq"}));

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

    parser.add_option(options.accToTaxMapFile, 'm', "acc-tax-map",
        "An NCBI or UniProt accession-to-taxid mapping file. Download from "
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ or "
        "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ .",
        seqan3::option_spec::DEFAULT, seqan3::file_ext_validator(taxExtensions));

    parser.add_option(options.taxDumpDir,'x', "tax-dump-dir",
        "A directory that contains nodes.dmp and names.dmp; unzipped from "
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz");

    parser.add_section("Output Options");

    // TODO Does this input directory structure work?
    parser.add_option(options.indexDir, 'i', "index",
        "The output directory for the index files (defaults to \"DATABASE.lambda\").",
        seqan3::option_spec::DEFAULT);

    std::string dbIndexTypeTmp = "fm";
    parser.add_option(dbIndexTypeTmp, '\0', "db-index-type", "Suffix array or full-text minute space.",
        seqan3::option_spec::ADVANCED);

    parser.add_option(options.truncateIDs, '\0', "truncate-ids",
        "Truncate IDs at first whitespace. This saves a lot of space and is irrelevant for all LAMBDA output formats "
        "other than BLAST Pairwise (.m0).");

    std::string inputAlphabetTmp = "auto";
    std::string alphabetReductionTmp = "murphy10";
    int geneticCodeTmp = 1;

    if (options.nucleotide_mode)
    {
        options.subjOrigAlphabet = AlphabetEnum::DNA5;
        options.transAlphabet    = AlphabetEnum::DNA5;
        options.reducedAlphabet  = AlphabetEnum::DNA5;
    }
    else
    {
        parser.add_section("Alphabet and Translation");

        parser.add_option(inputAlphabetTmp,'a', "input-alphabet",
            "Alphabet of the database sequences (specify to override auto-detection); "
            "if input is Dna, it will be translated.",
            seqan3::option_spec::ADVANCED, seqan3::value_list_validator({"auto", "dna5", "aminoacid"}));

        parser.add_option(geneticCodeTmp, 'g', "genetic-code",
            "The translation table to use if input is Dna. See "
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
            " for ids. (default is generic)", seqan3::option_spec::ADVANCED,
            seqan3::value_list_validator({0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25}));

        parser.add_option(alphabetReductionTmp, 'r', "alphabet-reduction", "Alphabet Reduction for seeding phase.",
            seqan3::option_spec::ADVANCED, seqan3::value_list_validator({"none", "murphy10"}));
    }

    std::string algorithmTmp = "radixsort";
    parser.add_option(algorithmTmp, '\0', "algorithm",
        "Algorithm for SA construction (also used for FM; see Memory Requirements below!).",
        seqan3::option_spec::ADVANCED,
        seqan3::value_list_validator({"mergesort", "quicksortbuckets", "quicksort", "radixsort", "skew7ext"}));

    parser.add_option(options.threads, 't', "threads",
        "Number of threads to run concurrently (ignored if a == skew7ext).", seqan3::option_spec::ADVANCED,
        seqan3::arithmetic_range_validator{1, static_cast<double>(options.threads)});

    // TODO change validator
    parser.add_option(options.tmpdir,'\0', "tmp-dir", "temporary directory used by skew, defaults to working directory.",
        seqan3::option_spec::ADVANCED, seqan3::path_existence_validator());

    parser.add_section("Remarks");
    parser.add_line("Please see the wiki (<https://github.com/seqan/lambda/wiki>) for more information on which indexes"
        " to chose and which algorithms to pick.", false);
    parser.add_line("Note that the indexes created are binary and not compatible between different CPU endiannesses. "
        "Also the on-disk format is still subject to change between Lambda versions.", false);

    // parse command line.
    parser.parse();

    // set db index type
    if (dbIndexTypeTmp == "sa")
        options.dbIndexType = DbIndexType::SUFFIX_ARRAY;
    else if (dbIndexTypeTmp == "bifm")
        options.dbIndexType = DbIndexType::BI_FM_INDEX;
    else
        options.dbIndexType = DbIndexType::FM_INDEX;

    // set options for protein alphabet, genetic code and alphabet reduction
    if (!options.nucleotide_mode)
    {
        if (inputAlphabetTmp == "auto")
            options.subjOrigAlphabet = AlphabetEnum::DNA4;
        else if (inputAlphabetTmp == "dna5")
            options.subjOrigAlphabet = AlphabetEnum::DNA5;
        else if (inputAlphabetTmp == "aminoacid")
            options.subjOrigAlphabet = AlphabetEnum::AMINO_ACID;
        else
            throw seqan3::parser_invalid_argument("ERROR: Invalid argument to --input-alphabet\n");

        if (alphabetReductionTmp == "murphy10")
        {
            options.reducedAlphabet = AlphabetEnum::MURPHY10;
            // TODO deprecate:
            options.alphReduction = 2;
        }
        else
        {
            options.reducedAlphabet = AlphabetEnum::AMINO_ACID;
            options.alphReduction = 0;
        }

        options.geneticCode = static_cast<seqan3::genetic_code>(geneticCodeTmp);
    }

    // set algorithm option
    if ((algorithmTmp == "mergesort") || (algorithmTmp == "quicksort") || (algorithmTmp == "quicksortbuckets"))
    {
        std::cerr << "WARNING: " << algorithmTmp << " tag is deprecated and superseded by \"radixsort\", please "
                  << "adapt your program calls.\n";
        algorithmTmp = "radixsort";
    }
    options.algo = algorithmTmp;

    setEnv("TMPDIR", options.tmpdir);

    // set hasSTaxIds based on taxonomy file
    options.hasSTaxIds = (options.accToTaxMapFile != "");

    if (options.indexDir == "")
        options.indexDir = options.dbFile + ".lambda";

    if (std::filesystem::exists(options.indexDir))
    {
        throw seqan3::parser_invalid_argument("ERROR: An output directory already exists at " + options.indexDir +
            "Remove it, or choose a different location.\n");
    }
    else
    {
        if (mkdir(options.indexDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
        {
            throw seqan3::parser_invalid_argument("ERROR: Cannot create output directory at " + options.indexDir + '\n');
        }
    }

    if (!options.taxDumpDir.empty())
    {
        if (!options.hasSTaxIds)
        {
            throw seqan3::parser_invalid_argument("ERROR: There is no point in inclduing a taxonomic tree in the index, if\n"
                                                  "       you don't also include taxonomic IDs for your sequences.\n");
        }
        //TODO check existance of directory
    }
}

#endif // header guard
