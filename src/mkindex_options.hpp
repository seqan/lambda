// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2017, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2017, Knut Reinert and Freie Universität Berlin
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

#include <seqan/basic.h>
#include <seqan/translation.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/blast.h>

using namespace seqan;

// --------------------------------------------------------------------------
// Class LambdaIndexerOptions
// --------------------------------------------------------------------------

struct LambdaIndexerOptions : public SharedOptions
{
    std::string     dbFile;
//     std::string     segFile = "";
    std::string     algo = "";
    std::string     accToTaxMapFile;
    std::string     taxDumpDir;

    bool            truncateIDs;

    int alphReduction;

    LambdaIndexerOptions()
        : SharedOptions()
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// INDEXER
ArgumentParser::ParseResult
parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    std::string programName = "lambda2 " + std::string(argv[0]);

    // this is important for option handling:
    if (std::string(argv[0]) == "mkindexn")
        options.blastProgram = BlastProgram::BLASTN;

    ArgumentParser parser(programName);

    setShortDescription(parser, "the Local Aligner for Massive Biological DatA");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-d DATABASE.fasta [-i INDEX.lambda]\\fP");

    sharedSetup(parser);

    // version checker initiated by top-level arg parser
    setDefaultValue(parser, "version-check", "0");
    hideOption(parser, "version-check");

    addOption(parser, ArgParseOption("v", "verbosity",
        "Display more/less diagnostic output during operation: 0 [only errors]; 1 [default]; 2 "
        "[+run-time, options and statistics].",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbosity", "1");
    setMinValue(parser, "verbosity", "0");
    setMaxValue(parser, "verbosity", "2");

    addDescription(parser, "This is the indexer_binary for creating lambda-compatible databases.");

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("d", "database",
        "Database sequences.",
        ArgParseArgument::INPUT_FILE,
        "IN"));
    setRequired(parser, "database");
    setValidValues(parser, "database", toCString(concat(getFileExtensions(SeqFileIn()), ' ')));

    addOption(parser, ArgParseOption("tx",
        "acc-tax-map",
        "An NCBI or UniProt accession-to-taxid mapping file. Download from "
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ or "
        "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ .",
        ArgParseArgument::INPUT_FILE));

    CharString taxExtensions = "accession2taxid dat";
#ifdef SEQAN_HAS_ZLIB
    taxExtensions+= " accession2taxid.gz";
    taxExtensions+= " accession2taxid.bgzf";
    taxExtensions+= " dat.gz";
    taxExtensions+= " dat.bgzf";
#endif
#ifdef SEQAN_HAS_BZIP2
    taxExtensions+= " accession2taxid.bz2";
    taxExtensions+= " dat.bz2";
#endif
    setValidValues(parser, "acc-tax-map", toCString(taxExtensions));

    addOption(parser, ArgParseOption("tt",
        "tax-dump-dir",
        "A directory that contains nodes.dmp and names.dmp; unzipped from "
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        ArgParseArgument::INPUT_DIRECTORY));

    addSection(parser, "Output Options");
    addOption(parser, ArgParseOption("i", "index",
        "The output directory for the index files (defaults to \"DATABASE.lambda\").",
        ArgParseArgument::OUTPUT_DIRECTORY,
        "OUT"));
    setValidValues(parser, "index", ".lambda");

    addOption(parser, ArgParseOption("di", "db-index-type",
        "Suffix array or full-text minute space.",
        ArgParseArgument::STRING,
        "type"));
#ifdef LAMBDA_LEGACY_PATHS
    setValidValues(parser, "db-index-type", "sa fm bifm");
#else
    setValidValues(parser, "db-index-type", "fm bifm");
#endif
    setDefaultValue(parser, "db-index-type", "fm");
    setAdvanced(parser, "db-index-type");

    addOption(parser, ArgParseOption("", "truncate-ids",
        "Truncate IDs at first whitespace. This saves a lot of space and is irrelevant for all LAMBDA output formats "
        "other than BLAST Pairwise (.m0).",
        ArgParseArgument::BOOL));
    setDefaultValue(parser, "truncate-ids", "on");

    if (options.blastProgram != BlastProgram::BLASTN)
    {
        addSection(parser, "Alphabets and Translation");

        addOption(parser, ArgParseOption("", "database-alphabet",
            "Alphabet of the database sequences (specify to override auto-detection);"
            "if input is Dna, it will be translated.",
            ArgParseArgument::STRING));
        setValidValues(parser, "database-alphabet", "auto dna5 aminoacid");
        setDefaultValue(parser, "database-alphabet", "auto");
        setAdvanced(parser, "database-alphabet");

        addOption(parser, ArgParseOption("g", "genetic-code",
            "The translation table to use if input is Dna. See "
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
            " for ids (default is generic).",
            ArgParseArgument::INTEGER));
        setDefaultValue(parser, "genetic-code", "1");
        setAdvanced(parser, "genetic-code");

        addOption(parser, ArgParseOption("ar", "alphabet-reduction",
            "Alphabet Reduction for seeding phase.",
            ArgParseArgument::STRING,
            "STR"));
        setValidValues(parser, "alphabet-reduction", "none murphy10");
        setDefaultValue(parser, "alphabet-reduction", "murphy10");
        setAdvanced(parser, "alphabet-reduction");
    }

    addSection(parser, "Algorithm");
    addOption(parser, ArgParseOption("a", "algorithm",
        "Algorithm for SA construction (also used for FM; see Memory "
        " Requirements below!).",
        ArgParseArgument::STRING,
        "STR"));
    setValidValues(parser, "algorithm", "mergesort quicksortbuckets quicksort radixsort skew7ext");
    setDefaultValue(parser, "algorithm", "radixsort");
    setAdvanced(parser, "algorithm");

#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads",
        "number of threads to run concurrently (ignored if a == skew7ext).",
        ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", std::to_string(omp_get_max_threads()));
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

    std::string tmpdir;
    getCwd(tmpdir);
    addOption(parser, ArgParseOption("td", "tmp-dir",
        "temporary directory used by skew, defaults to working directory.",
        ArgParseArgument::OUTPUT_DIRECTORY,
        "STR"));
    setDefaultValue(parser, "tmp-dir", tmpdir);
    setAdvanced(parser, "tmp-dir");

    //TODO move manual / auto-detect
//     addTextSection(parser, "Memory requirements and Speed");
//     addText(parser, "\033[1mmergesort [RAM]:\033[0m"
//                     "\t14 * size(dbSeqs)");
//     addText(parser, "\033[1mmergesort [speed]:\033[0m"
//                     "\tup to t threads");
//     addText(parser, "\033[1mquicksort and quicksortbuckets [RAM]:\033[0m"
//                     "\t7 * size(dbSeqs)");
//     addText(parser, "\033[1mquicksort [speed]:\033[0m"
//                     "\t1-2 threads");
//     addText(parser, "\033[1mquicksortbuckets [speed]:\033[0m"
//                     "\t1-2 threads for initial sort, up to t for buckets");
//     addText(parser, "\033[1mskew7ext [RAM]:\033[0m"
//                     "\t2 * size(dbSeqs)");
//     addText(parser, "\033[1mskew7ext [DISK]:\033[0m"
//                     "\t30 * size(dbSeqs)");
//     addText(parser, "\033[1mskew7ext [speed]:\033[0m"
//                     "\tnot parallelized");
//     addText(parser, "size(dbSeqs) refers to the total "
//                     "sequence length and does not include IDs (which can "
//                     "account for >50% of the file size for protein databases). "
//                     "The space is the maximum obseverved factor, for many "
//                     "databases the factor is smaller." );
//     addText(parser, "Use mergesort if you have enough memory! If not, you will "
//                     "probably want to use skew. For small databases and only a "
//                     "few cores the quicksorts might be a good tradeoff. "
//                     "mergesort and quicksortbuckets provide a rough progress "
//                     "estimate.");
// //     addText(parser, "Disk space required is in TMPDIR which you can set as "
// //                     "an environment variable.");

    addTextSection(parser, "Remarks");
    addText(parser, "Please see the wiki (<https://github.com/seqan/lambda/wiki>) for more information on which indexes"
        " to chose and which algorithms to pick.");

    addText(parser, "Note that the indexes created are binary and not compatible between different CPU endiannesses. "
        "Also the on-disk format is still subject to change between Lambda versions.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Options shared by lambda and its indexer
    res = parseCommandLineShared(options, parser);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    std::string buffer;
    int buf = 0;
    getOptionValue(buffer, parser, "db-index-type");
    if (buffer == "sa")
        options.dbIndexType = DbIndexType::SUFFIX_ARRAY;
    else if (buffer == "bifm")
        options.dbIndexType = DbIndexType::BI_FM_INDEX;
    else
        options.dbIndexType = DbIndexType::FM_INDEX;


    if (options.blastProgram == BlastProgram::BLASTN)
    {
        options.subjOrigAlphabet = AlphabetEnum::DNA5;
        options.transAlphabet    = AlphabetEnum::DNA5;
        options.reducedAlphabet  = AlphabetEnum::DNA5;
    }
    else
    {
        getOptionValue(buffer, parser, "database-alphabet");
        if (buffer == "auto")
            options.subjOrigAlphabet = AlphabetEnum::DNA4;
        else if (buffer == "dna5")
            options.subjOrigAlphabet = AlphabetEnum::DNA5;
        else if (buffer == "aminoacid")
            options.subjOrigAlphabet = AlphabetEnum::AMINO_ACID;
        else
            throw std::invalid_argument("ERROR: Invalid argument to --database-alphabet\n");


        getOptionValue(buffer, parser, "alphabet-reduction");
        if (buffer == "murphy10")
        {
            options.reducedAlphabet = AlphabetEnum::MURPHY10;
            //TODO deprecate:
            options.alphReduction = 2;
        }
        else
        {
            options.reducedAlphabet = AlphabetEnum::AMINO_ACID;
            options.alphReduction = 0;
        }

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

    }
//     getOptionValue(options.segFile, parser, "segfile");
    getOptionValue(options.algo, parser, "algorithm");
    if ((options.algo == "mergesort") || (options.algo == "quicksort") || (options.algo == "quicksortbuckets"))
    {
        std::cerr << "WARNING: " << options.algo << " tag is deprecated and superseded by \"radixsort\", please "
                  << "adapt your program calls.\n";
        options.algo = "radixsort";
    }

    getOptionValue(tmpdir, parser, "tmp-dir");
    setEnv("TMPDIR", tmpdir);

    getOptionValue(options.truncateIDs, parser, "truncate-ids");

    getOptionValue(options.dbFile, parser, "database");
    if (isSet(parser, "index"))
        getOptionValue(options.indexDir, parser, "index");
    else
        options.indexDir = options.dbFile + ".lambda";


    if (fileExists(options.indexDir.c_str()))
    {
        std::cerr << "ERROR: An output directory already exists at " << options.indexDir << '\n'
                  << "Remove it, or choose a different location.\n";
        return ArgumentParser::PARSE_ERROR;
    }
    else
    {
        if (mkdir(options.indexDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
        {
            std::cerr << "ERROR: Cannot create output directory at " << options.indexDir << '\n';;
            return ArgumentParser::PARSE_ERROR;
        }
    }

    getOptionValue(options.accToTaxMapFile, parser, "acc-tax-map");
    options.hasSTaxIds = (options.accToTaxMapFile != "");

    getOptionValue(options.taxDumpDir, parser, "tax-dump-dir");
    if (!options.taxDumpDir.empty())
    {
        if (!options.hasSTaxIds)
        {
            std::cerr << "ERROR: There is no point in inclduing a taxonomic tree in the index, if\n"
                         "       you don't also include taxonomic IDs for your sequences.\n";
            return ArgumentParser::PARSE_ERROR;
        }

        //TODO check existance of directory
    }

    return ArgumentParser::PARSE_OK;
}

#endif // header guard
