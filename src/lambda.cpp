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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#include <sharg/all.hpp>

#include "mkindex.hpp"
#include "search.hpp"
#include "shared_options.hpp"

void parseCommandLineMain(int argc, char const ** argv);

int main(int argc, char const ** argv)
{
    if (std::string(CMAKE_BUILD_TYPE) != "Release")
        std::cerr << "WARNING: This binary is not built in release mode and will be much slower than it should be!\n";

    int  until    = argc;
    bool skipNext = false;
    for (int i = 1; i < argc; ++i)
    {
        // version check expects a parameter
        if (std::string(argv[i]) == "--version-check")
            skipNext = true;

        if (argv[i][0] != '-')
        {
            if (skipNext)
            {
                skipNext = false;
            }
            else
            {
                until = i + 1;
                break;
            }
        }
    }

    try
    {
        parseCommandLineMain(until, argv);
    }
    catch (std::exception const & ext)
    {
        std::cerr << ext.what() << "\n";
        return -1;
    }

    --until; // undo the "+ 1" above

    std::string const subcommand_actual = argv[until];

    if (subcommand_actual.starts_with("search"))
    {
        return searchMain(argc - until, argv + until);
    }
    else if (subcommand_actual.starts_with("mkindex"))
    {
        return mkindexMain(argc - until, argv + until);
    }
    else
    {
        // shouldn't be reached
        std::cerr << "WRONG ARGUMENTS!\n";
        return -1;
    }
    return 0;
}

void parseCommandLineMain(int argc, char const ** argv)
{
    sharg::parser parser("lambda3", argc, argv, sharg::update_notifications::off);

    parser.info.short_description = "Lambda, the Local Aligner for Massive Biological DatA.";
    parser.info.synopsis.push_back("lambda3 [\\fIOPTIONS\\fP] COMMAND [\\fICOMMAND-OPTIONS\\fP]");

    sharedSetup(parser);

    std::string command{};
    parser.add_positional_option(
      command,
      sharg::config{
        .description = "The sub-program to execute. See above.",
        .validator =
          sharg::value_list_validator{"searchp", "searchn", "searchbs", "mkindexp", "mkindexn", "mkindexbs"}
    });

    parser.info.description.push_back("Available commands");
    parser.info.description.push_back(
      "\\fBsearchp  \\fP– Perform a protein search (BLASTP, BLASTX, TBLASTN, TBLASTX).");
    parser.info.description.push_back("\\fBsearchn  \\fP– Perform a nucleotide search (BLASTN, MEGABLAST).");
    parser.info.description.push_back("\\fBsearchbs \\fP– Perform a bisulfite search.");
    parser.info.description.push_back("\\fBmkindexp \\fP– Create an index for protein searches.");
    parser.info.description.push_back("\\fBmkindexn \\fP– Create an index for nucleotide searches.");
    parser.info.description.push_back("\\fBmkindexbs\\fP– Create an index for bisulfite searches.");
    parser.info.description.push_back(
      "To view the help page for a specific command, simply run 'lambda command --help'.");

    parser.parse();
}
