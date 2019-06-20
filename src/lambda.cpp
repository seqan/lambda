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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#include <seqan3/argument_parser/all.hpp>

#include "search.hpp"

#if 1
#include "mkindex.hpp"
#endif
// using namespace seqan;

void parseCommandLineMain(int argc, char const ** argv);

int main(int argc, char const ** argv)
{
    if (std::string(CMAKE_BUILD_TYPE) != "Release")
        std::cerr << "WARNING: This binary is not built in release mode and will be much slower than it should be!\n";

    int until = argc;
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
            } else
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

    // TODO change return values
    if ((std::string(argv[until]) == "searchp") || (std::string(argv[until]) == "searchn"))
    {

        try
        {
            searchMain(argc - until, argv + until);
        }
        catch (std::exception const & ext)
        {
            std::cerr << ext.what() << "\n";
            return -1;
        }


    }
    else

    if ((std::string(argv[until]) == "mkindexp") || (std::string(argv[until]) == "mkindexn"))
    {
    #if 0
        try
        {
            mkindexMain(argc - until, argv + until);
        }
        catch (std::exception const & ext)
        {
            std::cerr << ext.what() << "\n";
            return -1;
        }
    #endif
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
    seqan3::argument_parser parser("lambda3", argc, argv);

    parser.info.short_description = "Lambda, the Local Aligner for Massive Biological DatA.";
    parser.info.synopsis.push_back("[\\fIOPTIONS\\fP] COMMAND [\\fICOMMAND-OPTIONS\\fP]");

    sharedSetup(parser);

    std::string command{};
    parser.add_positional_option(command, "The sub-program to execute. See below.",
        seqan3::value_list_validator({"searchp", "searchn", "mkindexp", "mkindexn"}));

    parser.info.description.push_back("Available commands");
    parser.info.description.push_back("\\fBsearchp  \\fP– Perform a protein search (BLASTP, BLASTX, TBLASTN, TBLASTX).");
    parser.info.description.push_back("\\fBsearchn  \\fP– Perform a nucleotide search (BLASTN, MEGABLAST).");
    parser.info.description.push_back("\\fBmkindexp \\fP– Create an index for protein searches.");
    parser.info.description.push_back("\\fBmkindexn \\fP– Create an index for nucleotide searches.");
    parser.info.description.push_back("To view the help page for a specific command, simply run 'lambda command --help'.");

    parser.parse();
}
