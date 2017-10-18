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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#include "search.hpp"
#include "mkindex.hpp"

using namespace seqan;

ArgumentParser::ParseResult parseCommandLineMain(int argc, char const ** argv);

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

    ArgumentParser::ParseResult res = parseCommandLineMain(until, argv);

    if (res == ArgumentParser::PARSE_ERROR)
        return ArgumentParser::PARSE_ERROR;
    else if (res != ArgumentParser::PARSE_OK)
        return 0;

    --until; // undo the "+ 1" above

    if ((std::string(argv[until]) == "searchp") || (std::string(argv[until]) == "searchn"))
    {
        return searchMain(argc - until, argv + until);
    }
    else if ((std::string(argv[until]) == "mkindexp") || (std::string(argv[until]) == "mkindexn"))
    {
        return mkindexMain(argc - until, argv + until);
    } else
    {
        // shouldn't be reached
        std::cerr << "WRONG ARGUMENTS!\n";
        return -1;
    }
}

ArgumentParser::ParseResult parseCommandLineMain(int argc, char const ** argv)
{
    ArgumentParser parser("lambda2");

    setShortDescription(parser, "Lambda, the Local Aligner for Massive Biological DataA");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] COMMAND [\\fICOMMAND-OPTIONS\\fP]");
    sharedSetup(parser);

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "COMMAND"));
    setHelpText(parser, 0, "The sub-program to execute. See below.");
    setValidValues(parser, 0, "searchp searchn mkindexp mkindexn");

    addTextSection(parser, "Available commands");
    addText(parser, "\\fBsearchp  \\fP– Perform a protein search (BLASTP, BLASTX, TBLASTN, TBLASTX).");
    addText(parser, "\\fBsearchn  \\fP– Perform a nucleotide search (BLASTN, MEGABLAST).");
    addText(parser, "\\fBmkindexp \\fP– Create an index for protein searches.");
    addText(parser, "\\fBmkindexn \\fP– Create an index for nucleotide searches.");
    addText(parser, "To view the help page for a specific command, simply run 'lambda command --help'.");

    return parse(parser, argc, argv);
}
