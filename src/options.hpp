// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// options.h: contains the options and argument parser
// ==========================================================================


#ifndef SEQAN_LAMBDA_OPTIONS_H_
#define SEQAN_LAMBDA_OPTIONS_H_

#include <seqan/basic.h>
#include <seqan/blast.h>
#include <seqan/arg_parse.h>

#include <seqan/blast.h>

using namespace seqan;

// ==========================================================================
// Metafunctions
// ==========================================================================

// using TCDSpec = Packed<>;
using TCDSpec = Alloc<>;

template <typename TAlph>
using TCDStringSet = StringSet<String<TAlph, TCDSpec>, Owner<ConcatDirect<> > >;


template <BlastFormatOptions::Program p>
struct UnreducedStringSet
{
    typedef TCDStringSet<AminoAcid> Type;
};

template <>
struct UnreducedStringSet<BlastFormatOptions::BlastN>
{
    typedef TCDStringSet<Dna5> Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LambdaOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.


struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    CharString  dbFile;

    // for indexer, the file format of database sequences
    // for main app, the file format of query sequences
    // 0 -- fasta, 1 -- fastq
    int         fileFormat; 

    int         alphReduction;

    BlastFormatOptions::Program blastProg;

    SharedOptions() :
        verbosity(2),fileFormat(0),alphReduction(0),
        blastProg(BlastFormatOptions::BlastX)
    {}
};


struct LambdaOptions : public SharedOptions
{

    CharString      queryFile;
    bool            revComp     = true;

    CharString      output;

    unsigned short  queryPart = 0;

//     bool            semiGlobal;

    unsigned char   seedLength  = 0;
    unsigned char   maxSeedDist = 1;
    bool            hammingOnly = false;

    char            seedGravity     = 0;
    unsigned char   minSeedLength   = 0;
    unsigned char   minSeedScore    = 0;
    unsigned int    minSeedEVal     = 0;
    double          minSeedBitS     = -1;

    // 0 = manual, positive X = blosumX, negative Y = pamY
    int             scoringMethod   = 62;
    // scores
    int             gapOpen         = -11;
    int             gapExtend       = -1;
    int             match           = 0; // only for manual
    int             misMatch        = 0; // only for manual

    int             xDropOff    = 0;
    int             band        = -1;
    double          eCutOff     = 0;





    LambdaOptions() :
        SharedOptions()
    {}
};

struct LambdaIndexerOptions : public SharedOptions
{
    CharString      segFile;

    LambdaIndexerOptions()
        : SharedOptions(), segFile("")
    {}
};



// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(LambdaOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("lambda");
    // Set short description, version, and date.
    setShortDescription(parser, "BLAST compatible local aligner optimized for "
                                "NGS and Metagenomics.");
    setVersion(parser, "0.2");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI-q QUERY.fasta\\fP "
                         "\\fI-d DATABASE.fasta\\fP "
                         "[\\fI-o output.m8\\fP]");
    addDescription(parser, "Lambda is a local aligner that is faster than "
                           "BLAST, optimized for many query sequences.");

    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("q",
                                            "query",
                                            "Query sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "q");
    setValidValues(parser, "query", "fasta fa fna faa fastq fq");



    addOption(parser, seqan::ArgParseOption("d",
                                            "database",
                                            "Database sequences (fasta), with "
                                             "precomputed index (.sa).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "d");
    setValidValues(parser, "database", "fasta fa fna faa");

    // TODO implement some way to transfer options from index builder
    // maybe write a spec file that is loaded by this app

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "File to hold reports on hits (.m8 "
                                            "is blastall -m8 et cetera)",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "m0 m8 m9");
    setDefaultValue(parser, "output", "output.m8");


    addSection(parser, "Program Options");
    addOption(parser, seqan::ArgParseOption("p",
                                            "program",
                                            "Blast Operation Mode.",
                                            seqan::ArgParseArgument::STRING,
                                            "OUT"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "Alphabet Reduction for AminoAcid "
                                            "Alphabet (0 -> off; 2 -> Murphy10; "
                                            "8,10,12 -> Lambda*).",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "alph", "2");

    addOption(parser, seqan::ArgParseOption("qp",
                                            "querypart",
                                            "Partition query into N blocks;"
                                            "defaults to number of CPUs; "
                                            "increase if you run out of memory "
                                            "to mulitple of NCPU",
                                            seqan::ArgParseArgument::INTEGER));
//     setValidValues(parser, "alph", "0 10");
    setDefaultValue(parser, "querypart", "0");

    addSection(parser, "Seeding / Filtration");
    addOption(parser, seqan::ArgParseOption("su",
                                            "ungappedseeds",
                                            "Deacitvate EditDistance-Seeding.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "ungappedseeds", "1");

    addOption(parser, seqan::ArgParseOption("sl",
                                            "seedlength",
                                            "Length of the seeds (0 -> choose "
                                            "automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seedlength", "0");

    addOption(parser, seqan::ArgParseOption("sd",
                                            "seeddelta",
                                            "maximum seed distance.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seeddelta", "1");

    addOption(parser, seqan::ArgParseOption("sg",
                                            "seedgravity",
                                            "Seeds closer than this are joined"
                                            "(-1 -> choose automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seedgravity", "-1");

    addOption(parser, seqan::ArgParseOption("sm",
                                            "seedminlength",
                                            "after postproc shorter seeds are "
                                            "discarded"
                                            "(0 -> choose automatically).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seedminlength", "0");

    addOption(parser, seqan::ArgParseOption("ss",
                                            "seedminscore",
                                            "after postproc worse seeds are "
                                            "discarded.",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "seedminscore", "32");

//     addOption(parser, seqan::ArgParseOption("se",
//                                             "seedminevalue",
//                                             "after postproc worse seeds are "
//                                             "discarded"
//                                             "(0 -> off).",
//                                             seqan::ArgParseArgument::INTEGER));
//     setDefaultValue(parser, "seedminevalue", "100000");

        addOption(parser, seqan::ArgParseOption("sb",
                                            "seedminbits",
                                            "after postproc worse seeds are "
                                            "discarded"
                                            "(-1 -> off).",
                                            seqan::ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "seedminbits", "-1");

    addSection(parser, "Extension");

    addOption(parser, seqan::ArgParseOption("x",
                                            "xdrop",
                                            "Stop Banded extension if score "
                                            "x below the maximum seen"
                                            "(-1 means no xdrop).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "xdrop", "30");

    addOption(parser, seqan::ArgParseOption("b",
                                            "band",
                                            "Size of the DP-band used in "
                                            "extension (-3 means log2 "
                                            "of query length; -2 means sqrt of "
                                            "query length; -1 means full dp; "
                                            "n means band of size 2n+1)",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "band", "-3");

    addOption(parser, seqan::ArgParseOption("e",
                                            "evalue",
                                            "Minimum E-Value for Results.",
                                            seqan::ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "evalue", "0.1");


    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");
    addListItem(parser, "\\fBOMP_NUM_THREADS\\fP",
                        "number of threads to use.");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    getOptionValue(options.dbFile, parser, "database");

    CharString blastMode;
    getOptionValue(blastMode, parser, "program");
    if (blastMode == "blastn")
    {
        options.blastProg   = BlastFormatOptions::BlastN;
        options.seedLength  = 20;
        options.alphReduction = 0;
    }
    else if (blastMode == "blastp")
    {
        options.blastProg = BlastFormatOptions::BlastP;
        options.seedLength  = 10;
    }
    else if (blastMode == "blastx")
    {
        options.blastProg = BlastFormatOptions::BlastX;
        options.seedLength  = 10;
    }
    else if (blastMode == "tblastn")
    {
        options.blastProg = BlastFormatOptions::TBlastN;
        options.seedLength  = 10;
    }
    else if (blastMode == "tblastx")
    {
        options.blastProg = BlastFormatOptions::TBlastX;
        options.seedLength  = 10;
    }
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    //TODO adapt seedlength if alphReduction

    int buf = 0;
    getOptionValue(buf, parser, "seedlength");
    if (buf != 0) // seed length was specified manually
        options.seedLength = buf;

    buf = 0;
    getOptionValue(buf, parser, "seeddelta");
    options.maxSeedDist = buf;
    if (buf == 0)
    {
        options.hammingOnly = true;
    } else
    {
        getOptionValue(buf, parser, "ungappedseeds");
        options.hammingOnly = (buf != 0);
    }


    getOptionValue(buf, parser, "seedgravity");
    if (buf == -1)
        options.seedGravity = options.seedLength;
    else
        options.seedGravity = buf;

    unsigned foo = 0;
    getOptionValue(foo, parser, "seedminlength");
    options.minSeedLength = foo;
    if (options.minSeedLength == 0)
        options.minSeedLength = options.seedLength;

    getOptionValue(foo, parser, "seedminscore");
    options.minSeedScore = foo;

//     getOptionValue(foo, parser, "seedminevalue");
//     options.minSeedEVal = foo;

    getOptionValue(options.minSeedBitS, parser, "seedminbits");

    getOptionValue(options.eCutOff, parser, "evalue");

    getOptionValue(options.xDropOff, parser, "xdrop");
    if (options.xDropOff < -1)
    {
        std::cerr << "xDropOff must be in [ -1 : MAXINT ]" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.band, parser, "band");
    if (options.band < -3)
    {
        std::cerr << "band must be in [ -3 : MAXINT ]" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    //verifyFileFormat(); TODO verify fileformats for program mode

    getOptionValue(options.output, parser, "output");

    getOptionValue(options.alphReduction, parser, "alph");
    switch (options.alphReduction)
    {
        case 1: case 2: case 8: case 10: case 12:
        if (options.blastProg == BlastFormatOptions::BlastN)
        {
            std::cerr << "Alphabet Reduction only makes sense"
                            " for proteins." << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        } break;
        case 0: break;
        break;
        default:
            std::cerr << "Alphabet Reduction must be one of 8, 10, 12 or off" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(foo, parser, "querypart");
    options.queryPart = foo;
    if (options.queryPart == 0)
        options.queryPart = omp_get_max_threads();

    return seqan::ArgumentParser::PARSE_OK;
}


seqan::ArgumentParser::ParseResult
parseCommandLine(LambdaIndexerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("lambda_indexer");
    // Set short description, version, and date.
    setShortDescription(parser, "Indexer for Lambda");
    setVersion(parser, "0.2");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-i DATABASE.fasta\\fP [\\-o INDEXFILE.sa\\fP]");
    addDescription(parser, "This is the application skelleton and you should modify this string.");


    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("i",
                                            "input",
                                            "Database sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "i");
    setValidValues(parser, "input", "fasta fa fna faa");

    addOption(parser, seqan::ArgParseOption("s",
                                            "segfile",
                                            "SEG intervals for database"
                                            "(optional).",
                                            seqan::ArgParseArgument::INPUTFILE));

    setValidValues(parser, "segfile", "seg");


    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "Index of database sequences",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "sa fm");

    addSection(parser, "Program Options");
    addOption(parser, seqan::ArgParseOption("p",
                                            "program",
                                            "Blast Operation Mode.",
                                            seqan::ArgParseArgument::STRING,
                                            "OUT"));
    setValidValues(parser, "program", "blastn blastp blastx tblastn tblastx");
    setDefaultValue(parser, "program", "blastx");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "Alphabet Reduction for AminoAcid "
                                            "Alphabet (0 -> off).",
                                            seqan::ArgParseArgument::INTEGER));
    setDefaultValue(parser, "alph", "2");

    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Remarks");
    addText(parser, "Note that the indeces created are binary and not "
                    "necessarily platform independant.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.

    seqan::getOptionValue(options.dbFile, parser, "input");

    CharString blastMode;
    getOptionValue(blastMode, parser, "program");
    if (blastMode == "blastn")
    {
        options.blastProg   = BlastFormatOptions::BlastN;
        options.alphReduction = 0;

    }
    else if (blastMode == "blastp")
    {
        options.blastProg = BlastFormatOptions::BlastP;

    }
    else if (blastMode == "blastx")
    {
        options.blastProg = BlastFormatOptions::BlastX;

    }
    else if (blastMode == "tblastn")
    {
        options.blastProg = BlastFormatOptions::TBlastN;

    }
    else if (blastMode == "tblastx")
    {
        options.blastProg = BlastFormatOptions::TBlastX;

    }
    else
        return seqan::ArgumentParser::PARSE_ERROR;

    //verifyFileFormat(); TODO

    getOptionValue(options.alphReduction, parser, "alph");
    switch (options.alphReduction)
    {
        case 0: case 1: case 2: case 8: case 10: case 12:
        if (options.blastProg == BlastFormatOptions::BlastN)
        {
            std::cerr << "Alphabet Reduction only makes sense"
                            " for proteins." << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }
        break;
        default:
            std::cerr << "Alphabet Reduction must be one of 8, 10, 12 or off" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
    }


    getOptionValue(options.segFile, parser, "segfile");

    return seqan::ArgumentParser::PARSE_OK;
}



inline void
printOptions(LambdaOptions const & options)
{
    std::cout << "Seed-Length:   " << uint(options.seedLength) << "\n"
              << "Seed-Delta:    " << uint(options.maxSeedDist) << "\n"
              << "Hammingonly:   " << uint(options.hammingOnly) << "\n"
              << "SeedGravity:   " << uint(options.seedGravity) << "\n"
              << "MinSeedLength: " << uint(options.minSeedLength) << "\n"
              << "MinSeedScore:  " << uint(options.minSeedScore) << "\n"
              << "MinSeedEval:   " << uint(options.minSeedEVal) << "\n"
              << "MinSeedBitS:   " << options.minSeedBitS << "\n"
              << "xDropOff:      " << options.xDropOff << "\n"
              << "band:          " << options.band << "\n"
              << "eCutOff:       " << options.eCutOff << "\n"
              << "alphReduction: " << uint(options.alphReduction) << "\n";
}


#endif // header guard
