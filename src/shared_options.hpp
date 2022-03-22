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
// shared_options.h: contains the options and argument parser
// ==========================================================================

#pragma once

#include <bitset>
#include <cstdio>
#include <thread>
#include <unistd.h>

#include <seqan3/argument_parser/all.hpp>
#include <filesystem>

#include "shared_definitions.hpp"

// --------------------------------------------------------------------------
// Class SharedOptions
// --------------------------------------------------------------------------

struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose
    int32_t     verbosity = 1;

    std::string commandLine;

    std::filesystem::path indexFilePath;

    index_file_options indexFileOptions{};

    bool        nucleotide_mode = false;
    bool        need_to_translate = false;

    bool        isTerm          = true;
    unsigned    terminalCols    = 80;

    uint64_t    threads         = std::thread::hardware_concurrency();

    bool        hasSTaxIds      = false;

    SharedOptions()
    {
//         isTerm = seqan::isTerminal();
//         if (isTerm)
//         {
//             unsigned _rows;
//             seqan::getTerminalSize(terminalCols, _rows);
//         }
    }
};


inline void sharedSetup(seqan3::argument_parser & parser)
{
    // Set short description, version, and date
    parser.info.version = SEQAN_APP_VERSION;
    parser.info.date = __DATE__;
    parser.info.citation = "Hauswedell et al (2014); doi: 10.1093/bioinformatics/btu439";
    parser.info.short_copyright = "2013-2020 Hannes Hauswedell, released under the GNU AGPL v3 (or later); "
                                  "2016-2020 Knut Reinert and Freie Universität Berlin, released under the 3-clause-BSDL";
    parser.info.long_copyright = " Copyright (c) 2013-2020, Hannes Hauswedell\n"
                                 " All rights reserved.\n"
                                 "\n"
                                 " This program is free software: you can redistribute it and/or modify\n"
                                 " it under the terms of the GNU Affero General Public License as\n"
                                 " published by the Free Software Foundation, either version 3 of the\n"
                                 " License, or (at your option) any later version.\n"
                                 "\n"
                                 " Lambda is distributed in the hope that it will be useful,\n"
                                 " but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
                                 " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
                                 " GNU General Public License for more details.\n"
                                 "\n"
                                 " You should have received a copy of the GNU Affero General Public License\n"
                                 " along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
                                 "\n"
                                 " Copyright (c) 2016-2020 Knut Reinert and Freie Universität Berlin\n"
                                 " All rights reserved.\n"
                                 "\n"
                                 " Redistribution and use in source and binary forms, with or without\n"
                                 " modification, are permitted provided that the following conditions are met:\n"
                                 "\n"
                                 " * Redistributions of source code must retain the above copyright\n"
                                 "   notice, this list of conditions and the following disclaimer.\n"
                                 " * Redistributions in binary form must reproduce the above copyright\n"
                                 "   notice, this list of conditions and the following disclaimer in the\n"
                                 "   documentation and/or other materials provided with the distribution.\n"
                                 " * Neither the name of Knut Reinert or the FU Berlin nor the names of\n"
                                 "   its contributors may be used to endorse or promote products derived\n"
                                 "   from this software without specific prior written permission.\n"
                                 "\n"
                                 " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n"
                                 " AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n"
                                 " IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n"
                                 " ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE\n"
                                 " FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\n"
                                 " DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
                                 " SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
                                 " CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
                                 " LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY\n"
                                 " OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n"
                                 " DAMAGE.\n";
    parser.info.description.push_back("Lambda is a local aligner optimized for many query "
                                      "sequences and searches in protein space. It is compatible to BLAST, but "
                                      "much faster than BLAST and many other comparable tools.");


    parser.info.description.push_back("Detailed information is available in the wiki: "
                                      "<https://github.com/seqan/lambda/wiki>");
}
