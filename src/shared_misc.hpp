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
// shared_misc.h: misc stuff for indexer and search
// ==========================================================================

#pragma once

#include <chrono>
#include <dirent.h>
#include <forward_list>
#include <locale>
#include <sys/sysctl.h>
#include <type_traits>
#include <unistd.h>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/to_upper.hpp>

// ============================================================================
// Functions for translation and retranslation
// ============================================================================

template <typename TPos>
inline bool
inRange(TPos const i, TPos const beg, TPos const end)
{
    return ((i >= beg) && (i < end));
}

inline int64_t
intervalOverlap(uint64_t const s1, uint64_t const e1,
                uint64_t const s2, uint64_t const e2)
{
    return std::min(e1, e2) - std::max(s1, s2);
}

inline void
printProgressBar(uint64_t & lastPercent, uint64_t curPerc)
{
    //round down to even
    curPerc = curPerc & ~1;
//     #pragma omp critical(stdout)
    if ((curPerc > lastPercent) && (curPerc <= 100))
    {
        for (uint64_t i = lastPercent + 2; i <= curPerc; i+=2)
        {
            if (i == 100)
                std::cout << "|" << std::flush;
            else if (i % 10 == 0)
                std::cout << ":" << std::flush;
            else
                std::cout << "." << std::flush;
        }
        lastPercent = curPerc;
    }
}

struct alphabet_detection_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;

    template <typename alph>
    using sequence_container = std::vector<alph>;
};

template <typename TAlph, typename TRange>
constexpr bool all_valid(TRange && r)
{
    for (auto const c : r)
        if (!seqan3::char_is_valid_for<TAlph>(c))
            return false;
    return true;
}

inline AlphabetEnum detectSeqFileAlphabet(std::string const & path)
{
    seqan3::sequence_file_input<alphabet_detection_traits, seqan3::fields<seqan3::field::seq>> f{path};

    auto & seq = std::get<0>(*f.begin());

    if (all_valid<seqan3::dna5>(seq))
    {
        return AlphabetEnum::DNA5;
    }
    else if (all_valid<seqan3::dna15>(seq))
    {
        std::cerr << "\nWARNING: You query file was detected as non-standard DNA, but it could be AminoAcid, too.\n"
                    "To explicitly read as AminoAcid, add '--query-alphabet aminoacid'.\n"
                    "To ignore and disable this warning, add '--query-alphabet dna5'.\n";
        return AlphabetEnum::DNA5;
    }
    else if (all_valid<seqan3::aa27>(seq))
    {
        return AlphabetEnum::AMINO_ACID;
    }

    throw std::runtime_error("Your query file contains illegal characters in the first sequence.");

    // unreachable
    return AlphabetEnum::AMINO_ACID;
}

// ----------------------------------------------------------------------------
// print if certain verbosity is set
// ----------------------------------------------------------------------------

template <typename T>
inline void
myPrintImpl(SharedOptions const & /**/,
            T const & first)
{
    std::cout << first;
}

template <typename T, typename ... Args>
inline void
myPrintImpl(SharedOptions const & options,
            T const & first,
            Args const & ... args)
{
    myPrintImpl(options, first);
    myPrintImpl(options, args...);
}

template <typename... Args>
inline void
myPrint(SharedOptions const & options, const int verbose, Args const &... args)
{
    if (options.verbosity >= verbose)
    {
        myPrintImpl(options, args...);
        std::cout << std::flush;
    }
}

// ----------------------------------------------------------------------------
// Function sysTime()
// ----------------------------------------------------------------------------

inline double sysTime()
{
//     return std::time(NULL);
    return static_cast<double>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count()) / 1000;
//     return std::chrono::system_clock::now().time_since_epoch().count() / 1000;
}

// ----------------------------------------------------------------------------
// Function fileSize()
// ----------------------------------------------------------------------------

inline uint64_t fileSize(char const * fileName)
{
    struct stat st;
    if (stat(fileName, &st) != 0)
        throw std::runtime_error{"Could not read File.\n"};
    return st.st_size;
}

// ----------------------------------------------------------------------------
// Function dirSize()
// ----------------------------------------------------------------------------

inline uint64_t dirSize(char const * dirName)
{
    DIR *d;
    struct dirent *de;
    struct stat buf;
    int exists;
    uint64_t total_size;

    d = opendir(dirName);
    if (d == NULL)
        throw std::runtime_error{"Could not read index directory.\n"};

    total_size = 0;

    for (de = readdir(d); de != NULL; de = readdir(d))
    {
        std::string curPath = dirName + std::string{"/"} + de->d_name;
        exists = stat(curPath.c_str(), &buf);
        if (exists < 0)
        {
            closedir(d);
            throw std::runtime_error{"Could not read index directory.\n"};
        } else
        {
            total_size += buf.st_size;
        }
    }
    closedir(d);
    return total_size;
}

// ----------------------------------------------------------------------------
// Function fileSize()
// ----------------------------------------------------------------------------

inline uint64_t getTotalSystemMemory()
{
#if defined(__APPLE__)
    uint64_t mem;
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
    return mem;
#elif defined(__unix__)
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#else
#   error "no way to get phys pages"
#endif
}
