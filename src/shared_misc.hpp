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
// store.h: contains types and definitions for storing sequences and indices
// ==========================================================================

#ifndef LAMBDA_SHARED_MISC_H_
#define LAMBDA_SHARED_MISC_H_

#include <chrono>
#include <dirent.h>
#include <forward_list>
#include <locale>
#include <sys/sysctl.h>
#include <type_traits>
#include <unistd.h>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/view/to_upper.hpp>

// #include <seqan/basic.h>
// #include <seqan/sequence.h>
// #include <seqan/seq_io.h>
// #include <seqan/index.h>
// #include <seqan/align.h>
// #include <seqan/blast.h>

// using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// makes partial function specialization convenient
template <bool condition>
using MyEnableIf = typename std::enable_if<condition, int>::type;

// ============================================================================
// Functions for translation and retranslation
// ============================================================================

// template <typename TAlph>
// inline std::basic_ostream<char> &
// operator<<(std::basic_ostream<char> & out,
//            const seqan::Iter<const seqan::String<seqan::SimpleType<unsigned char,TAlph>,
//                                     seqan::Packed<> >,
//                       seqan::Packed<> > it)
// {
//     out << *it;
//     return out;
// }

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
    using sequence_container = std::basic_string<alph>;
};

AlphabetEnum detectSeqFileAlphabet(std::string const & path)
{
    seqan3::sequence_file_input<alphabet_detection_traits, seqan3::fields<seqan3::field::SEQ>> f{path};

    std::string & seq = std::get<0>(*f.begin());

    if (std::ranges::equal(seq | seqan3::view::char_to<seqan3::dna5> | seqan3::view::to_char,
                           seq | seqan3::view::to_upper))
    {
        return AlphabetEnum::DNA5;
    }
    else if (std::ranges::equal(seq | seqan3::view::char_to<seqan3::dna15> | seqan3::view::to_char,
                                seq | seqan3::view::to_upper))
    {
        std::cerr << "\nWARNING: You query file was detected as non-standard DNA, but it could be AminoAcid, too.\n"
                    "To explicitly read as AminoAcid, add '--query-alphabet aminoacid'.\n"
                    "To ignore and disable this warning, add '--query-alphabet dna5'.\n";
        return AlphabetEnum::DNA5;
    }
    else if (std::ranges::equal(seq | seqan3::view::char_to<seqan3::aa27> | seqan3::view::to_char,
                                seq | seqan3::view::to_upper))
    {
        return AlphabetEnum::AMINO_ACID;
    }

    throw std::runtime_error("Your query file contains illegal characters in the first sequence.");

    // unreachable
    return AlphabetEnum::AMINO_ACID;
}

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); an overload that truncates Ids at first Whitespace
// ----------------------------------------------------------------------------

// template <typename TSeqStringSet, typename TSpec, typename TRunnable>
// inline void
// _myReadRecordsImpl(TCDStringSet<seqan::String<char>> & meta,
//                    TSeqStringSet & seq,
//                    FormattedFile<Fastq, Input, TSpec> & file,
//                    TRunnable && runnable)
// {
//     typedef typename SeqFileBuffer_<TSeqStringSet, TSpec>::Type TSeqBuffer;
//
//     TSeqBuffer seqBuffer;
//
//     // reuse the memory of context(file).buffer for seqBuffer (which has a different type but same sizeof(Alphabet))
//     swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
//     swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
//     seqBuffer.data_capacity = context(file).buffer[1].data_capacity;
//
//     for (uint64_t count = 0; !atEnd(file); ++count) // count not used for abort condition
//     {
//         readRecord(context(file).buffer[0], seqBuffer, file);
//
//         // run whatever magic we are pushing in:
//         runnable(context(file).buffer[0], count);
//
//         appendValue(meta, context(file).buffer[0]);
//         appendValue(seq, seqBuffer);
//     }
//
//     swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
//     swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
//     context(file).buffer[1].data_capacity = seqBuffer.data_capacity;
//     seqBuffer.data_capacity = 0;
// }

// ----------------------------------------------------------------------------
// Generic Sequence loading
// ----------------------------------------------------------------------------

// template <typename TSpec1,
//           typename TSpec2,
//           typename TFile,
//           typename TRunnable>
// void
// myReadRecords(TCDStringSet<seqan::String<char, TSpec1>> & ids,
//               TCDStringSet<seqan::String<seqan::Dna5, TSpec2>> & seqs,
//               TFile                              & file,
//               TRunnable                         && runnable)
// {
//     TCDStringSet<seqan::String<seqan::Iupac>> tmpSeqs; // all IUPAC nucleic acid characters are valid input
//     try
//     {
//         _myReadRecordsImpl(ids, tmpSeqs, file, std::forward<TRunnable>(runnable));
//     }
//     catch(ParseError const & e)
//     {
//         std::string err;
//         err += "\nParseError thrown: ";
//         err += e.what();
//         err += "\nMake sure that the file is standards compliant. If you get an unexpected character warning "
//                "make sure you have set the right program parameter (-p), i.e. "
//                "Lambda expected nucleic acid alphabet, maybe the file was protein?\n";
//         throw std::runtime_error(err);
//     }
//
//     seqs = tmpSeqs; // convert IUPAC alphabet to Dna5
// }
//
// template <typename TSpec1,
//           typename TSpec2,
//           typename TFile,
//           typename TRunnable>
// void
// myReadRecords(TCDStringSet<seqan::String<char, TSpec1>>       & ids,
//               TCDStringSet<seqan::String<seqan::AminoAcid, TSpec2>>  & seqs,
//               TFile                                    & file,
//               TRunnable                               && runnable)
// {
//     try
//     {
//         _myReadRecordsImpl(ids, seqs, file, std::forward<TRunnable>(runnable));
//     }
//     catch(seqan::ParseError const & e)
//     {
//         std::string err;
//         err += "\nParseError thrown: ";
//         err += e.what();
//         err += "\nMake sure that the file is standards compliant.\n";
//         throw std::runtime_error(err);
//     }
//
//     if (seqan::length(seqs) > 0)
//     {
//         // warn if sequences look like DNA
//         if (seqan::CharString(seqan::String<seqan::Dna5>(seqan::CharString(seqs[0]))) == seqan::CharString(seqs[0]))
//             std::cout << "\nWarning: The first query sequence looks like nucleic acid, but amino acid is expected.\n"
//                          "           Make sure you have set the right program parameter (-p).\n";
//     }
// }
//
// template <typename TSpec1,
//           typename TCharSpec,
//           typename TSpec2,
//           typename TFile>
// void
// myReadRecords(TCDStringSet<seqan::String<char, TSpec1>>       & ids,
//               TCDStringSet<seqan::String<TCharSpec, TSpec2>>  & seqs,
//               TFile                                    & file)
// {
//     myReadRecords(ids, seqs, file, [] (auto const &, uint64_t const) {});
// }

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

inline void
myPrintImpl(SharedOptions const & options,
            std::stringstream const & first)
{
    std::string str = first.str();
//     std::cerr << "terminal cols: " << options.terminalCols
//               << " str.size() " << str.size() << "\n";
    if (options.isTerm && (str.size() >= (options.terminalCols -12)))
        std::cout << str.substr(str.size()-options.terminalCols+12,
                                options.terminalCols);
    else
        std::cout << str;
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

// template <typename ... Args>
// inline void
// myPrintImplThread(SharedOptions const & options,
// //                   T const & first,
//                   Args const & ... args)
// {
//     SEQAN_OMP_PRAGMA(critical(stdout))
//     {
// //                 std::cout << "\033[" << omp_get_thread_num() << "B";
// //                 std::cout << "\033E";
//         if (options.isTerm)
//         {
//             for (unsigned char i=0; i< omp_get_thread_num(); ++i)
//                 std::cout << std::endl;
//             std::cout << "\033[K";
//         }
//         std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << "| ";
//
//         myPrintImpl(options, args...);
//         std::cout << "\n" << std::flush;
//         if (options.isTerm)
//             std::cout << "\033[" << omp_get_thread_num()+1 << "A";
//     }
// }

template <typename... Args>
inline void
myPrint(SharedOptions const & options, const int verbose, Args const &... args)
{
    if (options.verbosity >= verbose)
    {
//         #if defined(_OPENMP)
//         if (omp_in_parallel())
//             myPrintImplThread(options, args...);
//         else
//         #endif
            myPrintImpl(options, args...);

        std::cout << std::flush;
    }
}

template <typename T>
inline void
appendToStatusImpl(std::stringstream & status,
                   T const & first)
{
    status << first;
}

template <typename T, typename ... Args>
inline void
appendToStatusImpl(std::stringstream & status,
                   T const & first,
                   Args const & ... args)
{
    appendToStatusImpl(status, first);
    appendToStatusImpl(status, args...);
}

template <typename... Args>
inline void
appendToStatus(std::stringstream & status,
               SharedOptions const & options,
               const int verbose,
               Args const & ... args)
{
    if (options.verbosity >= verbose)
        appendToStatusImpl(status, args...);
}

// ----------------------------------------------------------------------------
// Function sysTime()
// ----------------------------------------------------------------------------

inline double sysTime()
{
    //TODO
    return 0;
//     return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
}

// ----------------------------------------------------------------------------
// Function fileSize()
// ----------------------------------------------------------------------------

uint64_t fileSize(char const * fileName)
{
    struct stat st;
    if (stat(fileName, &st) != 0)
        throw std::runtime_error{"Could not read File.\n"};
    return st.st_size;
}

// ----------------------------------------------------------------------------
// Function dirSize()
// ----------------------------------------------------------------------------

uint64_t dirSize(char const * dirName)
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

uint64_t getTotalSystemMemory()
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

#endif // header guard
