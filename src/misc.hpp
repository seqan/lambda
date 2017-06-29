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
// store.h: contains types and definitions for storing sequences and indices
// ==========================================================================

#ifndef SEQAN_LAMBDA_MISC_H_
#define SEQAN_LAMBDA_MISC_H_

#include <unistd.h>
#include <type_traits>
#include <forward_list>
#include <sys/sysctl.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/blast.h>

#include "options.hpp"

using namespace seqan;

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

// template <typename TString, typename TSpec>
// inline void
// assign(StringSet<ModifiedString<TString, TSpec>, Owner<ConcatDirect<>> > & target,
//        StringSet<TString, Owner<ConcatDirect<>> > & source)
// {
//     target.limits = source.limits;
//     target.concat._host = &source.concat;
// }

/*
template <typename T>
inline uint64_t
length(std::deque<T> const & list)
{
    return list.size();
}
template <typename T>
inline uint64_t
length(std::forward_list<T> const & list)
{
    return std::distance(list.begin(), list.end());
}*/

template <typename TAlph>
inline std::basic_ostream<char> &
operator<<(std::basic_ostream<char> & out,
           const Iter<const String<SimpleType<unsigned char,TAlph>,
                                    seqan::Packed<> >,
                      seqan::Packed<> > it)
{
    out << *it;
    return out;
}


template <typename T1, typename T2>
inline uint64_t
quickHamming(T1 const & s1, T2 const & s2)
{
    SEQAN_ASSERT_EQ(length(s1), length(s2));

    uint64_t ret = 0;

    for (uint64_t i = 0; i < length(s1); ++i)
        if (s1[i] != s2[i])
            ++ret;

    return ret;
}

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
                std::cout << "|\n" << std::flush;
            else if (i % 10 == 0)
                std::cout << ":" << std::flush;
            else
                std::cout << "." << std::flush;
        }
        lastPercent = curPerc;
    }
}


template <typename TSource0, typename TGapsSpec0,
          typename TSource1, typename TGapsSpec1,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignContext>
inline TScoreValue
localAlignment2(Gaps<TSource0, TGapsSpec0> & row0,
                Gaps<TSource1, TGapsSpec1> & row1,
                Score<TScoreValue, TScoreSpec> const & scoringScheme,
                int const lowerDiag,
                int const upperDiag,
                TAlignContext & alignContext)
{
    clear(alignContext.traceSegment);

    typedef FreeEndGaps_<True, True, True, True> TFreeEndGaps;
    typedef AlignConfig2<LocalAlignment_<>,
                         DPBand,
                         TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace,
                                                      GapsLeft> > > TAlignConfig;

    TScoreValue score;
    DPScoutState_<Default> scoutState;
    score = _setUpAndRunAlignment(alignContext.dpContext,
                                  alignContext.traceSegment,
                                  scoutState,
                                  row0,
                                  row1,
                                  scoringScheme,
                                  TAlignConfig(lowerDiag, upperDiag));

    _adaptTraceSegmentsTo(row0, row1, alignContext.traceSegment);
    return score;
}

// ----------------------------------------------------------------------------
// Tag Truncate_ (private tag to signify truncating of IDs while reading)
// ----------------------------------------------------------------------------

struct Truncate_;

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); an overload that truncates Ids at first Whitespace
// ----------------------------------------------------------------------------

template <typename TSeqStringSet, typename TSpec, typename TRunnable>
inline void
_myReadRecordsImpl(TCDStringSet<String<char>> & meta,
                   TSeqStringSet & seq,
                   FormattedFile<Fastq, Input, TSpec> & file,
                   TRunnable && runnable)
{
    typedef typename SeqFileBuffer_<TSeqStringSet, TSpec>::Type TSeqBuffer;

    TSeqBuffer seqBuffer;

    // reuse the memory of context(file).buffer for seqBuffer (which has a different type but same sizeof(Alphabet))
    swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
    seqBuffer.data_capacity = context(file).buffer[1].data_capacity;

    for (uint64_t count = 0; !atEnd(file); ++count) // count not used for abort condition
    {
        readRecord(context(file).buffer[0], seqBuffer, file);

        // run whatever magic we are pushing in:
        runnable(context(file).buffer[0], count);

        appendValue(meta, context(file).buffer[0]);
        appendValue(seq, seqBuffer);
    }

    swapPtr(seqBuffer.data_begin, context(file).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(file).buffer[1].data_end);
    context(file).buffer[1].data_capacity = seqBuffer.data_capacity;
    seqBuffer.data_capacity = 0;
}

// ----------------------------------------------------------------------------
// Generic Sequence loading
// ----------------------------------------------------------------------------

template <typename TSpec1,
          typename TSpec2,
          typename TFile,
          typename TRunnable>
inline int
myReadRecords(TCDStringSet<String<char, TSpec1>> & ids,
              TCDStringSet<String<Dna5, TSpec2>> & seqs,
              TFile                              & file,
              TRunnable                         && runnable)
{
    TCDStringSet<String<Iupac>> tmpSeqs; // all IUPAC nucleic acid characters are valid input
    try
    {
        _myReadRecordsImpl(ids, tmpSeqs, file, std::forward<TRunnable>(runnable));
    }
    catch(ParseError const & e)
    {
        std::cerr << "\nParseError thrown: " << e.what() << '\n'
                  << "Make sure that the file is standards compliant. If you get an unexpected character warning "
                  << "make sure you have set the right program parameter (-p), i.e. "
                  << "Lambda expected nucleic acid alphabet, maybe the file was protein?\n";
        return -1;
    }

    seqs = tmpSeqs; // convert IUPAC alphabet to Dna5

    return 0;
}

template <typename TSpec1,
          typename TSpec2,
          typename TFile,
          typename TRunnable>
inline int
myReadRecords(TCDStringSet<String<char, TSpec1>>       & ids,
              TCDStringSet<String<AminoAcid, TSpec2>>  & seqs,
              TFile                                    & file,
              TRunnable                               && runnable)
{
    try
    {
        _myReadRecordsImpl(ids, seqs, file, std::forward<TRunnable>(runnable));
    }
    catch(ParseError const & e)
    {
        std::cerr << "\nParseError thrown: " << e.what() << '\n'
                  << "Make sure that the file is standards compliant.\n";
        return -1;
    }

    if (length(seqs) > 0)
    {
        // warn if sequences look like DNA
        if (CharString(String<Dna5>(CharString(seqs[0]))) == CharString(seqs[0]))
            std::cout << "\nWarning: The first query sequence looks like nucleic acid, but amino acid is expected.\n"
                         "           Make sure you have set the right program parameter (-p).\n";
    }

    return 0;
}

template <typename TSpec1,
          typename TCharSpec,
          typename TSpec2,
          typename TFile>
inline int
myReadRecords(TCDStringSet<String<char, TSpec1>>       & ids,
              TCDStringSet<String<TCharSpec, TSpec2>>  & seqs,
              TFile                                    & file)
{
    return  myReadRecords(ids, seqs, file, [] (auto const &, uint64_t const) {});
}

template <typename TLocalHolder>
inline int
_bandSize(uint64_t const seqLength, TLocalHolder & lH)
{
    switch (lH.options.band)
    {
        case -3: case -2:
        {
            int ret = 0;
            auto fit = lH.bandTable.find(seqLength);
            if (fit != lH.bandTable.end())
            {
                ret = fit->second;
            } else
            {
                if (lH.options.band == -3)
                    ret = ceil(std::log2(seqLength));
                else
                    ret = floor(sqrt(seqLength));
            }
            lH.bandTable[seqLength] = ret;
            return ret;
        } break;
        case -1:
            return std::numeric_limits<int>::max();
        default:
            return lH.options.band;
    }
}

// ----------------------------------------------------------------------------
// truncate sequences
// ----------------------------------------------------------------------------

// template <typename TString, typename TSpec>
// inline void
// _debug_shorten(StringSet<TString, TSpec > & seqs, unsigned const len)
// {
//     StringSet<TString, TSpec > copySeqs;
//     reserve(copySeqs.concat, length(seqs)*len, Exact());
// 
//     for (TString const & s : seqs)
//         if (length(s) >= len)
//             appendValue(copySeqs, prefix(s, len), Exact());
// 
//     clear(seqs);
//     reserve(seqs.concat, length(copySeqs)*len, Exact());
//     for (TString const & s : copySeqs)
//         appendValue(seqs, s);
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

template <typename ... Args>
inline void
myPrintImplThread(SharedOptions const & options,
//                   T const & first,
                  Args const & ... args)
{
    SEQAN_OMP_PRAGMA(critical(stdout))
    {
//                 std::cout << "\033[" << omp_get_thread_num() << "B";
//                 std::cout << "\033E";
        if (options.isTerm)
        {
            for (unsigned char i=0; i< omp_get_thread_num(); ++i)
                std::cout << std::endl;
            std::cout << "\033[K";
        }
        std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << "| ";

        myPrintImpl(options, args...);
        std::cout << "\n" << std::flush;
        if (options.isTerm)
            std::cout << "\033[" << omp_get_thread_num()+1 << "A";
    }
}

template <typename... Args>
inline void
myPrint(SharedOptions const & options, const int verbose, Args const &... args)
{
    if (options.verbosity >= verbose)
    {
        #if defined(_OPENMP)
        if (omp_in_parallel())
            myPrintImplThread(options, args...);
        else
        #endif
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
               LambdaOptions const & options,
               const int verbose,
               Args const & ... args)
{
    if (options.verbosity >= verbose)
        appendToStatusImpl(status, args...);
}

// ----------------------------------------------------------------------------
// Function computeEValueThreadSafe
// ----------------------------------------------------------------------------

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeEValueThreadSafe(TBlastMatch & match,
                        uint64_t ql,
                        BlastIOContext<TScore, p, h> & context)
{
#if defined(__FreeBSD__)
    // && version < 11 && defined(STDLIB_LLVM) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=192320
    // || version >= 11 && defined(STDLIB_GNU) because of https://bugs.freebsd.org/bugzilla/show_bug.cgi?id=215709
    static std::vector<std::unordered_map<uint64_t, uint64_t>> _cachedLengthAdjustmentsArray(omp_get_num_threads());
    std::unordered_map<uint64_t, uint64_t> & _cachedLengthAdjustments = _cachedLengthAdjustmentsArray[omp_get_thread_num()];
#else
    static thread_local std::unordered_map<uint64_t, uint64_t> _cachedLengthAdjustments;
#endif

    // convert to 64bit and divide for translated sequences
    ql = ql / (qIsTranslated(context.blastProgram) ? 3 : 1);
    // length adjustment not yet computed
    if (_cachedLengthAdjustments.find(ql) == _cachedLengthAdjustments.end())
        _cachedLengthAdjustments[ql] = _lengthAdjustment(context.dbTotalLength, ql, context.scoringScheme);

    uint64_t adj = _cachedLengthAdjustments[ql];

    match.eValue = _computeEValue(match.alignStats.alignmentScore,
                                  ql - adj,
                                  context.dbTotalLength - adj,
                                  context.scoringScheme);
    return match.eValue;
}

// ----------------------------------------------------------------------------
// compute LCA
// ----------------------------------------------------------------------------

template <typename T, typename T2>
T computeLCA(String<T> const & taxParents, String<T2> const & taxHeights, T n1, T n2)
{
    if (n1 == n2)
        return n1;

    // move up so that nodes are on same height
    for (auto i = taxHeights[n1]; i > taxHeights[n2]; --i)
        n1 = taxParents[n1];

    for (auto i = taxHeights[n2]; i > taxHeights[n1]; --i)
        n2 = taxParents[n2];

    while ((n1 != 0) && ( n2 != 0))
    {
        // common ancestor
        if (n1 == n2)
            return n1;

        // move up in parallel
        n1 = taxParents[n1];
        n2 = taxParents[n2];
    }

    SEQAN_FAIL("One of the paths didn't lead to root.");
    return 0; // avoid warnings on clang
}

// ----------------------------------------------------------------------------
// Function fileSize()
// ----------------------------------------------------------------------------

auto fileSize(char const * fileName)
{
    struct stat st;
    if (stat(fileName, &st) != 0)
        throw std::runtime_error{"ERROR: Could not read File!"};
    return st.st_size;
}

// ----------------------------------------------------------------------------
// Function dirSize()
// ----------------------------------------------------------------------------

unsigned long long dirSize(char const * dirName)
{
    DIR *d;
    struct dirent *de;
    struct stat buf;
    int exists;
    unsigned long long total_size;

    d = opendir(dirName);
    if (d == NULL)
        throw std::runtime_error{"ERROR: Could not read index directory!"};

    total_size = 0;

    for (de = readdir(d); de != NULL; de = readdir(d))
    {
        std::string curPath = dirName + std::string{"/"} + de->d_name;
        exists = stat(curPath.c_str(), &buf);
        if (exists < 0)
        {
            closedir(d);
            throw std::runtime_error{"ERROR: Could not read index directory!"};
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

unsigned long long getTotalSystemMemory()
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


// ----------------------------------------------------------------------------
// remove tag type
// ----------------------------------------------------------------------------

// template <typename T>
// T unTag(Tag<T> const & /**/)
// {
//     return T();
// }

// ----------------------------------------------------------------------------
// get plus-minus-range with bounds-checking for unsigned types
// ----------------------------------------------------------------------------

// template <typename TNum, typename TNum2>
// inline TNum
// _protectUnderflow(const TNum n, const TNum2 s)
// {
//     const TNum r = n -s;
//     return std::min(r, n);
// }
// 
// template <typename TNum, typename TNum2>
// inline TNum
// _protectOverflow(const TNum n, const TNum2 s)
// {
//     const TNum r = n + s;
//     return std::max(r, n);
// }
// 
// template <typename TGaps>
// inline bool
// _startsWithGap(TGaps const & gaps)
// {
//     SEQAN_ASSERT(length(gaps._array) > 0);
//     return (gaps._array[0] != 0);
// }
// 
// template <typename TGaps>
// inline int
// _endsWithGap(TGaps const & gaps)
// {
//     SEQAN_ASSERT(length(gaps._array) > 0);
//     if ((length(gaps._array)-1) % 2 == 1)
//         return -1;
//     return ((gaps._array[length(gaps._array)-1] != 0) ? 1 : 0);
// }
// 
// template <typename TGaps, typename TSeq>
// inline void
// _prependNonGaps(TGaps & gaps, TSeq const & seq)
// {
//     if (_startsWithGap(gaps))
//     {
//         insertValue(gaps._array, 0, length(seq)); // new non-gap column
//         insertValue(gaps._array, 0, 0); // empty gaps column
//     }
//     else
//     {
//         gaps._array[1] += length(seq);
//     }
// 
//     insert(value(gaps._source), 0, seq);
//     setBeginPosition(gaps, 0);
// }
// 
// template <typename TGaps, typename TSeq>
// inline void
// _appendNonGaps(TGaps & gaps, TSeq const & seq)
// {
//     switch (_endsWithGap(gaps))
//     {
//         case -1:
//         case  1:
//             appendValue(gaps._array, length(seq)); // new non-gap column
//             break;
//         case 0:
//             gaps._array[1] += length(seq);
//     }
//     append(value(gaps._source), seq);
//     setEndPosition(gaps, length(value(gaps._source)));
// }

#endif // header guard
