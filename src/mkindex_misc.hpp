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
// lambda_indexer_misc.hpp: misc stuff for indexer
// ==========================================================================

#ifndef LAMBDA_INDEXER_MISC_HPP_
#define LAMBDA_INDEXER_MISC_HPP_

template <typename TString>
void getCwd(TString & string)
{
    char cwd[1000];

#ifdef PLATFORM_WINDOWS
    _getcwd(cwd, 1000);
#else
    getcwd(cwd, 1000);
#endif

    assign(string, cwd);
}

template <typename TString, typename TValue>
bool setEnv(TString const & key, TValue & value)
{
#ifdef PLATFORM_WINDOWS
    return !_putenv_s(toCString(key), toCString(value));
#else
    return !setenv(toCString(key), toCString(value), true);
#endif
}

// ----------------------------------------------------------------------------
// Class ComparisonCounter
// ----------------------------------------------------------------------------

#if 0 // why is this still here?

template <typename TText, typename TSpec>
struct ComparisonCounter;

// no counting
template <typename TText>
struct ComparisonCounter<TText, Nothing>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
    uint64_t _lastPercent = 0;
    ComparisonCounter(TText const &,
                      uint64_t expectedComparisons = 0u)
    {
        (void)expectedComparisons;
    }

    // may be constexpr in c++14
    inline void inc() const
    {}
};

// every thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::false_type>
#else
template <typename TText, typename TSpec>
struct ComparisonCounter
#endif
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.2 * double(l) * std::log(l) / std::log(2);
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
        _comparisons = 0;
        _lastPercent = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        uint64_t comp = ++_comparisons;
        // it is not important that the henceforth _comparisons be actually
        // the same value (might not be due to SMP)

        // progress reporting
        if (comp & _checkEveryNHits)
        {
            uint64_t curPerc = comp * 100 / _expectedComparisons;
            if (curPerc < 100)
                printProgressBar(_lastPercent, curPerc);
        }
    }
};

// only one thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::true_type>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.2 * double(l) * std::log(l) / std::log(2) /
                                   omp_get_max_threads();
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
//         _comparisons = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        if (omp_get_thread_num() == 0) // only one thread counts
        {
            uint64_t comp = ++_comparisons;

            // progress reporting
            if (comp & _checkEveryNHits)
            {
                uint64_t curPerc = comp * 100 / _expectedComparisons;
                if (curPerc < 100)
                    printProgressBar(_lastPercent, curPerc);
            }
        }
    }
};
#endif

#endif

// ----------------------------------------------------------------------------
// function _readMappingFileNCBI
// ----------------------------------------------------------------------------

template <typename TInputIterator,
          typename TStaxIDs>
void
_readMappingFileUniProt(TInputIterator                                        & fit,
                        TStaxIDs                                              & sTaxIds,
                        std::vector<bool>                                     & taxIdIsPresent,
                        std::unordered_map<std::string, uint64_t>       const & accToIdRank)
{
    // skip line with headers
    skipLine(fit);

    //TODO this is too slow, investigate whether its the lookup or the allocs
    std::string acc;
    std::string nextColumn;
    while (!atEnd(fit))
    {
        clear(acc);
        clear(nextColumn);

        // read accession number
        readUntil(acc, fit, IsBlank());
        // skip whitespace
        skipUntil(fit, IsAlphaNum());
        // read accession number
        readUntil(nextColumn, fit, IsBlank());

        if ((nextColumn == "NCBI_TaxID") && (accToIdRank.count(acc) == 1))
        {
            auto & sTaxIdV = sTaxIds[accToIdRank.at(acc)];
            // skip whitespace
            skipUntil(fit, IsAlphaNum());
            // read tax id
            clear(nextColumn);
            readUntil(nextColumn, fit, IsWhitespace());
            uint32_t idNum = 0;
            try
            {
                idNum = lexicalCast<uint32_t>(nextColumn);
            }
            catch (BadLexicalCast const & badCast)
            {
                throw std::runtime_error(
                    std::string("Error: Expected taxonomical ID, but got something I couldn't read: ") +
                    std::string(badCast.what()) + "\n");
            }
            appendValue(sTaxIdV, idNum);
            if (taxIdIsPresent.size() < idNum + 1)
                taxIdIsPresent.resize(idNum + 1);
            taxIdIsPresent[idNum] = true;
        }

        skipLine(fit);
    }
}

template <typename TInputIterator,
          typename TStaxIDs>
void
_readMappingFileNCBI(TInputIterator                                        & fit,
                     TStaxIDs                                              & sTaxIds,
                     std::vector<bool>                                     & taxIdIsPresent,
                     std::unordered_map<std::string, uint64_t>       const & accToIdRank)
{
    // skip line with headers
    skipLine(fit);

    //TODO this is too slow, investigate whether its the lookup or the allocs
    std::string buf;
    while (!atEnd(fit))
    {
        clear(buf);
        // read accession number
        readUntil(buf, fit, IsBlank());
        // we have a sequence with this ID in our database
        if (accToIdRank.count(buf) == 1)
        {
            auto & sTaxIdV = sTaxIds[accToIdRank.at(buf)];
            // skip whitespace
            skipUntil(fit, IsAlphaNum());
            // skip versioned acc
            skipUntil(fit, IsBlank());
            // skip whitespace
            skipUntil(fit, IsAlphaNum());
            // read tax id
            clear(buf);
            readUntil(buf, fit, IsBlank());
            uint32_t idNum = 0;
            try
            {
                idNum = lexicalCast<uint32_t>(buf);
            }
            catch (BadLexicalCast const & badCast)
            {
                throw std::runtime_error(
                    std::string("Error: Expected taxonomical ID, but got something I couldn't read: ") +
                    std::string(badCast.what()) + "\n");
            }
            appendValue(sTaxIdV, idNum);
            if (taxIdIsPresent.size() < idNum + 1)
                taxIdIsPresent.resize(idNum + 1);
            taxIdIsPresent[idNum] = true;
        }

        skipLine(fit);
    }
}

/// REGEX version is 5x slower, but verifies file format correctness
// template <typename TInputIterator,
//           typename TStaxIDs>
// inline int
// _readMappingFileNCBI2(TInputIterator                                        & fit,
//                      TStaxIDs                                              & sTaxIds,
//                      std::vector<bool>                                     & taxIdIsPresent,
//                      std::unordered_map<std::string, uint64_t>       const & accToIdRank)
// {
//     // skip line with headers
//     skipLine(fit);
//
//     std::string buf;
//     //                        ACC       ACC.VERSION    taxid     gi
//     std::regex const lineRE{"(\\w+)\\s(\\w+\\.?\\d*)\\s(\\d+)\\s(\\d+)"};
//     std::smatch baseMatch;
//     uint32_t idNum = 0;
//     StringSet<std::string>
//
//     while (!atEnd(fit))
//     {
//         clear(buf);
//         // read line
//         readLine(buf, fit);
//
//         if (std::regex_match(buf, baseMatch, lineRE) && (baseMatch.size() == 5))
//         {
//             // we have a sequence with this ID in our database
//             if (accToIdRank.count(baseMatch[1]) == 1)
//             {
//                 auto & sTaxIdV = sTaxIds[accToIdRank.at(baseMatch[1])];
//                 idNum = 0;
//                 try
//                 {
//                     idNum = lexicalCast<uint32_t>(std::string(baseMatch[3]));
//                 }
//                 catch (BadLexicalCast const & badCast)
//                 {
//                     std::cerr << "Error: Expected taxonomical ID, but got something I couldn't read: "
//                             << badCast.what() << "\n";
//                     return -1;
//                 }
//
//                 appendValue(sTaxIdV, idNum);
//                 if (taxIdIsPresent.size() < idNum + 1)
//                     taxIdIsPresent.resize(idNum + 1);
//                 taxIdIsPresent[idNum] = true;
//             }
//
//         } else
//         {
//             std::cerr << "ERROR: The following line in the mapping file, did not satisfy the regex!\n"
//                       << "       " << buf << "\n\n";
//             return -1;
//         }
//     }
//
//     return 0;
// }

// ============================================================================
// Parallel BWT construction
// ============================================================================

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig,
          typename TOtherText,
          typename TSA,
          typename TCallback>
void
createRankDictionaryProgress(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> & lf,
                             TOtherText const & text,
                             TSA const & sa,
                             TCallback && progress)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;

    // Resize the RankDictionary.
    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);
    resize(lf.sentinels, seqNum + totalLen, Exact());
    resize(lf.bwt, seqNum + totalLen, Exact());

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 0; i < seqNum; ++i)
    {
        if (length(text[seqNum - (i + 1)]) > 0)
        {
            setValue(lf.bwt, i, back(text[seqNum - (i + 1)]));
            setValue(lf.sentinels, i, false);
        }
    }

    /* Compute the rest of the bwt.*/

    // align the chunk_size to underlying word boundaries to prevent parallel write to word spanning chunk boundary
    uint64_t chunkSize = _max((length(sa) / omp_get_max_threads() / 64) * 64, 1ull);
    uint64_t twoPercent = chunkSize / 50;
    // the 0th thread might get an additional chunk because of the above alignment so we count from the 1st instead
    uint32_t countThreadID = omp_get_max_threads() > 1 ? 1 : 0;

    SEQAN_OMP_PRAGMA(parallel for schedule(static, chunkSize))
    for (TSize i = 0; i < length(sa); ++i)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, sa[i], stringSetLimits(text));

        if (getSeqOffset(pos) != 0)
        {
            setValue(lf.bwt, i + seqNum, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
            setValue(lf.sentinels, i + seqNum, false);
        }
        else
        {
            setValue(lf.bwt, i + seqNum, lf.sentinelSubstitute);
            setValue(lf.sentinels, i + seqNum, true);
        }

        if (((static_cast<uint32_t>(omp_get_thread_num()) == countThreadID) && ((i % chunkSize) % twoPercent == 0)))
            progress(((i % chunkSize) / twoPercent) * 2);
    }

   // Update all ranks.
   updateRanks(lf.bwt);
   // Update the auxiliary RankDictionary of sentinel positions.
   updateRanks(lf.sentinels);
}

template <typename TText, typename TSpec, typename TConfig, typename TOtherText, typename TSA, typename TCallback>
void
createLFProgress(LF<TText, TSpec, TConfig> & lf, TOtherText const & text, TSA const & sa, TCallback && progress)
{
    typedef LF<TText, TSpec, TConfig>                          TLF;
    typedef typename Value<TLF>::Type                          TValue;
    typedef typename Size<TLF>::Type                           TSize;

    // Clear assuming undefined state.
    clear(lf);

    // Compute prefix sum.
    prefixSums<TValue>(lf.sums, text);

    // Choose the sentinel substitute.
    _setSentinelSubstitute(lf);

    // Create and index BWT bwt for rank queries.
    createRankDictionaryProgress(lf, text, sa, progress);

    // Add sentinels to prefix sum.
    TSize sentinelsCount = countSequences(text);
    for (TSize i = 0; i < length(lf.sums); ++i)
        lf.sums[i] += sentinelsCount;

    progress(100);
}

#endif // LAMBDA_INDEXER_MISC_HPP_
