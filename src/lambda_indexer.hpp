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
// lambda_indexer.hpp: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_LAMBDA_LAMBDA_INDEXER_H_
#define SEQAN_LAMBDA_LAMBDA_INDEXER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/seq_io.h>
#include <seqan/index.h>
// #include <seqan/index_extras.h>
#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include "misc.hpp"
#include "options.hpp"
#include "trans.hpp"
#include "alph.hpp"



using namespace seqan;


// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TOrigAlph>
inline int
loadSubjSeqsAndIds(TCDStringSet<TOrigAlph> & originalSeqs,
                   LambdaIndexerOptions const & options)
{
    int ret = 0;
    StringSet<CharString, Owner<ConcatDirect<>>> ids;

    double start = sysTime();
    myPrint(options, 1, "Loading Subject Sequences and Ids…");


    if (options.fileFormat)
        ret = loadSeqsAndIds(ids,
                             originalSeqs,
                             options.dbFile,
                             Fastq());
    else
        ret = loadSeqsAndIds(ids,
                             originalSeqs,
                             options.dbFile,
                             Fasta());
    if (ret)
        return ret;

    myPrint(options, 1,  " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n");

    unsigned long maxLen = 0ul;
    for (auto const & s : originalSeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    myPrint(options, 2, "Number of sequences read: ", length(originalSeqs),
            "\nLongest sequence read: ", maxLen, "\n\n");

    myPrint(options, 1, "Dumping Subj Ids...");

    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".ids");
    save(ids, toCString(_path));

    myPrint(options, 1, " done.\n");
    finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TLimits>
inline void
_saveOriginalSeqLengths(TLimits limits, // we want copy!
                       LambdaIndexerOptions const & options,
                       True const & /*SHasFrames*/)
{
    for (uint32_t i = 0; i < (length(limits) - 1); ++i)
        limits[i] = limits[i+1] - limits[i];
    // last entry not overwritten, should be the sum of all lengths

    myPrint(options, 1, " dumping untranslated subject lengths...");
    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".untranslengths");
    save(limits, toCString(_path));
}

template <typename TLimits>
inline void
_saveOriginalSeqLengths(TLimits const &/**/,
                       LambdaIndexerOptions const & /**/,
                       False const & /*SHasFrames*/)
{
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TOrigAlph>
inline void
translateOrSwap(TCDStringSet<TTransAlph> & out,
                TCDStringSet<TOrigAlph> & in,
                LambdaIndexerOptions const & options)
{
    //TODO more output
    myPrint(options, 1, "translating…");
    translate(out,
              in,
              SIX_FRAME,
              options.geneticCode);
}

template <typename TSameAlph>
inline void
translateOrSwap(TCDStringSet<TSameAlph> & out,
                TCDStringSet<TSameAlph> & in,
                LambdaIndexerOptions const & /**/)
{
    swap(out, in);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph>
inline void
dumpTranslatedSeqs(TCDStringSet<TTransAlph> const & translatedSeqs,
                   LambdaIndexerOptions const & options)
{
    double start = sysTime();
    myPrint(options, 1, "Dumping unreduced Subj Sequences...");

    //TODO save to TMPDIR instead
    std::string _path = options.dbFile + '.' + std::string(_alphName(TTransAlph()));
    save(translatedSeqs, _path.c_str());

    myPrint(options, 1, " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

// template <typename TTransAlph, typename TRedAlph>
// inline void
// reduceOrSwap(TCDStringSet<TRedAlph> & out,
//              TCDStringSet<TTransAlph> & in)
// {
//     //TODO more output
//     // reduce implicitly
//     myPrint(options, 1, "Reducing…");
//     out.concat = in.concat;
//     out.limits = in.limits;
// }
// 
// template <typename TSameAlph>
// inline void
// reduceOrSwap(TCDStringSet<TSameAlph> & out,
//              TCDStringSet<TSameAlph> & in)
// {
//     swap(out, in);
// }

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TRedAlph>
inline bool
checkIndexSize(TCDStringSet<TRedAlph> const & seqs)
{
    using SAV = typename SAValue<TCDStringSet<TRedAlph>>::Type;
    uint64_t curNumSeq = length(seqs);
    uint64_t maxNumSeq = std::numeric_limits<typename Value<SAV, 1>::Type>::max();

    if (curNumSeq >= maxNumSeq)
    {
        std::cerr << "Too many sequences to be indexed:\n  "
                  << length(seqs) << " in file, but only "
                  << maxNumSeq << " supported by index.\n";
        return false;
    }

    uint64_t maxLenSeq = std::numeric_limits<typename Value<SAV, 2>::Type>::max();
    uint64_t maxLen = 0ul;
    for (auto const & s : seqs)
        if (length(s) > maxLen)
            maxLen = length(s);

    if (maxLen >= maxLenSeq)
    {
        std::cerr << "Too long sequences to be indexed:\n  "
                  << "length" << maxLen << " present in file, but only "
                  << maxLenSeq << " supported by index.\n";
        return false;
    }
    return true;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

inline int
convertMaskingFile(uint64_t numberOfSeqs,
                   LambdaIndexerOptions const & options)

{
    int ret = 0;
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntStarts;
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntEnds;
//     resize(segIntervals, numberOfSeqs, Exact());

    if (options.segFile != "")
    {
        myPrint(options, 1, "Constructing binary seqan masking from seg-file...");

        std::ifstream stream;
        stream.open(toCString(options.segFile));
        if (!stream.is_open())
            return -1;

        typedef RecordReader<std::ifstream, SinglePass<> > TReader;
        TReader reader(stream);

//         StringSet<String<Tuple<unsigned, 2>>> _segIntervals;
//         auto & _segIntervals = segIntervals;
//         resize(_segIntervals, numberOfSeqs, Exact());
        StringSet<String<unsigned>> _segIntStarts;
        StringSet<String<unsigned>> _segIntEnds;
        resize(_segIntStarts, numberOfSeqs, Exact());
        resize(_segIntEnds, numberOfSeqs, Exact());
        CharString buf;
//         std::tuple<unsigned, unsigned> tup;

//         auto curSeq = begin(_segIntervals);
        unsigned curSeq = 0;
        while ((!atEnd(reader)) && (value(reader) == '>'))
        {
//             if (curSeq == end(_segIntervals))
//                 return -7;
            if (curSeq == numberOfSeqs)
                return -7;
            ret = skipLine(reader);
            if ((ret) && (ret != EOF_BEFORE_SUCCESS))
                return ret;

            unsigned curInt = 0;
            while ((!atEnd(reader)) && (value(reader) != '>'))
            {
                resize(_segIntStarts[curSeq], length(_segIntStarts[curSeq])+1);
                resize(_segIntEnds[curSeq], length(_segIntEnds[curSeq])+1);
                clear(buf);
                ret = readDigits(buf, reader);
                if (ret)
                    return ret;

//                 std::get<0>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntStarts[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);
                ret = skipNChars(reader, 3);
                if (ret)
                    return ret;

                clear(buf);
                ret = readDigits(buf, reader);
                if (ret)
                    return ret;

//                 std::get<1>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntEnds[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);

//                 appendValue(*curSeq, tup);

                ret = skipLine(reader);
                if ((ret) && (ret != EOF_BEFORE_SUCCESS))
                    return ret;
                curInt++;
            }
            curSeq++;
        }
//         if (curSeq != end(_segIntervals))
//             return -9;
        if (curSeq != numberOfSeqs)
            return -9;

        segIntStarts.concat = concat(_segIntStarts);
        segIntStarts.limits = stringSetLimits(_segIntStarts);
        segIntEnds.concat = concat(_segIntEnds);
        segIntEnds.limits = stringSetLimits(_segIntEnds);
//         segIntEnds = _segIntEnds;
//         segIntervals = _segIntervals; // non-concatdirect to concatdirect

        stream.close();

    } else
    {
        myPrint(options, 1, "No Seg-File specified, no masking will take place.\n");
//         resize(segIntervals, numberOfSeqs, Exact());
        resize(segIntStarts, numberOfSeqs, Exact());
        resize(segIntEnds, numberOfSeqs, Exact());
    }

//     for (unsigned u = 0; u < length(segIntStarts); ++u)
//     {
//         myPrint(options, 1,u, ": ";
//         for (unsigned v = 0; v < length(segIntStarts[u]); ++v)
//         {
//             myPrint(options, 1,'(', segIntStarts[u][v], ", ", segIntEnds[u][v], ")  ";
//         }
//         myPrint(options, 1,'\n';
//     }
    myPrint(options, 1, "Dumping binary seqan mask file...");
    CharString _path = options.dbFile;
    append(_path, ".binseg_s");
    save(segIntStarts, toCString(_path));
    _path = options.dbFile;
    append(_path, ".binseg_e");
    save(segIntEnds, toCString(_path));
    myPrint(options, 1, " done.\n\n");

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

#ifdef _OPENMP
#define TID omp_get_thread_num()
#else
#define TID 0
#endif

template <typename TIndex, typename TText, typename TFibre>
inline void
createIndexActual(TIndex & index,
                      TText const & text,
                      TFibre const &,
                      SaAdvancedSort<MergeSortTag> const &)
{
    ComparisonCounter<TText, std::true_type> counter(text);
    indexCreate(index, text, TFibre(), [&counter] () { counter.inc(); });
    printProgressBar(counter._lastPercent, 100);
}

template <typename TIndex, typename TText, typename TFibre>
inline void
createIndexActual(TIndex & index,
                      TText const & text,
                      TFibre const &,
                      SaAdvancedSort<QuickSortBucketTag> const &)
{
    uint64_t _lastPercent = 0;
    indexCreate(index, text, TFibre(),
        [&_lastPercent] (uint64_t curPerc)
        {
            if (TID == 0)
                printProgressBar(_lastPercent, curPerc);
        });
    printProgressBar(_lastPercent, 100);
}

template <typename TIndex, typename TText, typename TFibre, typename TAlgo>
inline void
createIndexActual(TIndex & index,
                      TText const & text,
                      TFibre const &,
                      TAlgo const &)
{
    indexCreate(index, text, TFibre(), [] () {});
}

template <typename TIndexSpec,
          typename TIndexSpecSpec,
          typename TString,
          typename TSpec,
          typename TRedAlph_,
          BlastFormatProgram p>
inline void
generateIndexAndDump(StringSet<TString, TSpec> & seqs,
                     LambdaIndexerOptions const & options,
                     TRedAlph_ const & /**/,
                     BlastFormat<BlastFormatFile::INVALID_File,p,
                       BlastFormatGeneration::INVALID_Generation> const & /**/)
{
//     using TFormat   = BlastFormat<BlastFormatFile::INVALID_File,
//                                   p,
//                                   BlastFormatGeneration::INVALID_Generation>;

    using TTransSeqs    = TCDStringSet<TransAlph<p>>;

    using TRedAlph      = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TRedSeqVirt   = ModifiedString<String<TransAlph<p>, PackSpec>,
                            ModView<FunctorConvert<TransAlph<p>,TRedAlph>>>;
    using TRedSeqsVirt  = StringSet<TRedSeqVirt, Owner<ConcatDirect<>>>;

    static bool constexpr
    indexIsFM           = std::is_same<TIndexSpec,
                                       TFMIndex<TIndexSpecSpec>
                                       >::value;
    static bool constexpr
    noReduction         = std::is_same<TransAlph<p>, TRedAlph>::value;

    using TRedSeqs      = typename std::conditional<
                            noReduction,
                            TTransSeqs,             // owner
                            TRedSeqsVirt>::type;    // modview
    using TRedSeqsACT   = typename std::conditional<
                            noReduction,
                            TTransSeqs &,           // reference to owner
                            TRedSeqsVirt>::type;    // modview

    using TDbIndex      = Index<TRedSeqs, TIndexSpec>;
    using TFullFibre    = typename std::conditional<indexIsFM,
                                                    FibreSALF,
                                                    FibreSA>::type;
    static bool constexpr
    hasProgress         = std::is_same<TIndexSpecSpec,
                                       SaAdvancedSort<QuickSortBucketTag>>::value
                          || std::is_same<TIndexSpecSpec,
                                       SaAdvancedSort<MergeSortTag>>::value;

//     using TCountPartial = std::is_same<TIndexSpecSpec,
//                                        SaAdvancedSort<MergeSortTag>>;
// TODO debug this
//     using TProgressCounter = typename std::conditional<
//                                 hasProgress,
//                                 ComparisonCounter<TRedSeqs, TCountPartial>,
//                                 ComparisonCounter<TRedSeqs, Nothing>>::type;
//     using TProgressCounter = ComparisonCounter<TRedSeqs, Nothing>;

    // Generate Index
    myPrint(options, 1, "Generating Index...");
    double s = sysTime();

//     std::cout << "indexIsFM: " << int(indexIsFM) << std::endl;

    // FM-Index needs reverse input
    if (indexIsFM)
        reverse(seqs);

    TRedSeqsACT redSubjSeqs(seqs);

//     TProgressCounter counter(redSubjSeqs, 0);
//     std::cout << "ExpectedNumComparisons: " << counter._expectedComparisons
//               << std::endl;
    if (hasProgress)
        myPrint(options, 1, "progress:\n"
                "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n|");
    TDbIndex dbIndex(redSubjSeqs);
    // instantiate SA

    // create SA with progressCallback function
    if (options.verbosity >= 1)
        createIndexActual(dbIndex, redSubjSeqs, TFullFibre(),
                          TIndexSpecSpec());
    else // don't print progress (independent of algo)
        createIndexActual(dbIndex, redSubjSeqs, TFullFibre(),
                         Nothing());

    // instantiate potential rest
//     std::cout << "\nActualNumComparisons: " << counter._comparisons
//               << std::endl;
//     indexRequire(dbIndex, TFullFibre()); // instantiate

    // since we dumped unreduced sequences before and reduced sequences are
    // only "virtual" we clear them before dump
    clear(seqs);
    if (!noReduction)
        clear(redSubjSeqs.limits);

    double e = sysTime() - s;
    if (!hasProgress)
        myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n\n");

    // Dump Index
    myPrint(options, 1, "Writing Index to disk...");
    s = sysTime();
    std::string path = toCString(options.dbFile);
    path += '.' + std::string(_alphName(TRedAlph()));
    if (indexIsFM)
        path += ".fm";
    else
        path += ".sa";
    save(dbIndex, path.c_str());
    e = sysTime() - s;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n");
}

#endif // header guard
