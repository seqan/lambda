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

#include "output.hpp" //TODO only needed because options are in one file, remove later
#include "misc.hpp"
#include "options.hpp"
#include "radix_inplace.h"
#include "lambda_indexer_misc.hpp"

using namespace seqan;

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TOrigAlph>
inline int
loadSubjSeqsAndIds(TCDStringSet<String<TOrigAlph>> & originalSeqs,
                   LambdaIndexerOptions const & options)
{
    typedef TCDStringSet<String<char, Alloc<>>>             TIDs;
    typedef TCDStringSet<String<char, Alloc<Truncate_>>>    TIDsTruncated;
    TIDs ids;
    // difference only in name, same layout in mem, so we can hack the type:
    TIDsTruncated* tIds = static_cast<TIDsTruncated*>((void*)&ids);

    double start = sysTime();
    myPrint(options, 1, "Loading Subject Sequences and Ids...");

    SeqFileIn infile(toCString(options.dbFile));
    int ret;
    if (options.truncateIDs)
        ret = myReadRecords(*tIds, originalSeqs, infile);
    else
        ret = myReadRecords(ids, originalSeqs, infile);

    if (ret)
        return ret;

    myPrint(options, 1,  " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n");

    if (length(originalSeqs) == 0)
    {
        std::cerr << "ERROR: No sequences in file. Aborting.\n";
        return -1;
    }
    unsigned long maxLen = 0ul;
    for (auto const & s : originalSeqs)
    {
        if (length(s) > maxLen)
        {
            maxLen = length(s);
        }
        else if (length(s) == 0)
        {
            std::cerr << "ERROR: Unexpectedly encountered a sequence of length 0 in the file."
                      << "Remove the entry and try again. Aborting.\n";
            return -1;
        }
    }
    myPrint(options, 2, "Number of sequences read: ", length(originalSeqs),
            "\nLongest sequence read: ", maxLen, "\n\n");

    if (length(originalSeqs) * 6 >= std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max())
    {
        std::cerr << "ERROR: Too many sequences submitted. The maximum (including frames) is "
                  << std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max()
                  << ".\n";
        return -1;
    }

    if (maxLen >= std::numeric_limits<SizeTypePos_<TOrigAlph>>::max())
    {
        std::cerr << "ERROR: one or more of your subject sequences are too long. "
                  << "The maximum length is " << std::numeric_limits<SizeTypePos_<TOrigAlph>>::max()
                  << ".\n";
        return -1;
    }


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
                       LambdaIndexerOptions const & options)
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

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TOrigAlph>
inline void
translateOrSwap(TCDStringSet<String<TTransAlph>> & out,
                TCDStringSet<String<TOrigAlph>> & in,
                LambdaIndexerOptions const & options)
{
    //TODO more output
    myPrint(options, 1, "translating...");
    translate(out,
              in,
              SIX_FRAME,
              options.geneticCode);
}

template <typename TSameAlph>
inline void
translateOrSwap(TCDStringSet<String<TSameAlph>> & out,
                TCDStringSet<String<TSameAlph>> & in,
                LambdaIndexerOptions const & /**/)
{
    swap(out, in);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph>
inline void
dumpTranslatedSeqs(TCDStringSet<String<TTransAlph>> const & translatedSeqs,
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
// reduceOrSwap(TCDStringSet<String<TRedAlph>> & out,
//              TCDStringSet<String<TTransAlph>> & in)
// {
//     //TODO more output
//     // reduce implicitly
//     myPrint(options, 1, "Reducing...");
//     out.concat = in.concat;
//     out.limits = in.limits;
// }
// 
// template <typename TSameAlph>
// inline void
// reduceOrSwap(TCDStringSet<String<TSameAlph>> & out,
//              TCDStringSet<String<TSameAlph>> & in)
// {
//     swap(out, in);
// }

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TRedAlph>
inline bool
checkIndexSize(TCDStringSet<String<TRedAlph>> const & seqs)
{
    using SAV = typename SAValue<TCDStringSet<String<TRedAlph>>>::Type;
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
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntStarts;
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntEnds;
//     resize(segIntervals, numberOfSeqs, Exact());

    if (options.segFile != "")
    {
        myPrint(options, 1, "Constructing binary seqan masking from seg-file...");

        std::ifstream stream;
        stream.open(toCString(options.segFile));
        if (!stream.is_open())
        {
            std::cerr << "ERROR: could not open seg file.\n";
            return -1;
        }

        auto reader = directionIterator(stream, Input());

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
        while (value(reader) == '>')
        {
//             if (curSeq == end(_segIntervals))
//                 return -7;
            if (curSeq == numberOfSeqs)
            {
                std::cerr << "ERROR: seg file has more entries then database.\n";
                return -7;
            }
            skipLine(reader);
            if (atEnd(reader))
                break;

            unsigned curInt = 0;
            while ((!atEnd(reader)) && (value(reader) != '>'))
            {
                resize(_segIntStarts[curSeq], length(_segIntStarts[curSeq])+1);
                resize(_segIntEnds[curSeq], length(_segIntEnds[curSeq])+1);
                clear(buf);
                readUntil(buf, reader, IsWhitespace());

//                 std::get<0>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntStarts[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);
                skipUntil(reader, IsDigit());

                clear(buf);
                readUntil(buf, reader, IsWhitespace());

//                 std::get<1>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntEnds[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);

//                 appendValue(*curSeq, tup);

                skipLine(reader);
                curInt++;
            }
            if (atEnd(reader))
                break;
            else
                curSeq++;
        }
//         if (curSeq != end(_segIntervals))
//             return -9;
        if (curSeq != (numberOfSeqs - 1))
        {
            std::cerr << "ERROR: seg file has less entries (" << curSeq + 1
                      << ") than database (" << numberOfSeqs << ").\n";
            return -9;
        }

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
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "\n");
    return 0;
}

// --------------------------------------------------------------------------
// Function createSuffixArray()
// --------------------------------------------------------------------------

// If there is no overload with progress function, then strip it
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgo,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  TAlgo const &,
                  TLambda const &)
{
    return createSuffixArray(SA, s, TAlgo());
}

// ----------------------------------------------------------------------------
// Function indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TLambda>
inline bool
indexCreateProgress(Index<TText, FMIndex<TSpec, TConfig> > & index,
                    FibreSALF const &,
                    TLambda const & progressCallback)
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >               TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type            TTempSA;
    typedef typename Size<TIndex>::Type                          TSize;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    TText const & text = indexText(index);

    if (empty(text))
        return false;

    TTempSA tempSA;

    std::cout << "Generating       0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n"
                 " (1) SuffixArray |" << std::flush;
    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    createSuffixArray(tempSA, text, TAlgo(), progressCallback);

    std::cout << " (2) FM-Index..." << std::flush;
    // Create the LF table.
    createLF(indexLF(index), text, tempSA);

    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());

    // Create the compressed SA.
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);
    std::cout << " done.\n" << std::flush;
    return true;
}

template <typename TText, typename TSpec, typename TLambda>
inline bool
indexCreateProgress(Index<TText, IndexSa<TSpec> > & index,
                    FibreSA const &,
                    TLambda const & progressCallback)
{
    typedef Index<TText, IndexSa<TSpec> >                        TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                TSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    TText const & text = indexText(index);

    if (empty(text))
        return false;

    TSA & sa = getFibre(index, FibreSA());

    std::cout << "Generating       0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n"
                 "  SuffixArray    |" << std::flush;
    // Create the full SA.
    resize(sa, lengthSum(text), Exact());
    createSuffixArray(sa, text, TAlgo(), progressCallback);

    return true;
}

// --------------------------------------------------------------------------
// Function generateIndexAndDump()
// --------------------------------------------------------------------------

#ifdef _OPENMP
#define TID omp_get_thread_num()
#else
#define TID 0
#endif

template <typename TIndexSpec,
          typename TIndexSpecSpec,
          typename TString,
          typename TSpec,
          typename TRedAlph_,
          BlastProgram p>
inline void
generateIndexAndDump(StringSet<TString, TSpec>        & seqs,
                     LambdaIndexerOptions       const & options,
                     BlastProgramSelector<p>    const &,
                     TRedAlph_                  const &)
{
    using TTransSeqs    = TCDStringSet<String<TransAlph<p>>>;

    using TRedAlph      = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TRedSeqVirt   = ModifiedString<String<TransAlph<p>, Alloc<>>,
                            ModView<FunctorConvert<TransAlph<p>,TRedAlph>>>;
    using TRedSeqsVirt  = StringSet<TRedSeqVirt, Owner<ConcatDirect<>>>;

    static bool constexpr
    indexIsFM           = std::is_same<TIndexSpec,
                                       TFMIndex<TIndexSpecSpec>
                                       >::value;
    static bool constexpr
    alphReduction       = !std::is_same<TransAlph<p>, TRedAlph>::value;

    using TRedSeqs      = typename std::conditional<
                            !alphReduction,
                            TTransSeqs,             // owner
                            TRedSeqsVirt>::type;    // modview
    using TRedSeqsACT   = typename std::conditional<
                            !alphReduction,
                            TTransSeqs &,           // reference to owner
                            TRedSeqsVirt>::type;    // modview

    using TDbIndex      = Index<TRedSeqs, TIndexSpec>;
    using TFullFibre    = typename std::conditional<indexIsFM,
                                                    FibreSALF,
                                                    FibreSA>::type;
    static bool constexpr
    hasProgress         = std::is_same<TIndexSpecSpec, RadixSortSACreateTag>::value;

    // Generate Index
    if (!hasProgress)
        myPrint(options, 1, "Generating Index...");

    double s = sysTime();

//     std::cout << "indexIsFM: " << int(indexIsFM) << std::endl;

    // FM-Index needs reverse input
    if (indexIsFM)
        reverse(seqs);

    TRedSeqsACT redSubjSeqs(seqs);

    TDbIndex dbIndex(redSubjSeqs);

    // instantiate SA
    if (hasProgress && (options.verbosity >= 1))
    {
        uint64_t _lastPercent = 0;
        indexCreateProgress(dbIndex, TFullFibre(),
                            [&_lastPercent] (uint64_t curPerc)
                            {
                                SEQAN_OMP_PRAGMA(critical(progressBar))
        //                         if (TID == 0)
                                printProgressBar(_lastPercent, curPerc);
                            });
    }
    else
    {
        indexCreate(dbIndex, TFullFibre());
    }

    // since we dumped unreduced sequences before and reduced sequences are
    // only "virtual" we clear them before dump
    if (alphReduction || indexIsFM)
        clear(seqs);
    if (alphReduction)
        clear(redSubjSeqs.limits); // limits part is not lightweight

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
