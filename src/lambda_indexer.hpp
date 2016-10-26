// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2016, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016, Knut Reinert and Freie Universität Berlin
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
                   std::unordered_map<std::string, uint64_t> & accToIdRank,
                   LambdaIndexerOptions const & options)
{
    typedef TCDStringSet<String<char, Alloc<>>>             TIDs;

    TIDs ids; // the IDs

    // see http://www.uniprot.org/help/accession_numbers
    std::regex const accRegEx{"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"};

    // lambda that truncates IDs at first whitespace
    auto truncateID = [] (auto && id, uint64_t const)
    {
        IsWhitespace isWhitespace;
        for (size_t i = 0; i < length(id); ++i)
        {
            if (isWhitespace(id[i]))
            {
                resize(id, i);
                break;
            }
        }
    };

    // lambda that extracts accession numbers and saves them in the map
    auto extractAccIds = [&accToIdRank, &accRegEx] (auto && id, uint64_t const rank)
    {
        // TODO avoid copying here by specializing regex_iterator
        std::string buf;
        assign(buf, id);

        for (auto it = std::sregex_iterator(buf.begin(), buf.end(), accRegEx), itEnd = std::sregex_iterator();
             it != itEnd;
             ++it)
        {
            SEQAN_ASSERT_MSG(accToIdRank.count(it->str()) == 0,
                             "An accession number appeared twice in the file, but they should be unique.");
            accToIdRank[it->str()] = rank;
        }
    };

    double start = sysTime();
    myPrint(options, 1, "Loading Subject Sequences and Ids...");

    SeqFileIn infile(toCString(options.dbFile));
    int ret;
    if (options.truncateIDs)
    {
        if (options.hasSTaxIds)
        {
            ret = myReadRecords(ids, originalSeqs, infile, [&] (auto && id, uint64_t const rank)
            {
                extractAccIds(std::forward<decltype(id)>(id), rank);
                truncateID(std::forward<decltype(id)>(id), rank);
            });
        } else
        {
            ret = myReadRecords(ids, originalSeqs, infile, truncateID);
        }
    } else
    {
        if (options.hasSTaxIds)
            ret = myReadRecords(ids, originalSeqs, infile, extractAccIds);
        else
            ret = myReadRecords(ids, originalSeqs, infile);
    }

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
    CharString _path = options.indexDir;
    append(_path, "/seq_ids");
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
    CharString _path = options.indexDir;
    append(_path, "/untranslated_seq_lengths");
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
    std::string _path = options.indexDir + "/translated_seqs";
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
// Function mapAndDumpTaxIDs()
// --------------------------------------------------------------------------

inline int
mapAndDumpTaxIDs(std::unordered_map<std::string, uint64_t>       const & accToIdRank,
                 uint64_t                                        const   numSubjects,
                 LambdaIndexerOptions                            const & options)

{
//     StringSet<String<uint64_t>> sTaxIds; // not concat because we resize inbetween
//     resize(sTaxIds, numSubjects);
//
//     // TODO
//     // load accToTaxId mapAndDumpTaxIDs
//     // on every line check if the col1 exists in accToIdRank
//     //      then appendValue(sTaxIds[accToIdRank[col1]], col3);
//
//     // Do something with the subjects that have no (valid) taxid?
//
//     // concat direct so that it's easier to read/write
//     StringSet<String<uint64_t>, Owner<ConcatDirect<>>> outSTaxIds = sTaxIds;
//     save(options.indexDir + "/staxids", outSTaxIds);
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
                  TLambda &&)
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
                    TLambda && progressCallback)
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

template <typename TText, typename TSpec, typename TConfig, typename TLambda>
inline bool
indexCreateProgress(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > & index,
                    FibreSALF const &,
                    TLambda && progressCallback)
{
    auto progressCallback2 = progressCallback; // need second lambda because counter internally increased

    std::cout << "Bi-Directional Index [forward]\n";
    bool ret = indexCreateProgress(index.fwd, FibreSALF(), progressCallback);
    if (!ret)
        return ret;

    std::cout << "Bi-Directional Index [backward]\n";
    return indexCreateProgress(index.rev, FibreSALF(), progressCallback2);
}

template <typename TText, typename TSpec, typename TLambda>
inline bool
indexCreateProgress(Index<TText, IndexSa<TSpec> > & index,
                    FibreSA const &,
                    TLambda && progressCallback)
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
          typename TDirection,
          BlastProgram p>
inline void
generateIndexAndDump(StringSet<TString, TSpec>        & seqs,
                     LambdaIndexerOptions       const & options,
                     BlastProgramSelector<p>    const &,
                     TRedAlph_                  const &,
                     Tag<TDirection>            const &)
{
    using TTransSeqs    = TCDStringSet<String<TransAlph<p>>>;

    using TRedAlph      = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TRedSeqVirt   = ModifiedString<String<TransAlph<p>, Alloc<>>,
                            ModView<FunctorConvert<TransAlph<p>,TRedAlph>>>;
    using TRedSeqsVirt  = StringSet<TRedSeqVirt, Owner<ConcatDirect<>>>;

    static bool constexpr
    indexIsFM           = std::is_same<TIndexSpec, TFMIndex<TIndexSpecSpec> >::value ||
                          std::is_same<TIndexSpec, BidirectionalIndex<TFMIndex<TIndexSpecSpec> > >::value;
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
    if (indexIsFM && std::is_same<Tag<TDirection>, Fwd>::value)
        reverse(seqs);

    TRedSeqsACT redSubjSeqs(seqs);

    TDbIndex dbIndex(redSubjSeqs);

    // instantiate SA
    if (hasProgress && (options.verbosity >= 1))
    {
        uint64_t _lastPercent = 0;
        indexCreateProgress(dbIndex, TFullFibre(),
                            [_lastPercent] (uint64_t curPerc) mutable
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
    std::decay_t<decltype(redSubjSeqs.limits)> tmpLimits;
    if (alphReduction || !indexIsFM) // fm indexes don't dump them anyways
    {
        if (indexIsFM && (std::is_same<Tag<TDirection>, Rev>::value))
        {
            // these makes redSubjSeqs appear empty and deactivates output
            swap(tmpLimits, redSubjSeqs.limits);
        } else
        {
            clear(seqs);
            clear(redSubjSeqs.limits); // limits part is not lightweight
        }
    }

    double e = sysTime() - s;
    if (!hasProgress)
        myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n\n");

    // Dump Index
    myPrint(options, 1, "Writing Index to disk...");
    s = sysTime();
    std::string path = options.indexDir + "/index";

    if (std::is_same<Tag<TDirection>, Rev>::value)
        path += ".rev";

    save(dbIndex, path.c_str());

    e = sysTime() - s;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n");

    if (alphReduction && indexIsFM && (std::is_same<Tag<TDirection>, Rev>::value))
    {
        // we swap back so that the sequences can be used for building the second index
        swap(tmpLimits, redSubjSeqs.limits);
//         redSubjSeqs.limits = tmpLimits;
    }
}

#endif // header guard
