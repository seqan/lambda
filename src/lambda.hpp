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
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// lambda.hpp: contains the main progam pipeline
// ==========================================================================


#ifndef SEQAN_LAMBDA_LAMBDA_H_
#define SEQAN_LAMBDA_LAMBDA_H_

#include <type_traits>
#include <iomanip>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
// #include <seqan/index_extras.h>

#include <seqan/blast.h>

#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include <seqan/align_extend.h>

#include "options.hpp"
#include "match.hpp"
#include "misc.hpp"
#include "trans.hpp"
#include "alph.hpp"
#include "holders.hpp"

using namespace seqan;


enum COMPUTERESULT_
{
    SUCCESS = 0,
    PREALIGNSCORE,
    PREALIGNEVAL,
    ALIGNSCORE,
    ALIGNEVAL
};


// comparison operator to sort SA-Values based on the strings in the SA they refer to
template <typename TSav, typename TStringSet>
struct Comp :
    public ::std::binary_function < TSav, TSav, bool >
{
//     TSeedIndex const & index;
//     typedef typename Fibre<TSeedIndex, EsaSA>::Type TSa;
//      typedef typename SAValue<TSa>::Type TSav;

    TStringSet const & stringSet;

    Comp(TStringSet const & _stringSet)
        : stringSet(_stringSet)
    {}

    inline bool operator() (TSav const & i, TSav const & j) const
    {
//        return false;
         return (value(stringSet,getSeqNo(i)) < value(stringSet,getSeqNo(j)));
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <typename TGlobalHolder>
inline int
loadDbIndexFromDisk(TGlobalHolder       & globalHolder,
                    LambdaOptions const & options)
{
    std::cout << "Loading Database Index from disk…" << std::flush;
    double start = sysTime();
    int ret = open(globalHolder.dbIndex, toCString(options.dbFile));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    double finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "No of Fibres: " << length(indexSA(globalHolder.dbIndex))
              << "\t no of Seqs in Db: " << length(indexText(globalHolder.dbIndex))
              << "\n\n";

    return 0;
}


// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------


// Generic, with translation and reduction
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          MyEnableIf<!std::is_same<TRedAlph,AminoAcid>::value &&
          !std::is_same<TRedAlph,Dna5>::value> = 0>
inline int
loadQueryImpl(GlobalDataHolder<TRedAlph,
                                TScoreScheme,
                                m, p, g>         & globalHolder,
              LambdaOptions                     const & options)
{
    StringSet<String<Dna5>, Owner<ConcatDirect<>>> untranslatedSeqs;

    int ret = 0;
    if (options.fileFormat)
        ret = loadSeqsAndIds(globalHolder.qryIds,
                             untranslatedSeqs,
                             options.queryFile,
                             Fastq());
    else
        ret = loadSeqsAndIds(globalHolder.qryIds,
                             untranslatedSeqs,
                             options.queryFile,
                             Fasta());
    if (ret)
        return ret;

    std::cout << "translating…" << std::flush;
    translate(globalHolder.qrySeqs,
              untranslatedSeqs,
              SIX_FRAME,
              options.geneticCode);

    // preserve lengths of untranslated sequences
    resize(globalHolder.untransQrySeqLengths,
           length(untranslatedSeqs.limits),
           Exact());
    for (uint32_t i = 0;
         i < (length(globalHolder.untransQrySeqLengths) - 1);
         ++i)
    {
        globalHolder.untransQrySeqLengths[i] =
            untranslatedSeqs.limits[i + 1] - untranslatedSeqs.limits[i];
    }
    // save sum of lengths (both strings have n + 1 elements
    back(untranslatedSeqs.limits) = length(untranslatedSeqs.concat);

    // reduce implicitly
    std::cout << "reducing…" << std::flush;
    globalHolder.redQrySeqs.concat = globalHolder.qrySeqs.concat;
    globalHolder.redQrySeqs.limits = globalHolder.qrySeqs.limits;

    return 0;
}

// only translation
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          MyEnableIf<std::is_same<TRedAlph,AminoAcid>::value> = 0 >
inline int
loadQueryImpl(GlobalDataHolder<TRedAlph,
                                TScoreScheme,
                                m, p, g>         & globalHolder,
              LambdaOptions                     const & options)
{
    StringSet<String<Dna5>, Owner<ConcatDirect<>>> untranslatedSeqs;

    int ret = 0;
    if (options.fileFormat)
        ret = loadSeqsAndIds(globalHolder.qryIds,
                             untranslatedSeqs,
                             options.queryFile,
                             Fastq());
    else
        ret = loadSeqsAndIds(globalHolder.qryIds,
                             untranslatedSeqs,
                             options.queryFile,
                             Fasta());
    if (ret)
        return ret;

    std::cout << "translating…" << std::flush;
    translate(globalHolder.qrySeqs,
              untranslatedSeqs,
              SIX_FRAME,
              options.geneticCode);

    return 0;
}

// none
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          MyEnableIf<std::is_same<TRedAlph,Dna5>::value> = 0>
inline int
loadQueryImpl(GlobalDataHolder<TRedAlph,
                                TScoreScheme,
                                m, p, g>         & globalHolder,
              LambdaOptions                     const & options)
{
    if (options.fileFormat)
        return loadSeqsAndIds(globalHolder.qryIds,
                              globalHolder.qrySeqs,
                              options.queryFile,
                              Fastq());
    else
        return loadSeqsAndIds(globalHolder.qryIds,
                              globalHolder.qrySeqs,
                              options.queryFile,
                              Fasta());
}



template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
loadQuery(GlobalDataHolder<TRedAlph,
                            TScoreScheme,
                            m, p, g>         & globalHolder,
          LambdaOptions                     const & options)
{
    double start = sysTime();
    std::cout << "Loading Query Sequences and Ids…" << std::flush;
    int ret = loadQueryImpl(globalHolder, options);
    if (ret)
        return ret;

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    unsigned long maxLen = 0ul;
    for (auto const & s : globalHolder.qrySeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(globalHolder.qrySeqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubjects()
// --------------------------------------------------------------------------


template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
loadSubjects(GlobalDataHolder<TRedAlph,
                              TScoreScheme,
                              m, p, g>         & globalHolder,
             LambdaOptions                    const & options)
{
//     typedef BlastFormat<m,p,g> TFormat;

    double start = sysTime();
    std::cout << "Loading Subj Sequences…" << std::flush;

    CharString _dbSeqs = options.dbFile;
    if (options.alphReduction > 0)
        append(_dbSeqs, ".unredsubj"); // get unreduced stringset
    else // stringset already dumped by index dump
        append(_dbSeqs, ".txt");
    int ret = open(globalHolder.subjSeqs, toCString(_dbSeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Amount: " << length(globalHolder.subjSeqs) << "\n\n"<< std::flush;


    start = sysTime();
    std::cout << "Loading Subj Ids…" << std::flush;
    _dbSeqs = options.dbFile;
    append(_dbSeqs, ".ids");
    ret = open(globalHolder.subjIds, toCString(_dbSeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    globalHolder.dbSpecs.dbName = options.dbFile;
    // if subjects where translated, we don't have the untranslated seqs at all
    // but we still need the data for statistics and position un-translation
    if (sHasFrames(p))
    {
        start = sysTime();
        std::cout << "Loading Lengths of untranslated Subj sequences…" << std::flush;
        _dbSeqs = options.dbFile;
        append(_dbSeqs, ".untranslengths");
        ret = open(globalHolder.untransSubjSeqLengths, toCString(_dbSeqs));
        if (ret != true)
        {
            std::cout << " failed.\n" << std::flush;
            return 1;
        }
        std::cout << " done.\n";
        finish = sysTime() - start;
        std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

        // last value has sum of lengths
        globalHolder.dbSpecs.dbTotalLength =
            back(globalHolder.untransSubjSeqLengths);
        globalHolder.dbSpecs.dbNumberOfSeqs =
            length(globalHolder.untransSubjSeqLengths) - 1;
    } else
    {
        globalHolder.dbSpecs.dbTotalLength =
            length(concat(globalHolder.subjSeqs));
        globalHolder.dbSpecs.dbNumberOfSeqs =
            length(globalHolder.subjSeqs);
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSegintervals()
// --------------------------------------------------------------------------


template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
loadSegintervals(GlobalDataHolder<TRedAlph,
                              TScoreScheme,
                              m, p, g>         & globalHolder,
             LambdaOptions                    const & options)
{

    double start = sysTime();
    std::cout << "Loading Database Masking file…" << std::flush;

    CharString segFileS = options.dbFile;
    append(segFileS, ".binseg_s.concat");
    CharString segFileE = options.dbFile;
    append(segFileE, ".binseg_e.concat");
    struct stat buffer;
    // file exists
    if ((stat(toCString(segFileS), &buffer) == 0) &&
        (stat(toCString(segFileE), &buffer) == 0))
    {
        //cut off ".concat" again
        resize(segFileS, length(segFileS) - 7);
        resize(segFileE, length(segFileE) - 7);

        int ret = open(globalHolder.segIntStarts, toCString(segFileS));
        if (ret != true)
            return -89;
        ret = open(globalHolder.segIntEnds, toCString(segFileE));
        if (ret != true)
            return -89;

    } else
    {
        std::cout << "ERROR: Masking files not found.\n";
        return -89;
    }

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    return 0;
}


// --------------------------------------------------------------------------
// Function prepareScoring()
// --------------------------------------------------------------------------


template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline void
prepareScoringMore(GlobalDataHolder<TRedAlph,
                                    TScoreScheme,
                                    m, p, g>         & globalHolder,
                   LambdaOptions                    const & options,
                   std::true_type                   const & /**/)
{
    setScoreMatch(globalHolder.scoreScheme, options.match);
    setScoreMismatch(globalHolder.scoreScheme, options.misMatch);
}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline void
prepareScoringMore(GlobalDataHolder<TRedAlph,
                                    TScoreScheme,
                                    m, p, g>              & /*globalHolder*/,
                   LambdaOptions                    const & /*options*/,
                   std::false_type                  const & /**/)
{
}


template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
prepareScoring(GlobalDataHolder<TRedAlph, TScoreScheme, m, p, g>
                                                      & globalHolder,
               LambdaOptions                    const & options)
{

    setScoreGapOpen  (globalHolder.scoreScheme, options.gapOpen);
    setScoreGapExtend(globalHolder.scoreScheme, options.gapExtend);
    blastScoringScheme2seqanScoringScheme(globalHolder.scoreScheme);

    prepareScoringMore(globalHolder, options,
                       std::is_same<TScoreScheme, Score<int, Simple>>());

    if (!assignScoreScheme(globalHolder.blastScoringAdapter,
                      globalHolder.scoreScheme))
    {
        ::std::cerr << "Could not computer Karlin-Altschul-Values for "
                    << "Scoring Scheme. Exiting.\n";
        return -1;
    }
    return 0;
}


/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function generateSeeds()
// --------------------------------------------------------------------------
/*
#define THREADLINE std::cout << "\0338" << std::endl << "Thread " << lH.i \
<< std::endl; \
for (unsigned char i=0; i< lH.i*20; ++i) std::cout << std::endl;*/

template <typename TLocalHolder>
inline int
generateSeeds(TLocalHolder & lH)
{
    lH.statusStr << "Generating Seeds…";
    if (lH.options.isTerminal)
        myPrint(lH.options, 2, lH.statusStr);

    double start = sysTime();
    for (unsigned long i = lH.indexBeginQry; i < lH.indexEndQry; ++i)
    {
        for (unsigned j = 0;
             (j* lH.options.seedOffset + lH.options.seedLength)
                <= length(value(lH.gH.redQrySeqs, i));
             ++j)
        {
            appendValue(lH.seeds, infix(value(lH.gH.redQrySeqs, i),
                                     j* lH.options.seedOffset,
                                     j* lH.options.seedOffset
                                     + lH.options.seedLength),
                        Generous());
            appendValue(lH.seedRefs,  i, Generous());
            appendValue(lH.seedRanks, j, Generous());

//             std::cout << "seed: " << back(lH.seeds) << "\n";
        }
    }
    double finish = sysTime() - start;
    lH.statusStr << " done. " ;
    if (lH.options.verbosity > 2)
        lH.statusStr << finish << "[s]. "
                     << length(lH.seeds) << " seeds created.";

    myPrint(lH.options, 2, lH.statusStr);
    return 0;
}


// --------------------------------------------------------------------------
// Function generateTrieOverSeeds()
// --------------------------------------------------------------------------


template <typename TLocalHolder>
inline int
generateTrieOverSeeds(TLocalHolder & lH)
{
    lH.statusStr << "Generating Query-Index…";
    if (lH.options.isTerminal)
        myPrint(lH.options, 2, lH.statusStr);

    double start = sysTime();
    // we only want full length seed sequences in index, build up manually
    typedef typename Fibre<typename TLocalHolder::TSeedIndex, EsaSA>::Type TSa;
    //TODO maybe swap here instead
    lH.seedIndex = decltype(lH.seedIndex)(lH.seeds);
    TSa & sa = indexSA(lH.seedIndex);
    resize(sa, length(lH.seeds));
    for (unsigned u = 0; u < length(lH.seeds); ++u)
    {
        assignValueI1(value(sa,u), u);
        assignValueI2(value(sa,u), 0);
    }
    Comp<typename Value<TSa>::Type, typename TLocalHolder::TSeeds const> comp(lH.seeds);
    std::sort(begin(sa, Standard()), end(sa, Standard()), comp);
    typename Iterator<typename TLocalHolder::TSeedIndex, TopDown<> >::Type it(lH.seedIndex); // instantiate
    double finish = sysTime() - start;

    lH.statusStr << " done. " ;
    if (lH.options.verbosity > 2)
        lH.statusStr << finish << "[s]. "
                     << length(sa) << " fibres in SeedIndex. ";
    myPrint(lH.options, 2, lH.statusStr);

    return 0;
}

// --------------------------------------------------------------------------
// Function search()
// --------------------------------------------------------------------------


template <typename TLocalHolder>
inline int
search(TLocalHolder & lH)
{
    // FIND
    lH.statusStr << "Seeding…";
    if (lH.options.isTerminal)
        myPrint(lH.options, 2, lH.statusStr);

    double start = sysTime();
    if (lH.options.hammingOnly)
    {
        typedef Backtracking<HammingDistance>           BackSpec;
        typedef Finder2<decltype(lH.gH.dbIndex),
                        decltype(lH.seedIndex),
                        BackSpec>                       LambdaFinder;

        LambdaFinder finder;

//         (void*)lH.gH.dbIndex;
        find(finder, lH.gH.dbIndex, lH.seedIndex, lH.options.maxSeedDist, lH);
    }
    else
    {
        typedef Backtracking<EditDistance>              BackSpec;
        typedef Finder2<decltype(lH.gH.dbIndex),
                        decltype(lH.seedIndex),
                        BackSpec>                       LambdaFinder;

        LambdaFinder finder;

        find(finder, lH.gH.dbIndex, lH.seedIndex, lH.options.maxSeedDist, lH);
    }

    double finish = sysTime() - start;

    lH.statusStr << " done. " << finish << "[s]. ";
    myPrint(lH.options, 2, lH.statusStr);

    return 0;
}


// --------------------------------------------------------------------------
// Function joinAndFilterMatches()
// --------------------------------------------------------------------------


template <typename TLocalHolder>
inline void
joinAndFilterMatches(TLocalHolder & lH)
{
//     auto const originalNum = length(lH.matches);
    lH.statusStr << "Sorting hits…";
    if (lH.options.isTerminal)
        myPrint(lH.options, 2, lH.statusStr);

    /// sort

    double start = sysTime();
//    std::sort(begin(lH.matches, Standard()), end(lH.matches, Standard()));
   std::sort(lH.matches.begin(), lH.matches.end());
//     lH.matches.sort();
    double finish = sysTime() - start;
    lH.statusStr << " done. " ;
    if (lH.options.verbosity > 2)
        lH.statusStr << finish << "[s]. ";
    myPrint(lH.options, 2, lH.statusStr);

//     // join and remove duplicates, filter too short matches
//     std::cout <<": Joining and filtering matches…" << std::flush;
//     start = sysTime();


//     std::cout << "Matches before DEBUG checks:" << lH.matches.size() << "\n";
//     // DEBUG
//     for (auto it = lH.matches.begin(),
//               itEnd = lH.matches.end(),
//               itN = std::next(it, 1);
//          (it != itEnd) && (itN != itEnd);
//         )
//     {
//         itN = std::next(it, 1);
//         auto const & redQryInfix = infix(lH.gH.redQrySeqs[it->qryId],
//                                          it->qryStart,
//                                          it->qryStart+lH.options.seedLength);
//         auto const & redSubjInfix = infix(value(indexText(lH.gH.dbIndex), it->subjId),
//                                           it->subjStart,
//                                           it->subjStart+lH.options.seedLength);
// 
//         /*if (((it->qryEnd - it->qryStart) != lH.options.seedLength) ||
//             ((it->subjEnd - it->subjStart) != lH.options.seedLength))
//         {
//             std::cout << "Match with UNLENGTH\n"
//                       << " qId: " << it->qryId << " sId: " << it->subjId
//                       << "\n " << redQryInfix
//                       << "\n " << redSubjInfix
//                       << "\n\n";
//             lH.matches.erase(it);
//         } else */
//         if (quickHamming(redQryInfix,redSubjInfix) >
//                    lH.options.maxSeedDist)
//         {
//             std::cout << "Match with UNMATCH\n"
//                       << " qId: " << it->qryId << " sId: " << it->subjId
//                       << "\n " << redQryInfix
//                       << "\n " << redSubjInfix
//                       << "\n\n";
//             lH.matches.erase(it);
//         }
// 
//         it = itN;
//     }
// 
//     std::cout << "Matches after DEBUG checks:" << lH.matches.size() << "\n";

    
//     // merge
//     for (auto it = lH.matches.begin(),
//               itEnd = lH.matches.end(),
//               itN = std::next(it, 1);
//          (it != itEnd) && (itN != itEnd);
//         )
//     {
//         if (// same sequence pair
//             (it->qryId == itN->qryId) &&
//             (it->subjId == itN->subjId) &&
//             // and matches overlap
//             (overlap(*it, *itN, lH.options.seedGravity)))
//         {
//             mergeUnto(*it, *itN);
//             lH.matches.erase(itN);
//             itN = std::next(it, 1);
//             ++lH.stats.hitsMerged;
//             continue; // goto to next match
//         }
// 
//         if ((itN->qryEnd - itN->qryStart) < lH.options.minSeedLength)
//         {
//             lH.matches.erase(itN);
//             ++lH.stats.hitsTooShort;
//         }
//         ++it;
//         itN = std::next(it, 1);
// 
//     }
//     // check length of first list element
//     auto it = lH.matches.begin();
//     if ((it->qryEnd - it->qryStart) < lH.options.minSeedLength)
//     {
//         lH.matches.erase(it);
//         ++lH.stats.hitsTooShort;
//     }

//     if (lH.options.verbosity >= 3)
//     {
//         finish = sysTime() - start;
//         THREADLINE
//         std::cout <<":  done.\n" << std::flush;
//         THREADLINE
//         std::cout <<": Runtime: " << finish << "s \n" << std::flush;
//         THREADLINE
//         std::cout <<": No of matches before joining and filtering: " << originalNum
//                 << ". After: " << length(lH.matches) << "\n\n" << std::flush;
//     }

}


template <typename TBlastMatch,
          typename TLocalHolder>
inline int
computeBlastMatch(TBlastMatch   & bm,
                  Match         const & m,
                  TLocalHolder  & lH)
{
    using TFormat = typename TLocalHolder::TGlobalHolder::TFormat;

    const unsigned long qryLength = length(value(lH.gH.qrySeqs, m.qryId));

    auto const &  curQry = value(lH.gH.qrySeqs, m.qryId);
    auto const & curSubj = value(lH.gH.subjSeqs, m.subjId);


    SEQAN_ASSERT_LEQ(bm.qStart, bm.qEnd);
    SEQAN_ASSERT_LEQ(bm.sStart, bm.sEnd);


//     auto qryInfix = infix(curQry,
//                                 bm.qStart,
//                                 bm.qEnd);
//     auto subjInfix  = infix(curSubj,
//                                 bm.sStart,
//                                 bm.sEnd);

//     std::cout << "Query Id: " << m.qryId
//               << "\t TrueQryId: " << getTrueQryId(bm.m, lH.options, TFormat())
//               << "\t length(qryIds): " << length(qryIds)
//               << "Subj Id: " << m.subjId
//               << "\t TrueSubjId: " << getTrueSubjId(bm.m, lH.options, TFormat())
//               << "\t length(subjIds): " << length(subjIds) << "\n\n";

    resize(rows(bm.align), 2);
    auto & row0  = row(bm.align, 0);
    auto & row1  = row(bm.align, 1);

    assignSource(row0, infix(curQry, bm.qStart, bm.qEnd));
    assignSource(row1, infix(curSubj,bm.sStart, bm.sEnd));

//     std::cout << "== Positions\n";
//     std::cout << "   " <<  bm.qStart << " - " << bm.qEnd << " [before ali]\n";
//     std::cout << bm.align << std::endl;

    int scr = 0;

//     unsigned short seedLeng = 0;
//     double         seedE = 0;
//     double         seedB = 0;

    unsigned row0len = bm.qEnd - bm.qStart;
    unsigned row1len = bm.sEnd - bm.sStart;
    unsigned short band = (!lH.options.hammingOnly) * (lH.options.maxSeedDist);

//     // TODO FIGURE THIS OUT
//     if ((row0len > (lH.options.seedLength + band)) ||
//         (row1len > (lH.options.seedLength + band)))
//     {
//         #pragma omp atomic
//         ++mergeCount;
//         std::cout << "qrId " << m.qryId << "\tsId: " << m.subjId << "\trow0len: " << row0len << "\trow1len: " << row1len << "\n";
//         std::cout << source(row0) << "\n";
//         std::cout << source(row1) << "\n";
//     }

    auto seedsInSeed = std::max(row0len, row1len) / lH.options.seedLength;

    unsigned short maxDist =  0;
    if (lH.options.maxSeedDist <= 1)
        maxDist = std::abs(int(row1len) - int(row0len));
    else
        maxDist = std::abs(int(row1len) - int(row0len)) + (seedsInSeed * band);

    // fast local alignment without DP-stuff
    if (maxDist == 0)
    {
        int scores[row0len]; // non-iso, waiting for std::dynarray
        scores[0] = 0;
        unsigned newEnd = 0;
        unsigned newBeg = 0;
        // score the diagonal
        for (unsigned i = 0; i < row0len; ++i)
        {
            scores[i] += score(lH.gH.scoreScheme, source(row0)[i], source(row1)[i]);
            if (scores[i] < 0)
            {
                scores[i] = 0;
            } else if (scores[i] >= scr)
            {
                scr = scores[i];
                newEnd = i + 1;
            }
            if (i <row0len -1)
                scores[i+1] = scores[i];
        }
        if (newEnd == 0)
        {
            return PREALIGNSCORE;
        }
        // backtrack
        for (unsigned i = newEnd - 1; i > 0; --i)
        {
            if (scores[i] == 0)
            {
                newBeg = i + 1;
                break;
            }
        }
        setEndPosition(row0, newEnd);
        setEndPosition(row1, newEnd);
        setBeginPosition(row0, newBeg);
        setBeginPosition(row1, newBeg);

    } else
    {
        // compute with DP-code
//         scr = localAlignment(bm.align, lH.gH.scoreScheme, -maxDist, +maxDist);
        scr = localAlignment2(bm.align, lH.gH.scoreScheme, -maxDist, +maxDist, lH.alignContext);
    }

    // save new bounds of alignment
    bm.qEnd    =  bm.qStart  + endPosition(row0);
    bm.qStart  += beginPosition(row0);

    bm.sEnd   =  bm.sStart + endPosition(row1);
    bm.sStart += beginPosition(row1);


    if (scr < lH.options.minSeedScore)
        return PREALIGNSCORE;

    // OLD WAY extension with birte's code
    if (false)
    {
    //     std::cout << "   " <<  bm.qStart << " - " << bm.qEnd << " [after ali]\n";
    //     std::cout << bm.align << std::endl;
        decltype(lH.gH.scoreScheme) extScheme(lH.gH.scoreScheme);
        setScoreGapOpen  (extScheme, -8);
        setScoreGapExtend(extScheme, -8);

        Seed<Simple>           seed(bm.sStart, bm.qStart,
                                    bm.sEnd, bm.qEnd);
        extendSeed(seed,
                curSubj,
                curQry,
                EXTEND_BOTH,
                extScheme,
//                    lH.gH.scoreScheme,
                int(lH.options.xDropOff),
                GappedXDrop());

        bm.sStart = beginPositionH(seed);
        bm.qStart = beginPositionV(seed);
        bm.sEnd   = endPositionH(seed);
        bm.qEnd   = endPositionV(seed);

        assignSource(row0, infix(curQry,
                                bm.qStart,
                                bm.qEnd));
        assignSource(row1, infix(curSubj,
                                bm.sStart,
                                bm.sEnd));

        //DEBUG
        auto oldscr = scr;

        scr = localAlignment(bm.align,
                            lH.gH.scoreScheme,
//                                 alignConfig,
                            lowerDiagonal(seed)-beginDiagonal(seed),
                            upperDiagonal(seed)-beginDiagonal(seed));
        // save new bounds of alignment
        bm.qEnd    =  bm.qStart  + endPosition(row0);
        bm.qStart  += beginPosition(row0);

        bm.sEnd   =  bm.sStart + endPosition(row1);
        bm.sStart += beginPosition(row1);

        if (scr < 0) // alignment got screwed up
        {
            std::cout << "SCREW UP\n";
            std::cout << "beginDiag: " << beginDiagonal(seed)
                << "\tlowDiag: " << lowerDiagonal(seed)
                << "\tupDiag: " << upperDiagonal(seed) << '\n';
            std::cout << "oldscore: " << oldscr
                    << "\tseedscore: " << score(seed)
                    << "\tscore: " << scr << '\n';
            std::cout << bm.align << '\n';
        }


    }
    //NEW WAY extension with dp
    if (true)
    {
        if (false) // ungapped second prealign
        {
            Tuple<decltype(bm.qStart), 4> positions =
                { { bm.qStart, bm.sStart, bm.qEnd, bm.sEnd} };

            decltype(lH.gH.scoreScheme) extScheme(lH.gH.scoreScheme);
            setScoreGapOpen  (extScheme, -100);
            setScoreGapExtend(extScheme, -100);
            scr = extendAlignment(bm.align,
                                  lH.alignContext,
                                    scr,
                                    curQry,
                                    curSubj,
                                    positions,
                                    EXTEND_BOTH,
                                    0, // band of 0 size
                                    0, // band of 0 size
                                    1, // xdrop of 1
                                    extScheme);
            bm.qStart  = beginPosition(row0);
            bm.qEnd    = endPosition(row0);

            bm.sStart =  beginPosition(row1);
            bm.sEnd   =  endPosition(row1);
        }

        if (((bm.qStart > 0) && (bm.sStart > 0)) ||
            ((bm.qEnd < qryLength - 1) && (bm.sEnd < length(curSubj) -1)))
        {
            // we want to allow more gaps in longer query sequences
            switch (lH.options.band)
            {
                case -3: maxDist = ceil(log2(qryLength)); break;
                case -2: maxDist = floor(sqrt(qryLength)); break;
                case -1: break;
                default: maxDist = lH.options.band; break;
            }

            Tuple<decltype(bm.qStart), 4> positions =
                    { { bm.qStart, bm.sStart, bm.qEnd, bm.sEnd} };

            if (lH.options.band != -1)
            {
                if (lH.options.xDropOff != -1)
                {
                    scr = extendAlignment(bm.align,
                                          lH.alignContext,
                                        scr,
                                        curQry,
                                        curSubj,
                                        positions,
                                        EXTEND_BOTH,
                                        -maxDist,
                                        +maxDist,
                                        lH.options.xDropOff,
                                        lH.gH.scoreScheme);
                } else
                {
                    //TODO add alignContext to other calls
                    scr = extendAlignment(bm.align,
                                        scr,
                                        curQry,
                                        curSubj,
                                        positions,
                                        EXTEND_BOTH,
                                        -maxDist,
                                        +maxDist,
                                        lH.gH.scoreScheme);
                }
            } else
            {
                if (lH.options.xDropOff != -1)
                {
                    scr = extendAlignment(bm.align,
                                        scr,
                                        curQry,
                                        curSubj,
                                        positions,
                                        EXTEND_BOTH,
                                        lH.options.xDropOff,
                                        lH.gH.scoreScheme);
                } else
                {
                    scr = extendAlignment(bm.align,
                                        scr,
                                        curQry,
                                        curSubj,
                                        positions,
                                        EXTEND_BOTH,
                                        lH.gH.scoreScheme);
                }
            }
            bm.sStart = beginPosition(row1);
            bm.qStart = beginPosition(row0);
            bm.sEnd   = endPosition(row1);
            bm.qEnd   = endPosition(row0);

    //         std::cout << "AFTER:\n" << bm.align << "\n";
        }
    }

//     std::cerr << "AFTEREXT:\n "<< bm.align << "\n";

    if (scr <= 0)
    {
        std::cout << "## LATE FAIL\n" << bm.align << '\n';

        return ALIGNSCORE;
    }
//     std::cout << "##LINE: " << __LINE__ << '\n';

//     std::cout << "ALIGN BEFORE STATS:\n" << bm.align << "\n";

    calcStatsAndScore(bm, lH.gH.scoreScheme);
//     const unsigned long qryLength = length(row0);
    bm.bitScore = calcBitScore(bm.score, lH.gH.blastScoringAdapter);
    // TODO possibly cache the lengthAdjustments
    auto const lengthAdj = _lengthAdjustment(lH.gH.dbSpecs.dbTotalLength,
                                             qryLength,
                                             lH.gH.blastScoringAdapter);
    bm.eValue = calcEValue(bm.score,
                         lH.gH.dbSpecs.dbTotalLength - lengthAdj,
                         qryLength - lengthAdj,
                         lH.gH.blastScoringAdapter);

    if (bm.eValue > lH.options.eCutOff)
    {

        return ALIGNEVAL;
    }

    // make a Match with updated positions
    // blast is 1-indexed, not 0-indexed, and last pos is lsat pos
    // ON the sequence. Seqan end positions are 1 behind this evens out for end
//     m.qryStart  = bm.qStart;
//     m.qryEnd    = bm.qEnd;
//     m.subjStart = bm.sStart;
//     m.subjEnd   = bm.sEnd;

    // UNTRANSLATE and add 1
//     bm.qStart  = getTrueQryStartPos (m.qryId, bm.qStart, bm.qEnd,
//                                      lH.options, TFormat());
//     bm.qEnd    = getTrueQryEndPos   (m.qryId, bm.qStart, bm.qEnd,
//                                      lH.options, TFormat());
//     bm.sStart  = getTrueSubjStartPos(m.subjId, bm.sStart, bm.sEnd,
//                                      lH.options, TFormat());
//     bm.sEnd    = getTrueSubjEndPos  (m.subjId, bm.sStart, bm.sEnd,
//                                      lH.options, TFormat());

    bm.qFrameShift = getQryFrameShift(m.qryId, lH.options, TFormat()) + 1;
    if (qryIsReverseComplemented(m.qryId, lH.options, TFormat()))
        bm.qFrameShift = -bm.qFrameShift;

    bm.sFrameShift = getSubjFrameShift(m.subjId, lH.options, TFormat()) + 1;
    if (qryIsReverseComplemented(m.subjId, lH.options, TFormat()))
        bm.sFrameShift = -bm.sFrameShift;

//     std::cout << "Successfull Hit hat origSeedLeng: " << seedLeng
//               << " and origE: " << seedE << "\n";
//     (void)seedLeng;
//     (void)seedE;
//     (void)seedB;
//     std::cout << "See Bit Score: " << seedB << "\n" << std::flush;
    return 0;
}

template <typename TStream, typename TLocalHolder>
inline int
iterateMatches(TStream & stream, TLocalHolder & lH)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TFormat       = typename TGlobalHolder::TFormat;
    using TPos          = typename Match::TPos;
    using TBlastRecord  = BlastRecord<
                           typename Value<typename TGlobalHolder::TIds>::Type,// const &,
                           typename Value<typename TGlobalHolder::TIds>::Type,// const &,
                           TPos,
                           typename TLocalHolder::TAlign>;

    constexpr TPos TPosMax = std::numeric_limits<TPos>::max();
//     constexpr uint8_t qFactor = qHasRevComp(TFormat()) ? 3 : 1;
//     constexpr uint8_t sFactor = sHasRevComp(TFormat()) ? 3 : 1;

    double start = sysTime();
    lH.statusStr << "Extending and writing hits…";
    myPrint(lH.options, 2, lH.statusStr);

    //DEBUG
//     std::cout << "Length of matches:   " << length(lH.matches);
//     for (auto const & m :  lH.matches)
//     {
//         std::cout << m.qryId << "\t" << getTrueQryId(m,lH.options, TFormat()) << "\n";
//     }

    // outer loop over records (all matches of one query)
    for (auto it = lH.matches.begin(),
              itN = std::next(it, 1),
              itEnd = lH.matches.end();
         it != itEnd;
         ++it)
    {
        itN = std::next(it,1);
        auto const trueQryId = getTrueQryId(it->qryId,lH.options, TFormat());

        TBlastRecord record(lH.gH.qryIds[trueQryId]);

        record.qLength = qHasFrames(TFormat())
                            ? lH.gH.untransQrySeqLengths[trueQryId]
                            : length(lH.gH.qrySeqs[it->qryId]);

        // inner loop over matches per record
        for (; it != itEnd; ++it)
        {
            auto const trueSubjId = getTrueSubjId(it->subjId,
                                                  lH.options,
                                                  TFormat());
            itN = std::next(it,1);
//             std::cout << "FOO\n" << std::flush;
//             std::cout << "QryStart: " << it->qryStart << "\n" << std::flush;
//             std::cout << "SubjStart: " << it->subjStart << "\n" << std::flush;
//             std::cout << "BAR\n" << std::flush;
            if (!((it->qryStart == TPosMax) && (it->subjStart == TPosMax)))
            {
//                 std::cout << "BAX\n" << std::flush;
                // create blastmatch in list without copy or move
                record.matches.emplace_back(
                lH.gH.qryIds [getTrueQryId(it->qryId, lH.options, TFormat())],
                lH.gH.subjIds[getTrueSubjId(it->subjId, lH.options, TFormat())]);

                auto & bm = back(record.matches);

                bm.qStart    = it->qryStart;
                bm.qEnd      = it->qryStart + lH.options.seedLength;
                bm.sStart    = it->subjStart;
                bm.sEnd      = it->subjStart + lH.options.seedLength;

                // merge putative siblings into this match
                for (auto it2 = itN;
                     (it2 != itEnd) &&
                     (trueQryId == getTrueQryId(it2->qryId,
                                                lH.options,
                                                TFormat())) &&
                     (trueSubjId == getTrueSubjId(it2->subjId,
                                                  lH.options,
                                                  TFormat()));
                    ++it2)
                {
                    // same frame
                    if (qryIsSameFrame(it->qryId, it2->qryId,
                                        lH.options, TFormat()) &&
                        subjIsSameFrame(it->subjId, it2->subjId,
                                        lH.options, TFormat()))
                    {

//                         TPos const qDist = (it2->qryStart >= bm.qEnd)
//                                             ? it2->qryStart - bm.qEnd // upstream
//                                             : 0; // overlap
// 
//                         TPos sDist = TPosMax; // subj match region downstream of *it
//                         if (it2->subjStart >= bm.sEnd) // upstream
//                             sDist = it2->subjStart - bm.sEnd;
//                         else if (it2->subjStart >= it->subjStart) // overlap
//                             sDist = 0;

                        // due to sorting it2->qryStart never <= it->qStart
                        // so subject sequences must have same order
                        if (it2->subjStart < it->subjStart)
                            continue;

                        long const qDist = it2->qryStart - bm.qEnd;
                        long const sDist = it2->subjStart - bm.sEnd;

                        if ((qDist == sDist) &&
                            (qDist <= (long)lH.options.seedGravity))
                        {
                            bm.qEnd = std::max(bm.qEnd,
                                               (TPos)(it2->qryStart
                                               + lH.options.seedLength));
                            bm.sEnd = std::max(bm.sEnd,
                                               (TPos)(it2->subjStart
                                               + lH.options.seedLength));
                            ++lH.stats.hitsMerged;

                            it2->qryStart = TPosMax;
                            it2->subjStart = TPosMax;
                        }
                    }
                }

                // do the extension and statistics
                int lret = computeBlastMatch(bm, *it, lH);

                switch (lret)
                {
                    case COMPUTERESULT_::SUCCESS:
                        // set remaining unset members of match
                        bm.sLength = sHasFrames(TFormat())
                                    ? lH.gH.untransSubjSeqLengths[trueSubjId]
                                    : length(lH.gH.subjSeqs[it->subjId]);
                        bm.qLength = record.qLength;
    //                     lastMatch = ma;
    //                     ++lH.stats.goodMatches;
                        break;
                    case ALIGNEVAL:
                        ++lH.stats.hitsFailedExtendAlignEValTest;
                        break;
                    case ALIGNSCORE:
                        ++lH.stats.hitsFailedExtendAlignScoreTest;
                        break;
                    case PREALIGNEVAL:
                        ++lH.stats.hitsFailedSeedAlignEValTest;
                        break;
                    case PREALIGNSCORE:
                        ++lH.stats.hitsFailedSeedAlignScoreTest;
                        break;
                    default:
                        return lret;
                        break;
                }

                if (lret != 0)// discard match
                {
                    record.matches.pop_back();
                } else
                {
                    // filter the following matches for duplicate-candidates
                    for (auto it2 = itN;
                         (it2 != itEnd) &&
                         (trueQryId == getTrueQryId(it2->qryId,
                                                    lH.options,
                                                    TFormat())) &&
                         (trueSubjId == getTrueSubjId(it2->subjId,
                                                      lH.options,
                                                      TFormat()));
                         ++it2)
                    {
                        // same frame and same range
                        if ((it->qryId == it2->qryId) &&
                            (it->subjId == it2->subjId) &&
                            (intervalOverlap(it2->qryStart,
                                             it2->qryStart + lH.options.seedLength,
                                             bm.qStart,
                                             bm.qEnd) > 0) &&
                            (intervalOverlap(it2->subjStart,
                                             it2->subjStart + lH.options.seedLength,
                                             bm.sStart,
                                             bm.sEnd) > 0))
                        {
                            // deactivated alignment check to get rid of
                            // duplicates early on
//                             auto const & row0 = row(bm.align, 0);
//                             auto const & row1 = row(bm.align, 1);
//                             // part of alignment
//                             if (toSourcePosition(row0,
//                                                  toViewPosition(row1,
//                                                                 it2->subjStart
//                                                                 - bm.sStart))
//                                 == TPos(it2->qryStart - bm.qStart))
//                             {
                                ++lH.stats.hitsPutativeDuplicate;
                                it2->qryStart = TPosMax;
                                it2->subjStart = TPosMax;
//                             }
                        }
                    }
                }
            }

            // last item or new TrueQryId
            if ((itN == itEnd) ||
                (trueQryId != getTrueQryId(itN->qryId, lH.options, TFormat())))
                break;
        }

        if (length(record.matches) > 0)
        {
            ++lH.stats.qrysWithHit;
            // sort and remove duplicates -> STL, yeah!
            auto const before = record.matches.size();
            record.matches.sort();
            record.matches.unique();
            lH.stats.hitsFinal += record.matches.size();
            lH.stats.hitsDuplicate += before - record.matches.size();

            int lret = 0;
            SEQAN_OMP_PRAGMA(critical(filewrite))
            {
                lret = writeRecord(stream, record, lH.gH.dbSpecs, TFormat());
            }
            if (lret)
                return lret;
        }

    }

    double finish = sysTime() - start;

    lH.statusStr << " done. " << finish << "[s]. ";
    myPrint(lH.options, 2, lH.statusStr);
    return 0;
}

void printStats(StatsHolder const & stats, LambdaOptions const & options)
{
    if (options.verbosity >= 2)
    {
        if (options.isTerminal)
            for (unsigned char i=0; i< options.threads+3; ++i)
                std::cout << std::endl;
        unsigned long rem = stats.hitsAfterSeeding;
        auto const w = _numberOfDigits(rem); // number of digits
        #define R  " " << std::setw(w)
        #define RR " = " << std::setw(w)
        #define BLANKS for (unsigned i = 0; i< w; ++i) std::cout << " ";
        std::cout << "\033[1m   HITS                         "; BLANKS; std::cout << "Remaining\033[0m"
                 << "\n   after Seeding               "; BLANKS; std::cout << R << rem;
        std::cout<< "\n - masked                   " << R << stats.hitsMasked       << RR << (rem-= stats.hitsMasked);
        std::cout<< "\n - merged                   " << R << stats.hitsMerged       << RR << (rem-= stats.hitsMerged);
//         std::cout<< "\n - tooShort                 " << R << stats.hitsTooShort     << RR << (rem-= stats.hitsTooShort);
        std::cout<< "\n - putative duplicates      " << R << stats.hitsPutativeDuplicate   << RR << (rem-= stats.hitsPutativeDuplicate);
        std::cout<< "\n - failed seed align test   " << R << stats.hitsFailedSeedAlignScoreTest  << RR << (rem-= stats.hitsFailedSeedAlignScoreTest);
        std::cout<< "\n - failed extend align test " << R << stats.hitsFailedExtendAlignEValTest << RR << (rem-= stats.hitsFailedExtendAlignEValTest);
        std::cout<< "\n - duplicates               " << R << stats.hitsDuplicate    << "\033[1m" << RR << (rem-= stats.hitsDuplicate)
                 << "\033[0m\n\n";

        if (rem != stats.hitsFinal)
            std::cout << "WARNING: hits dont add up\n";
    }

    if (options.verbosity >= 1)
    {
        auto const w = _numberOfDigits(stats.hitsFinal);
        std::cout << "Number of valid hits:                           "
                  << std::setw(w) << stats.hitsFinal
                  << "\nNumber of Queries with at least one valid hit:  "
                  << std::setw(w) << stats.qrysWithHit
                  << "\n\n";
    }

}

#endif // HEADER GUARD
