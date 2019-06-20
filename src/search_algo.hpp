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
// lambda.hpp: contains the main progam pipeline
// ==========================================================================

#pragma once

#include <type_traits>
#include <iomanip>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/misc/terminal.h>

#include <seqan/index.h>

#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include <seqan/align_extend.h>

// ============================================================================
// Forwards
// ============================================================================

// template <typename TRedAlph_,
//           typename TIndexSpec_,
//           typename TFileFormat,
//           seqan::BlastProgram p,
//           BlastTabularSpec h>
// class GlobalDataHolder;
//
// template <typename TGlobalHolder_,
//           typename TScoreExtension>
// class LocalDataHolder;

// ============================================================================
// Classes, structs, enums
// ============================================================================

enum COMPUTERESULT_
{
    SUCCESS = 0,
    PREEXTEND,
    PERCENTIDENT,
    EVALUE,
    OTHER_FAIL
};

//TODO replace with lambda
// comparison operator to sort SA-Values based on the strings in the SA they refer to
template <typename TSav, typename TStringSet>
struct Comp :
    public ::std::binary_function < TSav, TSav, bool >
{
    TStringSet const & stringSet;

    Comp(TStringSet const & _stringSet)
        : stringSet(_stringSet)
    {}

    inline bool operator() (TSav const & i, TSav const & j) const
    {
         return (value(stringSet,getSeqNo(i)) < value(stringSet,getSeqNo(j)));
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function readIndexOption()
// --------------------------------------------------------------------------

// inline void
// readIndexOption(std::string            & optionString,
//                 std::string      const & optionIdentifier,
//                 LambdaOptions    const & options)
// {
//     std::ifstream f{(options.indexDir + "/option:" + optionIdentifier).c_str(),
//                     std::ios_base::in | std::ios_base::binary};
//     if (f.is_open())
//     {
//         auto fit = directionIterator(f, Input());
//         readLine(optionString, fit);
//         f.close();
//     }
//     else
//     {
//         throw IndexException("Expected option specifier:\n" + options.indexDir + "/option:" + optionIdentifier);
//     }
// }

// --------------------------------------------------------------------------
// Function readIndexOptions()
// --------------------------------------------------------------------------

void readIndexOptions(LambdaOptions & options)
{
    // Check that directory exists and is readable
//     struct stat path_stat;
// //     stat(indexFilePath.c_str(), &path_stat);
//     if (stat(indexFilePath.c_str(), &path_stat) || !S_ISDIR(path_stat.st_mode))
//         throw IndexException("Index directory does not exist or is not readable.\n");

    fake_index_file f{options.indexFileOptions};

    if (options.indexFilePath.extension() == ".lba")
    {
        std::ifstream is{options.indexFilePath.c_str(), std::ios::binary};
        cereal::BinaryInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", f));
    } else if (options.indexFilePath.extension() == ".lta")
    {
        std::ifstream is{options.indexFilePath.c_str(), std::ios::binary};
        cereal::JSONInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", f));
    } else
    {
        throw 59;
    }
}

// --------------------------------------------------------------------------
// Function checkRAM()
// --------------------------------------------------------------------------

void
checkRAM(LambdaOptions const & options)
{
    myPrint(options, 1, "Checking memory requirements... ");
    uint64_t ram = getTotalSystemMemory();
    uint64_t sizeIndex = 0;
    uint64_t sizeQuery = 0;

    sizeIndex = fileSize(options.indexFilePath.c_str());

    sizeQuery = fileSize(options.queryFile.c_str());

    uint64_t requiredRAM = ((sizeIndex + sizeQuery) * 11) / 10; // give it +10% TODO verify

    if (requiredRAM >= ram)
    {
        myPrint(options, 1, "done.\n");
        std::cerr << "WARNING: You need approximately " << requiredRAM / 1024 / 1024 << "MB of memory, "
                  << "but you have only " << ram / 1024 / 1024
                  << " :'(\nYou should abort this run and try on a machine with more memory!";
    }

    myPrint(options, 1, "met.\n");
    myPrint(options, 2, "Detected: ", ram / 1024 / 1024, "MB, Estimated: ", requiredRAM / 1024 / 1024, "MB\n\n");
}

// --------------------------------------------------------------------------
// Function prepareScoring()
// --------------------------------------------------------------------------

template <typename TGlobalHolder>
void
prepareScoring(TGlobalHolder       & globalHolder,
               LambdaOptions const & options)
{
    //TODO Seqan3 scoring scheme

//     using TGlobalHolder = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;
    if constexpr (std::is_same_v<typename TGlobalHolder::TScoreScheme, seqan::Score<int, seqan::Simple>>)
    {
        seqan::setScoreMatch(context(globalHolder.outfile).scoringScheme, options.match);
        seqan::setScoreMismatch(context(globalHolder.outfile).scoringScheme, options.misMatch);
    } else
    {
        switch (options.scoringMethod)
        {
            case 45:
                seqan::setScoreMatrixById(seqan::context(globalHolder.outfile).scoringScheme._internalScheme,
                                        seqan::AminoAcidScoreMatrixID::BLOSUM45);
                break;
            case 62:
                seqan::setScoreMatrixById(seqan::context(globalHolder.outfile).scoringScheme._internalScheme,
                                        seqan::AminoAcidScoreMatrixID::BLOSUM62);
                break;
            case 80:
                seqan::setScoreMatrixById(seqan::context(globalHolder.outfile).scoringScheme._internalScheme,
                                        seqan::AminoAcidScoreMatrixID::BLOSUM80);
                break;
            default:
                break;
        }
    }

    setScoreGapOpenBlast(context(globalHolder.outfile).scoringScheme, options.gapOpen);
    setScoreGapExtend(context(globalHolder.outfile).scoringScheme, options.gapExtend);

    if (!isValid(context(globalHolder.outfile).scoringScheme))
        throw std::runtime_error{"Could not computer Karlin-Altschul-Values for Scoring Scheme.\n"};
}

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <typename TGlobalHolder>
void
loadDbIndexFromDisk(TGlobalHolder       & globalHolder,
                    LambdaOptions const & options)
{
    std::string strIdent = "Loading Database Index...";
    myPrint(options, 1, strIdent);
    double start = sysTime();

    if (options.indexFilePath.extension() == ".lba")
    {
        std::ifstream is{options.indexFilePath.c_str(), std::ios::binary};
        cereal::BinaryInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", globalHolder.indexFile));
    } else if (options.indexFilePath.extension() == ".lta")
    {
        std::ifstream is{options.indexFilePath.c_str(), std::ios::binary};
        cereal::JSONInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", globalHolder.indexFile));
    } else
    {
        throw 88;
    }

    if constexpr (!std::is_lvalue_reference_v<typename TGlobalHolder::TRedSbjSeqs>) // is view
    {
        // update view, now that underlying range is available
        globalHolder.redSubjSeqs = typename TGlobalHolder::TRedSbjSeqs{globalHolder.indexFile.transSeqs};
    }

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    // this is actually part of prepareScoring(), but the values are just available now
    if constexpr (sIsTranslated(TGlobalHolder::blastProgram))
    {
        // last value has sum of lengths
        seqan::context(globalHolder.outfileBlastTab).dbTotalLength  = back(globalHolder.untransSubjSeqLengths);
        seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = length(globalHolder.untransSubjSeqLengths) - 1;

        seqan::context(globalHolder.outfileBlastRep).dbTotalLength  = back(globalHolder.untransSubjSeqLengths);
        seqan::context(globalHolder.outfileBlastRep).dbNumberOfSeqs = length(globalHolder.untransSubjSeqLengths) - 1;
    } else
    {
        seqan::context(globalHolder.outfileBlastTab).dbTotalLength  = length(concat(globalHolder.subjSeqs));
        seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = length(globalHolder.subjSeqs);

        seqan::context(globalHolder.outfileBlastRep).dbTotalLength  = length(concat(globalHolder.subjSeqs));
        seqan::context(globalHolder.outfileBlastRep).dbNumberOfSeqs = length(globalHolder.subjSeqs);
    }
}
#if 0
// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------

// BLASTX, TBLASTX
template <typename TSourceAlph, typename TSpec1,
          typename TTargetAlph, typename TSpec2,
          typename TUntransLengths,
          MyEnableIf<!std::is_same<TSourceAlph, TTargetAlph>::value> = 0>
inline void
loadQueryImplTrans(TCDStringSet<seqan::String<TTargetAlph, TSpec1>> & target,
                   TCDStringSet<seqan::String<TSourceAlph, TSpec2>> & source,
                   TUntransLengths                            & untransQrySeqLengths,
                   LambdaOptions                        const & options)
{
    myPrint(options, 1, "translating...");
    // translate
    translate(target,
              source,
              SIX_FRAME,
              options.geneticCode);

    // preserve lengths of untranslated sequences
    resize(untransQrySeqLengths,
           length(source.limits),
           Exact());

#ifdef __clang__
    SEQAN_OMP_PRAGMA(parallel for)
#else
    SEQAN_OMP_PRAGMA(parallel for simd)
#endif
    for (uint32_t i = 0; i < (length(untransQrySeqLengths) - 1); ++i)
        untransQrySeqLengths[i] = source.limits[i + 1] - source.limits[i];

    // save sum of lengths (both strings have n + 1 elements
    back(source.limits) = length(source.concat);
}

// BLASTN
template <typename TSpec1,
          typename TSpec2,
          typename TUntransLengths>
inline void
loadQueryImplTrans(TCDStringSet<String<TransAlph<BlastProgram::BLASTN>, TSpec1>> & target,
                   TCDStringSet<String<TransAlph<BlastProgram::BLASTN>, TSpec2>> & source,
                   TUntransLengths                                                & /**/,
                   LambdaOptions                                            const & options)
{
    using TAlph = TransAlph<BlastProgram::BLASTN>;
//     using TReverseCompl =  ModifiedString<ModifiedString<String<TAlph>,
//                             ModView<FunctorComplement<TAlph>>>, ModReverse>;
    myPrint(options, 1, " generating reverse complements...");
    // no need for translation, but we need reverse complements
    resize(target.concat, length(source.concat) * 2);
    resize(target.limits, length(source) * 2 + 1);

    target.limits[0] = 0;
    uint64_t const l = length(target.limits) - 1;
    for (uint64_t i = 1; i < l; i+=2)
    {
        target.limits[i] = target.limits[i-1] + length(source[i/2]);
        target.limits[i+1] = target.limits[i] + length(source[i/2]);
    }

    FunctorComplement<TAlph> functor;
    uint64_t tBeg, tBegNext, len, sBeg;
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic) private(tBeg, tBegNext, len, sBeg))
    for (uint64_t i = 0; i < (l - 1); i+=2)
    {
        tBeg = target.limits[i];
        tBegNext = target.limits[i+1];
        len = tBegNext - tBeg;
        sBeg = source.limits[i/2];
        // avoid senseless copying and iterate manually
        for (uint32_t j = 0; j < len; ++j)
        {
            target.concat[tBeg + j] = source.concat[sBeg + j];
            target.concat[tBegNext + j] = functor(source.concat[sBeg+len-j-1]);
        }
    }
}

// BLASTP, TBLASTN
template <typename TSpec1,
          typename TSpec2,
          typename TUntransLengths>
inline void
loadQueryImplTrans(TCDStringSet<String<TransAlph<BlastProgram::BLASTP>, TSpec1>> & target,
                   TCDStringSet<String<TransAlph<BlastProgram::BLASTP>, TSpec2>> & source,
                   TUntransLengths           & /**/,
                   LambdaOptions       const & /**/)
{
    // no need for translation, but sequences have to be in right place
    swap(target, source);
}

template <BlastTabularSpec h,
          seqan::BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
void
loadQuery(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>      & globalHolder,
          LambdaOptions                                                 & options)
{
    using TGH = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;
    double start = sysTime();

    std::string strIdent = "Loading Query Sequences and Ids...";
    myPrint(options, 1, strIdent);

    TCDStringSet<String<OrigQryAlph<p>, typename TGH::TQryTag>> origSeqs;

    try
    {
        SeqFileIn infile(toCString(options.queryFile));
        myReadRecords(globalHolder.qryIds, origSeqs, infile);
    }
    catch(std::exception const & e)
    {
        throw QueryException{"There was an file system or format error."};
    }

    if (length(origSeqs) == 0)
    {
        throw QueryException{"Zero sequences submitted."};
    }

    // translate
    loadQueryImplTrans(globalHolder.qrySeqs,
                       origSeqs,
                       globalHolder.untransQrySeqLengths,
                       options);

    // sam and bam need original sequences if translation happened
    if (qIsTranslated(TGH::blastProgram) && (options.outFileFormat > 0) &&
        (options.samBamSeq > 0))
        std::swap(origSeqs, globalHolder.untranslatedQrySeqs);

    if (TGH::alphReduction)
        globalHolder.redQrySeqs.limits = globalHolder.qrySeqs.limits;

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    unsigned long maxLen = 0ul;
    for (auto const & s : globalHolder.qrySeqs)
        if (length(s) > maxLen)
            maxLen = length(s);

    myPrint(options, 2, "Runtime: ", finish, "s \n",
            "Number of effective query sequences: ",
            length(globalHolder.qrySeqs), "\nLongest query sequence: ",
            maxLen, "\n\n");

    if (length(globalHolder.qrySeqs) >= std::numeric_limits<typename TGH::TMatch::TQId>::max())
    {
        throw QueryException{"Too many sequences submitted. The maximum (including frames) is " +
                             std::to_string(std::numeric_limits<typename TGH::TMatch::TQId>::max()) +
                             "."};
    }

    if (maxLen >= std::numeric_limits<typename TGH::TMatch::TPos>::max())
    {
        throw QueryException{"One or more of your query sequences are too long. The maximum length is " +
                             std::to_string(std::numeric_limits<typename TGH::TMatch::TPos>::max()) +
                             "."};
    }

    // TODO: after changing this, make options const again
    if (options.extensionMode == LambdaOptions::ExtensionMode::AUTO)
    {
        if (maxLen <= 100)
        {
        #if defined(SEQAN_SIMD_ENABLED)
            options.extensionMode = LambdaOptions::ExtensionMode::FULL_SIMD;
        #else
            options.extensionMode = LambdaOptions::ExtensionMode::FULL_SERIAL;
        #endif
            options.xDropOff = -1;
            options.filterPutativeAbundant = false;
            options.filterPutativeDuplicates = false;
            options.mergePutativeSiblings = false;
        }
        else
        {
            options.extensionMode = LambdaOptions::ExtensionMode::XDROP;
        }
    }
}

/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function generateSeeds()
// --------------------------------------------------------------------------

// perform a fast local alignment score calculation on the seed and see if we
// reach above threshold
// WARNING the following function only works for hammingdistanced seeds
template <typename TGlobalHolder,
          typename TScoreExtension>
inline bool
seedLooksPromising(LocalDataHolder<TGlobalHolder, TScoreExtension> const & lH,
                   typename TGlobalHolder::TMatch const & m)
{
    int64_t effectiveQBegin = m.qryStart;
    int64_t effectiveSBegin = m.subjStart;
    uint64_t actualLength = m.qryEnd - m.qryStart;
    uint64_t effectiveLength = std::max(static_cast<uint64_t>(lH.options.seedLength * lH.options.preScoring),
                                        actualLength);

    if (effectiveLength > actualLength)
    {
        effectiveQBegin -= (effectiveLength - actualLength) / 2;
        effectiveSBegin -= (effectiveLength - actualLength) / 2;

        int64_t min = std::min(effectiveQBegin, effectiveSBegin);
        if (min < 0)
        {
            effectiveQBegin -= min;
            effectiveSBegin -= min;
            effectiveLength += min;
        }

        effectiveLength = std::min({
                            static_cast<uint64_t>(std::ranges::size(lH.gH.qrySeqs[m.qryId]) - effectiveQBegin),
                            static_cast<uint64_t>(std::ranges::size(lH.gH.subjSeqs[m.subjId]) - effectiveSBegin),
                            effectiveLength});
    }

    auto const & qSeq = lH.gH.qrySeqs[m.qryId]
                      | seqan3::view::slice(effectiveQBegin, effectiveQBegin + effectiveLength);
    auto const & sSeq = lH.gH.subjSeqs[m.subjId]
                      | seqan3::view::slice(effectiveSBegin, effectiveSBegin + effectiveLength);

    int             s = 0;
    int      maxScore = 0;
    int const thresh  = lH.options.preScoringThresh * effectiveLength;

    // score the diagonal
    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        s += score(seqanScheme(context(lH.gH.outfile).scoringScheme), qSeq[i], sSeq[i]);
        if (s < 0)
            s = 0;
        else if (s > maxScore)
            maxScore = s;

        if (maxScore >= thresh)
            return true;
    }

    return false;
}

// --------------------------------------------------------------------------
// Function onFind()
// --------------------------------------------------------------------------
#if 0
template <typename TGlobalHolder,
          typename TScoreExtension,
          typename TSeedId,
          typename TSubjOcc>
inline void
onFind(LocalDataHolder<TGlobalHolder, TScoreExtension> & lH,
       TSeedId const & seedId,
       TSubjOcc subjOcc)
{
    using TMatch = typename TGlobalHolder::TMatch;
    SEQAN_ASSERT_LEQ_MSG(getSeqOffset(subjOcc) + lH.options.seedLength,
                         length(lH.gH.subjSeqs[getSeqNo(subjOcc)]),
                         "Seed reaches beyond end of subject sequence! Please report a bug with your files at "
                         "http://www.seqan.de/lambda !");

    if (TGlobalHolder::indexIsFM) // positions are reversed
        setSeqOffset(subjOcc,
                     length(lH.gH.subjSeqs[getSeqNo(subjOcc)])
                     - getSeqOffset(subjOcc)
                     - lH.options.seedLength);

    TMatch m {static_cast<typename TMatch::TQId>(lH.seedRefs[seedId]),
              static_cast<typename TMatch::TSId>(getSeqNo(subjOcc)),
              static_cast<typename TMatch::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset),
              static_cast<typename TMatch::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset + lH.options.seedLength),
              static_cast<typename TMatch::TPos>(getSeqOffset(subjOcc)),
              static_cast<typename TMatch::TPos>(getSeqOffset(subjOcc) + lH.options.seedLength)};

    bool discarded = false;

     if (!seedLooksPromising(lH, m))
     {
         discarded = true;
         ++lH.stats.hitsFailedPreExtendTest;
     }

    if (!discarded)
        lH.matches.emplace_back(m);
}

template <typename TGlobalHolder,
          typename TScoreExtension,
          typename TSubjOcc>
inline void
onFindVariable(LocalDataHolder<TGlobalHolder, TScoreExtension> & lH,
               TSubjOcc subjOcc,
               typename TGlobalHolder::TMatch::TQId const seedId,
               typename TGlobalHolder::TMatch::TPos const seedBegin,
               typename TGlobalHolder::TMatch::TPos const seedLength)
{
    using TMatch = typename TGlobalHolder::TMatch;
    if (TGlobalHolder::indexIsFM) // positions are reversed
        setSeqOffset(subjOcc,
                     length(lH.gH.subjSeqs[getSeqNo(subjOcc)])
                     - getSeqOffset(subjOcc)
                     - seedLength);

    TMatch m {seedId,
              static_cast<typename TGlobalHolder::TMatch::TSId>(getSeqNo(subjOcc)),
              seedBegin,
              static_cast<typename TGlobalHolder::TMatch::TPos>(seedBegin + seedLength),
              static_cast<typename TGlobalHolder::TMatch::TPos>(getSeqOffset(subjOcc)),
              static_cast<typename TGlobalHolder::TMatch::TPos>(getSeqOffset(subjOcc) + seedLength)};

    SEQAN_ASSERT_LT(m.qryStart,  m.qryEnd);
    SEQAN_ASSERT_LT(m.subjStart, m.subjEnd);

    if (!seedLooksPromising(lH, m))
        ++lH.stats.hitsFailedPreExtendTest;
    else
        lH.matches.emplace_back(m);
}

// --------------------------------------------------------------------------
// Function search()
// --------------------------------------------------------------------------

template <typename TIndexIt, typename TGoDownTag, typename TNeedleIt, typename TLambda, typename TLambda2>
inline void
__goDownNoErrors(TIndexIt & indexIt,
                 TGoDownTag const &,
                 TNeedleIt needleIt,
                 TNeedleIt const & needleItEnd,
                 TLambda & continRunnable,
                 TLambda2 & reportRunnable)
{
    TIndexIt prevIndexIt;

    do
        prevIndexIt = indexIt;
    while ((needleIt != needleItEnd) &&
           goDown(indexIt, *(needleIt++), TGoDownTag()) &&
           continRunnable(prevIndexIt, indexIt, true));

    reportRunnable(prevIndexIt, true);
}

//TODO make number of errors configurable
template <typename TIndexIt, typename TGoDownTag, typename TNeedleIt, typename TLambda, typename TLambda2>
inline void
__goDownErrors(TIndexIt const & indexIt,
               TGoDownTag const &,
               TNeedleIt const & needleIt,
               TNeedleIt const & needleItEnd,
               TLambda & continRunnable,
               TLambda2 & reportRunnable)
{
    using TAlph = typename Value<TNeedleIt>::Type;

    unsigned contin = 0;

    if (needleIt != needleItEnd)
    {
        for (unsigned i = 0; i < ValueSize<TAlph>::VALUE; ++i)
        {
            TIndexIt nextIndexIt(indexIt);
            if (goDown(nextIndexIt, static_cast<TAlph>(i), TGoDownTag()) &&
                continRunnable(indexIt, nextIndexIt, ordValue(*needleIt) != i))
            {
                ++contin;
                if (ordValue(*needleIt) == i)
                    __goDownErrors(nextIndexIt, TGoDownTag(), needleIt + 1, needleItEnd, continRunnable, reportRunnable);
                else
                    __goDownNoErrors(nextIndexIt, TGoDownTag(), needleIt + 1, needleItEnd, continRunnable, reportRunnable);
            }
        }
    }

    if (contin == 0)
        reportRunnable(indexIt, false);
}

template <typename TGlobalHolder,
          typename TScoreExtension>
inline void
_searchSingleIndex(LocalDataHolder<TGlobalHolder, TScoreExtension> & lH)
{
    typedef typename Iterator<typename TGlobalHolder::TDbIndex, TopDown<> >::Type TIndexIt;

    // TODO optionize
    size_t constexpr seedHeurFactor = /*TGlobalHolder::indexIsBiFM ? 5 :*/ 10;
    size_t constexpr minResults = 1;

    size_t needlesSum = 0;
    size_t needlesPos = 0;
    size_t oldTotalMatches = 0;

    TIndexIt root(lH.gH.dbIndex);
    TIndexIt indexIt = root;

    std::function<bool(TIndexIt const &, TIndexIt const &, bool const)> continRunnable;

    /* It is important to note some option dependencies:
     * lH.options.maxSeedDist == 0 -> lH.options.seedHalfExact == true
     * lH.options.maxSeedDist == 0 -> TGlobalHolder::indexIsBiFM == false
     * TGlobalHolder::indexIsBiFM  -> lH.options.seedHalfExact == true
     * [these are enforced in options.hpp and save us some comparisons here
     */
    SEQAN_ASSERT((lH.options.maxSeedDist != 0) || lH.options.seedHalfExact);
    SEQAN_ASSERT((lH.options.maxSeedDist != 0) || TGlobalHolder::indexIsBiFM);
    SEQAN_ASSERT((!TGlobalHolder::indexIsBiFM) || lH.options.seedHalfExact);

    size_t const goExactLength = lH.options.seedHalfExact ? (lH.options.seedLength / 2) : 0;

    for (size_t i = lH.indexBeginQry; i < lH.indexEndQry; ++i)
    {
        if (length(lH.gH.redQrySeqs[i]) < lH.options.seedLength)
            continue;

        size_t desiredOccs      = 0;

        // the next sequences belong to a new set of query sequences
        if ((i % qNumFrames(TGlobalHolder::blastProgram)) == 0)
        {
            needlesSum = lH.gH.redQrySeqs.limits[i + qNumFrames(TGlobalHolder::blastProgram)] - lH.gH.redQrySeqs.limits[i];
            // BROKEN:lengthSum(infix(lH.gH.redQrySeqs, lH.indexBeginQry, lH.indexEndQry));
            // the above is faster anyway (but only works on concatdirect sets)

            needlesPos = 0;
            oldTotalMatches = length(lH.matches); // need to subtract matchcount from other queries
        }

        if (lH.options.adaptiveSeeding)
        {
            continRunnable = [&lH, &desiredOccs] (auto const & prevIndexIt, auto const & indexIt, bool const/*hasError*/)
            {
                // ADAPTIVE SEEDING:

                // always continue if minimum seed length not reached
                // TODO currently unclear why considering hasError provides no benefit, and why +1 does
                if (repLength(indexIt) <= (lH.options.seedLength + /*hasError* */ lH.options.seedDeltaIncreasesLength))
                    return true;
                else if (repLength(indexIt) > 2000) // maximum recursion depth
                    return false;

                // always continue if it means not loosing hits
                if (countOccurrences(indexIt) == countOccurrences(prevIndexIt))
                    return true;

                // do vodoo heuristics to see if this hit is to frequent
                if (countOccurrences(indexIt) < desiredOccs)
                    return false;

                return true;
            };
        } else
        {
            continRunnable = [&lH] (auto const &, auto const & indexIt, bool const /*hasError*/)
            {
                // NON-ADAPTIVE
                return (repLength(indexIt) <= (lH.options.seedLength +
                                               /*hasError* */ lH.options.seedDeltaIncreasesLength));
            };
        }

        /* FORWARD SEARCH */
        for (size_t seedBegin = 0; /* below */; seedBegin += lH.options.seedOffset)
        {
            // skip proteine 'X' or Dna 'N'
            while ((lH.gH.qrySeqs[i][seedBegin] == unknownValue<TransAlph<TGlobalHolder::blastProgram>>()) &&
                   (seedBegin <= length(lH.gH.redQrySeqs[i]) - lH.options.seedLength))
                ++seedBegin;

            // termination criterium
            if (seedBegin > length(lH.gH.redQrySeqs[i]) - lH.options.seedLength)
                break;

            indexIt = root;

            if (lH.options.adaptiveSeeding)
            {
                desiredOccs = (length(lH.matches) - oldTotalMatches) >= lH.options.maxMatches
                            ? minResults
                            : (lH.options.maxMatches - (length(lH.matches) - oldTotalMatches)) * seedHeurFactor /
                                std::max((needlesSum - needlesPos - seedBegin) / lH.options.seedOffset, static_cast<size_t>(1));

                if (desiredOccs == 0)
                    desiredOccs = minResults;
            }

            // go down some characters without errors if bidirectional or halfExact
            for (size_t k = 0; k < goExactLength; ++k)
                if (!goDown(indexIt, lH.gH.redQrySeqs[i][seedBegin + k]))
                    break;
            // if unsuccessful, move to next seed
            if (repLength(indexIt) != goExactLength)
                continue;

            auto reportRunnable = [&lH, &i, &seedBegin] (auto const & indexIt, bool const hasError)
            {
                if (repLength(indexIt) >= lH.options.seedLength + hasError * lH.options.seedDeltaIncreasesLength)
                {
                #ifdef LAMBDA_MICRO_STATS
                    appendValue(lH.stats.seedLengths, repLength(indexIt));
                #endif
                    lH.stats.hitsAfterSeeding += countOccurrences(indexIt);
                    for (auto occ : getOccurrences(indexIt))
                        onFindVariable(lH, occ, i, seedBegin, repLength(indexIt));
                }
            };

            if (lH.options.maxSeedDist)
            {
                __goDownErrors(indexIt,
                               Fwd(),
                               begin(lH.gH.redQrySeqs[i], Standard()) + seedBegin + goExactLength,
                               end(lH.gH.redQrySeqs[i], Standard()),
                               continRunnable,
                               reportRunnable);
            }
            else
                __goDownNoErrors(indexIt,
                                 Fwd(),
                                 begin(lH.gH.redQrySeqs[i], Standard()) + seedBegin + goExactLength,
                                 end(lH.gH.redQrySeqs[i], Standard()),
                                 continRunnable,
                                 reportRunnable);
        }

        /* REVERSE SEARCH on BIDIRECTIONAL INDEXES */
        if (TGlobalHolder::indexIsBiFM)
        {
            using   TRevNeedle      = ModifiedString<decltype(lH.gH.redQrySeqs[0]), ModReverse>;
            TRevNeedle revNeedle{lH.gH.redQrySeqs[i]};
            for (size_t seedBegin = lH.options.seedLength - 1; /* below */; seedBegin += lH.options.seedOffset)
            {

                // skip proteine 'X' or Dna 'N'
                while ((lH.gH.qrySeqs[i][seedBegin] == unknownValue<TransAlph<TGlobalHolder::blastProgram>>()) &&
                    seedBegin < length(lH.gH.redQrySeqs[i]))                 // [different abort condition than above]
                    ++seedBegin;

                // termination criterium
                if (seedBegin >= length(lH.gH.redQrySeqs[i]))                // [different abort condition than above]
                    break;

                indexIt = root;

                if (lH.options.adaptiveSeeding)
                {
                    desiredOccs = (length(lH.matches) - oldTotalMatches) >= lH.options.maxMatches
                                ? minResults
                                : (lH.options.maxMatches - (length(lH.matches) - oldTotalMatches)) * seedHeurFactor /
                                    std::max((needlesSum - needlesPos - seedBegin) / lH.options.seedOffset, static_cast<size_t>(1));

                    if (desiredOccs == 0)
                        desiredOccs = minResults;
                }

                // go down seedOffset number of characters without errors
                for (size_t k = 0; k < (lH.options.seedLength - goExactLength); ++k)
                    if (!goDown(indexIt, lH.gH.redQrySeqs[i][seedBegin - k], Rev())) // [rev and  - instead of fwd]
                        break;
                // if unsuccessful, move to next seed
                if (repLength(indexIt) != (lH.options.seedLength - goExactLength))
                    continue;

                auto reportRunnable = [&lH, &i, &seedBegin] (auto const & indexIt, bool const hasOneError)
                {
                    if ((repLength(indexIt) >= lH.options.seedLength) && (hasOneError))        // [must have one error for rev]
                    {
                        //TODO remove debug stuff
                    #ifdef LAMBDA_MICRO_STATS
                        appendValue(lH.stats.seedLengths, repLength(indexIt));
                    #endif
                        lH.stats.hitsAfterSeeding += countOccurrences(indexIt);
                        for (auto occ : getOccurrences(indexIt))                    // [different start pos]
                            onFindVariable(lH, occ, i, seedBegin - repLength(indexIt) + 1,  repLength(indexIt));
                    }
                };

                // [rev and reverse needle]
                __goDownErrors(indexIt,
                               Rev(),
                               end(revNeedle, Standard()) - seedBegin + lH.options.seedLength - goExactLength - 1,
                               end(revNeedle, Standard()),
                               continRunnable,
                               reportRunnable);

            }
        }

        needlesPos += length(lH.gH.redQrySeqs[i]);
    }
}
#endif
template <typename TLocalHolder>
inline void
search(TLocalHolder & lH)
{
//     _searchSingleIndex(lH);

    //TODO new search here
}

// --------------------------------------------------------------------------
// Function joinAndFilterMatches()
// --------------------------------------------------------------------------

template <typename TLocalHolder>
inline void
sortMatches(TLocalHolder & lH)
{
//    std::sort(begin(lH.matches, Standard()), end(lH.matches, Standard()));
//     std::sort(lH.matches.begin(), lH.matches.end());

//     if (lH.matches.size() > lH.options.maxMatches)
//     {
//         MatchSortComp   comp(lH.matches);
//         std::sort(lH.matches.begin(), lH.matches.end(), comp);
//     } else

    if ((lH.options.filterPutativeAbundant) &&
        (lH.matches.size() > lH.options.maxMatches))
        // more expensive sort to get likely targets to front
        myHyperSortSingleIndex(lH.matches, lH.gH);
    else
        std::sort(lH.matches.begin(), lH.matches.end());
}

// --------------------------------------------------------------------------
// Function _setFrames()
// --------------------------------------------------------------------------

template <typename TBlastMatch,
          typename TLocalHolder>
inline void
_setFrames(TBlastMatch                          & bm,
           typename TLocalHolder::TMatch  const & m,
           TLocalHolder                   const &)
{
    if constexpr (seqan::qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = (m.qryId % 3) + 1;
        if (m.qryId % 6 > 2)
            bm.qFrameShift = -bm.qFrameShift;
    } else if constexpr (seqan::qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = 1;
        if (m.qryId % 2)
            bm.qFrameShift = -bm.qFrameShift;
    } else
    {
        bm.qFrameShift = 0;
    }

    if constexpr (seqan::sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = (m.subjId % 3) + 1;
        if (m.subjId % 6 > 2)
            bm.sFrameShift = -bm.sFrameShift;
    } else if constexpr (seqan::sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = 1;
        if (m.subjId % 2)
            bm.sFrameShift = -bm.sFrameShift;
    } else
    {
        bm.sFrameShift = 0;
    }
}

// --------------------------------------------------------------------------
// Function _writeMatches()
// --------------------------------------------------------------------------

template <typename TBlastRecord,
          typename TLocalHolder>
inline void
_writeRecord(TBlastRecord & record,
             TLocalHolder & lH)
{
    if (seqan::length(record.matches) > 0)
    {
        ++lH.stats.qrysWithHit;
        // sort and remove duplicates -> STL, yeah!
        auto const before = record.matches.size();

        if (!lH.options.filterPutativeDuplicates)
        {
            record.matches.sort([] (auto const & m1, auto const & m2)
            {
                return std::tie(m1._n_sId,
                                m1.qStart,
                                m1.qEnd,
                                m1.sStart,
                                m1.sEnd,
                                m1.qFrameShift,
                                m1.sFrameShift) <
                       std::tie(m2._n_sId,
                                m2.qStart,
                                m2.qEnd,
                                m2.sStart,
                                m2.sEnd,
                                m2.qFrameShift,
                                m2.sFrameShift);
            });
            record.matches.unique([] (auto const & m1, auto const & m2)
            {
                return std::tie(m1._n_sId,
                                m1.qStart,
                                m1.qEnd,
                                m1.sStart,
                                m1.sEnd,
                                m1.qFrameShift,
                                m1.sFrameShift) ==
                       std::tie(m2._n_sId,
                                m2.qStart,
                                m2.qEnd,
                                m2.sStart,
                                m2.sEnd,
                                m2.qFrameShift,
                                m2.sFrameShift);
            });
            lH.stats.hitsDuplicate += before - record.matches.size();
        }

        // sort by evalue before writing
        record.matches.sort([] (auto const & m1, auto const & m2)
        {
            return m1.bitScore > m2.bitScore;
        });

        // cutoff abundant
        if (record.matches.size() > lH.options.maxMatches)
        {
            lH.stats.hitsAbundant += record.matches.size() -
                                        lH.options.maxMatches;
            record.matches.resize(lH.options.maxMatches);
        }
        lH.stats.hitsFinal += record.matches.size();

        // compute LCA
        if (lH.options.computeLCA)
        {
            record.lcaTaxId = 0;
            for (auto const & bm : record.matches)
            {
                if ((length(lH.gH.sTaxIds[bm._n_sId]) > 0) && (lH.gH.taxParents[lH.gH.sTaxIds[bm._n_sId][0]] != 0))
                {
                    record.lcaTaxId = lH.gH.sTaxIds[bm._n_sId][0];
                    break;
                }
            }

            if (record.lcaTaxId != 0)
                for (auto const & bm : record.matches)
                    for (uint32_t const sTaxId : lH.gH.sTaxIds[bm._n_sId])
                        if (lH.gH.taxParents[sTaxId] != 0) // TODO do we want to skip unassigned subjects
                            record.lcaTaxId = computeLCA(lH.gH.taxParents, lH.gH.taxHeights, sTaxId, record.lcaTaxId);

            record.lcaId = lH.gH.taxNames[record.lcaTaxId];
        }

        myWriteRecord(lH, record);
    }
}

// --------------------------------------------------------------------------
// Function computeBlastMatch()
// --------------------------------------------------------------------------

template <typename TBlastMatch,
          typename TLocalHolder>
inline void
_setupAlignInfix(TBlastMatch & bm,
                 typename TLocalHolder::TMatch const & m,
                 TLocalHolder & lH)
{
    int64_t startMod = (int64_t)m.subjStart - (int64_t)m.qryStart;

    bm.qEnd = length(lH.gH.qrySeqs[m.qryId]);
    decltype(bm.qEnd) band = _bandSize(bm.qEnd , lH);
    if (startMod >= 0)
    {
        bm.sStart = startMod;
        bm.qStart = 0;
    }
    else
    {
        bm.sStart = 0;
        bm.qStart = -startMod;
    }
    bm.sEnd = _min(bm.sStart + bm.qEnd - bm.qStart + band, length(lH.gH.subjSeqs[m.subjId]));

    if (bm.sStart >= band)
        bm.sStart -= band;
    else
        bm.sStart = 0;

    seqan::assignSource(bm.alignRow0, seqan::infix(lH.gH.qrySeqs[m.qryId],   bm.qStart, bm.qEnd));
    seqan::assignSource(bm.alignRow1, seqan::infix(lH.gH.subjSeqs[m.subjId], bm.sStart, bm.sEnd));
}

template <typename TBlastMatch,
          typename TLocalHolder>
inline auto
_untrueQryId(TBlastMatch const & bm,
             TLocalHolder const &)
{
    using namespace seqan;

    if (qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 6 + bm.qFrameShift - 1;
        else
            return bm._n_qId * 6 - bm.qFrameShift + 2;
    } else if (qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 2;
        else
            return bm._n_qId * 2 + 1;
    } else
    {
        return bm._n_qId;
    }
}

template <typename TBlastMatch,
          typename TLocalHolder>
inline auto
_untrueSubjId(TBlastMatch const & bm,
              TLocalHolder const &)
{
    using namespace seqan;

    if (sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.sFrameShift > 0)
            return bm._n_sId * 6 + bm.sFrameShift - 1;
        else
            return bm._n_sId * 6 - bm.sFrameShift + 2;
    } else if (sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.sFrameShift > 0)
            return bm._n_sId * 2;
        else
            return bm._n_sId * 2 + 1;
    } else
    {
        return bm._n_sId;
    }
}

template <typename TBlastMatch,
          typename TLocalHolder>
inline void
_expandAlign(TBlastMatch & bm,
             TLocalHolder const & lH)
{
    using namespace seqan;

    auto oldQLen = length(source(bm.alignRow0));
    auto oldSLen = length(source(bm.alignRow1));

    // replace source from underneath without triggereng reset
    value(bm.alignRow0._source) = lH.gH.qrySeqs[_untrueQryId(bm, lH)];
    value(bm.alignRow1._source) = lH.gH.subjSeqs[_untrueSubjId(bm, lH)];

    // insert fields into array gaps
    if (bm.alignRow0._array[0] == 0)
        bm.alignRow0._array[1] += bm.qStart;
    else
        insert(bm.alignRow0._array, 0, std::vector<uint64_t>{0, bm.qStart});
    if (bm.alignRow0._array[length(bm.alignRow0._array) - 1] == 0)
        bm.alignRow0._array[length(bm.alignRow0._array) - 2] += length(source(bm.alignRow0)) - oldQLen;
    else
        append(bm.alignRow0._array, std::vector<uint64_t>{length(source(bm.alignRow0)) - oldQLen, 0});

    if (bm.alignRow1._array[0] == 0)
        bm.alignRow1._array[1] += bm.sStart;
    else
        insert(bm.alignRow1._array, 0, std::vector<uint64_t>{0, bm.sStart});
    if (bm.alignRow1._array[length(bm.alignRow1._array) - 1] == 0)
        bm.alignRow1._array[length(bm.alignRow1._array) - 2] += length(source(bm.alignRow1)) - oldSLen;
    else
        append(bm.alignRow1._array, std::vector<uint64_t>{length(source(bm.alignRow1)) - oldSLen, 0});

    // the begin positions from the align object are relative to the infix created above
    bm.qEnd   = bm.qStart + endPosition(bm.alignRow0);
    bm.qStart = bm.qStart + beginPosition(bm.alignRow0);
    bm.sEnd   = bm.sStart + endPosition(bm.alignRow1);
    bm.sStart = bm.sStart + beginPosition(bm.alignRow1);

    // set clipping positions on new gaps objects
    setBeginPosition(bm.alignRow0, bm.qStart);
    setEndPosition(bm.alignRow0, bm.qEnd);
    setBeginPosition(bm.alignRow1, bm.sStart);
    setEndPosition(bm.alignRow1, bm.sEnd);
}

#ifdef SEQAN_SIMD_ENABLED

template <typename TDepSetH,
          typename TDepSetV,
          typename TBlastMatches>
inline void
_setupDepSets(TDepSetH & depSetH, TDepSetV & depSetV, TBlastMatches const & blastMatches)
{
    using namespace seqan;

    using TSimdAlign    = typename SimdVector<int16_t>::Type;
    unsigned constexpr sizeBatch = LENGTH<TSimdAlign>::VALUE;
    unsigned const      fullSize = sizeBatch * ((length(blastMatches) + sizeBatch - 1) / sizeBatch);

    clear(depSetH);
    clear(depSetV);
    reserve(depSetH, fullSize);
    reserve(depSetV, fullSize);

    for (auto const & bm : blastMatches)
    {
        appendValue(depSetH, source(bm.alignRow0));
        appendValue(depSetV, source(bm.alignRow1));
    }

    // fill up last batch
    for (size_t i = length(blastMatches); i < fullSize; ++i)
    {
        appendValue(depSetH, source(back(blastMatches).alignRow0));
        appendValue(depSetV, source(back(blastMatches).alignRow1));
    }
}

template <typename TDepSetH,
          typename TDepSetV,
          typename TBlastMatches,
          typename TLocalHolder,
          bool withTrace>
inline void
_performAlignment(TDepSetH & depSetH,
                  TDepSetV & depSetV,
                  TBlastMatches & blastMatches,
                  TLocalHolder & lH,
                  std::integral_constant<bool, withTrace> const &)
{
    using namespace seqan;

    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TAlignConfig  = AlignConfig2<LocalAlignment_<>,
                                       DPBandConfig<BandOff>,
                                       FreeEndGaps_<True, True, True, True>,
                                       std::conditional_t<withTrace,
                                                          TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> >,
                                                          TracebackOff> >;
    using TSimdAlign    = typename SimdVector<int16_t>::Type;
    using TSimdScore    = Score<TSimdAlign, ScoreSimdWrapper<typename TGlobalHolder::TScoreScheme> >;
    using TSize         = typename Size<typename TLocalHolder::TAlignRow0>::Type;
    using TMatch        = typename TGlobalHolder::TMatch;
    using TPos          = typename TMatch::TPos;
    using TTraceSegment = TraceSegment_<TPos, TSize>;

    unsigned constexpr sizeBatch = LENGTH<TSimdAlign>::VALUE;
    unsigned const      fullSize = sizeBatch * ((length(blastMatches) + sizeBatch - 1) / sizeBatch);

    TSimdScore simdScoringScheme(seqanScheme(context(lH.gH.outfile).scoringScheme));
    StringSet<String<TTraceSegment> > trace;

    // TODO when band is available, create inside block with band
    TAlignConfig config;//(0, 2*band)

    auto matchIt = blastMatches.begin();
    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        auto infSetH = infixWithLength(depSetH, pos, sizeBatch);
        auto infSetV = infixWithLength(depSetV, pos, sizeBatch);

        TSimdAlign resultsBatch;

        clear(trace);
        resize(trace, sizeBatch, Exact());

        // TODO pass in lH.dpSIMDContext
        _prepareAndRunSimdAlignment(resultsBatch,
                                    trace,
                                    infSetH,
                                    infSetV,
                                    simdScoringScheme,
                                    config,
                                    typename TLocalHolder::TScoreExtension());

        for(auto x = pos; x < pos + sizeBatch && x < length(blastMatches); ++x)
        {
            if constexpr (withTrace)
                _adaptTraceSegmentsTo(matchIt->alignRow0, matchIt->alignRow1, trace[x - pos]);
            else
                matchIt->alignStats.alignmentScore = resultsBatch[x - pos];

            ++matchIt;
        }
    }

}

template <typename TLocalHolder>
inline int
iterateMatchesFullSimd(TLocalHolder & lH)
{
    using namespace seqan;

    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
//     using TMatch        = typename TGlobalHolder::TMatch;
//     using TPos          = typename TMatch::TPos;
    using TBlastPos     = uint32_t; //TODO why can't this be == TPos
    using TBlastMatch   = BlastMatch<
                           typename TLocalHolder::TAlignRow,
                           typename TLocalHolder::TAlignRow,
                           TBlastPos,
                           decltype(lH.gH.qryIds[0]),
                           decltype(lH.gH.indexFile.subjIds[0]),
                           >;
    using TBlastRecord  = BlastRecord<TBlastMatch,
                                      decltype(lH.gH.qryIds[0]),
                                      std::vector<std::string>,
                                      decltype(lH.gH.indexFile.taxNames[0]),
                                      uint32_t>;
    // statistics
#ifdef LAMBDA_MICRO_STATS
    ++lH.stats.numQueryWithExt;
    lH.stats.numExtScore += length(lH.matches);

    double start = sysTime();
#endif

    // Prepare string sets with sequences.
    StringSet<typename Source<typename TLocalHolder::TAlignRow0>::Type> depSetH;
    StringSet<typename Source<typename TLocalHolder::TAlignRow1>::Type> depSetV;

    // container of blastMatches (possibly from multiple queries
    decltype(TBlastRecord().matches) blastMatches;

    // create blast matches
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // create blastmatch in list without copy or move
        blastMatches.emplace_back(lH.gH.qryIds [it->qryId / qNumFrames(TGlobalHolder::blastProgram)],
                                    lH.gH.subjIds[it->subjId / sNumFrames(TGlobalHolder::blastProgram)]);

        auto & bm = back(blastMatches);

        bm._n_qId = it->qryId / qNumFrames(TGlobalHolder::blastProgram);
        bm._n_sId = it->subjId / sNumFrames(TGlobalHolder::blastProgram);

        bm.sLength = sIsTranslated(TGlobalHolder::blastProgram)
                        ? lH.gH.untransSubjSeqLengths[bm._n_sId]
                        : length(lH.gH.subjSeqs[it->subjId]);

        _setupAlignInfix(bm, *it, lH);

        _setFrames(bm, *it, lH);

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.sTaxIds[bm._n_sId];
    }
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtend      += sysTime() - start;
    lH.stats.timeExtendTrace += sysTime() - start; //TODO remove this line!

    // filter out duplicates
    start = sysTime();
#endif
    auto before = length(blastMatches);
    blastMatches.sort([] (auto const & l, auto const & r)
    {
        return std::tie(l._n_qId, l._n_sId, l.sStart, l.sEnd, l.qStart, l.qEnd, l.qFrameShift, l.sFrameShift) <
               std::tie(r._n_qId, r._n_sId, r.sStart, r.sEnd, r.qStart, r.qEnd, r.qFrameShift, r.sFrameShift);
    });
    blastMatches.unique([] (auto const & l, auto const & r)
    {
        return std::tie(l._n_qId, l._n_sId, l.sStart, l.sEnd, l.qStart, l.qEnd, l.qFrameShift, l.sFrameShift) ==
               std::tie(r._n_qId, r._n_sId, r.sStart, r.sEnd, r.qStart, r.qEnd, r.qFrameShift, r.sFrameShift);
    });
    lH.stats.hitsDuplicate += (before - length(blastMatches));

    // sort by lengths to minimize padding in SIMD
    blastMatches.sort([] (auto const & l, auto const & r)
    {
        return std::make_tuple(length(source(l.alignRow0)), length(source(l.alignRow1))) <
               std::make_tuple(length(source(r.alignRow0)), length(source(r.alignRow1)));
    });
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeSort += sysTime() - start;

    start = sysTime();
#endif
    // fill batches
    _setupDepSets(depSetH, depSetV, blastMatches);

    // Run extensions WITHOUT ALIGNMENT
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::false_type());

    // copmute evalues and filter based on evalue
    for (auto it = blastMatches.begin(), itEnd = blastMatches.end(); it != itEnd; /*below*/)
    {
        TBlastMatch & bm = *it;

        computeEValueThreadSafe(bm,
                                qIsTranslated(TGlobalHolder::blastProgram)
                                    ? lH.gH.untransQrySeqLengths[bm._n_qId]
                                    : length(lH.gH.qrySeqs[bm._n_qId]),
                                seqan::context(lH.gH.outfile));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            it = blastMatches.erase(it);
            continue;
        }

        ++it;
    }
    if (length(blastMatches) == 0)
        return 0;

    // statistics
#ifdef LAMBDA_MICRO_STATS
    lH.stats.numExtAli += length(blastMatches);
    lH.stats.timeExtend += sysTime() - start;
    start = sysTime();
#endif

    // reset and fill batches
    _setupDepSets(depSetH, depSetV, blastMatches);

    // Run extensions WITH ALIGNMENT
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::true_type());

    // sort by query
    blastMatches.sort([] (auto const & lhs, auto const & rhs)
    {
        return lhs._n_qId < rhs._n_qId;
    });

    // compute the rest of the match properties
    for (auto it = blastMatches.begin(), itEnd = blastMatches.end(); it != itEnd; /*below*/)
    {
        TBlastMatch & bm = *it;

        _expandAlign(bm, lH);

        computeAlignmentStats(bm, seqan::context(lH.gH.outfile));

        if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        {
            ++lH.stats.hitsFailedExtendPercentIdentTest;
            it = blastMatches.erase(it);
            continue;
        }

        computeBitScore(bm, seqan::context(lH.gH.outfile));

        // evalue computed previously

        ++it;
    }
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif

    if (length(blastMatches) == 0)
        return 0;

    // devide matches into records (per query) and write
    for (auto it = blastMatches.begin(), itLast = blastMatches.begin();
         length(blastMatches) > 0;
         /*below*/)
    {
        if ((it == blastMatches.end()) || ((it != blastMatches.begin()) && (it->_n_qId != itLast->_n_qId)))
        {
            // create a record for each query
            TBlastRecord record(lH.gH.qryIds[itLast->_n_qId]);
            record.qLength = (qIsTranslated(TGlobalHolder::blastProgram)
                                ? lH.gH.untransQrySeqLengths[itLast->_n_qId]
                                : length(lH.gH.qrySeqs[itLast->_n_qId]));
            // move the matches into the record
            record.matches.splice(record.matches.begin(),
                                  blastMatches,
                                  blastMatches.begin(),
                                  it);
            // write to file
            _writeRecord(record, lH);

            it = blastMatches.begin();
            itLast = blastMatches.begin();
        } else
        {
            itLast = it;
            ++it;
        }
    }

    return 0;
}

#endif // SEQAN_SIMD_ENABLED

template <typename TLocalHolder>
inline int
iterateMatchesFullSerial(TLocalHolder & lH)
{
    using namespace seqan;

    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
//     using TMatch        = typename TGlobalHolder::TMatch;
//     using TPos          = typename TMatch::TPos;
    using TBlastPos     = uint32_t; //TODO why can't this be == TPos
    using TBlastMatch   = BlastMatch<
                           typename TLocalHolder::TAlignRow0,
                           typename TLocalHolder::TAlignRow1,
                           TBlastPos,
                           typename Value<typename TGlobalHolder::TQryIds>::Type,// const &,
                           typename Value<typename TGlobalHolder::TSubjIds>::Type// const &,
                           >;
    using TBlastRecord  = BlastRecord<TBlastMatch,
                                      typename Value<typename TGlobalHolder::TQryIds>::Type,
                                      std::vector<std::string>,
                                      typename Value<typename TGlobalHolder::TTaxNames>::Type,
                                      uint32_t>;

    auto const trueQryId = lH.matches[0].qryId / qNumFrames(TGlobalHolder::blastProgram);

    TBlastRecord record(lH.gH.qryIds[trueQryId]);
    record.qLength = (qIsTranslated(TGlobalHolder::blastProgram)
                        ? lH.gH.untransQrySeqLengths[trueQryId]
                        : length(lH.gH.qrySeqs[lH.matches[0].qryId]));

    unsigned band = _bandSize(record.qLength, lH);

#ifdef LAMBDA_MICRO_STATS
    double start = sysTime();
#endif

    // create blast matches
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // create blastmatch in list without copy or move
        record.matches.emplace_back(lH.gH.qryIds [it->qryId / qNumFrames(TGlobalHolder::blastProgram)],
                                    lH.gH.subjIds[it->subjId / sNumFrames(TGlobalHolder::blastProgram)]);

        auto & bm = back(record.matches);
        auto &  m = *it;

        bm._n_qId = it->qryId / qNumFrames(TGlobalHolder::blastProgram);
        bm._n_sId = it->subjId / sNumFrames(TGlobalHolder::blastProgram);

        bm.sLength = sIsTranslated(TGlobalHolder::blastProgram)
                        ? lH.gH.untransSubjSeqLengths[bm._n_sId]
                        : length(lH.gH.subjSeqs[it->subjId]);

        _setupAlignInfix(bm, *it, lH);

        _setFrames(bm, m, lH);

        // Run extension WITHOUT TRACEBACK
        typedef AlignConfig2<LocalAlignment_<>,
                            DPBandConfig<BandOn>,
                            FreeEndGaps_<True, True, True, True>,
                            TracebackOff> TAlignConfig;

        DPScoutState_<Default> scoutState;

        bm.alignStats.alignmentScore = _setUpAndRunAlignment(lH.alignContext.dpContext,
                                                             lH.alignContext.traceSegment,
                                                             scoutState,
                                                             source(bm.alignRow0),
                                                             source(bm.alignRow1),
                                                             seqanScheme(context(lH.gH.outfile).scoringScheme),
                                                             TAlignConfig(-band, +band));

        computeEValueThreadSafe(bm, record.qLength, seqan::context(lH.gH.outfile));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            record.matches.pop_back();
            continue;
        }

        // Run extension WITH TRACEBACK
        localAlignment(bm.alignRow0,
                       bm.alignRow1,
                       seqanScheme(context(lH.gH.outfile).scoringScheme),
                       -band,
                       +band);

        _expandAlign(bm, lH);

        computeAlignmentStats(bm, seqan::context(lH.gH.outfile));

        if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        {
            ++lH.stats.hitsFailedExtendPercentIdentTest;
            record.matches.pop_back();
            continue;
        }

        computeBitScore(bm, seqan::context(lH.gH.outfile));

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.sTaxIds[bm._n_sId];
    }

#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif

    _writeRecord(record, lH);

    return 0;
}

template <typename TLocalHolder>
inline int
iterateMatches(TLocalHolder & lH)
{
#ifdef SEQAN_SIMD_ENABLED
    if (lH.options.extensionMode == LambdaOptions::ExtensionMode::FULL_SIMD)
        return iterateMatchesFullSimd(lH);
    else
#endif

    return iterateMatchesFullSerial(lH);

}
#endif
