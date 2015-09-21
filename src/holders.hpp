// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2015, Hannes Hauswedell, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
// ==========================================================================
// holders.hpp: Data container structs
// ==========================================================================


#ifndef SEQAN_LAMBDA_HOLDERS_H_
#define SEQAN_LAMBDA_HOLDERS_H_


#include "match.hpp"
#include "options.hpp"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct StatsHolder
// ----------------------------------------------------------------------------

struct StatsHolder
{
// seeding
    uint64_t hitsAfterSeeding;
    uint64_t hitsMerged;
    uint64_t hitsTooShort;
    uint64_t hitsMasked;

// pre-extension
    uint64_t hitsFailedPreExtendTest;
    uint64_t hitsPutativeDuplicate;
    uint64_t hitsPutativeAbundant;

// post-extension
    uint64_t hitsFailedExtendPercentIdentTest;
    uint64_t hitsFailedExtendEValueTest;
    uint64_t hitsAbundant;
    uint64_t hitsDuplicate;

// final
    uint64_t hitsFinal;
    uint64_t qrysWithHit;

    StatsHolder()
    {
        clear();
    }

    void clear()
    {
        hitsAfterSeeding = 0;
        hitsMerged = 0;
        hitsTooShort = 0;
        hitsMasked = 0;

        hitsFailedPreExtendTest = 0;
        hitsPutativeDuplicate = 0;
        hitsPutativeAbundant = 0;

        hitsFailedExtendPercentIdentTest = 0;
        hitsFailedExtendEValueTest = 0;
        hitsAbundant = 0;
        hitsDuplicate = 0;

        hitsFinal = 0;
        qrysWithHit = 0;
    }

    StatsHolder plus(StatsHolder const & rhs)
    {
        hitsAfterSeeding += rhs.hitsAfterSeeding;
        hitsMerged += rhs.hitsMerged;
        hitsTooShort += rhs.hitsTooShort;
        hitsMasked += rhs.hitsMasked;

        hitsFailedPreExtendTest += rhs.hitsFailedPreExtendTest;
        hitsPutativeDuplicate += rhs.hitsPutativeDuplicate;
        hitsPutativeAbundant += rhs.hitsPutativeAbundant;

        hitsFailedExtendPercentIdentTest += rhs.hitsFailedExtendPercentIdentTest;
        hitsFailedExtendEValueTest += rhs.hitsFailedExtendEValueTest;
        hitsAbundant += rhs.hitsAbundant;
        hitsDuplicate += rhs.hitsDuplicate;

        hitsFinal += rhs.hitsFinal;
        qrysWithHit += rhs.qrysWithHit;
        return *this;
    }

    StatsHolder operator+(StatsHolder const& rhs)
    {
        StatsHolder tmp(*this);
        return tmp.plus(rhs);
    }

    StatsHolder operator+=(StatsHolder const& rhs)
    {
        this->plus(rhs);
        return *this;
    }
};


// ----------------------------------------------------------------------------
// struct GlobalDataHolder  -- one object per program
// ----------------------------------------------------------------------------

template <typename TRedAlph_,
          typename TScoreScheme_,
          typename TIndexSpec_,
          typename TFileFormat,
          BlastProgram p,
          BlastTabularSpec h>
class GlobalDataHolder
{
public:
    static constexpr BlastProgram blastProgram = p;

    /* Sequence storage types */
    using TStringTag    = Alloc<>;
#if defined(LAMBDA_MMAPPED_DB)
    using TDirectStringTag = MMap<>;
#else
    using TDirectStringTag = TStringTag;
#endif
    using TQryTag  = TStringTag;//typename std::conditional<qNumFrames(p) == 1, TStringTag, TDirectStringTag>::type;
    using TSubjTag = TDirectStringTag; // even if subjects were translated they are now loaded from disk

    /* Possibly translated but yet unreduced sequences */
    template <typename TSpec>
    using TTransSeqs     = StringSet<String<TransAlph<p>, TSpec>, Owner<ConcatDirect<>>>;
    using TTransQrySeqs  = TTransSeqs<TQryTag>;
    using TTransSubjSeqs = TTransSeqs<TSubjTag>;

    TTransQrySeqs       qrySeqs;
    TTransSubjSeqs      subjSeqs;

    /* Reduced sequence objects, either as modstrings or as references to trans-strings */
    using TRedAlph       = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    template <typename TSpec>
    using TRedAlphModString = ModifiedString<String<TransAlph<p>, TSpec>,
                                ModView<FunctorConvert<TransAlph<p>, TRedAlph>>>;

    static constexpr bool alphReduction = !std::is_same<TransAlph<p>, TRedAlph>::value;

    using TRedQrySeqs   = typename std::conditional<
                            alphReduction,
                            StringSet<TRedAlphModString<TQryTag>, Owner<ConcatDirect<>>>, // modview
                            TTransQrySeqs &>::type;                                       // reference to owner
    using TRedSubjSeqs  = typename std::conditional<
                            alphReduction,
                            StringSet<TRedAlphModString<TSubjTag>, Owner<ConcatDirect<>>>, // modview
                            TTransSubjSeqs &>::type;                                       // reference to owner

    TRedQrySeqs         redQrySeqs;
    TRedSubjSeqs        redSubjSeqs;

    /* INDECES AND THEIR TYPE */
    static constexpr bool
    indexIsFM           = std::is_same<TIndexSpec_, TFMIndex<>>::value;
    using TIndexSpec    = TIndexSpec_;
    using TDbIndex      = Index<typename std::remove_reference<TRedSubjSeqs>::type, TIndexSpec>;

    TDbIndex            dbIndex;

    // TODO maybe remove these for other specs?
    using TPositions    = typename StringSetLimits<TTransQrySeqs>::Type;
    TPositions          untransQrySeqLengths; // used iff qIsTranslated(p)
    TPositions          untransSubjSeqLengths; // used iff sIsTranslated(p)

    using TMasking      = StringSet<String<unsigned>, Owner<ConcatDirect<>>>;
    TMasking            segIntStarts;
    TMasking            segIntEnds;

    template <typename TSpec>
    using TIds          = StringSet<String<char, TSpec>, Owner<ConcatDirect<>>>;
    using TQryIds       = TIds<TQryTag>;
    using TSubjIds      = TIds<TSubjTag>;
    TQryIds             qryIds;
    TSubjIds            subjIds;

    // OUTPUT FILE //
    using TScoreScheme  = TScoreScheme_;
    using TIOContext    = BlastIOContext<TScoreScheme, p, h>;
    using TFile         = typename std::conditional<std::is_same<TFileFormat, BlastTabular>::value,
                                                    BlastTabularFileOut<TIOContext>,
                                                    BlastReportFileOut<TIOContext>>::type;
    TFile               outfile;

    StatsHolder                 stats;

    GlobalDataHolder() :
        redQrySeqs(qrySeqs), redSubjSeqs(subjSeqs), stats()
    {}
};

// ----------------------------------------------------------------------------
// struct LocalDataHolder  -- one object per thread
// ----------------------------------------------------------------------------

template <typename TMatch,
          typename TGlobalHolder_,
          typename TScoreExtension>
class LocalDataHolder
{
public:
    using TGlobalHolder = TGlobalHolder_;
    using TRedQrySeq    = typename Value<typename std::remove_reference<typename TGlobalHolder::TRedQrySeqs>::type>::Type;
    using TSeeds        = StringSet<typename Infix<TRedQrySeq const>::Type>;
    using TSeedIndex    = Index<TSeeds, IndexSa<>>;


    // references to global stuff
    LambdaOptions     const & options;
    TGlobalHolder /*const*/ & gH;
    static constexpr BlastProgram blastProgram = TGlobalHolder::blastProgram;

    // this is the localHolder for the i-th part of the queries
    uint64_t            i;

    // regarding range of queries
    uint64_t            indexBeginQry;
    uint64_t            indexEndQry;

    // regarding seedingp
    TSeeds              seeds;
    TSeedIndex          seedIndex;
//     std::forward_list<TMatch>   matches;
    std::vector<TMatch>   matches;
    std::vector<typename Match::TQId>   seedRefs;  // mapping seed -> query
    std::vector<uint16_t>               seedRanks; // mapping seed -> relative rank

    // regarding extension
    using TAlignRow0 = Gaps<typename Infix<typename Value<typename TGlobalHolder::TTransQrySeqs>::Type>::Type,
                            ArrayGaps>;
    using TAlignRow1 = Gaps<typename Infix<typename Value<typename TGlobalHolder::TTransSubjSeqs>::Type>::Type,
                            ArrayGaps>;
    using TDPContext = DPContext<typename Value<typename TGlobalHolder::TScoreScheme>::Type, TScoreExtension>;
    using TAliExtContext = AliExtContext_<TAlignRow0, TAlignRow1, TDPContext>;

    TAliExtContext      alignContext;


    // regarding the gathering of stats
    StatsHolder         stats;

    // progress string
    std::stringstream   statusStr;

    // constructor
    LocalDataHolder(LambdaOptions     const & _options,
                    TGlobalHolder     /*const*/ & _globalHolder) :
        options(_options), gH(_globalHolder), stats()
    {}

    // copy constructor SHALLOW COPY ONLY, REQUIRED FOR firsprivate()
    LocalDataHolder(LocalDataHolder const & rhs) :
        options(rhs.options), gH(rhs.gH), stats()
    {}

    void init(uint64_t const _i)
    {
        i = _i;

        if (options.doubleIndexing)
        {
            indexBeginQry = (length(gH.qrySeqs) / options.queryPart) * i;
            indexEndQry = (i+1 == options.queryPart) // last interval
                            ? length(gH.qrySeqs) // reach until end
                            : (length(gH.qrySeqs) / options.queryPart) * (i+1);
            // make sure different frames of one sequence in same interval
            indexBeginQry -= (indexBeginQry % qNumFrames(blastProgram));
            indexEndQry -= (indexEndQry % qNumFrames(blastProgram));
        } else
        {
            indexBeginQry = qNumFrames(blastProgram) * i;
            indexEndQry = qNumFrames(blastProgram) * (i+1);
        }

        clear(seeds);
        clear(seedIndex);
        matches.clear();
        seedRefs.clear();
        seedRanks.clear();
//         stats.clear();
        statusStr.clear();
        statusStr.precision(2);
    }

};

// perform a fast local alignment score calculation on the seed and see if we
// reach above threshold
// WARNING the following function only works for hammingdistanced seeds
template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension>
inline bool
seedLooksPromising(
            LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> const & lH,
            TMatch const & m)
{
    // no pre-scoring, but still filter out XXX and NNN hits
//     if (!lH.options.preScoring))
//     {
//         for (unsigned i = m.qryStart, count = 0; i < m.qryStart + lH.options.seedLength; ++i)
//             if (lH.gH.qrySeqs[m.qryId][i] == unkownValue<typename TGlobalHolder::TRedAlph>())
//                 if (++count > lH.options.maxSeedDist)
//                     return false;
//         return true;
//     }

    int64_t effectiveQBegin = m.qryStart;
    int64_t effectiveSBegin = m.subjStart;
    uint64_t effectiveLength = lH.options.seedLength * lH.options.preScoring;
    if (lH.options.preScoring > 1)
    {
        effectiveQBegin -= (lH.options.preScoring - 1) *
                           lH.options.seedLength / 2;
        effectiveSBegin -= (lH.options.preScoring - 1) *
                           lH.options.seedLength / 2;
//         std::cout << effectiveQBegin << "\t" << effectiveSBegin << "\n";
        int64_t min = std::min(effectiveQBegin, effectiveSBegin);
        if (min < 0)
        {
            effectiveQBegin -= min;
            effectiveSBegin -= min;
            effectiveLength += min;
        }

        effectiveLength = std::min({
                            length(lH.gH.qrySeqs[m.qryId]) - effectiveQBegin,
                            length(lH.gH.subjSeqs[m.subjId]) - effectiveSBegin,
                            effectiveLength});
//         std::cout << effectiveQBegin << "\t" << effectiveSBegin << "\t"
//                   << effectiveLength << "\n";
    }

    auto const & qSeq = infix(lH.gH.qrySeqs[m.qryId],
                              effectiveQBegin,
                              effectiveQBegin + effectiveLength);
    auto const & sSeq = infix(lH.gH.subjSeqs[m.subjId],
                              effectiveSBegin,
                              effectiveSBegin + effectiveLength);
    int maxScore = 0;

    int scores[effectiveLength+1]; // C99, C++14, -Wno-vla before that
    scores[0] = 0;

    // score the diagonal
    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        scores[i] += score(seqanScheme(context(lH.gH.outfile).scoringScheme), qSeq[i], sSeq[i]);
        if (scores[i] < 0)
            scores[i] = 0;
        else if (scores[i] > maxScore)
            maxScore = scores[i];

//         if (i < static_cast<uint64_t>(effectiveLength - 1)) // TODO remove
        scores[i+1] = scores[i];
    }

    return (maxScore >= int(lH.options.preScoringThresh * effectiveLength));
}

template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension,
          typename TSeedId,
          typename TSubjOcc>
inline void
onFindImpl(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
           TSeedId const & seedId,
           TSubjOcc subjOcc)
{
    if (TGlobalHolder::indexIsFM) // positions are reversed
        setSeqOffset(subjOcc,
                     length(lH.gH.subjSeqs[getSeqNo(subjOcc)])
                     - getSeqOffset(subjOcc)
                     - lH.options.seedLength);

    Match m {static_cast<Match::TQId>(lH.seedRefs[seedId]),
             static_cast<Match::TSId>(getSeqNo(subjOcc)),
             static_cast<Match::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset),
             static_cast<Match::TPos>(getSeqOffset(subjOcc))};

    bool discarded = false;
    auto const halfSubjL = lH.options.seedLength /  2;

    if (!sIsTranslated(lH.gH.blastProgram))
    {
        for (unsigned k = 0; k < length(lH.gH.segIntStarts[m.subjId]); ++k)
        {
            // more than half of the seed falls into masked interval
            if (intervalOverlap(m.subjStart,
                                m.subjStart + lH.options.seedLength,
                                lH.gH.segIntStarts[m.subjId][k],
                                lH.gH.segIntEnds[m.subjId][k])
                    >= halfSubjL)
            {
                ++lH.stats.hitsMasked;
                discarded = true;
                break;
            }
        }
    }

     if ((!discarded) && (!seedLooksPromising(lH, m)))
     {
         discarded = true;
         ++lH.stats.hitsFailedPreExtendTest;
     }

    if (!discarded)
        lH.matches.emplace_back(m);
}

template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension,
          typename TFinder>
inline void
onFindDoubleIndex(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
                  TFinder const & finder)
{
    auto qryOccs = getOccurrences(back(finder.patternStack));
    auto subjOccs = getOccurrences(back(finder.textStack));

    lH.stats.hitsAfterSeeding += length(qryOccs) * length(subjOccs);

    for (unsigned i = 0; i < length(qryOccs); ++i)
        for (unsigned j = 0; j < length(subjOccs); ++j)
            onFindImpl(lH, getSeqNo(qryOccs[i]), subjOccs[j]);
}

template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension,
          typename TSeedListIterator,
          typename TIndexIterator>
inline void
onFindSingleIndex(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
                  TSeedListIterator const & seedIt,
                  TIndexIterator const & indexIt)
{
    auto qryOcc = position(seedIt);
    auto subjOccs = getOccurrences(indexIt);

    lH.stats.hitsAfterSeeding += length(subjOccs);

    for (unsigned j = 0; j < length(subjOccs); ++j)
        onFindImpl(lH, qryOcc, subjOccs[j]);
}

#endif // SEQAN_LAMBDA_HOLDERS_H_
