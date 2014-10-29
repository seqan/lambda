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
        hitsAfterSeeding = rhs.hitsAfterSeeding;
        hitsMerged = rhs.hitsMerged;
        hitsTooShort = rhs.hitsTooShort;
        hitsMasked = rhs.hitsMasked;

        hitsFailedPreExtendTest = rhs.hitsFailedPreExtendTest;
        hitsPutativeDuplicate = rhs.hitsPutativeDuplicate;
        hitsPutativeAbundant = rhs.hitsPutativeAbundant;

        hitsFailedExtendPercentIdentTest = rhs.hitsFailedExtendPercentIdentTest;
        hitsFailedExtendEValueTest = rhs.hitsFailedExtendEValueTest;
        hitsAbundant = rhs.hitsAbundant;
        hitsDuplicate = rhs.hitsDuplicate;

        hitsFinal = rhs.hitsFinal;
        qrysWithHit = rhs.qrysWithHit;
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
          typename TScoreScheme,
          typename TIndexSpec_,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
class GlobalDataHolder
{
public:
    using TFormat       = BlastFormat<m,p,g>;

    // SEQUENCES AND THEIR TYPE //
    using TTransSeqs    = TCDStringSet<TransAlph<p>>;

    using TRedAlph      = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TRedSeqVirt   = ModifiedString<String<TransAlph<p>, PackSpec>,
                                ModView<FunctorConvert<TransAlph<p>,TRedAlph>>>;
    using TRedSeqsVirt  = StringSet<TRedSeqVirt, Owner<ConcatDirect<>>>;

    static bool constexpr
    indexIsFM           = std::is_same<TIndexSpec_, FMIndex<>>::value;
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

    // these are the sequences in *unreduced space*
    TTransSeqs          qrySeqs;
    TTransSeqs          subjSeqs;

    // reduced query sequences if using reduction, otherwise & = transSeqs
    TRedSeqsACT         redQrySeqs = qrySeqs;
    TRedSeqsACT         redSubjSeqs = subjSeqs;

    // INDECES AND THEIR TYPE //
    using TIndexSpec    = TIndexSpec_;
    using TDbIndex      = Index<TRedSeqs, TIndexSpec>;

    TDbIndex            dbIndex;
    BlastDbSpecs<>      dbSpecs;

    // TODO maybe remove these for other specs?
    using TPositions    = typename StringSetLimits<TTransSeqs>::Type;
    TPositions          untransQrySeqLengths; // used iff qHasFrames(p)
    TPositions          untransSubjSeqLengths; // used iff sHasFrames(p)

    using TMasking      = StringSet<String<unsigned>, Owner<ConcatDirect<>>>;
    TMasking            segIntStarts;
    TMasking            segIntEnds;

    using TIds          = StringSet<CharString, Owner<ConcatDirect<>>>;
    TIds                qryIds;
    TIds                subjIds;

    using TBlastScoringAdapter = BlastScoringAdapter<TScoreScheme>;
    TScoreScheme                scoreScheme;
    TBlastScoringAdapter        blastScoringAdapter;

    StatsHolder                 stats = StatsHolder();
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
    using TFormat       = typename TGlobalHolder::TFormat;
    using TRedQrySeq    = typename Value<typename TGlobalHolder::TRedSeqs>::Type;
    using TSeeds        = StringSet<typename Infix<TRedQrySeq const>::Type>;
    using TSeedIndex    = Index<TSeeds,IndexSa<>>;

    // references to global stuff
    LambdaOptions     const & options;
    TGlobalHolder /*const*/ & gH;

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
    typedef Align<
        typename Infix<
            typename Value<
                typename TGlobalHolder::TTransSeqs>::Type>::Type,
                ArrayGaps> TAlign;

    typedef DPContext<typename Value<decltype(gH.scoreScheme)>::Type,
                      TScoreExtension> TDPContext;
    typedef AliExtContext_<TAlign,
                          TDPContext> TAliExtContext;

    TAliExtContext      alignContext;


    // regarding the gathering of stats
    StatsHolder         stats;

    // progress string
    std::stringstream         statusStr;

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
        } else
        {
            indexBeginQry = qNumFrames(TFormat()) * i;
            indexEndQry = qNumFrames(TFormat()) * (i+1);
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
    auto const & qSeq = infix(lH.gH.qrySeqs[m.qryId],
                              m.qryStart,
                              m.qryStart + lH.options.seedLength);
    auto const & sSeq = infix(lH.gH.subjSeqs[m.subjId],
                              m.subjStart,
                              m.subjStart + lH.options.seedLength);
    int maxScore = 0;

    int scores[lH.options.seedLength]; // C99, C++14, -Wno-vla before that
    scores[0] = 0;

    // score the diagonal
    for (unsigned i = 0; i < lH.options.seedLength; ++i)
    {
        scores[i] += score(lH.gH.scoreScheme, qSeq[i], sSeq[i]);
        if (scores[i] < 0)
            scores[i] = 0;
        else if (scores[i] >= maxScore)
            maxScore = scores[i];

        if (i < static_cast<unsigned>(lH.options.seedLength - 1))
            scores[i+1] = scores[i];
    }

    return (maxScore >= lH.options.minSeedScore);
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

     if ((!TGlobalHolder::noReduction) && (!discarded) &&
         (!seedLooksPromising(lH, m)))
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
