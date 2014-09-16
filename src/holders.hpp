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
    unsigned long hitsAfterSeeding;
    unsigned long hitsMerged;
    unsigned long hitsTooShort;
    unsigned long hitsMasked;

// extension
    unsigned long hitsFailedSeedAlignScoreTest;
    unsigned long hitsFailedSeedAlignEValTest;
    unsigned long hitsFailedExtendAlignScoreTest;
    unsigned long hitsFailedExtendAlignEValTest;
    unsigned long hitsPutativeDuplicate;
    unsigned long hitsDuplicate;

// final
    unsigned long hitsFinal;
    unsigned long qrysWithHit;


//     void hitsFinal()
//     {
//         return hitsAfterSeeding - hitsMerged - 

    void clear()
    {
        hitsFinal       = 0;
        hitsAfterSeeding      = 0;
        hitsMerged     = 0;
        hitsTooShort   = 0;
        hitsMasked      = 0;

        qrysWithHit    = 0;
        hitsPutativeDuplicate = 0;
        hitsDuplicate     = 0;

        hitsFailedExtendAlignEValTest   = 0;
        hitsFailedExtendAlignScoreTest  = 0;
        hitsFailedSeedAlignEValTest = 0;
        hitsFailedSeedAlignScoreTest = 0;
    }

    StatsHolder plus(StatsHolder const & rhs)
    {
        hitsFinal       += rhs.hitsFinal;
        hitsAfterSeeding       +=  rhs.hitsAfterSeeding      ;
        hitsMerged      +=  rhs.hitsMerged     ;
        hitsTooShort    +=  rhs.hitsTooShort   ;
        hitsMasked       +=  rhs.hitsMasked      ;

        qrysWithHit     +=  rhs.qrysWithHit    ;
        hitsPutativeDuplicate += rhs.hitsPutativeDuplicate;
        hitsDuplicate      +=  rhs.hitsDuplicate     ;

        hitsFailedExtendAlignEValTest    +=  rhs.hitsFailedExtendAlignEValTest   ;
        hitsFailedExtendAlignScoreTest   +=  rhs.hitsFailedExtendAlignScoreTest  ;
        hitsFailedSeedAlignEValTest +=  rhs.hitsFailedSeedAlignEValTest;
        hitsFailedSeedAlignScoreTest+=  rhs.hitsFailedSeedAlignScoreTest;
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


template <typename TSet1,
          typename TSet2,
          MyEnableIf<std::is_same<TSet1, TSet2>::value> = 0>
inline TSet1 &
initHelper(TSet1 && /**/, TSet2 & set2)
{
    return set2;
}

template <typename TSet1,
          typename TSet2,
          MyEnableIf<!std::is_same<TSet1, TSet2>::value> = 0>
inline TSet1 &&
initHelper(TSet1 && set1, TSet2 & /**/)
{
    return std::move(set1);
}


template <typename TRedAlph_,
          typename TScoreScheme,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
class GlobalDataHolder
{
public:
    using TRedAlph      = RedAlph<p, TRedAlph_>;
    using TFormat       = BlastFormat<m,p,g>;
    using TTransSeqs    = TCDStringSet<TransAlph<p>>;
    using TRedSeqs      = TCDStringSet<TRedAlph>;
    using TRedSeqsACT   = typename std::conditional<
                            std::is_same<TTransSeqs, TRedSeqs>::value,
                            TRedSeqs &,
                            TRedSeqs>::type;
    using TSubjSeqs     = typename std::conditional<
                            std::is_same<TTransSeqs, TRedSeqs>::value,
                            TTransSeqs &,
                            TTransSeqs>::type;
    using TDbSAIndex    = Index<TRedSeqs, IndexSa<> >;
    using TDbFMIndex    = Index<TRedSeqs, FMIndex<> >;
    using TPositions    = typename StringSetLimits<TTransSeqs>::Type;
    using TIds          = StringSet<CharString, Owner<ConcatDirect<>>>;
    using TMasking      = StringSet<String<unsigned>, Owner<ConcatDirect<>>>;
    using TBlastScoringAdapter = BlastScoringAdapter<TScoreScheme>;

    TDbSAIndex                  dbSAIndex;
    TDbFMIndex                  dbFMIndex;
    bool                        dbIndexIsFM = false;
    BlastDbSpecs<>              dbSpecs;

    // these are always translated, except for BLASTN mode, but *never reduced*
    TTransSeqs                  qrySeqs;
    TSubjSeqs                   subjSeqs = initHelper(TTransSeqs(),
                                                      indexText(dbSAIndex));
    // reduced query sequences if using reduction, otherwise const & = qrySeqs
    TRedSeqsACT                 redQrySeqs = initHelper(TRedSeqs(), qrySeqs);


    // TODO maybe remove these for other specs?
    TPositions                  untransQrySeqLengths; // used iff qHasFrames(p)
    TPositions                  untransSubjSeqLengths; // used iff sHasFrames(p)

    TMasking                    segIntStarts;
    TMasking                    segIntEnds;

    TIds                        qryIds;
    TIds                        subjIds;

    TScoreScheme                scoreScheme;
    TBlastScoringAdapter        blastScoringAdapter;

    StatsHolder                 stats;
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
    using TRedQrySeq    = String<typename TGlobalHolder::TRedAlph, TCDSpec>;
    using TSeeds        = StringSet<typename Infix<TRedQrySeq const>::Type>;
    using TSeedIndex    = Index<TSeeds,IndexSa<>>;

    // references to global stuff
    LambdaOptions     const & options;
    TGlobalHolder /*const*/ & gH;

    // this is the localHolder for the i-th part of the queries
    unsigned short      i;

    // regarding range of queries
    unsigned long       indexBeginQry;
    unsigned long       indexEndQry;

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

    LocalDataHolder(LambdaOptions     const & _options,
                    TGlobalHolder     /*const*/ & _globalHolder) :
        options(_options), gH(_globalHolder)
    {}

    void init(const int _i)
    {
        i = _i;
//         indexBeginQry = length(gH.qrySeqs) / options.queryPart * i;
//         indexEndQry = length(gH.qrySeqs) / options.queryPart * (i+1) + 1;
        indexBeginQry = (length(gH.qrySeqs) / options.queryPart) * i;
        indexEndQry = (i+1 == options.queryPart) // last interval
                        ? length(gH.qrySeqs) // reach until end
                        : (length(gH.qrySeqs) / options.queryPart) * (i+1);

//         std::cout << "Thread " << i << " beg: " << indexBeginQry
//                   << " end: " << indexEndQry << "\n" << std::flush;
        clear(seedIndex);
        matches.clear();
        clear(seedRefs);
        clear(seedRanks);
        stats.clear();
        statusStr.clear();
        statusStr.precision(2);
    }

};

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
    if (lH.gH.dbIndexIsFM)
        setSeqOffset(subjOcc,
                     suffixLength(subjOcc, lH.gH.dbFMIndex)
                     - lH.options.seedLength);

    Match m {static_cast<Match::TQId>(lH.seedRefs[seedId]),
             static_cast<Match::TSId>(getSeqNo(subjOcc)),
             static_cast<Match::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset),
             static_cast<Match::TPos>(getSeqOffset(subjOcc))};

    bool masked = false;
    auto const halfSubjL = lH.options.seedLength /  2;

    for (unsigned k = 0;
            k < length(lH.gH.segIntStarts[m.subjId]);
            ++k)
    {
        // more than half of the seed falls into masked interval
        if (intervalOverlap(m.subjStart,
                            m.subjStart + lH.options.seedLength,
                            lH.gH.segIntStarts[m.subjId][k],
                            lH.gH.segIntEnds[m.subjId][k])
                >= halfSubjL)
        {
            ++lH.stats.hitsMasked;
            masked = true;
            break;
        }
    }

    if (!masked)
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

    for (unsigned i = 0; i< length(qryOccs); ++i)
        for (unsigned j = 0; j< length(subjOccs); ++j)
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

    for (unsigned j = 0; j< length(subjOccs); ++j)
        onFindImpl(lH, qryOcc, subjOccs[j]);
}

#endif // SEQAN_LAMBDA_HOLDERS_H_
