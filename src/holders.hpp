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



// base
template <typename TAlph_,
          typename TRedAlph_,
          typename TScoreScheme,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
class GlobalDataHolderBase
{
public:
    using TAlph         = TAlph_;
    using TRedAlph      = TRedAlph_;
//     using TScoreScheme  = TScoreScheme;
    using TFormat       = BlastFormat<m,p,g>;
    using TUnredSeqs    = TCDStringSet<TAlph>;
    using TRedSeqs      = TCDStringSet<TRedAlph>;
    using TDbIndex      = Index<TRedSeqs,IndexSa<> >;
    using TIds          = StringSet<CharString, Owner<ConcatDirect<>>>;
    using TMasking      = StringSet<String<unsigned>, Owner<ConcatDirect<>>>;
    using TBlastParams  = typename BlastStatisticalParameters<TScoreScheme>::Type;

    TUnredSeqs                  qrySeqs;
    TUnredSeqs                  subjSeqs;
    TDbIndex                    dbIndex;
    unsigned long long          dbTotalLength;
    unsigned long               dbNumberOfSeqs;

    TMasking                    segIntStarts;
    TMasking                    segIntEnds;

    TIds                        qryIds;
    TIds                        subjIds;

    TScoreScheme                scoreScheme;
    TBlastParams                blastParams;

    StatsHolder                 stats;
};


// protein
template <typename TRedAlph_,
          typename TScoreScheme,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
class GlobalDataHolder :
    public GlobalDataHolderBase<AminoAcid, TRedAlph_, TScoreScheme, m, p, g>
{
public:
    TCDStringSet<TRedAlph_>     redQrySeqs;
};

// unreduced protein
template <typename TScoreScheme,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
class GlobalDataHolder<AminoAcid, TScoreScheme, m, p, g> :
    public GlobalDataHolderBase<AminoAcid, AminoAcid, TScoreScheme, m, p, g>
{
public:
    using Base = GlobalDataHolderBase<AminoAcid, AminoAcid, TScoreScheme, m, p, g>;
    TCDStringSet<AminoAcid> const & redQrySeqs = Base::qrySeqs;
};

// blastN
template <typename TScoreScheme,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
class GlobalDataHolder<Dna5, TScoreScheme, m, p, g> :
    public GlobalDataHolderBase<Dna5, Dna5, TScoreScheme, m, p, g>
{
public:
    using Base = GlobalDataHolderBase<Dna5, Dna5, TScoreScheme, m, p, g>;
    TCDStringSet<Dna5> const & redQrySeqs = Base::qrySeqs;
};


// TODO add aliExtContext

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

    // regarding seeding
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
                typename TGlobalHolder::TUnredSeqs>::Type>::Type,
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
          typename TFinder>
inline void
onFind(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
       TFinder const & finder)
{
//     unsigned maxLen = lH.options.seedLength +
//                             ((!lH.options.hammingOnly) * lH.options.maxSeedDist);
//     TODO fix this in finder:
//     if ((length(finder.patternStack) - 1 > maxLen) ||
//         (length(finder.textStack) - 1 > maxLen))
//     {
// //         std::cout << __LINE__ << "\n";
//         return;
//     }

    auto qryOccs = getOccurrences(back(finder.patternStack));
    auto subjOccs = getOccurrences(back(finder.textStack));

//     std::cout << __LINE__ << "\n";
//     std::cout << length(qryOccs) << "\t" << length(subjOccs) << "\t"
//               << length(finder.patternStack) << std::endl;
//     for (auto const & s : finder.scoreStack)
//     {
//         for (auto const & ss : s)
//             std::cout << unsigned(ss) << ' ';
//         std::cout << std::endl;
//     }

    auto const seedNum = length(qryOccs) * length(subjOccs);
//     reserve(lH.matches, length(lH.matches)+seedNum, Generous());
    lH.stats.hitsAfterSeeding += seedNum;
    for (unsigned i = 0; i< length(qryOccs); ++i)
    {
        for (unsigned j = 0; j< length(subjOccs); ++j)
        {
            auto const & qryOcc  = qryOccs[i];
            auto const & subjOcc = subjOccs[j];

            auto const seedId = getSeqNo(qryOcc);
            Match m {static_cast<Match::TQId>(lH.seedRefs[seedId]),
                     static_cast<Match::TSId>(getSeqNo(subjOcc)),
                     static_cast<Match::TPos>(lH.seedRanks[seedId] * lH.options.seedLength),
                     static_cast<Match::TPos>(getSeqOffset(subjOcc))};

            bool masked = false;

            for (unsigned k = 0;
                 k < length(lH.gH.segIntStarts[m.subjId]);
                 ++k)
            {
                auto const halfSubjL = lH.options.seedLength /  2;
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
    }
}



#endif // SEQAN_LAMBDA_HOLDERS_H_
