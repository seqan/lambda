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
// holders.hpp: Data container structs
// ==========================================================================

#ifndef LAMBDA_SEARCH_DATASTRUCTURES_H_
#define LAMBDA_SEARCH_DATASTRUCTURES_H_

#include <seqan/align_extend.h>

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct Match
// ----------------------------------------------------------------------------

template<typename TAlph>
struct Match
{
    typedef SizeTypeNum_<TAlph>    TQId;
    typedef SizeTypeNum_<TAlph>    TSId;
    typedef SizeTypePos_<TAlph>    TPos;

    TQId qryId;
    TSId subjId;
    TPos qryStart;
    TPos qryEnd;

    TPos subjStart;
    TPos subjEnd;

    inline bool operator== (Match const & m2) const
    {
         return std::tie(qryId, subjId, qryStart, subjStart, qryEnd, subjEnd)
             == std::tie(m2.qryId, m2.subjId, m2.qryStart, m2.subjStart, m2.qryEnd, m2.subjEnd);
    }
    inline bool operator< (Match const & m2) const
    {
         return std::tie(qryId, subjId, qryStart, subjStart, qryEnd, subjEnd)
              < std::tie(m2.qryId, m2.subjId, m2.qryStart, m2.subjStart, m2.qryEnd, m2.subjEnd);
    }
};

template <typename TAlph>
inline void
setToSkip(Match<TAlph> & m)
{
    using TPos          = typename Match<TAlph>::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    m.qryStart = posMax;
    m.subjStart = posMax;
}

template <typename TAlph>
inline bool
isSetToSkip(Match<TAlph> const & m)
{
    using TPos          = typename Match<TAlph>::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    return (m.qryStart == posMax) && (m.subjStart == posMax);
}

template <typename TAlph>
inline void
_printMatch(Match<TAlph> const & m)
{
    std::cout << "MATCH  Query " << m.qryId
              << "(" << m.qryStart << ", " << m.qryEnd
              << ")   on Subject "<< m.subjId
              << "(" << m.subjStart << ", " << m.subjEnd
              << ")" <<  std::endl << std::flush;
}

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
#ifdef LAMBDA_MICRO_STATS
    std::vector<uint16_t> seedLengths;
#endif

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

#ifdef LAMBDA_MICRO_STATS
// times
    double timeGenSeeds;
    double timeSearch;
    double timeSort;
    double timeExtend;
    double timeExtendTrace;

// extension counters
    uint64_t numQueryWithExt;
    uint64_t numExtScore;
    uint64_t numExtAli;
#endif

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

    #ifdef LAMBDA_MICRO_STATS
        seedLengths.clear();
        timeGenSeeds = 0;
        timeSearch = 0;
        timeSort = 0;
        timeExtend = 0;
        timeExtendTrace = 0;

        numQueryWithExt = 0;
        numExtScore = 0;
        numExtAli = 0;
    #endif
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

    #ifdef LAMBDA_MICRO_STATS
        append(seedLengths, rhs.seedLengths);
        timeGenSeeds += rhs.timeGenSeeds;
        timeSearch   += rhs.timeSearch;
        timeSort     += rhs.timeSort;
        timeExtend   += rhs.timeExtend;
        timeExtendTrace += rhs.timeExtendTrace;

        numQueryWithExt += rhs.numQueryWithExt;
        numExtScore += rhs.numExtScore;
        numExtAli += rhs.numExtAli;
    #endif
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

void printStats(StatsHolder const & stats, LambdaOptions const & options)
{
    if (options.verbosity >= 2)
    {
        unsigned long rem = stats.hitsAfterSeeding;
        auto const w = _numberOfDigits(rem); // number of digits
        #define R  " " << std::setw(w)
        #define RR " = " << std::setw(w)
        #define BLANKS for (unsigned i = 0; i< w; ++i) std::cout << " ";
        std::cout << "\033[1m   HITS                         "; BLANKS;
        std::cout << "Remaining\033[0m"
                  << "\n   after Seeding               "; BLANKS;
        std::cout << R << rem;
        if (stats.hitsMasked)
            std::cout << "\n - masked                   " << R << stats.hitsMasked
                      << RR << (rem -= stats.hitsMasked);
        if (options.mergePutativeSiblings)
            std::cout << "\n - merged                   " << R << stats.hitsMerged
                      << RR << (rem -= stats.hitsMerged);
        if (options.filterPutativeDuplicates)
            std::cout << "\n - putative duplicates      " << R
                      << stats.hitsPutativeDuplicate << RR
                      << (rem -= stats.hitsPutativeDuplicate);
        if (options.filterPutativeAbundant)
            std::cout << "\n - putative abundant        " << R
                      << stats.hitsPutativeAbundant   << RR
                      << (rem -= stats.hitsPutativeAbundant);
        if (options.preScoring)
            std::cout << "\n - failed pre-extend test   " << R
                      << stats.hitsFailedPreExtendTest  << RR
                      << (rem -= stats.hitsFailedPreExtendTest);
        std::cout << "\n - failed %-identity test   " << R
                  << stats.hitsFailedExtendPercentIdentTest << RR
                  << (rem -= stats.hitsFailedExtendPercentIdentTest);
        std::cout << "\n - failed e-value test      " << R
                  << stats.hitsFailedExtendEValueTest << RR
                  << (rem -= stats.hitsFailedExtendEValueTest);
        std::cout << "\n - duplicates               " << R
                  << stats.hitsDuplicate              << RR
                  << (rem -= stats.hitsDuplicate);
        std::cout << "\n - abundant                 " << R
                  << stats.hitsAbundant << "\033[1m"  << RR
                  << (rem -= stats.hitsAbundant)
                  << "\033[0m\n\n";

        if (rem != stats.hitsFinal)
            std::cout << "WARNING: hits don't add up\n";

    #ifdef LAMBDA_MICRO_STATS
        std::cout << "Detailed Non-Wall-Clock times:\n"
                  << " genSeeds:    " << stats.timeGenSeeds << "\n"
                  << " search:      " << stats.timeSearch << "\n"
                  << " sort:        " << stats.timeSort << "\n"
                  << " extend:      " << stats.timeExtend << "\n"
                  << " extendTrace: " << stats.timeExtendTrace << "\n\n";

        if (length(stats.seedLengths))
        {
            double _seedLengthSum       = std::accumulate(stats.seedLengths.begin(), stats.seedLengths.end(), 0.0);
            double seedLengthMean       = _seedLengthSum / stats.seedLengths.size();

            double _seedLengthMeanSqSum = std::inner_product(stats.seedLengths.begin(),
                                                            stats.seedLengths.end(),
                                                            stats.seedLengths.begin(),
                                                            0.0);
            double seedLengthStdDev     = std::sqrt(_seedLengthMeanSqSum / stats.seedLengths.size() -
                                                    seedLengthMean * seedLengthMean);
            uint16_t seedLengthMax      = *std::max_element(stats.seedLengths.begin(), stats.seedLengths.end());

            std::cout << "SeedStats:\n"
                    << " avgLength:   " << seedLengthMean << "\n"
                    << " stddev:      " << seedLengthStdDev << "\n"
                    << " max:         " << seedLengthMax << "\n\n";
        }
    #ifdef SEQAN_SIMD_ENABLED
        if (stats.numQueryWithExt > 0)
            std::cout << "Number of Extensions stats:\n"
                    << " # queries with Extensions:    " << stats.numQueryWithExt << "\n"
                    << " avg # extensions without Ali: " << stats.numExtScore / stats.numQueryWithExt << "\n"
                    << " avg # extensions with    Ali: " << stats.numExtAli / stats.numQueryWithExt << "\n\n";
    #endif
    #endif
    }

    if (options.verbosity >= 1)
    {
        auto const w = _numberOfDigits(stats.hitsFinal);
        std::cout << "Number of valid hits:                           "
                  << std::setw(w) << stats.hitsFinal
                  << "\nNumber of Queries with at least one valid hit:  "
                  << std::setw(w) << stats.qrysWithHit
                  << "\n";
    }

}

// ----------------------------------------------------------------------------
// struct GlobalDataHolder  -- one object per program
// ----------------------------------------------------------------------------

template <typename T>
inline T &
_initHelper(T & t1, T &&)
{
//     std::cout << "FOO\n";
    return t1;
}

template <typename T, typename T2>
inline T2 &&
_initHelper(T &, T2 && t2)
{
//     std::cout << "BAR\n";
    return std::move(t2);
}

template <typename TRedAlph_,
          typename TIndexSpec_,
          typename TFileFormat,
          BlastProgram p,
          BlastTabularSpec h>
class GlobalDataHolder
{
public:
    using TRedAlph       = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TMatch         = Match<TRedAlph>;

    static constexpr BlastProgram blastProgram  = p;
    static constexpr bool indexIsBiFM           = std::is_same<TIndexSpec_, BidirectionalIndex<TFMIndexInBi<>>>::value;
    static constexpr bool indexIsFM             = std::is_same<TIndexSpec_, TFMIndex<>>::value || indexIsBiFM;
    static constexpr bool alphReduction         = !std::is_same<TransAlph<p>, TRedAlph>::value;

    /* Sequence storage types */
    using TStringTag    = Alloc<>;
#if defined(LAMBDA_MMAPPED_DB)
    using TDirectStringTag = MMap<>;
#else
    using TDirectStringTag = TStringTag;
#endif
    using TQryTag  = TStringTag;
    using TSubjTag = TDirectStringTag; // even if subjects were translated they are now loaded from disk

    /* untranslated query sequences (ONLY USED FOR SAM/BAM OUTPUT) */
    using TUntransQrySeqs = StringSet<String<OrigQryAlph<p>, TQryTag>, Owner<ConcatDirect<>>>;

    /* Possibly translated but yet unreduced sequences */
    template <typename TSpec>
    using TTransSeqs     = StringSet<String<TransAlph<p>, TSpec>, Owner<ConcatDirect<>>>;
    using TTransQrySeqs  = TTransSeqs<TQryTag>;
    using TTransSubjSeqs = TTransSeqs<TSubjTag>;
    using TTransSubjReal = typename std::conditional<
                            alphReduction || indexIsFM,
                            TTransSubjSeqs,                     // real type
                            TTransSubjSeqs &>::type;            // will be initialized in constructor

    /* Reduced sequence objects, either as modstrings or as references to trans-strings */
    template <typename TSpec>
    using TRedAlphModString = ModifiedString<String<TransAlph<p>, TSpec>,
                                ModView<FunctorConvert<TransAlph<p>, TRedAlph>>>;

    using TRedQrySeqs   = typename std::conditional<
                            alphReduction,
                            StringSet<TRedAlphModString<TQryTag>, Owner<ConcatDirect<>>>, // modview
                            TTransQrySeqs &>::type;                                       // reference to owner
    using TRedSubjSeqs  = typename std::conditional<
                            alphReduction,
                            StringSet<TRedAlphModString<TSubjTag>, Owner<ConcatDirect<>>>, // modview
                            TTransSubjSeqs &>::type;                                       // reference to owner

    /* sequence ID strings */
    template <typename TSpec>
    using TIds          = StringSet<String<char, TSpec>, Owner<ConcatDirect<>>>;
    using TQryIds       = TIds<TQryTag>;
    using TSubjIds      = TIds<TSubjTag>;

    /* indeces and their type */
    using TIndexSpec    = TIndexSpec_;
    using TDbIndex      = Index<typename std::remove_reference<TRedSubjSeqs>::type, TIndexSpec>;

    /* output file */
    using TScoreScheme  = std::conditional_t<std::is_same<TRedAlph, Dna5>::value,
                                             Score<int, Simple>,
                                             Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>>;
//     using TScoreScheme  = TScoreScheme_;
    using TIOContext    = BlastIOContext<TScoreScheme, p, h>;
    using TFile         = FormattedFile<TFileFormat, Output, TIOContext>;
    using TBamFile      = FormattedFile<Bam, Output, BlastTabular>;

    /* misc types */
    using TPositions    = typename StringSetLimits<TTransQrySeqs>::Type;
    using TMasking      = StringSet<String<unsigned>, Owner<ConcatDirect<>>>;
    using TTaxIDs       = StringSet<String<uint32_t>, Owner<ConcatDirect<>>>;
    using TTaxParents   = String<uint32_t>;
    using TTaxHeights   = String<uint8_t>;
    using TTaxNames     = StringSet<CharString, Owner<ConcatDirect<>>>;

    /* the actual members */
    TDbIndex            dbIndex;

    TUntransQrySeqs     untranslatedQrySeqs;    // used iff outformat is sam or bam

    TTransQrySeqs       qrySeqs;
    TTransSubjReal      subjSeqs;

    TRedQrySeqs         redQrySeqs;
    TRedSubjSeqs        redSubjSeqs;

    TQryIds             qryIds;
    TSubjIds            subjIds;

    TFile               outfile;
    TBamFile            outfileBam;

    TPositions          untransQrySeqLengths;   // used iff qIsTranslated(p)
    TPositions          untransSubjSeqLengths;  // used iff sIsTranslated(p)

    TTaxIDs             sTaxIds;
    TTaxParents         taxParents;
    TTaxHeights         taxHeights;
    TTaxNames           taxNames;

    StatsHolder         stats;

    GlobalDataHolder() :
        subjSeqs(_initHelper(indexText(dbIndex), TTransSubjSeqs())),//std::integral_constant<bool, alphReduction || indexIsFM>())),// : TTransSubjSeqs()),
        redQrySeqs(qrySeqs),
        redSubjSeqs(subjSeqs),
        stats()
    {}
};

/* Documentation on the confusing type resolution used in the above class:
 *
 * !alphReduction && !indexIsFM  e.g. BLASTN and SA-Index
 *
 *   subjSeqs           is & and initialized with indexText()
 *   redSubjSeqs        is & and initialized with subjSeqs
 *   indexText(dbIndex) is non-ref owner StringSet assigned by loadDbIndexFromDisk()
 *
 * !alphReduction && indexIsFM  e.g. BLASTN and FM-Index
 *
 *   subjSeqs           is non-ref owner StringSet and assigned in loadSubjects()
 *   redSubjSeqs        is & and initialized with subjSeqs
 *   indexText(dbIndex) is non-ref owner StringSet, but never set (fmIndex doesnt need it)
 *
 * alphReduction && indexIsFM  e.g. default
 *
 *   subjSeqs           is non-ref owner StringSet and assigned in loadSubjects()
 *   redSubjSeqs        is lightweight reduced StringSet and initialized with subjSeqs
 *   indexText(dbIndex) is lightweight reduced StringSet, but never set (fmIndex doesnt need it)
 *
 * alphReduction && !indexIsFM  e.g. default
 *
 *   subjSeqs           is non-ref owner StringSet and assigned in loadSubjects()
 *   redSubjSeqs        is lightweight reduced StringSet and initialized with subjSeqs
 *   indexText(dbIndex) is lightweight reduced StringSet and assigned redSubjSeqs in loadDbIndexFromDisk
 */

// ----------------------------------------------------------------------------
// struct LocalDataHolder  -- one object per thread
// ----------------------------------------------------------------------------

template <typename TGlobalHolder_,
          typename TScoreExtension_>
class LocalDataHolder
{
public:
    using TGlobalHolder = TGlobalHolder_;
    using TRedQrySeq    = typename Value<typename std::remove_reference<typename TGlobalHolder::TRedQrySeqs>::type>::Type;
    using TSeeds        = StringSet<typename Infix<TRedQrySeq const>::Type>;
    using TSeedIndex    = Index<TSeeds, IndexSa<>>;
    using TMatch        = typename TGlobalHolder::TMatch;
    using TScoreExtension = TScoreExtension_;


    // references to global stuff
    LambdaOptions     const & options;
    TGlobalHolder /*const*/ & gH;
    static constexpr BlastProgram blastProgram = TGlobalHolder::blastProgram;

    // this is the localHolder for the i-th part of the queries
    uint64_t            i;
    uint64_t            nBlocks;

    // regarding range of queries
    uint64_t            indexBeginQry;
    uint64_t            indexEndQry;

    // regarding seedingp
    TSeeds              seeds;
    TSeedIndex          seedIndex;
//     std::forward_list<TMatch>   matches;
    std::vector<TMatch>   matches;
    std::vector<typename TMatch::TQId> seedRefs;  // mapping seed -> query
    std::vector<typename TMatch::TPos> seedRanks; // mapping seed -> relative rank

    // regarding extension
    using TAlignRow0 = Gaps<typename Infix<typename Value<typename TGlobalHolder::TTransQrySeqs>::Type>::Type,
                            ArrayGaps>;
    using TAlignRow1 = Gaps<typename Infix<typename Value<typename TGlobalHolder::TTransSubjSeqs>::Type>::Type,
                            ArrayGaps>;

#if (SEQAN_VERSION_MINOR < 4)
    using TDPContextNoSIMD = DPContext<typename Value<typename TGlobalHolder::TScoreScheme>::Type, TScoreExtension>;
#else
//     #if defined(SEQAN_SIMD_ENABLED)
//     using TCellValueSIMD  = typename SimdVector<int16_t>::TYPE;
//     using TDPCellSIMD     = DPCell_<TCellValueSIMD, TScoreExtension_>;
//     using TTraceValueSIMD = typename TraceBitMap_<TCellValueSIMD>::Type;
//     using TScoreHostSIMD  = String<TDPCellSIMD, Alloc<OverAligned> >;
//     using TTraceHostSIMD  = String<TTraceValueSIMD, Alloc<OverAligned> >;
//     using TDPContextSIMD  = DPContext<TDPCellSIMD, TTraceValueSIMD, TScoreHostSIMD, TTraceHostSIMD>;
//     #endif

    using TCellValueNoSIMD  = int16_t;
    using TDPCellNoSIMD     = DPCell_<TCellValueNoSIMD, TScoreExtension_>;
    using TTraceValueNoSIMD = typename TraceBitMap_<TCellValueNoSIMD>::Type;
    using TScoreHostNoSIMD  = String<TDPCellNoSIMD, Alloc<OverAligned> >;
    using TTraceHostNoSIMD  = String<TTraceValueNoSIMD, Alloc<OverAligned> >;
    using TDPContextNoSIMD  = DPContext<TDPCellNoSIMD, TTraceValueNoSIMD, TScoreHostNoSIMD, TTraceHostNoSIMD>;
#endif

    using TAliExtContext = AliExtContext_<TAlignRow0, TAlignRow1, TDPContextNoSIMD>;

    TAliExtContext      alignContext;
// #if defined(SEQAN_SIMD_ENABLED)
//     TDPContextSIMD      alignSIMDContext;
// #endif

    // map from sequence length to band size
    std::unordered_map<uint64_t, int> bandTable;

    // regarding the gathering of stats
    StatsHolder         stats;

    // progress string
    std::stringstream   statusStr;

    // constructor
    LocalDataHolder(LambdaOptions     const & _options,
                    TGlobalHolder     /*const*/ & _globalHolder) :
        options(_options), gH(_globalHolder), stats()
    {
        if (options.extensionMode == LambdaOptions::ExtensionMode::FULL_SIMD)
        {
            // division with rounding up
            nBlocks = (length(gH.redQrySeqs) + qNumFrames(blastProgram) * 10 - 1) / (qNumFrames(blastProgram) * 10);
        } else
        {
            nBlocks = length(gH.redQrySeqs) / qNumFrames(blastProgram);
        }
    }

    // copy constructor SHALLOW COPY ONLY, REQUIRED FOR firsprivate()
    LocalDataHolder(LocalDataHolder const & rhs) :
        options(rhs.options), gH(rhs.gH), stats()
    {
    }

    void init(uint64_t const _i)
    {
        i = _i;

        if (options.extensionMode == LambdaOptions::ExtensionMode::FULL_SIMD)
        {
            indexBeginQry = qNumFrames(blastProgram) * i * 10;
            indexEndQry = _min(qNumFrames(blastProgram) * (i+1) * 10, length(gH.qrySeqs));
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

#endif // LAMBDA_SEARCH_DATASTRUCTURES_H_
