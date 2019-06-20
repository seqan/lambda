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

#pragma once

#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/deep.hpp>

#include <seqan/align_extend.h>

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct Match
// ----------------------------------------------------------------------------

// template<typename TAlph>
struct Match
{
//     typedef SizeTypeNum_<TAlph>    TQId;
//     typedef SizeTypeNum_<TAlph>    TSId;
//     typedef SizeTypePos_<TAlph>    TPos;
    using TQId = uint64_t;
    using TSId = uint64_t;
    using TPos = uint64_t;

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

// template <typename TAlph>
inline void
setToSkip(Match & m)
{
    using TPos          = typename Match::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    m.qryStart = posMax;
    m.subjStart = posMax;
}

// template <typename TAlph>
inline bool
isSetToSkip(Match const & m)
{
    using TPos          = typename Match::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    return (m.qryStart == posMax) && (m.subjStart == posMax);
}


inline void
_printMatch(Match const & m)
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
        auto const w = seqan::_numberOfDigits(rem); // number of digits
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
        auto const w = seqan::_numberOfDigits(stats.hitsFinal);
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

template <DbIndexType   c_dbIndexType,
          AlphabetEnum  c_origQryAlph,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph>
class GlobalDataHolder
{
public:
    using TOrigQryAlph   = _alphabetEnumToType<c_origQryAlph>;
    using TOrigSbjAlph   = _alphabetEnumToType<c_origSbjAlph>;
    using TTransAlph     = _alphabetEnumToType<c_transAlph>;
    using TRedAlph       = _alphabetEnumToType<c_redAlph>;
    using TMatch         = Match;

    static constexpr seqan::BlastProgram blastProgram =
        (c_transAlph != AlphabetEnum::AMINO_ACID)                                       ? seqan::BlastProgram::BLASTN  :
        (c_origQryAlph == AlphabetEnum::AMINO_ACID && c_origSbjAlph == c_origQryAlph)   ? seqan::BlastProgram::BLASTP  :
        (c_origQryAlph == AlphabetEnum::AMINO_ACID && c_origSbjAlph != c_origQryAlph)   ? seqan::BlastProgram::TBLASTN :
        (c_origQryAlph != AlphabetEnum::AMINO_ACID && c_origSbjAlph == c_origQryAlph)   ? seqan::BlastProgram::TBLASTX :
      /*(c_origQryAlph != AlphabetEnum::AMINO_ACID && c_origSbjAlph != c_origQryAlph) ?*/ seqan::BlastProgram::BLASTX;



//     static constexpr bool indexIsBiFM           = std::is_same<TIndexSpec_, BidirectionalIndex<TFMIndexInBi<>>>::value;
//     static constexpr bool indexIsFM             = std::is_same<TIndexSpec_, TFMIndex<>>::value || indexIsBiFM;
//     static constexpr bool alphReduction         = !std::is_same<TransAlph<p>, TRedAlph>::value;

    /* Sequence storage types */
//     using TStringTag    = Alloc<>;
// #if defined(LAMBDA_MMAPPED_DB)
//     using TDirectStringTag = MMap<>;
// #else
//     using TDirectStringTag = TStringTag;
// #endif
//     using TQryTag  = TStringTag;
//     using TSubjTag = TDirectStringTag; // even if subjects were translated they are now loaded from disk

    /* untranslated query sequences (ONLY USED FOR SAM/BAM OUTPUT) */
    using TUntransQrySeqs = TCDStringSet<std::vector<TOrigQryAlph>>;

    /* Possibly translated but yet unreduced sequences */
    using TTransQrySeqs  = TCDStringSet<std::vector<TTransAlph>>;
    using TTransSubjSeqs = TCDStringSet<std::vector<TTransAlph>>;

    using TSeqAn2TransSeq = seqan::String<std::conditional_t<c_transAlph == AlphabetEnum::AMINO_ACID,
                                                      seqan::AminoAcid,
                                                      seqan::Dna5>>;

    /* Reduced sequence objects, either as modstrings or as references to trans-strings */
    template <typename TSpec>
    using TRedAlphModString =
        decltype(std::declval<TSpec const &> | seqan3::view::deep{seqan3::view::convert<TRedAlph>});

    using TRedQrySeqs   = std::conditional_t<c_transAlph == c_redAlph,
                                             TTransQrySeqs &,                          // reference to owner
                                             TRedAlphModString<TTransQrySeqs>>;        // modview
    using TRedSbjSeqs   = std::conditional_t<c_transAlph == c_redAlph,
                                             TTransSubjSeqs &,                          // reference to owner
                                             TRedAlphModString<TTransSubjSeqs>>;        // modview

    /* sequence ID strings */
    using TIds          = TCDStringSet<std::string>;
    using TQryIds       = TIds;
    using TSubjIds      = TIds;

    /* indeces and their type */
    using TDbIndex      = std::conditional_t<c_dbIndexType == DbIndexType::BI_FM_INDEX,
                                             seqan3::bi_fm_index<true>,
                                             seqan3::fm_index<true>>;

    /* output file */
//     using TScoreScheme3  = std::conditional_t<seqan3::NucleotideAlphabet<TRedAlph>,
//                                               seqan3::nucleotide_scoring_scheme<>,
//                                               seqan3::aminoacid_scoring_scheme<>>;

    using TScoreScheme  =
        std::conditional_t<seqan3::NucleotideAlphabet<TRedAlph>,
                           seqan::Score<int, seqan::Simple>,
                           seqan::Score<int, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecSelectable>>>;
    using TIOContext    = seqan::BlastIOContext<TScoreScheme, blastProgram>;
    using TBlastTabFile = seqan::FormattedFile<seqan::BlastTabular, seqan::Output, TIOContext>;
    using TBlastRepFile = seqan::FormattedFile<seqan::BlastReport, seqan::Output, TIOContext>;
    using TBamFile      = seqan::FormattedFile<seqan::Bam, seqan::Output, seqan::BlastTabular>;

    /* misc types */
    using TPositions    = std::vector<size_t>;

    /* the actual members */
    index_file<c_dbIndexType, c_origSbjAlph, c_transAlph> indexFile;

    TRedSbjSeqs         redSubjSeqs;

    TUntransQrySeqs     untranslatedQrySeqs;    // used iff outformat is sam or bam
    TPositions          untransQrySeqLengths;   // used iff qIsTranslated(p)
    TTransQrySeqs       qrySeqs;
    TRedQrySeqs         redQrySeqs;
    TQryIds             qryIds;

    TBlastTabFile       outfileBlastTab;
    TBlastRepFile       outfileBlastRep;
    TBamFile            outfileBam;

    StatsHolder         stats;

    GlobalDataHolder() :
        redQrySeqs(qrySeqs),
        redSubjSeqs(indexFile.transSeqs),
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
    using TMatch        = typename TGlobalHolder::TMatch;
    using TScoreExtension = TScoreExtension_;


    // references to global stuff
    LambdaOptions     const & options;
    TGlobalHolder /*const*/ & gH;
    static constexpr seqan::BlastProgram blastProgram = TGlobalHolder::blastProgram;

    // this is the localHolder for the i-th part of the queries
    uint64_t            i;
    uint64_t            nBlocks;

    // regarding range of queries
    uint64_t            indexBeginQry;
    uint64_t            indexEndQry;

    // regarding seedingp
    std::vector<TMatch>   matches;
    std::vector<typename TMatch::TQId> seedRefs;  // mapping seed -> query
    std::vector<typename TMatch::TPos> seedRanks; // mapping seed -> relative rank

    // regarding extension
    using TAlignRow = seqan::Gaps<seqan::Infix<typename TGlobalHolder::TSeqAn2TransSeq>,
                                               seqan::ArrayGaps>;

#if (SEQAN_VERSION_MINOR < 4)
    using TDPContextNoSIMD = seqan::DPContext<typename seqan::Value<typename TGlobalHolder::TScoreScheme>::Type, TScoreExtension>;
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
    using TDPCellNoSIMD     = seqan::DPCell_<TCellValueNoSIMD, TScoreExtension_>;
    using TTraceValueNoSIMD = typename seqan::TraceBitMap_<TCellValueNoSIMD>::Type;
    using TScoreHostNoSIMD  = seqan::String<TDPCellNoSIMD, seqan::Alloc<seqan::OverAligned> >;
    using TTraceHostNoSIMD  = seqan::String<TTraceValueNoSIMD, seqan::Alloc<seqan::OverAligned> >;
    using TDPContextNoSIMD  = seqan::DPContext<TDPCellNoSIMD, TTraceValueNoSIMD, TScoreHostNoSIMD, TTraceHostNoSIMD>;
#endif

//     using TAliExtContext = AliExtContext_<TAlignRow0, TAlignRow1, TDPContextNoSIMD>;

//     TAliExtContext      alignContext;
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

        matches.clear();
        seedRefs.clear();
        seedRanks.clear();
//         stats.clear();
        statusStr.clear();
        statusStr.precision(2);
    }
};

