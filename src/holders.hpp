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
// holders.hpp: Data container structs
// ==========================================================================

#ifndef SEQAN_LAMBDA_HOLDERS_H_
#define SEQAN_LAMBDA_HOLDERS_H_

#include <seqan/align_extend.h>

#include "options.hpp"
#include "match.hpp"

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

void printStats(StatsHolder const & stats, LambdaOptions const & options)
{
    if ((options.verbosity >= 1) && options.isTerm && options.doubleIndexing)
        for (unsigned char i=0; i < options.threads + 3; ++i)
            std::cout << std::endl;

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
        std::cout << "\n - masked                   " << R << stats.hitsMasked
                  << RR << (rem -= stats.hitsMasked);
        std::cout << "\n - merged                   " << R << stats.hitsMerged
                  << RR << (rem -= stats.hitsMerged);
        std::cout << "\n - putative duplicates      " << R
                  << stats.hitsPutativeDuplicate << RR
                  << (rem -= stats.hitsPutativeDuplicate);
        std::cout << "\n - putative abundant        " << R
                  << stats.hitsPutativeAbundant   << RR
                  << (rem -= stats.hitsPutativeAbundant);
        std::cout << "\n - failed pre-extend test   " << R
                  << stats.hitsFailedPreExtendTest  << RR
                  << (rem -= stats.hitsFailedPreExtendTest);
        std::cout << "\n - failed %-identity test   " << R
                  << stats.hitsFailedExtendPercentIdentTest << RR
                  << (rem -= stats.hitsFailedExtendPercentIdentTest);
        std::cout << "\n - failed e-value test      " << R
                  << stats.hitsFailedExtendEValueTest << RR
                  << (rem -= stats.hitsFailedExtendEValueTest);
        std::cout << "\n - abundant                 " << R
                  << stats.hitsAbundant << RR
                  << (rem -= stats.hitsAbundant);
        std::cout << "\n - duplicates               " << R
                  << stats.hitsDuplicate << "\033[1m" << RR
                  << (rem -= stats.hitsDuplicate)
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
    static constexpr bool indexIsFM             = std::is_same<TIndexSpec_, TFMIndex<>>::value;
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

    TMasking            segIntStarts;
    TMasking            segIntEnds;

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
          typename TScoreExtension>
class LocalDataHolder
{
public:
    using TGlobalHolder = TGlobalHolder_;
    using TRedQrySeq    = typename Value<typename std::remove_reference<typename TGlobalHolder::TRedQrySeqs>::type>::Type;
    using TSeeds        = StringSet<typename Infix<TRedQrySeq const>::Type>;
    using TSeedIndex    = Index<TSeeds, IndexSa<>>;
    using TMatch        = typename TGlobalHolder::TMatch;


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
    std::vector<typename TMatch::TQId> seedRefs;  // mapping seed -> query
    std::vector<typename TMatch::TPos> seedRanks; // mapping seed -> relative rank

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

#endif // SEQAN_LAMBDA_HOLDERS_H_
