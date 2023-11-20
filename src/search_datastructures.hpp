// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2020, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2020, Knut Reinert and Freie Universität Berlin
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
// search_datastructures.hpp: Data container structs
// ==========================================================================

#pragma once

#include <fmindex-collection/search/all.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>

#include <bio/ranges/views/convert.hpp>
#include <bio/ranges/views/deep.hpp>
#include <bio/ranges/views/translate_join.hpp>
#include <bio/ranges/views/type_reduce.hpp>

#include <seqan/align_extend.h>

#include "bisulfite_scoring.hpp"
#include "search_options.hpp"

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct Match
// ----------------------------------------------------------------------------

struct Match
{
    using TQId = uint64_t;
    using TSId = uint64_t;
    using TPos = uint64_t;

    TQId qryId;
    TSId subjId;
    TPos qryStart;
    TPos qryEnd;

    TPos subjStart;
    TPos subjEnd;

    constexpr friend auto operator<=>(Match const &, Match const &) = default;
};

// template <typename TAlph>
inline void setToSkip(Match & m)
{
    using TPos            = typename Match::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    m.qryStart            = posMax;
    m.subjStart           = posMax;
}

// template <typename TAlph>
inline bool isSetToSkip(Match const & m)
{
    using TPos            = typename Match::TPos;
    constexpr TPos posMax = std::numeric_limits<TPos>::max();
    return (m.qryStart == posMax) && (m.subjStart == posMax);
}

inline void _printMatch(Match const & m)
{
    std::cout << "MATCH  Query " << m.qryId << "(" << m.qryStart << ", " << m.qryEnd << ")   on Subject " << m.subjId
              << "(" << m.subjStart << ", " << m.subjEnd << ")" << std::endl
              << std::flush;
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

    // post-extension
    uint64_t hitsFailedExtendPercentIdentTest;
    uint64_t hitsFailedExtendBitScoreTest;
    uint64_t hitsFailedExtendEValueTest;
    uint64_t hitsAbundant;
    uint64_t hitsDuplicate;
    uint64_t hitsDuplicate2;

    // final
    uint64_t hitsFinal;
    uint64_t qrysWithHit;
    uint64_t pairs;

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
        hitsMerged       = 0;
        hitsTooShort     = 0;
        hitsMasked       = 0;

        hitsFailedPreExtendTest = 0;

        hitsFailedExtendPercentIdentTest = 0;
        hitsFailedExtendBitScoreTest     = 0;
        hitsFailedExtendEValueTest       = 0;
        hitsAbundant                     = 0;
        hitsDuplicate                    = 0;
        hitsDuplicate2                   = 0;

        hitsFinal   = 0;
        qrysWithHit = 0;
        pairs       = 0;

#ifdef LAMBDA_MICRO_STATS
        seedLengths.clear();
        timeGenSeeds    = 0;
        timeSearch      = 0;
        timeSort        = 0;
        timeExtend      = 0;
        timeExtendTrace = 0;

        numQueryWithExt = 0;
        numExtScore     = 0;
        numExtAli       = 0;
#endif
    }

    StatsHolder plus(StatsHolder const & rhs)
    {
        hitsAfterSeeding += rhs.hitsAfterSeeding;
        hitsMerged += rhs.hitsMerged;
        hitsTooShort += rhs.hitsTooShort;
        hitsMasked += rhs.hitsMasked;

        hitsFailedPreExtendTest += rhs.hitsFailedPreExtendTest;

        hitsFailedExtendPercentIdentTest += rhs.hitsFailedExtendPercentIdentTest;
        hitsFailedExtendBitScoreTest += rhs.hitsFailedExtendBitScoreTest;
        hitsFailedExtendEValueTest += rhs.hitsFailedExtendEValueTest;
        hitsAbundant += rhs.hitsAbundant;
        hitsDuplicate += rhs.hitsDuplicate;
        hitsDuplicate2 += rhs.hitsDuplicate2;

        hitsFinal += rhs.hitsFinal;
        qrysWithHit += rhs.qrysWithHit;
        pairs += rhs.pairs;

#ifdef LAMBDA_MICRO_STATS
        seqan::append(seedLengths, rhs.seedLengths);
        timeGenSeeds += rhs.timeGenSeeds;
        timeSearch += rhs.timeSearch;
        timeSort += rhs.timeSort;
        timeExtend += rhs.timeExtend;
        timeExtendTrace += rhs.timeExtendTrace;

        numQueryWithExt += rhs.numQueryWithExt;
        numExtScore += rhs.numExtScore;
        numExtAli += rhs.numExtAli;
#endif
        return *this;
    }

    StatsHolder operator+(StatsHolder const & rhs)
    {
        StatsHolder tmp(*this);
        return tmp.plus(rhs);
    }

    StatsHolder operator+=(StatsHolder const & rhs)
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
        auto const    w   = seqan::_numberOfDigits(rem); // number of digits
#define R  " " << std::setw(w)
#define RR " = " << std::setw(w)
#define BLANKS                                                                                                         \
    for (unsigned i = 0; i < w; ++i)                                                                                   \
        std::cout << " ";
        std::cout << "\033[1m   HITS                         ";
        BLANKS;
        std::cout << "Remaining\033[0m"
                  << "\n   after Seeding               ";
        BLANKS;
        std::cout << R << rem;
        if (stats.hitsMasked)
            std::cout << "\n - masked                   " << R << stats.hitsMasked << RR << (rem -= stats.hitsMasked);
        if (options.preScoring)
            std::cout << "\n - failed pre-extend test   " << R << stats.hitsFailedPreExtendTest << RR
                      << (rem -= stats.hitsFailedPreExtendTest);
        std::cout << "\n - failed e-value test      " << R << stats.hitsFailedExtendEValueTest << RR
                  << (rem -= stats.hitsFailedExtendEValueTest);
        std::cout << "\n - failed bitScore test     " << R << stats.hitsFailedExtendBitScoreTest << RR
                  << (rem -= stats.hitsFailedExtendBitScoreTest);
        std::cout << "\n - failed %-identity test   " << R << stats.hitsFailedExtendPercentIdentTest << RR
                  << (rem -= stats.hitsFailedExtendPercentIdentTest);
        std::cout << "\n - duplicates               " << R << stats.hitsDuplicate << RR << (rem -= stats.hitsDuplicate);
        std::cout << "\n - late duplicates          " << R << stats.hitsDuplicate2 << RR
                  << (rem -= stats.hitsDuplicate2);
        std::cout << "\n - abundant                 " << R << stats.hitsAbundant << "\033[1m" << RR
                  << (rem -= stats.hitsAbundant) << "\033[0m\n\n";

        if (rem != stats.hitsFinal)
            std::cout << "WARNING: hits don't add up\n";

#ifdef LAMBDA_MICRO_STATS
        std::cout << "Detailed Non-Wall-Clock times:\n"
                  << " genSeeds:    " << stats.timeGenSeeds << "\n"
                  << " search:      " << stats.timeSearch << "\n"
                  << " sort:        " << stats.timeSort << "\n"
                  << " extend:      " << stats.timeExtend << "\n"
                  << " extendTrace: " << stats.timeExtendTrace << "\n\n";

        if (seqan::length(stats.seedLengths))
        {
            double _seedLengthSum = std::accumulate(stats.seedLengths.begin(), stats.seedLengths.end(), 0.0);
            double seedLengthMean = _seedLengthSum / stats.seedLengths.size();

            double _seedLengthMeanSqSum =
              std::inner_product(stats.seedLengths.begin(), stats.seedLengths.end(), stats.seedLengths.begin(), 0.0);
            double seedLengthStdDev =
              std::sqrt(_seedLengthMeanSqSum / stats.seedLengths.size() - seedLengthMean * seedLengthMean);
            uint16_t seedLengthMax = *std::max_element(stats.seedLengths.begin(), stats.seedLengths.end());

            std::cout << "SeedStats:\n"
                      << " avgLength:   " << seedLengthMean << "\n"
                      << " stddev:      " << seedLengthStdDev << "\n"
                      << " max:         " << seedLengthMax << "\n\n";
        }
        if (stats.numQueryWithExt > 0)
            std::cout << "Number of Extensions stats:\n"
                      << " # queries with Extensions:    " << stats.numQueryWithExt << "\n"
                      << " avg # extensions without Ali: " << stats.numExtScore / stats.numQueryWithExt << "\n"
                      << " avg # extensions with    Ali: " << stats.numExtAli / stats.numQueryWithExt << "\n\n";
#endif
    }

    if (options.verbosity >= 1)
    {
        auto const w = seqan::_numberOfDigits(stats.hitsFinal);
        std::cout << "Number of total hits:                           " << std::setw(w) << stats.hitsFinal
                  << "\nNumber of Query-Subject pairs:                  " << std::setw(w) << stats.pairs
                  << "\nNumber of Queries with at least one valid hit:  " << std::setw(w) << stats.qrysWithHit << "\n";
    }
}

enum class IterativeSearchMode
{
    OFF,
    PHASE1,
    PHASE2
};

// ----------------------------------------------------------------------------
// struct GlobalDataHolder  -- one object per program
// ----------------------------------------------------------------------------

template <DbIndexType  c_dbIndexType_,
          AlphabetEnum c_origSbjAlph_,
          AlphabetEnum c_transAlph_,
          AlphabetEnum c_redAlph_,
          AlphabetEnum c_origQryAlph_>
class GlobalDataHolder
{
public:
    static constexpr DbIndexType  c_dbIndexType = c_dbIndexType_;
    static constexpr AlphabetEnum c_origSbjAlph = c_origSbjAlph_;
    static constexpr AlphabetEnum c_transAlph   = c_transAlph_;
    static constexpr AlphabetEnum c_redAlph     = c_redAlph_;
    static constexpr AlphabetEnum c_origQryAlph = c_origQryAlph_;

    using TOrigQryAlph = _alphabetEnumToType<c_origQryAlph>;
    using TOrigSbjAlph = _alphabetEnumToType<c_origSbjAlph>;
    using TTransAlph   = _alphabetEnumToType<c_transAlph>;
    using TRedAlph     = _alphabetEnumToType<c_redAlph>;
    using TMatch       = Match;

    // clang-format off
    static constexpr seqan::BlastProgram blastProgram =
        (c_transAlph != AlphabetEnum::AMINO_ACID)                                       ? seqan::BlastProgram::BLASTN  :
        (c_origQryAlph == AlphabetEnum::AMINO_ACID && c_origSbjAlph == c_origQryAlph)   ? seqan::BlastProgram::BLASTP  :
        (c_origQryAlph == AlphabetEnum::AMINO_ACID && c_origSbjAlph != c_origQryAlph)   ? seqan::BlastProgram::TBLASTN :
        (c_origQryAlph != AlphabetEnum::AMINO_ACID && c_origSbjAlph == c_origQryAlph)   ? seqan::BlastProgram::TBLASTX :
      /*(c_origQryAlph != AlphabetEnum::AMINO_ACID && c_origSbjAlph != c_origQryAlph) ?*/ seqan::BlastProgram::BLASTX;
    // clang-format on

    /* untranslated query sequences (ONLY USED FOR SAM/BAM OUTPUT) */
    using TQrySeqs = std::vector<std::vector<TOrigQryAlph>>;
    using TSbjSeqs = TCDStringSet<std::vector<TOrigSbjAlph>>;

    /* Translated sequence objects, either as modstrings or as references to original strings */
    using TTransQrySeqs = decltype(std::declval<TQrySeqs &>() | qryTransView<c_origQryAlph, c_transAlph, c_redAlph>);
    using TTransSbjSeqs = decltype(std::declval<TSbjSeqs &>() | sbjTransView<c_origSbjAlph, c_transAlph, c_redAlph>);

    using TRedQrySeqs = decltype(std::declval<TTransQrySeqs &>() | redView<c_transAlph, c_redAlph>);
    using TRedSbjSeqs = decltype(std::declval<TTransSbjSeqs &>() | redView<c_transAlph, c_redAlph>);

    /* sequence ID strings */
    using TIds     = TCDStringSet<std::string>;
    using TQryIds  = std::vector<std::string>;
    using TSubjIds = TIds;

    /* index */
    using TIndexFile   = index_file<c_dbIndexType, c_origSbjAlph, c_redAlph>;
    using TIndex       = typename TIndexFile::TIndex;
    using TIndexCursor = fmindex_collection::select_cursor_t<TIndex>;

    /* output file */
    // SeqAn2 scoring scheme for local alignment of extended seeds. This can be adapted for bisulfite scoring.
    using TScoreSchemeAlign =
      std::conditional_t<bio::alphabet::nucleotide<TTransAlph>,
                         std::conditional_t<c_redAlph == AlphabetEnum::DNA3BS,
                                            seqan::Score<int, seqan::ScoreMatrix<seqan::Dna5, seqan::BisulfiteMatrix>>,
                                            seqan::Score<int, seqan::Simple>>,
                         seqan::Score<int, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecSelectable>>>;

    // SeqAn2 scoring scheme for blast statistics (does not work with bisulfite scoring scheme)
    using TScoreSchemeStats =
      std::conditional_t<bio::alphabet::nucleotide<TTransAlph>,
                         seqan::Score<int, seqan::Simple>,
                         seqan::Score<int, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecSelectable>>>;

    using TIOContext    = seqan::BlastIOContext<TScoreSchemeStats, blastProgram>;
    using TBlastTabFile = seqan::FormattedFile<seqan::BlastTabular, seqan::Output, TIOContext>;
    using TBlastRepFile = seqan::FormattedFile<seqan::BlastReport, seqan::Output, TIOContext>;
    using TBamFile      = seqan::FormattedFile<seqan::Bam, seqan::Output, seqan::BlastTabular>;

    // Number of frames for subject and query sequences
    // clang-format off
    static constexpr uint8_t qryNumFrames = (c_redAlph == AlphabetEnum::DNA3BS)  ? 4 :
                                            seqan::qIsTranslated(blastProgram)   ? 6 :
                                            seqan::qHasRevComp(blastProgram)     ? 2 : 1;
    static constexpr uint8_t sbjNumFrames = (c_redAlph == AlphabetEnum::DNA3BS)  ? 2 :
                                            seqan::sIsTranslated(blastProgram)   ? 6 :
                                            seqan::sHasRevComp(blastProgram)     ? 2 : 1;
    // clang-format on

    /* misc types */
    using TPositions = std::vector<size_t>;

    /* the actual members */
    TIndexFile indexFile;

    // origin sbjSeqs in indexFile
    TTransSbjSeqs transSbjSeqs = indexFile.seqs | sbjTransView<c_origSbjAlph, c_transAlph, c_redAlph>;
    TRedSbjSeqs   redSbjSeqs   = transSbjSeqs | redView<c_transAlph, c_redAlph>;

    /* the following are not used when the query is loaded on-demand / lazily */
    TQryIds       qryIds;
    TQrySeqs      qrySeqs;
    TTransQrySeqs transQrySeqs = qrySeqs | qryTransView<c_origQryAlph, c_transAlph, c_redAlph>;
    TRedQrySeqs   redQrySeqs   = transQrySeqs | redView<c_transAlph, c_redAlph>;

    TBlastTabFile outfileBlastTab;
    TBlastRepFile outfileBlastRep;
    TBamFile      outfileBam;

    TScoreSchemeAlign scoringSchemeAlign;
    TScoreSchemeAlign scoringSchemeAlignBSRev;

    StatsHolder stats{};

    size_t              queryTotal        = 0;
    std::atomic<size_t> queryCount        = 0;
    size_t              records_per_batch = -1; // this is always overwritten

    GlobalDataHolder()                                     = default;
    GlobalDataHolder(GlobalDataHolder const &)             = delete;
    GlobalDataHolder(GlobalDataHolder &&)                  = delete;
    GlobalDataHolder & operator=(GlobalDataHolder const &) = delete;
    GlobalDataHolder & operator=(GlobalDataHolder &&)      = delete;
};

// ----------------------------------------------------------------------------
// struct LocalDataHolder  -- one object per thread
// ----------------------------------------------------------------------------

template <typename TGlobalHolder_>
class LocalDataHolder
{
public:
    using TGlobalHolder   = TGlobalHolder_;
    using TMatch          = typename TGlobalHolder::TMatch;
    using TScoreExtension = seqan::AffineGaps;
    using TSeqInfix0 = decltype(std::declval<std::ranges::range_reference_t<typename TGlobalHolder::TTransQrySeqs>>() |
                                bio::views::slice(0, 1));
    using TSeqInfix1 = decltype(std::declval<std::ranges::range_reference_t<typename TGlobalHolder::TTransSbjSeqs>>() |
                                bio::views::slice(0, 1));

    // references to global stuff
    LambdaOptions const &                options;
    TGlobalHolder /*const*/ &            gH;
    static constexpr seqan::BlastProgram blastProgram = TGlobalHolder::blastProgram;

    // temporary storage used by lazy-loading
    typename TGlobalHolder::TQryIds                qryIdsTmp;
    typename TGlobalHolder::TQrySeqs               qrySeqsTmp;
    // refers to temporary storage (lazy-loading) or subset of qry in globalHolder (eager-loading)
    decltype(qryIdsTmp | bio::views::slice(0, 0))  qryIds;
    decltype(qrySeqsTmp | bio::views::slice(0, 0)) qrySeqs;
    typename TGlobalHolder::TTransQrySeqs          transQrySeqs =
      qrySeqs | qryTransView<TGlobalHolder::c_origQryAlph, TGlobalHolder::c_transAlph, TGlobalHolder::c_redAlph>;
    typename TGlobalHolder::TRedQrySeqs redQrySeqs =
      transQrySeqs | redView<TGlobalHolder::c_transAlph, TGlobalHolder::c_redAlph>;
    size_t           queryCount = 0;
    std::vector<int> successfulQueries; // vector<bool> -.-

    // regarding seeding
    std::vector<typename TGlobalHolder::TIndexCursor>                    cursor_buffer;
    std::vector<std::pair<typename TGlobalHolder::TIndexCursor, size_t>> cursor_tmp_buffer;
    std::vector<std::pair<typename TGlobalHolder::TIndexCursor, size_t>> cursor_tmp_buffer2;
    std::vector<size_t>                                                  offset_modifier_buffer;
    std::vector<std::pair<size_t, size_t>>                               matches_buffer;
    std::vector<TMatch>                                                  matches;
    // allows local override of options for iterative search
    IterativeSearchMode                                                  iterativeSearch = IterativeSearchMode::OFF;
    LambdaOptions::SearchOpts                                            searchOpts;

    // regarding extension
    using TAlignRow0   = seqan::Gaps<TSeqInfix0, seqan::ArrayGaps>;
    using TAlignRow1   = seqan::Gaps<TSeqInfix1, seqan::ArrayGaps>;
    using TBlastMatch  = seqan::BlastMatch<TAlignRow0,
                                          TAlignRow1,
                                          uint32_t,
                                          std::string_view,
                                          std::string_view,
                                          std::vector<std::string>, // not used
                                          std::span<uint32_t>>;
    using TBlastRecord = seqan::BlastRecord<TBlastMatch,
                                            std::string_view,
                                            std::vector<std::string>, // not used
                                            std::string_view,
                                            uint32_t>;
    std::list<TBlastMatch>       blastMatches;
    std::unordered_set<uint64_t> uniqSubjIds;

    // regarding the gathering of stats
    StatsHolder stats{};

    // currently used search scheme
    search_schemes::Scheme searchScheme0; // iterative phase 1 scheme
    search_schemes::Scheme searchScheme1; // iterative phase 2 scheme
    search_schemes::Scheme searchScheme;  // the one currently in use

    LocalDataHolder()                                    = delete;
    LocalDataHolder(LocalDataHolder const &)             = delete;
    LocalDataHolder(LocalDataHolder &&)                  = delete;
    LocalDataHolder & operator=(LocalDataHolder const &) = delete;
    LocalDataHolder & operator=(LocalDataHolder &&)      = delete;

    LocalDataHolder(LambdaOptions const & _options, TGlobalHolder & _gH) :
      options{_options},
      gH{_gH},
      searchOpts{options.searchOpts},
      stats{},
      searchScheme0{search_schemes::expand(search_schemes::generator::pigeon_opt(0, options.searchOpts0.maxSeedDist),
                                           options.searchOpts0.seedLength)},
      searchScheme1{search_schemes::expand(search_schemes::generator::pigeon_opt(0, options.searchOpts.maxSeedDist),
                                           options.searchOpts.seedLength)},
      searchScheme{searchScheme1}
    {
        if (_options.iterativeSearch)
        {
            iterativeSearch = IterativeSearchMode::PHASE1;
            searchOpts      = options.searchOpts0;
            searchScheme    = searchScheme0;
        }
    }

    void reset()
    {
        // clear storage
        queryCount = 0;
        matches.clear();
        qryIdsTmp.clear();
        qrySeqsTmp.clear();
        blastMatches.clear();

        // stats explicitly not cleared, because accumulated
    }

    void resetViews()
    {
        transQrySeqs =
          qrySeqs | qryTransView<TGlobalHolder::c_origQryAlph, TGlobalHolder::c_transAlph, TGlobalHolder::c_redAlph>;
        redQrySeqs = transQrySeqs | redView<TGlobalHolder::c_transAlph, TGlobalHolder::c_redAlph>;
    }
};
