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
// search_algo.hpp: functions to set up, perform and evaluate the search
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

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate_join.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

#include "bisulfite_scoring.hpp"
#include "evaluate_bisulfite_alignment.hpp"
#include "view_dna_n_to_random.hpp"
#include "view_reduce_to_bisulfite.hpp"

#include <fmindex-collection/DenseCSA.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/locate.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/expand.h>

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

template <DbIndexType   c_dbIndexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void
prepareScoring(GlobalDataHolder<c_dbIndexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>  & globalHolder,
               LambdaOptions const & options)
{
    if constexpr (c_transAlph != AlphabetEnum::AMINO_ACID)
    {
        // Seqan2
        seqan::setScoreMatch(context(globalHolder.outfileBlastTab).scoringScheme, options.match);
        seqan::setScoreMismatch(context(globalHolder.outfileBlastTab).scoringScheme, options.misMatch);

        if constexpr (c_redAlph == AlphabetEnum::DNA3BS)
        {
            // Seqan2
            seqan::setScoreBisulfiteMatrix(globalHolder.scoringSchemeAlign, options.match, options.misMatch, bsDirection::fwd);
            seqan::setScoreBisulfiteMatrix(globalHolder.scoringSchemeAlignBSRev, options.match, options.misMatch, bsDirection::rev);

            // Seqan3
            globalHolder.scoringSchemePreScoring.set_bisulfite_scheme(seqan3::match_score{options.match},
                                                         seqan3::mismatch_score{options.misMatch});
            globalHolder.scoringSchemePreScoringBSRev.set_bisulfite_scheme(seqan3::match_score{options.match},
                                                         seqan3::mismatch_score{options.misMatch}, bsDirection::rev);
        }
        else
        {
            // Seqan2
            globalHolder.scoringSchemeAlign = seqan::seqanScheme(context(globalHolder.outfileBlastTab).scoringScheme);
            globalHolder.scoringSchemeAlignBSRev = seqan::seqanScheme(context(globalHolder.outfileBlastTab).scoringScheme);

            // Seqan3
            globalHolder.scoringSchemePreScoring.set_simple_scheme(seqan3::match_score{options.match},
                                                         seqan3::mismatch_score{options.misMatch});
            globalHolder.scoringSchemePreScoringBSRev.set_simple_scheme(seqan3::match_score{options.match},
                                                         seqan3::mismatch_score{options.misMatch});
        }
    }
    else
    {
        seqan::AminoAcidScoreMatrixID       seqan2_matrix_id{};
        seqan3::aminoacid_similarity_matrix seqan3_matrix_id{};

        switch (options.scoringMethod)
        {
            case 45:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM45;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::blosum45;
                break;
            case 62:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM62;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::blosum62;
                break;
            case 80:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM80;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::blosum80;
                break;
            default:
                break;
        }

        // seqan2
        seqan::setScoreMatrixById(seqan::context(globalHolder.outfileBlastTab).scoringScheme._internalScheme,
                                  seqan2_matrix_id);
        seqan::setScoreMatrixById(globalHolder.scoringSchemeAlign,
                                  seqan2_matrix_id);
        seqan::setScoreMatrixById(globalHolder.scoringSchemeAlignBSRev,
                                  seqan2_matrix_id);
        // Seqan3
        globalHolder.scoringSchemePreScoring.set_similarity_matrix(seqan3_matrix_id);
        globalHolder.scoringSchemePreScoringBSRev.set_similarity_matrix(seqan3_matrix_id);
    }

    // seqan2
    seqan::setScoreGapOpenBlast(seqan::context(globalHolder.outfileBlastTab).scoringScheme, options.gapOpen);
    seqan::setScoreGapExtend(seqan::context(globalHolder.outfileBlastTab).scoringScheme, options.gapExtend);

    seqan::setScoreGapOpen(globalHolder.scoringSchemeAlign, options.gapOpen + options.gapExtend);
    seqan::setScoreGapExtend(globalHolder.scoringSchemeAlign, options.gapExtend);

    seqan::setScoreGapOpen(globalHolder.scoringSchemeAlignBSRev, options.gapOpen + options.gapExtend);
    seqan::setScoreGapExtend(globalHolder.scoringSchemeAlignBSRev, options.gapExtend);

    if (!seqan::isValid(seqan::context(globalHolder.outfileBlastTab).scoringScheme))
        throw std::runtime_error{"Could not compute Karlin-Altschul-Values for Scoring Scheme.\n"};
}

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void
loadDbIndexFromDisk(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
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

    globalHolder.transSbjSeqs = globalHolder.indexFile.seqs | sbjTransView<c_origSbjAlph, c_transAlph, c_redAlph>;
    globalHolder.redSbjSeqs =   globalHolder.transSbjSeqs | redView<c_transAlph, c_redAlph>;

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    // this is actually part of prepareScoring(), but the values are just available now
    if constexpr (c_origSbjAlph != c_transAlph)
    {
        // last value has sum of lengths
        // seqan::context(globalHolder.outfileBlastTab).dbTotalLength  = globalHolder.indexFile.origSeqLengths.back();
        // seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = globalHolder.indexFile.origSeqLengths.size() - 1;
        seqan::context(globalHolder.outfileBlastTab).dbTotalLength  =
            std::accumulate(globalHolder.indexFile.seqs.begin(), globalHolder.indexFile.seqs.end(), 0,
                [](size_t sum, auto const & a)
                {
                    return sum + std::ranges::size(a);
                });
        seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = std::ranges::size(globalHolder.indexFile.seqs);
    } else
    {
        seqan::context(globalHolder.outfileBlastTab).dbTotalLength  =
            std::ranges::distance(globalHolder.redSbjSeqs | std::views::join);
        seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = std::ranges::size(globalHolder.redSbjSeqs);
    }

    seqan::context(globalHolder.outfileBlastTab).dbName = options.indexFilePath;
    //TODO did we forget outfileBlastRep here? Or is that copied in output?
}

// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void
countQuery(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>       & globalHolder,
           LambdaOptions                                                                       const & options)
{
    double start = sysTime();

    std::string strIdent = "Counting Query Sequences ...";
    myPrint(options, 1, strIdent);

    // TODO potentially optimise this for fasta/fastq with simple 'grep -c'
    seqan3::sequence_file_input<QueryFileTraits<c_origQryAlph>, seqan3::fields<>> infile{options.queryFile};

    // parse the file completely and get count in one line:
    globalHolder.queryTotal = std::ranges::distance(infile);

    // batch-size as set in options (unless too few sequences)
    globalHolder.records_per_batch = std::max<size_t>(std::min<size_t>(globalHolder.queryTotal / (options.threads * 10),
                                                                       options.maximumQueryBlockSize),
                                                      1);
    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function seedLooksPromising()
// --------------------------------------------------------------------------

// perform a fast local alignment score calculation on the seed and see if we
// reach above threshold
// WARNING the following function only works for hammingdistanced seeds
template <typename TGlobalHolder>
inline bool
seedLooksPromising(LocalDataHolder<TGlobalHolder> const & lH,
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
                            static_cast<uint64_t>(std::ranges::size(lH.transQrySeqs[m.qryId]) - effectiveQBegin),
                            static_cast<uint64_t>(std::ranges::size(lH.gH.transSbjSeqs[m.subjId]) - effectiveSBegin),
                            effectiveLength});
    }

    auto const & qSeq = lH.transQrySeqs[m.qryId]
                      | seqan3::views::slice(effectiveQBegin, static_cast<int64_t>(effectiveQBegin + effectiveLength));
    auto const & sSeq = lH.gH.transSbjSeqs[m.subjId]
                      | seqan3::views::slice(effectiveSBegin, static_cast<int64_t>(effectiveSBegin + effectiveLength));

    int             s = 0;
    int      maxScore = 0;
    int const thresh  = lH.options.preScoringThresh * effectiveLength;

    // score the diagonal
    auto & currentScoringScheme = TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS && m.subjId % 2 ?
                                  lH.gH.scoringSchemePreScoringBSRev :  lH.gH.scoringSchemePreScoring;

    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        s += currentScoringScheme.score(qSeq[i], sSeq[i]);

        if (s < 0)
            s = 0;
        else if (s > maxScore)
            maxScore = s;

        if (maxScore >= thresh)
            return true;
    }
    return false;
}

template <DbIndexType   c_indexType,
          typename      TGlobalHolder,
          typename      TSeed>
inline void
search_impl(LocalDataHolder<TGlobalHolder> & lH, TSeed && seed)
{
    assert(lH.queries.size() == 1);
    assert(lH.queries[0].size() == seed.size());

    if constexpr (c_indexType == DbIndexType::FM_INDEX)
    {
        // prepare query for fmindex_collection, by shifting by one and reversing the string
        for (size_t i{0}; i < seed.size(); ++i)
        {
            lH.queries[0][seed.size()-i-1] = seed[i].to_rank()+1;
        }

        fmindex_collection::search_backtracking::search(
            lH.gH.indexFile.index,
            lH.queries,
            lH.options.maxSeedDist,
            [&](size_t /*queryId*/, auto cursor, size_t /*errors*/)
            {
                lH.cursor_buffer.push_back(cursor);
            }
        );
    }
    else if constexpr (c_indexType == DbIndexType::BI_FM_INDEX)
    {
        // prepare query for fmindex_collection, by shifting by one value
        for (size_t i{0}; i < seed.size(); ++i)
        {
            lH.queries[0][i] = seed[i].to_rank()+1;
        }

        fmindex_collection::search_pseudo::search</*editdistance=*/false>(
            lH.gH.indexFile.index,
            lH.queries,
            lH.searchScheme,
            [&](size_t /*queryId*/, auto cursor, size_t /*errors*/)
            {
                lH.cursor_buffer.push_back(cursor);
            }
        );
    }
}

template <typename      TGlobalHolder,
          typename      TSeed>
inline void
searchHalfExactImpl(LocalDataHolder<TGlobalHolder> & lH, TSeed && seed)
{

    auto extendRight = [&](auto & cursor, auto c)
    {
        auto newCursor = cursor.extendRight(c.to_rank()+1);
        if (newCursor.empty())
        {
            return false;
        }
        cursor = newCursor;
        return true;
    };

    using alph_t = std::ranges::range_value_t<TSeed>;

    lH.cursor_tmp_buffer.clear();
    lH.cursor_tmp_buffer2.clear();
    size_t const seedFirstHalfLength = lH.options.seedLength / 2;
    size_t const seedSecondHalfLength = lH.options.seedLength - seedFirstHalfLength;

    lH.cursor_tmp_buffer.emplace_back(lH.gH.indexFile.index, 0);
    auto & c = lH.cursor_tmp_buffer.back().first;

    // extend by half exactly
    for (size_t i = 0; i < seedFirstHalfLength; ++i)
    {
        auto newCursor = c.extendRight(seed[i].to_rank()+1);
        if (newCursor.empty())
        {
            return;
        }
        c = newCursor;
    }


    // manual backtracking
    for (size_t i = 0; i < seedSecondHalfLength; ++i)
    {
        auto seed_at_i = seed[seedFirstHalfLength + i];

        for (auto & [ cursor, error_count ] : lH.cursor_tmp_buffer)
        {
            if (error_count < lH.options.maxSeedDist)
            {
                for (size_t r = 0; r < seqan3::alphabet_size<alph_t>; ++r)
                {
                    alph_t cur_letter = seqan3::assign_rank_to(r, alph_t{});

                    lH.cursor_tmp_buffer2.emplace_back(cursor, error_count + (cur_letter != seed_at_i));
                    if (!extendRight(lH.cursor_tmp_buffer2.back().first, cur_letter))
                        lH.cursor_tmp_buffer2.pop_back();
                }
            }
            else
            {
                lH.cursor_tmp_buffer2.emplace_back(cursor, error_count);
                if (!extendRight(lH.cursor_tmp_buffer2.back().first, seed_at_i))
                    lH.cursor_tmp_buffer2.pop_back();
            }
        }
        lH.cursor_tmp_buffer.clear();
        std::swap(lH.cursor_tmp_buffer, lH.cursor_tmp_buffer2);
    }

    lH.cursor_buffer.reserve(lH.cursor_buffer.size() + lH.cursor_tmp_buffer.size());
    std::ranges::copy(lH.cursor_tmp_buffer | std::views::elements<0>, std::back_inserter(lH.cursor_buffer));
    lH.cursor_tmp_buffer.clear();
}


template <typename TGlobalHolder>
inline void
search(LocalDataHolder<TGlobalHolder> & lH)
{
    using TTransAlph = typename TGlobalHolder::TTransAlph;
    using TMatch = typename TGlobalHolder::TMatch;

    constexpr auto dbIndexType = TGlobalHolder::c_dbIndexType;

    /* The flag changes seeding to take into account how successful previous seeding operations were.
     * It basically changes the "heuristicFactor" below to be evidence-based (although the formula is also
     * slightly different).
     * This makes seeding more robust to different kinds of datasets and options, but due to multi-threading
     * it also makes the results becoming non-deterministic (i.e. repeated runs of the program may produce
     * slightly different amounts of results).
     */
#ifdef LAMBDA_NONDETERMINISTIC_SEEDS
    [[maybe_unused]] size_t const hitPerFinalHit = (double)lH.stats.hitsAfterSeeding / std::max<double>(lH.stats.hitsFinal, 1.0);
#endif

    size_t hitsThisSeq = 0;                     // hits per untranslated sequence (multiple frames counted together)
    size_t needlesSum = 0;                      // cumulative size of all frames of a sequence
    size_t needlesPos = 0;                      // current position in the cumulative length
    constexpr size_t heuristicFactor = 10;      // a vague approximation of the success rate of seeds

    /* These values are used to track how many hits we already have per untranslated sequence,
     * how many seeds we still have after the current one, and based on that, how many hits we
     * expect that the current seed needs to have.
     * With this information, the current seed can be elongated to get better (but fewer) hits.
     */

    for (size_t i = 0; i < std::ranges::size(lH.redQrySeqs); ++i)
    {
        if (lH.redQrySeqs[i].size() < lH.options.seedLength)
            continue;

        if (i % TGlobalHolder::qryNumFrames == 0) // reset on every "real" new read
        {
            hitsThisSeq = 0;
            needlesSum = 0;
            needlesPos = 0;
            for (size_t j = 0; j < TGlobalHolder::qryNumFrames; ++j)
                needlesSum += lH.redQrySeqs[i + j].size();
        }

        for (size_t seedBegin = 0; /* below */; seedBegin += lH.options.seedOffset)
        {
            // skip proteine 'X' or Dna 'N'
            while ((seedBegin <= (lH.redQrySeqs[i].size() - lH.options.seedLength)) &&
                   (lH.transQrySeqs[i][seedBegin] == seqan3::assign_char_to('`', TTransAlph{}))) // assume that '°' gets converted to UNKNOWN
                ++seedBegin;

            // termination criterium
            if (seedBegin > (lH.redQrySeqs[i].size() - lH.options.seedLength))
                break;

            // results are in cursor_buffer
            lH.cursor_buffer.clear();
            if (lH.options.seedHalfExact)
                searchHalfExactImpl(lH, lH.redQrySeqs[i] | seqan3::views::slice(seedBegin, seedBegin + lH.options.seedLength));
            else
                search_impl<dbIndexType>(lH, lH.redQrySeqs[i] | seqan3::views::slice(seedBegin, seedBegin + lH.options.seedLength));

            if (lH.options.adaptiveSeeding)
                lH.offset_modifier_buffer.clear();

            for (auto & cursor : lH.cursor_buffer)
            {
                size_t seedLength = lH.options.seedLength;

                // elongate seeds
                if (lH.options.adaptiveSeeding)
                {
#ifdef LAMBDA_NONDETERMINISTIC_SEEDS
                    size_t desiredOccs = hitsThisSeq >= lH.options.maxMatches * hitPerFinalHit
                                       ? 1
                                       : (lH.options.maxMatches * hitPerFinalHit - hitsThisSeq) /
                                std::max<size_t>((needlesSum - needlesPos - seedBegin) / lH.options.seedOffset, 1ul);

#else // lambda2 mode BUT NOT QUIET
                    // desiredOccs == the number of seed hits we estimate that we need to reach lH.options.maxMatches
                    // hitsThisSeq >= lH.options.maxMatches → if we have more than we need already, only look for one
                    // (lH.options.maxMatches - hitsThisSeq) * heuristicFactor → total desired hits
                    // ((needlesSum - needlesPos - seedBegin) / lH.options.seedOffset → number of remaining seeds
                    // dividing the last two yields desired hits FOR THE CURRENT SEED
                    size_t desiredOccs = hitsThisSeq >= lH.options.maxMatches
                                       ? 1
                                       : (lH.options.maxMatches - hitsThisSeq) * heuristicFactor /
                                std::max<size_t>((needlesSum - needlesPos - seedBegin) / lH.options.seedOffset, 1ul);
#endif

                    if (desiredOccs == 0)
                        desiredOccs = 1;

                    // This aborts when we fall under the threshold
                    auto old_cursor = cursor;
                    size_t old_count = cursor.count();
                    while (seedBegin + seedLength < lH.redQrySeqs[i].size())
                    {
                        cursor = cursor.extendRight((lH.redQrySeqs[i][seedBegin + seedLength]).to_rank()+1);

                        size_t new_count = cursor.count();

                        if (new_count < desiredOccs && new_count < old_count) // we always continue to extend if we don't loose anything
                        {
                            // revert last extension
                            cursor = old_cursor;
                            break;
                        }

                        ++seedLength;
                        old_count = new_count;
                        old_cursor = cursor;
                    }
                }

                // discard over-abundant seeds (typically seeds at end of sequence that couldn't be elongated)
                if (cursor.count() > heuristicFactor * lH.options.maxMatches)
                    continue;

                // locate hits
                for (auto [subjNo, subjOffset] : fmindex_collection::LocateLinear{lH.gH.indexFile.index, cursor})
                {
                    // !TODO Should this be handled by the cursor?
                    if (dbIndexType == DbIndexType::FM_INDEX)
                    {
                        subjOffset -= seedLength;
                    }
                    TMatch m {static_cast<typename TMatch::TQId>(i),
                              static_cast<typename TMatch::TSId>(subjNo),
                              static_cast<typename TMatch::TPos>(seedBegin),
                              static_cast<typename TMatch::TPos>(seedBegin + seedLength),
                              static_cast<typename TMatch::TPos>(subjOffset),
                              static_cast<typename TMatch::TPos>(subjOffset + seedLength)};

                    ++lH.stats.hitsAfterSeeding;

                    if (!seedLooksPromising(lH, m))
                    {
                        ++lH.stats.hitsFailedPreExtendTest;
                    }
                    else
                    {
                        lH.matches.push_back(m);
                        ++hitsThisSeq;
#ifdef LAMBDA_MICRO_STATS
                        lH.stats.seedLengths.push_back(seedLength);
#endif

                    }
                }
            }
        }

        needlesPos += lH.redQrySeqs[i].size();
    }
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
    } else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        bm.qFrameShift = (m.qryId % 2) + 1;
        if (m.qryId % 4 > 1)
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
    } else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        bm.sFrameShift = (m.subjId % 2) + 1;
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
    auto const & const_gH = lH.gH;

    if (record.matches.size() > 0)
    {
        ++lH.stats.qrysWithHit;
        // sort and remove duplicates -> STL, yeah!
        auto const before = record.matches.size();

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
                if ((lH.gH.indexFile.sTaxIds[bm._n_sId].size() > 0) && (const_gH.indexFile.taxonParentIDs[lH.gH.indexFile.sTaxIds[bm._n_sId][0]] != 0))
                {
                    record.lcaTaxId = lH.gH.indexFile.sTaxIds[bm._n_sId][0];
                    break;
                }
            }

            if (record.lcaTaxId != 0)
                for (auto const & bm : record.matches)
                    for (uint32_t const sTaxId : lH.gH.indexFile.sTaxIds[bm._n_sId])
                        if (const_gH.indexFile.taxonParentIDs[sTaxId] != 0) // TODO do we want to skip unassigned subjects
                            record.lcaTaxId = computeLCA(const_gH.indexFile.taxonParentIDs,
                                                         const_gH.indexFile.taxonHeights,
                                                         sTaxId,
                                                         record.lcaTaxId);

            record.lcaId = const_gH.indexFile.taxonNames[record.lcaTaxId];
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

    bm.qEnd = lH.transQrySeqs[m.qryId].size();
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
    bm.sEnd = std::min<size_t>(bm.sStart + bm.qEnd - bm.qStart + band, lH.gH.transSbjSeqs[m.subjId].size());

    if (bm.sStart >= band)
        bm.sStart -= band;
    else
        bm.sStart = 0;

    seqan::assignSource(bm.alignRow0, lH.transQrySeqs[m.qryId] | seqan3::views::slice(bm.qStart, bm.qEnd));
    seqan::assignSource(bm.alignRow1, lH.gH.transSbjSeqs[m.subjId] | seqan3::views::slice(bm.sStart, bm.sEnd));
}

template <typename TBlastMatch,
          typename TLocalHolder>
inline auto
_untrueQryId(TBlastMatch const & bm,
             TLocalHolder const &)
{
    using namespace seqan;

    if constexpr (seqan::qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 6 + bm.qFrameShift - 1;
        else
            return bm._n_qId * 6 - bm.qFrameShift + 2;
    } else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 4;
        else
            return bm._n_qId * 4 + 2;
    } else if constexpr (seqan::qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
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

    if constexpr (seqan::sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.sFrameShift > 0)
            return bm._n_sId * 6 + bm.sFrameShift - 1;
        else
            return bm._n_sId * 6 - bm.sFrameShift + 2;
    } else if constexpr (seqan::sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram) ||
                         TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
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

    auto oldQLen = seqan::length(source(bm.alignRow0));
    auto oldSLen = seqan::length(source(bm.alignRow1));

    // replace source from underneath without triggering reset
    bm.alignRow0._source = lH.transQrySeqs[_untrueQryId(bm, lH)] | seqan3::views::slice(0, std::ranges::size(lH.transQrySeqs[_untrueQryId(bm, lH)]));
    bm.alignRow1._source = lH.gH.transSbjSeqs[_untrueSubjId(bm, lH)] | seqan3::views::slice(0, std::ranges::size(lH.gH.transSbjSeqs[_untrueSubjId(bm, lH)]));

    // insert fields into array gaps
    if (bm.alignRow0._array[0] == 0)
        bm.alignRow0._array[1] += bm.qStart;
    else
        insert(bm.alignRow0._array, 0, std::vector<uint64_t>{0, bm.qStart});
    if (bm.alignRow0._array[seqan::length(bm.alignRow0._array) - 1] == 0)
        bm.alignRow0._array[seqan::length(bm.alignRow0._array) - 2] += length(source(bm.alignRow0)) - oldQLen;
    else
        append(bm.alignRow0._array, std::vector<uint64_t>{seqan::length(source(bm.alignRow0)) - oldQLen, 0});

    if (bm.alignRow1._array[0] == 0)
        bm.alignRow1._array[1] += bm.sStart;
    else
        insert(bm.alignRow1._array, 0, std::vector<uint64_t>{0, bm.sStart});
    if (bm.alignRow1._array[length(bm.alignRow1._array) - 1] == 0)
        bm.alignRow1._array[length(bm.alignRow1._array) - 2] += seqan::length(source(bm.alignRow1)) - oldSLen;
    else
        append(bm.alignRow1._array, std::vector<uint64_t>{seqan::length(source(bm.alignRow1)) - oldSLen, 0});

    // the begin positions from the align object are relative to the infix created above
    bm.qEnd   = bm.qStart + seqan::endPosition(bm.alignRow0);
    bm.qStart = bm.qStart + seqan::beginPosition(bm.alignRow0);
    bm.sEnd   = bm.sStart + seqan::endPosition(bm.alignRow1);
    bm.sStart = bm.sStart + seqan::beginPosition(bm.alignRow1);

    // set clipping positions on new gaps objects
    seqan::setBeginPosition(bm.alignRow0, bm.qStart);
    seqan::setEndPosition(bm.alignRow0, bm.qEnd);
    seqan::setBeginPosition(bm.alignRow1, bm.sStart);
    seqan::setEndPosition(bm.alignRow1, bm.sEnd);
}

#ifdef SEQAN_SIMD_ENABLED

template <typename TDepSetH,
          typename TDepSetV,
          typename TBlastMatches>
inline void
_setupDepSets(TDepSetH & depSetH, TDepSetV & depSetV, TBlastMatches const & blastMatches)
{
    using TSimdAlign    = typename seqan::SimdVector<int16_t>::Type;
    unsigned constexpr sizeBatch = seqan::LENGTH<TSimdAlign>::VALUE;
    unsigned const      fullSize = sizeBatch * ((seqan::length(blastMatches) + sizeBatch - 1) / sizeBatch);

    seqan::clear(depSetH);
    seqan::clear(depSetV);
    seqan::reserve(depSetH, fullSize);
    seqan::reserve(depSetV, fullSize);

    for (auto const & bm : blastMatches)
    {
        seqan::appendValue(depSetH, seqan::source(bm.alignRow0));
        seqan::appendValue(depSetV, seqan::source(bm.alignRow1));
    }

    // fill up last batch
    for (size_t i = seqan::length(blastMatches); i < fullSize; ++i)
    {
        seqan::appendValue(depSetH, seqan::source(seqan::back(blastMatches).alignRow0));
        seqan::appendValue(depSetV, seqan::source(seqan::back(blastMatches).alignRow1));
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
                  std::integral_constant<bool, withTrace> const &,
                  bsDirection const dir = bsDirection::fwd)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TAlignConfig  = seqan::AlignConfig2<seqan::LocalAlignment_<>,
                                              seqan::DPBandConfig<seqan::BandOff>,
                                              seqan::FreeEndGaps_<seqan::True, seqan::True, seqan::True, seqan::True>,
                                              std::conditional_t<withTrace,
                                                                 seqan::TracebackOn<
                                                                     seqan::TracebackConfig_<seqan::CompleteTrace,
                                                                                             seqan::GapsLeft> >,
                                                                 seqan::TracebackOff> >;

    using TSimdAlign    = typename seqan::SimdVector<int16_t>::Type;
    using TSimdScore    = seqan::Score<TSimdAlign, seqan::ScoreSimdWrapper<typename TGlobalHolder::TScoreSchemeAlign> >;
    using TSize         = typename seqan::Size<typename TLocalHolder::TAlignRow0>::Type;
    using TMatch        = typename TGlobalHolder::TMatch;
    using TPos          = typename TMatch::TPos;
    using TTraceSegment = seqan::TraceSegment_<TPos, TSize>;

    unsigned constexpr sizeBatch = seqan::LENGTH<TSimdAlign>::VALUE;
    unsigned const      fullSize = sizeBatch * ((seqan::length(blastMatches) + sizeBatch - 1) / sizeBatch);

    auto const & currentScoringScheme = dir == bsDirection::fwd ? lH.gH.scoringSchemeAlign : lH.gH.scoringSchemeAlignBSRev;
    TSimdScore simdScoringScheme(currentScoringScheme);
    seqan::StringSet<seqan::String<TTraceSegment> > trace;

    // TODO when band is available, create inside block with band
    TAlignConfig config;//(0, 2*band)

    auto matchIt = blastMatches.begin();
    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        auto infSetH = infixWithLength(depSetH, pos, sizeBatch);
        auto infSetV = infixWithLength(depSetV, pos, sizeBatch);

        TSimdAlign resultsBatch;

        seqan::clear(trace);
        seqan::resize(trace, sizeBatch, seqan::Exact());

        // TODO pass in lH.dpSIMDContext
        seqan::_prepareAndRunSimdAlignment(resultsBatch,
                                           trace,
                                           infSetH,
                                           infSetV,
                                           simdScoringScheme,
                                           config,
                                           typename TLocalHolder::TScoreExtension());

        for (auto x = pos; x < pos + sizeBatch && x < seqan::length(blastMatches); ++x)
        {
            if constexpr (withTrace)
                seqan::_adaptTraceSegmentsTo(matchIt->alignRow0, matchIt->alignRow1, trace[x - pos]);
            else
                matchIt->alignStats.alignmentScore = resultsBatch[x - pos];

            ++matchIt;
        }
    }
}

template <typename TLocalHolder>
inline void
iterateMatchesFullSimd(TLocalHolder & lH, bsDirection const dir = bsDirection::fwd)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TBlastMatch   = typename TLocalHolder::TBlastMatch;

    auto const & const_gH = lH.gH;

    // statistics
#ifdef LAMBDA_MICRO_STATS
    ++lH.stats.numQueryWithExt;
    lH.stats.numExtScore += seqan::length(lH.matches);

    double start = sysTime();
#endif

    // Prepare string sets with sequences.
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow0>::Type> depSetH;
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow1>::Type> depSetV;

    // create blast matches
    std::list<TBlastMatch> blastMatches;
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // In BS-mode, skip those results that have wrong orientation
        if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
        {
            if ((dir == bsDirection::fwd && (it->subjId % 2)) || (dir == bsDirection::rev && !(it->subjId % 2)))
                continue;
        }
        // create blastmatch in list without copy or move
        blastMatches.emplace_back(lH.qryIds [it->qryId / TGlobalHolder::qryNumFrames],
                                  const_gH.indexFile.ids[it->subjId / TGlobalHolder::sbjNumFrames]);

        TBlastMatch & bm = blastMatches.back();

        bm._n_qId = it->qryId / TGlobalHolder::qryNumFrames;
        bm._n_sId = it->subjId / TGlobalHolder::sbjNumFrames;

        bm.qLength = //std::ranges::size(lH.transQrySeqs[it->qryId]);
                    std::ranges::size(lH.qrySeqs[bm._n_qId]);

        bm.sLength = // std::ranges::size(lH.gH.transSbjSeqs[it->subjId]);
                     std::ranges::size(lH.gH.indexFile.seqs[bm._n_sId]);

        _setupAlignInfix(bm, *it, lH);

        _setFrames(bm, *it, lH);

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.indexFile.sTaxIds[bm._n_sId];
    }
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtend      += sysTime() - start;
    lH.stats.timeExtendTrace += sysTime() - start; //TODO remove this line!

    // filter out duplicates
    start = sysTime();
#endif
    auto before = seqan::length(blastMatches);
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
    lH.stats.hitsDuplicate += (before - seqan::length(blastMatches));

    // sort by lengths to minimize padding in SIMD
    blastMatches.sort([] (auto const & l, auto const & r)
    {
        return std::make_tuple(seqan::length(seqan::source(l.alignRow0)), seqan::length(seqan::source(l.alignRow1))) <
               std::make_tuple(seqan::length(seqan::source(r.alignRow0)), seqan::length(seqan::source(r.alignRow1)));
    });
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeSort += sysTime() - start;

    start = sysTime();
#endif
    // fill batches
    _setupDepSets(depSetH, depSetV, blastMatches);

    // Run extensions WITHOUT ALIGNMENT
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::false_type(), dir);

    auto const & currentScoringScheme = dir == bsDirection::fwd ? lH.gH.scoringSchemeAlign : lH.gH.scoringSchemeAlignBSRev;

    // copmute evalues and filter based on evalue
    for (auto it = blastMatches.begin(), itEnd = blastMatches.end(); it != itEnd; /*below*/)
    {
        TBlastMatch & bm = *it;

        computeEValueThreadSafe(bm, bm.qLength, seqan::context(lH.gH.outfileBlastTab));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            it = blastMatches.erase(it);
            continue;
        }

        ++it;
    }
    if (seqan::length(blastMatches) == 0)
        return;

    // statistics
#ifdef LAMBDA_MICRO_STATS
    lH.stats.numExtAli += seqan::length(blastMatches);
    lH.stats.timeExtend += sysTime() - start;
    start = sysTime();
#endif

    // reset and fill batches
    _setupDepSets(depSetH, depSetV, blastMatches);

    // Run extensions WITH ALIGNMENT
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::true_type(), dir);

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

        seqan::computeAlignmentStats(bm.alignStats, bm.alignRow0, bm.alignRow1, currentScoringScheme);

        if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        {
            ++lH.stats.hitsFailedExtendPercentIdentTest;
            it = blastMatches.erase(it);
            continue;
        }

        seqan::computeBitScore(bm, seqan::context(lH.gH.outfileBlastTab));

        // evalue computed previously

        ++it;
    }
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif

    // move the blastMatches into localHolder's cache
    lH.blastMatches.splice(lH.blastMatches.end(), blastMatches);
}

#endif // SEQAN_SIMD_ENABLED

template <typename TLocalHolder>
inline void
iterateMatchesFullSerial(TLocalHolder & lH)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;

    size_t band = _bandSize(lH.transQrySeqs[lH.matches[0].qryId].size(), lH);

#ifdef LAMBDA_MICRO_STATS
    double start = sysTime();
#endif

    auto const & const_gH = lH.gH;

    // create blast matches
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // create blastmatch in list without copy or move
        lH.blastMatches.emplace_back(lH.qryIds [it->qryId / TGlobalHolder::qryNumFrames],
                                     const_gH.indexFile.ids[it->subjId / TGlobalHolder::sbjNumFrames]);

        auto & bm = lH.blastMatches.back();
        auto &  m = *it;

        bm._n_qId = it->qryId / TGlobalHolder::qryNumFrames;
        bm._n_sId = it->subjId / TGlobalHolder::sbjNumFrames;

        bm.qLength = //std::ranges::size(lH.transQrySeqs[it->qryId]);
                     std::ranges::size(lH.qrySeqs[bm._n_qId]);
        bm.sLength = //std::ranges::size(lH.gH.transSbjSeqs[it->subjId]);
                     std::ranges::size(lH.gH.indexFile.seqs[bm._n_sId]);

        _setupAlignInfix(bm, *it, lH);

        _setFrames(bm, m, lH);

        auto & currentScoringScheme = TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS && (it->subjId % 2) ?
                                      lH.gH.scoringSchemeAlignBSRev :  lH.gH.scoringSchemeAlign;

        // Run extension WITHOUT TRACEBACK
        bm.alignStats.alignmentScore = localAlignmentScore(bm.alignRow0,
                                                           bm.alignRow1,
                                                           currentScoringScheme,
                                                           -band,
                                                           +band);

        computeEValueThreadSafe(bm, bm.qLength, seqan::context(lH.gH.outfileBlastTab));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            lH.blastMatches.pop_back();
            continue;
        }

        // Run extension WITH TRACEBACK
        seqan::localAlignment(bm.alignRow0,
                              bm.alignRow1,
                              currentScoringScheme,
                              -band,
                              +band);

        _expandAlign(bm, lH);

        seqan::computeAlignmentStats(bm.alignStats, bm.alignRow0, bm.alignRow1, currentScoringScheme);

        if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        {
            ++lH.stats.hitsFailedExtendPercentIdentTest;
            lH.blastMatches.pop_back();
            continue;
        }

        seqan::computeBitScore(bm, seqan::context(lH.gH.outfileBlastTab));

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.indexFile.sTaxIds[bm._n_sId];

    }

#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif
}

template <typename TLocalHolder>
inline void
writeRecords(TLocalHolder & lH)
{
    using TBlastRecord  = typename TLocalHolder::TBlastRecord;

    // divide matches into records (per query) and write
    for (auto it = lH.blastMatches.begin(), itLast = lH.blastMatches.begin();
         seqan::length(lH.blastMatches) > 0;
         /*below*/)
    {
        if ((it == lH.blastMatches.end()) || ((it != lH.blastMatches.begin()) && (it->_n_qId != itLast->_n_qId)))
        {
            // create a record for each query
            TBlastRecord record(lH.qryIds[itLast->_n_qId]);
            record.qLength = seqan::length(lH.qrySeqs[itLast->_n_qId]);
            // move the matches into the record
            record.matches.splice(record.matches.begin(),
                                  lH.blastMatches,
                                  lH.blastMatches.begin(),
                                  it);
            // write to file
            _writeRecord(record, lH);

            it = lH.blastMatches.begin();
            itLast = lH.blastMatches.begin();
        } else
        {
            itLast = it;
            ++it;
        }
    }
}

template <typename TLocalHolder>
inline void
iterateMatches(TLocalHolder & lH)
{
#ifdef SEQAN_SIMD_ENABLED
    if (lH.options.extensionMode == LambdaOptions::ExtensionMode::FULL_SIMD)
    {
        iterateMatchesFullSimd(lH, bsDirection::fwd);
        if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
        {
            iterateMatchesFullSimd(lH, bsDirection::rev);
            lH.blastMatches.sort([] (auto const & lhs, auto const & rhs)
            {
                return lhs._n_qId < rhs._n_qId;
            });
        }
    }
    else
#endif

    iterateMatchesFullSerial(lH);

}
