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

#include <cstdint>
#include <iomanip>
#include <type_traits>

#include <seqan/align_extend.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/blast.h>
#include <seqan/misc/terminal.h>
#include <seqan/reduced_aminoacid.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/translation.h>

#include <bio/alphabet/aminoacid/aa27.hpp>
#include <bio/io/seq/reader.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/complement.hpp>
#include <bio/ranges/views/translate_join.hpp>
#if __cpp_lib_ranges <= 202106L
#    include <bio/ranges/views/persist.hpp>
#endif

#include <fmindex-collection/DenseCSA.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>

#include "bisulfite_scoring.hpp"
#include "evaluate_bisulfite_alignment.hpp"
#include "search_datastructures.hpp"
#include "search_misc.hpp"
#include "search_options.hpp"
#include "view_dna_n_to_random.hpp"
#include "view_reduce_to_bisulfite.hpp"

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function readIndexOptions()
// --------------------------------------------------------------------------

void readIndexOptions(LambdaOptions & options)
{
    std::string filename(options.indexFilePath);
    size_t      t = std::min<size_t>(options.threads, 4); // more than four threads is harmful

    /* verify the correct index_generation */
    uint64_t indexGeneration = -1;
    if (filename.ends_with(".lba") || filename.ends_with(".lba.gz"))
    {
        bio::io::transparent_istream is{filename, {.threads = t}};
        cereal::BinaryInputArchive   iarchive(is);
        iarchive(cereal::make_nvp("generation", indexGeneration));
    }
    else if (filename.ends_with(".lta") || filename.ends_with(".lta.gz"))
    {
        bio::io::transparent_istream is{filename, {.threads = t}};
        cereal::JSONInputArchive     iarchive(is);
        iarchive(cereal::make_nvp("generation", indexGeneration));
    }
    else
    {
        throw 59;
    }

    if (indexGeneration != supportedIndexGeneration)
    {
        std::string error = "ERROR: this version of Lambda only supports INDEXES of generation " +
                            std::to_string(supportedIndexGeneration) +
                            ", but the provided index was of generation: " + std::to_string(indexGeneration) +
                            ". PLEASE RECREATE THE INDEX!\n";
        throw std::runtime_error{error};
    }

    /* load index "options" */
    fake_index_file              f{options.indexFileOptions};
    bio::io::transparent_istream is{filename, {.threads = t}};

    if (filename.ends_with(".lba") || filename.ends_with(".lba.gz"))
    {
        cereal::BinaryInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", f));
    }
    else if (filename.ends_with(".lta") || filename.ends_with(".lta.gz"))
    {
        cereal::JSONInputArchive iarchive(is);
        iarchive(cereal::make_nvp("lambda index", f));
    }
    else
    {
        throw 59;
    }
}

// --------------------------------------------------------------------------
// Function checkRAM()
// --------------------------------------------------------------------------

void checkRAM(LambdaOptions const & options)
{
    myPrint(options, 1, "Checking memory requirements... ");
    uint64_t ram       = getTotalSystemMemory();
    uint64_t sizeIndex = 0;
    uint64_t sizeQuery = 0;

    sizeIndex = fileSize(options.indexFilePath.c_str());

    sizeQuery = fileSize(options.queryFile.c_str());

    uint64_t requiredRAM = 0;

    //TODO I have verified this for 2-16 threads on uniprot, but we should also double-check on Uniref50
    if (options.lazyQryFile)
        requiredRAM = (sizeIndex * 12) / 10; // give it +20%
    else
        requiredRAM = ((sizeIndex + sizeQuery) * 12) / 10; // give it +20%

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

template <DbIndexType  c_dbIndexType,
          AlphabetEnum c_origSbjAlph,
          AlphabetEnum c_transAlph,
          AlphabetEnum c_redAlph,
          AlphabetEnum c_origQryAlph>
void prepareScoring(
  GlobalDataHolder<c_dbIndexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
  LambdaOptions const &                                                                   options)
{
    if constexpr (c_transAlph != AlphabetEnum::AMINO_ACID)
    {
        // Seqan2
        seqan::setScoreMatch(context(globalHolder.outfileBlastTab).scoringScheme, options.match);
        seqan::setScoreMismatch(context(globalHolder.outfileBlastTab).scoringScheme, options.misMatch);

        if constexpr (c_redAlph == AlphabetEnum::DNA3BS)
        {
            // Seqan2
            seqan::setScoreBisulfiteMatrix(globalHolder.scoringSchemeAlign,
                                           options.match,
                                           options.misMatch,
                                           bsDirection::fwd);
            seqan::setScoreBisulfiteMatrix(globalHolder.scoringSchemeAlignBSRev,
                                           options.match,
                                           options.misMatch,
                                           bsDirection::rev);
        }
        else
        {
            // Seqan2
            globalHolder.scoringSchemeAlign = seqan::seqanScheme(context(globalHolder.outfileBlastTab).scoringScheme);
            globalHolder.scoringSchemeAlignBSRev =
              seqan::seqanScheme(context(globalHolder.outfileBlastTab).scoringScheme);
        }
    }
    else
    {
        seqan::AminoAcidScoreMatrixID seqan2_matrix_id{};

        switch (options.scoringMethod)
        {
            case 45:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM45;
                break;
            case 62:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM62;
                break;
            case 80:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM80;
                break;
            default:
                break;
        }

        // seqan2
        seqan::setScoreMatrixById(seqan::context(globalHolder.outfileBlastTab).scoringScheme._internalScheme,
                                  seqan2_matrix_id);
        seqan::setScoreMatrixById(globalHolder.scoringSchemeAlign, seqan2_matrix_id);
        seqan::setScoreMatrixById(globalHolder.scoringSchemeAlignBSRev, seqan2_matrix_id);
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

template <DbIndexType  c_indexType,
          AlphabetEnum c_origSbjAlph,
          AlphabetEnum c_transAlph,
          AlphabetEnum c_redAlph,
          AlphabetEnum c_origQryAlph>
void loadDbIndexFromDisk(
  GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
  LambdaOptions const &                                                                 options)
{
    std::string strIdent = "Loading Database Index...";
    myPrint(options, 1, strIdent);
    double start = sysTime();

    {
        std::string filename(options.indexFilePath);
        size_t      threads = std::min<size_t>(options.threads, 4); // more than four threads is harmful
        bio::io::transparent_istream is{filename, {.threads = threads}};

        if (filename.ends_with(".lba") || filename.ends_with(".lba.gz"))
        {
            cereal::BinaryInputArchive iarchive(is);
            iarchive(cereal::make_nvp("lambda index", globalHolder.indexFile));
        }
        else if (filename.ends_with(".lta") || filename.ends_with(".lta.gz"))
        {
            cereal::JSONInputArchive iarchive(is);
            iarchive(cereal::make_nvp("lambda index", globalHolder.indexFile));
        }
        else
        {
            throw 88;
        }
    }

    globalHolder.transSbjSeqs = globalHolder.indexFile.seqs | sbjTransView<c_origSbjAlph, c_transAlph, c_redAlph>;
    globalHolder.redSbjSeqs   = globalHolder.transSbjSeqs | redView<c_transAlph, c_redAlph>;

    size_t searchSpaceSize = 0ull;

    if (options.verbosity == 2)
    {
        searchSpaceSize = bio::meta::overloaded{[]<typename T>(bio::ranges::concatenated_sequences<T> const & seqs)
                                                { return seqs.concat_size(); },
                                                [](auto const & seqs)
                                                {
                                                    auto v = seqs | std::views::transform(std::ranges::size);
                                                    return std::reduce(v.begin(), v.end(), 0ull);
                                                }}(globalHolder.transSbjSeqs);
    }

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    myPrint(options, 2, "    # original subjects:   ", globalHolder.indexFile.seqs.size(), "\n");
    myPrint(options, 2, "    # translated subjects: ", globalHolder.transSbjSeqs.size(), "\n");
    myPrint(options, 2, "    # reduced subjects:    ", globalHolder.redSbjSeqs.size(), "\n");
    myPrint(options, 2, "    size of search space:  ", searchSpaceSize, "\n");
    bool const indexHasSTaxIDs = globalHolder.indexFile.sTaxIds.size() == globalHolder.indexFile.seqs.size();
    myPrint(options, 2, "    has taxonomic IDs:     ", indexHasSTaxIDs, "\n");
    bool const indexHasTaxTree = globalHolder.indexFile.taxonNames.size() >= globalHolder.indexFile.seqs.size();
    myPrint(options, 2, "    has taxonomic tree:    ", indexHasTaxTree, "\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    if (options.hasSTaxIds && !indexHasSTaxIDs)
    {
        throw std::runtime_error{
          "You requested printing of taxonomic IDs and/or taxonomic binning, but the index "
          "does not contain taxonomic information. Recreate it and provide --acc-tax-map ."};
    }
    if (options.computeLCA && !indexHasTaxTree)
    {
        throw std::runtime_error{
          "You requested taxonomic binning, but the index "
          "does not contain a taxonomic tree. Recreate it and provide --tax-dump-dir ."};
    }

    /* this is actually part of prepareScoring(), but the values are just available now */
    auto sizes = globalHolder.redSbjSeqs | std::views::transform(std::ranges::size);
    seqan::context(globalHolder.outfileBlastTab).dbTotalLength  = std::accumulate(sizes.begin(), sizes.end(), 0ull);
    seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = globalHolder.redSbjSeqs.size();
    seqan::context(globalHolder.outfileBlastTab).dbName         = options.indexFilePath;
}

// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------

template <DbIndexType  c_indexType,
          AlphabetEnum c_origSbjAlph,
          AlphabetEnum c_transAlph,
          AlphabetEnum c_redAlph,
          AlphabetEnum c_origQryAlph>
void loadQuery(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
               LambdaOptions const &                                                                 options)
{
    using TGH = GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>;

    double start = sysTime();

    std::string strIdent = "Loading Query Sequences...";
    myPrint(options, 1, strIdent);

    bio::io::seq::record r{.id = std::string{}, .seq = std::vector<typename TGH::TOrigQryAlph>{}, .qual = std::ignore};
    bio::io::seq::reader reader{options.queryFile, bio::io::seq::reader_options{.record = r}};
    for (auto & rec : reader)
    {
        globalHolder.qryIds.push_back(std::move(rec.id));
        globalHolder.qrySeqs.push_back(std::move(rec.seq));
    }

    // parse the file completely and get count in one line:
    globalHolder.queryTotal = globalHolder.qrySeqs.size();

    // batch-size as set in options (unless too few sequences)
    globalHolder.records_per_batch = std::max<size_t>(
      std::min<size_t>(globalHolder.queryTotal / (options.threads * 10), options.maximumQueryBlockSize),
      1);
    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

template <DbIndexType  c_indexType,
          AlphabetEnum c_origSbjAlph,
          AlphabetEnum c_transAlph,
          AlphabetEnum c_redAlph,
          AlphabetEnum c_origQryAlph>
void countQuery(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
                LambdaOptions const &                                                                 options)
{
    double start = sysTime();

    std::string strIdent = "Counting Query Sequences...";
    myPrint(options, 1, strIdent);

#if __GNUC__ >= 11
    bio::io::seq::record r{.id = std::ignore, .seq = std::ignore, .qual = std::ignore};
    bio::io::seq::reader reader{options.queryFile, bio::io::seq::reader_options{.record = r}};
#else // it is absolutely unclear why this workaround is needed here––CTAD with member-init works everywhere else
    bio::io::seq::record<bio::meta::ignore_t, bio::meta::ignore_t, bio::meta::ignore_t> r{};
    bio::io::seq::reader reader{options.queryFile, bio::io::seq::reader_options<decltype(r)>{}};
#endif

    // parse the file completely and get count in one line:
    globalHolder.queryTotal = std::ranges::distance(reader);

    // batch-size as set in options (unless too few sequences)
    globalHolder.records_per_batch = std::max<size_t>(
      std::min<size_t>(globalHolder.queryTotal / (options.threads * 10), options.maximumQueryBlockSize),
      1);
    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

template <DbIndexType  c_indexType,
          AlphabetEnum c_origSbjAlph,
          AlphabetEnum c_transAlph,
          AlphabetEnum c_redAlph,
          AlphabetEnum c_origQryAlph>
auto createQryView(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph> & globalHolder,
                   LambdaOptions const &                                                                 options)
{
    using TGH = GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>;
    bio::io::seq::record r{.id = std::string{}, .seq = std::vector<typename TGH::TOrigQryAlph>{}, .qual = std::ignore};
    bio::io::seq::reader reader{options.queryFile, bio::io::seq::reader_options{.record = r}};

#if __cpp_lib_ranges <= 202106L
    return std::move(reader) | bio::views::persist |
           views::async_input_buffer(globalHolder.records_per_batch * options.threads);
#else
    return std::move(reader) | views::async_input_buffer(globalHolder.records_per_batch * options.threads);
#endif
}

/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function seedLooksPromising()
// --------------------------------------------------------------------------

// perform a fast local alignment score calculation on the seed and see if we
// reach above threshold
// WARNING the following function only works for hammingdistanced seeds
template <typename TGlobalHolder>
inline bool seedLooksPromising(LocalDataHolder<TGlobalHolder> const & lH, typename TGlobalHolder::TMatch const & m)
{
    int64_t  effectiveQBegin = m.qryStart;
    int64_t  effectiveSBegin = m.subjStart;
    uint64_t actualLength    = m.qryEnd - m.qryStart;
    uint64_t effectiveLength =
      std::max(static_cast<uint64_t>(lH.searchOpts.seedLength * lH.options.preScoring), actualLength);

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

        effectiveLength =
          std::min({static_cast<uint64_t>(std::ranges::size(lH.transQrySeqs[m.qryId]) - effectiveQBegin),
                    static_cast<uint64_t>(std::ranges::size(lH.gH.transSbjSeqs[m.subjId]) - effectiveSBegin),
                    effectiveLength});
    }

    auto const & qSeq = lH.transQrySeqs[m.qryId] |
                        bio::views::slice(effectiveQBegin, static_cast<int64_t>(effectiveQBegin + effectiveLength));
    auto const & sSeq = lH.gH.transSbjSeqs[m.subjId] |
                        bio::views::slice(effectiveSBegin, static_cast<int64_t>(effectiveSBegin + effectiveLength));

    int       s        = 0;
    int       maxScore = 0;
    int const thresh   = lH.options.preScoringThresh * effectiveLength;

    // score the diagonal
    auto & currentScoringScheme = TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS && m.subjId % 2
                                    ? lH.gH.scoringSchemeAlignBSRev
                                    : lH.gH.scoringSchemeAlign;

    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        s += score(currentScoringScheme, qSeq[i], sSeq[i]);

        if (s < 0)
            s = 0;
        else if (s > maxScore)
            maxScore = s;

        if (maxScore >= thresh)
            return true;
    }
    return false;
}

template <typename TGlobalHolder, typename TSeed>
inline void search_impl(LocalDataHolder<TGlobalHolder> & lH, TSeed && seed)
{
    if constexpr (TGlobalHolder::c_dbIndexType == DbIndexType::FM_INDEX)
    {
        fmindex_collection::search_backtracking_with_buffers::search(
          lH.gH.indexFile.index,
          seed | bio::views::to_rank | fmindex_collection::add_sentinel,
          lH.searchOpts.maxSeedDist,
          lH.cursor_tmp_buffer,
          lH.cursor_tmp_buffer2,
          [&](auto cursor, size_t /*errors*/) { lH.cursor_buffer.push_back(cursor); });
    }
    else if constexpr (TGlobalHolder::c_dbIndexType == DbIndexType::BI_FM_INDEX)
    {
        if (lH.searchOpts.maxSeedDist == 0)
        {
            [&]()
            {
                using cursor_t = TGlobalHolder::TIndexCursor;

                auto query = seed | bio::views::to_rank | fmindex_collection::add_sentinel;

                auto cur = cursor_t{lH.gH.indexFile.index};
                for (size_t i{0}; i < query.size(); ++i)
                {
                    auto r = query[query.size() - i - 1];
                    cur    = cur.extendLeft(r);
                    if (cur.empty())
                    {
                        return;
                    }
                }
                lH.cursor_buffer.push_back(cur);
            }();
        }
        else if (lH.searchOpts.maxSeedDist == 1)
        {
            fmindex_collection::search_one_error::search(lH.gH.indexFile.index,
                                                         seed | bio::views::to_rank | fmindex_collection::add_sentinel,
                                                         [&](auto cursor, size_t /*errors*/)
                                                         { lH.cursor_buffer.push_back(cursor); });
        }
        else
        {
            fmindex_collection::search_pseudo::search</*editdistance=*/false>(
              lH.gH.indexFile.index,
              seed | bio::views::to_rank | fmindex_collection::add_sentinel,
              lH.searchScheme,
              [&](auto cursor, size_t /*errors*/) { lH.cursor_buffer.push_back(cursor); });
        }
    }
}

template <typename TGlobalHolder, typename TSeed>
inline void searchHalfExactImpl(LocalDataHolder<TGlobalHolder> & lH, TSeed && seed)
{
    auto extendRight = [&](auto & cursor, auto c)
    {
        auto newCursor = cursor.extendRight(c.to_rank() + 1);
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
    size_t const seedFirstHalfLength  = lH.searchOpts.seedLength / 2;
    size_t const seedSecondHalfLength = lH.searchOpts.seedLength - seedFirstHalfLength;

    lH.cursor_tmp_buffer.emplace_back(lH.gH.indexFile.index, 0);
    auto & c = lH.cursor_tmp_buffer.back().first;

    // extend by half exactly
    for (size_t i = 0; i < seedFirstHalfLength; ++i)
    {
        auto newCursor = c.extendRight(seed[i].to_rank() + 1);
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

        for (auto & [cursor, error_count] : lH.cursor_tmp_buffer)
        {
            if (error_count < lH.searchOpts.maxSeedDist)
            {
                for (size_t r = 0; r < bio::alphabet::size<alph_t>; ++r)
                {
                    alph_t cur_letter = bio::alphabet::assign_rank_to(r, alph_t{});

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
inline void search(LocalDataHolder<TGlobalHolder> & lH)
{
    using TTransAlph = typename TGlobalHolder::TTransAlph;
    using TMatch     = typename TGlobalHolder::TMatch;

    /* The flag changes seeding to take into account how successful previous seeding operations were.
     * It basically changes the "heuristicFactor" below to be evidence-based (although the formula is also
     * slightly different).
     * This makes seeding more robust to different kinds of datasets and options, but due to multi-threading
     * it also makes the results becoming non-deterministic (i.e. repeated runs of the program may produce
     * slightly different amounts of results).
     */
#ifdef LAMBDA_NONDETERMINISTIC_SEEDS
    [[maybe_unused]] size_t const hitPerFinalHit =
      (double)lH.stats.hitsAfterSeeding / std::max<double>(lH.stats.hitsFinal, 1.0);
#endif

    size_t           hitsThisSeq     = 0;  // hits per untranslated sequence (multiple frames counted together)
    size_t           needlesSum      = 0;  // cumulative size of all frames of a sequence
    size_t           needlesPos      = 0;  // current position in the cumulative length
    constexpr size_t heuristicFactor = 10; // a vague approximation of the success rate of seeds

    /* These values are used to track how many hits we already have per untranslated sequence,
     * how many seeds we still have after the current one, and based on that, how many hits we
     * expect that the current seed needs to have.
     * With this information, the current seed can be elongated to get better (but fewer) hits.
     */

    for (size_t i = 0; i < std::ranges::size(lH.redQrySeqs); ++i)
    {
        if (lH.redQrySeqs[i].size() < lH.searchOpts.seedLength)
            continue;

        if (i % TGlobalHolder::qryNumFrames == 0) // reset on every "real" new read
        {
            hitsThisSeq = 0;
            needlesSum  = 0;
            needlesPos  = 0;
            for (size_t j = 0; j < TGlobalHolder::qryNumFrames; ++j)
                needlesSum += lH.redQrySeqs[i + j].size();
        }

        for (size_t seedBegin = 0; /* below */; seedBegin += lH.searchOpts.seedOffset)
        {
            // skip proteine 'X' or Dna 'N', skip letter if next letter is the same
            while ((seedBegin < (lH.redQrySeqs[i].size() - lH.searchOpts.seedLength)) &&
                   ((lH.transQrySeqs[i][seedBegin] ==
                     bio::alphabet::assign_char_to('`', TTransAlph{})) || // assume that '°' gets converted to UNKNOWN
                    (lH.transQrySeqs[i][seedBegin] == lH.transQrySeqs[i][seedBegin + 1])))
                ++seedBegin;

            // termination criterium
            if (seedBegin > (lH.redQrySeqs[i].size() - lH.searchOpts.seedLength))
                break;

            // results are in cursor_buffer
            lH.cursor_buffer.clear();
            if (lH.options.seedHalfExact && lH.searchOpts.maxSeedDist != 0)
                searchHalfExactImpl(lH,
                                    lH.redQrySeqs[i] |
                                      bio::views::slice(seedBegin, seedBegin + lH.searchOpts.seedLength));
            else
                search_impl(lH, lH.redQrySeqs[i] | bio::views::slice(seedBegin, seedBegin + lH.searchOpts.seedLength));

            if (lH.options.adaptiveSeeding)
                lH.offset_modifier_buffer.clear();

            for (auto & cursor : lH.cursor_buffer)
            {
                size_t seedLength = lH.searchOpts.seedLength;

                // elongate seeds
                if (lH.options.adaptiveSeeding)
                {
#ifdef LAMBDA_NONDETERMINISTIC_SEEDS
                    size_t desiredOccs =
                      hitsThisSeq >= lH.options.maxMatches * hitPerFinalHit
                        ? 1
                        : (lH.options.maxMatches * hitPerFinalHit - hitsThisSeq) /
                            std::max<size_t>((needlesSum - needlesPos - seedBegin) / lH.searchOpts.seedOffset, 1ul);

#else
                    // lambda2 mode BUT NOT QUIET
                    // desiredOccs == the number of seed hits we estimate that we need to reach lH.options.maxMatches
                    // hitsThisSeq >= lH.options.maxMatches → if we have more than we need already, only look for one
                    // (lH.options.maxMatches - hitsThisSeq) * heuristicFactor → total desired hits
                    // ((needlesSum - needlesPos - seedBegin) / lH.searchOpts.seedOffset → number of remaining seeds
                    // dividing the last two yields desired hits FOR THE CURRENT SEED
                    size_t desiredOccs =
                      hitsThisSeq >= lH.options.maxMatches
                        ? 1
                        : (lH.options.maxMatches - hitsThisSeq) * heuristicFactor /
                            std::max<size_t>((needlesSum - needlesPos - seedBegin) / lH.searchOpts.seedOffset, 1ul);
#endif

                    if (desiredOccs == 0)
                        desiredOccs = 1;

                    // This aborts when we fall under the threshold
                    auto   old_cursor = cursor;
                    size_t old_count  = cursor.count();
                    while (seedBegin + seedLength < lH.redQrySeqs[i].size())
                    {
                        cursor = cursor.extendRight((lH.redQrySeqs[i][seedBegin + seedLength]).to_rank() + 1);

                        size_t new_count = cursor.count();

                        if (new_count < desiredOccs &&
                            new_count < old_count) // we always continue to extend if we don't loose anything
                        {
                            // revert last extension
                            cursor = old_cursor;
                            break;
                        }

                        ++seedLength;
                        old_count  = new_count;
                        old_cursor = cursor;
                    }
                }

                // discard over-abundant seeds (typically seeds at end of sequence that couldn't be elongated)
                if (cursor.count() > heuristicFactor * lH.options.maxMatches)
                    continue;

                // locate hits
                for (auto [subjNo, subjOffset] : fmindex_collection::LocateLinear{lH.gH.indexFile.index, cursor})
                {
                    TMatch m{static_cast<typename TMatch::TQId>(i),
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

template <typename TBlastMatch, typename TLocalHolder>
inline void _setFrames(TBlastMatch & bm, typename TLocalHolder::TMatch const & m, TLocalHolder const &)
{
    if constexpr (seqan::qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = (m.qryId % 3) + 1;
        if (m.qryId % 6 > 2)
            bm.qFrameShift = -bm.qFrameShift;
    }
    else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        bm.qFrameShift = (m.qryId % 2) + 1;
        if (m.qryId % 4 > 1)
            bm.qFrameShift = -bm.qFrameShift;
    }
    else if constexpr (seqan::qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = 1;
        if (m.qryId % 2)
            bm.qFrameShift = -bm.qFrameShift;
    }
    else
    {
        bm.qFrameShift = 0;
    }

    if constexpr (seqan::sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = (m.subjId % 3) + 1;
        if (m.subjId % 6 > 2)
            bm.sFrameShift = -bm.sFrameShift;
    }
    else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        bm.sFrameShift = (m.subjId % 2) + 1;
    }
    else if constexpr (seqan::sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = 1;
        if (m.subjId % 2)
            bm.sFrameShift = -bm.sFrameShift;
    }
    else
    {
        bm.sFrameShift = 0;
    }
}

// --------------------------------------------------------------------------
// Function _writeMatches()
// --------------------------------------------------------------------------

template <typename TBlastRecord, typename TLocalHolder>
inline void _writeRecord(TBlastRecord & record, TLocalHolder & lH)
{
    auto const & const_gH = lH.gH;

    if (record.matches.size() > 0)
    {
        ++lH.stats.qrysWithHit;
        // sort and remove duplicates -> STL, yeah!
        auto const before = record.matches.size();

        // sort matches, using an inverted bitScore to have the highest score first
        record.matches.sort(
          [](auto const & m1, auto const & m2)
          {
              // bitscores explicitly switched so larger scores are sorted first
              // clang-format off
            return std::tie(m1._n_sId,
                            m1.qStart,
                            m1.qEnd,
                            m1.sStart,
                            m1.sEnd,
                            m1.qFrameShift,
                            m1.sFrameShift,
                            m2.bitScore) <
                    std::tie(m2._n_sId,
                            m2.qStart,
                            m2.qEnd,
                            m2.sStart,
                            m2.sEnd,
                            m2.qFrameShift,
                            m2.sFrameShift,
                            m1.bitScore);
              // clang-format on
          });

        // removes duplicates and keeping the ones with the greatest score
        record.matches.unique(
          [](auto const & m1, auto const & m2)
          {
              return std::tie(m1._n_sId, m1.qStart, m1.qEnd, m1.sStart, m1.sEnd, m1.qFrameShift, m1.sFrameShift) ==
                     std::tie(m2._n_sId, m2.qStart, m2.qEnd, m2.sStart, m2.sEnd, m2.qFrameShift, m2.sFrameShift);
          });
        lH.stats.hitsDuplicate2 += before - record.matches.size();

        // sort by evalue before writing
        record.matches.sort([](auto const & m1, auto const & m2) { return m1.bitScore > m2.bitScore; });

        // cutoff abundant
        if (record.matches.size() > lH.options.maxMatches)
        {
            lH.stats.hitsAbundant += record.matches.size() - lH.options.maxMatches;
            record.matches.resize(lH.options.maxMatches);
        }
        lH.stats.hitsFinal += record.matches.size();

        /* count uniq qry-subj-pairs */
        lH.uniqSubjIds.clear();
        lH.uniqSubjIds.reserve(record.matches.size());
        for (auto const & bm : record.matches)
            lH.uniqSubjIds.insert(bm._n_sId);

        lH.stats.pairs += lH.uniqSubjIds.size();

        // compute LCA
        if (lH.options.computeLCA)
        {
            record.lcaTaxId = 0;
            for (auto const & bm : record.matches)
            {
                if ((lH.gH.indexFile.sTaxIds[bm._n_sId].size() > 0) &&
                    (const_gH.indexFile.taxonParentIDs[lH.gH.indexFile.sTaxIds[bm._n_sId][0]] != 0))
                {
                    record.lcaTaxId = lH.gH.indexFile.sTaxIds[bm._n_sId][0];
                    break;
                }
            }

            if (record.lcaTaxId != 0)
                for (auto const & bm : record.matches)
                    for (uint32_t const sTaxId : lH.gH.indexFile.sTaxIds[bm._n_sId])
                        // Unassigned subjects are simply ignored
                        if (const_gH.indexFile.taxonParentIDs[sTaxId] != 0)
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

template <typename TLocalHolder>
inline void _widenMatch(Match & m, TLocalHolder const & lH)
{
    // move sStart as far left as needed to cover the part of query before qryStart
    m.subjStart = (m.subjStart < m.qryStart) ? 0 : m.subjStart - m.qryStart;

    /* always align full query independent of hit-region */
    m.qryStart = 0;
    m.qryEnd   = lH.transQrySeqs[m.qryId].size();

    // there is no band in computation but this value extends begin and end of Subj to account for gaps
    uint64_t band = _bandSize(lH.transQrySeqs[m.qryId].size());

    // end on subject is beginning plus full query length plus band
    m.subjEnd =
      std::min<size_t>(m.subjStart + lH.transQrySeqs[m.qryId].size() + band, lH.gH.transSbjSeqs[m.subjId].size());

    // account for band in subj start
    m.subjStart = (band < m.subjStart) ? m.subjStart - band : 0;
}

template <typename TBlastMatch, typename TLocalHolder>
inline auto _untrueQryId(TBlastMatch const & bm, TLocalHolder const &)
{
    using namespace seqan;

    if constexpr (seqan::qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 6 + bm.qFrameShift - 1;
        else
            return bm._n_qId * 6 - bm.qFrameShift + 2;
    }
    else if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 4;
        else
            return bm._n_qId * 4 + 2;
    }
    else if constexpr (seqan::qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.qFrameShift > 0)
            return bm._n_qId * 2;
        else
            return bm._n_qId * 2 + 1;
    }
    else
    {
        return bm._n_qId;
    }
}

template <typename TBlastMatch, typename TLocalHolder>
inline auto _untrueSubjId(TBlastMatch const & bm, TLocalHolder const &)
{
    using namespace seqan;

    if constexpr (seqan::sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        if (bm.sFrameShift > 0)
            return bm._n_sId * 6 + bm.sFrameShift - 1;
        else
            return bm._n_sId * 6 - bm.sFrameShift + 2;
    }
    else if constexpr (seqan::sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram) ||
                       TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        if (bm.sFrameShift > 0)
            return bm._n_sId * 2;
        else
            return bm._n_sId * 2 + 1;
    }
    else
    {
        return bm._n_sId;
    }
}

template <typename TBlastMatch, typename TLocalHolder>
inline void _expandAlign(TBlastMatch & bm, TLocalHolder const & lH)
{
    using namespace seqan;

    auto oldQLen = seqan::length(source(bm.alignRow0));
    auto oldSLen = seqan::length(source(bm.alignRow1));

    // replace source from underneath without triggering reset
    bm.alignRow0._source = lH.transQrySeqs[_untrueQryId(bm, lH)] |
                           bio::views::slice(0, std::ranges::size(lH.transQrySeqs[_untrueQryId(bm, lH)]));
    bm.alignRow1._source = lH.gH.transSbjSeqs[_untrueSubjId(bm, lH)] |
                           bio::views::slice(0, std::ranges::size(lH.gH.transSbjSeqs[_untrueSubjId(bm, lH)]));

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

template <typename TDepSetH, typename TDepSetV, typename TBlastMatches>
inline void _setupDepSets(TDepSetH & depSetH, TDepSetV & depSetV, TBlastMatches const & blastMatches)
{
    using TSimdAlign             = typename seqan::SimdVector<int16_t>::Type;
    constexpr unsigned sizeBatch = seqan::LENGTH<TSimdAlign>::VALUE;
    unsigned const     fullSize  = sizeBatch * ((seqan::length(blastMatches) + sizeBatch - 1) / sizeBatch);

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

template <typename TDepSetH, typename TDepSetV, typename TBlastMatches, typename TLocalHolder, bool withTrace>
inline void _performAlignment(TDepSetH &      depSetH,
                              TDepSetV &      depSetV,
                              TBlastMatches & blastMatches,
                              TLocalHolder &  lH,
                              std::integral_constant<bool, withTrace> const &,
                              bsDirection const dir = bsDirection::fwd)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TAlignConfig  = seqan::AlignConfig2<
      seqan::LocalAlignment_<>,
      seqan::DPBandConfig<seqan::BandOff>,
      seqan::FreeEndGaps_<seqan::True, seqan::True, seqan::True, seqan::True>,
      std::conditional_t<withTrace,
                         seqan::TracebackOn<seqan::TracebackConfig_<seqan::CompleteTrace, seqan::GapsLeft>>,
                         seqan::TracebackOff>>;

    using TSimdAlign    = typename seqan::SimdVector<int16_t>::Type;
    using TSimdScore    = seqan::Score<TSimdAlign, seqan::ScoreSimdWrapper<typename TGlobalHolder::TScoreSchemeAlign>>;
    using TSize         = typename seqan::Size<typename TLocalHolder::TAlignRow0>::Type;
    using TMatch        = typename TGlobalHolder::TMatch;
    using TPos          = typename TMatch::TPos;
    using TTraceSegment = seqan::TraceSegment_<TPos, TSize>;

    constexpr unsigned sizeBatch = seqan::LENGTH<TSimdAlign>::VALUE;
    unsigned const     fullSize  = sizeBatch * ((seqan::length(blastMatches) + sizeBatch - 1) / sizeBatch);

    auto const & currentScoringScheme =
      dir == bsDirection::fwd ? lH.gH.scoringSchemeAlign : lH.gH.scoringSchemeAlignBSRev;
    TSimdScore                                     simdScoringScheme(currentScoringScheme);
    seqan::StringSet<seqan::String<TTraceSegment>> trace;

    // Alignment is not banded, because SeqAn2 doesn't support it for SIMD-fied alignment
    TAlignConfig config; //(0, 2*band)

    auto matchIt = blastMatches.begin();
    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        auto infSetH = infixWithLength(depSetH, pos, sizeBatch);
        auto infSetV = infixWithLength(depSetV, pos, sizeBatch);

        TSimdAlign resultsBatch;

        seqan::clear(trace);
        seqan::resize(trace, sizeBatch, seqan::Exact());

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
inline void _widenAndPreprocessMatches(std::span<Match> & matches, TLocalHolder & lH)
{
    auto before = matches.size();

    for (Match & m : matches)
        _widenMatch<TLocalHolder>(m, lH);

    std::ranges::sort(matches);

    if (matches.size() > 1)
    {
        // pairwise merge from left to right
        for (auto it = matches.begin(); it < matches.end() - 1; ++it)
        {
            Match & l = *it;
            Match & r = *(it + 1);
            if ((std::tie(l.qryId, l.subjId) == std::tie(r.qryId, r.subjId)) && (l.subjEnd >= r.subjStart))
            {
                l.subjEnd   = r.subjEnd;
                r.subjStart = l.subjStart;
            }
        }

        // pairwise "swallow" from right to left
        for (auto it = matches.rbegin(); it < matches.rend() - 1; ++it)
        {
            Match & r = *it;
            Match & l = *(it + 1);
            if ((std::tie(r.qryId, r.subjId) == std::tie(l.qryId, l.subjId)) && (r.subjStart < l.subjEnd))
            {
                l = r;
            }
        }

        auto [new_end, old_end] = std::ranges::unique(matches);               // move non-uniq to the end
        matches                 = std::span<Match>{matches.begin(), new_end}; // "resize" of the span
        lH.stats.hitsDuplicate += (before - matches.size());
    }
}

template <typename TLocalHolder>
inline void iterateMatchesFullSimd(std::span<Match> lambdaMatches, TLocalHolder & lH, bsDirection const dir)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TBlastMatch   = typename TLocalHolder::TBlastMatch;

    auto const & const_gH = lH.gH;

    // statistics
#ifdef LAMBDA_MICRO_STATS
    ++lH.stats.numQueryWithExt;
    lH.stats.numExtScore += seqan::length(lambdaMatches);

    double start = sysTime();
#endif

    // Prepare string sets with sequences.
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow0>::Type> depSetH;
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow1>::Type> depSetV;

    // pre-sort and filter
    _widenAndPreprocessMatches(lambdaMatches, lH);

    // create blast matches from Lambda matches
    std::list<TBlastMatch> blastMatches;
    for (Match const & m : lambdaMatches)
    {
        // create blastmatch in list without copy or move
        blastMatches.emplace_back(lH.qryIds[m.qryId / TGlobalHolder::qryNumFrames],
                                  const_gH.indexFile.ids[m.subjId / TGlobalHolder::sbjNumFrames]);

        TBlastMatch & bm = blastMatches.back();

        bm._n_qId = m.qryId / TGlobalHolder::qryNumFrames;
        bm._n_sId = m.subjId / TGlobalHolder::sbjNumFrames;

        bm.qLength = std::ranges::size(lH.qrySeqs[bm._n_qId]);
        bm.sLength = std::ranges::size(lH.gH.indexFile.seqs[bm._n_sId]);

        bm.qStart = m.qryStart;
        bm.qEnd   = m.qryEnd;
        bm.sStart = m.subjStart;
        bm.sEnd   = m.subjEnd;
        seqan::assignSource(bm.alignRow0, lH.transQrySeqs[m.qryId] | bio::views::slice(bm.qStart, bm.qEnd));
        seqan::assignSource(bm.alignRow1, lH.gH.transSbjSeqs[m.subjId] | bio::views::slice(bm.sStart, bm.sEnd));

        _setFrames(bm, m, lH);

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.indexFile.sTaxIds[bm._n_sId];
    }

    // sort by lengths to minimize padding in SIMD
    blastMatches.sort(
      [](auto const & l, auto const & r)
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

    auto const & currentScoringScheme =
      dir == bsDirection::fwd ? lH.gH.scoringSchemeAlign : lH.gH.scoringSchemeAlignBSRev;

    // copmute evalues and filter based on evalue
    for (auto it = blastMatches.begin(), itEnd = blastMatches.end(); it != itEnd; /*below*/)
    {
        TBlastMatch & bm = *it;

        if (lH.options.minBitScore >= 0)
        {
            seqan::computeBitScore(bm, seqan::context(lH.gH.outfileBlastTab));

            if (bm.bitScore < lH.options.minBitScore)
            {
                ++lH.stats.hitsFailedExtendBitScoreTest;
                it = blastMatches.erase(it);
                continue;
            }
        }

        if (lH.options.maxEValue >= 0)
        {
            computeEValueThreadSafe(bm, bm.qLength, seqan::context(lH.gH.outfileBlastTab));

            if (bm.eValue > lH.options.maxEValue)
            {
                ++lH.stats.hitsFailedExtendEValueTest;
                it = blastMatches.erase(it);
                continue;
            }
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
    blastMatches.sort([](auto const & lhs, auto const & rhs) { return lhs._n_qId < rhs._n_qId; });

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

        // not computed previously
        if (lH.options.minBitScore < 0)
            seqan::computeBitScore(bm, seqan::context(lH.gH.outfileBlastTab));

        if (lH.options.maxEValue < 0)
            computeEValueThreadSafe(bm, bm.qLength, seqan::context(lH.gH.outfileBlastTab));

        ++it;
    }
#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif

    // move the blastMatches into localHolder's cache
    lH.blastMatches.splice(lH.blastMatches.end(), blastMatches);
}

template <typename TLocalHolder>
inline void writeRecords(TLocalHolder & lH)
{
    using TBlastRecord = typename TLocalHolder::TBlastRecord;

    // divide matches into records (per query) and write
    for (auto it = lH.blastMatches.begin(), itLast = lH.blastMatches.begin(); seqan::length(lH.blastMatches) > 0;
         /*below*/)
    {
        if ((it == lH.blastMatches.end()) || ((it != lH.blastMatches.begin()) && (it->_n_qId != itLast->_n_qId)))
        {
            // create a record for each query
            TBlastRecord record(lH.qryIds[itLast->_n_qId]);
            record.qLength = seqan::length(lH.qrySeqs[itLast->_n_qId]);
            // move the matches into the record
            record.matches.splice(record.matches.begin(), lH.blastMatches, lH.blastMatches.begin(), it);
            // write to file
            _writeRecord(record, lH);

            it     = lH.blastMatches.begin();
            itLast = lH.blastMatches.begin();
        }
        else
        {
            itLast = it;
            ++it;
        }
    }
}

template <typename TLocalHolder>
inline void iterateMatches(TLocalHolder & lH)
{
    if constexpr (TLocalHolder::TGlobalHolder::c_redAlph == AlphabetEnum::DNA3BS)
    {
        std::ranges::sort(lH.matches,
                          [](Match const & l, Match const & r) {
                              return std::tuple<bool, Match const &>{l.subjId % 2, l} <
                                     std::tuple<bool, Match const &>{r.subjId % 2, r};
                          });

        auto it = std::ranges::find_if(lH.matches, [](Match const & m) { return m.subjId % 2; });

        iterateMatchesFullSimd(std::span{lH.matches.begin(), it}, lH, bsDirection::fwd);
        iterateMatchesFullSimd(std::span{it, lH.matches.end()}, lH, bsDirection::rev);
        lH.blastMatches.sort([](auto const & lhs, auto const & rhs) { return lhs._n_qId < rhs._n_qId; });
    }
    else
    {
        iterateMatchesFullSimd(lH.matches, lH, bsDirection::fwd);
    }
}

//-----------------------------------------------------------------------
// iterativeSearch
//-----------------------------------------------------------------------

void iterativeSearchPre(auto & lH)
{
    switch (lH.iterativeSearch)
    {
        case IterativeSearchMode::PHASE1:
            lH.successfulQueries.clear();
            lH.successfulQueries.resize(lH.qryIds.size());
            for (auto const & bm : lH.blastMatches)
            {
                assert(bm._n_qId < lH.successfulQueries.size());
                lH.successfulQueries[bm._n_qId] = true;
            }
            break;
        default:
            break;
    }
}

void iterativeSearchPost(auto & lH)
{
    switch (lH.iterativeSearch)
    {
        case IterativeSearchMode::PHASE1:
            {
                /* remove those queries from set that are already successful */
                size_t successfulCount =
                  std::accumulate(lH.successfulQueries.begin(), lH.successfulQueries.end(), 0ull);

                if (successfulCount != lH.qryIds.size() && successfulCount != 0)
                {
                    size_t                index = 0;
                    [[maybe_unused]] auto subr1 = std::ranges::remove_if(
                      lH.qryIds,
                      [&](size_t pos) { return lH.successfulQueries[pos] == true; },
                      [&](auto &&) { return index++; }); // converts element to index
                    assert(subr1.size() == successfulCount);

                    index                       = 0;
                    [[maybe_unused]] auto subr2 = std::ranges::remove_if(
                      lH.qrySeqs,
                      [&](size_t pos) { return lH.successfulQueries[pos] == true; },
                      [&](auto &&) { return index++; });
                    assert(subr2.size() == successfulCount);
                }

                lH.qryIds  = lH.qryIds | std::views::take(lH.qryIds.size() - successfulCount);
                lH.qrySeqs = lH.qrySeqs | std::views::take(lH.qrySeqs.size() - successfulCount);

                /* only switch to PHASE2 if there are any left */
                if (!lH.qryIds.empty())
                {
                    lH.iterativeSearch = IterativeSearchMode::PHASE2;
                    lH.searchOpts      = lH.options.searchOpts; // default
                    lH.searchScheme    = lH.searchScheme1;

                    // clear some caches (since we will not reset all of the holder)
                    lH.matches.clear();
                    lH.blastMatches.clear();
                }
                break;
            }
        case IterativeSearchMode::PHASE2:
            lH.iterativeSearch = IterativeSearchMode::PHASE1;
            lH.searchOpts      = lH.options.searchOpts0; // alternative scores
            lH.searchScheme    = lH.searchScheme0;
            break;
        default:
            break;
    }
}
