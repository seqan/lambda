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
// lambda.hpp: contains the main progam pipeline
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

#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/translate_join.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

// #include "search_datastructures.hpp"

// ============================================================================
// Forwards
// ============================================================================

// template <typename TRedAlph_,
//           typename TIndexSpec_,
//           typename TFileFormat,
//           seqan::BlastProgram p,
//           BlastTabularSpec h>
// class GlobalDataHolder;
//
// template <typename TGlobalHolder_,
//           typename TScoreExtension>
// class LocalDataHolder;

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
// Function readIndexOption()
// --------------------------------------------------------------------------

// inline void
// readIndexOption(std::string            & optionString,
//                 std::string      const & optionIdentifier,
//                 LambdaOptions    const & options)
// {
//     std::ifstream f{(options.indexDir + "/option:" + optionIdentifier).c_str(),
//                     std::ios_base::in | std::ios_base::binary};
//     if (f.is_open())
//     {
//         auto fit = directionIterator(f, Input());
//         readLine(optionString, fit);
//         f.close();
//     }
//     else
//     {
//         throw IndexException("Expected option specifier:\n" + options.indexDir + "/option:" + optionIdentifier);
//     }
// }

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

        // Seqan3
        globalHolder.scoringScheme.set_simple_scheme(seqan3::match_score{options.match},
                                                     seqan3::mismatch_score{options.misMatch});
    }
    else
    {
        seqan::AminoAcidScoreMatrixID       seqan2_matrix_id{};
        seqan3::aminoacid_similarity_matrix seqan3_matrix_id{};

        switch (options.scoringMethod)
        {
            case 45:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM45;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::BLOSUM45;
                break;
            case 62:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM62;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::BLOSUM62;
                break;
            case 80:
                seqan2_matrix_id = seqan::AminoAcidScoreMatrixID::BLOSUM80;
                seqan3_matrix_id = seqan3::aminoacid_similarity_matrix::BLOSUM80;
                break;
            default:
                break;
        }

        // seqan2
        seqan::setScoreMatrixById(seqan::context(globalHolder.outfileBlastTab).scoringScheme._internalScheme,
                                  seqan2_matrix_id);
        // Seqan3
        globalHolder.scoringScheme.set_similarity_matrix(seqan3_matrix_id);
    }

    // seqan2
    seqan::setScoreGapOpenBlast(seqan::context(globalHolder.outfileBlastTab).scoringScheme, options.gapOpen);
    seqan::setScoreGapExtend(seqan::context(globalHolder.outfileBlastTab).scoringScheme, options.gapExtend);

    if (!seqan::isValid(seqan::context(globalHolder.outfileBlastTab).scoringScheme))
        throw std::runtime_error{"Could not computer Karlin-Altschul-Values for Scoring Scheme.\n"};

    // seqan3
    globalHolder.gapScheme.set_affine(seqan3::gap_score{options.gapExtend},
                                      seqan3::gap_open_score{options.gapOpen});
}

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <typename TGlobalHolder>
void
loadDbIndexFromDisk(TGlobalHolder       & globalHolder,
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

    if constexpr (!std::is_lvalue_reference_v<typename TGlobalHolder::TTransSbjSeqs>) // is view
    {
        // update view, now that underlying range is available
        globalHolder.transSbjSeqs = globalHolder.indexFile.seqs | seqan3::views::translate_join;
    }

    if constexpr (!std::is_lvalue_reference_v<typename TGlobalHolder::TRedSbjSeqs>) // is view
    {
        // update view, now that underlying range is available
        if constexpr (TGlobalHolder::blastProgram == seqan::BlastProgram::BLASTN)
        {
            globalHolder.redSbjSeqs = globalHolder.transSbjSeqs
                                    | seqan3::views::dna_n_to_random;
        }
        else
        {
            globalHolder.redSbjSeqs = globalHolder.transSbjSeqs
                                    | seqan3::views::deep{seqan3::views::convert<typename TGlobalHolder::TRedAlph>};
        }
    }

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    // this is actually part of prepareScoring(), but the values are just available now
    if constexpr (seqan::sIsTranslated(TGlobalHolder::blastProgram))
    {
        // last value has sum of lengths
        // seqan::context(globalHolder.outfileBlastTab).dbTotalLength  = globalHolder.indexFile.origSeqLengths.back();
        // seqan::context(globalHolder.outfileBlastTab).dbNumberOfSeqs = globalHolder.indexFile.origSeqLengths.size() - 1;
        seqan::context(globalHolder.outfileBlastTab).dbTotalLength  =
            std::accumulate(globalHolder.indexFile.seqs.begin(), globalHolder.indexFile.seqs.end(), 0,
                [](size_t sum, std::vector<typename TGlobalHolder::TOrigSbjAlph> const & a)
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
}

// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------

template <AlphabetEnum c_origQryAlph>
struct QueryFileTraits : std::conditional_t<c_origQryAlph == AlphabetEnum::DNA5,
                                            seqan3::sequence_file_input_default_traits_dna,
                                            seqan3::sequence_file_input_default_traits_aa>
{
    template <typename _sequence_container>
    using sequence_container_container      = TCDStringSet<_sequence_container>;

    template <typename _id_container>
    using id_container_container            = TCDStringSet<_id_container>;

    template <typename _quality_container>
    using quality_container_container       = TCDStringSet<_quality_container>;
};

template <DbIndexType   c_indexType,
          AlphabetEnum  c_origSbjAlph,
          AlphabetEnum  c_transAlph,
          AlphabetEnum  c_redAlph,
          AlphabetEnum  c_origQryAlph>
void
loadQuery(GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>       & globalHolder,
          LambdaOptions                                                                       const & options)
{
    using TGH = GlobalDataHolder<c_indexType, c_origSbjAlph, c_transAlph, c_redAlph, c_origQryAlph>;
    // using TOrigAlph = typename TGH::TOrigQryAlph;
    using TRedAlph  = typename TGH::TRedAlph;
    double start = sysTime();

    std::string strIdent = "Loading Query Sequences and Ids...";
    myPrint(options, 1, strIdent);

    // TCDStringSet<std::vector<TOrigAlph>> origSeqs;

    // load
    seqan3::sequence_file_input<QueryFileTraits<c_origQryAlph>, seqan3::fields<seqan3::field::id, seqan3::field::seq>>
        infile{options.queryFile};

// TODO change to this once https://github.com/seqan/seqan3/issues/1135 is fixed
//     globalHolder.qryIds                 = std::move(seqan3::get<0>(infile));
//     globalHolder.qrySeqs    = std::move(seqan3::get<1>(infile));

    for (auto & [ id, seq ] : infile)
    {
        globalHolder.qryIds.push_back(std::move(id));
        globalHolder.qrySeqs.push_back(std::move(seq));
    }


    if (globalHolder.qrySeqs.size() == 0)
    {
        throw QueryException{"Zero sequences submitted."};
    }

    if constexpr (TGH::blastProgram == seqan::BlastProgram::BLASTN)
    {
        //TODO view-ify this at some point. Better yet, move to index.
        // ATTENTION, this won't work well if qryseqs is concatenated_sequences

        decltype(globalHolder.qrySeqs) tmp;

        tmp.resize(std::ranges::size(globalHolder.qrySeqs) * 2);
        for (size_t i = 0; i < std::ranges::size(tmp); ++i)
        {
            if (i % 2 == 0)
                tmp[i] = std::move(globalHolder.qrySeqs[i / 2]);
            else
                tmp[i] = tmp[i-1] | std::views::reverse
                                  | seqan3::views::complement
                                  | seqan3::views::to<std::ranges::range_value_t<decltype(tmp)>>;

        }

        std::swap(globalHolder.qrySeqs, tmp);
    }

    if constexpr (c_origQryAlph != c_transAlph)
    {
        globalHolder.transQrySeqs = globalHolder.qrySeqs
                                  | seqan3::views::translate_join;
    }

    if constexpr (c_transAlph != c_redAlph)
    {
        if constexpr (c_transAlph != AlphabetEnum::AMINO_ACID)
        {
            globalHolder.redQrySeqs = globalHolder.transQrySeqs
                                    | seqan3::views::dna_n_to_random;
        }
        else
        {
            globalHolder.redQrySeqs = globalHolder.transQrySeqs
                                    | seqan3::views::deep{seqan3::views::convert<TRedAlph>};
        }
    }


    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    unsigned long maxLen = 0ul;
    for (auto const & s : globalHolder.transQrySeqs)
        if (s.size() > maxLen)
            maxLen = s.size();

    myPrint(options, 2, "Runtime: ", finish, "s \n",
            "Number of effective query sequences: ",
            globalHolder.transQrySeqs.size(), "\nLongest query sequence: ",
            maxLen, "\n\n");
}

/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function generateSeeds()
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
                            static_cast<uint64_t>(std::ranges::size(lH.gH.transQrySeqs[m.qryId]) - effectiveQBegin),
                            static_cast<uint64_t>(std::ranges::size(lH.gH.transSbjSeqs[m.subjId]) - effectiveSBegin),
                            effectiveLength});
    }

    auto const & qSeq = lH.gH.transQrySeqs[m.qryId]
                      | seqan3::views::slice(effectiveQBegin, static_cast<int64_t>(effectiveQBegin + effectiveLength));
    auto const & sSeq = lH.gH.transSbjSeqs[m.subjId]
                      | seqan3::views::slice(effectiveSBegin, static_cast<int64_t>(effectiveSBegin + effectiveLength));

    int             s = 0;
    int      maxScore = 0;
    int const thresh  = lH.options.preScoringThresh * effectiveLength;

    // score the diagonal
    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        s += lH.gH.scoringScheme.score(qSeq[i], sSeq[i]);
        if (s < 0)
            s = 0;
        else if (s > maxScore)
            maxScore = s;

        if (maxScore >= thresh)
            return true;
    }

    return false;
}

template <typename TGlobalHolder>
inline void
search(LocalDataHolder<TGlobalHolder> & lH)
{
    using TTransAlph = typename TGlobalHolder::TTransAlph;
    using TMatch = typename TGlobalHolder::TMatch;

    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{static_cast<uint8_t>(lH.options.maxSeedDist)}};

    for (size_t i = lH.indexBeginQry; i < lH.indexEndQry; ++i)
    {
        if (lH.gH.redQrySeqs[i].size() < lH.options.seedLength)
            continue;

        for (size_t seedBegin = 0; /* below */; seedBegin += lH.options.seedOffset)
        {
            // skip proteine 'X' or Dna 'N'
            while ((lH.gH.transQrySeqs[i][seedBegin] == seqan3::assign_char_to('`', TTransAlph{})) && // assume that '°' gets converted to UNKNOWN
                   (seedBegin <= (lH.gH.redQrySeqs[i].size() - lH.options.seedLength)))
                ++seedBegin;

            // termination criterium
            if (seedBegin > (lH.gH.redQrySeqs[i].size() - lH.options.seedLength))
                break;

            auto results = seqan3::search(lH.gH.redQrySeqs[i] | seqan3::views::slice(seedBegin, seedBegin + lH.options.seedLength),
                                          lH.gH.indexFile.index,
                                          cfg);

            for (auto [ subjNo, subjOffset ] : results)
            {
                TMatch m {static_cast<typename TMatch::TQId>(i),
                          static_cast<typename TMatch::TSId>(subjNo),
                          static_cast<typename TMatch::TPos>(seedBegin),
                          static_cast<typename TMatch::TPos>(seedBegin + lH.options.seedLength),
                          static_cast<typename TMatch::TPos>(subjOffset),
                          static_cast<typename TMatch::TPos>(subjOffset + lH.options.seedLength)};

                ++lH.stats.hitsAfterSeeding;

                if (!seedLooksPromising(lH, m))
                    ++lH.stats.hitsFailedPreExtendTest;
                else
                    lH.matches.emplace_back(m);
            }
        }
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
    if (record.matches.size() > 0)
    {
        ++lH.stats.qrysWithHit;
        // sort and remove duplicates -> STL, yeah!
        auto const before = record.matches.size();

        if (!lH.options.filterPutativeDuplicates)
        {
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
        }

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
                if ((lH.gH.indexFile.sTaxIds[bm._n_sId].size() > 0) && (lH.gH.indexFile.taxonParentIDs[lH.gH.indexFile.sTaxIds[bm._n_sId][0]] != 0))
                {
                    record.lcaTaxId = lH.gH.indexFile.sTaxIds[bm._n_sId][0];
                    break;
                }
            }

            if (record.lcaTaxId != 0)
                for (auto const & bm : record.matches)
                    for (uint32_t const sTaxId : lH.gH.indexFile.sTaxIds[bm._n_sId])
                        if (lH.gH.indexFile.taxonParentIDs[sTaxId] != 0) // TODO do we want to skip unassigned subjects
                            record.lcaTaxId = computeLCA(lH.gH.indexFile.taxonParentIDs,
                                                         lH.gH.indexFile.taxonHeights,
                                                         sTaxId,
                                                         record.lcaTaxId);

            record.lcaId = lH.gH.indexFile.taxonNames[record.lcaTaxId];
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

    bm.qEnd = lH.gH.transQrySeqs[m.qryId].size();
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

    seqan::assignSource(bm.alignRow0, lH.gH.transQrySeqs[m.qryId] | seqan3::views::slice(bm.qStart, bm.qEnd));
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
    } else if constexpr (seqan::sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
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

    // replace source from underneath without triggereng reset
    bm.alignRow0._source = lH.gH.transQrySeqs[_untrueQryId(bm, lH)] | seqan3::views::slice(0, std::ranges::size(lH.gH.transQrySeqs[_untrueQryId(bm, lH)]));
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
                  std::integral_constant<bool, withTrace> const &)
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
    using TSimdScore    = seqan::Score<TSimdAlign, seqan::ScoreSimdWrapper<typename TGlobalHolder::TScoreScheme> >;
    using TSize         = typename seqan::Size<typename TLocalHolder::TAlignRow0>::Type;
    using TMatch        = typename TGlobalHolder::TMatch;
    using TPos          = typename TMatch::TPos;
    using TTraceSegment = seqan::TraceSegment_<TPos, TSize>;

    unsigned constexpr sizeBatch = seqan::LENGTH<TSimdAlign>::VALUE;
    unsigned const      fullSize = sizeBatch * ((seqan::length(blastMatches) + sizeBatch - 1) / sizeBatch);

    TSimdScore simdScoringScheme(seqan::seqanScheme(seqan::context(lH.gH.outfileBlastTab).scoringScheme));
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

        for(auto x = pos; x < pos + sizeBatch && x < seqan::length(blastMatches); ++x)
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
inline int
iterateMatchesFullSimd(TLocalHolder & lH)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TBlastPos     = uint32_t;
    using TBlastMatch   = seqan::BlastMatch<
                           typename TLocalHolder::TAlignRow0,
                           typename TLocalHolder::TAlignRow1,
                           TBlastPos,
                           std::string_view,
                           std::string_view,
                           std::vector<std::string>, // not used
                           std::span<uint32_t>>;

    using TBlastRecord  = seqan::BlastRecord<TBlastMatch,
                                             std::string_view,
                                             std::vector<std::string>, // not used
                                             std::string_view,
                                             uint32_t>;
    // statistics
#ifdef LAMBDA_MICRO_STATS
    ++lH.stats.numQueryWithExt;
    lH.stats.numExtScore += seqan::length(lH.matches);

    double start = sysTime();
#endif

    // Prepare string sets with sequences.
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow0>::Type> depSetH;
    seqan::StringSet<typename seqan::Source<typename TLocalHolder::TAlignRow1>::Type> depSetV;

    // container of blastMatches (possibly from multiple queries
    decltype(TBlastRecord().matches) blastMatches;

    // create blast matches
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // create blastmatch in list without copy or move
        blastMatches.emplace_back(lH.gH.qryIds [it->qryId / seqan::qNumFrames(TGlobalHolder::blastProgram)],
                                  lH.gH.indexFile.ids[it->subjId / seqan::sNumFrames(TGlobalHolder::blastProgram)]);

        auto & bm = seqan::back(blastMatches);

        bm._n_qId = it->qryId / seqan::qNumFrames(TGlobalHolder::blastProgram);
        bm._n_sId = it->subjId / seqan::sNumFrames(TGlobalHolder::blastProgram);

        bm.sLength = seqan::length(lH.gH.indexFile.seqs[bm._n_sId]);

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
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::false_type());

    // copmute evalues and filter based on evalue
    for (auto it = blastMatches.begin(), itEnd = blastMatches.end(); it != itEnd; /*below*/)
    {
        TBlastMatch & bm = *it;

        computeEValueThreadSafe(bm, seqan::length(lH.gH.qrySeqs[bm._n_qId]), seqan::context(lH.gH.outfileBlastTab));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            it = blastMatches.erase(it);
            continue;
        }

        ++it;
    }
    if (seqan::length(blastMatches) == 0)
        return 0;

    // statistics
#ifdef LAMBDA_MICRO_STATS
    lH.stats.numExtAli += seqan::length(blastMatches);
    lH.stats.timeExtend += sysTime() - start;
    start = sysTime();
#endif

    // reset and fill batches
    _setupDepSets(depSetH, depSetV, blastMatches);

    // Run extensions WITH ALIGNMENT
    _performAlignment(depSetH, depSetV, blastMatches, lH, std::true_type());

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

        seqan::computeAlignmentStats(bm, seqan::context(lH.gH.outfileBlastTab));

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

    if (seqan::length(blastMatches) == 0)
        return 0;

    // devide matches into records (per query) and write
    for (auto it = blastMatches.begin(), itLast = blastMatches.begin();
         seqan::length(blastMatches) > 0;
         /*below*/)
    {
        if ((it == blastMatches.end()) || ((it != blastMatches.begin()) && (it->_n_qId != itLast->_n_qId)))
        {
            // create a record for each query
            TBlastRecord record(lH.gH.qryIds[itLast->_n_qId]);
            record.qLength = seqan::length(lH.gH.qrySeqs[itLast->_n_qId]);
            // move the matches into the record
            record.matches.splice(record.matches.begin(),
                                  blastMatches,
                                  blastMatches.begin(),
                                  it);
            // write to file
            _writeRecord(record, lH);

            it = blastMatches.begin();
            itLast = blastMatches.begin();
        } else
        {
            itLast = it;
            ++it;
        }
    }

    return 0;
}

#endif // SEQAN_SIMD_ENABLED

template <typename TLocalHolder>
inline int
iterateMatchesFullSerial(TLocalHolder & lH)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
    using TBlastPos     = uint32_t;

    using TBlastMatch   = seqan::BlastMatch<
                           typename TLocalHolder::TAlignRow0,
                           typename TLocalHolder::TAlignRow1,
                           TBlastPos,
                           std::string_view,
                           std::string_view,
                           std::vector<std::string>, // not used
                           std::span<uint32_t>>;

    using TBlastRecord  = seqan::BlastRecord<TBlastMatch,
                                             std::string_view,
                                             std::vector<std::string>, // not used
                                             std::string_view,
                                             uint32_t>;

    auto const trueQryId = lH.matches[0].qryId / qNumFrames(TGlobalHolder::blastProgram);

    TBlastRecord record(lH.gH.qryIds[trueQryId]);
    record.qLength = qIsTranslated(TGlobalHolder::blastProgram) ?
                     lH.gH.qrySeqs[trueQryId].size() :
                     lH.gH.transQrySeqs[trueQryId].size();


    unsigned band = _bandSize(record.qLength, lH);

#ifdef LAMBDA_MICRO_STATS
    double start = sysTime();
#endif

    // create blast matches
    for (auto it = lH.matches.begin(), itEnd = lH.matches.end(); it != itEnd; ++it)
    {
        // create blastmatch in list without copy or move
        record.matches.emplace_back(lH.gH.qryIds [it->qryId / seqan::qNumFrames(TGlobalHolder::blastProgram)],
                                    lH.gH.indexFile.ids[it->subjId / seqan::sNumFrames(TGlobalHolder::blastProgram)]);

        auto & bm = record.matches.back();
        auto &  m = *it;

        bm._n_qId = it->qryId / seqan::qNumFrames(TGlobalHolder::blastProgram);
        bm._n_sId = it->subjId / seqan::sNumFrames(TGlobalHolder::blastProgram);

        bm.sLength = seqan::sIsTranslated(TGlobalHolder::blastProgram)
                        ? lH.gH.indexFile.seqs[it->subjId].size()
                        : lH.gH.transSbjSeqs[it->subjId].size();

        _setupAlignInfix(bm, *it, lH);

        _setFrames(bm, m, lH);

        // Run extension WITHOUT TRACEBACK
        typedef seqan::AlignConfig2<seqan::LocalAlignment_<>,
                                    seqan::DPBandConfig<seqan::BandOn>,
                                    seqan::FreeEndGaps_<seqan::True, seqan::True, seqan::True, seqan::True>,
                                    seqan::TracebackOff> TAlignConfig;

        seqan::DPScoutState_<seqan::Default> scoutState;

        bm.alignStats.alignmentScore =
        seqan::_setUpAndRunAlignment(lH.alignContext.dpContext,
                                     lH.alignContext.traceSegment,
                                     scoutState,
                                     seqan::source(bm.alignRow0),
                                     seqan::source(bm.alignRow1),
                                     seqan::seqanScheme(seqan::context(lH.gH.outfileBlastTab).scoringScheme),
                                     TAlignConfig(-band, +band));

        computeEValueThreadSafe(bm, record.qLength, seqan::context(lH.gH.outfileBlastTab));

        if (bm.eValue > lH.options.eCutOff)
        {
            ++lH.stats.hitsFailedExtendEValueTest;
            record.matches.pop_back();
            continue;
        }

        // Run extension WITH TRACEBACK
        seqan::localAlignment(bm.alignRow0,
                              bm.alignRow1,
                              seqan::seqanScheme(seqan::context(lH.gH.outfileBlastTab).scoringScheme),
                              -band,
                              +band);

        _expandAlign(bm, lH);

        seqan::computeAlignmentStats(bm, seqan::context(lH.gH.outfileBlastTab));

        if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        {
            ++lH.stats.hitsFailedExtendPercentIdentTest;
            record.matches.pop_back();
            continue;
        }

        seqan::computeBitScore(bm, seqan::context(lH.gH.outfileBlastTab));

        if (lH.options.hasSTaxIds)
            bm.sTaxIds = lH.gH.indexFile.sTaxIds[bm._n_sId];

    }

#ifdef LAMBDA_MICRO_STATS
    lH.stats.timeExtendTrace += sysTime() - start;
#endif

    _writeRecord(record, lH);

    return 0;
}

template <typename TLocalHolder>
inline int
iterateMatches(TLocalHolder & lH)
{
#ifdef SEQAN_SIMD_ENABLED
    if (lH.options.extensionMode == LambdaOptions::ExtensionMode::FULL_SIMD)
        return iterateMatchesFullSimd(lH);
    else
#endif

    return iterateMatchesFullSerial(lH);

}
