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
// search_output.hpp: contains routines for file-writing
// ==========================================================================

#pragma once

#include <seqan/bam_io.h>
#include <seqan/blast.h>
#include <seqan/modifier.h>

template <typename TVoidSpec = void>
struct SamBamExtraTags
{
    enum Enum
    {
        //         Q_START,
        //         S_START,
        BIT_SCORE,
        Q_AA_CIGAR,
        EDIT_DISTANCE,
        MATCH_COUNT,
        SCORE,
        E_VALUE,
        P_IDENT,
        P_POS,
        Q_FRAME,
        Q_AA_SEQ,
        S_FRAME,
        S_TAX_IDS,
        LCA_ID,
        LCA_TAX_ID
    };

    // clang-format off
    static constexpr const std::array<std::pair<const char*, const char*>, 14> keyDescPairs
    {
        {
//             { "ZS", "query start (in DNA if original was DNA)" },       //  Q_START,
//             { "YS", "subject start (in DNA if original was DNA)" },  //  S_START,
            { "AS", "bit score" },                                      //  BIT_SCORE,
            { "OC", "query protein cigar (* for BLASTN)"},              //  Q_AA_CIGAR,
            { "NM", "edit distance (in protein space unless BLASTN)"},  //  EDIT_DISTANCE
            { "IH", "number of matches this query has"},                //  MATCH_COUNT
            { "ar", "raw score" },                                      //  SCORE,
            { "ae", "expect value" },                                   //  E_VALUE,
            { "ai", "% identity (in protein space unless BLASTN) " },   //  P_IDENT,
            { "ap", "% positive (in protein space unless BLASTN)"},     //  P_POS,
            { "qf", "query frame" },                                    //  Q_FRAME,
            { "qs", "query protein sequence (* for BLASTN)"},           //  Q_AA_SEQ,
            { "sf", "subject frame" },                                  //  S_FRAME,
            { "st", "subject taxonomy IDs (* if n/a)" },                //  S_TAX_IDS,
            { "ls", "lowest common ancestor scientific name" },         //  LCA_ID,
            { "lt", "lowest common ancestor taxonomy ID" },             //  LCA_TAX_ID,
        }
    };
    // clang-format on
};

template <typename TVoidSpec>
constexpr const std::array<std::pair<char const *, char const *>, 14> SamBamExtraTags<TVoidSpec>::keyDescPairs;

// ----------------------------------------------------------------------------
// Function _untranslatedClipPositions()
// ----------------------------------------------------------------------------

// similar to _untranslatePositions() from the blast module
template <typename TSequence1, typename TSequence2, typename TNum>
inline void _untranslateSequence(TSequence1 &       target,
                                 TSequence2 const & source,
                                 TNum const         qStart,
                                 TNum const         qEnd,
                                 int const          qFrameShift)
{
    if (qFrameShift >= 0)
    {
        seqan::copy_range(
          source | seqan3::views::slice(3 * qStart + std::abs(qFrameShift) - 1, 3 * qEnd + std::abs(qFrameShift) - 1) |
            seqan3::views::to_char,
          target);
    }
    else
    {
        seqan::copy_range(source |
                            seqan3::views::slice(seqan::length(source) - (3 * qEnd + std::abs(qFrameShift) - 1),
                                                 seqan::length(source) - (3 * qStart + std::abs(qFrameShift) - 1)) |
                            seqan3::views::to_char,
                          target);

        reverseComplement(target);
    }
}

// ----------------------------------------------------------------------------
// Function blastMatchToCigar() convert seqan align to cigar
// ----------------------------------------------------------------------------

template <typename TCigar, typename TBlastMatch, typename TBlastRecord, typename TLocalHolder>
inline void blastMatchOneCigar(TCigar & cigar, TBlastMatch const & m, TBlastRecord const & r, TLocalHolder const & lH)
{
    using TCElem        = typename seqan::Value<TCigar>::Type;
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;

    SEQAN_ASSERT_EQ(seqan::length(m.alignRow0), seqan::length(m.alignRow1));

    // translate positions into dna space
    unsigned const transFac       = qIsTranslated(TGlobalHolder::blastProgram) ? 3 : 1;
    // clips resulting from translation / frameshift are always hard clips
    unsigned const leftFrameClip  = std::abs(m.qFrameShift) - 1;
    unsigned const rightFrameClip = qIsTranslated(TGlobalHolder::blastProgram) ? (r.qLength - leftFrameClip) % 3 : 0;

    // regular clipping from local alignment (regions outside match) can be hard or soft
    unsigned const leftClip  = m.qStart * transFac;
    unsigned const rightClip = (seqan::length(seqan::source(m.alignRow0)) - m.qEnd) * transFac;

    if (lH.options.samBamHardClip)
    {
        if (leftFrameClip + leftClip > 0)
            seqan::appendValue(cigar, TCElem('H', leftFrameClip + leftClip));
    }
    else
    {
        if (leftFrameClip > 0)
            seqan::appendValue(cigar, TCElem('H', leftFrameClip));
        if (leftClip > 0)
            seqan::appendValue(cigar, TCElem('S', leftClip));
    }

    for (unsigned i = 0, count = 0; i < seqan::length(m.alignRow0); /* incremented below */)
    {
        // deletion in query
        count = 0;
        while (seqan::isGap(m.alignRow0, i) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            seqan::appendValue(cigar, TCElem('D', count * transFac));

        // insertion in query
        count = 0;
        while (seqan::isGap(m.alignRow1, i) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            seqan::appendValue(cigar, TCElem('I', count * transFac));

        // match or mismatch
        count = 0;
        while ((!seqan::isGap(m.alignRow0, i)) && (!seqan::isGap(m.alignRow1, i)) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            seqan::appendValue(cigar, TCElem('M', count * transFac));
    }

    if (lH.options.samBamHardClip)
    {
        if (rightFrameClip + rightClip > 0)
            seqan::appendValue(cigar, TCElem('H', rightFrameClip + rightClip));
    }
    else
    {
        if (rightClip > 0)
            seqan::appendValue(cigar, TCElem('S', rightClip));
        if (rightFrameClip > 0)
            seqan::appendValue(cigar, TCElem('H', rightFrameClip));
    }

    if (m.qFrameShift < 0)
        reverse(cigar);
}

// translation happened and we want both cigars
template <typename TCigar, typename TBlastMatch, typename TBlastRecord, typename TLocalHolder>
inline void blastMatchTwoCigar(TCigar &             dnaCigar,
                               TCigar &             protCigar,
                               TBlastMatch const &  m,
                               TBlastRecord const & r,
                               TLocalHolder const & lH)
{
    using TCElem = typename seqan::Value<TCigar>::Type;

    SEQAN_ASSERT_EQ(seqan::length(m.alignRow0), seqan::length(m.alignRow1));

    // clips resulting from translation / frameshift are always hard clips
    unsigned const leftFrameClip  = std::abs(m.qFrameShift) - 1;     // in dna space
    unsigned const rightFrameClip = (r.qLength - leftFrameClip) % 3; // in dna space
    // regular clipping from local alignment (regions outside match) can be hard or soft
    unsigned const leftClip       = m.qStart;                                           // in protein space
    unsigned const rightClip      = seqan::length(seqan::source(m.alignRow0)) - m.qEnd; // in protein space

    if (lH.options.samBamHardClip)
    {
        if (leftFrameClip + leftClip > 0)
            seqan::appendValue(dnaCigar, TCElem('H', leftFrameClip + 3 * leftClip));
        if (leftClip > 0)
            seqan::appendValue(protCigar, TCElem('H', leftClip));
    }
    else
    {
        if (leftFrameClip > 0)
            seqan::appendValue(dnaCigar, TCElem('H', leftFrameClip));

        if (leftClip > 0)
        {
            seqan::appendValue(dnaCigar, TCElem('S', 3 * leftClip));
            seqan::appendValue(protCigar, TCElem('S', leftClip));
        }
    }

    for (unsigned i = 0, count = 0; i < seqan::length(m.alignRow0); /* incremented below */)
    {
        // deletion in query
        count = 0;
        while (seqan::isGap(m.alignRow0, i) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            seqan::appendValue(dnaCigar, TCElem('D', count * 3));
            seqan::appendValue(protCigar, TCElem('D', count));
        }

        // insertion in query
        count = 0;
        while (seqan::isGap(m.alignRow1, i) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            seqan::appendValue(dnaCigar, TCElem('I', count * 3));
            seqan::appendValue(protCigar, TCElem('I', count));
        }

        // match or mismatch
        count = 0;
        while ((!seqan::isGap(m.alignRow0, i)) && (!seqan::isGap(m.alignRow1, i)) && (i < seqan::length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            seqan::appendValue(dnaCigar, TCElem('M', count * 3));
            seqan::appendValue(protCigar, TCElem('M', count));
        }
    }

    if (lH.options.samBamHardClip)
    {
        if (rightFrameClip + rightClip > 0)
            seqan::appendValue(dnaCigar, TCElem('H', rightFrameClip + 3 * rightClip));
        if (rightClip > 0)
            seqan::appendValue(protCigar, TCElem('H', rightClip));
    }
    else
    {
        if (rightClip > 0)
        {
            seqan::appendValue(dnaCigar, TCElem('S', 3 * rightClip));
            seqan::appendValue(protCigar, TCElem('S', rightClip));
        }

        if (rightFrameClip > 0)
            seqan::appendValue(dnaCigar, TCElem('H', rightFrameClip));
    }

    if (m.qFrameShift < 0)
        seqan::reverse(dnaCigar);
    // protCigar never reversed
}

// ----------------------------------------------------------------------------
// Function myWriteHeader()
// ----------------------------------------------------------------------------

template <typename TGH, typename TLambdaOptions>
inline void myWriteHeader(TGH & globalHolder, TLambdaOptions const & options)
{
    if (options.outFileFormat <= 0) // Blast
    {
        std::string versionString;
        seqan::append(versionString, _programTagToString(TGH::blastProgram));
        seqan::append(versionString, " 2.2.26+ [created by LAMBDA");
        if (options.versionInformationToOutputFile)
        {
            seqan::append(versionString, "-");
            seqan::append(versionString, SEQAN_APP_VERSION);
        }
        seqan::append(versionString, ", see http://seqan.de/lambda and please cite correctly in your academic work]");

        if (options.outFileFormat == -1) // BLAST-rep
        {
            // copy details from tabular file
            seqan::context(globalHolder.outfileBlastRep) = seqan::context(globalHolder.outfileBlastTab);

            if (options.versionInformationToOutputFile)
                seqan::context(globalHolder.outfileBlastRep).versionString += versionString;

            seqan::open(globalHolder.outfileBlastRep, options.output.c_str());
            seqan::context(globalHolder.outfileBlastRep).fields = options.columns;
            seqan::writeHeader(globalHolder.outfileBlastRep);
        }
        else // BLAST-tab
        {
            if (options.blastTabularWithComments)
                seqan::context(globalHolder.outfileBlastTab).tabularSpec = seqan::BlastTabularSpec::COMMENTS;
            else
                seqan::context(globalHolder.outfileBlastTab).tabularSpec = seqan::BlastTabularSpec::NO_COMMENTS;

            if (options.versionInformationToOutputFile)
                seqan::context(globalHolder.outfileBlastTab).versionString += versionString;

            seqan::open(globalHolder.outfileBlastTab, options.output.c_str());
            seqan::context(globalHolder.outfileBlastTab).fields = options.columns;
            seqan::writeHeader(globalHolder.outfileBlastTab);
        }
    }
    else // SAM or BAM
    {
        seqan::open(globalHolder.outfileBam, options.output.c_str());
        auto & context        = seqan::context(globalHolder.outfileBam);
        auto & subjSeqLengths = seqan::contigLengths(context);
        auto & subjIds        = seqan::contigNames(context);

        // compute seqan::lengths ultra-fast
        seqan::resize(subjSeqLengths, globalHolder.indexFile.seqs.size());
        SEQAN_OMP_PRAGMA(parallel for simd)
        for (size_t i = 0; i < globalHolder.indexFile.seqs.size(); ++i)
            subjSeqLengths[i] = globalHolder.indexFile.seqs[i].size();

        // set namestore
        resize(subjIds, globalHolder.indexFile.ids.size());
        SEQAN_OMP_PRAGMA(parallel for)
        for (unsigned i = 0; i < globalHolder.indexFile.ids.size(); ++i)
        {
            //TODO replace with assign algo
            subjIds[i] = globalHolder.indexFile.ids[i] | seqan3::detail::take_until(seqan3::is_space) |
                         seqan3::ranges::to<std::string>();
        }

        typedef seqan::BamHeaderRecord::TTag TTag;

        // CREATE HEADER
        seqan::BamHeader       header;
        // Fill first header line.
        seqan::BamHeaderRecord firstRecord;
        firstRecord.type = seqan::BAM_HEADER_FIRST;
        seqan::appendValue(firstRecord.tags, TTag("VN", "1.4"));
        //         seqan::appendValue(firstRecord.tags, TTag("SO", "unsorted"));
        seqan::appendValue(firstRecord.tags, TTag("GO", "query"));
        seqan::appendValue(header, firstRecord);

        // Fill program header line.
        if (options.versionInformationToOutputFile)
        {
            seqan::BamHeaderRecord pgRecord;
            pgRecord.type = seqan::BAM_HEADER_PROGRAM;
            seqan::appendValue(pgRecord.tags, TTag("ID", "lambda"));
            seqan::appendValue(pgRecord.tags, TTag("PN", "lambda"));
            seqan::appendValue(pgRecord.tags, TTag("VN", SEQAN_APP_VERSION));
            seqan::appendValue(pgRecord.tags, TTag("CL", options.commandLine));
            seqan::appendValue(header, pgRecord);
        }

        // Fill homepage header line.
        seqan::BamHeaderRecord hpRecord0;
        hpRecord0.type = seqan::BAM_HEADER_COMMENT;
        seqan::appendValue(hpRecord0.tags,
                           TTag("CO",
                                "Lambda is a high performance BLAST compatible local aligner, "
                                "please see http://seqan.de/lambda for more information."));
        seqan::appendValue(header, hpRecord0);
        seqan::BamHeaderRecord hpRecord1;
        hpRecord1.type = seqan::BAM_HEADER_COMMENT;
        seqan::appendValue(hpRecord1.tags,
                           TTag("CO",
                                "SAM/BAM dialect documentation is available here: "
                                "https://github.com/seqan/lambda/wiki/Output-Formats"));
        seqan::appendValue(header, hpRecord1);
        seqan::BamHeaderRecord hpRecord2;
        hpRecord2.type = seqan::BAM_HEADER_COMMENT;
        seqan::appendValue(hpRecord2.tags,
                           TTag("CO",
                                "If you use any results found by Lambda, please cite "
                                "Hauswedell et al. (2014) doi: 10.1093/bioinformatics/btu439"));
        seqan::appendValue(header, hpRecord2);

        // Fill extra tags header line.
        seqan::BamHeaderRecord tagRecord;
        tagRecord.type            = seqan::BAM_HEADER_COMMENT;
        std::string columnHeaders = "Optional tags as follow";
        for (unsigned i = 0; i < seqan::length(SamBamExtraTags<>::keyDescPairs); ++i)
        {
            if (options.samBamTags[i])
            {
                columnHeaders += '\t';
                columnHeaders += std::get<0>(SamBamExtraTags<>::keyDescPairs[i]);
                columnHeaders += ':';
                columnHeaders += std::get<1>(SamBamExtraTags<>::keyDescPairs[i]);
            }
        }
        seqan::appendValue(tagRecord.tags, TTag("CO", columnHeaders));
        seqan::appendValue(header, tagRecord);

        // sam and we don't want the headers
        if (!options.samWithRefHeader && (options.outFileFormat == 1))
        {
            // we only write the header records that we actually created ourselves
            for (unsigned i = 0; i < seqan::length(header); ++i)
                seqan::write(globalHolder.outfileBam.iter,
                             header[i],
                             seqan::context(globalHolder.outfileBam),
                             seqan::Sam());
            std::cout << __FILE__ << ": " << __LINE__ << '\n';
        }
        else
        {
            // ref header records are automatically added with default writeHeader()
            seqan::writeHeader(globalHolder.outfileBam, header);
            std::cout << __FILE__ << ": " << __LINE__ << '\n';
        }
    }
}

// ----------------------------------------------------------------------------
// Function myWriteRecord()
// ----------------------------------------------------------------------------

template <typename TLH, typename TRecord>
inline void myWriteRecord(TLH & lH, TRecord const & record)
{
    using TGH = typename TLH::TGlobalHolder;
    if (lH.options.outFileFormat == 0) // BLAST-Tab
    {
        SEQAN_OMP_PRAGMA(critical(filewrite))
        {
            seqan::writeRecord(lH.gH.outfileBlastTab, record);
        }
    }
    else if (lH.options.outFileFormat == -1) // BLAST
    {
        SEQAN_OMP_PRAGMA(critical(filewrite))
        {
            seqan::writeRecord(lH.gH.outfileBlastRep, record);
        }
    }
    else // SAM or BAM
    {
        // convert multi-match blast-record to multiple SAM/BAM-Records
        std::vector<seqan::BamAlignmentRecord> bamRecords;
        bamRecords.resize(record.matches.size());

        seqan::String<seqan::CigarElement<>> protCigar;
        std::string                          protCigarString = "*";

        auto mIt = seqan::begin(record.matches, seqan::Standard());
        for (auto & bamR : bamRecords)
        {
            // untranslate for sIsTranslated
            if constexpr (seqan::sIsTranslated(TGH::blastProgram))
            {
                bamR.beginPos = mIt->sStart * 3 + std::abs(mIt->sFrameShift) - 1;
                if (mIt->sFrameShift < 0)
                    bamR.beginPos = record.qLength - bamR.beginPos;
            }
            else
            {
                bamR.beginPos = mIt->sStart;
            }

            bamR.flag = seqan::BAM_FLAG_SECONDARY; // all are secondary for now
            if (mIt->qFrameShift < 0)
                bamR.flag |= seqan::BAM_FLAG_RC;
            // truncated query name+
            for (char c : record.qId | seqan3::detail::take_until(seqan3::is_space))
                seqan::appendValue(bamR.qName, c);
            // reference ID
            bamR.rID = mIt->_n_sId;

            // compute cigar
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_CIGAR]) // amino acid cigar, too?
            {
                seqan::clear(protCigar);
                // native protein
                if constexpr ((TGH::blastProgram == seqan::BlastProgram::BLASTP) ||
                              (TGH::blastProgram == seqan::BlastProgram::TBLASTN))
                    blastMatchOneCigar(protCigar, *mIt, record, lH);
                else if constexpr (qIsTranslated(TGH::blastProgram)) // translated
                    blastMatchTwoCigar(bamR.cigar, protCigar, *mIt, record, lH);
                else // BLASTN can't have protein sequence
                    blastMatchOneCigar(bamR.cigar, *mIt, record, lH);
            }
            else
            {
                if constexpr ((TGH::blastProgram != seqan::BlastProgram::BLASTP) &&
                              (TGH::blastProgram != seqan::BlastProgram::TBLASTN))
                    blastMatchOneCigar(bamR.cigar, *mIt, record, lH);
            }
            // we want to include the seq
            bool writeSeq = false;
            if (lH.options.samBamSeq > 1)
            {
                writeSeq = true;
            }
            else if (lH.options.samBamSeq == 1) // only uniq sequences
            {
                if (mIt == begin(record.matches, seqan::Standard()))
                {
                    writeSeq = true;
                }
                else
                {
                    decltype(mIt) mPrevIt = mIt - 1;
                    writeSeq              = ((mIt->qFrameShift != mPrevIt->qFrameShift) ||
                                (seqan::beginPosition(mIt->alignRow0) != seqan::beginPosition(mPrevIt->alignRow0)) ||
                                (seqan::endPosition(mIt->alignRow0) != seqan::endPosition(mPrevIt->alignRow0)));
                }
            }

            if constexpr (TGH::blastProgram == seqan::BlastProgram::BLASTN)
            {
                if (lH.options.samBamHardClip)
                {
                    if (writeSeq)
                    {
                        seqan::copy_range(
                          seqan::source(mIt->alignRow0) |
                            seqan3::views::slice(seqan::beginPosition(mIt->alignRow0),
                                                 seqan::endPosition(mIt->alignRow0)) |
                            std::views::transform([](seqan3::dna5 a) { return seqan::Iupac{seqan3::to_char(a)}; }),
                          bamR.seq);
                    }
                }
                else
                {
                    if (writeSeq)
                    {
                        seqan::copy_range(
                          seqan::source(mIt->alignRow0) |
                            std::views::transform([](seqan3::dna5 a) { return seqan::Iupac{seqan3::to_char(a)}; }),
                          bamR.seq);
                    }
                }
            }
            else if constexpr (qIsTranslated(TGH::blastProgram))
            {
                if (lH.options.samBamHardClip)
                {
                    if (writeSeq)
                        _untranslateSequence(bamR.seq,
                                             lH.qrySeqs[mIt->_n_qId],
                                             mIt->qStart,
                                             mIt->qEnd,
                                             mIt->qFrameShift);
                }
                else
                {
                    if (writeSeq)
                        _untranslateSequence(bamR.seq,
                                             lH.qrySeqs[mIt->_n_qId],
                                             decltype(seqan::length(seqan::source(mIt->alignRow0)))(0u),
                                             seqan::length(seqan::source(mIt->alignRow0)),
                                             mIt->qFrameShift);
                }
            } // else original query is protein and cannot be printed

            // custom tags
            //TODO untranslate?
            //             if (lH.options.samBamTags[SamBamExtraTags<>::Q_START])
            //                 seqan::appendTagValue(bamR.tags,
            //                                std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_START]),
            //                                uint32_t(mIt->qStart), 'I');
            //      case    S_START:
            if (lH.options.samBamTags[SamBamExtraTags<>::E_VALUE])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::E_VALUE]),
                                      float(mIt->eValue),
                                      'f');
            if (lH.options.samBamTags[SamBamExtraTags<>::BIT_SCORE])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::BIT_SCORE]),
                                      uint16_t(mIt->bitScore),
                                      'S');
            if (lH.options.samBamTags[SamBamExtraTags<>::SCORE])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::SCORE]),
                                      uint8_t(mIt->alignStats.alignmentScore),
                                      'C');
            if (lH.options.samBamTags[SamBamExtraTags<>::P_IDENT])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::P_IDENT]),
                                      uint8_t(mIt->alignStats.alignmentIdentity),
                                      'C');
            if (lH.options.samBamTags[SamBamExtraTags<>::P_POS])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::P_POS]),
                                      uint16_t(mIt->alignStats.alignmentSimilarity),
                                      'S');
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_FRAME])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_FRAME]),
                                      int8_t(mIt->qFrameShift),
                                      'c');
            if (lH.options.samBamTags[SamBamExtraTags<>::S_FRAME])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::S_FRAME]),
                                      int8_t(mIt->sFrameShift),
                                      'c');
            if (lH.options.samBamTags[SamBamExtraTags<>::S_TAX_IDS])
            {
                //TODO append integer array, instead of transforming to string
                seqan::CharString buf;
                auto              it = begin(buf);
                if (seqan::length(mIt->sTaxIds) == 0)
                {
                    buf = "*";
                }
                else
                {
                    seqan::appendNumber(it, mIt->sTaxIds[0]);
                    for (unsigned i = 1; i < seqan::length(mIt->sTaxIds); ++i)
                    {
                        write(it, ";");
                        seqan::appendNumber(it, mIt->sTaxIds[i]);
                    }
                }
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::S_TAX_IDS]),
                                      buf,
                                      'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::LCA_ID])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::LCA_ID]),
                                      record.lcaId,
                                      'Z');
            if (lH.options.samBamTags[SamBamExtraTags<>::LCA_TAX_ID])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::LCA_TAX_ID]),
                                      uint32_t(record.lcaTaxId),
                                      'I');
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_SEQ])
            {
                if ((TGH::blastProgram == seqan::BlastProgram::BLASTN) || (!writeSeq))
                    seqan::appendTagValue(bamR.tags,
                                          std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                                          "*",
                                          'Z');
                else if (lH.options.samBamHardClip)
                    seqan::appendTagValue(
                      bamR.tags,
                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                      seqan::source(mIt->alignRow0) |
                        seqan3::views::slice(seqan::beginPosition(mIt->alignRow0), seqan::endPosition(mIt->alignRow0)) |
                        seqan3::views::to_char,
                      'Z');
                else // full prot sequence
                    seqan::appendTagValue(bamR.tags,
                                          std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                                          seqan::source(mIt->alignRow0) | seqan3::views::to_char,
                                          'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_CIGAR])
            {
                if (seqan::empty(protCigar))
                {
                    protCigarString = "*";
                }
                else
                {
                    seqan::clear(protCigarString);
                    for (unsigned i = 0; i < seqan::length(protCigar); ++i)
                    {
                        seqan::appendNumber(protCigarString, protCigar[i].count);
                        seqan::appendValue(protCigarString, protCigar[i].operation);
                    }
                }
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_CIGAR]),
                                      protCigarString,
                                      'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::EDIT_DISTANCE])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::EDIT_DISTANCE]),
                                      uint32_t(mIt->alignStats.alignmentLength - mIt->alignStats.numMatches),
                                      'I');
            if (lH.options.samBamTags[SamBamExtraTags<>::MATCH_COUNT])
                seqan::appendTagValue(bamR.tags,
                                      std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::MATCH_COUNT]),
                                      uint32_t(seqan::length(record.matches)),
                                      'I');

            // goto next match
            ++mIt;
        }

        bamRecords.front().flag -= seqan::BAM_FLAG_SECONDARY; // remove BAM_FLAG_SECONDARY for first

        SEQAN_OMP_PRAGMA(critical(filewrite))
        {
            for (auto & r : bamRecords)
                seqan::writeRecord(lH.gH.outfileBam, r);
        }
    }
}

// ----------------------------------------------------------------------------
// Function myWriteFooter()
// ----------------------------------------------------------------------------

template <typename TGH, typename TLambdaOptions>
inline void myWriteFooter(TGH & globalHolder, TLambdaOptions const & options)
{
    if (options.outFileFormat == 0) // BLAST
    {
        seqan::writeFooter(globalHolder.outfileBlastTab);
    }
    else if (options.outFileFormat == -1) // BLAST
    {
        seqan::writeFooter(globalHolder.outfileBlastRep);
    }
}
