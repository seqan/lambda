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
// output.hpp: contains routines for file-writing
// ==========================================================================

#ifndef SEQAN_LAMBDA_OUTPUT_H_
#define SEQAN_LAMBDA_OUTPUT_H_

#include <seqan/blast.h>
#include <seqan/bam_io.h>

using namespace seqan;

template <typename TVoidSpec = void>
struct SamBamExtraTags
{
    enum Enum
    {
//         Q_START,
//         S_START,
        E_VALUE,
        BIT_SCORE,
        SCORE,
        P_IDENT,
        P_POS,
        Q_FRAME,
        S_FRAME,
        S_TAX_IDS,
        Q_AA_SEQ,
        Q_AA_CIGAR,
        EDIT_DISTANCE,
        MATCH_COUNT
    };

    static constexpr const std::array<std::pair<const char*, const char*>, 12> keyDescPairs
    {
        {
//             { "ZS", "query start (in DNA if original was DNA)" },       //  Q_START,
//             { "YS", "subject start (in DNA if original was DNA)" },  //  S_START,
            { "ZE", "expect value" },                                   //  E_VALUE,
            { "AS", "bit score" },                                      //  BIT_SCORE,
            { "ZR", "raw score" },                                      //  SCORE,
            { "ZI", "% identity (in protein space unless BLASTN) " },  //  P_IDENT,
            { "ZP", "% positive (in protein space unless BLASTN)"},     //  P_POS,
            { "ZF", "query frame" },                                    //  Q_FRAME,
            { "YF", "subject frame" },                                  //  S_FRAME,
            { "YT", "subject taxonomy IDs (* if n/a)" },                //  S_TAX_IDS,
            { "ZQ", "query protein sequence (* for BLASTN)"},          //  Q_AA_SEQ,
            { "OC", "query protein cigar (* for BLASTN)"},             //  Q_AA_CIGAR,
            { "NM", "edit distance (in protein space unless BLASTN)"}, //  EDIT_DISTANCE
            { "IH", "number of matches this query has"},                //  MATCH_COUNT
        }
    };

};

template <typename TVoidSpec>
constexpr const std::array<std::pair<const char*, const char*>, 12> SamBamExtraTags<TVoidSpec>::keyDescPairs;

// ----------------------------------------------------------------------------
// Function _untranslatedClipPositions()
// ----------------------------------------------------------------------------

// similar to _untranslatePositions() from the blast module
template <typename TSequence1, typename TSequence2, typename TNum>
inline void
_untranslateSequence(TSequence1                     & target,
                     TSequence2               const & source,
                     TNum                     const   qStart,
                     TNum                     const   qEnd,
                     int                      const   qFrameShift)
{
    if (qFrameShift >= 0)
    {
        target = infix(source,
                       3 * qStart + std::abs(qFrameShift) - 1,
                       3 * qEnd + std::abs(qFrameShift) - 1);
    }
    else
    {
        target = infix(source,
                       length(source) - (3 * qEnd + std::abs(qFrameShift) - 1),
                       length(source) - (3 * qStart + std::abs(qFrameShift) - 1));
        reverseComplement(target);
    }
}

// ----------------------------------------------------------------------------
// Function blastMatchToCigar() convert seqan align to cigar
// ----------------------------------------------------------------------------

template <typename TCigar, typename TBlastMatch, typename TLocalHolder>
inline void
blastMatchOneCigar(TCigar & cigar,
                   TBlastMatch const & m,
                   TLocalHolder const & lH)
{
    using TCElem = typename Value<TCigar>::Type;
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;

    SEQAN_ASSERT_EQ(length(m.alignRow0), length(m.alignRow1));

    // translate positions into dna space
    unsigned const transFac       = qIsTranslated(TGlobalHolder::blastProgram) ? 3 : 1;
    // clips resulting from translation / frameshift are always hard clips
    unsigned const leftFrameClip  = std::abs(m.qFrameShift) - 1;
    unsigned const rightFrameClip = qIsTranslated(TGlobalHolder::blastProgram) ? (m.qLength - leftFrameClip) % 3 : 0;
    // regular clipping from local alignment (regions outside match) can be hard or soft
    unsigned const leftClip       = m.qStart * transFac;
    unsigned const rightClip      = (length(source(m.alignRow0)) - m.qEnd) * transFac;

    if (lH.options.samBamHardClip)
    {
        if (leftFrameClip + leftClip > 0)
            appendValue(cigar, TCElem('H', leftFrameClip + leftClip));
    } else
    {
        if (leftFrameClip > 0)
            appendValue(cigar, TCElem('H', leftFrameClip));
        if (leftClip > 0)
            appendValue(cigar, TCElem('S', leftClip));
    }

    for (unsigned i = 0, count = 0; i < length(m.alignRow0); /* incremented below */)
    {
        // deletion in query
        count = 0;
        while (isGap(m.alignRow0, i) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            appendValue(cigar, TCElem('D', count * transFac));

        // insertion in query
        count = 0;
        while (isGap(m.alignRow1, i) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            appendValue(cigar, TCElem('I', count * transFac));

        // match or mismatch
        count = 0;
        while ((!isGap(m.alignRow0, i)) && (!isGap(m.alignRow1, i)) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            appendValue(cigar, TCElem('M', count * transFac));
    }

    if (lH.options.samBamHardClip)
    {
        if (rightFrameClip + rightClip > 0)
            appendValue(cigar, TCElem('H', rightFrameClip + rightClip));
    } else
    {
        if (rightClip > 0)
            appendValue(cigar, TCElem('S', rightClip));
        if (rightFrameClip > 0)
            appendValue(cigar, TCElem('H', rightFrameClip));
    }

    if (m.qFrameShift < 0)
        reverse(cigar);
}

// translation happened and we want both cigars
template <typename TCigar, typename TBlastMatch, typename TLocalHolder>
inline void
blastMatchTwoCigar(TCigar & dnaCigar,
                   TCigar & protCigar,
                   TBlastMatch const & m,
                   TLocalHolder const & lH)
{
    using TCElem = typename Value<TCigar>::Type;
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;

    SEQAN_ASSERT_EQ(length(m.alignRow0), length(m.alignRow1));

    // clips resulting from translation / frameshift are always hard clips
    unsigned const leftFrameClip  = std::abs(m.qFrameShift) - 1;            // in dna space
    unsigned const rightFrameClip = (m.qLength - leftFrameClip) % 3;          // in dna space
    // regular clipping from local alignment (regions outside match) can be hard or soft
    unsigned const leftClip       = m.qStart;                               // in protein space
    unsigned const rightClip      = length(source(m.alignRow0)) - m.qEnd;   // in protein space

    if (lH.options.samBamHardClip)
    {
        if (leftFrameClip + leftClip > 0)
            appendValue(dnaCigar, TCElem('H', leftFrameClip + 3 * leftClip));
        if (leftClip > 0)
            appendValue(protCigar, TCElem('H', leftClip));

    } else
    {
        if (leftFrameClip > 0)
            appendValue(dnaCigar, TCElem('H', leftFrameClip));

        if (leftClip > 0)
        {
            appendValue(dnaCigar, TCElem('S', 3 * leftClip));
            appendValue(protCigar, TCElem('S', leftClip));
        }
    }

    for (unsigned i = 0, count = 0; i < length(m.alignRow0); /* incremented below */)
    {
        // deletion in query
        count = 0;
        while (isGap(m.alignRow0, i) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            appendValue(dnaCigar, TCElem('D', count * 3));
            appendValue(protCigar, TCElem('D', count));
        }

        // insertion in query
        count = 0;
        while (isGap(m.alignRow1, i) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            appendValue(dnaCigar, TCElem('I', count * 3));
            appendValue(protCigar, TCElem('I', count));
        }

        // match or mismatch
        count = 0;
        while ((!isGap(m.alignRow0, i)) && (!isGap(m.alignRow1, i)) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
        {
            appendValue(dnaCigar, TCElem('M', count * 3));
            appendValue(protCigar, TCElem('M', count));
        }
    }

    if (lH.options.samBamHardClip)
    {
        if (rightFrameClip + rightClip > 0)
            appendValue(dnaCigar, TCElem('H', rightFrameClip + 3 * rightClip));
        if (rightClip > 0)
            appendValue(protCigar, TCElem('H', rightClip));

    } else
    {
        if (rightClip > 0)
        {
            appendValue(dnaCigar, TCElem('S', 3 * rightClip));
            appendValue(protCigar, TCElem('S', rightClip));
        }

        if (rightFrameClip > 0)
            appendValue(dnaCigar, TCElem('H', rightFrameClip));
    }

    if (m.qFrameShift < 0)
        reverse(dnaCigar);
    // protCigar never reversed
}

// ----------------------------------------------------------------------------
// Function myWriteHeader()
// ----------------------------------------------------------------------------

template <typename TGH, typename TLambdaOptions>
inline void
myWriteHeader(TGH & globalHolder, TLambdaOptions const & options)
{
    if (options.outFileFormat == 0) // BLAST
    {
        open(globalHolder.outfile, toCString(options.output));
        context(globalHolder.outfile).fields = options.columns;
        auto & versionString = context(globalHolder.outfile).versionString;
        clear(versionString);
        append(versionString, _programTagToString(TGH::blastProgram));
        append(versionString, " 2.2.26+ [created by LAMBDA");
        if (options.versionInformationToOutputFile)
        {
            append(versionString, "-");
            append(versionString, SEQAN_APP_VERSION);
        }
        append(versionString, ", see http://seqan.de/lambda and please cite correctly in your academic work]");
        writeHeader(globalHolder.outfile);
    } else // SAM or BAM
    {
        open(globalHolder.outfileBam, toCString(options.output));
        auto & context          = seqan::context(globalHolder.outfileBam);
        auto & subjSeqLengths   = contigLengths(context);
        auto & subjIds          = contigNames(context);

        // set sequence lengths
        if (sIsTranslated(TGH::blastProgram))
        {
            //TODO can we get around a copy?
            subjSeqLengths = globalHolder.untransSubjSeqLengths;
        } else
        {
            // compute lengths ultra-fast
            resize(subjSeqLengths, length(globalHolder.subjSeqs));
#ifdef __clang__
            SEQAN_OMP_PRAGMA(parallel for)
#else
            SEQAN_OMP_PRAGMA(parallel for simd)
#endif
            for (unsigned i = 0; i < length(subjSeqLengths); ++i)
                subjSeqLengths[i] = globalHolder.subjSeqs.limits[i+1] - globalHolder.subjSeqs.limits[i];
        }
        // set namestore
        resize(subjIds, length(globalHolder.subjIds));
        SEQAN_OMP_PRAGMA(parallel for)
        for (unsigned i = 0; i < length(globalHolder.subjIds); ++i)
            subjIds[i] = prefix(globalHolder.subjIds[i],
                                std::find(begin(globalHolder.subjIds[i], Standard()),
                                          end(globalHolder.subjIds[i], Standard()),
                                          ' ')
                                - begin(globalHolder.subjIds[i], Standard()));

        typedef BamHeaderRecord::TTag   TTag;

        // CREATE HEADER
        BamHeader header;
        // Fill first header line.
        BamHeaderRecord firstRecord;
        firstRecord.type = BAM_HEADER_FIRST;
        appendValue(firstRecord.tags, TTag("VN", "1.4"));
//         appendValue(firstRecord.tags, TTag("SO", "unsorted"));
        appendValue(firstRecord.tags, TTag("GO", "query"));
        appendValue(header, firstRecord);

        // Fill program header line.
        if (options.versionInformationToOutputFile)
        {
            BamHeaderRecord pgRecord;
            pgRecord.type = BAM_HEADER_PROGRAM;
            appendValue(pgRecord.tags, TTag("ID", "lambda"));
            appendValue(pgRecord.tags, TTag("PN", "lambda"));
            appendValue(pgRecord.tags, TTag("VN", SEQAN_APP_VERSION));
            appendValue(pgRecord.tags, TTag("CL", options.commandLine));
            appendValue(header, pgRecord);
        }

        // Fill homepage header line.
        BamHeaderRecord hpRecord0;
        hpRecord0.type = BAM_HEADER_COMMENT;
        appendValue(hpRecord0.tags, TTag("CO", "Lambda is a high performance BLAST compatible local aligner, "
                                         "please see http://seqan.de/lambda for more information."));
        appendValue(header, hpRecord0);
        BamHeaderRecord hpRecord1;
        hpRecord1.type = BAM_HEADER_COMMENT;
        appendValue(hpRecord1.tags, TTag("CO", "SAM/BAM dialect documentation is available here: "
                                         "https://github.com/seqan/lambda/wiki/Output-Formats"));
        appendValue(header, hpRecord1);
        BamHeaderRecord hpRecord2;
        hpRecord2.type = BAM_HEADER_COMMENT;
        appendValue(hpRecord2.tags, TTag("CO", "If you use any results found by Lambda, please cite "
                                         "Hauswedell et al. (2014) doi: 10.1093/bioinformatics/btu439"));
        appendValue(header, hpRecord2);

        // Fill extra tags header line.
        BamHeaderRecord tagRecord;
        tagRecord.type = BAM_HEADER_COMMENT;
        std::string columnHeaders = "Optional tags as follow";
        for (unsigned i = 0; i < length(SamBamExtraTags<>::keyDescPairs); ++i)
        {
            if (options.samBamTags[i])
            {
                columnHeaders += '\t';
                columnHeaders += std::get<0>(SamBamExtraTags<>::keyDescPairs[i]);
                columnHeaders += ':';
                columnHeaders += std::get<1>(SamBamExtraTags<>::keyDescPairs[i]);
            }
        }
        appendValue(tagRecord.tags, TTag("CO", columnHeaders));
        appendValue(header, tagRecord);

        // sam and we don't want the headers
        if (!options.samWithRefHeader && (options.outFileFormat == 1))
        {
            // we only write the header records that we actually created ourselves
            for (unsigned i = 0; i < length(header); ++i)
                write(globalHolder.outfileBam.iter, header[i], seqan::context(globalHolder.outfileBam), Sam());
        }
        else
        {
            // ref header records are automatically added with default writeHeader()
            writeHeader(globalHolder.outfileBam, header);
        }
    }
}

// ----------------------------------------------------------------------------
// Function myWriteRecord()
// ----------------------------------------------------------------------------

template <typename TLH, typename TRecord>
inline void
myWriteRecord(TLH & lH, TRecord const & record)
{
    using TGH = typename TLH::TGlobalHolder;
    if (lH.options.outFileFormat == 0) // BLAST
    {
        SEQAN_OMP_PRAGMA(critical(filewrite))
        {
            writeRecord(lH.gH.outfile, record);
        }
    } else // SAM or BAM
    {
        // convert multi-match blast-record to multiple SAM/BAM-Records

        std::vector<BamAlignmentRecord> bamRecords;
        bamRecords.resize(record.matches.size());

        String<CigarElement<>> protCigar;
        std::string protCigarString = "*";

        auto mIt = begin(record.matches, Standard());
        for (auto & bamR : bamRecords)
        {
            // untranslate for sIsTranslated
            if (sIsTranslated(TGH::blastProgram))
            {
                bamR.beginPos = mIt->sStart * 3 + std::abs(mIt->sFrameShift) - 1;
                if (mIt->sFrameShift < 0)
                    bamR.beginPos = mIt->qLength - bamR.beginPos;
            } else
            {
                bamR.beginPos   = mIt->sStart;
            }

            bamR.flag       = BAM_FLAG_SECONDARY; // all are secondary for now
            if (mIt->qFrameShift < 0)
                bamR.flag   |= BAM_FLAG_RC;
            // truncated query name
            bamR.qName      = prefix(mIt->qId,
                                     std::find(begin(mIt->qId, Standard()),
                                               end(mIt->qId, Standard()),
                                               ' ')
                                     - begin(mIt->qId, Standard()));
            // reference ID
            bamR.rID        = mIt->_n_sId;

            // compute cigar
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_CIGAR]) // amino acid cigar, too?
            {
                clear(protCigar);
                // native protein
                if ((TGH::blastProgram == BlastProgram::BLASTP) || (TGH::blastProgram == BlastProgram::TBLASTN))
                    blastMatchOneCigar(protCigar, *mIt, lH);
                else if (qIsTranslated(TGH::blastProgram)) // translated
                    blastMatchTwoCigar(bamR.cigar, protCigar, *mIt, lH);
                else // BLASTN can't have protein sequence
                    blastMatchOneCigar(bamR.cigar, *mIt, lH);
            }
            else
            {
                if ((TGH::blastProgram != BlastProgram::BLASTP) && (TGH::blastProgram != BlastProgram::TBLASTN))
                    blastMatchOneCigar(bamR.cigar, *mIt, lH);
            }
            // we want to include the seq
            bool writeSeq = false;
            if ((lH.options.samBamSeq > 1) || (mIt == begin(record.matches, Standard())))
            {
                writeSeq = true;
            }
            else if ((lH.options.samBamSeq == 1))// && lH.options.samBamHardClip)// only uniq sequences
            {
                decltype(mIt) mPrevIt = mIt - 1;
                writeSeq = ((mIt->qFrameShift              != mPrevIt->qFrameShift) ||
                            (beginPosition(mIt->alignRow0) != beginPosition(mPrevIt->alignRow0)) ||
                            (endPosition(mIt->alignRow0)   != endPosition(mPrevIt->alignRow0)));
            }

            if (TGH::blastProgram == BlastProgram::BLASTN)
            {
                if (lH.options.samBamHardClip)
                {
                    if (writeSeq)
                        bamR.seq = infix(source(mIt->alignRow0),
                                         beginPosition(mIt->alignRow0),
                                         endPosition(mIt->alignRow0));
                } else
                {
                    if (writeSeq)
                        bamR.seq = source(mIt->alignRow0);
                }
            }
            else if (qIsTranslated(TGH::blastProgram))
            {
                if (lH.options.samBamHardClip)
                {
                    if (writeSeq)
                        _untranslateSequence(bamR.seq,
                                             lH.gH.untranslatedQrySeqs[mIt->_n_qId],
                                             mIt->qStart,
                                             mIt->qEnd,
                                             mIt->qFrameShift);
                } else
                {
                    if (writeSeq)
                        _untranslateSequence(bamR.seq,
                                             lH.gH.untranslatedQrySeqs[mIt->_n_qId],
                                             decltype(length(source(mIt->alignRow0)))(0u),
                                             length(source(mIt->alignRow0)),
                                             mIt->qFrameShift);
                }
            } // else original query is protein and cannot be printed


            // custom tags
            //TODO untranslate?
//             if (lH.options.samBamTags[SamBamExtraTags<>::Q_START])
//                 appendTagValue(bamR.tags,
//                                std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_START]),
//                                uint32_t(mIt->qStart), 'I');
            //      case    S_START:
            if (lH.options.samBamTags[SamBamExtraTags<>::E_VALUE])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::E_VALUE]),
                               float(mIt->eValue), 'f');
            if (lH.options.samBamTags[SamBamExtraTags<>::BIT_SCORE])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::BIT_SCORE]),
                               uint16_t(mIt->bitScore), 'S');
            if (lH.options.samBamTags[SamBamExtraTags<>::SCORE])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::SCORE]),
                               uint8_t(mIt->alignStats.alignmentScore), 'C');
            if (lH.options.samBamTags[SamBamExtraTags<>::P_IDENT])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::P_IDENT]),
                               uint8_t(mIt->alignStats.alignmentIdentity), 'C');
            if (lH.options.samBamTags[SamBamExtraTags<>::P_POS])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::P_POS]),
                               uint16_t(mIt->alignStats.alignmentSimilarity), 'S');
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_FRAME])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_FRAME]),
                               int8_t(mIt->qFrameShift), 'c');
            if (lH.options.samBamTags[SamBamExtraTags<>::S_FRAME])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::S_FRAME]),
                               int8_t(mIt->sFrameShift), 'c');
            if (lH.options.samBamTags[SamBamExtraTags<>::S_TAX_IDS])
            {
                //TODO append integer array, instead of transforming to string
                CharString buf;
                auto it = begin(buf);
                if (length(mIt->sTaxIds) == 0)
                {
                    buf = "*";
                } else
                {
                    appendNumber(it, mIt->sTaxIds[0]);
                    for (unsigned i = 1; i < length(mIt->sTaxIds); ++i)
                    {
                        write(it, ";");
                        appendNumber(it, mIt->sTaxIds[i]);
                    }
                }
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                               buf, 'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_SEQ])
            {
                if ((TGH::blastProgram == BlastProgram::BLASTN) || (!writeSeq))
                    appendTagValue(bamR.tags,
                                   std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                                   "*", 'Z');
                else if (lH.options.samBamHardClip)
                    appendTagValue(bamR.tags,
                                   std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                                   infix(source(mIt->alignRow0),
                                         beginPosition(mIt->alignRow0),
                                         endPosition(mIt->alignRow0)),
                                   'Z');
                else // full prot sequence
                    appendTagValue(bamR.tags,
                                   std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_SEQ]),
                                   source(mIt->alignRow0),
                                   'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::Q_AA_CIGAR])
            {
                if (empty(protCigar))
                {
                    protCigarString = "*";
                }
                else
                {
                    clear(protCigarString);
                    for (unsigned i = 0; i < length(protCigar); ++i)
                    {
                        appendNumber(protCigarString, protCigar[i].count);
                        appendValue(protCigarString, protCigar[i].operation);
                    }

                }
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::Q_AA_CIGAR]),
                               protCigarString, 'Z');
            }
            if (lH.options.samBamTags[SamBamExtraTags<>::EDIT_DISTANCE])
                appendTagValue(bamR.tags,
                                std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::EDIT_DISTANCE]),
                                uint32_t(mIt->alignStats.alignmentLength - mIt->alignStats.numMatches), 'I');
            if (lH.options.samBamTags[SamBamExtraTags<>::MATCH_COUNT])
                appendTagValue(bamR.tags,
                               std::get<0>(SamBamExtraTags<>::keyDescPairs[SamBamExtraTags<>::MATCH_COUNT]),
                               uint32_t(length(record.matches)), 'I');

            // goto next match
            ++mIt;
        }

        bamRecords.front().flag -= BAM_FLAG_SECONDARY; // remove BAM_FLAG_SECONDARY for first

        SEQAN_OMP_PRAGMA(critical(filewrite))
        {
            for (auto & r : bamRecords)
                writeRecord(lH.gH.outfileBam, r);
        }
    }
}

// ----------------------------------------------------------------------------
// Function myWriteFooter()
// ----------------------------------------------------------------------------

template <typename TGH, typename TLambdaOptions>
inline void
myWriteFooter(TGH & globalHolder, TLambdaOptions const & options)
{
    if (options.outFileFormat == 0) // BLAST
    {
        writeFooter(globalHolder.outfile);
    }
}

#endif // SEQAN_LAMBDA_OUTPUT_H_
