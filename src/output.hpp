// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2015, Hannes Hauswedell, FU Berlin
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
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.
// ==========================================================================
// Author: Hannes Hauswedell <hannes.hauswedell @ fu-berlin.de>
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
    enum class Enum : uint8_t
    {
        Q_START,
//         S_START,
        E_VALUE,
        BIT_SCORE,
        SCORE,
        P_IDENT,
        P_POS,
        Q_FRAME,
        S_FRAME,
        Q_AA_SEQ,
        Q_AA_CIGAR,
        EDIT_DISTANCE,
        MATCH_COUNT
    };

    static constexpr const std::array<std::pair<const char*, const char*>, 12> keyDescPairs
    {
        {
            { "ZS", "query start (in DNA if original was DNA)" },       //  Q_START,
//             { "YS", "subject start (in DNA if original was DNA)" },  //  S_START,
            { "ZE", "expect value" },                                   //  E_VALUE,
            { "AS", "bit score" },                                      //  BIT_SCORE,
            { "ZR", "raw score" },                                      //  SCORE,
            { "ZI", "% identity (in proteine space unless BLASTN) " },  //  P_IDENT,
            { "ZP", "% positive (in proteine space unless BLASTN)"},     //  P_POS,
            { "ZF", "query frame" },                                    //  Q_FRAME,
            { "YF", "subject frame" },                                  //  S_FRAME,
            { "ZQ", "query proteine sequence (* for BLASTN)"},          //  Q_AA_SEQ,
            { "ZC", "query proteine cigar (* for BLASTN)"},             //  Q_AA_CIGAR,
            { "NM", "edit distance (in proteine space unless BLASTN)"}, //  EDIT_DISTANCE
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
template <typename TSequence1, typename TSequence2, typename TBlastMatch>
inline void
_untranslateSequence(TSequence1                     & target,
                     TSequence2               const & source,
                     TBlastMatch              const & m)
{
    if (m.qFrameShift >= 0)
    {
        target = infix(source,
                       3 * m.qStart + std::abs(m.qFrameShift) - 1,
                       3 * m.qEnd + std::abs(m.qFrameShift) - 1);
    }
    else
    {
        static thread_local Dna5String buf;
        buf = source;
        reverseComplement(buf);
        target = infix(buf,
                       3 * m.qStart + std::abs(m.qFrameShift) - 1,
                       3 * m.qEnd + std::abs(m.qFrameShift) - 1);
    }
}

// ----------------------------------------------------------------------------
// Function _untranslatedClipPositions()
// ----------------------------------------------------------------------------

// similar to _untranslatePositions() from the blast module
template <BlastProgram p>
inline void
_untranslatedClipPositions(unsigned & leftClip,         // should start with qStart
                           unsigned & rightClip,        // should start with qEnd
                           int const frameShift,
                           unsigned const realLength,   // untranslated Length if translated
                           BlastProgramSelector<p> const & selector)
{
    if (qIsTranslated(selector))
    {
        leftClip  = leftClip * 3 + std::abs(frameShift) - 1;
        rightClip = rightClip * 3 + std::abs(frameShift) - 1;
    }

    rightClip = realLength - rightClip;

    if (frameShift < 0)
        std::swap(leftClip, rightClip);
}

// ----------------------------------------------------------------------------
// Function blastMatchToCigar() convert seqan align to cigar
// ----------------------------------------------------------------------------

//TODO this could be done nicer, I guess
template <typename TCigar, typename TBlastMatch, typename TLocalHolder>
inline void
blastMatchToDnaCigar(TCigar & cigar, TBlastMatch const & m, unsigned const untransMatchLength, TLocalHolder const & lH)
{
    using TCElem = typename Value<TCigar>::Type;

    SEQAN_ASSERT_EQ(length(m.alignRow0), length(m.alignRow1));

    // hard clipping
    unsigned leftClip   = m.qStart;
    unsigned rightClip  = m.qEnd;
    _untranslatedClipPositions(leftClip,
                               rightClip,
                               m.qFrameShift,
                               untransMatchLength,
                               context(lH.gH.outfile).blastProgram);
    if (leftClip > 0)
        appendValue(cigar, TCElem('H', leftClip));

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
            appendValue(cigar, TCElem('D', (qIsTranslated(lH.gH.blastProgram) ? count * 3 : count)));

        // insertion in query
        count = 0;
        while (isGap(m.alignRow1, i) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            appendValue(cigar, TCElem('I', (qIsTranslated(lH.gH.blastProgram) ? count * 3 : count)));

        // match or mismatch
        count = 0;
        while ((!isGap(m.alignRow0, i)) && (!isGap(m.alignRow1, i)) && (i < length(m.alignRow0)))
        {
            ++count;
            ++i;
        }
        if (count > 0)
            appendValue(cigar, TCElem('M', (qIsTranslated(lH.gH.blastProgram) ? count * 3 : count)));
    }

    if (rightClip > 0)
        appendValue(cigar, TCElem('H', rightClip));
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
        writeHeader(globalHolder.outfile);
    } else // SAM or BAM
    {
        open(globalHolder.outfileBam, toCString(options.output));
        auto & context          = seqan::context(globalHolder.outfileBam);
        auto & subjSeqLengths   = contigLengths(context);
        auto & subjIds          = contigNames(context);

        // set sequence lengths
        if (sIsTranslated(globalHolder.blastProgram))
        {
            //TODO can we get around a copy?
            subjSeqLengths = globalHolder.untransSubjSeqLengths;
        } else
        {
            // compute lengths ultra-fast
            resize(subjSeqLengths, length(globalHolder.subjSeqs));
            SEQAN_OMP_PRAGMA(parallel for simd)
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
        BamHeaderRecord pgRecord;
        pgRecord.type = BAM_HEADER_PROGRAM;
        appendValue(pgRecord.tags, TTag("ID", "lambda"));
        appendValue(pgRecord.tags, TTag("PN", "lambda"));
        appendValue(pgRecord.tags, TTag("VN", SEQAN_APP_VERSION));
        appendValue(pgRecord.tags, TTag("CL", options.commandLine));
        appendValue(header, pgRecord);

        // TODO add comment line describing extra fields and possibly link to homepage

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

        auto mIt = begin(record.matches, Standard());
        for (auto & bamR : bamRecords)
        {
            bamR.beginPos   = mIt->sStart;

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
            blastMatchToDnaCigar(bamR.cigar, *mIt, record.qLength, lH);
            // we want to include the seq
            if (lH.options.samBamSeq)
            {
                // only dna sequences supported
                if (lH.gH.blastProgram == BlastProgram::BLASTN)
                    bamR.seq = infix(source(mIt->alignRow0),
                                     beginPosition(mIt->alignRow0),
                                     endPosition(mIt->alignRow0));
                // untranslation is ok, too
                else if (qIsTranslated(lH.gH.blastProgram))
                    _untranslateSequence(bamR.seq,
                                         lH.gH.untranslatedQrySeqs[mIt->_n_qId],
                                        * mIt);
                // else no sequence is available
            }

            // custom tags

            std::vector<SamBamExtraTags<>::Enum> tags;
            appendValue(tags, SamBamExtraTags<>::Enum::BIT_SCORE);
            appendValue(tags, SamBamExtraTags<>::Enum::EDIT_DISTANCE);
            appendValue(tags, SamBamExtraTags<>::Enum::E_VALUE);
            appendValue(tags, SamBamExtraTags<>::Enum::P_IDENT);
            appendValue(tags, SamBamExtraTags<>::Enum::P_POS);
            appendValue(tags, SamBamExtraTags<>::Enum::Q_FRAME);

            for (auto const & tag : lH.options.samBamColumns)
            {
                switch (tag)
                {
                    case SamBamExtraTags<>::Enum::Q_START:
                        // TODO untranslate
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint32_t(mIt->qLength), 'I');
                        break;
            //      case    S_START:
                    case SamBamExtraTags<>::Enum::E_VALUE:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       float(mIt->eValue), 'f');
                        break;
                    case SamBamExtraTags<>::Enum::BIT_SCORE:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint16_t(mIt->bitScore), 'S');
                        break;
                    case SamBamExtraTags<>::Enum::SCORE:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint8_t(mIt->alignStats.alignmentScore), 'C');
                        break;
                    case SamBamExtraTags<>::Enum::P_IDENT:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint8_t(mIt->alignStats.alignmentIdentity), 'C');
                        break;
                    case SamBamExtraTags<>::Enum::P_POS:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint16_t(mIt->alignStats.alignmentSimilarity), 'S');
                        break;
                    case SamBamExtraTags<>::Enum::Q_FRAME:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       int8_t(mIt->qFrameShift), 'c');
                        break;
                    case SamBamExtraTags<>::Enum::S_FRAME:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       int8_t(mIt->sFrameShift), 'c');
                        break;
                    case SamBamExtraTags<>::Enum::Q_AA_SEQ:
                        if (lH.gH.blastProgram == BlastProgram::BLASTN)
                            appendTagValue(bamR.tags,
                                           std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                           "*", 'Z');
                        else
                            appendTagValue(bamR.tags,
                                           std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                           infix(source(mIt->alignRow0),
                                                 beginPosition(mIt->alignRow0),
                                                 endPosition(mIt->alignRow0)),
                                           'Z');
                        break;
                    case SamBamExtraTags<>::Enum::Q_AA_CIGAR:
                        //TODO
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       "*", 'Z');
                        break;
                    case SamBamExtraTags<>::Enum::EDIT_DISTANCE:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint32_t(mIt->alignStats.alignmentLength - mIt->alignStats.numMatches), 'I');
                        break;
                    case SamBamExtraTags<>::Enum::MATCH_COUNT:
                        appendTagValue(bamR.tags,
                                       std::get<0>(SamBamExtraTags<>::keyDescPairs[(uint8_t)tag]),
                                       uint32_t(length(record.matches)), 'I');
                        break;
                }
            }
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
