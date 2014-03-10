// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// store.h: contains types and definitions for storing sequences and indices
// ==========================================================================

#ifndef SEQAN_LAMBDA_MISC_H_
#define SEQAN_LAMBDA_MISC_H_

#include <type_traits>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
// #include <seqan/index_extras.h>

#include <seqan/blast.h>

// #include "options.hpp"
// #include "finder.hpp"


using namespace seqan;

// typedef StringSet<CharString, Owner<ConcatDirect<> > > TDbSeqs;
// typedef StringSet<CharString, Owner<ConcatDirect<> > > TQuerySeqs;
// 
// typedef Index<TDbSeqs,      IndexSa<> > TDbIndex;
// typedef Index<TQuerySeqs,   IndexSa<> > TQueryIndex;


// ============================================================================
// Forwards
// ============================================================================






// ============================================================================
// Metafunctions
// ============================================================================

// makes partial function specialization possible

namespace ns_enabler {
    enum class enabler {};
}

template <bool Condition>
using MyEnableIf = typename std::enable_if<Condition, ns_enabler::enabler>::type;


// ============================================================================
// Functions for translation and retranslation
// ============================================================================

template <class T>
int wdth(T number)
{
    int digits = 0;
    if (number < 0) digits = 1; // remove this line if '-' counts as a digit
    while (number) {
        number /= 10;
        digits++;
    }
    return digits;
}



template <typename TAlph>
inline std::basic_ostream<char> &
operator<<(std::basic_ostream<char> & out,
           const Iter<const String<SimpleType<unsigned char,TAlph>,
                                    seqan::Packed<> >,
                      seqan::Packed<> > it)
{
    out << *it;
    return out;
}







int64_t intervalOverlap(uint64_t const s1, uint64_t const e1,
                        uint64_t const s2, uint64_t const e2)
{
    return std::min(e1, e2) - std::max(s1, s2);
}


template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec, typename TAlignContext>
TScoreValue localAlignment2(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlignContext & alignContext)
{
//     typedef Align<TSequence, TAlignSpec> TAlign;
//     typedef typename Size<TAlign>::Type TSize;
//     typedef typename Position<TAlign>::Type TPosition;
//     typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    SEQAN_ASSERT_EQ(length(rows(align)), 2u);

    clear(alignContext.traceSegment);

    TScoreValue score;

    score = _setUpAndRunAlignment(alignContext.dpContext,
                                  alignContext.traceSegment,
                                  source(row(align, 0)),
                                  source(row(align, 1)),
                                  scoringScheme,
                                  lowerDiag, upperDiag,
                                  SmithWaterman());

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), alignContext.traceSegment);
    return score;
}



// ----------------------------------------------------------------------------
// Generic Sequence loading
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TFormat>
inline int
loadSequences(StringSet<TString, TSpec > & seqs,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > ids;

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}

template <typename TString, typename TSpec, typename TFormat>
inline int
loadIds(StringSet<TString, TSpec > & ids,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > seqs;

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}

template <typename TSeqs, typename TIds, typename TFormat>
inline int
loadSeqsAndIds(TIds             & ids,
               TSeqs            & seqs,
               CharString const & path,
               TFormat    const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}


template <typename TString, typename TSpec, typename TFormat>
inline int
loadIds2(StringSet<TString, TSpec > & ids,
         CharString const & path,
         TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, SinglePass<> > TReader;
    TReader reader(stream);

    int res = 0;
    CharString seq;
    for (unsigned i = 0; i < length(ids); ++i)
    {
//         res = readRecord(ids[i], seq, reader, Fasta());
        if (value(reader) != '>')
            return 9;

        res = goNext(reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = skipBlanks(reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = readLine(ids[i], reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = skipUntilLineBeginsWithChar(reader, '>');
        if ((res) && (res != EOF_BEFORE_SUCCESS))
            std::cerr << "Error : " << res << "\n";
    }
    stream.close();
    return 0;
}


// ----------------------------------------------------------------------------
// truncate sequences
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void
_debug_shorten(StringSet<TString, TSpec > & seqs, unsigned const len)
{
    StringSet<TString, TSpec > copySeqs;
    reserve(copySeqs.concat, length(seqs)*len, Exact());

    for (TString const & s : seqs)
        if (length(s) >= len)
            appendValue(copySeqs, prefix(s, len), Exact());

    clear(seqs);
    reserve(seqs.concat, length(copySeqs)*len, Exact());
    for (TString const & s : copySeqs)
        appendValue(seqs, s);
}

// ----------------------------------------------------------------------------
// get plus-minus-range with bounds-checking for unsigned types
// ----------------------------------------------------------------------------

template <typename TNum, typename TNum2>
inline TNum
_protectUnderflow(const TNum n, const TNum2 s)
{
    const TNum r = n -s;
    return std::min(r, n);
}

template <typename TNum, typename TNum2>
inline TNum
_protectOverflow(const TNum n, const TNum2 s)
{
    const TNum r = n + s;
    return std::max(r, n);
}

template <typename TGaps>
inline bool
_startsWithGap(TGaps const & gaps)
{
    SEQAN_ASSERT(length(gaps._array) > 0);
    return (gaps._array[0] != 0);
}

template <typename TGaps>
inline int
_endsWithGap(TGaps const & gaps)
{
    SEQAN_ASSERT(length(gaps._array) > 0);
    if ((length(gaps._array)-1) % 2 == 1)
        return -1;
    return ((gaps._array[length(gaps._array)-1] != 0) ? 1 : 0);
}

template <typename TGaps, typename TSeq>
inline void
_prependNonGaps(TGaps & gaps, TSeq const & seq)
{
    if (_startsWithGap(gaps))
    {
        insertValue(gaps._array, 0, length(seq)); // new non-gap column
        insertValue(gaps._array, 0, 0); // empty gaps column
    }
    else
    {
        gaps._array[1] += length(seq);
    }

    insert(value(gaps._source), 0, seq);
    setBeginPosition(gaps, 0);
}

template <typename TGaps, typename TSeq>
inline void
_appendNonGaps(TGaps & gaps, TSeq const & seq)
{
    switch (_endsWithGap(gaps))
    {
        case -1:
        case  1:
            appendValue(gaps._array, length(seq)); // new non-gap column
            break;
        case 0:
            gaps._array[1] += length(seq);
    }
    append(value(gaps._source), seq);
    setEndPosition(gaps, length(value(gaps._source)));
}

#endif // header guard