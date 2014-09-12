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
// lambda_indexer.hpp: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_LAMBDA_LAMBDA_INDEXER_H_
#define SEQAN_LAMBDA_LAMBDA_INDEXER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/seq_io.h>
#include <seqan/index.h>
// #include <seqan/index_extras.h>
#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include "misc.hpp"
#include "options.hpp"
#include "trans.hpp"
#include "alph.hpp"



using namespace seqan;


// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TOrigAlph>
inline int
loadSubjSeqsAndIds(TCDStringSet<TOrigAlph> & originalSeqs,
                   LambdaIndexerOptions const & options)
{
    int ret = 0;
    StringSet<CharString, Owner<ConcatDirect<>>> ids;

    double start = sysTime();
    std::cout << "Loading Subject Sequences and Ids…" << std::flush;


    if (options.fileFormat)
        ret = loadSeqsAndIds(ids,
                             originalSeqs,
                             options.dbFile,
                             Fastq());
    else
        ret = loadSeqsAndIds(ids,
                             originalSeqs,
                             options.dbFile,
                             Fasta());
    if (ret)
        return ret;

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    unsigned long maxLen = 0ul;
    for (auto const & s : originalSeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(originalSeqs)
            << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    std::cout << "Dumping Subj Ids..." << std::flush;

    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".ids");
    save(ids, toCString(_path));

    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TLimits>
inline void
_saveOriginalSeqLengths(TLimits limits, // we want copy!
                       LambdaIndexerOptions const & options,
                       True const & /*SHasFrames*/)
{
    for (uint32_t i = 0; i < (length(limits) - 1); ++i)
        limits[i] = limits[i+1] - limits[i];
    // last entry not overwritten, should be the sum of all lengths

    std::cout << " dumping untranslated subject lengths..." << std::flush;
    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".untranslengths");
    save(limits, toCString(_path));
}

template <typename TLimits>
inline void
_saveOriginalSeqLengths(TLimits const &/**/,
                       LambdaIndexerOptions const & /**/,
                       False const & /*SHasFrames*/)
{
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TOrigAlph>
inline void
translateOrSwap(TCDStringSet<TTransAlph> & out,
                TCDStringSet<TOrigAlph> & in,
                LambdaIndexerOptions const & options)
{
    //TODO more output
    std::cout << "translating…" << std::flush;
    translate(out,
              in,
              SIX_FRAME,
              options.geneticCode);
}

template <typename TSameAlph>
inline void
translateOrSwap(TCDStringSet<TSameAlph> & out,
                TCDStringSet<TSameAlph> & in,
                LambdaIndexerOptions const & /**/)
{
    swap(out, in);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph>
inline void
_dumpTranslatedSeqs(TCDStringSet<TTransAlph> const & translatedSeqs,
                    LambdaIndexerOptions const & options)
{
    double start = sysTime();
    std::cout << "Dumping unreduced Subj Sequences..." << std::flush;

    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".unredsubj");
    save(translatedSeqs, toCString(_path));

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;
}

template <typename TTransAlph>
inline void
dumpTranslatedSeqs(TCDStringSet<TTransAlph> const & translatedSeqs,
                   LambdaIndexerOptions const & options)
{
    if (options.alphReduction > 0)
        _dumpTranslatedSeqs(translatedSeqs, options);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TRedAlph>
inline void
reduceOrSwap(TCDStringSet<TRedAlph> & out,
             TCDStringSet<TTransAlph> & in)
{
    //TODO more output
    // reduce implicitly
    std::cout << "reducing…" << std::flush;
    out.concat = in.concat;
    out.limits = in.limits;
}

template <typename TSameAlph>
inline void
reduceOrSwap(TCDStringSet<TSameAlph> & out,
             TCDStringSet<TSameAlph> & in)
{
    swap(out, in);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TRedAlph>
inline bool
checkIndexSize(TCDStringSet<TRedAlph> const & seqs)
{
    using SAV = typename SAValue<TCDStringSet<TRedAlph>>::Type;
    uint64_t curNumSeq = length(seqs);
    uint64_t maxNumSeq = std::numeric_limits<typename Value<SAV, 1>::Type>::max();

    if (curNumSeq >= maxNumSeq)
    {
        std::cerr << "Too many sequences to be indexed:\n  "
                  << length(seqs) << " in file, but only "
                  << maxNumSeq << " supported by index.\n";
        return false;
    }

    uint64_t maxLenSeq = std::numeric_limits<typename Value<SAV, 2>::Type>::max();
    uint64_t maxLen = 0ul;
    for (auto const & s : seqs)
        if (length(s) > maxLen)
            maxLen = length(s);

    if (maxLen >= maxLenSeq)
    {
        std::cerr << "Too long sequences to be indexed:\n  "
                  << "length" << maxLen << " present in file, but only "
                  << maxLenSeq << " supported by index.\n";
        return false;
    }
    return true;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

inline int
convertMaskingFile(uint64_t numberOfSeqs,
                   LambdaIndexerOptions const & options)

{
    int ret = 0;
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntStarts;
    StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntEnds;
//     resize(segIntervals, numberOfSeqs, Exact());

    if (options.segFile != "")
    {
        std::cout << "Constructing binary seqan masking from seg-file...\n"
                << std::flush;

        std::ifstream stream;
        stream.open(toCString(options.segFile));
        if (!stream.is_open())
            return -1;

        typedef RecordReader<std::ifstream, SinglePass<> > TReader;
        TReader reader(stream);

//         StringSet<String<Tuple<unsigned, 2>>> _segIntervals;
//         auto & _segIntervals = segIntervals;
//         resize(_segIntervals, numberOfSeqs, Exact());
        StringSet<String<unsigned>> _segIntStarts;
        StringSet<String<unsigned>> _segIntEnds;
        resize(_segIntStarts, numberOfSeqs, Exact());
        resize(_segIntEnds, numberOfSeqs, Exact());
        CharString buf;
//         std::tuple<unsigned, unsigned> tup;

//         auto curSeq = begin(_segIntervals);
        unsigned curSeq = 0;
        while ((!atEnd(reader)) && (value(reader) == '>'))
        {
//             if (curSeq == end(_segIntervals))
//                 return -7;
            if (curSeq == numberOfSeqs)
                return -7;
            ret = skipLine(reader);
            if ((ret) && (ret != EOF_BEFORE_SUCCESS))
                return ret;

            unsigned curInt = 0;
            while ((!atEnd(reader)) && (value(reader) != '>'))
            {
                resize(_segIntStarts[curSeq], length(_segIntStarts[curSeq])+1);
                resize(_segIntEnds[curSeq], length(_segIntEnds[curSeq])+1);
                clear(buf);
                ret = readDigits(buf, reader);
                if (ret)
                    return ret;

//                 std::get<0>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntStarts[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);
                ret = skipNChars(reader, 3);
                if (ret)
                    return ret;

                clear(buf);
                ret = readDigits(buf, reader);
                if (ret)
                    return ret;

//                 std::get<1>(tup) = strtoumax(toCString(buf), 0, 10);
                _segIntEnds[curSeq][curInt] = strtoumax(toCString(buf), 0, 10);

//                 appendValue(*curSeq, tup);

                ret = skipLine(reader);
                if ((ret) && (ret != EOF_BEFORE_SUCCESS))
                    return ret;
                curInt++;
            }
            curSeq++;
        }
//         if (curSeq != end(_segIntervals))
//             return -9;
        if (curSeq != numberOfSeqs)
            return -9;

        segIntStarts.concat = concat(_segIntStarts);
        segIntStarts.limits = stringSetLimits(_segIntStarts);
        segIntEnds.concat = concat(_segIntEnds);
        segIntEnds.limits = stringSetLimits(_segIntEnds);
//         segIntEnds = _segIntEnds;
//         segIntervals = _segIntervals; // non-concatdirect to concatdirect

        stream.close();

    } else
    {
        std::cout << "No Seg-File specified, no masking will take place,\n"
                << std::flush;
//         resize(segIntervals, numberOfSeqs, Exact());
        resize(segIntStarts, numberOfSeqs, Exact());
        resize(segIntEnds, numberOfSeqs, Exact());
    }

//     for (unsigned u = 0; u < length(segIntStarts); ++u)
//     {
//         std::cout << u << ": ";
//         for (unsigned v = 0; v < length(segIntStarts[u]); ++v)
//         {
//             std::cout << '(' << segIntStarts[u][v] << ", " << segIntEnds[u][v] << ")  ";
//         }
//         std::cout << '\n';
//     }
    std::cout << "Dumping binary seqan mask file...\n"
                << std::flush;
    CharString _path = options.dbFile;
    append(_path, ".binseg_s");
    save(segIntStarts, toCString(_path));
    _path = options.dbFile;
    append(_path, ".binseg_e");
    save(segIntEnds, toCString(_path));
    std::cout << "Done.\n\n" << std::flush;

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TIndexSpec,typename TString, typename TSpec>
inline void
generateIndexAndDump(StringSet<TString, TSpec> & seqs,
                     LambdaIndexerOptions const & options)
{
    bool constexpr isFM = std::is_same<TIndexSpec, FMIndex<> >::value;
    using TDBIndex =  Index<StringSet<TString, TSpec>, TIndexSpec >;
    using TFibre = typename std::conditional<isFM,
                                             FibreSALF,
                                             FibreSA>::type;
    // Generate Index
    std::cout << "Generating Index..." << std::flush;
    double s = sysTime();

    if (isFM)
        reverse(seqs);

    TDBIndex dbIndex(seqs);
    indexRequire(dbIndex, TFibre());// instantiate

//     if (isFM)
//         reverse(seqs);

    double e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n\n" << std::flush;


    // Dump Index
    std::cout << "Writing Index to disk..." << std::flush;
    s = sysTime();
    save(dbIndex, toCString(options.dbFile));
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n" << std::flush;
}

#endif // header guard