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

template <typename TString, typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step03_generateIndexAndDump(StringSet<TString, TSpec> & seqs,
                            LambdaIndexerOptions const & options,
                            BlastFormat<m, p, g> const & /*tag*/)
{
    typedef Index<StringSet<TString, TSpec>, IndexSa<> > TDbIndex;

    // Generate Index
    std::cout << "Generating Index..." << std::flush;
    double s = sysTime();
    TDbIndex dbIndex(seqs);
    typename Iterator<TDbIndex, TopDown<> >::Type it(dbIndex); // instantiate
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

    return 0;
}

template <typename TAlph, typename TSpec1, typename TSpec2,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step02_reduceAlphabet(StringSet<String<TAlph, TSpec1>, TSpec2> & oldSubj,
                      LambdaIndexerOptions const & options,
                      BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    switch (options.alphReduction)
    {
        case 1:
        {
            typedef AminoAcid10 TAA;
            StringSet<String<TAA, TSpec1>,Owner<ConcatDirect<>>> newSubj;
            newSubj.concat = oldSubj.concat; //implicit conversion
            newSubj.limits = oldSubj.limits;
            return step03_generateIndexAndDump(newSubj, options, TFormat());
        } break;
        case 2:
        {
            typedef ReducedAminoAcid<Murphy10> TAA;
            StringSet<String<TAA, TSpec1>,Owner<ConcatDirect<>>> newSubj;
            newSubj.concat = oldSubj.concat; //implicit conversion
            newSubj.limits = oldSubj.limits;
            return step03_generateIndexAndDump(newSubj, options, TFormat());
        } break;
        case 8:
        {
            typedef ReducedAminoAcid<ClusterReduction<8>> TAA;
            StringSet<String<TAA, TSpec1>,Owner<ConcatDirect<>>> newSubj;
            newSubj.concat = oldSubj.concat; //implicit conversion
            newSubj.limits = oldSubj.limits;
            return step03_generateIndexAndDump(newSubj, options, TFormat());
        } break;

        case 10:
        {
            typedef ReducedAminoAcid<ClusterReduction<10>> TAA;
            StringSet<String<TAA, TSpec1>,Owner<ConcatDirect<>>> newSubj;
            newSubj.concat = oldSubj.concat; //implicit conversion
            newSubj.limits = oldSubj.limits;
            return step03_generateIndexAndDump(newSubj, options, TFormat());
        } break;

        case 12:
        {
            typedef ReducedAminoAcid<ClusterReduction<12>> TAA;
            StringSet<String<TAA, TSpec1>,Owner<ConcatDirect<>>> newSubj;
            newSubj.concat = oldSubj.concat; //implicit conversion
            newSubj.limits = oldSubj.limits;
            return step03_generateIndexAndDump(newSubj, options, TFormat());
        } break;

        default:
            break;
    }

    return step03_generateIndexAndDump(oldSubj, options, TFormat());
}

// --------------------------------------------------------------------------
// Function preprocessSubjSeqs()
// --------------------------------------------------------------------------




template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                   TOldSubj & oldSubj,
                                   LambdaIndexerOptions const & /**/,
                                   BlastFormat<m,
                                               BlastFormatOptions::BlastN,
                                               g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    LambdaIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastP,
                                                g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;

}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    LambdaIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastX,
                                                g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    LambdaIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastN,
                                                g> const & /*tag*/)
{
    translate(newSubj, oldSubj, SIX_FRAME);
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    LambdaIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastX,
                                                g> const & /*tag*/)
{
    translate(newSubj, oldSubj, SIX_FRAME);
}

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step01_preprocessSubjSeqs(
                       StringSet<CharString, Owner<ConcatDirect<> > >  & oldSubj,
                       LambdaIndexerOptions       const & options,
                       BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    typename UnreducedStringSet<p>::Type newSubj;

    double start = sysTime();
    std::cout << "Preprocessing Subj Sequences..." << std::flush;

    // depending on blastProgram
    step01a_preprocessSubj(newSubj, oldSubj, options, TFormat());

    clear(oldSubj);
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Number of queries after preproc: " << length(newSubj)
                << "\n\n" << std::flush;

    if (options.alphReduction > 0)
    {
        start = sysTime();
        std::cout << "Dumping unreduced Subj Sequences..." << std::flush;

        //TODO save to TMPDIR instead
        CharString _path = options.dbFile;
        append(_path, ".unredsubj");
        save(newSubj, toCString(_path));

        std::cout << " done.\n";
        finish = sysTime() - start;
        std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

        return step02_reduceAlphabet(newSubj, options, TFormat());
    }

    return step03_generateIndexAndDump(newSubj, options, TFormat());
}

// --------------------------------------------------------------------------
// Function loadSubjSequences()
// --------------------------------------------------------------------------



template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step00_loadSubjSequences(LambdaIndexerOptions const & options,
                         BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    // load sequences as CharString first to be error tolerant
    StringSet<CharString, Owner<ConcatDirect<> > > seqs;
    StringSet<CharString, Owner<ConcatDirect<> > > ids;
    double start = sysTime();
    std::cout << "Loading Database Sequences..." << std::flush;
    int ret = 0;
    if (options.fileFormat == 1)
        ret = loadSeqsAndIds(ids, seqs, options.dbFile, Fastq());
    else
        ret = loadSeqsAndIds(ids, seqs, options.dbFile, Fasta());
    if (ret)
        return ret;
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    unsigned long maxLen = 0ul;
    for (auto const & s : seqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(seqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    start = sysTime();
    std::cout << "Dumping Subj Ids..." << std::flush;

    //TODO save to TMPDIR instead
    CharString _path = options.dbFile;
    append(_path, ".ids");
    save(ids, toCString(_path));

    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    // MASKING
//     String<std::forward_list<std::tuple<unsigned, unsigned>>> segIntervals;
//     StringSet<String<Tuple<unsigned, 2>>, Owner<ConcatDirect<>>> segIntervals;
    {
        StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntStarts;
        StringSet<String<unsigned>, Owner<ConcatDirect<>>> segIntEnds;
    //     resize(segIntervals, length(seqs));

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
    //         resize(_segIntervals, length(seqs));
            StringSet<String<unsigned>> _segIntStarts;
            StringSet<String<unsigned>> _segIntEnds;
            resize(_segIntStarts, length(seqs));
            resize(_segIntEnds, length(seqs));
            CharString buf;
    //         std::tuple<unsigned, unsigned> tup;

    //         auto curSeq = begin(_segIntervals);
            unsigned curSeq = 0;
            while ((!atEnd(reader)) && (value(reader) == '>'))
            {
    //             if (curSeq == end(_segIntervals))
    //                 return -7;
                if (curSeq == length(seqs))
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
            if (curSeq != length(seqs))
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
    //         resize(segIntervals, length(seqs));
            resize(segIntStarts, length(seqs));
            resize(segIntEnds, length(seqs));
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
        _path = options.dbFile;
        append(_path, ".binseg_s");
        save(segIntStarts, toCString(_path));
        _path = options.dbFile;
        append(_path, ".binseg_e");
        save(segIntEnds, toCString(_path));
        std::cout << "Done.\n\n" << std::flush;
    }
    return step01_preprocessSubjSeqs(seqs, options, TFormat());
}

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
beginPipeline(LambdaIndexerOptions const & options,
              BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;
    return step00_loadSubjSequences(options, TFormat());
}

#endif // header guard