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
// trans.hpp: Functions and Metafunctions for translation
// ==========================================================================

#ifndef SEQAN_LAMBDA_TRANS_H_
#define SEQAN_LAMBDA_TRANS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
// #include <seqan/index_extras.h>

#include <seqan/blast.h>

#include "options.hpp"
#include "match.hpp"

using namespace seqan;

/*
char translationTable[4][4][4] =
{
    { // a??
        { 'K', 'N', 'K', 'N' }, // aa?
        { 'T', 'T', 'T', 'T' }, // ac?
        { 'R', 'S', 'R', 'S' }, // ag?
        { 'I', 'I', 'M', 'I' }  // au?
    },
    { // c??
        { 'Q', 'H', 'Q', 'H' }, // ca?
        { 'P', 'P', 'P', 'P' }, // cc?
        { 'R', 'R', 'R', 'R' }, // cg?
        { 'L', 'L', 'L', 'L' }  // cu?
    },
    { // g??
        { 'E', 'D', 'E', 'D' }, // ga?
        { 'A', 'A', 'A', 'A' }, // gc?
        { 'G', 'G', 'G', 'G' }, // gg?
        { 'V', 'V', 'V', 'V' }  // gu?
    },
    { // u??
        { '*', 'Y', '*', 'Y' }, // ua?
        { 'S', 'S', 'S', 'S' }, // uc?
        { '*', 'C', 'W', 'C' }, // ug?
        { 'L', 'F', 'L', 'F' }  // uu?
    }
};


// ----------------------------------------------------------------------------
// append6FrameTranslation()
// ----------------------------------------------------------------------------

template <typename TSString>
inline void
appendSixFrameTranslation(StringSet<String<AminoAcid>,
                                    Owner<ConcatDirect<>>> & target,
                          TSString                   const & source)
{
    Dna5String dna = source;
    for (unsigned char c : { 0, 1} )
    {
        (void)c;
        for (unsigned char cc : { 0, 1, 2 } )
        {
            String<AminoAcid> t;
//             reserve(back(target), (length(source) - cc) / 3, Exact());
            for (unsigned int i = cc; i < length(dna) - 3; i+=3)
            {
                if ((ordValue(value(dna, i  ) == 4)) ||
                    (ordValue(value(dna, i+1) == 4)) || // 'N' encountered
                    (ordValue(value(dna, i+2) == 4)) )
                    appendValue(t, 'X', Exact());
                else
                    appendValue(t,
                                translationTable[ordValue(value(dna, i  ))]
                                                [ordValue(value(dna, i+1))]
                                                [ordValue(value(dna, i+2))],
                                Exact());
            }
            appendValue(target, t, Exact());
        }
        reverseComplement(dna);
    }
}

template <typename TAAString, typename TDnaString>
inline void
translateInPlace(TAAString & target, TDnaString const & source)
{
//     SEQAN_ASSERT_EQ(length(target) -1, (length(source)-1)/3);
    for (unsigned int i = 0; i+2 < length(source); i+=3)
    {
        if ((ordValue(value(source, i  ) == 4)) ||
            (ordValue(value(source, i+1) == 4)) ||
            (ordValue(value(source, i+2) == 4)) )   // 'N' encountered
            target[i/3] = 'X';
        else
            target[i/3] = translationTable[ordValue(value(source, i  ))]
                                          [ordValue(value(source, i+1))]
                                          [ordValue(value(source, i+2))];
    }
}

template <typename TStringSet>
inline void
sixFrameTranslateImpl(StringSet<String<AminoAcid>,
                            Owner<ConcatDirect<>>> & target,
                      TStringSet             const & source)
{

    clear(target);
    resize(target.limits, length(source) * 6 + 1, Exact());
    target.limits[0] = 0;
    for (unsigned i = 0; i+1 < length(target.limits); ++i)
        target.limits[i+1] =    target.limits[i]
                            +   (length(source[(i)/6]) / 3)
                            - (((length(source[(i)/6]) % 3) < (i%3))
                                ? 1 : 0);
        //explenation:
        //                     start with last item's start pos
        //                  add  current dnastring's length / 3 (3DNA -> 1 AA)
        //                  shorten for shifted frames, if source has not
        //                  enough characters


    resize(target.concat, back(target.limits), Exact());

     #pragma omp parallel for schedule(guided)
    for (unsigned i = 0; i < length(target); ++i)
    {
        auto trgt = infix(target.concat, target.limits[i], target.limits[i+1]);

        if (i % 6 > 2)
        {
            String<typename Value<typename Value<TStringSet>::Type>::Type> src = source[i/6];
            reverseComplement(src);
            translateInPlace(trgt,
                             suffix(src, i % 3));
        } else
            translateInPlace(trgt,
                             suffix(source[i/6], i % 3));
    }
}

template <typename TSpec1, typename TSpec2>
inline void
sixFrameTranslate(StringSet<String<AminoAcid>,
                            Owner<ConcatDirect<>>> & target,
                  StringSet<String<Dna5, TSpec1>,
                            TSpec2>          const & source)
{
    return sixFrameTranslateImpl(target, source);
}

template <typename TSpec1, typename TSpec2>
inline void
sixFrameTranslate(StringSet<String<AminoAcid>,
                            Owner<ConcatDirect<>>> & target,
                  StringSet<String<Dna, TSpec1>,
                            TSpec2>          const & source)
{
    return sixFrameTranslateImpl(target, source);
}

inline void
sixFrameTranslate(StringSet<String<AminoAcid>,
                            Owner<ConcatDirect<>>>       & target,
                  StringSet<CharString,
                            Owner<ConcatDirect<>>> const & source)
{
    StringSet<String<Dna5>, Owner<ConcatDirect<>>> new_source;
    new_source.concat = source.concat; //implicit conversion to Dna5
    new_source.limits = source.limits;
    return sixFrameTranslateImpl(target, new_source);
}
*/

// ----------------------------------------------------------------------------
// getTrueQryId()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr TId
getTrueQryId(TId const & qryId,
             LambdaOptions const & /**/,
             BlastProgram const blastProgram)
{
    return qryId / qNumFrames(blastProgram);
}

// ----------------------------------------------------------------------------
// getTrueSubjId()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr TId
getTrueSubjId(TId                   const & subjId,
              LambdaOptions         const & /**/,
              BlastProgram const blastProgram)
{
    return subjId / sNumFrames(blastProgram);
}

// ----------------------------------------------------------------------------
// qryIsReverseComplemented()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr bool
qryIsReverseComplemented(TId const & qryId,
                         LambdaOptions         const & /**/,
                         BlastProgram const blastProgram)
{
    return (qNumFrames(blastProgram) == 2
                ? (qryId % 2 == 1)
                : (qNumFrames(blastProgram) == 6
                    ? (qryId % 6 > 2)
                    : false));
}

// ----------------------------------------------------------------------------
// subjIsReverseComplemented()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr bool
subjIsReverseComplemented(TId                   const & subjId,
                          LambdaOptions         const & /**/,
                          BlastProgram const blastProgram)
{
    return (sNumFrames(blastProgram) == 2
                ? (subjId % 2 == 1)
                : (sNumFrames(blastProgram) == 6
                    ? (subjId % 6 > 2)
                    : false));
}

// ----------------------------------------------------------------------------
// getQryFrameShift()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastFormatGeneration g>
constexpr unsigned char
getQryFrameShift(TId const & qryId,
             LambdaOptions const & /**/,
             BlastProgram const blastProgram)
{
    return (qIsTranslated(blastProgram) ? qryId % 3 : 0);
}

// ----------------------------------------------------------------------------
// getSubjFrameShift()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr unsigned char
getSubjFrameShift(TId                   const & /**/,
                  LambdaOptions         const & /**/,
                  BlastProgram const blastProgram)
{
    return (sIsTranslated(blastProgram) ? subjId % 3 : 0);
}


// ----------------------------------------------------------------------------
// getTrueQryStartPos()
// ----------------------------------------------------------------------------

// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastProgram p,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryStartPos(TId                                  const & /**/,
//                    TPos                                 const & qryStart,
//                    TPos                                 const & /**/,
//                    LambdaOptions                        const & /**/,
//                    BlastProgram const blastProgram)
// {
//     return qryStart + 1;
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryStartPos(TId                                  const & qryId,
//                    TPos                                 const & qryStart,
//                    TPos                                 const & qryEnd,
//                    LambdaOptions                        const & options,
//                    BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::BLASTN,g> TFormat;
//     return (qryIsReverseComplemented(qryId, options, TFormat()))
//                 ? qryEnd
//                 : qryStart + 1;
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryStartPos(TId                                  const & qryId,
//                    TPos                                 const & qryStart,
//                    TPos                                 const & qryEnd,
//                    LambdaOptions                        const & options,
//                    BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::BLASTX,g> TFormat;
//     return (qryIsReverseComplemented(qryId, options, TFormat()))
//             ? (qryEnd   * 3 + getQryFrameShift(qryId, options, TFormat()))
//             : (qryStart * 3 + getQryFrameShift(qryId, options, TFormat()) + 1);
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryStartPos(TId                                  const & qryId,
//                    TPos                                 const & qryStart,
//                    TPos                                 const & qryEnd,
//                    LambdaOptions                        const & options,
//                    BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTX,g> TFormat;
//     return (qryIsReverseComplemented(qryId, options, TFormat()))
//             ? (qryEnd   * 3 + getQryFrameShift(qryId, options, TFormat()))
//             : (qryStart * 3 + getQryFrameShift(qryId, options, TFormat()) + 1);
// }


// ----------------------------------------------------------------------------
// getTrueQryEndPos()
// ----------------------------------------------------------------------------

// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastProgram p,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryEndPos(TId                                  const & /**/,
//                  TPos                                 const & /**/,
//                  TPos                                 const & qryEnd,
//                  LambdaOptions                        const & /**/,
//                  BlastProgram const blastProgram)
// {
//     return qryEnd + 1;
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryEndPos(TId                                  const & qryId,
//                  TPos                                 const & qryStart,
//                  TPos                                 const & qryEnd,
//                  LambdaOptions                        const & options,
//                  BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::BLASTN,g> TFormat;
//     return (qryIsReverseComplemented(qryId, options, TFormat()))
//                 ? qryStart + 1
//                 : qryEnd;
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryEndPos(TId                                  const & qryId,
//                  TPos                                 const & qryStart,
//                  TPos                                 const & qryEnd,
//                  LambdaOptions                        const & options,
//                  BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::BLASTX,g> TFormat;
//     return (!qryIsReverseComplemented(qryId, options, TFormat()))
//             ? (qryEnd   * 3 + getQryFrameShift(qryId, options, TFormat()))
//             : (qryStart * 3 + getQryFrameShift(qryId, options, TFormat()) + 1);
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueQryEndPos(TId                                  const & qryId,
//                  TPos                                 const & qryStart,
//                  TPos                                 const & qryEnd,
//                  LambdaOptions const & options,
//                  BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTX,g> TFormat;
//     return (!qryIsReverseComplemented(qryId, options, TFormat()))
//             ? (qryEnd   * 3 + getQryFrameShift(qryId, options, TFormat()))
//             : (qryStart * 3 + getQryFrameShift(qryId, options, TFormat()) + 1);
// }

// ----------------------------------------------------------------------------
// getTrueSubjStartPos()
// ----------------------------------------------------------------------------

// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastProgram p,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjStartPos(TId                                       const & /**/,
//                     TPos                                      const & subjStart,
//                     TPos                                      const & /**/,
//                     LambdaOptions const & /**/,
//                     BlastProgram const blastProgram)
// {
//     return subjStart + 1;
// }
// 
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjStartPos(TId                                     const & subjId,
//                     TPos                                    const & subjStart,
//                     TPos                                    const & subjEnd,
//                     LambdaOptions                           const & options,
//                     BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTN,g> TFormat;
//     return (subjIsReverseComplemented(subjId, options, TFormat()))
//             ? (subjEnd   * 3 + getQryFrameShift(subjId, options, TFormat()))
//             : (subjStart * 3 + getQryFrameShift(subjId, options, TFormat()) + 1);
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjStartPos(TId                                     const & subjId,
//                     TPos                                    const & subjStart,
//                     TPos                                    const & subjEnd,
//                     LambdaOptions                           const & options,
//                     BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTX,g> TFormat;
//     return (subjIsReverseComplemented(subjId, options, TFormat()))
//             ? (subjEnd   * 3 + getQryFrameShift(subjId, options, TFormat()))
//             : (subjStart * 3 + getQryFrameShift(subjId, options, TFormat()) + 1);
// }
// 

// ----------------------------------------------------------------------------
// getTrueSubjEndPos()
// ----------------------------------------------------------------------------

// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastProgram p,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjEndPos(TId                                       const & /**/,
//                   TPos                                      const & /**/,
//                   TPos                                      const & subjEnd,
//                   LambdaOptions                             const & /**/,
//                   BlastProgram const blastProgram)
// {
//     return subjEnd;
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjEndPos(TId                                       const & subjId,
//                   TPos                                      const & subjStart,
//                   TPos                                      const & subjEnd,
//                   LambdaOptions                             const & options,
//                   BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTN,g> TFormat;
//     return (!subjIsReverseComplemented(subjId, options, TFormat()))
//             ? (subjEnd   * 3 + getQryFrameShift(subjId, options, TFormat()))
//             : (subjStart * 3 + getQryFrameShift(subjId, options, TFormat()) + 1);
// }
// 
// template <typename TId,
//           typename TPos,
//           BlastFormatFile mf,
//           BlastFormatGeneration g>
// constexpr TPos
// getTrueSubjEndPos(TId                                       const & subjId,
//                   TPos                                      const & subjStart,
//                   TPos                                      const & subjEnd,
//                   LambdaOptions                             const & options,
//                   BlastProgram const blastProgram)
// {
//     typedef BlastFormat<mf,BlastProgram::TBLASTX,g> TFormat;
//     return (!subjIsReverseComplemented(subjId, options, TFormat()))
//             ? (subjEnd   * 3 + getQryFrameShift(subjId, options, TFormat()))
//             : (subjStart * 3 + getQryFrameShift(subjId, options, TFormat()) + 1);
// }
// 

// ----------------------------------------------------------------------------
// getTrueSubjEndPos()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr bool
qryIsSameFrame(TId                                       const & qryId1,
               TId                                       const & qryId2,
               LambdaOptions                             const & options,
               BlastProgram const blastProgram)
{
    typedef BlastFormat<mf,p,g> TFormat;
    return ((getQryFrameShift(qryId1, options, TFormat()) ==
             getQryFrameShift(qryId2, options, TFormat())) &&
            (qryIsReverseComplemented(qryId1, options, TFormat()) ==
             qryIsReverseComplemented(qryId2, options, TFormat())));
}


// ----------------------------------------------------------------------------
// getTrueSubjEndPos()
// ----------------------------------------------------------------------------

template <typename TId,
          BlastFormatFile mf,
          BlastProgram p,
          BlastFormatGeneration g>
constexpr bool
subjIsSameFrame(TId                                       const & subjId1,
                TId                                       const & subjId2,
                LambdaOptions                             const & options,
                BlastProgram const blastProgram)
{
    typedef BlastFormat<mf,p,g> TFormat;
    return ((getSubjFrameShift(subjId1, options, TFormat()) ==
             getSubjFrameShift(subjId2, options, TFormat())) &&
            (subjIsReverseComplemented(subjId1, options, TFormat()) ==
             subjIsReverseComplemented(subjId2, options, TFormat())));
}



#endif // header guard