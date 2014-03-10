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

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return (options.revComp) ? (m.qryId / 2) : m.qryId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId / 6;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId / 6;
}

// ----------------------------------------------------------------------------
// getTrueSubjId()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match                 const & m,
              LambdaOptions        const & /**/,
              BlastFormat<mf,p,g>   const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match                                         const & m,
              LambdaOptions                                const & /**/,
              BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId / 6;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match                                         const & m,
              LambdaOptions                                const & /**/,
              BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId / 6;
}

// ----------------------------------------------------------------------------
// qryIsReverseComplemented()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr bool
qryIsReverseComplemented(Match                  const & /**/,
                         LambdaOptions         const & /**/,
                         BlastFormat<mf,p,g>    const & /*tag*/)
{
    return false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match                  const & m,
                         LambdaOptions         const & options,
                         BlastFormat<mf,
                                     BlastFormatOptions::BlastN,
                                     g>         const & /*tag*/)
{
    return (options.revComp) ? (m.qryId % 2 == 0) : false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match                  const & m,
                         LambdaOptions         const & /**/,
                         BlastFormat<mf,
                                     BlastFormatOptions::BlastX,
                                     g>         const & /*tag*/)
{
    return (m.qryId % 6 > 2);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match                  const & m,
                         LambdaOptions         const & /**/,
                         BlastFormat<mf,
                                  BlastFormatOptions::TBlastX,
                                  g>            const & /*tag*/)
{
    return (m.qryId % 6 > 2);
}

// ----------------------------------------------------------------------------
// subjIsReverseComplemented()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr bool
subjIsReverseComplemented(Match const & m,
                      LambdaOptions const & /**/,
                      BlastFormat<mf,p,g> const & /*tag*/)
{
    return false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
subjIsReverseComplemented(Match const & m,
                      LambdaOptions const & /**/,
                      BlastFormat<mf,
                                  BlastFormatOptions::TBlastN,
                                  g> const & /*tag*/)
{
    return (m.subjId % 6 > 2);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
subjIsReverseComplemented(Match const & m,
                      LambdaOptions const & /**/,
                      BlastFormat<mf,
                                  BlastFormatOptions::TBlastX,
                                  g> const & /*tag*/)
{
    return (m.subjId % 6 > 2);
}

// ----------------------------------------------------------------------------
// getQryFrameShift()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & /**/,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
{
    return 0;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
{
    return m.qryId % 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
{
    return m.qryId % 3;
}

// ----------------------------------------------------------------------------
// getSubjFrameShift()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & /**/,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
{
    return 0;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
{
    return m.subjId % 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
{
    return m.subjId % 3;
}


// ----------------------------------------------------------------------------
// getTrueQryStartPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    return m.qryStart + 1;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastN,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
                ? m.qryEnd
                : m.qryStart + 1;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastX,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}


// ----------------------------------------------------------------------------
// getTrueQryEndPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    return m.qryEnd + 1;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastN,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
                ? m.qryStart + 1
                : m.qryEnd;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastX,g> TFormat;
    return (!qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (!qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}

// ----------------------------------------------------------------------------
// getTrueSubjStartPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    return m.subjStart + 1;
}


template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastN,g> TFormat;
    return (subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}


// ----------------------------------------------------------------------------
// getTrueSubjEndPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             LambdaOptions const & /**/,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    return m.subjEnd;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastN,g> TFormat;
    return (!subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             LambdaOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (!subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 + getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()) + 1);
}



#endif // header guard