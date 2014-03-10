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
// alph.hpp: Functions and Metafunctions for Alphabet reduction
// ==========================================================================

#ifndef SEQAN_LAMBDA_ALPH_H_
#define SEQAN_LAMBDA_ALPH_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

// using namespace seqan;

namespace seqan {

// ----------------------------------------------------------------------------
// Specialization AminoAcid10
// ----------------------------------------------------------------------------

/**
.Spec.AminoAcid10:
..cat:Alphabets
..summary:Iupac code for amino acids.
..general:Class.SimpleType
..signature:AminoAcid10
..remarks:
...text:The @Metafunction.ValueSize@ of $AminoAcid10$ is 24.
...text:The amino acids are enumerated from 0 to 15 in this order:
...text:'A'=0, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=19.
...text:The remaining 4 symbols are:
...text: 'B'=20 (Aspartic Acid, Asparagine), 'Z'=21 (Glutamic Acid, Glutamine), 'X'=22 (unknown), '*'=23 (terminator)
...text:Objects of type $AminoAcid10$ can be converted to $char$ and vice versa.
Unkown values are converted to $'X'$.
...text:$AminoAcid10$ is typedef for $SimpleType<char,AminoAcid10_>$, while $AminoAcid10_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/

struct AminoAcid10_ {};
typedef SimpleType<unsigned char, AminoAcid10_> AminoAcid10;

template <> struct ValueSize<AminoAcid10>
{
    typedef __uint8 Type;
    static const Type VALUE = 10;
};

template <> struct BitsPerValue<AminoAcid10>
{
    typedef __uint8 Type;
    static const Type VALUE = 4;
};

inline AminoAcid10
unknownValueImpl(AminoAcid10 *)
{
    static const AminoAcid10 _result = AminoAcid10('X');
    return _result;
}

// --------------------------------------------------------------------------
// Amino Acid
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableAA10ToAscii_
{
    static char const VALUE[10];
};

template <typename T>
char const TranslateTableAA10ToAscii_<T>::VALUE[10] =
{
    'A', // 0 Acidic or Proto-Acidic
         // 'N' Asn Asparagine, 'D' Asp Apartic Acid,  'B' Asn or Asp
         // 'Q' Gln Glutamine, 'E' Glu Glutamic Acid, 'Z' Glu or Gln
    'B', // 1 BASIC
         // 'R' Arg Arginine, 'K' Lys Lysine
    'C', // 2 Cys Cystine
    'H', // 3 His Histidine
    'S', // 4 SMALL UNPOLAR
         // 'A' Ala Alanine, 'G' Gly Glycine
    'U', // 5 UNPOLAR
         // 'I' Ile Isoleucine, 'L' Leu Leucine, 'V' Val Valine, 'M' Met Methionine
    'O', // 6 POLAR
         // 'S' Ser Serine, 'T' Thr Threonine
    'P', // 7 Pro Proline
    'R', // 8 AROMATIC
         //'F' Phe Phenylalanine, 'W' Trp Tryptophan, 'Y' Tyr Tyrosine
    'X', // 9 'X' Unknown, '*' Terminator
};

template <typename T = void>
struct TranslateTableAsciiToAA10_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableAsciiToAA10_<T>::VALUE[256] =
{
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //0
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //1
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //2
//                                                     *
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //3
     9,   4,   0,   2,   0,   0,   8,   4,   3,   5,   9,   1,   5,   5,   0,   9, //4
//    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,

     7,   0,   1,   6,   6,   9,   5,   8,   9,   8,   0,   9,   9,   9,   9,   9, //5
//   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,

     9,   4,   0,   2,   0,   0,   8,   4,   3,   5,   9,   1,   5,   5,   0,   9, //6
//    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

     7,   0,   1,   6,   6,   9,   5,   8,   9,   8,   0,   9,   9,   9,   9,   9, //7
//   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,

     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //8
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //9
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //10
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //11
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //12
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //13
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //14
     9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9  //15
};

template <typename T = void>
struct TranslateTableAAToAA10_
{
    static char const VALUE[24];
};

template <typename T>
char const TranslateTableAAToAA10_<T>::VALUE[24] =
{
    4,  1,  0,  0,  0,  2,  0,  4,  3,  5,  5,  1,  5,  8,  7,  6,  6,  8,  8,
//  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y
    5,  0,  0,  9,  9
//  V   B   Z   X   *
};


template <typename T = void>
struct TranslateTableByteToAA10_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToAA10_<T>::VALUE[256] =
{
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  9,  9,  9,  9,  9,  9, //0
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //1
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //2
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //3
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //4
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //5
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //6
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //7
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //8
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //9
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //10
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //11
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //12
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //13
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, //14
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9  //15
};



inline void assign(char & c_target, AminoAcid10 const & source)
{
    c_target = TranslateTableAA10ToAscii_<>::VALUE[source.value];
}


template <>
struct CompareType<AminoAcid10, __uint8>
{
    typedef AminoAcid10 Type;
};

inline void assign(AminoAcid10 & target, __uint8 c_source)
{
    target.value = TranslateTableByteToAA10_<>::VALUE[c_source];
}

template <>
struct CompareType<AminoAcid10, char>
{
    typedef AminoAcid10 Type;
};

inline void assign(AminoAcid10 & target, char c_source)
{
    target.value = TranslateTableAsciiToAA10_<>::VALUE[(unsigned char) c_source];
}

template <>
struct CompareType<AminoAcid10, AminoAcid>
{
    typedef AminoAcid10 Type;
};

inline void assign(AminoAcid10 & target, AminoAcid c_source)
{
    target.value = TranslateTableAAToAA10_<>::VALUE[c_source.value];
}


template <>
struct CompareType<AminoAcid10, Unicode>
{
    typedef AminoAcid10 Type;
};

inline void assign(AminoAcid10 & target, Unicode c_source)
{
    target.value = TranslateTableAsciiToAA10_<>::VALUE[(unsigned char) c_source];
}

} // namespace seqan

#endif // header guard