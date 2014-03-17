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
// match.h: Main File for the match class
// ==========================================================================

#ifndef SEQAN_LAMBDA_FINDER_H_
#define SEQAN_LAMBDA_FINDER_H_

#include <forward_list>
#include <vector>

#include "misc.hpp"



using namespace seqan;

//-----------------------------------------------------------------------------
//                  Finder
//-----------------------------------------------------------------------------

struct Match
{
    typedef uint32_t    TQId;
    typedef uint32_t    TSId;// many suffixes in subject-index
    typedef uint16_t    TPos;
    TQId qryId;
    TSId subjId;
    TPos qryStart;
//     TPos qryEnd;

    TPos subjStart;
//     TPos subjEnd;

//     Match()
//     :
//         qryId(0), qryStart(0), /*qryEnd(0),*/ subjId(0), subjStart(0)/*, subjEnd(0)*/
//     {
//     }
// 
//     Match(Match const & m2)
//     {
//         qryId       = m2.qryId;
//         qryStart    = m2.qryStart;
//         /*qryEnd      = m2.qryEnd;*/
//         subjId      = m2.subjId;
//         subjStart   = m2.subjStart;
//         /*subjEnd     = m2.subjEnd;*/
//     }

    inline bool operator== (Match const & m2) const
    {
         return std::tie(qryId, subjId, qryStart, subjStart/*, qryEnd, subjEnd*/)
             == std::tie(m2.qryId, m2.subjId, m2.qryStart, m2.subjStart/*, m2.qryEnd, m2.subjEnd*/);
    }
    inline bool operator< (Match const & m2) const
    {
         return std::tie(qryId, subjId, qryStart, subjStart/*, qryEnd, subjEnd*/)
              < std::tie(m2.qryId, m2.subjId, m2.qryStart, m2.subjStart/*, m2.qryEnd, m2.subjEnd*/);
    }
};

// inline bool
// overlap(Match const & m1, Match const & m2, unsigned char d = 0)
// {
//     SEQAN_ASSERT_EQ(m1.qryId,  m2.qryId);
//     SEQAN_ASSERT_EQ(m1.subjId, m2.subjId);
// 
// //     if (     //match beginnings overlap
// //         ((m1.qryStart  >= _protectUnderflow(m2.qryStart, d)) &&
// //          (m1.qryStart  <= _protectOverflow( m2.qryEnd,   d)) &&
// //          (m1.subjStart >= _protectUnderflow(m2.subjStart,d)) &&
// //          (m1.subjStart <= _protectOverflow( m2.subjEnd,  d))
// //         ) || // OR match endings overlap
// //         ((m1.qryEnd  >= _protectUnderflow(m2.qryStart, d)) &&
// //          (m1.qryEnd  <= _protectOverflow( m2.qryEnd,   d)) &&
// //          (m1.subjEnd >= _protectUnderflow(m2.subjStart,d)) &&
// //          (m1.subjEnd <= _protectOverflow( m2.subjEnd,  d))
// //         )
// //        )
// //         return true;
// 
//     //DEBUG DEBUG DEBUG TODO TODO TODO WHEN ALIGNMENT IS FIXED
//     // also add check if diffference of overlaps is  <= maximum number of gaps
//     if ((intervalOverlap(m1.qryStart, m1.qryEnd, m2.qryStart, m2.qryEnd) > -d) &&
//         (intervalOverlap(m1.subjStart, m1.subjEnd, m2.subjStart, m2.subjEnd) > -d) &&
//         ((m1.qryStart < m2.qryStart ) == (m1.subjStart < m2.subjStart))) // same order
//         return true;
// 
//     return false;
// }

// inline bool contains(Match const & m1, Match const & m2)
// {
// //     if (m1.qryId != m2.qryId)
// //         return false;
//     SEQAN_ASSERT_EQ(m1.qryId,  m2.qryId);
//     if (m1.qryStart != m2.qryStart)
//         return false;
// 
//     if (m1.qryStart > m2.qryStart)
//         return false;
//     if (m1.subjStart > m2.subjStart)
//         return false;
//     if (m1.qryEnd < m2.qryEnd)
//         return false;
//     if (m1.subjEnd < m2.subjEnd)
//         return false;
//     return true;
// }
// 
// inline void
// mergeUnto(Match & m1, Match const & m2)
// {
//     SEQAN_ASSERT_EQ(m1.qryId,  m2.qryId);
//     SEQAN_ASSERT_EQ(m1.subjId, m2.subjId);
//     m1.qryStart  = std::min(m1.qryStart, m2.qryStart);
//     m1.qryEnd    = std::max(m1.qryEnd,   m2.qryEnd);
//     m1.subjStart = std::min(m1.subjStart,m2.subjStart);
//     m1.subjEnd   = std::max(m1.subjEnd,  m2.subjEnd);
// }


// inline void
// _printMatch(Match const & m)
// {
//     std::cout << "MATCH  Query " << m.qryId
//               << "(" << m.qryStart << ", " << m.qryEnd
//               << ")   on Subject "<< m.subjId
//               << "(" << m.subjStart << ", " << m.subjEnd
//               << ")" <<  std::endl << std::flush;
// }





#endif // header guard
