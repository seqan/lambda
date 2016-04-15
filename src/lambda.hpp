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
// lambda.hpp: contains the main progam pipeline
// ==========================================================================


#ifndef SEQAN_LAMBDA_LAMBDA_H_
#define SEQAN_LAMBDA_LAMBDA_H_

#include <type_traits>
#include <iomanip>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/misc/terminal.h>

#include <seqan/index.h>

#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include <seqan/align_extend.h>

#include "options.hpp"
#include "match.hpp"
#include "misc.hpp"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

template <typename TRedAlph_,
          typename TIndexSpec_,
          typename TFileFormat,
          BlastProgram p,
          BlastTabularSpec h>
class GlobalDataHolder;

template <typename TGlobalHolder_,
          typename TScoreExtension>
class LocalDataHolder;

// ============================================================================
// Classes, structs, enums
// ============================================================================

enum COMPUTERESULT_
{
    SUCCESS = 0,
    PREEXTEND,
    PERCENTIDENT,
    EVALUE,
    OTHER_FAIL
};

//TODO replace with lambda
// comparison operator to sort SA-Values based on the strings in the SA they refer to
template <typename TSav, typename TStringSet>
struct Comp :
    public ::std::binary_function < TSav, TSav, bool >
{
    TStringSet const & stringSet;

    Comp(TStringSet const & _stringSet)
        : stringSet(_stringSet)
    {}

    inline bool operator() (TSav const & i, TSav const & j) const
    {
         return (value(stringSet,getSeqNo(i)) < value(stringSet,getSeqNo(j)));
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function prepareScoring()
// --------------------------------------------------------------------------

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline void
prepareScoringMore(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>  & globalHolder,
                   LambdaOptions                                                    const & options,
                   std::true_type                                                   const & /**/)
{
    setScoreMatch(context(globalHolder.outfile).scoringScheme, options.match);
    setScoreMismatch(context(globalHolder.outfile).scoringScheme, options.misMatch);
}

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline void
prepareScoringMore(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h> & globalHolder,
                   LambdaOptions                                                    const & options,
                   std::false_type                                                  const & /**/)
{
    switch (options.scoringMethod)
    {
//         case 0:
//             return argConv3(options, TOutFormat(), Th(), Tp(), TRedAlph(), Score<int, Simple>());
        case 45:
            setScoreMatrixById(context(globalHolder.outfile).scoringScheme._internalScheme,
                               AminoAcidScoreMatrixID::BLOSUM45);
            break;
        case 62:
            setScoreMatrixById(context(globalHolder.outfile).scoringScheme._internalScheme,
                               AminoAcidScoreMatrixID::BLOSUM62);
            break;
        case 80:
            setScoreMatrixById(context(globalHolder.outfile).scoringScheme._internalScheme,
                               AminoAcidScoreMatrixID::BLOSUM80);
            break;
        default:
            break;
    }
}

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline int
prepareScoring(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h> & globalHolder,
               LambdaOptions                                                    const & options)
{
    using TGlobalHolder = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;
    prepareScoringMore(globalHolder,
                       options,
                       std::is_same<typename TGlobalHolder::TScoreScheme, Score<int, Simple>>());

    setScoreGapOpenBlast(context(globalHolder.outfile).scoringScheme, options.gapOpen);
    setScoreGapExtend(context(globalHolder.outfile).scoringScheme, options.gapExtend);

    if (!isValid(context(globalHolder.outfile).scoringScheme))
    {
        std::cerr << "Could not computer Karlin-Altschul-Values for "
                  << "Scoring Scheme. Exiting.\n";
        return -1;
    }
    return 0;
}

// --------------------------------------------------------------------------
// Function loadSubjects()
// --------------------------------------------------------------------------

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline int
loadSubjects(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h> & globalHolder,
             LambdaOptions                                                    const & options)
{
    using TGH = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;

    double start, finish;
    std::string strIdent;
    int ret;
    CharString _dbSeqs;

    if (TGH::indexIsFM || TGH::alphReduction) // otherwise sequences are loaded as part of index
    {
        start = sysTime();
        strIdent = "Loading Subj Sequences...";
        myPrint(options, 1, strIdent);

        _dbSeqs = options.dbFile;
        append(_dbSeqs, ".");
        append(_dbSeqs, _alphName(TransAlph<p>()));

        ret = open(globalHolder.subjSeqs, toCString(_dbSeqs), OPEN_RDONLY);
        if (ret != true)
        {
            std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                        << " failed.\n";
            return 1;
        }

        if (TGH::alphReduction)
            globalHolder.redSubjSeqs.limits = globalHolder.subjSeqs.limits;

        finish = sysTime() - start;
        myPrint(options, 1, " done.\n");
        myPrint(options, 2, "Runtime: ", finish, "s \n", "Amount: ",
                length(globalHolder.subjSeqs), "\n\n");
    }

    start = sysTime();
    strIdent = "Loading Subj Ids...";
    myPrint(options, 1, strIdent);

    _dbSeqs = options.dbFile;
    append(_dbSeqs, ".ids");
    ret = open(globalHolder.subjIds, toCString(_dbSeqs), OPEN_RDONLY);
    if (ret != true)
    {
        std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                  << " failed.\n";
        return 1;
    }
    finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    context(globalHolder.outfile).dbName = options.dbFile;

    // if subjects where translated, we don't have the untranslated seqs at all
    // but we still need the data for statistics and position un-translation
    if (sIsTranslated(p))
    {
        start = sysTime();
        std::string strIdent = "Loading Lengths of untranslated Subj sequences...";
        myPrint(options, 1, strIdent);

        _dbSeqs = options.dbFile;
        append(_dbSeqs, ".untranslengths");
        ret = open(globalHolder.untransSubjSeqLengths, toCString(_dbSeqs), OPEN_RDONLY);
        if (ret != true)
        {
            std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                      << " failed.\n";
            return 1;
        }

        finish = sysTime() - start;
        myPrint(options, 1, " done.\n");
        myPrint(options, 2, "Runtime: ", finish, "s \n\n");
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <typename TGlobalHolder>
inline int
loadDbIndexFromDisk(TGlobalHolder       & globalHolder,
                    LambdaOptions const & options)
{
    std::string strIdent = "Loading Database Index...";
    myPrint(options, 1, strIdent);
    double start = sysTime();
    std::string path = toCString(options.dbFile);
    path += '.' + std::string(_alphName(typename TGlobalHolder::TRedAlph()));
    if (TGlobalHolder::indexIsFM)
        path += ".fm";
    else
        path += ".sa";

    // Check if the index is of the old format (pre 0.9.0) by looking for different files
    if ((globalHolder.blastProgram != BlastProgram::BLASTN) && // BLASTN indexes are compatible
        ((TGlobalHolder::alphReduction && fileExists(toCString(path + ".txt.concat"))) ||
        (!TGlobalHolder::alphReduction && TGlobalHolder::indexIsFM && !fileExists(toCString(path + ".lf.drv.wtc.24")))))
    {
        std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                  << " failed.\n"
                  << "It appears you tried to open an old index (created before 0.9.0) which "
                  << "is not supported. Please remove the old files and create a new index with lambda_indexer!\n";
        return 200;
    }

    int ret = open(globalHolder.dbIndex, path.c_str(), OPEN_RDONLY);
    if (ret != true)
    {
        std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                  << " failed. "
                  << "Did you use the same options as with lambda_indexer?\n";
        return 1;
    }

    // assign previously loaded sub sequences (possibly modifier-wrapped
    // to the text-member of our new index (unless isFM, which doesnt need text)
    if ((!TGlobalHolder::indexIsFM) && (TGlobalHolder::alphReduction))
        indexText(globalHolder.dbIndex) = globalHolder.redSubjSeqs;

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n", "No of Fibres: ",
            length(indexSA(globalHolder.dbIndex)), "\n\n");

    // this is actually part of prepareScoring(), but the values are just available now
    if (sIsTranslated(globalHolder.blastProgram ))
    {
        // last value has sum of lengths
        context(globalHolder.outfile).dbTotalLength  = back(globalHolder.untransSubjSeqLengths);
        context(globalHolder.outfile).dbNumberOfSeqs = length(globalHolder.untransSubjSeqLengths) - 1;
    } else
    {
        context(globalHolder.outfile).dbTotalLength  = length(concat(globalHolder.subjSeqs));
        context(globalHolder.outfile).dbNumberOfSeqs = length(globalHolder.subjSeqs);
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function loadSegintervals()
// --------------------------------------------------------------------------

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline int
loadSegintervals(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>     & globalHolder,
                 LambdaOptions                                            const & options)
{

    double start = sysTime();
    std::string strIdent = "Loading Database Masking file...";
    myPrint(options, 1, strIdent);

    CharString segFileS = options.dbFile;
    append(segFileS, ".binseg_s.concat");
    CharString segFileE = options.dbFile;
    append(segFileE, ".binseg_e.concat");
    bool fail = false;
    struct stat buffer;
    // file exists
    if ((stat(toCString(segFileS), &buffer) == 0) &&
        (stat(toCString(segFileE), &buffer) == 0))
    {
        //cut off ".concat" again
        resize(segFileS, length(segFileS) - 7);
        resize(segFileE, length(segFileE) - 7);

        fail = !open(globalHolder.segIntStarts, toCString(segFileS), OPEN_RDONLY);
        if (!fail)
            fail = !open(globalHolder.segIntEnds, toCString(segFileE), OPEN_RDONLY);
    } else
    {
        fail = true;
    }

    if (fail)
    {
        std::cerr << ((options.verbosity == 0) ? strIdent : std::string())
                  << " failed.\n";
        return 1;
    }

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
    return 0;
}

// --------------------------------------------------------------------------
// Function loadQuery()
// --------------------------------------------------------------------------

// BLASTX, TBLASTX
template <typename TSourceAlph, typename TSpec1,
          typename TTargetAlph, typename TSpec2,
          typename TUntransLengths,
          MyEnableIf<!std::is_same<TSourceAlph, TTargetAlph>::value> = 0>
inline void
loadQueryImplTrans(TCDStringSet<String<TTargetAlph, TSpec1>> & target,
                   TCDStringSet<String<TSourceAlph, TSpec2>> & source,
                   TUntransLengths                            & untransQrySeqLengths,
                   LambdaOptions                        const & options)
{
    myPrint(options, 1, "translating...");
    // translate
    translate(target,
              source,
              SIX_FRAME,
              options.geneticCode);

    // preserve lengths of untranslated sequences
    resize(untransQrySeqLengths,
           length(source.limits),
           Exact());

#ifdef __clang__
    SEQAN_OMP_PRAGMA(parallel for)
#else
    SEQAN_OMP_PRAGMA(parallel for simd)
#endif
    for (uint32_t i = 0; i < (length(untransQrySeqLengths) - 1); ++i)
        untransQrySeqLengths[i] = source.limits[i + 1] - source.limits[i];

    // save sum of lengths (both strings have n + 1 elements
    back(source.limits) = length(source.concat);
}

// BLASTN
template <typename TSpec1,
          typename TSpec2,
          typename TUntransLengths>
inline void
loadQueryImplTrans(TCDStringSet<String<TransAlph<BlastProgram::BLASTN>, TSpec1>> & target,
                   TCDStringSet<String<TransAlph<BlastProgram::BLASTN>, TSpec2>> & source,
                   TUntransLengths                                                & /**/,
                   LambdaOptions                                            const & options)
{
    using TAlph = TransAlph<BlastProgram::BLASTN>;
//     using TReverseCompl =  ModifiedString<ModifiedString<String<TAlph>,
//                             ModView<FunctorComplement<TAlph>>>, ModReverse>;
    myPrint(options, 1, " generating reverse complements...");
    // no need for translation, but we need reverse complements
    resize(target.concat, length(source.concat) * 2);
    resize(target.limits, length(source) * 2 + 1);

    target.limits[0] = 0;
    uint64_t const l = length(target.limits) - 1;
    for (uint64_t i = 1; i < l; i+=2)
    {
        target.limits[i] = target.limits[i-1] + length(source[i/2]);
        target.limits[i+1] = target.limits[i] + length(source[i/2]);
    }

    FunctorComplement<TAlph> functor;
    uint64_t tBeg, tBegNext, len, sBeg;
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic) private(tBeg, tBegNext, len, sBeg))
    for (uint64_t i = 0; i < (l - 1); i+=2)
    {
        tBeg = target.limits[i];
        tBegNext = target.limits[i+1];
        len = tBegNext - tBeg;
        sBeg = source.limits[i/2];
        // avoid senseless copying and iterate manually
        for (uint32_t j = 0; j < len; ++j)
        {
            target.concat[tBeg + j] = source.concat[sBeg + j];
            target.concat[tBegNext + j] = functor(source.concat[sBeg+len-j-1]);
        }
    }
}

// BLASTP, TBLASTN
template <typename TSpec1,
          typename TSpec2,
          typename TUntransLengths>
inline void
loadQueryImplTrans(TCDStringSet<String<TransAlph<BlastProgram::BLASTP>, TSpec1>> & target,
                   TCDStringSet<String<TransAlph<BlastProgram::BLASTP>, TSpec2>> & source,
                   TUntransLengths           & /**/,
                   LambdaOptions       const & /**/)
{
    // no need for translation, but sequences have to be in right place
    swap(target, source);
}

template <BlastTabularSpec h,
          BlastProgram p,
          typename TRedAlph,
          typename TIndexSpec,
          typename TOutFormat>
inline int
loadQuery(GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h> & globalHolder,
          LambdaOptions                                                    const & options)
{
    using TGH = GlobalDataHolder<TRedAlph, TIndexSpec, TOutFormat, p, h>;
    double start = sysTime();

    std::string strIdent = "Loading Query Sequences and Ids...";
    myPrint(options, 1, strIdent);

    TCDStringSet<String<OrigQryAlph<p>, typename TGH::TQryTag>> origSeqs;

//     std::cout << "FOO " <<  toCString(options.queryFile) << " BAR" << std::endl;
    try
    {
        SeqFileIn infile(toCString(options.queryFile));
        int ret = myReadRecords(globalHolder.qryIds, origSeqs, infile);
        if (ret)
            return ret;
    }
    catch(IOError const & e)
    {
        std::cerr << "\nIOError thrown: " << e.what() << '\n'
                  << "Could not read the query file, make sure it exists and is readable.\n";
        return -1;
    }

    // translate
    loadQueryImplTrans(globalHolder.qrySeqs,
                       origSeqs,
                       globalHolder.untransQrySeqLengths,
                       options);

    // sam and bam need original sequences if translation happened
    if (qIsTranslated(globalHolder.blastProgram) && (options.outFileFormat > 0) &&
        (options.samBamSeq > 0))
        std::swap(origSeqs, globalHolder.untranslatedQrySeqs);

    if (TGH::alphReduction)
        globalHolder.redQrySeqs.limits = globalHolder.qrySeqs.limits;

    double finish = sysTime() - start;
    myPrint(options, 1, " done.\n");

    unsigned long maxLen = 0ul;
    for (auto const & s : globalHolder.qrySeqs)
        if (length(s) > maxLen)
            maxLen = length(s);

    myPrint(options, 2, "Runtime: ", finish, "s \n",
            "Number of effective query sequences: ",
            length(globalHolder.qrySeqs), "\nLongest query sequence: ",
            maxLen, "\n\n");

    if (length(globalHolder.qrySeqs) >= std::numeric_limits<typename TGH::TMatch::TQId>::max())
    {
        std::cerr << "ERROR: Too many sequences submitted. The maximum (including frames) is "
                  << std::numeric_limits<typename TGH::TMatch::TQId>::max()
                  << ".\n";
        return -1;
    }

    if (maxLen >= std::numeric_limits<typename TGH::TMatch::TPos>::max())
    {
        std::cerr << "ERROR: one or more of your query sequences are too long. "
                  << "The maximum length is " << std::numeric_limits<typename TGH::TMatch::TPos>::max()
                  << ".\n";
        return -1;
    }
    return 0;
}

/// THREAD LOCAL STUFF

// --------------------------------------------------------------------------
// Function generateSeeds()
// --------------------------------------------------------------------------
/*
#define THREADLINE std::cout << "\0338" << std::endl << "Thread " << lH.i \
<< std::endl; \
for (unsigned char i=0; i< lH.i*20; ++i) std::cout << std::endl;*/

template <typename TLocalHolder>
inline int
generateSeeds(TLocalHolder & lH)
{
    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, "Block ", std::setw(4), 
                       lH.i, ": Generating Seeds...");
        if (lH.options.isTerm)
            myPrint(lH.options, 1, lH.statusStr);
    }

    double start = sysTime();
    for (unsigned long i = lH.indexBeginQry; i < lH.indexEndQry; ++i)
    {
        for (unsigned j = 0;
             (j* lH.options.seedOffset + lH.options.seedLength)
                <= length(value(lH.gH.redQrySeqs, i));
             ++j)
        {
            appendValue(lH.seeds, infix(value(lH.gH.redQrySeqs, i),
                                     j* lH.options.seedOffset,
                                     j* lH.options.seedOffset
                                     + lH.options.seedLength),
                        Generous());
            appendValue(lH.seedRefs,  i, Generous());
            appendValue(lH.seedRanks, j, Generous());

//             std::cout << "seed: " << back(lH.seeds) << "\n";
        }
    }
    double finish = sysTime() - start;
    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, " done. ");
        appendToStatus(lH.statusStr, lH.options, 2, finish, "s. ",
                       length(lH.seeds), " seeds created.");

        myPrint(lH.options, 1, lH.statusStr);
    }
    return 0;
}


// --------------------------------------------------------------------------
// Function generateTrieOverSeeds()
// --------------------------------------------------------------------------


template <typename TLocalHolder>
inline int
generateTrieOverSeeds(TLocalHolder & lH)
{
    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, "Generating Query-Index...");
        if (lH.options.isTerm)
            myPrint(lH.options, 1, lH.statusStr);
    }

    double start = sysTime();
    // we only want full length seed sequences in index, build up manually
    typedef typename Fibre<typename TLocalHolder::TSeedIndex, EsaSA>::Type TSa;
    //TODO maybe swap here instead
    lH.seedIndex = decltype(lH.seedIndex)(lH.seeds);
    TSa & sa = indexSA(lH.seedIndex);
    resize(sa, length(lH.seeds));
    for (unsigned u = 0; u < length(lH.seeds); ++u)
    {
        assignValueI1(value(sa,u), u);
        assignValueI2(value(sa,u), 0);
    }
    Comp<typename Value<TSa>::Type, typename TLocalHolder::TSeeds const> comp(lH.seeds);
    std::sort(begin(sa, Standard()), end(sa, Standard()), comp);
    typename Iterator<typename TLocalHolder::TSeedIndex, TopDown<> >::Type it(lH.seedIndex); // instantiate
    double finish = sysTime() - start;

    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, " done. ");
        appendToStatus(lH.statusStr, lH.options, 2, finish, "s. ",
                       length(sa), " fibres in SeedIndex. ");;
        myPrint(lH.options, 1, lH.statusStr);
    }

    return 0;
}

// perform a fast local alignment score calculation on the seed and see if we
// reach above threshold
// WARNING the following function only works for hammingdistanced seeds
template <typename TGlobalHolder,
          typename TScoreExtension>
inline bool
seedLooksPromising(LocalDataHolder<TGlobalHolder, TScoreExtension> const & lH,
                   typename TGlobalHolder::TMatch const & m)
{
    // no pre-scoring, but still filter out XXX and NNN hits
//     if (!lH.options.preScoring))
//     {
//         for (unsigned i = m.qryStart, count = 0; i < m.qryStart + lH.options.seedLength; ++i)
//             if (lH.gH.qrySeqs[m.qryId][i] == unkownValue<typename TGlobalHolder::TRedAlph>())
//                 if (++count > lH.options.maxSeedDist)
//                     return false;
//         return true;
//     }

    int64_t effectiveQBegin = m.qryStart;
    int64_t effectiveSBegin = m.subjStart;
    uint64_t actualLength = m.qryEnd - m.qryStart;
    uint64_t effectiveLength = std::max(static_cast<uint64_t>(lH.options.seedLength * lH.options.preScoring),
                                        actualLength);

    if (effectiveLength > actualLength)
    {
        effectiveQBegin -= (lH.options.preScoring - 1) *
                           actualLength / 2;
        effectiveSBegin -= (lH.options.preScoring - 1) *
                           actualLength / 2;
//         std::cout << effectiveQBegin << "\t" << effectiveSBegin << "\n";
        int64_t min = std::min(effectiveQBegin, effectiveSBegin);
        if (min < 0)
        {
            effectiveQBegin -= min;
            effectiveSBegin -= min;
            effectiveLength += min;
        }

        effectiveLength = std::min({
                            length(lH.gH.qrySeqs[m.qryId]) - effectiveQBegin,
                            length(lH.gH.subjSeqs[m.subjId]) - effectiveSBegin,
                            effectiveLength});
//         std::cout << effectiveQBegin << "\t" << effectiveSBegin << "\t"
//                   << effectiveLength << "\n";
    }

    auto const & qSeq = infix(lH.gH.qrySeqs[m.qryId],
                              effectiveQBegin,
                              effectiveQBegin + effectiveLength);
    auto const & sSeq = infix(lH.gH.subjSeqs[m.subjId],
                              effectiveSBegin,
                              effectiveSBegin + effectiveLength);
    int maxScore = 0;

    int scores[effectiveLength+1]; // C99, C++14, -Wno-vla before that
    scores[0] = 0;

    // score the diagonal
    for (uint64_t i = 0; i < effectiveLength; ++i)
    {
        scores[i] += score(seqanScheme(context(lH.gH.outfile).scoringScheme), qSeq[i], sSeq[i]);
        if (scores[i] < 0)
            scores[i] = 0;
        else if (scores[i] > maxScore)
            maxScore = scores[i];

//         if (i < static_cast<uint64_t>(effectiveLength - 1)) // TODO remove
        scores[i+1] = scores[i];
    }

    return (maxScore >= int(lH.options.preScoringThresh * effectiveLength));
}

// --------------------------------------------------------------------------
// Function onFind()
// --------------------------------------------------------------------------

template <typename TGlobalHolder,
          typename TScoreExtension,
          typename TSeedId,
          typename TSubjOcc>
inline void
onFind(LocalDataHolder<TGlobalHolder, TScoreExtension> & lH,
       TSeedId const & seedId,
       TSubjOcc subjOcc)
{
    using TMatch = typename TGlobalHolder::TMatch;
    SEQAN_ASSERT_LEQ_MSG(getSeqOffset(subjOcc) + lH.options.seedLength,
                         length(lH.gH.subjSeqs[getSeqNo(subjOcc)]),
                         "ERROR: Seed reaches beyond end of subject sequence! Please report a bug with your files at "
                         "http://www.seqan.de/lambda !");

    if (TGlobalHolder::indexIsFM) // positions are reversed
        setSeqOffset(subjOcc,
                     length(lH.gH.subjSeqs[getSeqNo(subjOcc)])
                     - getSeqOffset(subjOcc)
                     - lH.options.seedLength);

    Match m {static_cast<Match::TQId>(lH.seedRefs[seedId]),
             static_cast<Match::TSId>(getSeqNo(subjOcc)),
             static_cast<Match::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset),
             static_cast<Match::TPos>(lH.seedRanks[seedId] * lH.options.seedOffset + lH.options.seedLength),
             static_cast<Match::TPos>(getSeqOffset(subjOcc)),
             static_cast<Match::TPos>(getSeqOffset(subjOcc) + lH.options.seedLength)};

    bool discarded = false;
    auto const halfSubjL = lH.options.seedLength /  2;

    if (!sIsTranslated(lH.gH.blastProgram))
    {
        for (unsigned k = 0; k < length(lH.gH.segIntStarts[m.subjId]); ++k)
        {
            // more than half of the seed falls into masked interval
            if (intervalOverlap(m.subjStart,
                                m.subjEnd,
                                lH.gH.segIntStarts[m.subjId][k],
                                lH.gH.segIntEnds[m.subjId][k])
                    >= halfSubjL)
            {
                ++lH.stats.hitsMasked;
                discarded = true;
                break;
            }
        }
    }

     if ((!discarded) && (!seedLooksPromising(lH, m)))
     {
         discarded = true;
         ++lH.stats.hitsFailedPreExtendTest;
     }

    if (!discarded)
        lH.matches.emplace_back(m);
}

template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension,
          typename TSubjOcc>
inline void
onFindVariable(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
               TSubjOcc subjOcc,
               Match::TQId const seedId,
               Match::TPos const seedBegin,
               Match::TPos const seedLength)
{
    if (TGlobalHolder::indexIsFM) // positions are reversed
        setSeqOffset(subjOcc,
                     length(lH.gH.subjSeqs[getSeqNo(subjOcc)])
                     - getSeqOffset(subjOcc)
                     - seedLength);

    Match m {seedId,
             static_cast<Match::TSId>(getSeqNo(subjOcc)),
             seedBegin,
             static_cast<Match::TPos>(seedBegin + seedLength),
             static_cast<Match::TPos>(getSeqOffset(subjOcc)),
             static_cast<Match::TPos>(getSeqOffset(subjOcc) + seedLength)};

     if (!seedLooksPromising(lH, m))
         ++lH.stats.hitsFailedPreExtendTest;
     else
        lH.matches.emplace_back(m);
}

// --------------------------------------------------------------------------
// Function search()
// --------------------------------------------------------------------------

template <typename TIndexIt, typename TNeedleIt, typename TLambda, typename TLambda2>
inline void
__goDownNoErrors(TIndexIt const & indexIt,
                 TNeedleIt const & needleIt,
                 TNeedleIt const & needleItEnd,
                 TLambda & continRunnable,
                 TLambda2 & reportRunnable)
{
    TIndexIt nextIndexIt(indexIt);
    if ((needleIt != needleItEnd) &&
        goDown(nextIndexIt, *needleIt) &&
        continRunnable(indexIt, nextIndexIt))
    {
        __goDownNoErrors(nextIndexIt, needleIt + 1, needleItEnd, continRunnable, reportRunnable);
    } else
    {
        reportRunnable(indexIt);
    }
}

template <typename TIndexIt, typename TNeedleIt, typename TLambda, typename TLambda2>
inline void
__goDownErrors(TIndexIt const & indexIt,
               TNeedleIt const & needleIt,
               TNeedleIt const & needleItEnd,
               TLambda & continRunnable,
               TLambda2 & reportRunnable)
{
    using TAlph = typename Value<TNeedleIt>::Type;

    unsigned contin = 0;

    if (needleIt != needleItEnd)
    {
        for (unsigned i = 0; i < ValueSize<TAlph>::VALUE; ++i)
        {
            TIndexIt nextIndexIt(indexIt);
            if (goDown(nextIndexIt, static_cast<TAlph>(i)) &&
                continRunnable(indexIt, nextIndexIt))
            {
                ++contin;
                if (ordValue(*needleIt) == i)
                    __goDownErrors(nextIndexIt, needleIt + 1, needleItEnd, continRunnable, reportRunnable);
                else
                    __goDownNoErrors(nextIndexIt, needleIt + 1, needleItEnd, continRunnable, reportRunnable);
            }
        }
    }

    if (contin == 0)
        reportRunnable(indexIt);
}

template <typename TMatch,
          typename TGlobalHolder,
          typename TScoreExtension>
inline void
__serachAdaptive(LocalDataHolder<TMatch, TGlobalHolder, TScoreExtension> & lH,
                 uint64_t const seedLength)
{
    typedef typename Iterator<typename TGlobalHolder::TDbIndex, TopDown<> >::Type TIndexIt;

    // TODO optionize
    size_t constexpr seedHeurFactor = 1;//2;
    size_t constexpr minResults = 1;

    size_t needlesSum = lengthSum(infix(lH.gH.redQrySeqs, lH.indexBeginQry, lH.indexEndQry));
    size_t needlesPos = 0;

    TIndexIt root(lH.gH.dbIndex);
    TIndexIt indexIt = root;

    for (size_t i = lH.indexBeginQry; i < lH.indexEndQry; ++i)
    {
        for (size_t seedBegin = 0; seedBegin <= length(lH.gH.redQrySeqs[i]) - seedLength;)
        {
            indexIt = root;
            uint64_t maxSeedExtension = 0;

            auto continRunnable = [&] (TIndexIt const & prevIndexIt, TIndexIt const & indexIt)
            {
                // NON-ADAPTIVE
//              return (repLength(indexIt) <= seedLength);
                // ADAPTIVE SEEDING:

                // always continue if minimum seed length not reached
                if (repLength(indexIt) <= seedLength)
                    return true;

                // always continue if it means not loosing hits
                if (countOccurrences(indexIt) == countOccurrences(prevIndexIt))
                    return true;

                // do vodoo heuristics to see if this hit is to frequent
                size_t desiredOccs = length(lH.matches) >= lH.options.maxMatches
                                        ? minResults
                                        : (lH.options.maxMatches - length(lH.matches)) * seedHeurFactor / ((needlesSum - needlesPos - repLength(indexIt)) / repLength(indexIt));
//                                      : lH.options.maxMatches / (length(lH.gH.redQrySeqs[i]) / repLength(indexIt));
                if (desiredOccs == 0)
                    desiredOccs = 1u;
                if (countOccurrences(indexIt) < desiredOccs)
                    return false;

                return true;
            };


            auto reportRunnable = [&] (TIndexIt const & indexIt)
            {
                if (repLength(indexIt) >= seedLength)
                {
                    maxSeedExtension = std::max(maxSeedExtension, repLength(indexIt));
                    lH.stats.hitsAfterSeeding += countOccurrences(indexIt);
                    for (auto const & occ : getOccurrences(indexIt))
                        onFindVariable(lH, occ, i, seedBegin, repLength(indexIt));// pos - seedBegin);
                }
            };

            __goDownErrors(indexIt, begin(lH.gH.redQrySeqs[i], Standard()) + seedBegin, end(lH.gH.redQrySeqs[i], Standard()), continRunnable, reportRunnable);

            // set beginning of next seed (-2 so we have some overlap)
            if (maxSeedExtension <= seedLength - lH.options.seedOffset)
//                     maxSeedExtension = lH.options.seedOffset;
                ++seedBegin;
            else
                seedBegin += maxSeedExtension + lH.options.seedOffset - lH.options.seedLength;
//                 seedBegin += lH.options.seedOffset;

        }

        needlesPos += length(lH.gH.redQrySeqs[i]);
    }
}

template <typename BackSpec, typename TLocalHolder>
inline void
__searchDoubleIndex(TLocalHolder & lH)
{
    appendToStatus(lH.statusStr, lH.options, 1, "Seeding...");
    if (lH.options.isTerm)
        myPrint(lH.options, 1, lH.statusStr);

    double start = sysTime();

    using LambdaFinder = Finder_<decltype(lH.gH.dbIndex),
                                 decltype(lH.seedIndex),
                                 typename std::conditional<std::is_same<BackSpec, Backtracking<Exact>>::value,
                                                           Backtracking<HammingDistance>,
                                                           BackSpec>::type >;

    LambdaFinder finder;

    auto delegate = [&lH] (LambdaFinder const & finder)
    {
        auto qryOccs = getOccurrences(back(finder.patternStack));
        auto subjOccs = getOccurrences(back(finder.textStack));

        lH.stats.hitsAfterSeeding += length(qryOccs) * length(subjOccs);

        for (unsigned i = 0; i < length(qryOccs); ++i)
            for (unsigned j = 0; j < length(subjOccs); ++j)
                onFind(lH, getSeqNo(qryOccs[i]), subjOccs[j]);
    };

    _find(finder, lH.gH.dbIndex, lH.seedIndex, lH.options.maxSeedDist, delegate);

    double finish = sysTime() - start;

    appendToStatus(lH.statusStr, lH.options, 1, " done. ");
    appendToStatus(lH.statusStr, lH.options, 2, finish, "s. #hits: ",
                   length(lH.matches), " ");
    myPrint(lH.options, 1, lH.statusStr);
}

template <typename BackSpec, typename TLocalHolder>
inline void
__searchSingleIndex(TLocalHolder & lH)
{
    typedef typename Iterator<decltype(lH.seeds) const, Rooted>::Type TSeedsIt;
    typedef typename Iterator<decltype(lH.gH.dbIndex),TopDown<>>::Type TIndexIt;

//     SEQAN_OMP_PRAGMA(critical(stdout))
//     {
//         std::cout << "ReadId: " << lH.i << std::endl;
//         for (auto const & seed : lH.seeds)
//             std::cout << "\"" << toCString(CharString(seed)) << "\" ";
//         std::cout << std::endl;
//     }
    auto delegate = [&lH] (TIndexIt & indexIt,
                           TSeedsIt const & seedsIt,
                           int /*score*/)
    {
        auto qryOcc = position(seedsIt);
        auto subjOccs = getOccurrences(indexIt);

        lH.stats.hitsAfterSeeding += length(subjOccs);

        for (unsigned j = 0; j < length(subjOccs); ++j)
            onFind(lH, qryOcc, subjOccs[j]);
    };

    find(lH.gH.dbIndex, lH.seeds, int(lH.options.maxSeedDist), delegate,
         Backtracking<BackSpec>());
}

template <typename BackSpec, typename TLocalHolder>
inline void
__search(TLocalHolder & lH)
{
    if (lH.options.doubleIndexing)
        __searchDoubleIndex<BackSpec>(lH);
    else
        __searchSingleIndex<BackSpec>(lH);
}

template <typename TLocalHolder>
inline void
search(TLocalHolder & lH)
{
    if (lH.options.adaptiveSeeding)
        __serachAdaptive(lH, lH.options.seedLength);
    else if (lH.options.maxSeedDist == 0)
        __search<Backtracking<Exact>>(lH);
    else if (lH.options.hammingOnly)
        __search<Backtracking<HammingDistance>>(lH);
    else
#if 0 // reactivate if edit-distance seeding is readded
        __search<Backtracking<EditDistance>>(lH);
#else
        return;
#endif
}

// --------------------------------------------------------------------------
// Function joinAndFilterMatches()
// --------------------------------------------------------------------------


template <typename TLocalHolder>
inline void
sortMatches(TLocalHolder & lH)
{
    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, "Sorting hits...");
        if (lH.options.isTerm)
            myPrint(lH.options, 1, lH.statusStr);
    }

    double start = sysTime();

//    std::sort(begin(lH.matches, Standard()), end(lH.matches, Standard()));
//     std::sort(lH.matches.begin(), lH.matches.end());

//     if (lH.matches.size() > lH.options.maxMatches)
//     {
//         MatchSortComp   comp(lH.matches);
//         std::sort(lH.matches.begin(), lH.matches.end(), comp);
//     } else

    if ((lH.options.filterPutativeAbundant) &&
        (lH.matches.size() > lH.options.maxMatches))
        // more expensive sort to get likely targets to front
        myHyperSortSingleIndex(lH.matches, lH.options.doubleIndexing, lH.gH);
    else
        std::sort(lH.matches.begin(), lH.matches.end());

    double finish = sysTime() - start;

    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1, " done. ");
        appendToStatus(lH.statusStr, lH.options, 2, finish, "s. ");
        myPrint(lH.options, 1, lH.statusStr);
    }
}

template <typename TBlastMatch,
          typename TLocalHolder>
inline int
computeBlastMatch(TBlastMatch         & bm,
                  typename TLocalHolder::TMatch const & m,
                  TLocalHolder        & lH)
{
    using TMatch = typename TLocalHolder::TMatch;
    using TPos   = typename TMatch::TPos;

    const unsigned long qryLength = length(value(lH.gH.qrySeqs, m.qryId));

    SEQAN_ASSERT_LEQ(bm.qStart, bm.qEnd);
    SEQAN_ASSERT_LEQ(bm.sStart, bm.sEnd);

//     auto qryInfix = infix(curQry,
//                                 bm.qStart,
//                                 bm.qEnd);
//     auto subjInfix  = infix(curSubj,
//                                 bm.sStart,
//                                 bm.sEnd);

//     std::cout << "Query Id: " << m.qryId
//               << "\t TrueQryId: " << getTrueQryId(bm.m, lH.options, lH.gH.blastProgram)
//               << "\t length(qryIds): " << length(qryIds)
//               << "Subj Id: " << m.subjId
//               << "\t TrueSubjId: " << getTrueSubjId(bm.m, lH.options, lH.gH.blastProgram)
//               << "\t length(subjIds): " << length(subjIds) << "\n\n";

    assignSource(bm.alignRow0, infix(lH.gH.qrySeqs[m.qryId], bm.qStart, bm.qEnd));
    assignSource(bm.alignRow1, infix(lH.gH.subjSeqs[m.subjId],bm.sStart, bm.sEnd));

//     std::cout << "== Positions\n";
//     std::cout << "   " <<  bm.qStart << " - " << bm.qEnd << " [before ali]\n";
//     std::cout << bm.align << std::endl;

    int scr = 0;

//     unsigned short seedLeng = 0;
//     double         seedE = 0;
//     double         seedB = 0;

    TPos row0len = bm.qEnd - bm.qStart;
    TPos row1len = bm.sEnd - bm.sStart;
    TPos band = (!lH.options.hammingOnly) * (lH.options.maxSeedDist);

//     // TODO FIGURE THIS OUT
//     if ((row0len > (lH.options.seedLength + band)) ||
//         (row1len > (lH.options.seedLength + band)))
//     {
//         #pragma omp atomic
//         ++mergeCount;
//         std::cout << "qrId " << m.qryId << "\tsId: " << m.subjId << "\trow0len: " << row0len << "\trow1len: " << row1len << "\n";
//         std::cout << source(row0) << "\n";
//         std::cout << source(row1) << "\n";
//     }

    auto seedsInSeed = std::max(row0len, row1len) / lH.options.seedLength;

    TPos  maxDist =  0;
    if (lH.options.maxSeedDist <= 1)
        maxDist = std::abs(int(row1len) - int(row0len));
    else
        maxDist = std::abs(int(row1len) - int(row0len)) + (seedsInSeed * band);

    // fast local alignment without DP-stuff
    if (maxDist == 0)
    {
        int scores[row0len+1]; // C99, C++14, -Wno-vla before that
        scores[0] = 0;
        unsigned newEnd = 0;
        unsigned newBeg = 0;
        // score the diagonal
        for (unsigned i = 0; i < row0len; ++i)
        {
            scores[i] += score(seqanScheme(context(lH.gH.outfile).scoringScheme),
                               source(bm.alignRow0)[i],
                               source(bm.alignRow1)[i]);
            if (scores[i] < 0)
            {
                scores[i] = 0;
            } else if (scores[i] >= scr)
            {
                scr = scores[i];
                newEnd = i + 1;
            }
//             if (i <row0len -1)
            scores[i+1] = scores[i];
        }
        if (newEnd == 0) // no local alignment
        {
            return OTHER_FAIL; // TODO change to PREEXTEND?
        }
        // backtrack
        for (unsigned i = newEnd - 1; i > 0; --i)
        {
            if (scores[i] == 0)
            {
                newBeg = i + 1;
                break;
            }
        }
        setEndPosition(bm.alignRow0, newEnd);
        setEndPosition(bm.alignRow1, newEnd);
        setBeginPosition(bm.alignRow0, newBeg);
        setBeginPosition(bm.alignRow1, newBeg);

    } else
    {
        // compute with DP-code
//         scr = localAlignment(bm.align, seqanScheme(context(lH.gH.outfile).scoringScheme), -maxDist, +maxDist);
        scr = localAlignment2(bm.alignRow0,
                              bm.alignRow1,
                              seqanScheme(context(lH.gH.outfile).scoringScheme),
                              -maxDist,
                              +maxDist,
                              lH.alignContext);
    }

    // save new bounds of alignment
    bm.qEnd    =  bm.qStart  + endPosition(bm.alignRow0);
    bm.qStart  += beginPosition(bm.alignRow0);

    bm.sEnd   =  bm.sStart + endPosition(bm.alignRow1);
    bm.sStart += beginPosition(bm.alignRow1);

//     if (scr < lH.options.minSeedScore)
//         return PREEXTEND;

#if 0
// OLD WAY extension with birte's code
    {
    //     std::cout << "   " <<  bm.qStart << " - " << bm.qEnd << " [after ali]\n";
    //     std::cout << bm.align << std::endl;
        decltype(seqanScheme(context(lH.gH.outfile).scoringScheme)) extScheme(seqanScheme(context(lH.gH.outfile).scoringScheme));
        setScoreGapOpen  (extScheme, -8);
        setScoreGapExtend(extScheme, -8);

        Seed<Simple>           seed(bm.sStart, bm.qStart,
                                    bm.sEnd, bm.qEnd);
        extendSeed(seed,
                curSubj,
                curQry,
                EXTEND_BOTH,
                extScheme,
//                    seqanScheme(context(lH.gH.outfile).scoringScheme),
                int(lH.options.xDropOff),
                GappedXDrop());

        bm.sStart = beginPositionH(seed);
        bm.qStart = beginPositionV(seed);
        bm.sEnd   = endPositionH(seed);
        bm.qEnd   = endPositionV(seed);

        assignSource(row0, infix(curQry,
                                bm.qStart,
                                bm.qEnd));
        assignSource(row1, infix(curSubj,
                                bm.sStart,
                                bm.sEnd));

        //DEBUG
        auto oldscr = scr;

        scr = localAlignment(bm.align,
                            seqanScheme(context(lH.gH.outfile).scoringScheme),
//                                 alignConfig,
                            lowerDiagonal(seed)-beginDiagonal(seed),
                            upperDiagonal(seed)-beginDiagonal(seed));
        // save new bounds of alignment
        bm.qEnd    =  bm.qStart  + endPosition(row0);
        bm.qStart  += beginPosition(row0);

        bm.sEnd   =  bm.sStart + endPosition(row1);
        bm.sStart += beginPosition(row1);

        if (scr < 0) // alignment got screwed up
        {
            std::cout << "SCREW UP\n";
            std::cout << "beginDiag: " << beginDiagonal(seed)
                << "\tlowDiag: " << lowerDiagonal(seed)
                << "\tupDiag: " << upperDiagonal(seed) << '\n';
            std::cout << "oldscore: " << oldscr
                    << "\tseedscore: " << score(seed)
                    << "\tscore: " << scr << '\n';
            std::cout << bm.align << '\n';
        }
    }
#endif

#if 0
    // ungapped second prealign
    {
        Tuple<decltype(bm.qStart), 4> positions =
            { { bm.qStart, bm.sStart, bm.qEnd, bm.sEnd} };

        decltype(seqanScheme(context(lH.gH.outfile).scoringScheme)) extScheme(seqanScheme(context(lH.gH.outfile).scoringScheme));
        setScoreGapOpen  (extScheme, -100);
        setScoreGapExtend(extScheme, -100);
        scr = extendAlignment(bm.align,
                                lH.alignContext,
                                scr,
                                curQry,
                                curSubj,
                                positions,
                                EXTEND_BOTH,
                                0, // band of 0 size
                                0, // band of 0 size
                                1, // xdrop of 1
                                extScheme);
        bm.qStart  = beginPosition(row0);
        bm.qEnd    = endPosition(row0);

        bm.sStart =  beginPosition(row1);
        bm.sEnd   =  endPosition(row1);
    }
#endif
    if (((bm.qStart > 0) && (bm.sStart > 0)) ||
        ((bm.qEnd < qryLength - 1) && (bm.sEnd < length(lH.gH.subjSeqs[m.subjId]) -1)))
    {
        // we want to allow more gaps in longer query sequences
        switch (lH.options.band)
        {
            case -3: maxDist = ceil(log2(qryLength)); break;
            case -2: maxDist = floor(sqrt(qryLength)); break;
            case -1: break;
            default: maxDist = lH.options.band; break;
        }

        Tuple<decltype(bm.qStart), 4> positions =
                { { bm.qStart, bm.sStart, bm.qEnd, bm.sEnd} };

        if (lH.options.band != -1)
        {
            if (lH.options.xDropOff != -1)
            {
                scr = _extendAlignmentImpl(bm.alignRow0,
                                           bm.alignRow1,
                                           scr,
                                           lH.gH.qrySeqs[m.qryId],
                                           lH.gH.subjSeqs[m.subjId],
                                           positions,
                                           EXTEND_BOTH,
                                           -maxDist,
                                           +maxDist,
                                           lH.options.xDropOff,
                                           seqanScheme(context(lH.gH.outfile).scoringScheme),
                                           True(),
                                           True(),
                                           lH.alignContext);
            } else
            {
                scr = _extendAlignmentImpl(bm.alignRow0,
                                           bm.alignRow1,
                                           scr,
                                           lH.gH.qrySeqs[m.qryId],
                                           lH.gH.subjSeqs[m.subjId],
                                           positions,
                                           EXTEND_BOTH,
                                           -maxDist,
                                           +maxDist,
                                           lH.options.xDropOff,
                                           seqanScheme(context(lH.gH.outfile).scoringScheme),
                                           True(),
                                           False(),
                                           lH.alignContext);
            }
        } else
        {
            if (lH.options.xDropOff != -1)
            {
                scr = _extendAlignmentImpl(bm.alignRow0,
                                           bm.alignRow1,
                                           scr,
                                           lH.gH.qrySeqs[m.qryId],
                                           lH.gH.subjSeqs[m.subjId],
                                           positions,
                                           EXTEND_BOTH,
                                           -maxDist,
                                           +maxDist,
                                           lH.options.xDropOff,
                                           seqanScheme(context(lH.gH.outfile).scoringScheme),
                                           False(),
                                           True(),
                                           lH.alignContext);
            } else
            {
                scr = _extendAlignmentImpl(bm.alignRow0,
                                           bm.alignRow1,
                                           scr,
                                           lH.gH.qrySeqs[m.qryId],
                                           lH.gH.subjSeqs[m.subjId],
                                           positions,
                                           EXTEND_BOTH,
                                           -maxDist,
                                           +maxDist,
                                           lH.options.xDropOff,
                                           seqanScheme(context(lH.gH.outfile).scoringScheme),
                                           False(),
                                           False(),
                                           lH.alignContext);
            }
        }
        bm.sStart = beginPosition(bm.alignRow1);
        bm.qStart = beginPosition(bm.alignRow0);
        bm.sEnd   = endPosition(bm.alignRow1);
        bm.qEnd   = endPosition(bm.alignRow0);

//         std::cout << "AFTER:\n" << bm.align << "\n";
    }

//     std::cerr << "AFTEREXT:\n "<< bm.align << "\n";

    if (scr <= 0)
    {
//         std::cout << "## LATE FAIL\n" << bm.align << '\n';

        return OTHER_FAIL;
    }
//     std::cout << "##LINE: " << __LINE__ << '\n';

//     std::cout << "ALIGN BEFORE STATS:\n" << bm.align << "\n";

    computeAlignmentStats(bm, context(lH.gH.outfile));

    if (bm.alignStats.alignmentIdentity < lH.options.idCutOff)
        return PERCENTIDENT;

//     const unsigned long qryLength = length(row0);
    computeBitScore(bm, context(lH.gH.outfile));

    // the length adjustment cache must no be written to by multiple threads
    SEQAN_OMP_PRAGMA(critical(evalue_length_adj_cache))
    {
        computeEValue(bm, context(lH.gH.outfile));
    }

    if (bm.eValue > lH.options.eCutOff)
    {
        return EVALUE;
    }

    if (qIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = (m.qryId % 3) + 1;
        if (m.qryId % 6 > 2)
            bm.qFrameShift = -bm.qFrameShift;
    } else if (qHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.qFrameShift = 1;
        if (m.qryId % 2)
            bm.qFrameShift = -bm.qFrameShift;
    } else
    {
        bm.qFrameShift = 0;
    }

    if (sIsTranslated(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = (m.subjId % 3) + 1;
        if (m.subjId % 6 > 2)
            bm.sFrameShift = -bm.sFrameShift;
    } else if (sHasRevComp(TLocalHolder::TGlobalHolder::blastProgram))
    {
        bm.sFrameShift = 1;
        if (m.subjId % 2)
            bm.sFrameShift = -bm.sFrameShift;
    } else
    {
        bm.sFrameShift = 0;
    }

    return 0;
}


template <typename TLocalHolder>
inline int
iterateMatches(TLocalHolder & lH)
{
    using TGlobalHolder = typename TLocalHolder::TGlobalHolder;
//     using TMatch        = typename TGlobalHolder::TMatch;
//     using TPos          = typename TMatch::TPos;
    using TBlastPos     = uint32_t; //TODO why can't this be == TPos
    using TBlastMatch   = BlastMatch<
                           typename TLocalHolder::TAlignRow0,
                           typename TLocalHolder::TAlignRow1,
                           TBlastPos,
                           typename Value<typename TGlobalHolder::TQryIds>::Type,// const &,
                           typename Value<typename TGlobalHolder::TSubjIds>::Type// const &,
                           >;
    using TBlastRecord  = BlastRecord<TBlastMatch>;

//     constexpr TPos TPosMax = std::numeric_limits<TPos>::max();
//     constexpr uint8_t qFactor = qHasRevComp(lH.gH.blastProgram) ? 3 : 1;
//     constexpr uint8_t sFactor = sHasRevComp(lH.gH.blastProgram) ? 3 : 1;

    double start = sysTime();
    if (lH.options.doubleIndexing)
    {
        appendToStatus(lH.statusStr, lH.options, 1,
                       "Extending and writing hits...");
        myPrint(lH.options, 1, lH.statusStr);
    }

    //DEBUG
//     std::cout << "Length of matches:   " << length(lH.matches);
//     for (auto const & m :  lH.matches)
//     {
//         std::cout << m.qryId << "\t" << getTrueQryId(m,lH.options, lH.gH.blastProgram) << "\n";
//     }

//     double topMaxMatchesMedianBitScore = 0;
    // outer loop over records
    // (only one iteration if single indexing is used)
    for (auto it = lH.matches.begin(),
              itN = std::next(it, 1),
              itEnd = lH.matches.end();
         it != itEnd;
         ++it)
    {
        itN = std::next(it,1);
        auto const trueQryId = it->qryId / qNumFrames(lH.gH.blastProgram);

        TBlastRecord record(lH.gH.qryIds[trueQryId]);

        record.qLength = (qIsTranslated(lH.gH.blastProgram)
                            ? lH.gH.untransQrySeqLengths[trueQryId]
                            : length(lH.gH.qrySeqs[it->qryId]));

//         topMaxMatchesMedianBitScore = 0;

        // inner loop over matches per record
        for (; it != itEnd; ++it)
        {
            auto const trueSubjId = it->subjId / sNumFrames(lH.gH.blastProgram);
            itN = std::next(it,1);
//             std::cout << "FOO\n" << std::flush;
//             std::cout << "QryStart: " << it->qryStart << "\n" << std::flush;
//             std::cout << "SubjStart: " << it->subjStart << "\n" << std::flush;
//             std::cout << "BAR\n" << std::flush;
            if (!isSetToSkip(*it))
            {
                // ABUNDANCY and PUTATIVE ABUNDANCY CHECKS
                if ((lH.options.filterPutativeAbundant) && (record.matches.size() % lH.options.maxMatches == 0))
                {
                    if (record.matches.size() / lH.options.maxMatches == 1)
                    {
                        // numMaxMatches found the first time
                        record.matches.sort();
                    }
                    else if (record.matches.size() / lH.options.maxMatches > 1)
                    {
                        double medianTopNMatchesBefore = 0.0;
//                         if (lH.options.filterPutativeAbundant)
                        {
                            medianTopNMatchesBefore =
                            (std::next(record.matches.begin(),
                                       lH.options.maxMatches / 2))->bitScore;
                        }

                        uint64_t before = record.matches.size();
                        record.matches.sort();
                        // if we filter putative duplicates we never need to check for real duplicates
                        if (!lH.options.filterPutativeDuplicates)
                        {
                            record.matches.unique();
                            lH.stats.hitsDuplicate += before - record.matches.size();
                            before = record.matches.size();
                        }
                        if (record.matches.size() > (lH.options.maxMatches + 1))
                            // +1 so as not to trigger % == 0 in the next run
                            record.matches.resize(lH.options.maxMatches + 1);

                        lH.stats.hitsAbundant += before - record.matches.size();

//                         if (lH.options.filterPutativeAbundant)
                        {
                            double medianTopNMatchesAfter =
                            (std::next(record.matches.begin(),
                                       lH.options.maxMatches / 2))->bitScore;
                            // no new matches in top n/2
                            if (int(medianTopNMatchesAfter) <=
                                int(medianTopNMatchesBefore))
                            {
                                // declare all the rest as putative abundant
                                while ((it != itEnd) &&
                                       (trueQryId == it->qryId / qNumFrames(lH.gH.blastProgram)))
                                {
                                    // not already marked as abundant, duplicate or merged
                                    if (!isSetToSkip(*it))
                                        ++lH.stats.hitsPutativeAbundant;
                                    ++it;
                                }
                                // move back so if-loop's increment still valid
                                std::advance(it, -1);
                                break;
                            }
                        }
                    }
                }
//                 std::cout << "BAX\n" << std::flush;
                // create blastmatch in list without copy or move
                record.matches.emplace_back(lH.gH.qryIds [trueQryId],
                                            lH.gH.subjIds[trueSubjId]);

                auto & bm = back(record.matches);

                bm.qStart    = it->qryStart;
                bm.qEnd      = it->qryEnd; // it->qryStart + lH.options.seedLength;
                bm.sStart    = it->subjStart;
                bm.sEnd      = it->subjEnd;//it->subjStart + lH.options.seedLength;

                bm.qLength = record.qLength;
                bm.sLength = sIsTranslated(lH.gH.blastProgram)
                                ? lH.gH.untransSubjSeqLengths[trueSubjId]
                                : length(lH.gH.subjSeqs[it->subjId]);

                // MERGE PUTATIVE SIBLINGS INTO THIS MATCH
                if (lH.options.mergePutativeSiblings)
                {
                    for (auto it2 = itN;
                        (it2 != itEnd) &&
                        (trueQryId == it2->qryId / qNumFrames(lH.gH.blastProgram)) &&
                        (trueSubjId == it2->subjId / sNumFrames(lH.gH.blastProgram));
                        ++it2)
                    {
                        // same frame
                        if ((it->qryId % qNumFrames(lH.gH.blastProgram) == it2->qryId % qNumFrames(lH.gH.blastProgram)) &&
                            (it->subjId % sNumFrames(lH.gH.blastProgram) == it2->subjId % sNumFrames(lH.gH.blastProgram)))
                        {

    //                         TPos const qDist = (it2->qryStart >= bm.qEnd)
    //                                             ? it2->qryStart - bm.qEnd // upstream
    //                                             : 0; // overlap
    //
    //                         TPos sDist = TPosMax; // subj match region downstream of *it
    //                         if (it2->subjStart >= bm.sEnd) // upstream
    //                             sDist = it2->subjStart - bm.sEnd;
    //                         else if (it2->subjStart >= it->subjStart) // overlap
    //                             sDist = 0;

                            // due to sorting it2->qryStart never <= it->qStart
                            // so subject sequences must have same order
                            if (it2->subjStart < it->subjStart)
                                continue;

                            long const qDist = it2->qryStart - bm.qEnd;
                            long const sDist = it2->subjStart - bm.sEnd;

                        if ((qDist == sDist) &&
                            (qDist <= (long)lH.options.seedGravity))
                        {
                            bm.qEnd = std::max(bm.qEnd,
                                               static_cast<TBlastPos>(it2->qryStart
                                               + lH.options.seedLength));
                            bm.sEnd = std::max(bm.sEnd,
                                               static_cast<TBlastPos>(it2->subjStart
                                               + lH.options.seedLength));
                            ++lH.stats.hitsMerged;

                                setToSkip(*it2);
                            }
                        }
                    }
                }

                // do the extension and statistics
                int lret = computeBlastMatch(bm, *it, lH);

                switch (lret)
                {
                    case COMPUTERESULT_::SUCCESS:
    //                     ++lH.stats.goodMatches;
                        if (lH.options.outFileFormat > 0)
                        {
                            bm._n_qId = it->qryId / qNumFrames(lH.gH.blastProgram);
                            bm._n_sId = it->subjId / sNumFrames(lH.gH.blastProgram);
                        }
                        break;
                    case EVALUE:
                        ++lH.stats.hitsFailedExtendEValueTest;
                        break;
                    case PERCENTIDENT:
                        ++lH.stats.hitsFailedExtendPercentIdentTest;
                        break;
                    case PREEXTEND:
                        ++lH.stats.hitsFailedPreExtendTest;
                        break;
                    default:
                        std::cerr << "Unexpected Extension Failure:\n"
                          << "qryId: " << it->qryId << "\t"
                          << "subjId: " << it->subjId << "\t"
                          << "seed    qry: " << infix(lH.gH.redQrySeqs,
                                                      it->qryStart,
                                                      it->qryEnd)
//                                                       it->qryStart + lH.options.seedLength)
                          << "\n       subj: " << infix(lH.gH.redSubjSeqs,
                                                      it->subjStart,
                                                      it->subjEnd)
//                                                       it->subjStart + lH.options.seedLength)
                          << "\nunred  qry: " << infix(lH.gH.qrySeqs,
                                                      it->qryStart,
                                                      it->qryEnd)
//                                                       it->qryStart + lH.options.seedLength)
                          << "\n       subj: " << infix(lH.gH.subjSeqs,
                                                      it->subjStart,
                                                      it->subjEnd)
//                                                       it->subjStart + lH.options.seedLength)
                          << "\nmatch    qry: " << infix(lH.gH.qrySeqs,
                                                      bm.qStart,
                                                      bm.qEnd)
                          << "\n       subj: " << infix(lH.gH.subjSeqs,
                                                      bm.sStart,
                                                      bm.sEnd)
                          << "\nalign: " << bm.alignRow0 << "\n       " << bm.alignRow1
                          << "\n";
                        return lret;
                        break;
                }

                if (lret != 0)// discard match
                {
                    record.matches.pop_back();
                } else if (lH.options.filterPutativeDuplicates)
                {
                    // PUTATIVE DUBLICATES CHECK
                    for (auto it2 = itN;
                         (it2 != itEnd) &&
                         (trueQryId == it2->qryId / qNumFrames(lH.gH.blastProgram)) &&
                         (trueSubjId == it2->subjId / sNumFrames(lH.gH.blastProgram));
                         ++it2)
                    {
                        // same frame and same range
                        if ((it->qryId == it2->qryId) &&
                            (it->subjId == it2->subjId) &&
                            (intervalOverlap(it2->qryStart,
                                             it2->qryEnd,
//                                              it2->qryStart + lH.options.seedLength,
                                             bm.qStart,
                                             bm.qEnd) > 0) &&
                            (intervalOverlap(it2->subjStart,
                                             it2->subjEnd,
//                                              it2->subjStart + lH.options.seedLength,
                                             bm.sStart,
                                             bm.sEnd) > 0))
                        {
                            // deactivated alignment check to get rid of
                            // duplicates early on
//                             auto const & row0 = row(bm.align, 0);
//                             auto const & row1 = row(bm.align, 1);
//                             // part of alignment
//                             if (toSourcePosition(row0,
//                                                  toViewPosition(row1,
//                                                                 it2->subjStart
//                                                                 - bm.sStart))
//                                 == TPos(it2->qryStart - bm.qStart))
//                             {
                                ++lH.stats.hitsPutativeDuplicate;
                                setToSkip(*it2);
//                             }
                        }
                    }
                }
            }

            // last item or new TrueQryId
            if ((itN == itEnd) ||
                (trueQryId != itN->qryId / qNumFrames(lH.gH.blastProgram)))
                break;
        }

        if (length(record.matches) > 0)
        {
            ++lH.stats.qrysWithHit;
            // sort and remove duplicates -> STL, yeah!
            auto const before = record.matches.size();
            record.matches.sort();
            if (!lH.options.filterPutativeDuplicates)
            {
                record.matches.unique();
                lH.stats.hitsDuplicate += before - record.matches.size();
            }
            if (record.matches.size() > lH.options.maxMatches)
            {
                lH.stats.hitsAbundant += record.matches.size() -
                                         lH.options.maxMatches;
                record.matches.resize(lH.options.maxMatches);
            }
            lH.stats.hitsFinal += record.matches.size();

            myWriteRecord(lH, record);
        }

    }

    if (lH.options.doubleIndexing)
    {
        double finish = sysTime() - start;

        appendToStatus(lH.statusStr, lH.options, 1, " done. ");
        appendToStatus(lH.statusStr, lH.options, 2, finish, "s. ");
        myPrint(lH.options, 1, lH.statusStr);
    }

    return 0;
}

#endif // HEADER GUARD
