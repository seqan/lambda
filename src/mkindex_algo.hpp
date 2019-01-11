// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2019, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2019, Knut Reinert and Freie Universität Berlin
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
// lambda_indexer.hpp: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_LAMBDA_LAMBDA_INDEXER_H_
#define SEQAN_LAMBDA_LAMBDA_INDEXER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/translation.h>
#include <seqan/reduced_aminoacid.h>

#include "mkindex_misc.hpp"
#include "mkindex_saca.hpp"
#include "shared_misc.hpp"
#include "shared_options.hpp"
#include "search_output.hpp" //TODO only needed because options are in one file, remove later

using namespace seqan;

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TOrigAlph>
void
loadSubjSeqsAndIds(TCDStringSet<String<TOrigAlph>> & originalSeqs,
                   std::unordered_map<std::string, uint64_t> & accToIdRank,
                   LambdaIndexerOptions const & options)
{
    // Make sure we have enough RAM to load the file
    auto ram = getTotalSystemMemory();
    auto fS = fileSize(toCString(options.dbFile));

    if (fS >= ram)
        std::cerr << "WARNING: Your sequence file is already larger than your physical memory!\n"
                  << "         This means you will likely encounter a crash with \"bad_alloc\".\n"
                  << "         Split you sequence file into many smaller ones or use a computer\n"
                  << "         with more memory!\n";


    typedef TCDStringSet<String<char, Alloc<>>>             TIDs;

    TIDs ids; // the IDs

    // see http://www.uniprot.org/help/accession_numbers
    // https://www.ncbi.nlm.nih.gov/Sequin/acc.html
    // https://www.ncbi.nlm.nih.gov/refseq/about/
    // TODO: make sure these don't trigger twice on one ID
    std::regex const accRegEx{"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}|" // UNIPROT
                              "[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|"                                       // NCBI nucl
                              "[A-Z]{3}[0-9]{5}|"                                                     // NCBI prot
                              "[A-Z]{4}[0-9]{8,10}|"                                                  // NCBI wgs
                              "[A-Z]{5}[0-9]{7}|"                                                     // NCBI mga
                              "(NC|AC|NG|NT|NW|NZ|NM|NR|XM|XR|NP|AP|XP|YP|ZP)_[0-9]+|"                // RefSeq
                              "UPI[A-F0-9]{10}"};                                                     // UniParc

    // lambda that truncates IDs at first whitespace
    auto truncateID = [] (auto && id, uint64_t const)
    {
        IsWhitespace isWhitespace;
        for (size_t i = 0; i < length(id); ++i)
        {
            if (isWhitespace(id[i]))
            {
                resize(id, i);
                break;
            }
        }
    };

    uint64_t noAcc = 0;
    uint64_t multiAcc = 0;

    // lambda that extracts accession numbers and saves them in the map
    auto extractAccIds = [&accToIdRank, &accRegEx, &noAcc, &multiAcc] (auto && id, uint64_t const rank)
    {
        // TODO avoid copying here by specializing regex_iterator
        std::string buf;
        assign(buf, id);

        uint64_t count = 0;
        for (auto it = std::sregex_iterator(buf.begin(), buf.end(), accRegEx), itEnd = std::sregex_iterator();
             it != itEnd;
             ++it, ++count)
        {
            SEQAN_ASSERT_MSG(accToIdRank.count(it->str()) == 0,
                             "An accession number appeared twice in the file, but they should be unique.");
            // TODO store acc outside as well
            accToIdRank[it->str()] = rank;
        }

        switch (count)
        {
            case 0: ++noAcc; break;
            case 1: break;
            default: ++multiAcc; break;
        }
    };

    double start = sysTime();
    myPrint(options, 1, "Loading Subject Sequences and Ids...");

    SeqFileIn infile(toCString(options.dbFile));

    if (options.truncateIDs)
    {
        if (options.hasSTaxIds)
        {
            myReadRecords(ids, originalSeqs, infile, [&] (auto && id, uint64_t const rank)
            {
                extractAccIds(std::forward<decltype(id)>(id), rank);
                truncateID(std::forward<decltype(id)>(id), rank);
            });
        } else
        {
            myReadRecords(ids, originalSeqs, infile, truncateID);
        }
    } else
    {
        if (options.hasSTaxIds)
            myReadRecords(ids, originalSeqs, infile, extractAccIds);
        else
            myReadRecords(ids, originalSeqs, infile);
    }

    myPrint(options, 1,  " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n");

    if (length(originalSeqs) == 0)
    {
        throw std::runtime_error("ERROR: No sequences in file. Aborting.\n");
    }
    unsigned long maxLen = 0ul;
    for (auto const & s : originalSeqs)
    {
        if (length(s) > maxLen)
        {
            maxLen = length(s);
        }
        else if (length(s) == 0)
        {
            throw std::runtime_error("ERROR: Unexpectedly encountered a sequence of length 0 in the file."
                                     "Remove the entry and try again. Aborting.\n");
        }
    }
    myPrint(options, 2, "Number of sequences read: ", length(originalSeqs),
            "\nLongest sequence read: ", maxLen, "\n");

    if (length(originalSeqs) * 6 >= std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max())
    {
        throw std::runtime_error(std::string("ERROR: Too many sequences submitted. The maximum (including frames) is ")
                                 + std::to_string(std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max()) + ".\n");
    }

    if (maxLen >= std::numeric_limits<SizeTypePos_<TOrigAlph>>::max())
    {
        throw std::runtime_error(std::string("ERROR: one or more of your subject sequences are too long. "
                  "The maximum length is ") + std::to_string(std::numeric_limits<SizeTypePos_<TOrigAlph>>::max()) +
                  ".\n");
    }

    if (options.hasSTaxIds)
    {
        myPrint(options, 2, "Subjects without acc numbers:             ", noAcc, '/', length(ids), "\n",
                            "Subjects with more than one acc number:   ", multiAcc, '/', length(ids), "\n");
    }

    myPrint(options, 2, "\n");


    myPrint(options, 1, "Dumping Subj Ids...");

    //TODO save to TMPDIR instead
    CharString _path = options.indexDir;
    append(_path, "/seq_ids");
    save(ids, toCString(_path));

    myPrint(options, 1, " done.\n");
    finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TLimits>
inline void
_saveOriginalSeqLengths(TLimits limits, // we want copy!
                       LambdaIndexerOptions const & options)
{
    double start = sysTime();
    myPrint(options, 1, "Dumping untranslated subject lengths...");

    for (uint32_t i = 0; i < (length(limits) - 1); ++i)
        limits[i] = limits[i+1] - limits[i];
    // last entry not overwritten, should be the sum of all lengths

    CharString _path = options.indexDir;
    append(_path, "/untranslated_seq_lengths");
    save(limits, toCString(_path));
    myPrint(options, 1, " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TOrigAlph>
inline void
translateOrSwap(TCDStringSet<String<TTransAlph>> & out,
                TCDStringSet<String<TOrigAlph>> & in,
                LambdaIndexerOptions const & options)
{
    double start = sysTime();
    myPrint(options, 1, "Translating Subj Sequences...");

    translate(out,
              in,
              SIX_FRAME,
              options.geneticCode);

    myPrint(options, 1, " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

template <typename TSameAlph>
inline void
translateOrSwap(TCDStringSet<String<TSameAlph>> & out,
                TCDStringSet<String<TSameAlph>> & in,
                LambdaIndexerOptions const & /**/)
{
    swap(out, in);
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph>
inline void
dumpTranslatedSeqs(TCDStringSet<String<TTransAlph>> const & translatedSeqs,
                   LambdaIndexerOptions const & options)
{
    double start = sysTime();
    myPrint(options, 1, "Dumping unreduced Subj Sequences...");

    //TODO save to TMPDIR instead
    std::string _path = options.indexDir + "/translated_seqs";
    save(translatedSeqs, _path.c_str());

    myPrint(options, 1, " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

// template <typename TTransAlph, typename TRedAlph>
// inline void
// reduceOrSwap(TCDStringSet<String<TRedAlph>> & out,
//              TCDStringSet<String<TTransAlph>> & in)
// {
//     //TODO more output
//     // reduce implicitly
//     myPrint(options, 1, "Reducing...");
//     out.concat = in.concat;
//     out.limits = in.limits;
// }
//
// template <typename TSameAlph>
// inline void
// reduceOrSwap(TCDStringSet<String<TSameAlph>> & out,
//              TCDStringSet<String<TSameAlph>> & in)
// {
//     swap(out, in);
// }

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TRedAlph, BlastProgram p>
void
checkIndexSize(TCDStringSet<String<TRedAlph>> const & seqs,
               LambdaIndexerOptions const & options,
               BlastProgramSelector<p> const &)
{
    myPrint(options, 1, "Checking parameters of to-be-built index...");

    // check number of sequences
    using SAV = typename SAValue<TCDStringSet<String<TRedAlph>>>::Type;
    uint64_t curNumSeq = length(seqs);
    uint64_t maxNumSeq = std::numeric_limits<typename Value<SAV, 1>::Type>::max();

    if (curNumSeq >= maxNumSeq)
    {
        throw std::invalid_argument(std::string("ERROR: Too many sequences to be indexed:\n  ") +
                                    std::to_string(length(seqs)) +
                                    std::string(" in file, but only ") +
                                    std::to_string(maxNumSeq) +
                                    std::string(" supported by index.\n"));
    }

    // check length of sequences
    uint64_t maxLenSeq = std::numeric_limits<typename Value<SAV, 2>::Type>::max();
    uint64_t maxLen = 0ul;
    for (auto const & s : seqs)
        if (length(s) > maxLen)
            maxLen = length(s);

    if (maxLen >= maxLenSeq)
    {
        std::string err;
        err += "Sequences too long to be indexed:\n  ";
        err += "length";
        err += std::to_string(maxLen);
        err += " present in file, but only ";
        err += std::to_string(maxLenSeq);
        err += " supported by index.\n";
        #ifndef LAMBDA_LONG_PROTEIN_SUBJ_SEQS
        if (p != BlastProgram::BLASTN)
            err += "You can recompile Lambda and add -DLAMBDA_LONG_PROTEIN_SUBJ_SEQS=1 to activate\n"
                   "support for longer protein sequences.\n";
        #endif

        throw std::invalid_argument(err);
    }

    // check available RAM
    auto ram = getTotalSystemMemory();
    auto lS = lengthSum(seqs);
    unsigned long long factor = 0;
    if (options.algo == "radixsort")
        factor = sizeof(SizeTypeNum_<TRedAlph>) + sizeof(SizeTypePos_<TRedAlph>) + 4; // 4 is good heuristic
    else if (options.algo == "skew7ext")
        factor = 6; // TODO do some tests!
    auto estimatedSize = lS * factor;

    myPrint(options, 1, "done.\n");
    if (estimatedSize >= ram)
    {
        std::cerr << "WARNING: Lambda estimates that it will need " << estimatedSize / 1024 / 1024 << "MB\n"
                  << "         of memory to index this file, but you have only " << ram / 1024 / 1024 << "MB\n"
                  << "         available on your system.\n"
                  << "         This means you will likely encounter a crash with \"bad_alloc\".\n"
                  << "         Split you sequence file into many smaller ones or use a computer\n"
                  << "         with more memory!\n";
    } else
    {
        myPrint(options, 2, "Detected RAM: ", ram / 1024 / 1024, "MB, Estimated RAM usage: ",
                estimatedSize / 1024 / 1024, "MB\n\n");
    }
}

// --------------------------------------------------------------------------
// Function mapAndDumpTaxIDs()
// --------------------------------------------------------------------------

void
mapAndDumpTaxIDs(std::vector<bool>                                     & taxIdIsPresent,
                 std::unordered_map<std::string, uint64_t>       const & accToIdRank,
                 uint64_t                                        const   numSubjects,
                 LambdaIndexerOptions                            const & options)

{
    StringSet<String<uint32_t>> sTaxIds; // not concat because we resize inbetween
    resize(sTaxIds, numSubjects);

    // c++ stream
    std::ifstream fin(toCString(options.accToTaxMapFile), std::ios_base::in | std::ios_base::binary);
    if (!fin.is_open())
    {
        throw std::invalid_argument(std::string("ERROR: Could not open acc-to-tax-map file at ") +
                                    options.accToTaxMapFile + "\n");
    }

    // transparent decompressor
    VirtualStream<char, Input> vfin {fin};
    // stream iterator
    auto fit = directionIterator(vfin, Input());

    myPrint(options, 1, "Parsing acc-to-tax-map file... ");

    double start = sysTime();

    if (std::regex_match(options.accToTaxMapFile, std::regex{R"raw(.*\.accession2taxid(\.(gz|bgzf|bz2))?)raw"}))
    {
        _readMappingFileNCBI(fit, sTaxIds, taxIdIsPresent, accToIdRank);
    } else if (std::regex_match(options.accToTaxMapFile, std::regex{R"raw(.*\.dat(\.(gz|bgzf|bz2))?)raw"}))
    {
        _readMappingFileUniProt(fit, sTaxIds, taxIdIsPresent, accToIdRank);
    } else
    {
        throw std::invalid_argument("ERROR: extension of acc-to-tax-map file not handled.\n");
    }

    // root node is always present
    if (taxIdIsPresent.size() < 2)
        taxIdIsPresent.resize(2);
    taxIdIsPresent[1] = true;

    myPrint(options, 1, "done.\n");

    myPrint(options, 2, "Runtime: ", sysTime() - start, "s \n");

    // TODO do something with the subjects that have no (valid) taxid?

    uint64_t nomap = 0;
    uint64_t multi = 0;

    for (auto const & s : sTaxIds)
    {
        if (length(s) == 0)
            ++nomap;
        else if (length(s) > 1)
            ++multi;
    }

    myPrint(options, 2, "Subjects without tax IDs:             ", nomap, '/', numSubjects, "\n",
                        "Subjects with more than one tax ID:   ", multi, '/', numSubjects, "\n\n");
    if ((nomap > 0) && ((numSubjects / nomap) < 5))
        myPrint(options, 1, "WARNING: ", double(nomap) * 100 / numSubjects, "% of subjects have no taxID.\n"
                            "         Maybe you specified the wrong map file?\n\n");

    myPrint(options, 1,"Dumping Subject Taxonomy IDs... ");
    start = sysTime();
    // concat direct so that it's easier to read/write
    StringSet<String<uint32_t>, Owner<ConcatDirect<>>> outSTaxIds = sTaxIds;
    save(outSTaxIds, std::string(options.indexDir + "/staxids").c_str());
    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");
}

// --------------------------------------------------------------------------
// Function mapAndDumpTaxIDs()
// --------------------------------------------------------------------------

void
parseAndDumpTaxTree(std::vector<bool>          & taxIdIsPresent,
                    LambdaIndexerOptions const & options)

{
    String<uint32_t> taxonParentIDs; // ever position has the index of its parent node
    reserve(taxonParentIDs, 2'000'000); // reserve 2million to save reallocs

    std::string path = options.taxDumpDir + "/nodes.dmp";

    std::ifstream fin(path.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin.is_open())
    {
        throw std::runtime_error(std::string("ERROR: Could not open ") + path + "\n");
    }

    // transparent decompressor
    VirtualStream<char, Input> vfin{fin};
    // stream iterator
    auto fit = directionIterator(vfin, Input());

    myPrint(options, 1, "Parsing nodes.dmp... ");

    double start = sysTime();

    std::string buf;
    std::regex const numRegEx{"\\b\\d+\\b"};

    while (!atEnd(fit))
    {
        clear(buf);
        // read line
        readLine(buf, fit);

        uint32_t n = 0;
        uint32_t parent = 0;
        unsigned i = 0;
        for (auto it = std::sregex_iterator(buf.begin(), buf.end(), numRegEx), itEnd = std::sregex_iterator();
             (it != itEnd) && (i < 2);
             ++it, ++i)
        {
            try
            {
                if (i == 0)
                    n = lexicalCast<uint32_t>(it->str());
                else
                    parent = lexicalCast<uint32_t>(it->str());
            }
            catch (BadLexicalCast const & badCast)
            {
                throw std::runtime_error(
                    std::string("Error: Expected taxonomical ID, but got something I couldn't read: ") +
                    std::string(badCast.what()) + "\n");
            }
        }
        if (length(taxonParentIDs) <= n)
            resize(taxonParentIDs, n +1, 0);
        taxonParentIDs[n] = parent;
    }
    // also resize these, since we get new, possibly higher cardinality nodes
    taxIdIsPresent.resize(length(taxonParentIDs), false);

    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n");

    if (options.verbosity >= 2)
    {
        uint32_t heightMax = 0;
        uint32_t numNodes = 0;
        for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
        {
            if (taxonParentIDs[i] > 0)
                ++numNodes;
            uint32_t height = 0;
            uint32_t curPar = taxonParentIDs[i];
            while (curPar > 1)
            {
                curPar = taxonParentIDs[curPar];
                ++height;
            }

            heightMax = std::max(heightMax, height);
        }
        myPrint(options, 2, "Number of nodes in tree: ", numNodes, "\n");
        myPrint(options, 2, "Maximum Tree Height: ", heightMax, "\n\n");
    }

    myPrint(options, 1, "Thinning and flattening Tree... ");
    start = sysTime();

    // taxIdIsPresent are all directly present taxIds
    // taxIdIsPresentOrParent are also the (recursive) parents of the above
    // we need to differentiate this later, because we will remove some intermediate nodes
    // but we may not remove any that are directly present AND parents of directly present ones
    std::vector<bool> taxIdIsPresentOrParent{taxIdIsPresent};
    // mark parents as present, too
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
    {
        if (taxIdIsPresent[i])
        {
            // get ancestors:
            uint32_t curPar = i;
            do
            {
                curPar = taxonParentIDs[curPar];
                taxIdIsPresentOrParent[curPar] = true;
            } while (curPar > 1);
        }
    }

    // set unpresent nodes to 0
    SEQAN_OMP_PRAGMA(parallel for)
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
        if (!taxIdIsPresentOrParent[i])
            taxonParentIDs[i] = 0;

    // count inDegrees
    String<uint32_t> inDegrees;
    resize(inDegrees, length(taxonParentIDs), 0);
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
    {
        // increase inDegree of parent
        uint32_t curPar = taxonParentIDs[i];
        ++inDegrees[curPar];
    }

    // skip parents with indegree 1 (flattening)
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
    {
        uint32_t curPar = taxonParentIDs[i];
        // those intermediate nodes that themselve represent sequences may not be skipped
        while ((curPar > 1) && (inDegrees[curPar] == 1) && (!taxIdIsPresent[curPar]))
            curPar = taxonParentIDs[curPar];

        taxonParentIDs[i] = curPar;
    }

    // remove nodes that are now disconnected
    SEQAN_OMP_PRAGMA(parallel for)
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
    {
        // those intermediate nodes that themselve represent sequences may not be skipped
        if ((inDegrees[i] == 1) && (!taxIdIsPresent[i]))
        {
            taxonParentIDs[i] = 0;
            taxIdIsPresentOrParent[i] = false;
        }
    }

    String<uint8_t> taxonHeights;
    resize(taxonHeights, length(taxonParentIDs), 0);

    {
        uint32_t heightMax = 0;
        uint32_t numNodes = 0;
        for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
        {
            if (taxonParentIDs[i] > 0)
                ++numNodes;
            uint32_t height = 0;
            uint32_t curPar = taxonParentIDs[i];
            while (curPar > 1)
            {
                curPar = taxonParentIDs[curPar];
                ++height;
            }
            taxonHeights[i] = height;
            heightMax = std::max(heightMax, height);
        }

        myPrint(options, 1, "done.\n");
        myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n");

        myPrint(options, 2, "Number of nodes in tree: ", numNodes, "\n");
        myPrint(options, 2, "Maximum Tree Height: ", heightMax, "\n\n");
    }

    myPrint(options, 1,"Dumping Taxonomy Tree... ");
    start = sysTime();
    save(taxonParentIDs, std::string(options.indexDir + "/tax_parents").c_str());
    save(taxonHeights,   std::string(options.indexDir + "/tax_heights").c_str());
    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");

    // DEBUG
    #ifndef NDEBUG
    for (uint32_t i = 0; i < length(taxonParentIDs); ++i)
    {
        if (!taxIdIsPresentOrParent[i] && (taxonParentIDs[i] != 0))
            std::cerr << "WARNING: TaxID " << i << " has parent, but shouldn't.\n";

        if (taxIdIsPresentOrParent[i] && (taxonParentIDs[i] == 0))
            std::cerr << "WARNING: TaxID " << i << " has no parent, but should.\n";
        if (taxIdIsPresent[i] && (taxonParentIDs[i] == 0))
            std::cerr << "WARNING: TaxID " << i << " has no parent, but should. 2\n";

        if (taxIdIsPresent[i] && !taxIdIsPresentOrParent[i])
            std::cerr << "WARNING: TaxID " << i << " disappeared, but shouldn't have.\n";

        if (!taxIdIsPresent[i] && taxIdIsPresentOrParent[i] && (inDegrees[i] == 1))
            std::cerr << "WARNING: TaxID " << i << " should have disappeared, but didn't.\n";
    }
    #endif

    /** read the names **/

    StringSet<CharString> taxonNames; // ever position has the index of its parent node
    resize(taxonNames, length(taxonParentIDs));

    path = options.taxDumpDir + "/names.dmp";

    std::ifstream fin2(path.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin2.is_open())
        throw std::runtime_error(std::string("ERROR: Could not open ") + path + "\n");

    // transparent decompressor
    VirtualStream<char, Input> vfin2{fin2};
    // stream iterator
    fit = directionIterator(vfin2, Input());

    myPrint(options, 1, "Parsing names.dmp... ");

    start = sysTime();

    std::regex const wordRegEx{R"([\w.,\"<> ]+)"};
    std::string name;

    while (!atEnd(fit))
    {
        clear(buf);
        // read line
        readLine(buf, fit);

        uint32_t taxId = 0;

        auto itWord = std::sregex_iterator(buf.begin(), buf.end(), wordRegEx);
        if (itWord == std::sregex_iterator())
        {
            throw std::runtime_error("Error: Expected taxonomical ID in first column, but couldn't find it.\n");
        } else
        {
            try
            {
                taxId = lexicalCast<uint64_t>(itWord->str());
            }
            catch (BadLexicalCast const & badCast)
            {
                throw std::runtime_error(std::string("Error: Expected taxonomical ID in first column, but got something"
                                                     " I couldn't read: ") + std::string(badCast.what()) + "\n");
            }

            if (taxId >= length(taxonNames))
            {
                throw std::runtime_error(std::string("Error: taxonomical ID is ") + std::to_string(taxId) +
                                         ", but no such taxon in tree.\n");
            }
        }

        // we don't need this name
        if (!taxIdIsPresentOrParent[taxId])
            continue;

        if (++itWord == std::sregex_iterator())
            throw std::runtime_error("Error: Expected name in second column, but couldn't find it.\n");
        else
            name = itWord->str();

        while (++itWord != std::sregex_iterator())
        {
            if (itWord->str() == "scientific name")
                taxonNames[taxId] = name;
        }
    }

    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n");

    taxonNames[0] = "invalid";
    for (uint32_t i = 0; i < length(taxonNames); ++i)
    {
        if (taxIdIsPresentOrParent[i] && empty(taxonNames[i]))
        {
            std::cerr << "Warning: Taxon with ID " << i << " has no name associated, defaulting to \"n/a\".\n";
            taxonNames[i] = "n/a";
        }
    }

    myPrint(options, 1,"Dumping Taxon names... ");
    start = sysTime();
    // concat direct so that it's easier to read/write
    StringSet<CharString, Owner<ConcatDirect<>>> outTaxonNames = taxonNames;
    save(outTaxonNames, std::string(options.indexDir + "/tax_names").c_str());
    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");
}

// --------------------------------------------------------------------------
// Function createSuffixArray()
// --------------------------------------------------------------------------

// If there is no overload with progress function, then strip it
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgo,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  TAlgo const &,
                  TLambda &&)
{
    return createSuffixArray(SA, s, TAlgo());
}

// ----------------------------------------------------------------------------
// Function indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
void
indexCreateProgress(Index<TText, FMIndex<TSpec, TConfig> > & index,
                    FibreSALF const &,
                    LambdaIndexerOptions const & options)
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >               TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type            TTempSA;
    typedef typename Size<TIndex>::Type                          TSize;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    TText const & text = indexText(index);

    if (empty(text))
        return;

    TTempSA tempSA;
    uint64_t lastPercent = 0;

    double s = sysTime();
    myPrint(options, 1, "Generating Index 0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n"
                        " Progress:       |");
    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    if (options.verbosity >= 1)
    {
        createSuffixArray(tempSA,
                          text,
                          TAlgo(),
                          [&lastPercent] (uint64_t curPerc)
                          {
                              // needs locking, because called from multiple threads
                              SEQAN_OMP_PRAGMA(critical(progressBar))
                              printProgressBar(lastPercent, curPerc * 0.85); // 85% of progress
                          });
    } else
    {
        createSuffixArray(tempSA,
                          text,
                          TAlgo());
    }
    double sacaTime = sysTime() - s;

    if (options.verbosity >= 1)
        printProgressBar(lastPercent, 85);

    // Create the LF table.
    s = sysTime();
    if (options.verbosity >= 1)
    {
        createLFProgress(indexLF(index),
                         text,
                         tempSA,
                         [&lastPercent] (uint64_t curPerc)
                         {
                             // doesn't need locking, only writes from one thread
                             printProgressBar(lastPercent, curPerc * 0.1); // 10% of progress
                         });
    } else
    {
        createLFProgress(indexLF(index),
                         text,
                         tempSA,
                         [] (uint64_t) {});
    }
    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());
    double bwtTime = sysTime() - s;

    if (options.verbosity >= 1)
        printProgressBar(lastPercent, 95);

    // Create the sampled SA.
    s = sysTime();
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);
    double sampleTime = sysTime() - s;

    if (options.verbosity >= 1)
        printProgressBar(lastPercent, 100);

    myPrint(options, 1, "\n");
    myPrint(options, 2, "SA  construction runtime: ", sacaTime, "s\n");
    myPrint(options, 2, "BWT construction runtime: ", bwtTime, "s\n");
    myPrint(options, 2, "SA  sampling runtime:     ", sampleTime, "s\n");
    myPrint(options, 1, "\n");
}

template <typename TText, typename TSpec, typename TConfig>
void
indexCreateProgress(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > & index,
                    FibreSALF const &,
                    LambdaIndexerOptions const & options)
{
    myPrint(options, 1, "Bi-Directional Index [forward]\n");
    indexCreateProgress(index.fwd, FibreSALF(), options);

    myPrint(options, 1, "Bi-Directional Index [backward]\n");
    indexCreateProgress(index.rev, FibreSALF(), options);
}

template <typename TText, typename TSpec>
void
indexCreateProgress(Index<TText, IndexSa<TSpec> > & index,
                    FibreSA const &,
                    LambdaIndexerOptions const & options)
{
    typedef Index<TText, IndexSa<TSpec> >                        TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                TSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    TText const & text = indexText(index);

    if (empty(text))
        return;

    TSA & sa = getFibre(index, FibreSA());

    myPrint(options, 1, "Generating Index 0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n"
                        "  Progress:      |");
    // Create the full SA.
    resize(sa, lengthSum(text), Exact());
    if (options.verbosity >= 1)
    {
        createSuffixArray(sa,
                          text,
                          TAlgo(),
                          [lastPercent = uint64_t{0ull}] (uint64_t curPerc) mutable
                          {
                              SEQAN_OMP_PRAGMA(critical(progressBar))
                              printProgressBar(lastPercent, curPerc); // 100% of progress
                          });
    } else
    {
        createSuffixArray(sa,
                          text,
                          TAlgo());
    }
}

template <typename T>
inline void
_clearSparseSuffixArray(T &, std::false_type const &)
{}

template <typename T>
inline void
_clearSparseSuffixArray(T & dbIndex, std::true_type const &)
{
    // reverse index does not require sampled suffix array, but its size :|
    clear(getFibre(getFibre(getFibre(dbIndex, FibreSA()), FibreSparseString()), FibreValues()));
    clear(getFibre(getFibre(getFibre(dbIndex, FibreSA()), FibreSparseString()), FibreIndicators()));
}


// --------------------------------------------------------------------------
// Function generateIndexAndDump()
// --------------------------------------------------------------------------

#ifdef _OPENMP
#define TID omp_get_thread_num()
#else
#define TID 0
#endif

template <typename TIndexSpec,
          typename TIndexSpecSpec,
          typename TString,
          typename TSpec,
          typename TRedAlph_,
          typename TDirection,
          BlastProgram p>
inline void
generateIndexAndDump(StringSet<TString, TSpec>        & seqs,
                     LambdaIndexerOptions       const & options,
                     BlastProgramSelector<p>    const &,
                     TRedAlph_                  const &,
                     Tag<TDirection>            const &)
{
    using TTransSeqs    = TCDStringSet<String<TransAlph<p>>>;

    using TRedAlph      = RedAlph<p, TRedAlph_>; // ensures == Dna5 for BlastN
    using TRedSeqVirt   = ModifiedString<String<TransAlph<p>, Alloc<>>,
                            ModView<FunctorConvert<TransAlph<p>,TRedAlph>>>;
    using TRedSeqsVirt  = StringSet<TRedSeqVirt, Owner<ConcatDirect<>>>;

    static bool constexpr
    indexIsFM           = std::is_same<TIndexSpec, TFMIndex<TIndexSpecSpec> >::value ||
                          std::is_same<TIndexSpec, TFMIndexInBi<TIndexSpecSpec> >::value;
    static bool constexpr
    alphReduction       = !std::is_same<TransAlph<p>, TRedAlph>::value;

    using TRedSeqs      = typename std::conditional<
                            !alphReduction,
                            TTransSeqs,             // owner
                            TRedSeqsVirt>::type;    // modview
    using TRedSeqsACT   = typename std::conditional<
                            !alphReduction,
                            TTransSeqs &,           // reference to owner
                            TRedSeqsVirt>::type;    // modview

    using TDbIndex      = Index<TRedSeqs, TIndexSpec>;
    using TFullFibre    = typename std::conditional<indexIsFM,
                                                    FibreSALF,
                                                    FibreSA>::type;
    static bool constexpr
    hasProgress         = std::is_same<TIndexSpecSpec, RadixSortSACreateTag>::value;

    // Generate Index
    if (!hasProgress)
        myPrint(options, 1, "Generating Index...");

    double s = sysTime();

//     std::cerr << "indexIsFM: " << int(indexIsFM) << std::endl;

    // FM-Index needs reverse input
    if (indexIsFM && std::is_same<Tag<TDirection>, Fwd>::value)
        reverse(seqs);

    TRedSeqsACT redSubjSeqs(seqs);

    TDbIndex dbIndex(redSubjSeqs);

    // instantiate SA
    indexCreateProgress(dbIndex, TFullFibre(),  options);

    // since we dumped unreduced sequences before and reduced sequences are
    // only "virtual" we clear them before dump
    std::decay_t<decltype(redSubjSeqs.limits)> tmpLimits;
    if (alphReduction || !indexIsFM) // fm indexes don't dump them anyways
    {
        if (indexIsFM && (std::is_same<Tag<TDirection>, Rev>::value))
        {
            // these makes redSubjSeqs appear empty and deactivates output
            swap(tmpLimits, redSubjSeqs.limits);

            _clearSparseSuffixArray(dbIndex, std::integral_constant<bool, indexIsFM>{});
        } else
        {
            clear(seqs);
            clear(redSubjSeqs.limits); // limits part is not lightweight
        }
    }

    double e = sysTime() - s;
    if (!hasProgress)
    {
        myPrint(options, 1, " done.\n");
        myPrint(options, 2, "Runtime: ", e, "s \n\n");
    }

    // Dump Index
    myPrint(options, 1, "Writing Index to disk...");
    s = sysTime();
    std::string path = options.indexDir + "/index";

    if (std::is_same<Tag<TDirection>, Rev>::value)
        path += ".rev";

    save(dbIndex, path.c_str());

    e = sysTime() - s;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n");

    if (alphReduction && indexIsFM && (std::is_same<Tag<TDirection>, Rev>::value))
    {
        // we swap back so that the sequences can be used for building the second index
        swap(tmpLimits, redSubjSeqs.limits);
//         redSubjSeqs.limits = tmpLimits;
    }
}

#endif // header guard
