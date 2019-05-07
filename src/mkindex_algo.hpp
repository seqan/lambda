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

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/translation.hpp>
#include <seqan3/std/charconv>
#include <seqan3/std/concepts>

#include "mkindex_misc.hpp"
// #include "mkindex_saca.hpp"
#include "shared_misc.hpp"
#include "shared_options.hpp"

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TOrigAlph>
auto
loadSubjSeqsAndIds(LambdaIndexerOptions const & options)
{
    using TIDs          = TCDStringSet<std::string>;
    using TOrigSeqs     = TCDStringSet<std::vector<TOrigAlph>>;
    using TAccToIdRank  = std::unordered_map<std::string, uint64_t>;

    std::tuple<TIDs, TOrigSeqs, TAccToIdRank> ret;

    auto & ids          = std::get<0>(ret);
    auto & originalSeqs = std::get<1>(ret);
    auto & accToIdRank  = std::get<2>(ret);

    // Make sure we have enough RAM to load the file
    auto ram = getTotalSystemMemory();
    auto fS = fileSize(options.dbFile.c_str());

    if (fS >= ram)
        std::cerr << "WARNING: Your sequence file is already larger than your physical memory!\n"
                  << "         This means you will likely encounter a crash with \"bad_alloc\".\n"
                  << "         Split you sequence file into many smaller ones or use a computer\n"
                  << "         with more memory!\n";






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

    uint64_t noAcc = 0;
    uint64_t multiAcc = 0;

    // lambda that extracts accession numbers and saves them in the map
    auto extractAccIds = [&accToIdRank, &accRegEx, &noAcc, &multiAcc] (auto && id, uint64_t const rank)
    {

        std::conditional_t<(bool)std::Constructible<std::string const &, decltype(id)>, std::string const &, std::string>
        buf{id};

        uint64_t count = 0;
        for (auto it = std::sregex_iterator(buf.begin(), buf.end(), accRegEx), itEnd = std::sregex_iterator();
             it != itEnd;
             ++it, ++count)
        {
            assert(accToIdRank.count(it->str()) == 0);
//                              "An accession number appeared twice in the file, but they should be unique.");

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

    using seq_traits = std::conditional_t<seqan3::NucleotideAlphabet<TOrigAlph>,
                                          seqan3::sequence_file_input_default_traits_dna,
                                          seqan3::sequence_file_input_default_traits_aa>;
    seqan3::sequence_file_input<seq_traits, seqan3::fields<seqan3::field::ID, seqan3::field::SEQ>>
        infile{options.dbFile};

    size_t count = 0;
    for (auto & [ id, seq ] : infile)
    {
        if (options.hasSTaxIds)
            extractAccIds(id, count);

        if (options.truncateIDs)
            ids.push_back(id | seqan3::view::take_until(seqan3::is_space));
        else
            ids.push_back(std::move(id));

        originalSeqs.push_back(std::move(seq));
        ++count;
    }


    myPrint(options, 1,  " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n");

    if (std::ranges::empty(originalSeqs))
    {
        throw std::runtime_error("ERROR: No sequences in file. Aborting.\n");
    }

    size_t maxLen = 0ul;
    for (auto const & s : originalSeqs)
    {
        if (std::ranges::size(s) > maxLen)
        {
            maxLen = std::ranges::size(s);
        }
        else if (std::ranges::size(s) == 0ul)
        {
            throw std::runtime_error("ERROR: Unexpectedly encountered a sequence of length 0 in the file."
                                     "Remove the entry and try again. Aborting.\n");
        }
    }
    myPrint(options, 2, "Number of sequences read: ", std::ranges::size(originalSeqs),
            "\nLongest sequence read: ", maxLen, "\n");

//     if (std::ranges::size(originalSeqs) * 6 >= std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max())
//     {
//         throw std::runtime_error(std::string("ERROR: Too many sequences submitted. The maximum (including frames) is ")
//                                  + std::to_string(std::numeric_limits<SizeTypeNum_<TOrigAlph>>::max()) + ".\n");
//     }
//
//     if (maxLen >= std::numeric_limits<SizeTypePos_<TOrigAlph>>::max())
//     {
//         throw std::runtime_error(std::string("ERROR: one or more of your subject sequences are too long. "
//                   "The maximum length is ") + std::to_string(std::numeric_limits<SizeTypePos_<TOrigAlph>>::max()) +
//                   ".\n");
//     }

    if (options.hasSTaxIds)
    {
        myPrint(options, 2, "Subjects without acc numbers:             ", noAcc, '/', std::ranges::size(ids), "\n",
                            "Subjects with more than one acc number:   ", multiAcc, '/', std::ranges::size(ids), "\n");
    }

    myPrint(options, 2, "\n");

    return ret;
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TSeqSet>
std::vector<uint64_t>
saveOriginalSeqLengths(TSeqSet const & seqSet,
                       LambdaIndexerOptions const & /**/)
{
//     double start = sysTime();
//     myPrint(options, 1, "Dumping untranslated subject lengths...");

    std::vector<uint64_t> limits;
    limits.resize(std::ranges::size(seqSet) + 1); // last holds sum

    uint64_t sum = 0;
    for (size_t i = 0; i < limits.size() - 1; ++i)
    {
        limits[i] = std::ranges::size(seqSet[i]);
        sum += std::ranges::size(seqSet[i]);
    }
    limits.back() = sum;

    return limits;

//     std::string _path = options.indexDir + "/untranslated_seq_lengths";
//
//     {
//         std::ofstream os{_path.c_str()};
//
//         cereal::BinaryOutputArchive oarchive(os); // Create an output archive
//         oarchive(limits);
//     }
//
//     myPrint(options, 1, " done.\n");
//     double finish = sysTime() - start;
//     myPrint(options, 2, "Runtime: ", finish, "s \n\n");
}

// --------------------------------------------------------------------------
// Function loadSubj()
// --------------------------------------------------------------------------

template <typename TTransAlph, typename TOrigAlph>
TCDStringSet<std::vector<TTransAlph>>
translateSeqs(TCDStringSet<std::vector<TOrigAlph>> & in,
              LambdaIndexerOptions const & options)
{
    TCDStringSet<std::vector<TTransAlph>> out;

    double start = sysTime();
    myPrint(options, 1, "Translating Subj Sequences...");

    out = in | seqan3::view::translate | std::view::join; //TODO geneticCode

    myPrint(options, 1, " done.\n");
    double finish = sysTime() - start;
    myPrint(options, 2, "Runtime: ", finish, "s \n\n");

    return out;
}

template <typename TTransAlph, typename TOrigAlph>
    requires std::Same<TTransAlph, TOrigAlph>
TCDStringSet<std::vector<TTransAlph>>
translateSeqs(TCDStringSet<std::vector<TOrigAlph>> & in,
              LambdaIndexerOptions const & /**/)
{
    return std::move(in);
}

// --------------------------------------------------------------------------
// Function checkIndexSize()
// --------------------------------------------------------------------------

template <typename TRedAlph>
void
checkIndexSize(TCDStringSet<std::vector<TRedAlph>> const & seqs,
               LambdaIndexerOptions const & options)
{
#if 0 // TODO: new index should handle arbitrary sizes, but we want the RAM check again
    myPrint(options, 1, "Checking parameters of to-be-built index...");

    // check number of sequences
    using SAV = typename SAValue<TCDStringSet<std::vector<TRedAlph>>>::Type;
    uint64_t curNumSeq = std::ranges::size(seqs);
    uint64_t maxNumSeq = std::numeric_limits<typename Value<SAV, 1>::Type>::max();

    if (curNumSeq >= maxNumSeq)
    {
        throw std::invalid_argument(std::string("ERROR: Too many sequences to be indexed:\n  ") +
                                    std::to_string(std::ranges::size(seqs)) +
                                    std::string(" in file, but only ") +
                                    std::to_string(maxNumSeq) +
                                    std::string(" supported by index.\n"));
    }

    // check length of sequences
    uint64_t maxLenSeq = std::numeric_limits<typename Value<SAV, 2>::Type>::max();
    uint64_t maxLen = 0ul;
    for (auto const & s : seqs)
        if (std::ranges::size(s) > maxLen)
            maxLen = std::ranges::size(s);

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
#else
    (void)seqs;
    (void)options;
#endif
}

// --------------------------------------------------------------------------
// Function mapAndDumpTaxIDs()
// --------------------------------------------------------------------------

auto mapTaxIDs(std::unordered_map<std::string, uint64_t>       const & accToIdRank,
               uint64_t                                        const   numSubjects,
               LambdaIndexerOptions                            const & options)

{
    using TTaxIds        = std::vector<std::vector<uint32_t>>;// not concat because we resize inbetween
    using TTaxIdsPresent = std::vector<bool>;

    std::tuple<TTaxIds, TTaxIdsPresent> ret;

    auto & sTaxIds          = std::get<0>(ret);
    auto & taxIdIsPresent   = std::get<1>(ret);

    taxIdIsPresent.reserve(2'000'000);
    sTaxIds.resize(numSubjects);

    // c++ stream
    std::ifstream fin(options.accToTaxMapFile.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin.is_open())
    {
        throw std::invalid_argument(std::string("ERROR: Could not open acc-to-tax-map file at ") +
                                    options.accToTaxMapFile + "\n");
    }

    // transparent decompressor
    auto vstream = seqan3::detail::make_secondary_istream(fin);

    auto file_view = std::ranges::subrange<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
    {
        std::istreambuf_iterator<char>{*vstream},
        std::istreambuf_iterator<char>{}
    };

    myPrint(options, 1, "Parsing acc-to-tax-map file... ");

    double start = sysTime();

    if (std::regex_match(options.accToTaxMapFile, std::regex{R"raw(.*\.accession2taxid(\.(gz|bgzf|bz2))?)raw"}))
    {
        _readMappingFileNCBI(file_view, sTaxIds, taxIdIsPresent, accToIdRank);
    } else if (std::regex_match(options.accToTaxMapFile, std::regex{R"raw(.*\.dat(\.(gz|bgzf|bz2))?)raw"}))
    {
        _readMappingFileUniProt(file_view, sTaxIds, taxIdIsPresent, accToIdRank);
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

    for (auto && s : sTaxIds)
    {
        if (std::ranges::size(s) == 0)
            ++nomap;
        else if (std::ranges::size(s) > 1)
            ++multi;
    }

    myPrint(options, 2, "Subjects without tax IDs:             ", nomap, '/', numSubjects, "\n",
                        "Subjects with more than one tax ID:   ", multi, '/', numSubjects, "\n\n");
    if ((nomap > 0) && ((numSubjects / nomap) < 5))
        myPrint(options, 1, "WARNING: ", double(nomap) * 100 / numSubjects, "% of subjects have no taxID.\n"
                            "         Maybe you specified the wrong map file?\n\n");

    return ret;

//     myPrint(options, 1,"Dumping Subject Taxonomy IDs... ");
//     start = sysTime();
//
//     std::string _path = options.indexDir + "/staxids";
//
//     {
//         std::ofstream os{_path.c_str()};
//
//         cereal::BinaryOutputArchive oarchive(os); // Create an output archive
//         oarchive(sTaxIds);
//     }
//     myPrint(options, 1, "done.\n");
//     myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");
}

// --------------------------------------------------------------------------
// Function mapAndDumpTaxIDs()
// --------------------------------------------------------------------------

auto parseAndStoreTaxTree(std::vector<bool>          & taxIdIsPresent,
                          LambdaIndexerOptions const & options)

{
    using TTaxonParentIDs   = std::vector<uint32_t>; // ever position has the index of its parent node
    using TTaxonHeights     = std::vector<uint8_t>;
    using TTaxonNames       = std::vector<std::string>;

    std::tuple<TTaxonParentIDs, TTaxonHeights, TTaxonNames> ret;

    auto & taxonParentIDs   = std::get<0>(ret);
    auto & taxonHeights     = std::get<1>(ret);
    auto & taxonNames       = std::get<2>(ret);

    taxonParentIDs.reserve(2'000'000); // reserve 2million to save reallocs

    std::string path = options.taxDumpDir + "/nodes.dmp";

    std::ifstream fin(path.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin.is_open())
    {
        throw std::runtime_error(std::string("ERROR: Could not open ") + path + "\n");
    }

    // transparent decompressor
    auto vstream = seqan3::detail::make_secondary_istream(fin);

    auto file_view = std::ranges::subrange<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
    {
        std::istreambuf_iterator<char>{*vstream},
        std::istreambuf_iterator<char>{}
    };

    myPrint(options, 1, "Parsing nodes.dmp... ");

    double start = sysTime();

    std::string buf;
    std::regex const numRegEx{"\\b\\d+\\b"};

    while (std::ranges::begin(file_view) != std::ranges::end(file_view))
    {
        // read line
        buf = file_view | seqan3::view::take_line;

        uint32_t n = 0;
        uint32_t parent = 0;
        unsigned i = 0;
        for (auto it = std::sregex_iterator(buf.begin(), buf.end(), numRegEx), itEnd = std::sregex_iterator();
             (it != itEnd) && (i < 2);
             ++it, ++i)
        {
            std::string strbuf = it->str();
            std::from_chars_result res;

            if (i == 0)
                res = std::from_chars(strbuf.data(), strbuf.data() + strbuf.size(), n);
            else
                res = std::from_chars(strbuf.data(), strbuf.data() + strbuf.size(), parent);

            if (res.ec != std::errc{})
            {
                throw std::runtime_error{
                    std::string{"Error: Expected taxonomical ID, but got something I couldn't read: "} +
                                strbuf + "\n"};
            }
        }
        if (std::ranges::size(taxonParentIDs) <= n)
            taxonParentIDs.resize(n +1, 0);
        taxonParentIDs[n] = parent;
    }
    // also resize these, since we get new, possibly higher cardinality nodes
    taxIdIsPresent.resize(std::ranges::size(taxonParentIDs), false);

    myPrint(options, 1, "done.\n");
    myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n");

    if (options.verbosity >= 2)
    {
        uint32_t heightMax = 0;
        uint32_t numNodes = 0;
        for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
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
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
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
//     SEQAN_OMP_PRAGMA(parallel for)
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
        if (!taxIdIsPresentOrParent[i])
            taxonParentIDs[i] = 0;

    // count inDegrees
    std::vector<uint32_t> inDegrees;
    inDegrees.resize(std::ranges::size(taxonParentIDs), 0);
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
    {
        // increase inDegree of parent
        uint32_t curPar = taxonParentIDs[i];
        ++inDegrees[curPar];
    }

    // skip parents with indegree 1 (flattening)
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
    {
        uint32_t curPar = taxonParentIDs[i];
        // those intermediate nodes that themselve represent sequences may not be skipped
        while ((curPar > 1) && (inDegrees[curPar] == 1) && (!taxIdIsPresent[curPar]))
            curPar = taxonParentIDs[curPar];

        taxonParentIDs[i] = curPar;
    }

    // remove nodes that are now disconnected
//     SEQAN_OMP_PRAGMA(parallel for)
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
    {
        // those intermediate nodes that themselve represent sequences may not be skipped
        if ((inDegrees[i] == 1) && (!taxIdIsPresent[i]))
        {
            taxonParentIDs[i] = 0;
            taxIdIsPresentOrParent[i] = false;
        }
    }

    taxonHeights.resize(std::ranges::size(taxonParentIDs), 0);

    {
        uint32_t heightMax = 0;
        uint32_t numNodes = 0;
        for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
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


//     myPrint(options, 1,"Dumping Taxonomy Tree... ");
//     start = sysTime();
//     std::string _path = options.indexDir + "/tax";
//
//     {
//         std::ofstream os{_path.c_str()};
//
//         cereal::BinaryOutputArchive oarchive(os);
//         oarchive(taxonParentIDs);
//         oarchive(taxonHeights);
//     }
//
//     myPrint(options, 1, "done.\n");
//     myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");

    // DEBUG
    #ifndef NDEBUG
    for (uint32_t i = 0; i < std::ranges::size(taxonParentIDs); ++i)
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
    taxonNames.resize(std::ranges::size(taxonParentIDs));

    path = options.taxDumpDir + "/names.dmp";

    std::ifstream fin2(path.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin2.is_open())
        throw std::runtime_error(std::string("ERROR: Could not open ") + path + "\n");

    // transparent decompressor
    auto vstream2 = seqan3::detail::make_secondary_istream(fin2);

    auto file_view2 = std::ranges::subrange<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
    {
        std::istreambuf_iterator<char>{*vstream2},
        std::istreambuf_iterator<char>{}
    };

    myPrint(options, 1, "Parsing names.dmp... ");

    start = sysTime();

    std::regex const wordRegEx{R"([\w.,\"<> ]+)"};
    std::string name;

    while (std::ranges::begin(file_view2) != std::ranges::end(file_view2))
    {
        // read line
        buf = file_view | seqan3::view::take_line;

        uint32_t taxId = 0;

        auto itWord = std::sregex_iterator(buf.begin(), buf.end(), wordRegEx);
        if (itWord == std::sregex_iterator())
        {
            throw std::runtime_error("Error: Expected taxonomical ID in first column, but couldn't find it.\n");
        } else
        {
            std::string strbuf = itWord->str();
            std::from_chars_result res;

            res = std::from_chars(strbuf.data(), strbuf.data() + strbuf.size(), taxId);

            if (res.ec != std::errc{})
            {
                throw std::runtime_error{
                    std::string{"Error: Expected taxonomical ID in first column, but got something I couldn't read: "} +
                                strbuf + "\n"};
            }

            if (taxId >= std::ranges::size(taxonNames))
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
    for (uint32_t i = 0; i < std::ranges::size(taxonNames); ++i)
    {
        if (taxIdIsPresentOrParent[i] && empty(taxonNames[i]))
        {
            std::cerr << "Warning: Taxon with ID " << i << " has no name associated, defaulting to \"n/a\".\n";
            taxonNames[i] = "n/a";
        }
    }

//     myPrint(options, 1,"Dumping Taxon names... ");
//     start = sysTime();
//     std::string _path2 = options.indexDir + "/tax_names";
//
//     {
//         std::ofstream os{_path2.c_str()};
//
//         cereal::BinaryOutputArchive oarchive(os); // Create an output archive
//         oarchive(taxonNames);
//     }
//
//     myPrint(options, 1, "done.\n");
//     myPrint(options, 2, "Runtime: ", sysTime() - start, "s\n\n");

    return ret;
}


template <bool is_bi,
          AlphabetEnum c_redAlph,
          typename TStringSet>
auto
generateIndex(TStringSet                       & seqs,
              LambdaIndexerOptions       const & options)
{
    // TODO reduced sequences should not be saved; change after index changes
    using TRedSeqs = TCDStringSet<std::vector<_alphabetEnumToType<c_redAlph>>>;
    using TIndex   = std::conditional_t<is_bi, seqan3::bi_fm_index<TRedSeqs>, seqan3::fm_index<TRedSeqs>>;

    std::tuple<TRedSeqs, TIndex> ret;

    auto & redSeqs = std::get<0>(ret);
    auto & index   = std::get<1>(ret);

    myPrint(options, 1, "Reducing sequences...");
    double s = sysTime();
    redSeqs = seqs | seqan3::view::deep{seqan3::view::convert<_alphabetEnumToType<c_redAlph>>};
    double e = sysTime() - s;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n\n");

    myPrint(options, 1, "Generating Index...");
    s = sysTime();
    index = TIndex{redSeqs};
    e = sysTime() - s;
    myPrint(options, 1, " done.\n");
    myPrint(options, 2, "Runtime: ", e, "s \n\n");

    return ret;
    // Dump Index
//     myPrint(options, 1, "Writing Index to disk...");
//     s = sysTime();
//     std::string _path = options.indexDir + "/index";
//
// //     index.store(_path);
//
// //TODO once sdsl-ceral-support in SeqAn3
//     {
//         std::ofstream os{_path.c_str()};
//
//         cereal::BinaryOutputArchive oarchive(os); // Create an output archive
//         oarchive(index);
//     }
//
//     e = sysTime() - s;
//     myPrint(options, 1, " done.\n");
//     myPrint(options, 2, "Runtime: ", e, "s \n");
}

#endif // header guard
