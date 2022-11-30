// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2020, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2020, Knut Reinert and Freie Universität Berlin
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
// mkindex_misc.hpp: misc stuff for indexer
// ==========================================================================

#pragma once

#include <charconv>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <bio/io/txt/reader.hpp>
#include <bio/ranges/to.hpp>

bool setEnv(std::string const & key, std::string const & value)
{
#ifdef PLATFORM_WINDOWS
    return !_putenv_s(key.c_str(), value.c_str());
#else
    return !setenv(key.c_str(), value.c_str(), true);
#endif
}

// ----------------------------------------------------------------------------
// hashmap type
// ----------------------------------------------------------------------------

struct hash_string
{
    using is_transparent = void;

    std::size_t operator()(std::string const & v) const { return std::hash<std::string>{}(v); }
    std::size_t operator()(char const * v) const { return std::hash<std::string_view>{}(v); }

    std::size_t operator()(std::string_view const & v) const { return std::hash<std::string_view>{}(v); }
};

using TaccToIdRank = std::unordered_map<std::string, uint64_t, hash_string, std::equal_to<void>>;

// ----------------------------------------------------------------------------
// function _readMappingFileNCBI
// ----------------------------------------------------------------------------

#if defined(__GNUC__) && __GNUC__ < 11
#    define LAMBDA_TO_STR(x) static_cast<std::string>(x)
#else
#    define LAMBDA_TO_STR(x) x
#endif

inline void _readMappingFileUniProt(std::filesystem::path const &        fileName,
                                    TaccToIdRank const &                 accToIdRank,
                                    std::vector<std::vector<uint32_t>> & sTaxIds,
                                    std::vector<bool> &                  taxIdIsPresent)
{
    bio::io::txt::reader reader{fileName, '\t'};

    for (auto & r : reader)
    {
        assert(r.fields.size() == 3);

        std::string_view acc      = r.fields[0];
        std::string_view category = r.fields[1];
        std::string_view tId      = r.fields[2];

        if (category == "NCBI_TaxID")
        {
            if (auto it = accToIdRank.find(LAMBDA_TO_STR(acc)); it != accToIdRank.end())
            {
                std::vector<uint32_t> & sTaxIdV = sTaxIds[it->second];

                uint32_t idNum = 0;
                auto [p, ec]   = std::from_chars(tId.data(), tId.data() + tId.size(), idNum);
                if (ec != std::errc{})
                {
                    std::string msg = "Error: Expected taxonomical ID, but got something I couldn't read: ";
                    msg += tId;
                    msg += '\n';
                    throw std::runtime_error(msg);
                }
                sTaxIdV.push_back(idNum);
                if (taxIdIsPresent.size() < idNum + 1)
                    taxIdIsPresent.resize(idNum + 1);
                taxIdIsPresent[idNum] = true;
            }
        }
    }
}

inline void _readMappingFileNCBI(std::filesystem::path const &        fileName,
                                 TaccToIdRank const &                 accToIdRank,
                                 std::vector<std::vector<uint32_t>> & sTaxIds,
                                 std::vector<bool> &                  taxIdIsPresent)
{
    bio::io::txt::reader reader{fileName, '\t', bio::io::txt::header_kind::first_line};

    if (reader.header() != "accession\taccession.version\ttaxid\tgi")
        throw std::runtime_error{"Unexpected first line in NCBI taxid file."};

    for (auto & r : reader)
    {
        assert(r.fields.size() == 4);

        std::string_view acc = r.fields[0];
        std::string_view tId = r.fields[2];

        if (auto it = accToIdRank.find(LAMBDA_TO_STR(acc)); it != accToIdRank.end())
        {
            std::vector<uint32_t> & sTaxIdV = sTaxIds[it->second];

            uint32_t idNum = 0;
            auto [p, ec]   = std::from_chars(tId.data(), tId.data() + tId.size(), idNum);
            if (ec != std::errc{})
            {
                std::string msg = "Error: Expected taxonomical ID, but got something I couldn't read: ";
                msg += tId;
                msg += '\n';
                throw std::runtime_error(msg);
            }
            sTaxIdV.push_back(idNum);
            if (taxIdIsPresent.size() < idNum + 1)
                taxIdIsPresent.resize(idNum + 1);
            taxIdIsPresent[idNum] = true;
        }
    }
}

#undef LAMBDA_TO_STR
