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

#include <seqan3/core/range/detail/misc.hpp> // seqan3::detail::consume
#include <seqan3/utility/char_operations/predicate.hpp>

inline constexpr auto not_eol = !(seqan3::is_char<'\r'> || seqan3::is_char<'\n'>);

bool setEnv(std::string const & key, std::string const & value)
{
#ifdef PLATFORM_WINDOWS
    return !_putenv_s(key.c_str(), value.c_str());
#else
    return !setenv(key.c_str(), value.c_str(), true);
#endif
}

// ----------------------------------------------------------------------------
// function _readMappingFileNCBI
// ----------------------------------------------------------------------------

template <typename TInputView, typename TStaxIDs>
void _readMappingFileUniProt(TInputView &                                      fiv,
                             TStaxIDs &                                        sTaxIds,
                             std::vector<bool> &                               taxIdIsPresent,
                             std::unordered_map<std::string, uint64_t> const & accToIdRank)
{
    // skip line with headers
    seqan3::detail::consume(fiv | std::views::take_while(not_eol));

    //TODO this is too slow, investigate whether its the lookup or the allocs
    std::string acc;
    std::string nextColumn;

    while (std::ranges::begin(fiv) != std::ranges::end(fiv))
    {
        // read accession number
        acc = fiv | std::views::take_while(!seqan3::is_blank) | seqan3::ranges::to<std::string>(); //TODO and_consume
        // skip whitespace
        seqan3::detail::consume(fiv | std::views::take_while(!seqan3::is_alnum));
        // read accession number
        nextColumn =
          fiv | std::views::take_while(!seqan3::is_blank) | seqan3::ranges::to<std::string>(); //TODO and_consume

        if ((nextColumn == "NCBI_TaxID") && (accToIdRank.count(acc) == 1))
        {
            auto & sTaxIdV = sTaxIds[accToIdRank.at(acc)];
            // skip whitespace
            seqan3::detail::consume(fiv | std::views::take_while(!seqan3::is_alnum));
            // read tax id
            nextColumn =
              fiv | std::views::take_while(!seqan3::is_space) | seqan3::ranges::to<std::string>(); //TODO and_consume

            uint32_t idNum = 0;
            auto [p, ec]   = std::from_chars(nextColumn.data(), nextColumn.data() + nextColumn.size(), idNum);
            (void)p;
            if (ec != std::errc{})
            {
                throw std::runtime_error(
                  std::string("Error: Expected taxonomical ID, but got something I couldn't read: ") + nextColumn +
                  "\n");
            }

            sTaxIdV.push_back(idNum);
            if (taxIdIsPresent.size() < idNum + 1)
                taxIdIsPresent.resize(idNum + 1);
            taxIdIsPresent[idNum] = true;
        }

        seqan3::detail::consume(fiv | std::views::take_while(not_eol));
    }
}

template <typename TInputView, typename TStaxIDs>
void _readMappingFileNCBI(TInputView &                                      fiv,
                          TStaxIDs &                                        sTaxIds,
                          std::vector<bool> &                               taxIdIsPresent,
                          std::unordered_map<std::string, uint64_t> const & accToIdRank)
{
    // skip line with headers
    seqan3::detail::consume(fiv | std::views::take_while(not_eol));

    //TODO this is too slow, investigate whether its the lookup or the allocs
    std::string buf;
    while (std::ranges::begin(fiv) != std::ranges::end(fiv))
    {
        // read accession number
        buf = fiv | std::views::take_while(!seqan3::is_blank) | seqan3::ranges::to<std::string>(); //TODO and_consume
        // we have a sequence with this ID in our database
        if (accToIdRank.count(buf) == 1)
        {
            auto & sTaxIdV = sTaxIds[accToIdRank.at(buf)];
            // skip whitespace
            seqan3::detail::consume(fiv | std::views::take_while(!seqan3::is_alnum));
            // skip versioned acc
            seqan3::detail::consume(fiv | std::views::take_while(!seqan3::is_blank));
            // skip whitespace
            seqan3::detail::consume(fiv | std::views::take_while(!seqan3::is_alnum));
            // read tax id
            buf =
              fiv | std::views::take_while(!seqan3::is_blank) | seqan3::ranges::to<std::string>(); //TODO and_consume

            uint32_t idNum = 0;
            auto [p, ec]   = std::from_chars(buf.data(), buf.data() + buf.size(), idNum);
            (void)p;
            if (ec != std::errc{})
            {
                throw std::runtime_error(
                  std::string("Error: Expected taxonomical ID, but got something I couldn't read: ") + buf + "\n");
            }
            sTaxIdV.push_back(idNum);
            if (taxIdIsPresent.size() < idNum + 1)
                taxIdIsPresent.resize(idNum + 1);
            taxIdIsPresent[idNum] = true;
        }

        seqan3::detail::consume(fiv | std::views::take_while(not_eol));
    }
}
