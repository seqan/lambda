#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

#pragma once

namespace seqan
{

template <std::ranges::InputRange T>
    requires !std::is_lvalue_reference_v<seqan3::reference_t<T>>
inline void const * getObjectId(T const & me)
{
//     return 0;
    return std::addressof(me);
}

}

#include <seqan/basic/metaprogramming_logic.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

namespace seqan
{

template <typename t>
inline constexpr bool is_new_range = false;
template <typename t, auto j>
inline constexpr bool is_new_range<std::span<t, j>> = true;
template <typename ...ts>
inline constexpr bool is_new_range<ranges::transform_view<ts...>> = true;

template <typename t>
SEQAN3_CONCEPT NonSeqAn2Range = is_new_range<seqan3::remove_cvref_t<t>>;

template <NonSeqAn2Range T>
struct Value<T>
{
    using Type = seqan3::value_type_t<T>;
};
template <NonSeqAn2Range T>
struct Value<T const>
{
    using Type = seqan3::value_type_t<T const>;
};
template <NonSeqAn2Range T>
struct Reference<T>
{
    using Type = seqan3::reference_t<T>;
};
template <NonSeqAn2Range T>
struct Reference<T const>
{
    using Type = seqan3::reference_t<T const>;
};
template <NonSeqAn2Range T>
struct GetValue<T>
{
    using Type = seqan3::reference_t<T const>;
};
template <NonSeqAn2Range T>
struct GetValue<T const>
{
    using Type = seqan3::reference_t<T const>;
};
template <NonSeqAn2Range T>
struct Position<T>
{
    using Type = size_t;
};
template <NonSeqAn2Range T>
struct Position<T const>
{
    using Type = size_t;
};
template <NonSeqAn2Range T>
struct Size<T>
{
    using Type = size_t;
};
template <NonSeqAn2Range T>
struct Size<T const>
{
    using Type = size_t;
};
template <NonSeqAn2Range T>
struct Iterator<T, Standard>
{
    typedef Iter<T, StdIteratorAdaptor> Type;
};
template <NonSeqAn2Range T>
struct Iterator<T const, Standard>
{
    typedef Iter<T const, StdIteratorAdaptor> Type;
};
template <NonSeqAn2Range T>
struct Iterator<T, Rooted>
{
    typedef Iter<T, AdaptorIterator<Iter<T, StdIteratorAdaptor> > > Type;
};
template <NonSeqAn2Range T>
struct Iterator<T const, Rooted>
{
    typedef Iter<T const, AdaptorIterator<Iter<T const, StdIteratorAdaptor> > > Type;
};
template <NonSeqAn2Range T>
struct IsSequence<T> : True
{};
template <NonSeqAn2Range T>
struct IsSequence<T const> : True
{};
template <NonSeqAn2Range T>
struct StdContainerIterator<T>
{
     using Type = std::ranges::iterator_t<T>;
};
template <NonSeqAn2Range T>
struct StdContainerIterator<T const>
{
     using Type = std::ranges::iterator_t<T const>;
};
template <NonSeqAn2Range T>
struct Reference<Iter<T, StdIteratorAdaptor> >
{
    using Type = seqan3::reference_t<T>;
};
template <NonSeqAn2Range T>
struct Reference<Iter<T const, StdIteratorAdaptor> >
{
    using Type = seqan3::reference_t<T const>;
};
template <NonSeqAn2Range T>
struct IsContiguous<T> :
    std::conditional_t<std::ranges::ContiguousRange<T>, True, False>
{};
template <NonSeqAn2Range T>
struct HasSubscriptOperator<T> :
    std::conditional_t<std::ranges::RandomAccessRange<T>, True, False>
{};

template <NonSeqAn2Range T>
SEQAN_CONCEPT_IMPL((T), (StlContainerConcept));



// –---------------------------------------------------------------------------
// range stuff
// –---------------------------------------------------------------------------

template <std::ranges::ForwardRange TSource, typename TVal, typename TSpec>
void copy_range(TSource && in, String<TVal, TSpec> & out)
{
    if constexpr (std::ranges::SizedRange<TSource>)
    {
        resize(out, std::ranges::size(in));
        size_t i = 0;
        for (auto && v : in)
            out[i++] = std::forward<decltype(v)>(v);
    }
    else
    {
        clear(out);
        for (auto && v : in)
            appendValue(out, std::forward<decltype(v)>(v));
    }
}

// // template <typename TValue, typename TSpec, std::ranges::InputRange T>
// // inline void append(String<TValue, TSpec> & out, T const & in)
// // {
// //     if constexpr (std::ranges::SizedRange<T>)
// //         reserve(out, length(out) + std::ranges::size(in));
// //
// //     for (auto && v : in)
// //         appendValue(out, std::forward<decltype(v)>(v));
// // }
// //
// // template <typename TValue, typename TSpec, std::ranges::InputRange T>
// // inline void append(String<TValue, TSpec> & out, T & in)
// // {
// //     if constexpr (std::ranges::SizedRange<T>)
// //         reserve(out, length(out) + std::ranges::size(in));
// //
// //     for (auto && v : in)
// //         appendValue(out, std::forward<decltype(v)>(v));
// // }
// //

// –---------------------------------------------------------------------------
// alphabet stuff
// –---------------------------------------------------------------------------

template <typename stream_t, seqan3::Alphabet alph_t>
    requires !std::Integral<alph_t>
inline void write(stream_t & s, alph_t const alph)
{
    write(s, seqan3::to_char(alph));
}



template <typename alph_t>
    requires !std::Integral<alph_t> && requires (alph_t & a) { { seqan3::to_char(a) }; }
inline bool operator==(alph_t alph, char c)
{
    return seqan3::to_char(alph) == c;
}

template <typename alph_t>
    requires !std::Integral<alph_t> && requires (alph_t & a) { { seqan3::to_char(a) }; }
inline bool operator==(char c, alph_t alph)
{
    return seqan3::to_char(alph) == c;
}

template <typename TValue, typename TSequenceValue, typename TSpec, seqan3::AminoacidAlphabet alph_t>
inline auto score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec>> const & scheme, alph_t const a1, alph_t const a2)
{
    return score(scheme, AminoAcid{seqan3::to_char(a1)}, AminoAcid{seqan3::to_char(a2)});
}

// conflict with seqan::complement
template <typename TValue, typename TSpec, seqan3::NucleotideAlphabet alph_t>
inline auto score(Score<TValue, TSpec> const & scheme, alph_t const a1, alph_t const a2)
{
    return score(scheme, Iupac{seqan3::to_char(a1)}, Iupac{seqan3::to_char(a2)});
}


// –---------------------------------------------------------------------------
// random fixes
// –---------------------------------------------------------------------------




} // namespace seqan

