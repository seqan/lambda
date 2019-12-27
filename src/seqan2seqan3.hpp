#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/range/views/translate_join.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

#pragma once

namespace seqan
{

template <std::ranges::input_range T>
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
template <typename char_t, typename traits_t>
inline constexpr bool is_new_range<std::basic_string_view<char_t, traits_t>> = true;
template <typename t, auto j>
inline constexpr bool is_new_range<std::span<t, j>> = true;
template <typename t>
inline constexpr bool is_new_range<std::span<t, -1>> = true;
template <typename ...ts>
inline constexpr bool is_new_range<ranges::transform_view<ts...>> = true;
template <typename ...ts>
inline constexpr bool is_new_range<ranges::subrange<ts...>> = true;
template <typename ...ts>
inline constexpr bool is_new_range<seqan3::detail::view_translate_join<ts...>> = true;
template <typename ...ts>
inline constexpr bool is_new_range<seqan3::detail::view_translate_single<ts...>> = true;
template <typename t, bool b1, bool b2>
inline constexpr bool is_new_range<seqan3::detail::view_take<t, b1, b2>> = true;

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
    std::conditional_t<std::ranges::contiguous_range<T>, True, False>
{};
template <NonSeqAn2Range T>
struct HasSubscriptOperator<T> :
    std::conditional_t<std::ranges::random_access_range<T>, True, False>
{};

template <NonSeqAn2Range T>
SEQAN_CONCEPT_IMPL((T), (StlContainerConcept));

template <NonSeqAn2Range TContainer>
struct GetValue<Iter<TContainer, StdIteratorAdaptor> >
{
    using Type = seqan3::reference_t<TContainer>;
};

template <NonSeqAn2Range TContainer>
struct GetValue<Iter<TContainer const, StdIteratorAdaptor> >
{
    using Type = seqan3::reference_t<TContainer const>;
};

// –---------------------------------------------------------------------------
// range stuff
// –---------------------------------------------------------------------------

template <std::ranges::forward_range TSource, typename TVal, typename TSpec>
void copy_range(TSource && in, String<TVal, TSpec> & out)
{
    if constexpr (std::ranges::sized_range<TSource>)
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

// // template <typename TValue, typename TSpec, std::ranges::input_range T>
// // inline void append(String<TValue, TSpec> & out, T const & in)
// // {
// //     if constexpr (std::ranges::sized_range<T>)
// //         reserve(out, length(out) + std::ranges::size(in));
// //
// //     for (auto && v : in)
// //         appendValue(out, std::forward<decltype(v)>(v));
// // }
// //
// // template <typename TValue, typename TSpec, std::ranges::input_range T>
// // inline void append(String<TValue, TSpec> & out, T & in)
// // {
// //     if constexpr (std::ranges::sized_range<T>)
// //         reserve(out, length(out) + std::ranges::size(in));
// //
// //     for (auto && v : in)
// //         appendValue(out, std::forward<decltype(v)>(v));
// // }
// //

// –---------------------------------------------------------------------------
// alphabet stuff
// –---------------------------------------------------------------------------

template <typename stream_t, seqan3::alphabet alph_t>
    requires !std::integral<alph_t>
inline void write(stream_t & s, alph_t const alph)
{
    write(s, seqan3::to_char(alph));
}



template <typename alph_t>
    requires !std::integral<alph_t> && requires (alph_t & a) { { seqan3::to_char(a) }; }
inline bool operator==(alph_t alph, char c)
{
    return seqan3::to_char(alph) == c;
}

template <typename alph_t>
    requires !std::integral<alph_t> && requires (alph_t & a) { { seqan3::to_char(a) }; }
inline bool operator==(char c, alph_t alph)
{
    return seqan3::to_char(alph) == c;
}

template <typename TValue, typename TSequenceValue, typename TSpec, seqan3::aminoacid_alphabet alph_t>
inline auto score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec>> const & scheme, alph_t const a1, alph_t const a2)
{
    return score(scheme, AminoAcid{seqan3::to_char(a1)}, AminoAcid{seqan3::to_char(a2)});
}

// conflict with seqan::complement
template <typename TValue, typename TSpec, seqan3::nucleotide_alphabet alph_t>
inline auto score(Score<TValue, TSpec> const & scheme, alph_t const a1, alph_t const a2)
{
    return score(scheme, Iupac{seqan3::to_char(a1)}, Iupac{seqan3::to_char(a2)});
}

// –---------------------------------------------------------------------------
// SIMD Support
// SeqAn2 converts input sequence values into SIMD lanes via implicit conversion.
// This works for SeqAn2 alphabets, because they are implicitly convertible
// to all integral types but fails for SeqAn3 alphabets that are not.
//
// Below is an overload of the SIMD representation conversion function in SeqAn2.
// It gets chosen for ranges-over-ranges-over-SeqAn3-alphabet.
// And it applies a two-dimensional SeqAn2 View on top of that range-of-range
// that lazily converts the SeqAn3 alphabet to its rank representation which
// can then be used to create the SIMD lanes.
//
// –---------------------------------------------------------------------------

#ifdef SEQAN_SIMD_ENABLED
struct seqan2_to_rank_inner
{
    using result_type = uint8_t;

    template <typename t>
    constexpr uint8_t operator()(t const c) const noexcept
    {
        return seqan3::to_rank(c);
    }
};

template <typename input_t>
struct seqan2_to_rank_outer
{
    using argument_type = typename GetValue<input_t const>::Type;
    using result_type = seqan::ModifiedString<std::remove_reference_t<argument_type> const,
                                              seqan::ModView<seqan2_to_rank_inner>>;

    constexpr result_type operator()(argument_type const & c) const noexcept
    {
        return result_type{c};
    }
};

template <typename TSimdVecs,
          typename TStrings>
    requires (!std::integral<decltype(std::declval<TStrings&>()[0][0])>) &&
             seqan3::alphabet<decltype(std::declval<TStrings&>()[0][0])>
inline void
_createSimdRepImpl(TSimdVecs & simdStr,
                   TStrings const & strings)
{
    using TModT   = seqan::ModifiedString<TStrings const, seqan::ModView<seqan2_to_rank_outer<TStrings const>>>;
    _createSimdRepImpl(simdStr, TModT(strings));
}
#endif // SEQAN_SIMD_ENABLED

} // namespace seqan
