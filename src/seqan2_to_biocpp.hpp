#pragma once

#include <ranges>
#include <span>

#include <bio/alphabet/aminoacid/all.hpp>
#include <bio/alphabet/concept.hpp>
#include <bio/alphabet/gap/all.hpp>
#include <bio/alphabet/nucleotide/concept.hpp>
#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/ranges/views/to_rank.hpp>
#include <bio/ranges/views/translate_join.hpp>

// Some things needs to be defined before including SeqAn2 headers.
// This does not make sense, but oh well ¯\_(ツ)_/¯

namespace seqan
{

template <std::ranges::input_range T>
    requires(!std::is_lvalue_reference_v<std::ranges::range_reference_t<T>>)
inline void const * getObjectId(T const & me)
{
    //     return 0;
    return std::addressof(me);
}

template <typename alph_t>
inline bio::alphabet::gapped<alph_t> gapValueImpl(bio::alphabet::gapped<alph_t> *)
{
    return bio::alphabet::gapped<alph_t>{bio::alphabet::gap{}};
}

template <typename alph_t>
inline bio::alphabet::gapped<alph_t> gapValueImpl(bio::alphabet::gapped<alph_t> const *)
{
    return bio::alphabet::gapped<alph_t>{bio::alphabet::gap{}};
}

} // namespace seqan

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

#include "holder_tristate_overload.h"
#include "bisulfite_scoring.hpp"

namespace seqan
{

// cpp20 span does not define const_iterator
#if defined(__cpp_lib_span)
template <typename TContainer>
    requires(!requires { typename TContainer::const_iterator; })
struct Value<Iter<TContainer const, StdIteratorAdaptor>>
{
    typedef TContainer const                                      TContainer_;
    typedef typename TContainer::iterator const                   TIterator_;
    typedef typename std::iterator_traits<TIterator_>::value_type Type;
};
#endif

template <typename t>
inline constexpr bool is_new_range = false;
template <typename char_t, typename traits_t>
inline constexpr bool is_new_range<std::basic_string_view<char_t, traits_t>> = true;
template <typename t, size_t j>
inline constexpr bool is_new_range<std::span<t, j>> = true;
template <typename t>
inline constexpr bool is_new_range<std::span<t, std::dynamic_extent>> = true;
template <typename... ts>
inline constexpr bool is_new_range<std::ranges::transform_view<ts...>> = true;
template <typename... ts>
inline constexpr bool is_new_range<std::ranges::subrange<ts...>> = true;
template <typename... ts>
inline constexpr bool is_new_range<bio::ranges::detail::view_transform_by_pos<ts...>> = true;
template <typename t>
inline constexpr bool is_new_range<std::ranges::take_view<t>> = true;

template <typename t>
concept NonSeqAn2Range = is_new_range<std::remove_cvref_t<t>>;

template <NonSeqAn2Range T>
struct Value<T>
{
    using Type = std::ranges::range_value_t<T>;
};
template <NonSeqAn2Range T>
struct Value<T const>
{
    using Type = std::ranges::range_value_t<T const>;
};
template <NonSeqAn2Range T>
struct Reference<T>
{
    using Type = std::ranges::range_reference_t<T>;
};
template <NonSeqAn2Range T>
struct Reference<T const>
{
    using Type = std::ranges::range_reference_t<T const>;
};
template <NonSeqAn2Range T>
struct GetValue<T>
{
    using Type = std::ranges::range_reference_t<T const>;
};
template <NonSeqAn2Range T>
struct GetValue<T const>
{
    using Type = std::ranges::range_reference_t<T const>;
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
    typedef Iter<T, AdaptorIterator<Iter<T, StdIteratorAdaptor>>> Type;
};
template <NonSeqAn2Range T>
struct Iterator<T const, Rooted>
{
    typedef Iter<T const, AdaptorIterator<Iter<T const, StdIteratorAdaptor>>> Type;
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
struct Reference<Iter<T, StdIteratorAdaptor>>
{
    using Type = std::ranges::range_reference_t<T>;
};
template <NonSeqAn2Range T>
struct Reference<Iter<T const, StdIteratorAdaptor>>
{
    using Type = std::ranges::range_reference_t<T const>;
};
template <NonSeqAn2Range T>
struct IsContiguous<T> : std::conditional_t<std::ranges::contiguous_range<T>, True, False>
{};
template <NonSeqAn2Range T>
struct HasSubscriptOperator<T> : std::conditional_t<std::ranges::random_access_range<T>, True, False>
{};

template <NonSeqAn2Range T>
SEQAN_CONCEPT_IMPL((T), (StlContainerConcept));

template <NonSeqAn2Range TContainer>
struct GetValue<Iter<TContainer, StdIteratorAdaptor>>
{
    using Type = std::ranges::range_reference_t<TContainer>;
};

template <NonSeqAn2Range TContainer>
struct GetValue<Iter<TContainer const, StdIteratorAdaptor>>
{
    using Type = std::ranges::range_reference_t<TContainer const>;
};

template <overload_c TSequence>
decltype(auto) source(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return value(gaps._source);
}

template <overload_c TSequence>
decltype(auto) source(Gaps<TSequence, ArrayGaps> & gaps)
{
    return value(gaps._source);
}

template <overload_c TSequence>
inline void assignSource(Gaps<TSequence, ArrayGaps> & gaps, TSequence const & source)
{
    setValue(gaps._source, source);
    _reinitArrayGaps(gaps);
}

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

// –---------------------------------------------------------------------------
// alphabet stuff
// –---------------------------------------------------------------------------

template <typename stream_t, bio::alphabet::alphabet alph_t>
    requires(!std::integral<alph_t>)
inline void write(stream_t & s, alph_t const alph)
{
    write(s, bio::alphabet::to_char(alph));
}

template <typename alph_t>
    requires(!std::integral<alph_t> && requires(alph_t & a) { {bio::alphabet::to_char(a)}; })
inline bool operator==(alph_t alph, char c)
{
    return bio::alphabet::to_char(alph) == c;
}

template <typename alph_t>
    requires(!std::integral<alph_t> && requires(alph_t & a) { {bio::alphabet::to_char(a)}; })
inline bool operator==(char c, alph_t alph)
{
    return bio::alphabet::to_char(alph) == c;
}

template <typename TValue, typename TSequenceValue, typename TSpec, bio::alphabet::aminoacid alph_t>
inline auto score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec>> const & scheme,
                  alph_t const                                              a1,
                  alph_t const                                              a2) noexcept
{
    return score(scheme, AminoAcid{bio::alphabet::to_char(a1)}, AminoAcid{bio::alphabet::to_char(a2)});
}

template <typename TValue, typename TSequenceValue, typename TSpec, bio::alphabet::alphabet alph_t>
inline auto score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec>> const & scheme,
                  bio::alphabet::gapped<alph_t> const                       a1,
                  bio::alphabet::gapped<alph_t> const                       a2) noexcept
{
    return score(scheme, a1.template convert_unsafely_to<alph_t>(), a2.template convert_unsafely_to<alph_t>());
}

template <typename TValue>
inline auto score(Score<TValue, ScoreMatrix<Dna5, BisulfiteMatrix>> const & scheme,
                  bio::alphabet::dna5 const                                 a1,
                  bio::alphabet::dna5 const                                 a2) noexcept
{
    return score(scheme, Dna5{bio::alphabet::to_char(a1)}, Dna5{bio::alphabet::to_char(a2)});
}

template <typename TValue, typename TSpec>
inline auto score(Score<TValue, TSpec> const & scheme,
                  bio::alphabet::dna5 const    a1,
                  bio::alphabet::dna5 const    a2) noexcept
{
    return score(scheme, Dna5{bio::alphabet::to_char(a1)}, Dna5{bio::alphabet::to_char(a2)});
}

template <typename TValue, typename TSpec>
inline auto score(Score<TValue, TSpec> const & scheme,
                  bio::alphabet::dna4 const    a1,
                  bio::alphabet::dna4 const    a2) noexcept
{
    return score(scheme, Dna{bio::alphabet::to_char(a1)}, Dna{bio::alphabet::to_char(a2)});
}

template <typename TValue, typename TSpec, bio::alphabet::alphabet alph_t>
inline auto score(Score<TValue, TSpec> const &        scheme,
                  bio::alphabet::gapped<alph_t> const a1,
                  bio::alphabet::gapped<alph_t> const a2) noexcept
{
    return score(scheme, a1.template convert_unsafely_to<alph_t>(), a2.template convert_unsafely_to<alph_t>());
}

template <bio::alphabet::alphabet alph_t>
    requires(!std::integral<alph_t>)
char convertImpl(seqan::Convert<char, alph_t>, alph_t const & a)
{
    return bio::alphabet::to_char(a);
}

template <std::integral int_t, bio::alphabet::alphabet alph_t>
    requires(!std::integral<alph_t>)
int_t convertImpl(seqan::Convert<int_t, alph_t>, alph_t const & a)
{
    return bio::alphabet::to_rank(a);
}

template <bio::alphabet::alphabet alph_t>
    requires(!std::integral<alph_t>)
alph_t convertImpl(seqan::Convert<alph_t, bio::alphabet::gapped<alph_t>>, bio::alphabet::gapped<alph_t> const & a)
{
    //     return a.template convert_unsafely_to<bio::alphabet::gap>();
    return a.template convert_unsafely_to<alph_t>();
}

template <bio::alphabet::alphabet alph_t>
    requires(!std::integral<alph_t>)
struct GappedValueType<alph_t>
{
    using Type = bio::alphabet::gapped<alph_t>;
};

// default is invalid; only used for alphabets that are scored with matrixes
template <typename seqan3_alph_t>
inline constexpr std::array<uint8_t, 1> biocpp_rank_to_seqan2rank{};

// aa27 are not:
template <>
inline constexpr std::array<uint8_t, 27> biocpp_rank_to_seqan2rank<bio::alphabet::aa27> =
  // clang-format off
//A  B  C  D  E  F  G  H  I  J  K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z   *       <- biocpp
{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25, 23, 24, 26 };
//A  B  C  D  E  F  G  H  I  J  K   L   M   N   O   P   Q   R   S   T   U   V   W   Y   Z   X   *       <- seqan2
//                                                                                  ↑   ↑   ↑
// dna5 is scored via matrix in bisulfite-mode
template <>
inline constexpr std::array<uint8_t, 5> biocpp_rank_to_seqan2rank<bio::alphabet::dna5> =
//A  C  G  N  T       <- biocpp
{ 0, 1, 2, 4, 3 };
//A  C  G  T  N       <- seqan2
//         ↑  ↑
// clang-format on

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

struct seqan2_to_rank_inner
{
    using result_type = uint8_t;

    template <typename t>
    constexpr uint8_t operator()(t const c) const noexcept
    {
        if constexpr (std::is_same_v<t, bio::alphabet::aa27> ||
                      std::is_same_v<t, bio::alphabet::dna5>) // dna5 used in bisulfite-mode
            return biocpp_rank_to_seqan2rank<t>[bio::alphabet::to_rank(c)];
        else // alphabets that are scored with seqan::SimpleScore don't need extra conversion
            return bio::alphabet::to_rank(c);
    }
};

template <typename input_t>
struct seqan2_to_rank_outer
{
    using argument_type = typename GetValue<input_t const>::Type;
    using result_type =
      seqan::ModifiedString<std::remove_reference_t<argument_type> const, seqan::ModView<seqan2_to_rank_inner>>;

    constexpr result_type operator()(argument_type const & c) const noexcept { return result_type{c}; }
};

template <typename TSimdVecs, typename TStrings>
    requires((!std::integral<bio::ranges::range_innermost_value_t<TStrings>>)&&bio::alphabet::alphabet<
             bio::ranges::range_innermost_value_t<TStrings>>)
inline void _createSimdRepImpl(TSimdVecs & simdStr, TStrings const & strings)
{
    using TModT = seqan::ModifiedString<TStrings const, seqan::ModView<seqan2_to_rank_outer<TStrings const>>>;
    _createSimdRepImpl(simdStr, TModT(strings));
}

} // namespace seqan
