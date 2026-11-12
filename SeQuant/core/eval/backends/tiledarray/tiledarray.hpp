#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_TILEDARRAY_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_TILEDARRAY_HPP

/// @file tiledarray.hpp
/// Includes TiledArray (and hence MADNESS) with default symbol visibility.
///
/// The SeQuant libraries are compiled with hidden visibility, which would also
/// apply to the TiledArray/MADNESS code they instantiate from the headers.
/// TiledArray and MADNESS rely on default visibility, though: they keep
/// process-wide state in statics of inline functions (e.g. the default World)
/// and compare type_info objects, which would then be duplicated per shared
/// library. All SeQuant code must include TiledArray via this header.

#ifdef SEQUANT_HAS_TILEDARRAY

#if defined(__GNUC__) && !(defined(_WIN32) || defined(__CYGWIN__))
#pragma GCC visibility push(default)
#endif

#include <TiledArray/einsum/tiledarray.h>
#include <tiledarray.h>

#if defined(__GNUC__) && !(defined(_WIN32) || defined(__CYGWIN__))
#pragma GCC visibility pop
#endif

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_TILEDARRAY_HPP
