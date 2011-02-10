#ifndef INDENTINGOSTREAM_HH
#define INDENTINGOSTREAM_HH

/*!
 * \file IndentingOStream.hh
 * \brief Typedefs for IndentingOBuf and IndentingOStream classes
 */

#include "io/FilteringOStream.hh"
#include "io/IndentFilter.hh"

//! Typedef for an indenting stream buffer
typedef FilteringOBuf<IndentFilter> IndentingOBuf;
//! Typedef for an indenting output stream
typedef FilteringOStream<IndentFilter> IndentingOStream;

#endif // INDENTINGOSTREAM_HH
