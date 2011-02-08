#ifndef INDENTINGOSTREAM_HH
#define INDENTINGOSTREAM_HH

/*!
 * \file IndentingOStream.hh
 * \brief Manipulators for an indenting output stream
 */

#include "io/FilteringOStream.hh"
#include "io/IndentFilter.hh"

//! Typedef for an indenting stream buffer
typedef FilteringOBuf<IndentFilter> IndentingOBuf;
//! Typedef for an indenting output stream
typedef FilteringOStream<IndentFilter> IndentingOStream;

namespace {

/*!
 * \brief Increase the indentation level
 *
 * When output stream \a os is using an IndentingOBuf for its stream
 * buffer, this stream manipulator will increase the indentation level by one.
 * \param os The output stream for which to increase indentation
 * \return The output stream with updated stream buffer
 */
inline std::ostream& indent(std::ostream& os)
{
	IndentingOBuf *buf = dynamic_cast<IndentingOBuf*>(os.rdbuf());
	if (buf) buf->filter()->indent();
	return os;
}

/*!
 * \brief Decrease the indentation level
 *
 * When output stream \a os is using an IndentingOBuf for its stream
 * buffer, this stream manipulator will decrease the indentation level by one.
 * \param os The output stream for which to decrease indentation
 * \return The output stream with updated stream buffer
 */
inline std::ostream& dedent(std::ostream& os)
{
	IndentingOBuf *buf = dynamic_cast<IndentingOBuf*>(os.rdbuf());
	if (buf) buf->filter()->dedent();
	return os;
}

} // namespace

#endif // INDENTINGOSTREAM_HH
