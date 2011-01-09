#ifndef INDENTINGSTREAMBUF_HH
#define INDENTINGSTREAMBUF_HH

/*!
 * \file IndentingStreambuf.hh
 * \brief Definition of the IndentingStreambuf class
 */

#include <iostream>

/*!
 * \brief Stream buffer for indenting output
 *
 * Class IndentingStreambuf is a filtering stream buffer that can be used
 * together with the indent() and dedent() manipulators to indent output
 * on a stream.
 */
class IndentingStreambuf: public std::streambuf
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new IndentingStreambuf that filters output for
		 * stream buffer \a buf, indenting lines with indentation
		 * string \a indent, starting with an indentation level
		 * \a level.
		 * \param buf    The output buffer to write to
		 * \param indent The string to indent with
		 * \param level  The initial indentation level
		 */
		IndentingStreambuf(std::streambuf* buf,
			const std::string& indent="\t", int level=0):
			std::streambuf(), _buf(buf), _indent(indent),
			_level(level) {}

		//! Increase the indentation level with one
		void indent() { _level++; }
		//! Decrease the indentation level with one
		void dedent() { if (_level > 0) _level--; }

	protected:
		/*!
		 * \brief Overflow handler
		 *
		 * This function overrides the std::streambuf::overflow()
		 * function which is called when the output buffer is full.
		 * Since an IndentingStreambuf does not buffer, this function
		 * is called for every character inserted into the filter. It
		 * checks if the character is a newline, and if so, flags the
		 * buffer to insert an indentation string on the next call
		 * to overflow().
		 * \param c The character to insert into the buffer
		 * \return eof() on error, \a c otherwise
		 */
		int_type overflow(int_type c);

	private:
		//! The stream buffer for which to filter
		std::streambuf* _buf;
		//! The indentation string
		std::string _indent;
		//! The current indentation level
		int _level;
		//! If \c true, insert indentation before the next character
		bool _indent_now;
};

namespace {

/*!
 * \brief Increase the indentation level
 *
 * When output stream \a os is using an IndentingStreambuf for its stream
 * buffer, this stream manipulator will increase the indentation level by one.
 * \param os The output stream for which to increase indentation
 * \return The output stream with updated stream buffer
 */
inline std::ostream& indent(std::ostream& os)
{
	IndentingStreambuf *buf = dynamic_cast<IndentingStreambuf*>(os.rdbuf());
	if (buf) buf->indent();
	return os;
}

/*!
 * \brief Decrease the indentation level
 *
 * When output stream \a os is using an IndentingStreambuf for its stream
 * buffer, this stream manipulator will decrease the indentation level by one.
 * \param os The output stream for which to decrease indentation
 * \return The output stream with updated stream buffer
 */
inline std::ostream& dedent(std::ostream& os)
{
	IndentingStreambuf *buf = dynamic_cast<IndentingStreambuf*>(os.rdbuf());
	if (buf) buf->dedent();
	return os;
}

} // namespace

#endif // INDENTINGSTREAMBUF_HH
