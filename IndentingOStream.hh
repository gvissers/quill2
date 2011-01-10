#ifndef INDENTINGOSTREAM_HH
#define INDENTINGOSTREAM_HH

/*!
 * \file IndentingOStream.hh
 * \brief Definition of the IndentingOStream class
 */

#include "IndentingStreambuf.hh"

/*!
 * \brief An indenting output stream
 *
 * Class IndentingOStream is an output stream which uses an IndentingStreambuf
 * filtering stream buffer to easily indent output. Example use:
 * [code]
 * IndentingOStream os(std::cout.rdbuf());
 * os << "level 0\n" << indent << "level 1\n" << "another level 1\n"
 * 	<< indent << "level 2\n" << dedent << "yet another level 1\n"
 * 	<< dedent << "final level 0\n";
 * [/code]
 */
class IndentingOStream: public std::ostream
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new IndentingOstream which will put its filtered
		 * output on stream buffer \a buf, indenting lines with
		 * indentation string \a indent, starting at indentation level
		 * \a level.
		 * \param buf    The stream buffer to write to
		 * \param indent The indentation string
		 * \param level  The initial indentation level
		 */
		IndentingOStream(std::streambuf *buf,
			const std::string& indent="\t", int level=0):
			std::ostream(&_buf), _buf(buf, indent, level) {}

	private:
		//! The filtering stream buffer
		IndentingStreambuf _buf;
};

#endif // INDENTINGOSTREAM_HH
