#ifndef INDENTFILTER_HH
#define INDENTFILTER_HH

/*!
 * \file IndentFilter.hh
 * \brief Definition of the IndentFilter class
 */

#include <iostream>

/*!
 * \brief Stream buffer filter for indenting output
 *
 * Class IndentFilter is a stream buffer filter that can be used as a template
 * argument to a FilteringOBuf stream buffer to indent output on an output
 * stream. Stream manipulators indent() and dedent() are defined for easy
 * indenting on the stream.
 */
class IndentFilter
{
	public:
		//! Local typedef for character or error code
		typedef std::streambuf::int_type int_type;
		//! Local typedef for the type of a character
		typedef std::streambuf::char_type char_type;
		//! Local typedef for character traits
		typedef std::streambuf::traits_type traits_type;

		/*!
		 * \brief Constructor
		 *
		 * Create a new IndentFilter that indents lines with indentation
		 * string \a indent, starting with an indentation level
		 * \a level.
		 * \param indent The string to indent with
		 * \param level  The initial indentation level
		 */
		IndentFilter(const std::string& indent="\t", int level=0):
			_indent(indent), _level(level), _indent_now(true) {}

		//! Increase the indentation level with one
		void indent() { _level++; }
		//! Decrease the indentation level with one
		void dedent() { if (_level > 0) _level--; }

		/*!
		 * \brief Write a character
		 *
		 * Write character \a c to output buffer \a sb. Depending
		 * on the current state of the indenter, one or more
		 * indentation strings are inserted before the characetr.
		 * \param sb The stream buffer to write to
		 * \param c The character to insert into the buffer
		 * \return eof() on error, \a c otherwise
		 */
		int_type operator()(std::streambuf *sb, int_type c);

	private:
		//! The indentation string
		std::string _indent;
		//! The current indentation level
		int _level;
		//! If \c true, insert indentation before the next character
		bool _indent_now;
};

#endif // INDENTFILTER_HH
