#ifndef LINEGETTER_HH
#define LINEGETTER_HH

/*!
 * \file LineGetter.hh
 * \brief Definition of the LineGetter class
 */

#include <iostream>
#include <Exception.hh>
#include "support.hh"

/*!
 * \brief Class for extracting lines
 *
 * Class LineGetter is used to read single lines from an input stream and
 * return them, ignoring empty or comment lines.
 */
class LineGetter
{
	NON_COPYABLE(LineGetter)

	public:
		struct ReadFailure;
		struct UnexpectedEOF;

		/*!
		 * \brief Constructor
		 *
		 * Create a new LineGetter extracting lines from input stream
		 * \a is, using comment characters \a cmt, starting with line
		 * number \a line, and having a text buffer of \a buf_size
		 * characters.
		 * \param is       The input stream to read from
		 * \param cmt      The comment characters. Lines starting with
		 *    one of these characters are ignored.
		 * \param line_nr  The initial line number
		 * \param buf_size Size of the line buffer
		 */
		LineGetter(std::istream& is, const std::string& cmt="",
			int line_nr=0, size_t buf_size=256):
			_is(is), _buf(new char[buf_size]),
			_buf_size(buf_size), _line_nr(line_nr), _cmt(cmt) {}
		//! Destructor
		~LineGetter() { delete[] _buf; }

		//! Return The input stream we are reading from
		std::istream& stream() { return _is; }
		//! Return the current line number
		int lineNumber() const { return _line_nr; }
		//! Return the current line
		const char* line() const { return _buf; }

		/*!
		 * \brief Read a line
		 *
		 * Read the next non-empty, non-comment line from the input
		 * stream, and return it. If the optional argument \a except
		 * is \c true (the default), an UnexpectedEOF exception is
		 * thrown when end of file is reached before an non-empty line
		 * was read. If \a except is \c false, 0 will be returned
		 * in this case.
		 * \param except Whether to throw an exception when
		 *    end of file is reached
		 * \exception UnexpectedEOF thrown when \a except is \c true,
		 *    and end of file is reached before a line was found.
		 * \note This function returns a pointer to an internal
		 * buffer is, which will be overwritten on the next
		 * call to next().
		 */
		const char* next(bool except=true);

		//! Set the comment characters to \a cmt
		void setComment(const std::string& cmt)
		{
			_cmt = cmt;
		}
		/*!
		 * \brief Change the buffer size
		 *
		 * Change the size of the text buffer to \a size characters.
		 * If \a size is less than the current buffer size, the
		 * buffer will be truncated.
		 * \param size The new size of the buffer.
		 */
		void setBufferSize(size_t size);

	private:
		//! The input stream to read from
		std::istream& _is;
		//! The data buffer holding the current line
		char *_buf;
		//! Size of the data buffer
		size_t _buf_size;
		//! The current line number
		int _line_nr;
		//! The comment character
		std::string _cmt;
};

//! %Exception thrown when we fail to read from an input stream
struct LineGetter::ReadFailure: public Li::Exception
{
	//! Constructor
	ReadFailure(): Exception("Read failure") {}
};

//! %Exception thrown when end of file was reached unexepectedly
struct LineGetter::UnexpectedEOF: public Li::Exception
{
	//! Constructor
	UnexpectedEOF(): Exception("unexpected end of file") {}
};

#endif // LINEGETTER_HH
