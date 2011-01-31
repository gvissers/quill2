#ifndef LINEISTREAM_HH
#define LINEISTREAM_HH

/*!
 * \file LineIStream.hh
 * \brief Definition of the LineIStream class
 */

#include <list>
#include <sstream>
#include <Exception.hh>
#include "IStreamUser.hh"

/*!
 * \brief Class for reading files line by line
 *
 * Class LineIStream provides a more or less standard C++ input stream
 * interface for reading (Quill job) files line by line. It provides a
 * line-based stream with the possibility of placing lines back on the
 * stream.
 */
class LineIStream: public IStreamUser<std::istringstream>, public std::istream
{
	public:
		struct ReadFailure;

		/*!
		 * \brief Constructor
		 *
		 * Create a new LineIstream object for reading from
		 * input stream \a is.
		 */
		explicit LineIStream(std::istream& is):
			IStreamUser<std::istringstream>(),
			std::istream(rdbuf()),
			_lines(), _input(is) {}

		//! Return the line currently being processed
		std::string line() const { return stream().str(); }

		//! Set the status flags of this stream to \a state
		void clear(iostate state=goodbit)
		{
			IStreamUser<std::istringstream>::clear(state);
			std::istream::clear(state);
		}

		/*!
		 * \brief Extract a line
		 *
		 * Extract a new line from the line buffer, or when this is
		 * empty, from the original input stream. The extracted line
		 * is used a string buffer for extraction. When no more
		 * lines are available, \c eofbit is set.
		 * \exception ReadFailure thrown when we fail to extract
		 *    a new line.
		 */
		void getLine();
		/*!
		 * \brief Place the last extracted line back in the stream
		 *
		 * Place the line that was last extracted back in the stream.
		 * the next call to getLine() will extract the same line again.
		 */
		void ungetLastLine()
		{
			ungetLine(line());
		}
		/*!
		 * \brief Place a line back in the stream
		 *
		 * Place \a line back in the stream. When extracting new line
		 * using getLine(), this line will be used before a new line
		 * is extracted from the original input stream.
		 * \param line The line to be placed back in the stream
		 */
		void ungetLine(const std::string& line)
		{
			_lines.push_back(line);
		}

		//! Return the string buffer associated with the current line
		std::stringbuf* rdbuf() { return IStreamUser<std::istringstream>::rdbuf(); }

	private:
		//! Buffer of lines
		std::list<std::string> _lines;
		//! The stream from which we read lines
		std::istream& _input;

		/*!
		 * \brief Extract a new line
		 *
		 * Get the next line from the original input stream, and
		 * return it. If no more lines are available, an empty
		 * string is returned.
		 * \exception ReadFailure thrown when we fail to read from
		 *    the input stream.
		 */
		std::string nextLine();
};

//! %Exception thrown when we fail to read a line from the input stream
struct LineIStream::ReadFailure: public Li::Exception
{
	//! Constructor
	ReadFailure(): Li::Exception("Failed to read a input line") {}
};

namespace {

/*!
 * \brief Apply a stream manipulator
 *
 * Apply stream manipulator \a f to input stream \a jis. This function is
 * a specialization for a LineIStream of the extraction operator for stream
 * manipulators.
 * \param jis The LineIStream to manipulate
 * \param f   The function to apply to \a jis
 * \return The result of \a f.
 */
inline LineIStream& operator>>(LineIStream& jis,
	LineIStream& (*f)(LineIStream&))
{
	return f(jis);
}

/*!
 * \brief Retrieve a new line from a line input stream
 *
 * Stream manipulator for line extraction for a LineIStream. In principle,
 * one could do things like (disregarding error checking)
 * \code
 * int i;
 * std::string s;
 * jis >> getline >> i >> s;
 * \endcode
 * to extract an integer and  string from the next line in input stream \a jis.
 * \param jis The LineIStream from which to extract the next line
 * \return The update input stream
 */
inline LineIStream& getline(LineIStream& jis)
{
	jis.getLine();
	return jis;
}

} //namespace

#endif // LINEISTREAM_HH
