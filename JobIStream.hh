#ifndef JOBISTREAM_HH
#define JOBISTREAM_HH

/*!
 * \file JobIStream.hh
 * \brief Definition of the JobIStream class
 */

#include <list>
#include <sstream>
#include <Exception.hh>
#include "FilteringIStream.hh"
#include "JobFilter.hh"

/*!
 * \brief Proxy class for std::istringstream
 *
 * Class IStringStreamUser can act as a base class for classes that require
 * a \c std::istringstream to be fully initialized before running the
 * contructor of another base class (e.g. when this second base class needs
 * access to the stream buffer of the embedded stream). It simply embeds a
 * \c std::istringstream object, and adds a few functions to access this
 * embedded stream.
 */
class IStringStreamUser
{
	protected:
		//! Constructor
		IStringStreamUser(): _iss() {}

		//! Return the embedded \c istringstream
		std::istringstream& stringstream() { return _iss; }

		//! Set the status flags of the embedded stream to \a state
		void clear(std::ios::iostate state=std::ios::goodbit)
		{
			_iss.clear(state);
		}

		//! Return the string currently associated with the stream
		std::string str() const { return _iss.str(); }
		//! Set the string currently associated with the stream
		void str(const std::string& s) { _iss.str(s); }

	private:
		//! The embedded stream object
		std::istringstream _iss;
};

/*!
 * \brief Class for reading job files
 *
 * Class JobIStream provides a more or less standard C++ input stream
 * interface for reading Quill job files. It provides a line-based stream
 * with the possibility of placing lines back on the stream, and automatically
 * filters out comments and empty lines.
 */
class JobIStream: public IStringStreamUser, public std::istream
{
	public:
		struct ReadFailure;

		/*!
		 * \brief Constructor
		 *
		 * Create a new JobIstream object for rading a job from
		 * input stream \a is.
		 */
		JobIStream(std::istream& is):
			IStringStreamUser(),
			std::istream(stringstream().rdbuf()),
			_lines(), _input(is) {}

		//! Return the line currently being processed
		std::string line() const { return str(); }

		//! Set the status flags of this stream to \a state
		void clear(iostate state=goodbit)
		{
			IStringStreamUser::clear(state);
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
		 */
		void ungetLine(const std::string& line)
		{
			_lines.push_back(line);
		}

	private:
		//! Buffer of lines
		std::list<std::string> _lines;
		//! The (filtered) stream from which we read lines
		FilteringIStream<JobFilter> _input;

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
struct JobIStream::ReadFailure: public Li::Exception
{
	//! Constructor
	ReadFailure(): Li::Exception("Failed to read a job line") {}
};

namespace {

/*!
 * \brief Apply a stream manipulator
 *
 * Apply stream manipulator \a f to input stream \a jis. This function is
 * a specialization for a JobIStream of the extraction operator for stream
 * manipulators.
 * \param jis The JobIStream to manipulate
 * \param f   The function to apply to \a jis
 * \return The result of \a f.
 */
inline JobIStream& operator>>(JobIStream& jis, JobIStream& (*f)(JobIStream&))
{
	return f(jis);
}

/*!
 * \brief Retrieve a new line from a job input stream
 *
 * Stream manipulator for line extraction for a JobIStream. In principle,
 * one could do things like (disregarding error checking)
 * \code
 * int i;
 * std::string s;
 * jis >> getline >> i >> s;
 * \endcode
 * to extract an integer and  string from the next line in input stream \a jis.
 * \param jis The JobIStream from which to extract the next line
 * \return The update input stream
 */
inline JobIStream& getline(JobIStream& jis)
{
	jis.getLine();
	return jis;
}

} //namespace

#endif // JOBISTREAM_HH
