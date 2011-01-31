#ifndef JOBISTREAM_HH
#define JOBISTREAM_HH

/*!
 * \file JobIStream.hh
 * \brief Definition of the JobIStream class
 */

#include "LineIStream.hh"
#include "FilteringIStream.hh"
#include "JobFilter.hh"

/*!
 * \brief Class for reading job files
 *
 * Class JobIStream provides a more or less standard C++ input stream
 * interface for reading Quill job files. It provides a line-based stream
 * with the possibility of placing lines back on the stream, and automatically
 * filters out comments and empty lines.
 */
struct JobIStream: public IStreamUser< FilteringIStream<JobFilter> >,
	public LineIStream
{
	/*!
	 * \brief Constructor
	 *
	 * Create a new JobIstream object for reading a job from
	 * input stream \a is.
	 */
	explicit JobIStream(std::istream& is):
		IStreamUser< FilteringIStream<JobFilter> >(is.rdbuf()),
		LineIStream(IStreamUser< FilteringIStream<JobFilter> >::stream())
		{}
};

#endif // JOBISTREAM_HH
