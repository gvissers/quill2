#ifndef FILTERINGLINEISTREAM_HH
#define FILTERINGLINEISTREAM_HH

/*!
 * \file FilteringLineIStream.hh
 * \brief Definition of the FilteringLineIStream class
 */

#include "LineIStream.hh"
#include "FilteringIStream.hh"

/*!
 * \brief Class for reading filtered lines
 *
 * Class FilteringLineIStream is an input stream that provides the ability
 * to read filtered input line by line. It provides a line-based stream
 * with the possibility of placing lines back on the stream, and can be used
 * to e.g. read job files line by line, while filtering out comments.
 * \tparam Filter The filter to apply to the input data
 */
template <typename Filter>
struct FilteringLineIStream: public IStreamUser< FilteringIStream<Filter> >,
	public LineIStream
{
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringLineIstream object for reading a job from
	 * input stream \a is.
	 * \param is The input stream to read from
	 */
	explicit FilteringLineIStream(std::istream& is):
		IStreamUser< FilteringIStream<Filter> >(is.rdbuf()),
		LineIStream(IStreamUser< FilteringIStream<Filter> >::stream())
		{}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringLineIstream object for reading a job from
	 * input stream \a is, using \a filter to filter the incoming data.
	 * \param is The input stream to read from
	 * \param filter The filter object to apply to the data from \a is.
	 */
	FilteringLineIStream(std::istream& is, Filter *filter):
		IStreamUser< FilteringIStream<Filter> >(is.rdbuf(), filter),
		LineIStream(IStreamUser< FilteringIStream<Filter> >::stream())
		{}
};

#endif // FILTERINGLINEISTREAM_HH
