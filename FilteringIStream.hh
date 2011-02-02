#ifndef FILTERINGISTREAM_HH
#define FILTERINGISTREAM_HH

/*!
 * \file FilteringIStream.hh
 * \brief Definition of the FilteringIStream class
 */

#include "FilteringIBuf.hh"

/*!
 * \brief Proxy class for FilteringIBuf
 *
 * Class FilteringIBufUser embeds a FilteringIBuf stream buffer object.
 * It acts as a parent class for FilteringIStream, allowing it to fully
 * initialize the filtering stream buffer before constructing its
 * std::istream interface.
 * \tparam Filter The filter function to apply to the stream data
 */
template <typename Filter>
class FilteringIBufUser
{
	public:
		//! Constructor
		FilteringIBufUser(std::streambuf *sb): _buf(sb) {}
		//! Constructor
		FilteringIBufUser(std::streambuf *sb, Filter *filter):
			_buf(sb, filter) {}

		//! Return the filtering stream buffer used
		FilteringIBuf<Filter> *rdbuf() { return &_buf; }

	private:
		//! The stream buffer object
		FilteringIBuf<Filter> _buf;
};

/*!
 * \brief Filtering input stream
 *
 * Class FilteringIStream is an input stream that uses applies a filter
 * object to the incoming data before it is parsed.
 * \tparam Filter The filter function to apply to the stream data
 */
template <typename Filter>
struct FilteringIStream: private FilteringIBufUser<Filter>, public std::istream
{
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringIStream which filters stream buffer \a sb
	 * using a default constructed Filter object.
	 * \param sb The stream buffer to filter
	 */
	FilteringIStream(std::streambuf *sb):
		FilteringIBufUser<Filter>(sb),
		std::istream(FilteringIBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringIStream which filters stream buffer \a sb
	 * using Filter \a filter.
	 * \param sb         The stream buffer to filter
	 * \param filter     The filter to apply to the stream buffer
	 */
	FilteringIStream(std::streambuf *sb, Filter *filter):
		FilteringIBufUser<Filter>(sb, filter),
		std::istream(FilteringIBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringIStream which filters the stream buffer
	 * associated with input stream \a is using a default constructed
	 * Filter object.
	 * \param is The input stream to filter
	 */
	FilteringIStream(std::istream& is):
		FilteringIBufUser<Filter>(is.rdbuf()),
		std::istream(FilteringIBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringIStream which filters the stream buffer
	 * associated with input stream \a is using Filter \a filter.
	 * \param is         The input stream to filter
	 * \param filter     The filter to apply to the stream
	 */
	FilteringIStream(std::istream& is, Filter *filter):
		FilteringIBufUser<Filter>(is.rdbuf(), filter),
		std::istream(FilteringIBufUser<Filter>::rdbuf()) {}

	//! Return the filtering stream buffer
	FilteringIBuf<Filter> *rdbuf()
	{
		return FilteringIBufUser<Filter>::rdbuf();
	}
};

#endif // FILTERINGISTREAM_HH
