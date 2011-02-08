#ifndef FILTERINGOSTREAM_HH
#define FILTERINGOSTREAM_HH

/*!
 * \file FilteringOStream.hh
 * \brief Definition of the FilteringOStream class
 */

#include "io/FilteringOBuf.hh"

/*!
 * \brief Proxy class for FilteringOBuf
 *
 * Class FilteringOBufUser embeds a FilteringOBuf stream buffer object.
 * It acts as a parent class for FilteringOStream, allowing it to fully
 * initialize the filtering stream buffer before constructing its
 * std::ostream interface.
 * \tparam Filter The filter function to apply to the stream data
 */
template <typename Filter>
class FilteringOBufUser
{
	public:
		//! Constructor
		FilteringOBufUser(std::streambuf *sb): _buf(sb) {}
		//! Constructor
		FilteringOBufUser(std::streambuf *sb, Filter *filter,
			bool own_filter):
			_buf(sb, filter, own_filter) {}

		//! Return the filtering stream buffer used
		FilteringOBuf<Filter> *rdbuf() { return &_buf; }

	private:
		//! The stream buffer object
		FilteringOBuf<Filter> _buf;
};

/*!
 * \brief Filtering output stream
 *
 * Class FilteringOStream is an output stream that uses applies a filter
 * object to the output data before it is being written.
 * \tparam Filter The filter function to apply to the stream data
 */
template <typename Filter>
struct FilteringOStream: private FilteringOBufUser<Filter>, public std::ostream
{
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringOStream which filters stream buffer \a sb
	 * using a default constructed Filter object.
	 * \param sb The stream buffer to filter
	 */
	FilteringOStream(std::streambuf *sb):
		FilteringOBufUser<Filter>(sb),
		std::ostream(FilteringOBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringOStream which filters stream buffer \a sb
	 * using Filter \a filter. If the optional parameter \a own_filter
	 * is \c true, ownership of \a filter is transferred to this
	 * object, and \a filter will be deleted when the FilteringOStream
	 * goes out of existence.
	 * \param sb         The stream buffer to filter
	 * \param filter     The filter to apply to the stream buffer
	 * \param own_filter When \c true, assume ownership of \a filter
	 */
	FilteringOStream(std::streambuf *sb, Filter *filter,
		bool own_filter=false):
		FilteringOBufUser<Filter>(sb, filter, own_filter),
		std::ostream(FilteringOBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringOStream which filters the stream buffer
	 * associated with output stream \a is using a default constructed
	 * Filter object.
	 * \param os The output stream to filter
	 */
	FilteringOStream(std::ostream& os):
		FilteringOBufUser<Filter>(os.rdbuf()),
		std::ostream(FilteringOBufUser<Filter>::rdbuf()) {}
	/*!
	 * \brief Constructor
	 *
	 * Create a new FilteringOStream which filters the stream buffer
	 * associated with output stream \a is using Filter \a filter.
	 * If the optional parameter \a own_filter is \c true, ownership of
	 * \a filter is transferred to this  object, and \a filter will
	 * be deleted when the FilteringOStream goes out of existence.
	 * \param os         The output stream to filter
	 * \param filter     The filter to apply to the stream
	 * \param own_filter When \c true, assume ownership of \a filter
	 */
	FilteringOStream(std::ostream& os, Filter *filter,
		bool own_filter=false):
		FilteringOBufUser<Filter>(os.rdbuf(), filter, own_filter),
		std::ostream(FilteringOBufUser<Filter>::rdbuf()) {}

	//! Return the filtering stream buffer
	FilteringOBuf<Filter> *rdbuf()
	{
		return FilteringOBufUser<Filter>::rdbuf();
	}
};

#endif // FILTERINGOSTREAM_HH
