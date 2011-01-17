#ifndef FILTERINGISTREAM_HH
#define FILTERINGISTREAM_HH

#include "FilteringIBuf.hh"

/*!
 * \brief Proxy class for FilteringIBuf
 *
 * Class FilteringIBufUser embeds a FilteringIBuf stream buffer object.
 * It acts as a parent class for FilteringIStream, aloowing it to fully
 * initialize the filtering stream buffer before constructing its
 * std::istream interface.
 */
template <typename Filter>
class FilteringIBufUser
{
	public:
		//! Constructor
		FilteringIBufUser(std::streambuf *sb): _buf(sb) {}
		//! Constructor
		FilteringIBufUser(std::streambuf *sb, Filter *filter,
			bool own_filter):
			_buf(sb, filter, own_filter) {}

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
	 * using Filter \a filter. If the optional parameter \a own_filter
	 * is \c true, ownership of \a filter is transferred to this
	 * object, and \a filter will be deleted when the FilteringIStream
	 * goes out of existence.
	 * \param sb         The stream buffer to filter
	 * \param filter     The filter to apply to the stream buffer
	 * \param own_filter When \c true, assume ownership of \a filter
	 */
	FilteringIStream(std::streambuf *sb, Filter *filter,
		bool own_filter=false):
		FilteringIBufUser<Filter>(sb, filter, own_filter),
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
	 * If the optional parameter \a own_filter is \c true, ownership of
	 * \a filter is transferred to this  object, and \a filter will 
	 * be deleted when the FilteringIStream goes out of existence.
	 * \param is         The input stream to filter
	 * \param filter     The filter to apply to the stream
	 * \param own_filter When \c true, assume ownership of \a filter
	 */
	FilteringIStream(std::istream& is, Filter *filter,
		bool own_filter=false):
		FilteringIBufUser<Filter>(is.rdbuf(), filter, own_filter),
		std::istream(FilteringIBufUser<Filter>::rdbuf()) {}

	//! Return the filtering stream buffer
	FilteringIBuf<Filter> *rdbuf()
	{
		return FilteringIBufUser<Filter>::rdbuf();
	}
};

#endif // FILTERINGISTREAM_HH
