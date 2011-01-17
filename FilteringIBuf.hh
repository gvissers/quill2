#ifndef FILTERINGIBUF_HH
#define FILTERINGIBUF_HH

/*!
 * \file FilteringIBuf.hh
 * \brief Definition of the FilteringIBuf clas
 */

/*!
 * \brief Class for filtering input streams
 *
 * Class FilteringIBuf is a template class for a filtering input stream
 * buffer. It is instantiated with the actual filtering function object,
 * which operates on the underlying stream buffer.
 * Based on the ideas of James Kanze, available from
 * http://lists.boost.org/Archives/boost/att-49459/fltrsbf1.htm
*/
template <typename Filter>
class FilteringIBuf: public std::streambuf
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new FilteringIBuf object for filtering stream
		 * buffer \a sb, using a default constructed \a Filter
		 * object as filtering function.
		 * \param sb The stream buffer to filter.
		 */
		FilteringIBuf(std::streambuf *sb):
			std::streambuf(), _sb(sb), _filter(new Filter()),
			_buf('\0'), _own_filter(true) {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new FilteringIBuf object for filtering stream
		 * buffer \a sb, using \a Filter as filtering function.
		 * \param sb     The stream buffer to filter.
		 * \param filter The filter function to apply to \a sb.
		 */
		FilteringIBuf(std::streambuf *sb, Filter *filter):
			std::streambuf(), _sb(sb), _filter(filter),
			_buf('\0'), _own_filter(false) {}
		//! Destructor
		virtual ~FilteringIBuf() { if (_own_filter) delete _filter; }

		/*!
		 * \brief Extract new characters
		 *
		 * Extract new characters from the source, store them in the
		 * buffer, and return the first extracted character. For this
		 * filtering buffer, this boils down to reading the next
		 * (filtered) characterd from the original stream buffer,
		 * and returning that.
		 * \return The next character extracted, or EOF when no more
		 *   characters are available
		 */
		virtual int_type underflow();

		//! Return the filter object
		Filter* filter() { return _filter; }

	private:
		//! The stream buffer we are filtering
		std::streambuf *_sb;
		//! The filter function
		Filter *_filter;
		//! Our own character buffer
		char_type _buf;
		//! Whether we own the filter object or not
		bool _own_filter;
};

template <typename Filter>
typename FilteringIBuf<Filter>::int_type FilteringIBuf<Filter>::underflow()
{
#if 0
	// Shouldn't happen, except on really broken compilers
	if (gptr() < egptr())
		return *gptr();
#endif

	int_type res = (*_filter)(_sb);
	if (res != std::char_traits<char_type>::eof())
	{
		_buf = std::char_traits<char_type>::to_char_type(res);
		setg(&_buf, &_buf, &_buf+1);
	}
	return res;
}

#endif // FILTERINGIBUF_HH
