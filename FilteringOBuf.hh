#ifndef FILTERINGOBUF_HH
#define FILTERINGOBUF_HH

/*!
 * \file FilteringOBuf.hh
 * \brief Definition of the FilteringOBuf clas
 */

/*!
 * \brief Class for filtering output streams
 *
 * Class FilteringOBuf is a template class for a filtering output stream
 * buffer. It is instantiated with the actual filtering function object,
 * which operates on the underlying stream buffer.
 * Based on the ideas of James Kanze, available from
 * http://lists.boost.org/Archives/boost/att-49459/fltrsbf1.htm
 * \tparam Filter The filter function to apply to the stream data
 */
template <typename Filter>
class FilteringOBuf: public std::streambuf
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new FilteringOBuf object for filtering stream
		 * buffer \a sb, using a default constructed \a Filter
		 * object as filtering function.
		 * \param sb The stream buffer to filter.
		 */
		FilteringOBuf(std::streambuf *sb):
			std::streambuf(), _sb(sb), _filter(new Filter()),
			_own_filter(true) {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new FilteringOBuf object for filtering stream
		 * buffer \a sb, using \a Filter as filtering function.
		 * When the optional parameter \a own_filter is \c true,
		 * ownership of \a filter is transferred to this object,
		 * and \a filter will be deleted when this FilteringOBuf
		 * goes out of existence.
		 * \param sb     The stream buffer to filter.
		 * \param filter The filter function to apply to \a sb.
		 * \param own_filter If \c true, transfer ownership of
		 *    \a filter to this object.
		 */
		FilteringOBuf(std::streambuf *sb, Filter *filter,
			bool own_filter=false):
			std::streambuf(), _sb(sb), _filter(filter),
			_own_filter(own_filter) {}
		//! Destructor
		virtual ~FilteringOBuf() { if (_own_filter) delete _filter; }

		/*!
		 * \brief Write new characters
		 *
		 * Pass character \a c onto the actual filtering object,
		 * which will determine what to do with it.
		 * \return The output of the filtering function. Normally
		 *    the character last written, or EOF on failure.
		 */
		virtual int_type overflow(int_type c)
		{
			return (*_filter)(_sb, c);
		}

		//! Return the filter object
		Filter* filter() { return _filter; }

	private:
		//! The stream buffer we are filtering
		std::streambuf *_sb;
		//! The filter function
		Filter *_filter;
		//! Whether we own the filter object or not
		bool _own_filter;
};

#endif // FILTERINGOBUF_HH
