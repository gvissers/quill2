#ifndef ISTREAMUSER_HH
#define ISTREAMUSER_HH

/*!
 * \file IStreamUser.hh
 * \brief Definition of the IStreamUser class
 */

#include <iostream>

/*!
 * \brief Proxy class for an input stream
 *
 * Class IStreamUser can act as a base class for classes that require
 * an embedded istream object of type \a EmbeddedIStream to be fully
 * initialized before running the contructor of another base class (e.g.
 * when this second base class needs access to the stream buffer of the
 * embedded stream). It simply embeds a \a EmbeddedIStream object, and adds
 * a few functions to access this stream.
 * \tparam EmbeddedIStream The type of the input stream to embed
 */
template <typename EmbeddedIStream>
class IStreamUser
{
	protected:
		/*!
		 * \brief Constructor
		 *
		 * Create a new IStreamUser object, embedding a default
		 * constructed EmbeddedIStream.
		 */
		IStreamUser(): _is() {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new IStreamUser object, embedding an
		 * EmbeddedIstream which uses \a sb for its stream buffer.
		 * \param sb The stream buffer of the embedded stream
		 */
		IStreamUser(std::streambuf* sb): _is(sb) {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new IStreamUser object, embedding a filtering
		 * EmbeddedIstream which uses \a sb for its stream buffer
		 * and filter object \a filter. The embedded type
		 * \a EmbeddedIStream should be a FilteringIStream.
		 * \param sb     The stream buffer of the embedded stream
		 * \param filter The filter object to apply to apply to
		 *    incoming data.
		 */
		template <typename Filter>
		IStreamUser(std::streambuf* sb, Filter *filter):
			_is(sb, filter) {}

		//! Return the embedded \c istringstream
		EmbeddedIStream& stream() { return _is; }
		//! Return the embedded \c istringstream
		const EmbeddedIStream& stream() const { return _is; }

		//! Set the status flags of the embedded stream to \a state
		void clear(std::ios::iostate state=std::ios::goodbit)
		{
			_is.clear(state);
		}

		//! Return the string buffer of the embedded stream
		std::stringbuf *rdbuf() { return _is.rdbuf(); }

	private:
		//! The embedded stream object
		EmbeddedIStream _is;
};

#endif // ISTREAMUSER_HH
