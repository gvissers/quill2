#ifndef IBUFFILTER_HH
#define IBUFFILTER_HH

/*!
 * \file IBufFilter.hh
 * \brief Definition of the IBufFilter class
 */

#include <iostream>

/*!
 * \brief Base class for stream buffer filters
 *
 * Class IBufFilter is a base class for input stream buffer filters. It
 * defines the interface that the filters should adhere to, as well as a few
 * common type definitions.
 */
struct IBufFilter
{
	//! Local typedef for character or error code
	typedef std::streambuf::int_type int_type;
	//! Local typedef for the type of a character
	typedef std::streambuf::char_type char_type;
	//! Local typedef for character traits
	typedef std::streambuf::traits_type traits_type;

	//! Destructor
	virtual ~IBufFilter() {}

	/*!
	 * \brief Apply the filter
	 *
	 * Apply this filter object to stream buffer \a sb, and return
	 * the first filtered character.
	 * available.
	 * \param sb The stream buffer to filter
	 * \return The first filtered character, or EOF when no more
	 *   data is available.
	 */
	virtual int_type operator()(std::streambuf *sb) = 0;
};

#endif // IBUFFILTER_HH
