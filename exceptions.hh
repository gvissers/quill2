#ifndef EXCEPTIONS_HH
#define EXCEPTIONS_HH

/*!
 * \file exceptions.hh
 * \brief Common exceptions in Quill
 */

#include <Exception.hh>

/*!
 * \brief %Exception thrown when an element data file could not be opened
 */
struct NoFile: public Li::Exception
{
	//! Constructor
	NoFile(const std::string& fname):
		Exception("Unable to open file \"" + fname + "\"") {}
};

#endif // EXCEPTIONS_HH
