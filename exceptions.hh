#ifndef EXCEPTIONS_HH
#define EXCEPTIONS_HH

/*!
 * \file exceptions.hh
 * \brief Common exceptions in Quill
 */

#include <sstream>
#include <Exception.hh>

//! %Exception thrown when an element data file could not be opened
struct NoFile: public Li::Exception
{
	//! Constructor
	NoFile(const std::string& fname):
		Exception("Unable to open file \"" + fname + "\"") {}
};

//! %Exception throw when the end of the input is reached unexpectedly
struct UnexpectedEOF: public Li::Exception
{
	//! Constructor
	UnexpectedEOF(): Exception("Unexpected end of file") {}
};

/*!
 * \brief %Exception thrown when trying to look up an element not listed
 * in the periodic table
 */
struct UnknownElement: public Li::Exception
{
	//! Constructor
	UnknownElement(const std::string& elem):
		Exception("Unknown element \"" + elem + "\"") {}
	//! Constructor
	UnknownElement(int number): Exception()
	{
		std::ostringstream os;
		os << number;
		setMsg("Unknown element " + number);
	}
};

//! %Exception thrown when an errors occurs while parsing input
struct ParseError: public Li::Exception
{
	//! Constructor
	ParseError(): Exception("Parse error") {}
	//! Constructor
	ParseError(const std::string& msg): Exception("Parse error: " + msg) {}
};

/*!
 * %Exception thrown when trying to access an element in an array or map
 * using an index not present in the container.
 */
struct InvalidIndex: public Li::Exception
{
	//! Constructor
	InvalidIndex(int i): Exception()
	{
		std::ostringstream os;
		os << "Invalid index " << i;
		setMsg(os.str());
	}
};

#endif // EXCEPTIONS_HH
