#ifndef JOBISTREAM_HH
#define JOBISTREAM_HH

/*!
 * \file JobIStream.hh
 * \brief Definition of the JobIStream class
 */

#include "io/FilteringLineIStream.hh"
#include "io/JobFilter.hh"

/*!
 * \brief Class for reading job files
 *
 * Class JobIStream provides a more or less standard C++ input stream
 * interface for reading Quill job files. It provides a line-based stream
 * with the possibility of placing lines back on the stream, and automatically
 * filters out comments and empty lines.
 */
typedef FilteringLineIStream<JobFilter> JobIStream;

#endif // JOBISTREAM_HH
