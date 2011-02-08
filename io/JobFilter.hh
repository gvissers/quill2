#ifndef JOBFILTER_HH
#define JOBFILTER_HH

/*!
 * \file JobFilter.hh
 * \brief Definition of the JobFilter class
 */

#include "io/IBufFilter.hh"

/*!
 * \brief Class for filtering job files
 *
 * Class JobFilter can be used as a template parameter to a FilteringIBuf
 * object, to filter job file data from an input stream. This function
 * will strip out comments (indicated by a '#' character), and will replace
 * comma characters by a space.
 */
class JobFilter: public IBufFilter
{
	public:
		//! Constructor
		JobFilter(): _string_open('\0'), _escaped(false) {}

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
		int_type operator()(std::streambuf *sb);

	private:
		//! If not zero, the enclosing character (\c ' or \c ") of the currently open string
		char_type _string_open;
		//! In a string, when this is \c true, escape the current character and don't interpret it.
		bool _escaped;
};

#endif // JOBFILTER_HH
