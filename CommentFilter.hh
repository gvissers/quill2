#ifndef COMMENTFILTER_HH
#define COMMENTFILTER_HH

/*!
 * \file CommentFilter.hh
 * \brief Definition of the CommentFilter class
 */

#include "IBufFilter.hh"

/*!
 * \brief Class for filtering out single-line comments
 *
 * Class CommentFilter can be used as a template parameter to a FilteringIBuf
 * object, to filter job file data from an input stream. This function
 * will strip out comments starting with a comment character that run till
 * the end of the line, and will replace commas by spaces.
 */
class CommentFilter: public IBufFilter
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new CommentFilter that filters comments starting
		 * with \a cmt_char.
		 * \param cmt_char The comment character
		 */
		CommentFilter(char cmt_char):
			IBufFilter(), _cmt_chars(1, cmt_char) {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new CommentFilter that filters comments starting
		 * with any of the characters in \a cmt_chars.
		 * \param cmt_chars The comment characters
		 */
		CommentFilter(const std::string& cmt_chars):
			IBufFilter(), _cmt_chars(cmt_chars) {}

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
		//! The comment characters
		std::string _cmt_chars;
};

#endif // COMMENTFILTER_HH
