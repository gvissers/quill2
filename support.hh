#ifndef SUPPORT_HH
#define SUPPORT_HH

/*!
 * \file support.hh
 * \brief Various support and macros for Quill
 */

#include <string>

#define NON_COPYABLE(class_name) \
	private:\
		class_name(const class_name&);\
		class_name& operator=(const class_name&);\

// NON_COPYABLE

/*!
 * \brief Create a string with only the first character upper case
 *
 * Create a copy of \a str, where the first character is capitalized, and
 * all following characters are lower case.
 * \param str The string to transform
 * \return string with only the first character in upper case
 */
std::string ucFirst(const std::string& str);

#endif // SUPPORT_HH
